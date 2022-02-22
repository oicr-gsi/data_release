# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 12:27:07 2020

@author: rjovelin
"""


import argparse
import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import matplotlib as mpl
from matplotlib.lines import Line2D
from datetime import datetime
import xhtml2pdf.pisa as pisa
import mistune
import base64
import copy
from PIL import Image
import math
import requests
import gzip
from itertools import zip_longest


def get_libraries(library_file):
    '''
    (str) -> dict
    
    Returns a dictionary with library alias as key and library aliquot ID as value if 
    provided in library_file or the empty string
    
    Parameters
    ----------
    - sample_file (str): Path to sample file. Sample file is a tab delimited file
                         that includes 1 or 2 columns. The first column is always 
                         the library alias (TGL17_0009_Ct_T_PE_307_CM).
                         The second and optional column is the library aliquot ID (eg. LDI32439)
    '''
    D = {}
    
    infile = open(library_file)
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) == 1:
                D[line[0]] = ''
            elif len(line) == 2:
                D[line[0]] = line[1]
    infile.close()
    return D
    

def get_workflow_id(file_path):
    '''
    (str) -> int
    
    Parameters
    ----------
    
    Returns the workflow id from the file path
           
    - file_path (str): File path
    '''
    
    workflow = -1
    
    k = file_path.split('/')
    for j in k:
        if j.isdigit():
            workflow = int(j)
    return workflow
    

def select_most_recent_workflow(L):
    '''
    (list) -> list
    
    Parameters
    ----------
    
    Returns a list of file paths corresponding to the most recent workflow ID for each file
    


    
    - L (list): List of file paths
    '''
    
    # create a dict to store workflow accession and file path for each file name
    d = {}
    
    for i in L:
        filename = os.path.basename(i)
        if filename not in d:
            d[filename] = []
        workflow = get_workflow_id(i)
        d[filename].append([workflow, i])
    # get the most recent workflow accession id and corresponding file
    for i in d:
        d[i].sort()
    # keep only the files corresponding to the most recent workflow run
    # when worflow id not in file path, eg FileImport, just select the file. assumes 1 record in FPR
    files = [d[i][-1][1] for i in d]
    return files


def extract_files(provenance, project, runs, workflow, nomiseq, library_aliases, files_release, exclude):
    '''
    (str, str, list, str, bool, dict, list, list) -> (dict, dict, dict)
  
    Returns a tuple with dictionaries with files extracted from FPR and their corresponding run
    respectively from release and withheld from release, and a dictionary with md5sums of the release files.
    Returns the files corresponding to the most recent workflow iteration if duplicate files
            
    Parameters
    ----------
    - provenance (str): Path to File Provenance Report
    - project (str): Project name as it appears in File Provenance Report or empty string.
                     Used to parse the FPR by project. Files are further filtered
                     by run is runs parameter if provided, or all files for
                     the project and workflow are used
    - runs (list): List of run IDs or empty list. Other runs are ignored if not empty    
    - workflow (str): Worflow used to generate the output files
    - nomiseq (bool): Exclude MiSeq runs if True
    - library_aliases (dict): Dictionary with library alias as key and library aliquot ID as value or the empty string.
                              Can be an empty dictionary
    - files_release (list): List of file names to release
    - exclude (list): List of samples or libraries to exclude from release
    '''
    
    # create a dict {run: [files]}
    D, K = {}, {}
    
    # create a dict {file: md5sum}
    M = {}
        
    # make a list of records
    records = []
    
    # create a list to store project runs
    project_runs = []
    
    
    # parse the file provenance record
    if project: 
        records.extend(get_FPR_records(project, provenance, 'project'))
    elif runs:
        for run in runs:
            records.extend(get_FPR_records(run, provenance, 'run'))

    # parse the records
    for i in records:
        i = i.rstrip().split('\t')
        # get file path
        file_path = i[46]
        # get sample name
        sample_name = i[7]
        # get parent sample name
        parent_sample = i[9].split(':')[0]
        # skip if project provided and not in record
        if project and project != i[1]:
            continue
        # skip miseq runs if miseq should be excluded
        if nomiseq and 'miseq' in i[22].lower():
            continue
        # get run id
        run_id = i[18]
        # parse files from all runs if runs not provided
        if len(runs) == 0:
            project_runs.append(run_id)
        else:
            project_runs = runs
        # skip runs that are not specified  
        if run_id not in project_runs:
            continue
        # get library aliases and aliquot
        aliquot = i[56].split('_')[-1]
        library = i[13]
        # skip if incorrect workflow
        if workflow == 'bcl2fastq':
            # skip if not casava or bcl2fastq workflow    
            if not (i[30].lower() == 'casava' or i[30].lower() == 'bcl2fastq'):
                continue
        else:
            # skip if not provided workflow
            if workflow.lower() != i[30].lower():
                continue
        # check if file list is provided
        if files_release:
            # skip files not in the allow list
            if os.path.basename(file_path) not in files_release:
                continue
        # check if library aliases are provided
        if library_aliases:
            # skip libraries not included in file
            if library not in library_aliases and parent_sample not in library_aliases:
                continue
            else: 
                # skip libraries if aliquot ID is provided and incorrect
                if library in library_aliases:
                    if library_aliases[library] and library_aliases[library] != aliquot:
                        continue
                elif parent_sample in library_aliases:
                    if library_aliases[parent_sample] and library_aliases[parent_sample] != aliquot:
                        continue
        # check if sample or library is excluded 
        if sample_name in exclude or library in exclude or parent_sample in exclude:
            if run_id not in K:
                K[run_id] = []
            K[run_id].append(file_path)
        else:
            # record files per run
            if run_id not in D:
                D[run_id] = []
            D[run_id].append(file_path)
            # record md5sum
            if run_id not in M:
                M[run_id] = []
            M[run_id].append([file_path, i[47]]) 

    # select the files corresponding to the most recent workflow ID for each run
    # get the corresponding md5sums of the file
    for i in D:
        # remove duplicates due to multiple records of full paths because of merging workflows
        D[i] = list(set(D[i]))    
        L = select_most_recent_workflow(D[i])
        # remove duplicate files that are not the most recent
        toremove = [j for j in D[i] if j not in L]
        for j in toremove:
            D[i].remove(j)
        # remove files not selected
        to_remove = [j for j in M[i] if j[0] not in L]
        for j in to_remove:
            M[i].remove(j)
    
    return D, K, M



def generate_links(D, K, project_name, projects_dir, suffix, **keywords):
    '''
    (dict, dict, str, str, str, dict) -> None
    
    Parameters
    ----------
    - D (dict): Dictionary with run, list of files key, value pairs to release
    - K (dict): Dictionary with run, list of files key, value pairs to withold from release
    - project_name (str): Name of the project directory in projects_dir
    - projects_dir (str): Path to the directory in gsi space containing project directories 
    - suffix (str): Indicates fastqs or datafiles
    - keywords (str): Optional run name parameter. Valid option: run_name
    '''
    
    working_dir = os.path.join(projects_dir, project_name)
    os.makedirs(working_dir, exist_ok=True)


    for run in D:
        if 'run_name' in keywords and keywords['run_name']:
            run_name = keywords['run_name']
        else:
            run_name = run + '.{0}.{1}'.format(project_name, suffix)
        run_dir = os.path.join(working_dir, run_name)
        os.makedirs(run_dir, exist_ok=True)
        print('Linking files to folder {0}'.format(run_dir))
        for file in D[run]:
            filename = os.path.basename(file)
            link = os.path.join(run_dir, filename)
            if os.path.isfile(link) == False:
                os.symlink(file, link)

    if len(K) != 0:
        for run in K:
            if 'run_name' in keywords and keywords['run_name']:
                run_name = keywords['run_name'] + '.withhold'
            else:
                run_name = run + '.{0}.{1}.withhold'.format(project_name, suffix)
            run_dir = os.path.join(working_dir, run_name)
            os.makedirs(run_dir, exist_ok=True)
            for file in K[run]:
                filename = os.path.basename(file)
                link = os.path.join(run_dir, filename)
                if os.path.isfile(link) == False:
                    os.symlink(file, link)
    

def link_files(args):
    '''
    (str | None, str | None, list | None, str | None, str, bool, str, str, str, str | list, str) -> None
    
    Parameters
    ----------
    - libraries (str | None): Path to 1 or 2 columns tab-delimited file with library IDs.
                              The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM).
                              The second and optional column is the library aliquot ID (eg. LDI32439).
                              Only the samples with these library aliases are used if provided'
    - files (str | None): Path to file with file names to be released 
    - runs (list | None): List of run IDs. Include one or more run Id separated by white space.
                          Other runs are ignored if provided
    - project (str | None): Project name as it appears in File Provenance Report.
                            Used to parse the FPR by project. Files are further filtered
                            by run is runs parameter if provided, or all files for
                            the project and workflow are used
    - workflow (str): Worflow used to generate the output files
    - nomiseq (bool): Exclude MiSeq runs if activated
    - project_name (str): Project name used to create the project directory in gsi space
    - projects_dir (str): Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/
    - run_name (str): Specifies the run folder name. Run Id or run.withhold as run folder name if not specified. 
    - exclude (str | list): File with sample name or libraries to exclude from the release,
                            or a list of sample name or libraries
    - suffix (str): Indicates map for fastqs or datafiles in the output file name
    - provenance (str): Path to File Provenance Report.
                        Default is '/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz'
    '''
    
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
    
    # get the list of samples/libraries to exclude
    try:
        exclude = list(get_libraries(args.exclude[0]).keys())
    except:
        if args.exclude:
            exclude = list(map(lambda x: x.strip(), args.exclude))
        else:
            exclude = []
     
    # parse library file if exists 
    if args.libraries:
        libraries = get_libraries(args.libraries)
    else:
        libraries = {}

    # get the list of allowed file paths if exists
    if args.files:
        infile = open(args.files)
        files_names = infile.read().rstrip().split('\n')
        files_names = list(map(lambda x: os.path.basename(x), files_names))
        infile.close()
    else:
        files_names = []
        
    # extract files from FPR
    runs = args.runs if args.runs else []
    project = args.project if args.project else ''
    files_release, files_withhold, md5sums = extract_files(provenance, project, runs, args.workflow, args.nomiseq, libraries, files_names, exclude)
    print('Extracted files from File Provenance Report')
    
    # link files to project dir
    if args.suffix == 'fastqs':
        assert args.workflow.lower() in ['bcl2fastq', 'casava', 'fileimport', 'fileimportforanalysis']
    generate_links(files_release, files_withhold, args.project_name, args.projects_dir, args.suffix, run_name = args.run_name)
        
    # write summary md5sums
    run_dir = os.path.join(args.projects_dir, args.project_name)
    os.makedirs(run_dir, exist_ok=True)
    for i in md5sums:
        print('Generating md5sums summary file for run {0}'.format(i))
        filename = i + '.{0}.{1}.md5sums'.format(args.project_name, args.suffix)
        newfile = open(os.path.join(run_dir, filename), 'w')
        for j in md5sums[i]:
            newfile.write('\t'.join([os.path.basename(j[0]), j[1]]) +'\n')
        newfile.close()    
    print('Files were extracted from FPR {0}'.format(provenance))
    
def map_file_ids(L):
    '''
    (list)- > dict
    
    Returns a dictionary mapping files to various IDs including library ID, run ID and external ID
    
    Parameters
    ----------
    - L (list): List of records including library IDs and external IDs
    '''
    
    D = {}
    
    for i in L:
        j = i.rstrip().split('\t')
        try:
            file = j[46]
            run = j[23]
            ID = j[7]
            geo = {k.split('=')[0]:k.split('=')[1] for k in j[12].split(';')}
            lid = j[13]
            barcode = j[27]
            externalid = geo['geo_external_name']
            if 'geo_group_id' in geo:
                groupid = geo['geo_group_id']
            else:
                groupid = '-'
            if 'geo_group_id_description' in geo:
                groupdesc = geo['geo_group_id_description']
            else:
                groupdesc = '-'
            if 'geo_tube_id' in geo:
                tubeid = geo['geo_tube_id']
            else:
                tubeid = '-'
            if 'geo_targeted_resequencing' in geo:
                panel = geo['geo_targeted_resequencing']
            else:
                panel = '-'
            D[file] = [ID, lid, run, barcode, externalid, groupid, groupdesc, tubeid, panel]
        except:
            continue
    return D        


def write_map_file(projects_dir, project_name, run, L, suffix, add_tube, add_panel):
    '''
    (str, str, str, list, str) -> None
    
    Write a table file with external IDs to map released files to user sample IDs
    
    Parameters
    ----------
    - projects_dir (str): Path to the directory in gsi space containing project directories
    - project_name (str): Name of the project directory in projects_dir
    - run (str): Run ID corresponding to the libraries with external Ids in L
    - L (list): List of records including library IDs and external IDs
    - suffix (str): Indicates map for fastqs or datafiles in the output file name
    - add_tube (bool): Add tube_id column to sample map if True
    - add_panel (bool): Add panel column to sample map if True
    '''

    working_dir = os.path.join(projects_dir, project_name)
    os.makedirs(working_dir, exist_ok=True)
    
    output_map = os.path.join(working_dir, '{0}.{1}.{2}.map.txt'.format(run, project_name, suffix))
    newfile = open(output_map, 'w')
    header = ['sample', 'library', 'run', 'barcode', 'external_id', 'group_id', 'group_description']
    if add_tube:
        header.append('tube_id')
    if add_panel:
        header.append('panel')
    newfile.write('\t'.join(header) + '\n')

    for i in L:
        line = i[:-2]
        if add_tube:
            line.append(i[-2])
        if add_panel:
            line.append(i[-1])
        newfile.write('\t'.join(line) + '\n')
    newfile.close()


def list_files_release_folder(directories):
    '''
    (list) -> list
    
    Returns a list of full file paths, links resolved, of the fastqs in each folder of the list of direcrories
    
    Parameters
    ----------
    - directories (list): A list of full paths to directories containining released fastqs 
    '''

    # make a list of files
    files = []
    for run in directories:
        # make a list of links for the directory
        links = [os.path.join(run, i) for i in os.listdir(run) if os.path.isfile(os.path.join(run, i)) and '.fastq.gz' in i]
        filenames = [os.path.realpath(i) for i in links]
        files.extend(filenames)
    return files



def map_filename_workflow_accession(files):
    '''
    (list) -> dict
    
    Returns a dictionary mapping the file name of the released fastqs to the
    workflow accession that generated them
    
    Parameters
    ---------
    - files (list): List of full paths to the released fastqs 
    '''
    
    # map file names to their workflow accession
    file_names = {}
    for file in files:
        filename = os.path.basename(file)
        workflow_accession = get_workflow_id(file)
        assert filename not in file_names
        file_names[filename] = workflow_accession
    return file_names


def is_gzipped(file):
    '''
    (str) -> bool

    Return True if file is gzipped

    Parameters
    ----------
    - file (str): Path to file
    '''
    
    # open file in rb mode
    infile = open(file, 'rb')
    header = infile.readline()
    infile.close()
    if header.startswith(b'\x1f\x8b\x08'):
        return True
    else:
        return False
  

def get_FPR_records(name, provenance, level):
    '''
    (str, str, str) -> list
    
    Returns a list with all the records from the File Provenance Report for a given project or run
    Each individual record in the list is a list of fields    
    
    Parameters
    ----------
    - name (str): Name of a project or run as it appears in File Provenance Report
    - provenance (str): Path to File Provenance Report.
    - level (str): Run or project: Parse FPR for a given run or project
    '''
        
    # get the records for a single project
    records = []
    # open provenance for reading. allow gzipped file or not
    if is_gzipped(provenance):
        infile = gzip.open(provenance, 'rt', errors='ignore')
    else:
        infile = open(provenance)
    for line in infile:
        if name in line:
            line = line.rstrip().split('\t')
            if level == 'project' and name == line[1]:
                records.append(line)
            elif level == 'run' and name == line[18]:
                records.append(line)
    infile.close()
    return records


def collect_info_fastqs(records):
    '''
    (list) -> dict
    
    Returns a dictionary with relevant information for fastqs by parsing the
    File Provenance Report for a given project 
    
    Parameters
    ----------
    - records (list): Project level records from File Provenance Report
    '''
    
    # initiate dict
    D = {}
    
    # loop through each record. each individual record is a list of fields 
    for i in records:
        # only consider workflows generating fastqs
        if 'casava' in i[30].lower() or 'bcl2fastq' in i[30].lower():
            # get file path
            file = i[46]
            # get workflow accession from file path
            workflow_accession = get_workflow_id(file)
            # get file name
            filename = os.path.basename(file)
            # get md5sum, file SWID (unique file identifier) and get the lane of the flow cell
            md5sum, file_swid, lane = i[47], i[44], i[24]
            # get the run ID
            run = i[23]
            # remove lane from run
            run_alias = run[:run.index('_lane')]
            # get read count
            read_count = {k.split('=')[0]:k.split('=')[1] for k in i[45].split(';')}
            read_count = int(read_count['read_count'])
            # get instrument, aliquit ID and barcode and ID
            instrument, lid, barcode, ID = i[22], i[13], i[27], i[7]
            # get information about external sample ID, group ID, description and tube ID 
            geo = {k.split('=')[0]:k.split('=')[1] for k in i[12].split(';')}
            externalid = geo['geo_external_name']
            if 'geo_group_id' in geo:
                groupid = geo['geo_group_id']
            else:
                groupid = 'NA'
            if 'geo_group_id_description' in geo:
                groupdesc = geo['geo_group_id_description']
            else:
                groupdesc = 'NA'
            if 'geo_tube_id' in geo:
                tubeid = geo['geo_tube_id']
            else:
                tubeid = 'NA'
            D[file] = {'filename': filename, 'workflow_id': workflow_accession, 'md5sum': md5sum, 'file_swid': file_swid, 'ID': ID, 'lid': lid,
                       'run': run, 'barcode': barcode, 'external_id': externalid, 'group_id': groupid, 'group_desc': groupdesc,
                       'tube_id': tubeid, 'instrument': instrument, 'read_count': read_count, 'lane': lane, 'run_alias': run_alias}       
    return D



def stich_back_etl_record(record):
    '''
    (list) -> list
    
    Recreate bamqc etl record when fields have internal commas
    Fields with internal commas are also split and need to be stiched back
    
    Parameters
    ----------
    - record (list): List representation of a single record from the comma-separated
    table of bamqc etl, split on commas. 
    '''
    
    # collect indices of fields with double quotes. 
    # these fields have comma-separated values
    L = []
    # count the number of fields with double quotes (and thus with inernal commas)
    total = 0
    for k in range(len(record)):
        if record[k]:
            # find start and end position of fields with double quotes
            if record[k].startswith('"'):
                total += 1
                # make sure only one field has double quotes (and internal commas)
                assert total <= 1
                # collect start positions of all fields with internal commas
                L.append([k])
            if record[k][-1] == '"':
                # if ending with double quotes, collect end positions of all fields with internal quotes
                L[-1].append(k)
    # recreate records by stiching back fields with internal commas
    record = record[:L[0][0]] + [','.join(record[L[0][0]:L[0][1]+1])] + record[L[0][-1]+1:]
    return record
    

def parse_qc_etl(bamqc_table, project):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the bamqc table of qc-etl
        
    Parameters
    ----------
    - bamqc_table (str): Comma-separated table from bamqc generated by qc-etl
    - project (str): Specific project of interest
    '''
            
    # get versions for bamqc4
    #version = subprocess.check_output('/.mounts/labs/gsi/modulator/sw/Ubuntu18.04/gsi-qc-etl-0.53.0/bin/gsiqcetl versions bamqc4', shell=True).decode('utf-8').rstrip()
    
    # get tables
    #tables = subprocess.check_output('/.mounts/labs/gsi/modulator/sw/Ubuntu18.04/gsi-qc-etl-0.53.0/bin/gsiqcetl tables bamqc4 {0}'.format(version), shell=True).decode('utf-8').rstrip()
    
    # get bamqc records
    #bamqc_records =  subprocess.check_output('/.mounts/labs/gsi/modulator/sw/Ubuntu18.04/gsi-qc-etl-0.53.0/bin/gsiqcetl dump -f /scratch2/groups/gsi/production/qcetl/bamqc4/latest bamqc4 5 bamqc4', shell=True).decode('utf-8').rstrip().split('\n')
    
    infile = open(bamqc_table)
    bamqc_records = infile.read().rstrip().split('\n')
    header = bamqc_records.pop(0).split(',')
    infile.close()
    
    # initiate dictionary, collect information for each file at the run level  
    B = {}
    
    # loop over bamqc records, collect relevant information
    for j in bamqc_records:
        if project in j:
            i = j.split(',')
            # some fields have internal commas. need to stitch back entry for that field
            if len(i) != len(header):
                i =  stich_back_etl_record(i)
            # check that fields have been stiched back correctly
            assert len(i) == len(header)
            # collect relevant information
            run_alias = i[header.index('Run Alias')]
            instrument = i[header.index('instrument')]
            bases_mapped = int(i[header.index('bases mapped')])
            coverage = float(i[header.index('coverage')])
            mapped_reads = int(i[header.index('mapped reads')])
            barcodes = i[header.index('Barcodes')]
            lane = i[header.index('Lane Number')]
            sample = i[header.index('sample')]
            percent_duplicate = float(i[header.index('mark duplicates_PERCENT_DUPLICATION')])
            read_length = float(i[header.index('average read length')])
            
            total_bases_on_target = int(i[header.index('total bases on target')])
            total_reads = int(i[header.index('total reads')])
            
            paired_reads = int(i[header.index('paired reads')])
                        
            total_target_size = int(i[header.index('total target size')])
            
            # map file information
            d = {'run_alias': run_alias, 'instrument': instrument,
                 'bases_mapped': bases_mapped, 'coverage': coverage,
                 'mapped_reads': mapped_reads, 'sample': sample,
                 'total_bases_on_target': total_bases_on_target,
                 'total_reads': total_reads, 'barcodes':barcodes, 'lane': lane,
                 'percent_duplicate': percent_duplicate, 'paired_reads': paired_reads,
                 'read_length': read_length, 'total_target_size': total_target_size}
            # collect information at the run level
            if run_alias in B:
                B[run_alias].append(d)
            else:
                B[run_alias] = [d]
    return B



def update_information_released_fastqs(FPR_info, bamqc_info):
    '''
    (dict, dict) -> dict
    
    Update the information obtained from File Provenance Report in place with information
    collected from bamqc and returns a dictionary of libraries, runs for which QC is missing from the bamqc table
    
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - bamqc_info (dict): Information for each file and run for a given project collected from bamqc table
    '''
    
    # collect file info for libraries with missing QC 
    D = {}
    
    # map files to qc-etl data
    for file in FPR_info:
        found_qc = False
        # bamqc info is reported for each run
        run_alias = FPR_info[file]['run_alias']
        # check that run information is present in bamqc
        if run_alias in bamqc_info:
            # map FPR file info with bamqc file info
            for d in bamqc_info[run_alias]:
                # check if same file
                if FPR_info[file]['instrument'].replace('_', ' ') == d['instrument'] \
                    and FPR_info[file]['run_alias'] == d['run_alias'] \
                    and FPR_info[file]['barcode'] == d['barcodes'] \
                    and FPR_info[file]['lane'] == d['lane']:
                        # get sample name. edit to match sample name format in FPR
                        sample_name = d['sample'].split('_')
                        sample_name = '_'.join(sample_name[:2])
                        # check that sample is in file name
                        if sample_name in file:
                            # matched file between FPR and qc-etl. update file information
                            # add coverage
                            FPR_info[file]['bases_mapped'] = d['bases_mapped']
                            FPR_info[file]['mapped_reads'] = d['mapped_reads']
                            FPR_info[file]['total_bases_on_target'] = d['total_bases_on_target']
                            FPR_info[file]['coverage'] = round(d['coverage'], 2)
                            # add percent duplicate
                            FPR_info[file]['percent_duplicate'] = round(d['percent_duplicate'], 2)
                            FPR_info[file]['read_length'] = d['read_length']
                            FPR_info[file]['total_target_size'] = d['total_target_size']
                            found_qc = True
        if found_qc == False:
            print('WARNING. Cannot find bamQC for {0}'.format(os.path.basename(file)))
            if FPR_info[file]['lid'] in D:
                D[FPR_info[file]['lid']].append(FPR_info[file]['run'])
            else:
                D[FPR_info[file]['lid']] = [FPR_info[file]['run']]
            FPR_info[file]['bases_mapped'] = 'NA'
            FPR_info[file]['mapped_reads'] = 'NA'
            FPR_info[file]['total_bases_on_target'] = 'NA'
            FPR_info[file]['coverage'] = 'NA'
            FPR_info[file]['percent_duplicate'] = 'NA'
            FPR_info[file]['read_length'] = 'NA'
            FPR_info[file]['total_target_size'] = 'NA'
    return D

            



def create_ax(row, col, pos, figure, Data1, Data2, YLabel1, YLabel2, color1, color2, title = None, XLabel = None):
    '''
    (int, int, int, matplotlib.figure.Figure, list, list, str, str, str, str, str | None, str | None)
    
    Parameters
    ----------
    
    - row (int): Row position of the plot in figure
    - col (int): Column position of the plot in figure
    - figure (matplotlib.figure.Figure): Matplotlib figure
    - Data1 (list): List of metrics to plot in graph
    - Data2 (list): List of metrics to plot on the same graph as Data1
    - YLabel1 (str): Label of the Y axis for Data1
    - YLabel2 (str): Label of the Y axis for Data2
    - color1 (str): Color of markers and text related to Data1
    - color2 (str): Color of markers and text related to Data2
    - title (str | None): Title of the graph
    - XLabel (str | None): Label of the X axis    
    '''
    
    # create ax in figure to plot data 1
    ax1 = figure.add_subplot(row, col, pos)
    # create ax 2 in figure to plot data 2 using a different y scale
    ax2 = ax1.twinx()

    # plot data 1 and median  
    xcoord = [i/10 for i in range(len(Data1))]
    ax1.plot(xcoord, Data1, clip_on=False, linestyle='', marker= 'o', markerfacecolor = color1, markeredgecolor = color1, markeredgewidth = 1, markersize = 10, alpha=0.5)
    # compute median the data
    median1 = np.median(Data1)
    # plot median and mean. use zorder to bring line to background
    ax1.axhline(y=median1, color=color1, linestyle='-', linewidth=1.5, alpha=0.5, zorder=1)
    
    # plot data 2 and median  
    ax2.plot(xcoord, Data2, clip_on=False, linestyle='', marker= 'o', markerfacecolor = color2, markeredgecolor = color2, markeredgewidth = 1, markersize = 10, alpha=0.5)
    # compute median the data
    median2 = np.median(Data2)
    # plot median and mean. use zorder to bring line to background
    ax2.axhline(y=median2, color=color2, linestyle='-', linewidth=1.5, alpha=0.5, zorder=1)
    
    # start y axis at 0
    for i in [ax1, ax2]:
        i.set_ylim(ymin=0)
       
    # write axis labels
    if XLabel is not None:
        ax1.set_xlabel(XLabel, color='black', size=18, ha='center', weight= 'normal')
    ax1.set_ylabel(YLabel1, color=color1, size=18, ha='center', weight='normal')
    ax2.set_ylabel(YLabel2, color=color2, size=18, ha='center', weight='normal')

    # add title 
    if title is not None:
        ax1.set_title(title, weight='bold', pad =20, fontdict={'fontsize':40})

    # add xticks 
    plt.xticks(xcoord, ['' for i in range(0, len(xcoord), 5)], ha='center', fontsize=12, rotation=0)

    # add splace bewteen axis and tick labels
    for i in [ax1, ax2]:
        i.yaxis.labelpad = 17
        i.xaxis.labelpad = 17
       
    # do not show frame lines  
    for i in [ax1, ax2]:
        i.spines["top"].set_visible(False)    
        i.spines["bottom"].set_visible(True)    
        i.spines["right"].set_visible(False)    
        i.spines["left"].set_visible(False)    
        
    # offset the x axis
    for i in [ax1, ax2]:
        for loc, spine in i.spines.items():
            spine.set_position(('outward', 5))
            spine.set_smart_bounds(True)
    
    # disable scientific notation
    ax1.ticklabel_format(style='plain', axis='y')
    
    return ax1, ax2


def plot_qc_metrics(outputfile, width, height, Data1, Data2, YLabel1, YLabel2, color1, color2, XLabel=None):
    '''
    (str, int, int, list, foat) -> None
    
    Plot 2 metrics of interest on the same graph
        
    Parameters
    ----------
    - outputfile (str): Path to the output figure
    - width (int): Width in inches of the figure
    - height (int): Height in inches of the figure
    - Data1 (list): List of metrics to plot in graph
    - Data2 (list): List of metrics to plot on the same graph as Data1
    - YLabel1 (str): Label of the Y axis for Data1
    - YLabel2 (str): Label of the Y axis for Data2
    - color1 (str): Color of markers and text related to Data1
    - color2 (str): Color of markers and text related to Data2
    - XLabel (str | None): Label of the X axis    
    '''
    
    #create figure
    figure = plt.figure()
    figure.set_size_inches(width, height)
    # plot data
    ax1, ax2 = create_ax(1, 1, 1, figure, Data1, Data2, YLabel1, YLabel2, color1, color2, title = None, XLabel = XLabel)
    # make sure axes do not overlap
    plt.tight_layout(pad = 5)
    # write figure to file  
    figure.savefig(outputfile, bbox_inches = 'tight')
    plt.close()



def generate_figures(project_dir, project_name,  sequencers, library_metrics, metric1, metric2, YLabel1, YLabel2, color1, color2):
    '''
    (project_dir, str, list, dict, str, str, str, str, str, str) -> list
    
    Returns a list of paths to figure files
    
    Parameters
    ----------
    - project_dir (str): Path to the folder where figure files are written
    - project_name (str): name of the project
    - sequencers (list): List of intruments
    - library_metrics (dict): Dictionary with QC metrics of interest from FPR and bamQC for each library 
    - metric1 (str): Metrics of interest 1
    - metric2 (str): Metrics of interest 2
    - YLabel1 (str): Label of the Y axis for Data1
    - YLabel2 (str): Label of the Y axis for Data2
    - color1 (str): Color of markers and text related to Data1
    - color2 (str): Color of markers and text related to Data2
    '''
    
    # generate plots for each instrument. keep track of figure file names
    figure_files = {}
    for i in sequencers:
        outputfile = os.path.join(project_dir, '{0}.{1}.{2}.{3}.QC_plots.png'.format(project_name, i, ''.join(metric1.split()).replace('(%)', ''), ''.join(metric2.split()).replace('(%)', '')))
        # sort read counts in ascending order and coverage according to read count order
        Q1, Q2 = sort_metrics(library_metrics, i, metric1, metric2)
        
        plot_qc_metrics(outputfile, 13, 8, Q1, Q2, YLabel1, YLabel2, color1, color2, 'Samples')
        if i not in figure_files:
            figure_files[i] = []
        figure_files[i].append(outputfile)
    return figure_files


def resize_figure(filename, scaling_factor):
    '''
    (str, float) -> (int, int)
    
    Returns new file size with same proportions as a tuple of height and width
    
    Parameters
    ----------
    - filename (str): Path to figure file
    - scaling_factor (float): The factor applied to resize figure 
    '''
    # extract the original figure size
    image = Image.open(filename)
    width, height = image.size
    # resize while keeping proportions between height and width
    height, width = list(map(lambda x: x * scaling_factor, [height, width]))
    return height, width        


def map_instrument_type(sequencer):
    '''
    (str) -> str

    Returns a generic intrument name for the sequencer model
    
    Parameters
    ----------
    - sequencer (str): Name of the instrument on which libraries are sequenced
    '''
    
    instrument = ''
    if 'miseq' in sequencer.lower():
        instrument = 'MiSeq' 
    elif 'nextseq' in sequencer.lower():
        instrument = 'NextSeq'
    elif 'hiseq' in sequencer.lower():
        instrument = 'HiSeq'
    elif 'novaseq' in sequencer.lower():
        instrument = 'NovaSeq'
    return instrument


def group_qc_metric_by_instrument(library_metrics, metric):
    '''
    (dict, str) - > dict 
    
    Returns metric for each library across runs by instrument
    Uses a generic intrument collapsing sequencer models into a generic name
    
    Parameters
    ----------
    - library_metrics (dict): Dictionary with library-level metrics from FPR and qc-etl 
    - metric (str): Metric of interest
    '''
    
    # collect read_counts across runs for each instrument and file
    D = {}
    for library in library_metrics:
        instrument = map_instrument_type(library_metrics[library]['instrument'])
        if instrument in D:
            D[instrument].append([library_metrics[library][metric], library])
        else:
            D[instrument] = [[library_metrics[library][metric], library]]
    return D


def sort_metrics(library_metrics, instrument, metric1, metric2):
    '''
    (dict, str, str, str) -> (list, list)
    
    Returns a tuple with lists of metrics 1 and 2 sorted according to metrics 1
    
    Parameters
    ----------
    
    - library_metrics (dict): QC metrics extracted from the File provenance Report and bamqc for each library
    - instrument (str): Sequencing intrument
    - metric1 (str): Name of QC metric 1
    - metric2 (str): Name of QC metric 2
    '''
    
    # group metric 1 by instrument
    D1 = group_qc_metric_by_instrument(library_metrics, metric1)
    # make a sorted list of metric1 in ascending order
    M = [i for i in D1[instrument] if i[0] != 'NA']
    M.sort(key = lambda x: x[0])
    Q1 = [i[0] for i in M]
    # make a list of libraries corresponding to the sorted metric1 values
    libraries = [i[1] for i in M]
    
    # group metric 2 by instrument
    D2 = group_qc_metric_by_instrument(library_metrics, metric2)
    # make a list of metric2 sorted according to the order of metric1
    Q2 = []
    for i in libraries:
        for j in D2[instrument]:
            if j[1] == i:
                Q2.append(j[0])
    # remove missing QC values for each list
    while 'NA' in Q1 or 'NA' in Q2:
        if 'NA' in Q1:
            del Q1[Q1.index('NA')]
            del Q2[Q1.index('NA')]
        if 'NA' in Q2:
            del Q1[Q2.index('NA')]
            del Q2[Q2.index('NA')]
    assert len(Q1) == len(Q2)

    return Q1, Q2


def count_released_fastqs_by_instrument(FPR_info):
    '''
    (dict) -> dict
    
    Returns the count of released fastqs for each run and instrument
    
    Parameters
    ----------
    - FPR_info (dict): Information about the released fastqs collected from File Provenance Report
    '''
        
    # count released fastqs by instrument and run
    D = {}
    for file in FPR_info:
        instrument = map_instrument_type(FPR_info[file]['instrument'])
        run = FPR_info[file]['run_alias']
        if instrument not in D:
            D[instrument] = {}
        if run not in D[instrument]:
            D[instrument][run] = 1
        else:
            D[instrument][run] += 1
    return D
    


def encode_image(filename):
    '''
    (str)- > str
    
    Returns a string representing the encoding of image filename in base64
    
    Parameters
    ----------
    - filename (str): Path to the image file
    '''
    
    # encode image in base64
    with open(filename, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return encoded_string



def generate_table(library_metrics, header, column_size):
    '''
    (dict, list, dict, str) -> str
    
    Returns a html string representing a table
    
    Parameters
    ----------
    - library_metrics (dict): Dictionary with library-level metrics
    - header (list):
    - column_size (dict):
    '''
    
    # count the expected number of cells (excluding header) in tables
    cells = len(list(library_metrics.keys()))
    
    # add padding around text in cells    
    padding = '3px'
    
    # set up counter to track odd and even lines
    counter = 0

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('<tr>')
    for i in header:
        table.append('<th style="width:{0}; border-top: 1px solid #000000; border-bottom: 1px solid #000000; border-collapse: collapse; padding: {1}; text-align: left">{2}</th>'.format(column_size[i], padding, i))
    table.append('</tr>')
    # add lines in table
    for library in library_metrics:
        if counter % 2 == 0:
            table.append('<tr style="background-color: #eee">')
        else:
            table.append('<tr style="background-color: #fff">')
        for i in header:
            if i == 'run':
                j = str(library_metrics[library]['run_alias'])
                if ';' in j:
                    j = j.replace(';', ';\n')
            elif i == 'barcode':
                j = str(library_metrics[library][i])
                if '-' in j:
                    j = j.replace('-', '-\n')
            else:
                j = str(library_metrics[library][i])
            if counter + 1 == cells:
                table.append('<td style="border-bottom: 1px solid #000000; padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
            else:
                table.append('<td style="padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
        table.append('</tr>')
        # update counter
        counter += 1
    table.append('</table>')
    return ''.join(table)



def generate_table_md5sum(library_metrics, header, column_size):
    '''
    (dict, list, dict, str) -> str
    
    Returns a html string representing a table
    
    Parameters
    ----------
    - library_metrics (dict): Dictionary with library-level metrics
    - header (list): List of table column names
    - column_size (dict): Dictionary of column name, column width key, value pairs
    '''
    
    # count the expected number of cells (excluding header) in tables
    cells = 0
    for library in library_metrics:
        cells += len(library_metrics[library]['files'])
        
    # add padding around text in cells    
    padding = '3px'
    
    # set up counter to track odd and even lines
    counter = 0

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('<tr>')
    for i in header:
        table.append('<th style="width:{0}; border-top: 1px solid #000000; border-bottom: 1px solid #000000; border-collapse: collapse; padding: {1}; text-align: left">{2}</th>'.format(column_size[i], padding, i))
    table.append('</tr>')
    # add lines in table
    # create a list of file and md5sum
    files = []
    for library in library_metrics:
        for i in library_metrics[library]['files']:
            files.append(i)
    
    for file in files:
        if counter % 2 == 0:
            table.append('<tr style="background-color: #eee">')
        else:
            table.append('<tr style="background-color: #fff">')
        for i in range(len(header)):
            j = str(file[i])
            if counter + 1 == cells:
                table.append('<td style="border-bottom: 1px solid #000000; padding: {0}; font-size: 9.3px; text-align: left;">{1}</td>'.format(padding, j))
            else:
                table.append('<td style="padding: {0}; font-size: 9.3px;; text-align: left;">{1}</td>'.format(padding, j))
        table.append('</tr>')
        # update counter
        counter += 1
    table.append('</table>')
    return ''.join(table)


def generate_project_table(project_name, project_code, current_date, name, email):
    '''
    (str, str, str, str, str) -> str

    Returns a html string representing a table with project information

    Parameters
    ----------
    - project_name (str): Name of the project
    - project_code (str): Project code (PROXXXX) available in MISO
    - current_date (str): Date of the release (Y-M-D)
    - name (str): Name of contact personn releasing the data
    - email (str): Email of the contact personn releasing the data     
    '''

    content = [project_name, project_code, current_date, name, email]
    column_width = [15, 15, 10, 30, 30]

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('<tr>')
    for i in ['project', 'code', 'date', 'Genome Sequence Informatics', 'contact']:
        table.append('<th style="padding: 3px; border: 0.75px solid grey; border-collapse:collapse; text-align: center; font-size:10px">{0}</th>'.format(i))
    table.append('</tr>')
    # add lines in table
    table.append('<tr>')
    for i in range(len(content)):
        if i != len(content) -1:
            table.append('<td style="width: {0}%; padding: 3px; border: 0.75px solid grey;  border-collapse:collapse; text-align: center; font-size:10px">{1}</td>'.format(column_width[i], content[i]))
        else:
            table.append('<td style="width: {0}%; padding: 3px; border: 0.75px solid grey;  border-collapse:collapse; text-align: center; font-size:10px">{1}</td>'.format(column_width[i], '<a href = "mailto: {0}">{0}</a>'.format(content[i])))
    table.append('</tr>')
    table.append('</table>')
    return ''.join(table)


def generate_header_table(logo, width, height, level):
    '''
    (str, int, int, str) -> str

    Returns a html string representing a table with logo and report title

    Parameters
    ----------
    - logo (str): Path to the OICR logo
    - width (int): Width of the logo figure
    - height (int): Height of the logo figure
    - level (str): Single release or cumulative project report. Values: single or cumulative 
    '''

    if level == 'single':
        title = 'Data Release Report'
    elif level == 'cumulative':
        title = 'Cumulative Project Report'

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('</tr>')
    table.append('<tr>')
    table.append('<td style="width: 40%; padding: 3px; text-align: left"><img src="{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(logo, 'logo', width, height))
    table.append('<td style="width: 60%; padding: 3px; text-align: left"><p style="text-align: left; color: black; font-size:30px; font-family: Arial, Verdana, sans-serif; font-weight:bold">  {0}</p></td>'.format(title))
    table.append('</tr>')
    table.append('</table>')
    
    return ''.join(table)



def generate_figure_table(file1, file2=None):
   
    '''
    (str, str | None) -> str

    Returns a html string representing a table with figure files

    Parameters
    ----------
    - file1 (str): Path to the figure file with metrics 1 and 2 
    - file2 (str | None): Path to the figure file with metrics 3 and 4 if exists
    '''

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('</tr>')
    table.append('<tr>')
    # resize figure 1
    height1, width1 = resize_figure(file1, 0.3)
    table.append('<td style="width: 50%; padding: 3px; text-align: left"><img src="{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(file1, os.path.basename(file1), width1, height1))
    if file2:
        # resize figure 
        height2, width2 = resize_figure(file2, 0.3)
        table.append('<td style="width: 50%; padding: 3px; text-align: left"><img src="{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(file2, os.path.basename(file2), width2, height2))
    table.append('</tr>')
    table.append('</table>')
    
    return ''.join(table)


def list_file_count(sequencers, fastq_counts, level):
    '''
    (list, dict, str) -> list
    
    Returns a list of htm strings with counts of released fastqs by instrument and run
    
    Parameters
    ----------
    - sequencers (list): List of sequencing instruments
    - fastq_counts (dict): Counts of released fastqs for each run and instrument
    - level (str): Single release or cumulative project report. Values: single or cumulative 
    '''
    
    # count all files
    c = 0
    for i in fastq_counts:
        for j in fastq_counts[i]:
            c += fastq_counts[i][j]
    # store html in list
    L = []
    L.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">1. File Count</p>')
    if level == 'single':
        L.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">This release includes {0} fastqs. File count is broken down by instrument and run as follow.</p>'.format(c))
    elif level == 'cumulative':
        L.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">{0} fastqs have been released. File count is broken down by instrument and run as follow.</p>'.format(c))
    # add file count broken down by instrument and run
    for instrument in sequencers:
        if instrument in fastq_counts:
            L.append('<p style="text-align: left; color: black; font-size: 12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal">{0}</p>'.format(instrument + ':'))
            #Text.append('### {0}'.format(instrument))
            for run in fastq_counts[instrument]:
                L.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size: 12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li>{0}: {1}<li/></ul>'.format(run, fastq_counts[instrument][run]))
    return L            
    


def get_QC_status_from_nabu(api, file_swid):
    '''
    (str, str) -> (str | None, str)
    
    Returns a tuple with the file QC status and release ticket if file is released
        
    Parameters
    ----------
    - api (str): URL of the nabu API
    - file_swid (str): File unique identifier
    '''
    
    try:
        response = requests.get(api + '/fileqc/{0}'.format(file_swid), {'accept': 'application/json',})
    except:
        qcstatus, ticket = None, 'NA'
            
    # check response code
    if response.status_code == 200:
        d = response.json()
        if d['fileqcs']:
            assert len(d['fileqcs']) == 1
            qcstatus = d['fileqcs'][0]['qcstatus']
            if 'comment' in d['fileqcs'][0]:
                ticket = d['fileqcs'][0]['comment']
            else:
                ticket = 'NA'
        else:
            qcstatus, ticket = None, 'NA'
    else:
        qcstatus, ticket = None, 'NA'

    return qcstatus, ticket
        


def list_released_fastqs_project(api, FPR_info):
    '''
    (str, dict) -> list
    
    Returns a list of fastqs for a given project that were previously released by interrogating QC status in Nabu
    Pre-condition: Released files need to be marked in Nabu
        
    Parameters
    ----------
    - api (str): URL of the nabu API
    - file_swid (str): File unique identifier
    '''
    
    L = []
    
    for file in FPR_info:
        # get file QC status
        qcstatus, ticket = get_QC_status_from_nabu(api, FPR_info[file]['file_swid'])
        ticket = os.path.basename(ticket).upper()
        if qcstatus == 'PASS' and ticket.startswith('GDR'):
            L.append(file)
        elif qcstatus is None:
            print('WARNING. Could not retrieve QC status from Nabu for {0}'.format(file))
    return L



def get_cumulative_metrics(FPR_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with the total number of reads for each sequenced library
        
    Parameters
    ----------
    - FPR_info (dict): Dictionary with information for released fastqs pulled from FPR and qc-etl
    '''
    
    # create a dictionary for each library with cumulative number of reads, cumulative coverage, release dates
    D = {}

    for file in FPR_info:
        ID = FPR_info[file]['ID']
        run = FPR_info[file]['run']
        library = FPR_info[file]['library']
        reads =  FPR_info[file]['reads']      
        if library in D:
            assert D[library]['ID'] == ID
            D[library]['reads'] += reads
            D[library]['run'].append(run)
        else:
            D[library] = {'reads': reads, 'ID': ID, 'run': [run]}
    
    return D



def get_read(fastq_file):
    """
    (file) -- > itertools.zip_longest
    :param fastq_file: a fastq file open for reading in plain text mode
    
    Returns an iterator slicing the fastq into 4-line reads.
    Each element of the iterator is a tuple containing read information
    """
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)


def get_read_length(fastq):
    '''
    
    '''

    # open file for reading
    if is_gzipped(fastq):
        infile = gzip.open(fastq, 'rt', errors='ignore')
    else:
        infile = open(fastq)

    # create an iterator with the reads
    reads = get_read(infile)
    
    for i in reads:
        read_length = len(i[1])
        # assume all the reads have same length
        break
    
    infile.close()    
    return read_length
    
 
 







def compute_average_read_length(fastq):
    '''
    (str) -> float
    
    Returns the mean read length in the fastq file
    
    Parameters
    - fastq (str): Path to fastq file 
    '''
    
    # open file for reading
    if is_gzipped(fastq):
        infile = gzip.open(fastq, 'rt', errors='ignore')
    else:
        infile = open(fastq)

    # create an iterator with the reads
    reads = get_read(infile)
    
    # count read length
    d = {}
    
    for i in reads:
        seq = i[1]
        if len(seq) in d:
            d[len(seq)] += 1
        else:
            d[len(seq)] = 1
    
    mean = sum([i*d[i] for i in d.keys()]) / sum(list(d.values()))
    
    infile.close()
    return mean


def transform_metrics(FPR_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with the total number of reads, and other information for each sequenced library
        
    Parameters
    ----------
    - FPR_info (dict): Dictionary with information for released fastqs pulled from FPR and qc-etl
    '''
    
    # create a dictionary for each library with cumulative number of reads, cumulative coverage, release dates
    D = {}

    for file in FPR_info:
        ID = FPR_info[file]['ID']
        run = FPR_info[file]['run']
        run_alias = FPR_info[file]['run_alias']
        library = FPR_info[file]['lid']
        reads =  FPR_info[file]['read_count']      
        md5sum = FPR_info[file]['md5sum']
        barcode = FPR_info[file]['barcode']
        external_id = FPR_info[file]['external_id']
        instrument =  FPR_info[file]['instrument']
        tube_id = FPR_info[file]['tube_id']
                
        bases_mapped = FPR_info[file]['bases_mapped']
        mapped_reads = FPR_info[file]['mapped_reads']
        total_bases_on_target = FPR_info[file]['total_bases_on_target']
        duplicate = FPR_info[file]['percent_duplicate']
        read_length = FPR_info[file]['read_length']
        total_target_size = FPR_info[file]['total_target_size']
                
        if read_length == 'NA':
            # compute read length
            read_length = get_read_length(file)

        if library in D:
            assert D[library]['ID'] == ID
            assert D[library]['library'] == library
            D[library]['reads'].append(reads)
            
            D[library]['bases_mapped'].append(bases_mapped) 
            D[library]['mapped_reads'].append(mapped_reads)
            D[library]['total_bases_on_target'].append(total_bases_on_target)
            D[library]['duplicate (%)'].append(duplicate) 
            D[library]['read_length'].append(read_length)
            
            assert total_target_size == D[library]['total_target_size']
            
            D[library]['run'].append(run)
            D[library]['run_alias'].append(run_alias)
            D[library]['files'].append([os.path.basename(file), md5sum])
            D[library]['barcode'].append(barcode)
            D[library]['external_id'].append(external_id)
            D[library]['instrument'].append(instrument)
            D[library]['tube_id'].append(tube_id)     
                
        else:
            D[library] = {'reads': [reads], 'ID': ID, 'library': library, 'run': [run], 'run_alias': [run_alias],
                          'files': [[os.path.basename(file), md5sum]], 'barcode': [barcode],
                          'external_id': [external_id], 'instrument': [instrument],
                          'tube_id': [tube_id], 'duplicate (%)': [duplicate],
                          'bases_mapped': [bases_mapped], 'mapped_reads': [mapped_reads],
                          'total_bases_on_target': [total_bases_on_target],
                          'read_length': [read_length], 'total_target_size': total_target_size}
    
    for library in D:
        # collpase these fields
        for i in ['run', 'run_alias', 'barcode', 'external_id', 'instrument', 'tube_id']:
            D[library][i] = ';'.join(list(set(D[library][i])))
        # add read counts
        for i in ['reads', 'bases_mapped', 'mapped_reads', 'total_bases_on_target']:
            if 'NA' not in D[library][i]:
                D[library][i] = sum(D[library][i])
            else:
                D[library][i] = 'NA'
        # compute average read length
        D[library]['read_length'] = sum(D[library]['read_length']) / len(D[library]['read_length'])
        # replace with missing qc 
        if 'NA' in D[library]['duplicate (%)']:
            D[library]['duplicate (%)'] = 'NA'
        else:
            D[library]['duplicate (%)'] = list(set(D[library]['duplicate (%)']))[0]
    return D


def compute_coverage(library_metrics):
    '''
    (dict) -> None
    
    Update dictionary with library-level metrics with average coverage
    
    Parameters
    ----------
    - library_metrics (dict): Dictionary with library-level metrics 
    '''
         
    
    for library in library_metrics:
        # get coverage
        if library_metrics[library]['total_target_size'] != 'NA':
            coverage = library_metrics[library]['read_length'] * library_metrics[library]['reads'] / library_metrics[library]['total_target_size']
            coverage = round(coverage, 2)
        else:
            coverage = 'NA'
        library_metrics[library]['coverage'] = coverage
        


def compute_on_target(library_metrics):
    '''
    (dict) -> None
    
    Update dictionary with library-level metrics with on target rate
    
    Parameters
    ----------
    - library_metrics (dict): Dictionary with library-level metrics 
    '''
    
    for library in library_metrics:
        if library_metrics[library]['total_bases_on_target'] != 'NA' and library_metrics[library]['bases_mapped'] != 'NA':
            total_bases_on_target = library_metrics[library]['total_bases_on_target']
            bases_mapped = library_metrics[library]['bases_mapped']
            if bases_mapped == 0:
                on_target = 'NA'
            else:
                on_target = total_bases_on_target / bases_mapped * 100
                if on_target:
                    on_target = 100
                else:
                    on_target = round(on_target, 2)
                    # fix floating point approximations
                    if math.ceil(on_target) == 100:
                        on_target = math.ceil(on_target)
        else:
            on_target = 'NA'
        library_metrics[library]['on_target'] = on_target



def list_sequencers(fastq_counts):
    '''
    (dict) -> list
    
    Returns a list of sequencer instruments used to sequence the released fastqs
    
    Parameters
    ----------
    - fastq_counts (dict): Dictionary with counts of released files per run and instrument
    '''
        
    # make a list of possible sequencers
    sequencers = ['MiSeq', 'NextSeq', 'HiSeq', 'NovaSeq']
    # remove sequencers if not part of release
    to_remove = [i for i in sequencers if i not in fastq_counts]
    for i in to_remove:
        sequencers.remove(i)
    return sequencers







def write_report(args):
    '''
    (str, str, str, str, str, str, str, list)

    Write a PDF report with QC metrics and released fastqs for a given project

    - project (str): Project name as it appears in File Provenance Report
    - working-dir (str): Path to the directory with project directories and links to fastqs 
    - project_name (str): Project name used to create the project directory in gsi space
    - project_code (str): Project code from MISO
    - bamqc_table (str): Path to the bamqc table of qc-etl
    - contact_name (str): Name of the analyst releasing the data
    - contact_email (str): Email of the analyst releasing the data
    - run_directories (list): List of directories with links to fastqs
    - provenance (str): Path to File Provenance Report.
    - level (str): Simgle release or cumulative project level report. Values: single or cumulative 
    '''
    
    # check that runs are specified for single data release report
    if args.level == 'single' and args.run_directories is None:
        raise ValueError('Please provide a list of run folders')
        
    # get the project directory with release run folders
    project_dir = os.path.join(args.working_dir, args.project_name)
    # get the records for the project of interest
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
    records = get_FPR_records(args.project, provenance, 'project')
    print('Information was extracted from FPR {0}'.format(provenance))
    # collect relevant information from File Provenance Report about fastqs for project 
    FPR_info = collect_info_fastqs(records)
    # keep only info about released fastqs
    if args.level == 'single':
        # make a list of full paths to the released fastqs resolving the links in the run directories
        files = list_files_release_folder(args.run_directories)
    elif args.level == 'cumulative':
        files = list_released_fastqs_project(args.api, FPR_info)
    # map file names to their workflow accession
    file_names = map_filename_workflow_accession(files)
    # keep information about the listed fastqs
    to_remove = [i for i in FPR_info if os.path.basename(i) not in file_names or file_names[os.path.basename(i)] != FPR_info[i]['workflow_id']]
    for i in to_remove:
        del FPR_info[i]
    
    # count the number of released fastqs for each run and instrument
    fastq_counts = count_released_fastqs_by_instrument(FPR_info)
    
    # make a list sequencers used to sequence released fastqs
    sequencers = list_sequencers(fastq_counts)
        
    # collect information from bamqc table
    bamqc_info = parse_qc_etl(args.bamqc_table, args.project)
    
    # update FPR info with QC info from bamqc table and list libraries with missing QC
    missing_qc = update_information_released_fastqs(FPR_info, bamqc_info)
    
    # transform metrics per library instead of files
    library_metrics = transform_metrics(FPR_info)
    # add coverage for each library
    compute_coverage(library_metrics)
    # add on_target rate
    compute_on_target(library_metrics)
    
    # generate figure files
    figure_files1 = generate_figures(project_dir, args.project_name,  sequencers, library_metrics, 'reads', 'coverage', 'Read counts', 'Coverage', '#00CD6C', '#AF58BA')
    if args.level == 'single':
        # plot on target and duplicate rate
        figure_files2= generate_figures(project_dir, args.project_name,  sequencers, library_metrics, 'duplicate (%)', 'on_target', 'Percent duplicate', 'On target', '#009ADE', '#FFC61E')
        figure_files = {i: j + figure_files2[i] for i, j in figure_files1.items()}
    elif args.level == 'cumulative':
        figure_files = figure_files1
    
    # get current date (year-month-day)
    current_date = datetime.today().strftime('%Y-%m-%d')
    
    # write md5sums to separate text file
    if args.level == 'single':
        md5_file = os.path.join(project_dir, '{0}_fastqs_release_{1}.md5'.format(args.project_name, current_date))
        newfile = open(md5_file, 'w')
        newfile.write('\t'.join(['filename', 'md5sum']) +'\n')
        for file in sorted(list(FPR_info.keys())):
            newfile.write('\t'.join([FPR_info[file]['filename'], FPR_info[file]['md5sum']]) + '\n')
        newfile.close()       
    
    
    # make a list to store report
    Text = []
    # add title and logo
    logo = 'OICR_Logo_RGB_ENGLISH.png'
    height, width = resize_figure(logo, 0.085)
    Text.append(generate_header_table(logo, width, height, args.level))
    Text.append('<br />' * 3)
    # add information about project and contact personn
    Text.append(generate_project_table(args.project_name, args.project_code, current_date, args.contact_name, args.contact_email))
    Text.append('<br />' * 2)           
    # list the file count            
    Text.extend(list_file_count(sequencers, fastq_counts, args.level))
    Text.append('<br />')           
    # add QC plots
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">2. QC plots</p>')
    Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">QC plots are reported by instrument. Lines are the median of each metric. <span style="font-style: italic">Read counts</span> and <span style="font-style: italic">percent duplicate</span> are plotted by ascending order. <span style="font-style: italic">Mean coverage</span> and <span style="font-style: italic">on target rate</span> are plotted respectively according to the order of <span style="font-style: italic">read counts</span> and <span style="font-style: italic">percent duplicate</span></p>')
    Text.append('<br />')
    for i in sequencers:
        Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size: 12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li>{0}<li/></ul>'.format(i))
        if args.level == 'single':
            Text.append(generate_figure_table(figure_files[i][0], figure_files[i][1]))
        elif args.level == 'cumulative':
            Text.append(generate_figure_table(figure_files[i][0]))
        Text.append('<br />')

    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
    
    # add table with sample Ids
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 1. Sample identifiers</p>')
    header = ['ID', 'library', 'run', 'barcode', 'external_id']       
    column_size = {'ID': '10%', 'library': '25%', 'run': '35%', 'barcode': '10%', 'external_id': '20%'}
    Text.append(generate_table(library_metrics, header, column_size))            
    
    
    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
                
    # add QC metrics table
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 2. QC metrics</p>')
    if args.level == 'single':
        header = ['ID', 'library', 'run', 'reads', 'coverage', 'on_target', 'duplicate (%)']       
        column_size = {'ID': '10%', 'library': '24%', 'run': '29%', 'reads': '9%', 'coverage': '9%', 'on_target': '8%', 'duplicate (%)': '11%'}
    elif args.level == 'cumulative':
        header = ['ID', 'library', 'run', 'reads', 'coverage']
        column_size = {'ID': '15%', 'library': '25%', 'run': '40%', 'reads': '10%', 'coverage': '10%'}
    Text.append(generate_table(library_metrics, header, column_size))            
    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
        
    # add md5sums
    if args.level == 'single':
        Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 3. List of md5sums</p>')
        Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">A list of md5sums is also available in the accompanying file: <span style="color: black; font-style: italic">{0}</span></p>'.format(os.path.basename(md5_file)))
        header = ['filename', 'md5sum']       
        column_size = {'filename': '70%', 'md5sum': '30%'}
        Text.append(generate_table_md5sum(library_metrics, header, column_size))
    # convert to html
    renderer = mistune.Markdown()
    Text = '\n'.join(Text)
    html_str = renderer(Text)
    
    # convert html to pdf    
    report_file = os.path.join(project_dir, '{0}_Fastq_data_release_report.{1}.pdf'.format(args.project_name, current_date))
    newfile = open(report_file, "wb")
    pisa.CreatePDF(html_str, newfile)
    newfile.close()

    # remove figure files from disk
    for i in figure_files:
        for j in figure_files[i]:
            os.remove(j)
    




def map_external_ids(args):
    '''
    (str | None, str | None, list | None, str | None, str, bool, str, str, str | list, str) -> None

    Parameters
    ----------    
    - libraries (str | None): Path to 1 or 2 columns tab-delimited file with library IDs.
                              The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM).
                              The second and optional column is the library aliquot ID (eg. LDI32439).
                              Only the samples with these library aliases are used if provided'
    - files (str | None): Path to file with file names to be released 
    - runs (list | None): List of run IDs. Include one or more run Id separated by white space.
                          Other runs are ignored if provided
    - project (str | None): Project name as it appears in File Provenance Report.
                            Used to parse the FPR by project. Files are further filtered
                            by run is runs parameter if provided, or all files for
                            the project and workflow are used
    - workflow (str): Worflow used to generate the output files
    - nomiseq (bool): Exclude MiSeq runs if activated
    - project_name (str): Project name used to create the project directory in gsi space
    - projects_dir (str): Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/
    - exclude (str | list): File with sample name or libraries to exclude from the release,
                            or a list of sample name or libraries
    - suffix (str): Indicates map for fastqs or datafiles in the output file name
    - add_tube (bool): Add tube_id column to sample map if True
    - add_panel (bool): Add panel column to sample map if True
    - provenance (str): Path to File Provenance Report.
                        Default is '/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz'
    '''

    try:
        exclude = list(get_libraries(args.exclude[0]).keys())
    except:
        if args.exclude:
            exclude = list(map(lambda x: x.strip(), args.exclude))
        else:
            exclude = []
     
    # parse library file if exists 
    if args.libraries:
        libraries = get_libraries(args.libraries)
    else:
        libraries = {}


    # get the list of allowed file paths if exists
    if args.files:
        infile = open(args.files)
        files_names = infile.read().rstrip().split('\n')
        files_names = list(map(lambda x: os.path.basename(x), files_names))
        infile.close()
    else:
        files_names = []

    # extract files from FPR
    runs = args.runs if args.runs else []
    project = args.project if args.project else ''
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
    files_release, files_withhold, md5sums = extract_files(provenance, project, runs, args.workflow, args.nomiseq, libraries, files_names, exclude)
    print('Extracted files from FPR')
    
    for run in files_release:
        # make a list of items to write to file
        # make a list of records items
        print('Generating map for run {0}'.format(run))
        L, recorded = [], []
        # grab all the FPR records for that run
        records = get_FPR_records(run, provenance, 'run')
        # map files in run to external IDs
        mapped_files = map_file_ids(records)
        for file in files_release[run]:
            if file in mapped_files:
                R = mapped_files[file]
                if R not in recorded:
                    L.append(R)
                recorded.append(R)
        
        write_map_file(args.projects_dir, args.project_name, run, L, args.suffix, args.add_tube, args.add_panel)
    print('Information was extracted from FPR {0}'.format(provenance))

def map_swid_file(L):
    '''
    (list, str) -> dict
    
    Returns a dictionary mapping file paths to their swid 
    
    Parameters
    ----------
    - L (list): List of string records from the File Provenance Report
    '''

    # create a dict {file path: swid}
    D = {} 

    for i in L:
        j = i.rstrip().split('\t')
        swid = j[44]
        # get file path and resolve any links
        file = os.path.realpath(j[46])
        D[file] = swid         
    return D



def change_nabu_status(api, file_swid, qc_status, user_name, comment=None):
    '''
    (str, str, str, str, str | None) -> None
    
    Modify the record of file with file_swid in Nabu with QC status updated to qc_status,
    username updated to user_name and comment updated to comment if used 
    
    Parameters
    ----------
    - api (str): URL of the nabu API
    - file_swid (str): File unique identifier
    - qc_status (str): File QC status: PASS, PENDING or FAIL
    - comment (str): Jira ticket of the release
    '''
    
    if qc_status not in ['PASS', 'PENDING', 'FAIL']:
        raise ValueError('QC status is PASS, FAIL or PENDING')
    
    if comment:
        response = requests.post(api + '/fileqcs?fileswid={0}&qcstatus={1}&username={2}&comment={3}'.format(file_swid, qc_status, user_name, comment))
    else:
        response = requests.post(api + '/fileqcs?fileswid={0}&qcstatus={1}&username={2}'.format(file_swid, qc_status, user_name))
    
    # check response code
    if response.status_code == 201:
        # record created
        print('Sucessfully updated {0} status to {1}'.format(file_swid, qc_status))
    else:
        print('Could not update {0} status. Nabu response code: {1}'.format(file_swid, response.status_code))


def mark_files_nabu(args):
    '''
    (str, str, str) -> None

    Mark released files with user name and PASS and withheld files with user name and FAIL in Nabu

    Parameters
    ----------    
    - directory (str): Directory with links organized by project and run in gsi space 
    - status (str): Mark files accordingly when released or withheld. Valid options:
                     - pass: files that are released
                     - fail: files that are withheld
    - user (str): User name to appear in Nabu for each released or whitheld file
    - comment (str): A comment to used to tag the file. For instance the Jira ticket 
    - provenance (str): Path to File Provenance Report
    '''
    
    # check directory
    if os.path.isdir(args.directory) == False:
        raise ValueError('{0} is not a valid directory'.format(args.directory))
        
    # make a list of files in directory
    # get directory name and run Id
    run = os.path.basename(args.directory)
    if '.' in run:
        run_id = run[:run.index('.')]
    else:
        run_id = run
    if '.withhold' in run:
        if args.status.lower() != 'fail':
            raise ValueError('this run should not be released. Expected release value is fail')
    else:
        if args.status.lower() == 'fail':
            raise ValueError('this run should be released. Expected release value is pass')
    
    # get the real path of the links in directory
    files = [os.path.realpath(os.path.join(args.directory, i)) for i in os.listdir(args.directory) if os.path.isfile(os.path.join(args.directory, i))]
    
    # make a list of swids
    swids = []
    # dereference FPR
    provenance = os.path.realpath(args.provenance)
    try:
        records = get_FPR_records(run_id, provenance, 'run')
    except:
        raise ValueError('Cannot find records for run {0} in FPR {1}'.format(run_id, provenance))
    else:
        mapped_files = map_swid_file(records)
        # get the swid Ids
        swids = [mapped_files[file] for file in files if file in mapped_files]
        # list the released files without any swid
        no_swids = [file for file in files if file not in mapped_files]
        if no_swids:
            for file in no_swids:
                print('File {0} in directory {1} does not have a swid'.format(os.path.basename(file), args.directory))
        # mark files il nabu
        for i in swids:
            change_nabu_status(args.api, i, args.status.upper(), args.user, comment=args.comment)
    print('Information was extracted from FPR {0}'.format(provenance))
            
    
if __name__ == '__main__':

    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'dare.py', description='A tool to manage data release')
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')
    
   	# link files in gsi space 
    l_parser = subparsers.add_parser('link', help="Link files extracted from FPR to gsi space")
    l_parser.add_argument('-l', '--libraries', dest='libraries', help='Path to 1 or 2 columns tab-delimited file with library IDs.\
                          The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM). The second and optional column is the library aliquot ID (eg. LDI32439).\
                          Only the samples with these library aliases are used if provided')
    l_parser.add_argument('-w', '--workflow', dest='workflow', help='Worflow used to generate the output files', required = True)
    l_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space', required=True)
    l_parser.add_argument('-p', '--parent', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    l_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report. Used to parse the FPR by project. Files are further filtered by run is runs parameter if provided, or all files for the project and workflow are used')
    l_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run IDs. Include one or more run Id separated by white space. Other runs are ignored if provided')
    l_parser.add_argument('--exclude_miseq', dest='nomiseq', action='store_true', help='Exclude MiSeq runs if activated')
    l_parser.add_argument('-rn', '--run_name', dest='run_name', help='Optional run name parameter. Replaces run ID and run.withhold folder names if used')
    l_parser.add_argument('-e', '--exclude', dest='exclude', nargs='*', help='File with sample name or libraries to exclude from the release, or a list of sample name or libraries')
    l_parser.add_argument('-s', '--suffix', dest='suffix', help='Indicates if fastqs or datafiles are released by adding suffix to the directory name. Use fastqs or workflow name.', required=True)
    l_parser.add_argument('-f', '--files', dest='files', help='File with file names to be released')
    l_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    l_parser.set_defaults(func=link_files)
    
   	# map external IDs 
    m_parser = subparsers.add_parser('map', help="Map files to external IDs")
    m_parser.add_argument('-l', '--libraries', dest='libraries', help='Path to 1 or 2 columns tab-delimited file with library IDs.\
                          The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM). The second and optional column is the library aliquot ID (eg. LDI32439).\
                          Only the samples with these library aliases are used if provided')
    m_parser.add_argument('-w', '--workflow', dest='workflow', help='Worflow used to generate the output files', required = True)
    m_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space', required=True)
    m_parser.add_argument('-p', '--parent', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    m_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report. Used to parse the FPR by project. Files are further filtered by run is runs parameter if provided, or all files for the project and workflow are used')
    m_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run IDs. Include one or more run Id separated by white space. Other runs are ignored if provided')
    m_parser.add_argument('--exclude_miseq', dest='nomiseq', action='store_true', help='Exclude MiSeq runs if activated')
    m_parser.add_argument('-e', '--exclude', dest='exclude', nargs='*', help='File with sample name or libraries to exclude from the release, or a list of sample name or libraries')
    m_parser.add_argument('-s', '--suffix', dest='suffix', help='Indicates if fastqs or datafiles are released by adding suffix to the directory name. Use fastqs or workflow name.', required=True)
    m_parser.add_argument('-f', '--files', dest='files', help='File with file names to be released')
    m_parser.add_argument('--tube', dest='add_tube', action='store_true', help='Add tube_id to sample map if option is used. By default, tube_id is not added.')
    m_parser.add_argument('--panel', dest='add_panel', action='store_true', help='Add panel to sample if option is used. By default, panel is not added.')
    m_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    m_parser.set_defaults(func=map_external_ids)

    # mark files in nabu 
    n_parser = subparsers.add_parser('mark', help="Mark released or withheld files in Nabu")
    n_parser.add_argument('-u', '--user', dest='user', help='User name to appear in Nabu for each released or whitheld file', required=True)
    n_parser.add_argument('-s', '--status', dest='status', choices = ['fail', 'pass'], help='Mark files accordingly when released or withheld', required = True)
    n_parser.add_argument('-d', '--directory', dest='directory', help='Directory with links organized by project and run in gsi space', required=True)
    n_parser.add_argument('-c', '--comment', dest='comment', help='Comment to be added to the released file')
    n_parser.add_argument('-a', '--api', dest='api', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
    n_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    n_parser.set_defaults(func=mark_files_nabu)
    
    # write a report
    r_parser = subparsers.add_parser('report', help="Write a PDF report for released FASTQs")
    r_parser.add_argument('-p', '--project', dest='project', help='Project name as it appears in File Provenance Report', required=True)
    r_parser.add_argument('-d', '--directory', dest='working_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Project name as it appears in File Provenance Report')
    r_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space', required = True)
    r_parser.add_argument('-c', '--code', dest='project_code', help='Project code from MISO', required = True)
    r_parser.add_argument('-r', '--runs', dest='run_directories', nargs='*', help='List of directories with released fastqs')
    r_parser.add_argument('-q', '--qc_table', dest='bamqc_table', help='Path to the bamqc table', required=True)
    r_parser.add_argument('-ct', '--contact', dest='contact_name', help='Name of the contact personn releasing the data', required=True)
    r_parser.add_argument('-e', '--email', dest='contact_email', help='Email of the contact personn releasing the data', required=True)
    r_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz')
    r_parser.add_argument('-a', '--api', dest='api', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
    r_parser.add_argument('-l', '--level', dest='level', choices=['single', 'cumulative'], help='Generates a single release report or a cumulative project report', required = True)
    r_parser.set_defaults(func=write_report)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)