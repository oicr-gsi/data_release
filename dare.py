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
import io


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
            if 'geo_library_source_template_type' in geo:
                library_source = geo['geo_library_source_template_type']
            else:
                library_source = 'NA'
            if 'geo_tissue_type' in geo:
                tissue_type = geo['geo_tissue_type']
            else:
                tissue_type = 'NA'
            if 'geo_tissue_origin' in geo:
                tissue_origin = geo['geo_tissue_origin']
            else:
                tissue_origin = 'NA'
            
            sample_name = ID + '_' + tissue_origin + '_' + tissue_type +'_' + library_source + '_' + groupid
            
            D[file] = {'filename': filename, 'workflow_id': workflow_accession, 'md5sum': md5sum, 'file_swid': file_swid, 'ID': ID, 'lid': lid,
                       'run': run, 'barcode': barcode, 'external_id': externalid, 'group_id': groupid, 'group_desc': groupdesc,
                       'tube_id': tubeid, 'instrument': instrument, 'read_count': read_count, 'lane': lane, 'run_alias': run_alias,
                       'library_source': library_source, 'tissue_type': tissue_type, 'tissue_origin': tissue_origin,
                       'sample_name': sample_name}
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
    


def compute_on_target_rate(bases_mapped, total_bases_on_target):
    '''
    (int, int) -> float
    
    Returns the percent on target rate
    
    Parameters
    ----------
    - bases_mapped (int): Number of bases mapping the reference
    - total_bases_on_target (int): Number of bases mapping the target
    '''
    
    try:
        on_target = round(total_bases_on_target / bases_mapped * 100, 2)
    except:
        on_target = 'NA'
    finally:
        if on_target != 'NA':
            if math.ceil(on_target) == 100:
                on_target = math.ceil(on_target)
    return on_target


def parse_bamqc(bamqc_table, project):
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
            coverage_dedup = float(i[header.index('coverage deduplicated')])
            library = i[header.index('library')]
            mapped_reads = int(i[header.index('mapped reads')])
            barcodes = i[header.index('Barcodes')]
            lane = i[header.index('Lane Number')]
            sample = i[header.index('sample')]
            percent_duplicate = float(i[header.index('mark duplicates_PERCENT_DUPLICATION')])
            
            total_bases_on_target = int(i[header.index('total bases on target')])
            total_reads = int(i[header.index('total reads')])
                      
            on_target = compute_on_target_rate(bases_mapped, total_bases_on_target)    
            
            # map file information
            d = {'run_alias': run_alias, 'instrument': instrument,
                 'bases_mapped': bases_mapped, 'coverage': coverage,
                 'coverage_dedup': coverage_dedup, 'library': library,
                 'mapped_reads': mapped_reads, 'sample': sample,
                 'total_bases_on_target': total_bases_on_target,
                 'total_reads': total_reads, 'on_target': on_target,
                 'barcodes':barcodes, 'lane': lane, 'percent_duplicate': percent_duplicate}
            # collect information at the run level
            if run_alias in B:
                B[run_alias].append(d)
            else:
                B[run_alias] = [d]
    return B


def parse_merged_bamqc(merged_bamqc_table, project):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the bamqc4merged table of qc-etl
        
    Parameters
    ----------
    - merged_bamqc_table (str): Comma-separated table from bamqc4merged generated by qc-etl
    - project (str): Specific project of interest
    '''
            
    infile = open(merged_bamqc_table)
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
            library = i[header.index('library')]
            instrument = i[header.index('instrument')]
            coverage = float(i[header.index('coverage')])
            coverage_dedup = float(i[header.index('coverage deduplicated')])
            sample = i[header.index('sample')]
            group_id = i[header.index('Group ID')]
            donor = i[header.index('Donor')]
            library_design = i[header.index('Library Design')]
            tissue_origin = i[header.index('Tissue Origin')]
            tissue_type = i[header.index('Tissue Type')]
            sample_name = donor + '_' + tissue_origin + '_' + tissue_type + '_' + library_design
            assert sample == sample_name + '_' + group_id
            assert project == i[header.index('Project')]
            
            # map file information
            d = {'instrument': instrument, 'coverage': coverage, 'coverage_dedup': coverage_dedup,
                 'sample': sample, 'group_id': group_id, 'sample_name': sample_name,
                 'library_design': library_design, 'library': library}
            
            # collect information for each platform
            if instrument in B:
                assert sample not in B[instrument]
                B[instrument][sample] = d
            else:
                B[instrument] = {sample: d}
                
    return B





def map_bamqc_info_to_fpr(FPR_info, bamqc_info):
    '''
    (dict, dict) -> None
    
    Update the information obtained from File Provenance Report in place with information
    collected from bamqc
    
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - bamqc_info (dict): QC information for each paired fastq from the bamqc table
    '''
    
    for file in FPR_info:
        qc_found = False
        run_alias = FPR_info[file]['run_alias']
        # check that run in recorded in bamqc
        if run_alias in bamqc_info:
            # map file info with bamqc info
            for d in bamqc_info[run_alias]:
                if FPR_info[file]['lane'] == d['lane'] \
                    and FPR_info[file]['sample_name'] == d['sample'] \
                    and FPR_info[file]['barcode'] == d['barcodes'] \
                    and FPR_info[file]['instrument'].replace('_', ' ') == d['instrument']:
                        assert FPR_info[file]['lid'] == d['library']    
                        qc_found = True
                        FPR_info[file]['coverage'] = round(d['coverage'], 2)
                        FPR_info[file]['coverage_dedup'] = round(d['coverage_dedup'], 2)
                        FPR_info[file]['on_target'] = round(d['on_target'], 2)                
                        FPR_info[file]['percent_duplicate'] = round(d['percent_duplicate'], 2)

        if qc_found == False:
            FPR_info[file]['coverage'] = 'NA'
            FPR_info[file]['coverage_dedup'] = 'NA'
            FPR_info[file]['on_target'] = 'NA'                
            FPR_info[file]['percent_duplicate'] = 'NA'


           
def map_merged_bamqc_info_to_fpr(FPR_info, bamqc_info):
    '''
    (dict, dict) -> None
    
    Update the information obtained from File Provenance Report in place with information
    collected from bamqcmerged
    
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - bamqc_info (dict): QC metrics for each sample from the bamqc4merged table
    '''
    
    for file in FPR_info:
        qc_found = False
        instrument = FPR_info[file]['instrument'].replace('_', ' ')
        if instrument in bamqc_info:
            sample_name = FPR_info[file]['sample_name']
            if sample_name in bamqc_info[instrument]:
                # map file info with bamqc info
                assert bamqc_info[instrument][sample_name]['sample_name'] == '_'.join([FPR_info[file]['ID'], FPR_info[file]['tissue_origin'], FPR_info[file]['tissue_type'], FPR_info[file]['library_source']])
                assert FPR_info[file]['lid'] in bamqc_info[instrument][sample_name]['library']
                qc_found = True
                FPR_info[file]['coverage'] = round(bamqc_info[instrument][sample_name]['coverage'], 2)
                FPR_info[file]['coverage_dedup'] = round(bamqc_info[instrument][sample_name]['coverage_dedup'], 2)
                # add run-level relevant metrics
                FPR_info[file]['on_target'] = 'NA'                
                FPR_info[file]['percent_duplicate'] = 'NA'
        if qc_found == False:
            FPR_info[file]['coverage'] = 'NA'
            FPR_info[file]['coverage_dedup'] = 'NA'
            FPR_info[file]['on_target'] = 'NA'                
            FPR_info[file]['percent_duplicate'] = 'NA'


def get_run_level_sample_metrics(FPR_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with run-level QC metrics and read counts for each sample across instrument
        
    Parameters
    ----------
    - FPR_info (dict): Run level bam QC metrics and file stats for each released fastq 
    '''
    
    D = {}
    
    for file in FPR_info:
        # set up boolean to match paired end fastqs
        found = False
        # get file info        
        sample = FPR_info[file]['sample_name']
        lane = FPR_info[file]['lane']
        run = FPR_info[file]['run_alias']
        library = FPR_info[file]['lid']
        instrument = FPR_info[file]['instrument']
        barcode = FPR_info[file]['barcode']
        ext_id = FPR_info[file]['external_id']
        donor = FPR_info[file]['ID']
        library_source = FPR_info[file]['library_source']
        read_count = FPR_info[file]['read_count']
        coverage = FPR_info[file]['coverage']
        coverage_dedup = FPR_info[file]['coverage_dedup']
        on_target = FPR_info[file]['on_target']
        duplicate = FPR_info[file]['percent_duplicate']

        if instrument not in D:
            D[instrument] = {}
        
        if sample not in D[instrument]:
            D[instrument][sample] = [{'sample': sample, 'lane': lane, 'run': run, 'library': library,
                         'instrument': instrument, 'barcode': barcode, 'ext_id': ext_id,
                         'donor': donor, 'library_source': library_source,
                         'reads': read_count, 'coverage': coverage,
                         'coverage_dedup': coverage_dedup, 'on_target': on_target,
                         'duplicate (%)': duplicate, 'files': [file]}]
        else:
            # find paired fastq
            for d in D[instrument][sample]:
                if run == d['run'] and lane == d['lane'] and library == d['library'] \
                     and instrument == d['instrument'] and barcode == d['barcode'] \
                     and ext_id == d['ext_id'] and donor == d['donor'] and library_source == d['library_source']:
                         d['reads'] += read_count
                         assert coverage == d['coverage']
                         assert coverage_dedup == d['coverage_dedup']
                         assert duplicate == d['duplicate (%)']
                         assert on_target == d['on_target']
                         d['files'].append(file)
                         # update boolean. paired end fastq found
                         found = True
            if found == False:
                # record file info
                D[instrument][sample].append({'sample': sample, 'lane': lane, 'run': run, 'library': library,
                             'instrument': instrument, 'barcode': barcode, 'ext_id': ext_id,
                             'donor': donor, 'library_source': library_source,
                             'reads': read_count, 'coverage': coverage,
                             'coverage_dedup': coverage_dedup, 'on_target': on_target,
                             'duplicate (%)': duplicate, 'files': [file]})
                    
                    
            
            # if run == D[instrument][sample]['run'] and lane == D[instrument][sample]['lane'] and library == D[instrument][sample]['library'] \
            #     and instrument == D[instrument][sample]['instrument'] and barcode == D[instrument][sample]['barcode'] \
            #     and ext_id == D[instrument][sample]['ext_id'] and donor == D[instrument][sample]['donor'] and library_source == D[instrument][sample]['library_source']:
            #         D[instrument][sample]['read_count'] += read_count
            #         assert coverage == D[instrument][sample]['coverage']
            #         assert coverage_dedup == D[instrument][sample]['coverage_dedup']
            #         assert duplicate == D[instrument][sample]['duplicate']
            #         assert on_target == D[instrument][sample]['on_target']
            #         D[instrument][sample]['files'].append(file)
    return D                         
                                           
            
def get_cumulative_level_sample_metrics(FPR_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with cumulative QC metrics and read counts for each sample across instrument
        
    Parameters
    ----------
    - FPR_info (dict): Cumulative bam QC metrics and file stats for each released fastq 
    '''
    
    D = {}
    
    
    for file in FPR_info:
        sample = FPR_info[file]['sample_name']
        lane = FPR_info[file]['lane']
        run = FPR_info[file]['run_alias']
        library = FPR_info[file]['lid']
        instrument = FPR_info[file]['instrument']
        barcode = FPR_info[file]['barcode']
        ext_id = FPR_info[file]['external_id']
        donor = FPR_info[file]['ID']
        library_source = FPR_info[file]['library_source']
        read_count = FPR_info[file]['read_count']
        coverage = FPR_info[file]['coverage']
        coverage_dedup = FPR_info[file]['coverage_dedup']
        on_target = FPR_info[file]['on_target']
        duplicate = FPR_info[file]['percent_duplicate']
        
        if instrument not in D:
            D[instrument] = {}
               
        if sample not in D[instrument]:
            D[instrument][sample] = {'sample': sample, 'lane': lane, 'run': run, 'library': [library],
                              'instrument': instrument, 'barcode': barcode, 'ext_id': ext_id,
                              'donor': donor, 'library_source': library_source,
                              'read_count': read_count, 'coverage': coverage,
                              'coverage_dedup': coverage_dedup, 'on_target': on_target,
                              'duplicate (%)': duplicate, 'files': [file]}
        else:
            assert ext_id == D[instrument][sample]['ext_id']
            assert donor == D[instrument][sample]['donor']
            assert library_source == D[instrument][sample]['library_source']
            D[instrument][sample]['library'].append(library)  
            D[instrument][sample]['read_count'] += read_count
            assert coverage == D[instrument][sample]['coverage']
            assert coverage_dedup == D[instrument][sample]['coverage_dedup']
            assert duplicate == D[instrument][sample]['duplicate (%)']
            assert on_target == D[instrument][sample]['on_target']
            
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

    # set xticks
    # get all the ticks and set labels to empty str
    plt.xticks(xcoord, ['' for i in range(len(xcoord))], ha='center', fontsize=12, rotation=0)
    # set every N ticks
    N = 3
    xticks_pos = ax1.get_xticks()
    xticks_labels = ax1.get_xticklabels()
    myticks = [j for i,j in enumerate(xticks_pos) if not i % N]  # index of selected ticks
    newlabels = [label for i,label in enumerate(xticks_labels) if not i % N]
    plt.xticks(myticks, newlabels, ha='center', fontsize=12, rotation=0)
    
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
            #spine.set_smart_bounds(True)
    
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



def generate_figures(project_dir, level, project_name,  sample_metrics, metric1, metric2, YLabel1, YLabel2, color1, color2):
    '''
    (project_dir, str, str, dict, str, str, str, str, str, str) -> list
    
    Returns a list of paths to figure files
    
    Parameters
    ----------
    - project_dir (str): Path to the folder where figure files are written
    - level (str): Single or cumulative
    - project_name (str): name of the project
    - sample_metrics (dict): Dictionary with QC metrics of interest from FPR and bamQC for each library 
    - metric1 (str): Metrics of interest 1
    - metric2 (str): Metrics of interest 2
    - YLabel1 (str): Label of the Y axis for Data1
    - YLabel2 (str): Label of the Y axis for Data2
    - color1 (str): Color of markers and text related to Data1
    - color2 (str): Color of markers and text related to Data2
    '''
    
    # make a list of instruments
    instruments = sorted(list(sample_metrics.keys()))
        
    # generate plots for each instrument. keep track of figure file names
    figure_files = {}
    for i in instruments:
        outputfile = os.path.join(project_dir, '{0}.{1}.{2}.{3}.{4}.QC_plots.png'.format(project_name, i, level, ''.join(metric1.split()).replace('(%)', ''), ''.join(metric2.split()).replace('(%)', '')))
        # sort read counts in ascending order and coverage according to read count order
        Q1, Q2 = sort_metrics(sample_metrics, i, metric1, metric2, level)
        plot_qc_metrics(outputfile, 13, 8, Q1, Q2, YLabel1, YLabel2, color1, color2, 'Samples')
        if i not in figure_files:
            figure_files[i] = []
        figure_files[i].append(outputfile)
    return figure_files


def resize_image(image, scaling_factor):
    '''
    (str, float) -> (int, int)
    
    Returns new file size with same proportions as a tuple of height and width
    
    Parameters
    ----------
    - image (str): Path to figure file or base64 encoded image
    - scaling_factor (float): The factor applied to resize figure 
    '''
    
    if os.path.isfile(image):
        # extract the original figure size
        img = Image.open(image)
    else:
        imgdata = base64.b64decode(image)
        img = Image.open(io.BytesIO(imgdata))
    width, height = img.size
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


def group_qc_metric_by_instrument(sample_metrics, metric, level):
    '''
    (dict, str, str) - > dict 
    
    Returns metric for each sample across instruments
    
    
    Parameters
    ----------
    - sample_metrics (dict): Dictionary run-level or cumulative level sample QC metrics  
    - metric (str): Metric of interest
    - level (str): Single or cumulative
    '''
    
    # collect metric for each sample and instrument
    D = {}
        
    # check level because data structure is different for run-level or cumulative metrics
    if level == 'single':
        for instrument in sample_metrics:
            for sample in sample_metrics[instrument]:
                for d in sample_metrics[instrument][sample]:
                    if instrument in D:
                        D[instrument].append([d[metric], sample + '_' + d['run'] + '_' + d['lane']])
                    else:
                        D[instrument] = [[d[metric], sample + '_' + d['run'] + '_' + d['lane']]]
    elif level == 'cumulative':
        for instrument in sample_metrics:
            for sample in sample_metrics[instrument]:
                if instrument in D:
                    D[instrument].append([sample_metrics[instrument][sample][metric], sample])
                else:
                    D[instrument] = [sample_metrics[instrument][sample][metric], sample]
    return D


def sort_metrics(sample_metrics, instrument, metric1, metric2, level):
    '''
    (dict, str, str, str, str) -> (list, list)
    
    Returns a tuple with lists of metrics 1 and 2 sorted according to metrics 1
    
    Parameters
    ----------
    
    - sample_metrics (dict): Run-level or cumulative QC metrics for all samples 
    - instrument (str): Sequencing intrument
    - metric1 (str): Name of QC metric 1
    - metric2 (str): Name of QC metric 2
    - level (str): Single or cumulative
    '''
    
    # group metric 1 by instrument
    D1 = group_qc_metric_by_instrument(sample_metrics, metric1, level)
    
    # make a sorted list of metric1 in ascending order
    M = [i for i in D1[instrument] if i[0] != 'NA']
    M.sort(key = lambda x: x[0])
    Q1 = [i[0] for i in M]
    # make a list of samples corresponding to the sorted metric1 values
    samples = [i[1] for i in M]
    
    # group metric 2 by instrument
    D2 = group_qc_metric_by_instrument(sample_metrics, metric2, level)
    # make a list of metric2 sorted according to the order of metric1
    Q2 = []
    for i in samples:
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
        instrument = FPR_info[file]['instrument']
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



def generate_table(sample_metrics, header, column_size):
    '''
    (dict, list, dict, str) -> str
    
    Returns a html string representing a table
    
    Parameters
    ----------
    - sample_metrics (dict): Dictionary with run-level and cumulative sample metrics
    - header (list):
    - column_size (dict):
    '''
    
    # count the expected number of cells (excluding header) in tables
    cells = 0
    for instrument in sample_metrics:
        for sample in sample_metrics[instrument]:
            cells += len(sample_metrics[instrument][sample])
    
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

    for instrument in sample_metrics:
        samples = sorted(list(sample_metrics[instrument].keys()))
        # add lines in table
        for sample in samples:
            for d in sample_metrics[instrument][sample]:
                if counter % 2 == 0:
                    table.append('<tr style="background-color: #eee">')
                else:
                    table.append('<tr style="background-color: #fff">')
                for i in header:
                    if i == 'external_id':
                        j = str(d['ext_id'])
                    elif i == 'run':
                        j = str(d['run'])
                        if ';' in j:
                            j = j.replace(';', ';\n')
                    elif i == 'barcode':
                        j = str(d[i])
                        if '-' in j:
                            j = j.replace('-', '-\n')
                    else:
                        j = str(d[i])
                    if counter + 1 == cells:
                        table.append('<td style="border-bottom: 1px solid #000000; padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
                    else:
                        table.append('<td style="padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
                table.append('</tr>')
                # update counter
                counter += 1
    table.append('</table>')
    return ''.join(table)






    
    # sort libraries
    # libraries = sorted(list(library_metrics.keys()))
    # # add lines in table
    # for library in libraries:
    #     if counter % 2 == 0:
    #         table.append('<tr style="background-color: #eee">')
    #     else:
    #         table.append('<tr style="background-color: #fff">')
    #     for i in header:
    #         if i == 'run':
    #             j = str(library_metrics[library]['run_alias'])
    #             if ';' in j:
    #                 j = j.replace(';', ';\n')
    #         elif i == 'barcode':
    #             j = str(library_metrics[library][i])
    #             if '-' in j:
    #                 j = j.replace('-', '-\n')
    #         else:
    #             j = str(library_metrics[library][i])
    #         if counter + 1 == cells:
    #             table.append('<td style="border-bottom: 1px solid #000000; padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
    #         else:
    #             table.append('<td style="padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
    #     table.append('</tr>')
    #     # update counter
    #     counter += 1
    # table.append('</table>')
    # return ''.join(table)



def generate_table_md5sum(FPR_info, header, column_size):
    '''
    (dict, list, dict, str) -> str
    
    Returns a html string representing a table
    
    Parameters
    ----------
    - FPR_info (dict): Dictionary with QC metrics and file info
    - header (list): List of table column names
    - column_size (dict): Dictionary of column name, column width key, value pairs
    '''
    
    # count the expected number of cells (excluding header) in tables
    cells = len(FPR_info)
            
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
    for file in FPR_info:
        files.append([FPR_info[file]['filename'], FPR_info[file]['md5sum']])
    files.sort()
    
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


def generate_project_table(project_name, project_full_name, current_date, name, email):
    '''
    (str, str, str, str, str) -> str

    Returns a html string representing a table with project information

    Parameters
    ----------
    - project_name (str): Acronym of the project
    - project_full_name (str): Full name of the project 
    - current_date (str): Date of the release (Y-M-D)
    - name (str): Name of contact personn releasing the data
    - email (str): Email of the contact personn releasing the data     
    '''


    content = [project_name, project_full_name, current_date, name, email]
    column_width = [15, 15, 10, 30, 30]

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('<tr>')
    for i in ['Project', 'Name', 'Date', 'Genome Sequence Informatics', 'Contact']:
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
    - logo (str): String representation of the OICR logo
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
    table.append('<td style="width: 40%; padding: 3px; text-align: left"><img src="data:image/png;base64,{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(logo, 'logo', width, height))
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
    height1, width1 = resize_image(file1, 0.3)
    table.append('<td style="width: 50%; padding: 3px; text-align: left"><img src="{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(file1, os.path.basename(file1), width1, height1))
    if file2:
        # resize figure 
        height2, width2 = resize_image(file2, 0.3)
        table.append('<td style="width: 50%; padding: 3px; text-align: left"><img src="{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(file2, os.path.basename(file2), width2, height2))
    table.append('</tr>')
    table.append('</table>')
    
    return ''.join(table)


def list_file_count(fastq_counts, level):
    '''
    (dict, str) -> list
    
    Returns a list of htm strings with counts of released fastqs by instrument and run
    
    Parameters
    ----------
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
    sequencers = sorted(list(fastq_counts.keys()))
    for instrument in sequencers:
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



def get_logo():
    '''
    (None) -> str
    
    Returns a string representation of the OICR logo
    Pre-condition: The logo was previously converted to base64 image
    '''
    
    logo_image = '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'
    return logo_image


def write_md5sums(project_dir, project_name, current_date, FPR_info):
    '''
    (str, str, str, dict) -> str
    
    Writes a table file with file names and md5sums and returns thefile path
    
    
    
    Parameters
    ----------
    
    - working-dir (str): Path to the directory with project directories and links to fastqs 
    - project_name (str): Project name used to create the project directory in gsi space
    - current_date (str): Current date (Y-m-d)
    - FPR_info (dict): File information with QC metrics
    '''

    md5_file = os.path.join(project_dir, '{0}_fastqs_release_{1}.md5'.format(project_name, current_date))
    newfile = open(md5_file, 'w')
    newfile.write('\t'.join(['filename', 'md5sum']) +'\n')
    for file in sorted(list(FPR_info.keys())):
        newfile.write('\t'.join([FPR_info[file]['filename'], FPR_info[file]['md5sum']]) + '\n')
    newfile.close()       
    return md5_file


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
    
    # collect information from bamqc table
    if args.level == 'single':
        bamqc_info = parse_bamqc(args.bamqc_table, args.project)
        # update FPR info with QC info from bamqc table
        map_bamqc_info_to_fpr(FPR_info, bamqc_info)
        # re-organize metrics per sample and instrument
        sample_metrics = get_run_level_sample_metrics(FPR_info)
    elif args.level == 'cumulative':
        bamqc_info = parse_merged_bamqc(args.bamqc_table, args.project)   
        # update FPR info with QC info from bamqc merged table
        map_merged_bamqc_info_to_fpr(FPR_info, bamqc_info)
        # re-organize metrics per sample and instrument
        sample_metrics = get_cumulative_level_sample_metrics(FPR_info)
    

    # generate figure files
    figure_files1 = generate_figures(project_dir, args.level, args.project_name, sample_metrics, 'reads', 'coverage', 'Read counts', 'Coverage', '#00CD6C', '#AF58BA')
    if args.level == 'single':
        # plot on target and duplicate rate
        figure_files2= generate_figures(project_dir, args.level, args.project_name, sample_metrics, 'duplicate (%)', 'on_target', 'Percent duplicate', 'On target', '#009ADE', '#FFC61E')
        figure_files = {i: j + figure_files2[i] for i, j in figure_files1.items()}
    elif args.level == 'cumulative':
        figure_files = figure_files1
    
    # get current date (year-month-day)
    current_date = datetime.today().strftime('%Y-%m-%d')
    
    # write md5sums to separate text file
    if args.level == 'single':
        md5_file = write_md5sums(project_dir, args.project_name, current_date, FPR_info)
        
    # make a list to store report
    Text = []
    # add title and logo
    logo = str(get_logo())
    height, width = resize_image(logo, 0.085)
    Text.append(generate_header_table(logo, width, height, args.level))
    Text.append('<br />' * 3)
    # add information about project and contact personn
    Text.append(generate_project_table(args.project_name, args.project_full_name, current_date, args.contact_name, args.contact_email))
    Text.append('<br />' * 2)           
    # list the file count            
    Text.extend(list_file_count(fastq_counts, args.level))
    Text.append('<br />')           
    # add QC plots
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">2. QC plots</p>')
    Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">QC plots are reported by instrument. Lines are the median of each metric. <span style="font-style: italic">Read counts</span> and <span style="font-style: italic">percent duplicate</span> are plotted by ascending order. <span style="font-style: italic">Mean coverage</span> and <span style="font-style: italic">on target rate</span> are plotted respectively according to the order of <span style="font-style: italic">read counts</span> and <span style="font-style: italic">percent duplicate</span></p>')
    Text.append('<br />')
    
    
    for i in sorted(list(sample_metrics.keys())):
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
    header = ['donor', 'library', 'run', 'barcode', 'external_id']       
    column_size = {'donor': '10%', 'library': '25%', 'run': '35%', 'barcode': '10%', 'external_id': '20%'}
    Text.append(generate_table(sample_metrics, header, column_size))            
    
    
    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
                
    # add QC metrics table
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 2. QC metrics</p>')
    if args.level == 'single':
        header = ['donor', 'library', 'run', 'reads', 'coverage', 'on_target', 'duplicate (%)']       
        column_size = {'donor': '10%', 'library': '24%', 'run': '29%', 'reads': '9%', 'coverage': '9%', 'on_target': '8%', 'duplicate (%)': '11%'}
    elif args.level == 'cumulative':
        header = ['donor', 'library', 'run', 'reads', 'coverage']
        column_size = {'donor': '15%', 'library': '25%', 'run': '40%', 'reads': '10%', 'coverage': '10%'}
    Text.append(generate_table(sample_metrics, header, column_size))            
    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
        
    # add md5sums
    if args.level == 'single':
        Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 3. List of md5sums</p>')
        Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">A list of md5sums is also available in the accompanying file: <span style="color: black; font-style: italic">{0}</span></p>'.format(os.path.basename(md5_file)))
        header = ['filename', 'md5sum']       
        column_size = {'filename': '70%', 'md5sum': '30%'}
        Text.append(generate_table_md5sum(FPR_info, header, column_size))
    # convert to html
    renderer = mistune.Markdown()
    Text = '\n'.join(Text)
    html_str = renderer(Text)
    
    # convert html to pdf    
    report_name = '{0}_run_level_data_release_report.{1}.pdf' if args.level == 'single' else '{0}_cumulative_data_release_report.{1}.pdf'
    report_file = os.path.join(project_dir, report_name.format(args.project_name, current_date))
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
    r_parser.add_argument('-c', '--code', dest='project_full_name', help='Full name of the project', required = True)
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