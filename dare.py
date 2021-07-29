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

# can mark files other than fastqs in nabu?


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
    assert workflow != -1
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
    files = [d[i][-1][1] for i in d]
    return files


def extract_files(project, runs, workflow, nomiseq, library_aliases, files_release, exclude):
    '''
    (str, list, str, bool, dict, list, list) -> (dict, dict, dict)
  
    Returns a tuple with dictionaries with files extracted from FPR and their corresponding run
    respectively from release and withheld from release, and a dictionary with md5sums of the release files.
    Returns the files corresponding to the most recent workflow iteration if duplicate files
            
    Parameters
    ----------
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
        records.extend(subprocess.check_output('zcat /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz | grep {0}'.format(project), shell=True).decode('utf-8').rstrip().split('\n'))
    elif runs:
        for run in runs:
            records.extend(subprocess.check_output('zcat /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz | grep {0}'.format(run), shell=True).decode('utf-8').rstrip().split('\n'))
        
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
        for file in D[run]:
            filename = os.path.basename(file)
            link = os.path.join(run_dir, filename)
            subprocess.call('ln -s {0} {1}'.format(file, link), shell=True)


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
            subprocess.call('ln -s {0} {1}'.format(file, link), shell=True)


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
    files_release, files_withhold, md5sums = extract_files(project, runs, args.workflow, args.nomiseq, libraries, files_names, exclude)
    
    # link files to project dir
    if args.suffix == 'fastqs':
        assert args.workflow == 'bcl2fastq' or args.workflow.lower() == 'casava'
    generate_links(files_release, files_withhold, args.project_name, args.projects_dir, args.suffix, run_name = args.run_name)
        
    # write summary md5sums
    run_dir = os.path.join(args.projects_dir, args.project_name)
    os.makedirs(run_dir, exist_ok=True)
    for i in md5sums:
        filename = i + '.{0}.{1}.md5sums'.format(args.project_name, args.suffix)
        newfile = open(os.path.join(run_dir, filename), 'w')
        for j in md5sums[i]:
            newfile.write('\t'.join([os.path.basename(j[0]), j[1]]) +'\n')
        newfile.close()    
    
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
            D[file] = [ID, lid, run, barcode, externalid, groupid, groupdesc, tubeid]
        except:
            continue
    return D        


def write_map_file(projects_dir, project_name, run, L, suffix):
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
    '''

    working_dir = os.path.join(projects_dir, project_name)
    os.makedirs(working_dir, exist_ok=True)
    
    output_map = os.path.join(working_dir, '{0}.{1}.{2}.map.txt'.format(run, project_name, suffix))
    newfile = open(output_map, 'w')    
    newfile.write("sample\tlibrary\trun\tbarcode\texternal_id\tgroup_id\tgroup_description\ttube_id\n")
    for i in L:
       newfile.write('\t'.join(i) + '\n')
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
    

def get_project_records(project):
    '''
    (str) -> list
    
    Returns a list with all the records from the File Provenance Report for a given project
    Each individual record in the list is a list of fields    
    
    Parameters
    ----------
    - project (str): Name of a project as it appears in File Provenance Report
    '''
        
    # get the records for a single project
    records = subprocess.check_output('zcat /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz | grep {0}'.format(project), shell=True).decode('utf-8').rstrip().split('\n')
    for i in range(len(records)):
        records[i] = records[i].rstrip().split('\t')
    return records


def collect_info_released_fastqs(records, file_names):
    '''
    (list, dict) -> dict
    
    Returns a dictionary with relevant information for each released fastqs
    by parsing the File Provenance Report for a given project 
    
    Parameters
    ----------
    - records (list): Project level records from File Provenance Report
    - file_names (dict): Map of the names of the released fastqs and the workflow accession that generated them
    '''
    
    # initiate dict
    D = {}
    
    # loop through each record. each individual record is a list of fields 
    for i in records:
        # get file path
        file = i[46]
        # get workflow accession from file path
        workflow_accession = get_workflow_id(file)
        # get file name
        filename = os.path.basename(file)
        # only keep information about released fastqs
        if filename in file_names and workflow_accession == file_names[filename]:
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
            D[file] = {'filename': filename, 'md5sum': md5sum, 'file_swid': file_swid, 'ID': ID, 'lid': lid,
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
            
            total_bases_on_target = int(i[header.index('total bases on target')])
            total_reads = int(i[header.index('total reads')])
            if bases_mapped == 0:
                on_target = 'NA'
            else:
                on_target = total_bases_on_target / bases_mapped * 100
            # map file information
            d = {'run_alias': run_alias, 'instrument': instrument,
                 'bases_mapped': bases_mapped, 'coverage': coverage,
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



def update_information_released_fastqs(FPR_info, bamqc_info):
    '''
    (dict, dict) -> None
    
    Update the information obtained from File Provenance Report in place with information
    collected from bamqc
    
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - bamqc_info (dict): Information for each file and run for a given project collected from bamqc table
    '''
    
    # map files to qc-etl data
    for file in FPR_info:
        # bamqc info is reported for each run
        run_alias = FPR_info[file]['run_alias']
        # check that run information is present in bamqc
        if run_alias in bamqc_info:
            # map FPR file info with bamqc file info
            for d in bamqc_info[run_alias]:
                # check if same file
                if FPR_info[file]['instrument'].replace('_', ' ') == d['instrument'] and FPR_info[file]['run_alias'] == d['run_alias']:
                    if FPR_info[file]['barcode'] == d['barcodes'] and FPR_info[file]['lane'] == d['lane']:
                        # get sample name. edit to match sample name format in FPR
                        sample_name = d['sample'].split('_')
                        sample_name = '_'.join(sample_name[:2])
                        # check that sample is in file name
                        if sample_name in file:
                            # matched file between FPR and qc-etl. update file information
                            # add coverage
                            FPR_info[file]['coverage'] = round(d['coverage'] * 100, 2)
                            # add percent duplicate
                            FPR_info[file]['percent_duplicate'] = round(d['percent_duplicate'], 2)
                            # add on_target
                            if d['on_target'] > 100:
                                FPR_info[file]['on_target'] = 100
                            else:
                                FPR_info[file]['on_target'] = round(d['on_target'], 2)
                            # fix floating point approximations
                            if math.ceil(FPR_info[file]['on_target']) == 100:
                                FPR_info[file]['on_target'] = math.ceil(FPR_info[file]['on_target'])


def rename_metrics_FPR(FPR_info):
    '''
    
    
    '''
    
    for file in FPR_info:
        FPR_info[file]['duplicate (%)'] =  FPR_info[file]['percent_duplicate']
        FPR_info[file]['library'] =  FPR_info[file]['lid']
        FPR_info[file]['reads'] =  FPR_info[file]['read_count']
        
        for  i in ['percent_duplicate', 'lid', 'read_count']:
            del FPR_info[file][i]
    return FPR_info                        
    




                    

#def create_ax(row, col, pos, figure, Data, YLabel, title = None, XLabel = None):
#    '''
#    
#    
#    
#    
#    
#    
#    '''
#    
#    # create ax in figure to plot data 1
#    ax = figure.add_subplot(row, col, pos)
#    
#    # plot data 
#    xcoord = [i/10 for i in range(len(Data))]
#    ax.plot(xcoord, Data, clip_on=False, linestyle='', marker= 'o', markerfacecolor = 'pink', markeredgecolor = 'grey', markeredgewidth = 1, markersize = 10)
#
#    # compute median and mean of the data
#    median, mean = np.median(Data), np.mean(Data)
#    
#    # plot median and mean. use zorder to bring line to background
#    ax.axhline(y=mean, color='r', linestyle='-', linewidth=1.5, alpha=0.5, zorder=1)
#    ax.axhline(y=median, color='blue', linestyle='-', linewidth=1.5, alpha=0.5, zorder=1)
#    
#    # write axis labels
#    if XLabel is not None:
#        ax.set_xlabel(XLabel, color='black', size=18, ha='center', weight= 'normal')
#    ax.set_ylabel(YLabel, color='black', size=18, ha='center', weight='normal')
#
#    # add title 
#    if title is not None:
#        ax.set_title(title, weight='bold', pad =20, fontdict={'fontsize':40})
#
#    # add xticks 
#    plt.xticks(xcoord, ['' for i in range(0, len(xcoord), 5)], ha='center', fontsize=12, rotation=0)
#
#    # add splace bewteen axis and tick labels
#    ax.yaxis.labelpad = 17
#    ax.xaxis.labelpad = 17
#    
#    # do not show frame lines  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)    
#    ax.spines["left"].set_visible(False)    
#        
#    # offset the x axis
#    for loc, spine in ax.spines.items():
#        spine.set_position(('outward', 5))
#        spine.set_smart_bounds(True)
#    
#    # keep ticks only along the x axis, edit font size and change tick direction
#    if XLabel is not None:
#        plt.tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=False,
#                    labelbottom=True, colors = 'black', labelsize = 12, direction = 'out')
#        plt.tick_params(axis='x', which='both', bottom=True, top=False, right=False, left=False,
#                    labelbottom=True, colors = 'black', labelsize = 12, direction = 'out', labelrotation=90)
#        plt.tick_params(axis='y', which='both', bottom=True, top=False, right=False, left=False,
#                    labelbottom=True, colors = 'black', labelsize = 14, direction = 'out')
#    else:
#        plt.tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=False,
#                    labelbottom=False, colors = 'black', labelsize = 14, direction = 'out')
#    
#    # disable scientific notation
#    ax.ticklabel_format(style='plain', axis='y')
#    
#    # add a light grey horizontal grid to the plot, semi-transparent, 
#    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, linewidth = 0.5)  
#    # hide these grids behind plot objects
#    ax.set_axisbelow(True)
#
#    return ax
#



def create_ax(row, col, pos, figure, Data1, Data2, YLabel1, YLabel2, color1, color2, adjust_yaxis, title = None, XLabel = None):
    '''
    
    
    
    
    
    
    '''
    
    # create ax in figure to plot data 1
    ax1 = figure.add_subplot(row, col, pos)
    # create ax 2 in figure to plot data 2 using a different y scale
    ax2 = ax1.twinx()

    # plot data 1 and median  
    xcoord = [i/10 for i in range(len(Data1))]
    ax1.plot(xcoord, Data1, clip_on=False, linestyle='', marker= 'o', markerfacecolor = color1, markeredgecolor = color1, markeredgewidth = 1, markersize = 10)
    # compute median the data
    median1 = np.median(Data1)
    # plot median and mean. use zorder to bring line to background
    ax1.axhline(y=median1, color=color1, linestyle='-', linewidth=1.5, alpha=0.5, zorder=1)
    
    # plot data 2 and median  
    ax2.plot(xcoord, Data2, clip_on=False, linestyle='', marker= 'o', markerfacecolor = color2, markeredgecolor = color2, markeredgewidth = 1, markersize = 10)
    # compute median the data
    median2 = np.median(Data2)
    # plot median and mean. use zorder to bring line to background
    ax2.axhline(y=median2, color=color2, linestyle='-', linewidth=1.5, alpha=0.5, zorder=1)
    
    # start y axis at 0
    for i in [ax1, ax2]:
        i.set_ylim(ymin=0)
    if adjust_yaxis:
        ax1.set_ylim(ymax=max(Data1))
        #ax2.set_ylim(ymax=max(Data2))
        #ax2.set_yticks([0,20,40,60,80,100]) 
        #ax2.set_yticklabels(['0',  '20',  '40',  '60',  '80', '100'])
        
        ax2.yaxis.set_ticks(np.arange(0, max(Data2)+20, 20))
        ax2.yaxis.set_ticklabels(['0',  '20',  '40',  '60',  '80', '100'])
    #ax2.yaxis.set_ticks(np.arange(0, max(Data2)+20, 20))        
    #ax2.yaxis.set_ticklabels(['0',  '20',  '40',  '60',  '80', '100'])
        
        #plt.setp(ax2.get_yticklabels(), visible=True, rotation=30, ha='right')
    
    
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
    
    # keep ticks only along the x axis, edit font size and change tick direction
    if XLabel is not None:
        plt.tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=False,
                    labelbottom=True, colors = 'black', labelsize = 12, direction = 'out')
        plt.tick_params(axis='x', which='both', bottom=True, top=False, right=False, left=False,
                    labelbottom=True, colors = 'black', labelsize = 12, direction = 'out', labelrotation=90)
        plt.tick_params(axis='y', which='both', bottom=True, top=False, right=False, left=False,
                    labelbottom=True, colors = 'black', labelsize = 14, direction = 'out')
    else:
        plt.tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=False,
                    labelbottom=False, colors = 'black', labelsize = 14, direction = 'out')
    
    # disable scientific notation
#    for i in [ax1, ax2]:
#        i.ticklabel_format(style='plain', axis='y')
#        i.ticklabel_format(style='plain', axis='y')
    ax1.ticklabel_format(style='plain', axis='y')
    
    
    
#    # add a light grey horizontal grid to the plot, semi-transparent, 
#    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, linewidth = 0.5)  
#    # hide these grids behind plot objects
#    ax1.set_axisbelow(True)

    return ax1, ax2


















def plot_qc_metrics(outputfile, width, height, read_counts, coverage, YLabel1, YLabel2, color1, color2, adjust_yaxis, XLabel=None):
    '''
    (str, int, int, list, foat) -> None
    
    
    Parameters
    ----------
    - outputfile (str): Path to the output figure
    - width (int): Width in inches of the figure
    - height (int): Height in inches of the figure
    - read_counts (list): Ordered list of read counts for each sample across runs
    - coverage (float): Ordered list of coverage for each sample across runs
   
    
    
    
    -YLabel1:
        
        
    -YLabel2:    
        
        
    -color1


    -  color2    
        
        
    - adjust_yaxis
    
    '''
    
    #create figure
    figure = plt.figure()
    figure.set_size_inches(width, height)
    
    #figure = plt.figure(1, figsize = (4, 2))
    
    
    
    # plot read count
    #ax1 = create_ax(2, 1, 1, figure, read_counts, 'Read counts', title = None, XLabel = None)
    # plot coverage
    #ax2 = create_ax(2, 1, 2, figure, coverage, 'Coverage', title = None, XLabel = None)
    
    
    ax1, ax2 = create_ax(1, 1, 1, figure, read_counts, coverage, YLabel1, YLabel2, color1, color2, adjust_yaxis, title = None, XLabel = XLabel)

    
    
    
    # make sure axes do not overlap
    plt.tight_layout(pad = 5)
      
    figure.savefig(outputfile, bbox_inches = 'tight')
    plt.close()


#def generate_figures(project_dir, project_name,  sequencers, read_counts, coverage, metric1, metric2):
#    '''
#    (str, str, list, dict, dict, str, str)
#    
#    
#    
#    
#    
#    '''
#    
#    
#    
#    
#    
#    # generate plots for each instrument. keep track of figure file names
#    figure_files = {}
#    for i in sequencers:
#        if i in read_counts and i in coverage:
#            outputfile = os.path.join(project_dir, '{0}.{1}.{2}.{3}.QC_plots.png'.format(project_name, i, metric1, metric2))
#            plot_qc_metrics(outputfile, 12, 10, read_counts[i], coverage[i])
#            if i not in figure_files:
#                figure_files[i] = []
#            figure_files[i].append(outputfile)
#    return figure_files


def generate_figures(project_dir, project_name,  sequencers, FPR_info, metric1, metric2, YLabel1, YLabel2, color1, color2, adjust_yaxis):
    '''
    (str, str, list, dict, str, str)
    
    
    
    
    
    '''
    
    
    
    
    
    # generate plots for each instrument. keep track of figure file names
    figure_files = {}
    for i in sequencers:
        outputfile = os.path.join(project_dir, '{0}.{1}.{2}.{3}.QC_plots.png'.format(project_name, i, ''.join(metric1.split()).replace('(%)', ''), ''.join(metric2.split()).replace('(%)', '')))
        # sort read counts in ascending order and coverage according to read count order
        Q1, Q2 = sort_metrics(FPR_info, i, metric1, metric2)
        
        #plot_qc_metrics(outputfile, 15, 10, Q1, Q2)
        
        plot_qc_metrics(outputfile, 13, 8, Q1, Q2, YLabel1, YLabel2, color1, color2, adjust_yaxis, 'Samples')
        if i not in figure_files:
            figure_files[i] = []
        figure_files[i].append(outputfile)
    return figure_files









def resize_figure(filename, scaling_factor):
    '''
    (str, float) -> (int, int)
    
    Parameters
    ----------
    - filename (str): Path to figure file
    - scaling_factor (float): The factor applied to resize figure 
    
    Return new file size with same proportions as a tuple of height and width
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


#def group_qc_metric_by_instrument(FPR_info, metric):
#    '''
#    (dict, str) - > dict 
#    
#    Returns sorted metric across runs by instrument
#    Uses a generic intrument collapsing sequencer models into a generic name
#    
#    Parameters
#    ----------
#    - FPR_info (dict): Dictionary with information from FPR updated with information from QC-ETL
#                       for each released fastq of a given project
#    - metric (str): Metric of interest.
#                    Valid values:
#                    - read_count
#                    - coverage
#    '''
#    
#    # collect read_counts across runs for each instrument and file
#    D = {}
#    for file in FPR_info:
#        instrument = map_instrument_type(FPR_info[file]['instrument'])
#        if instrument in D:
#            D[instrument].append(FPR_info[file][metric])
#        else:
#            D[instrument] = [FPR_info[file][metric]]
#    # sort metric for lowest to highest value
#    for i in D:
#        D[i].sort()
#    return D
  

def group_qc_metric_by_instrument(FPR_info, metric):
    '''
    (dict, str) - > dict 
    
    Returns sorted metric across runs by instrument
    Uses a generic intrument collapsing sequencer models into a generic name
    
    Parameters
    ----------
    - FPR_info (dict): Dictionary with information from FPR updated with information from QC-ETL
                       for each released fastq of a given project
    - metric (str): Metric of interest.
                    Valid values:
                    - reads
                    - coverage
    '''
    
    # collect read_counts across runs for each instrument and file
    D = {}
    for file in FPR_info:
        instrument = map_instrument_type(FPR_info[file]['instrument'])
        if instrument in D:
            D[instrument].append([FPR_info[file][metric], file])
        else:
            D[instrument] = [[FPR_info[file][metric], file]]
    # sort metric for lowest to highest value
    for i in D:
        D[i].sort(key=lambda x: x[0])
    return D


def sort_metrics(FPR_info, instrument, metric1, metric2):
    '''
    
    
    '''
    
    # group metric 1 by instrument
    D1 = group_qc_metric_by_instrument(FPR_info, metric1)
    # make a sorted list of metric1 in ascending order
    Q1 = [i[0] for i in D1[instrument] if 'R1' in i[1]]
    # make a list of files corresponding to the sorted metric1 values
    files = [i[1] for i in D1[instrument] if 'R1' in i[1]]
    
    # group metric 2 by instrument
    D2 = group_qc_metric_by_instrument(FPR_info, metric2)
    # make a list of metric 2sorted according to the order of metric1
    Q2 = []
    for i in files:
        for j in D2[instrument]:
            if j[1] == i:
                Q2.append(j[0])
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
        run = FPR_info[file]['run']
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



def generate_table(FPR_info, header, column_size, table_type=None):
    '''
    (dict, list, dict, str)

    '''
    
    # count the expected number of cells (excluding header) in tables
    if table_type == 'md5sum':
        cells = len(list(FPR_info.keys()))
    else:
        cells = len(list(FPR_info.keys())) / 2
    
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
    for file in FPR_info:
        add_cells = False
        if table_type == 'md5sum':
            add_cells = True
        else:
            if 'R1' in file:
                add_cells = True
        if add_cells:
            if counter % 2 == 0:
                table.append('<tr style="background-color: #eee">')
            else:
                table.append('<tr style="background-color: #fff">')
            for i in header:
                if i == 'run' and 'run_alias' in FPR_info[file]:
                    j = str(FPR_info[file]['run_alias'])
                else:
                    j = str(FPR_info[file][i])
                if counter + 1 == cells:
                    table.append('<td style="border-bottom: 1px solid #000000; padding: {0}; text-align: left;">{1}</td>'.format(padding, j))
                else:
                    table.append('<td style="padding: {0}; text-align: left;">{1}</td>'.format(padding, j))
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


def generate_header_table(logo, width, height):
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

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('</tr>')
    table.append('<tr>')
    table.append('<td style="width: 40%; padding: 3px; text-align: left"><img src="{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(logo, 'logo', width, height))
    table.append('<td style="width: 60%; padding: 3px; text-align: left"><p style="text-align: left; color: black; font-size:30px; font-family: Arial, Verdana, sans-serif; font-weight:bold">  Data Release Report</p></td>')
    table.append('</tr>')
    table.append('</table>')
    
    return ''.join(table)



def generate_figure_table(file1, file2):
   
    '''
    (str, str, int, int) -> str

    Returns a html string representing a table with figure files

    Parameters
    ----------
    - file1 (str)
    
    
    - project_name (str): Name of the project
    
    
    
    - project_code (str): Project code (PROXXXX) available in MISO
    - current_date (str): Date of the release (Y-M-D)
    - name (str): Name of contact personn releasing the data
    - email (str): Email of the contact personn releasing the data     
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
    # resize figure 1
    height2, width2 = resize_figure(file2, 0.3)
    table.append('<td style="width: 50%; padding: 3px; text-align: left"><img src="{0}" alt="{1}" title="{1}" style="padding-right: 0px; padding-left:0px; width:{2}; height:{3}"></td>'.format(file2, os.path.basename(file2), width2, height2))
    table.append('</tr>')
    table.append('</table>')
    
    return ''.join(table)






def list_file_count(sequencers, fastq_counts):
    '''
    (list, dict) -> list
    
    Returns a list of htm strings with counts of released fastqs by instrument and run
    
    Parameters
    ----------
    - sequencers (list): List of sequencing instruments
    - fastq_counts (dict): Counts of released fastqs for each run and instrument
    '''
    
    # count all files
    c = 0
    for i in fastq_counts:
        for j in fastq_counts[i]:
            c += fastq_counts[i][j]
    # store html in list
    L = []
    L.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">1. File Count</p>')
    L.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">This release includes {0} fastqs. File count is broken down by instrument and run as follow.</p>'.format(c))
                 
    # add file count broken down by instrument and run
    for instrument in sequencers:
        if instrument in fastq_counts:
            L.append('<p style="text-align: left; color: black; font-size: 12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal">{0}</p>'.format(instrument + ':'))
            #Text.append('### {0}'.format(instrument))
            for run in fastq_counts[instrument]:
                # remove lane
                run_name = run[:run.index('_lane')]
                L.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size: 12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li>{0}: {1}<li/></ul>'.format(run_name, fastq_counts[instrument][run]))
    return L            
    

def write_report(args):
    '''

    - project: Project name as it appears in File Provenance Report
    - working-dir: default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Project name as it appears in File Provenance Report')
    - project_name: Project name used to create the project directory in gsi space
    - project_code: Project code from MISO
    - bamqc_table: Path to the bamqc table of qc-etl
    '''
    
    # get the project directory with release run folders
    project_dir = os.path.join(args.working_dir, args.project_name)
    
    # make a list of full paths to the released fastqs resolving the links in the run directories
    files = list_files_release_folder(args.run_directories)
    # map file names to their workflow accession
    file_names = map_filename_workflow_accession(files)
    # get the records for the project of interest
    records = get_project_records(args.project)
    # collect relevant information from File Provenance Report about each released fastq 
    FPR_info = collect_info_released_fastqs(records, file_names)
    # collect information from bamqc table
    bamqc_info = parse_qc_etl(args.bamqc_table, args.project)
    # update FPR info with QC info from bamqc table
    update_information_released_fastqs(FPR_info, bamqc_info)
    # rename QC metrics for tables
    FPR_info = rename_metrics_FPR(FPR_info)
    
    # count the number of released fastqs for each run and instrument
    fastq_counts = count_released_fastqs_by_instrument(FPR_info)
    
    # generate plots for coverage, read count by instrument
    # group read_counts across runs for each instrument
    
    #read_counts = group_qc_metric_by_instrument(FPR_info, 'read_count')
    #read_counts = group_qc_metric_by_instrument(FPR_info, 'reads')
    # group coverage across runs for each instrument
    #coverage = group_qc_metric_by_instrument(FPR_info, 'coverage')
    
    # make a list of possible sequencers
    sequencers = ['MiSeq', 'NextSeq', 'HiSeq', 'NovaSeq']
    # remove sequencers if not part of release
    to_remove = [i for i in sequencers if i not in fastq_counts]
    for i in to_remove:
        sequencers.remove(i)
    # generate plots for each instrument. keep track of figure file names
    #figure_files = generate_figures(project_dir, args.project_name,  sequencers, read_counts, coverage)
          

    figure_files1 = generate_figures(project_dir, args.project_name,  sequencers, FPR_info, 'reads', 'coverage', 'Read counts', 'Coverage', '#00CD6C', '#AF58BA', False)
    
    figure_files2= generate_figures(project_dir, args.project_name,  sequencers, FPR_info, 'duplicate (%)', 'on_target', 'Percent duplicate', 'On target', '#009ADE', '#FFC61E', False)
    
    figure_files = {i: j + figure_files2[i] for i, j in figure_files1.items()}
    
    # get current date (year-month-day)
    current_date = datetime.today().strftime('%Y-%m-%d')
    
    # write md5sums to file
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
    Text.append(generate_header_table(logo, width, height))
    #Text.append('<p style="text-align: left; color: black; font-size:30px; font-family: Arial, Verdana, sans-serif; font-weight:bold"> <img src="{0}" alt="{1}" title="{1}" style="float: left; padding-right: 300px; padding-left:0px; width:{2}; height:{3}"> Data Release Report</p>'.format(logo, 'logo', width, height))
    Text.append('<br />' * 3)
    # add information about project and contact personn
    Text.append(generate_project_table(args.project_name, args.project_code, current_date, args.contact_name, args.contact_email))
    Text.append('<br />' * 2)           
    # list the file count            
    Text.extend(list_file_count(sequencers, fastq_counts))
    Text.append('<br />')           
    # add QC plots
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">2. QC plots</p>')
    Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">QC plots are reported by instrument. Lines are the median of each metric. <span style="font-style: italic">Read counts</span> and <span style="font-style: italic">percent duplicate</span> are plotted by ascending order. <span style="font-style: italic">Mean coverage</span> and <span style="font-style: italic">on target rate</span> are plotted respectively according to the order of <span style="font-style: italic">read counts</span> and <span style="font-style: italic">percent duplicate</span></p>')
    Text.append('<br />')
    for i in sequencers:
        Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size: 12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li>{0}<li/></ul>'.format(i))
#        for plot in figure_files[i]:
#            # resize image
#            height, width = resize_figure(plot, 0.3)
#            Text.append('<img src="{0}" alt="{1}" title="{1}" style="display:block; padding-right: 100px; padding-left:30px; width:{2}; height:{3}">'.format(plot, 'qc_plot', width, height))
#            Text.append('<br />')
#    
        #Text.append(generate_figure_table(figure_files[i][0], figure_files[i][1], width, height))
        Text.append(generate_figure_table(figure_files[i][0], figure_files[i][1]))
        
        
        Text.append('<br />')





    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
    
    # add table with sample Ids
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 1. Sample identifiers</p>')
    #header = ['ID', 'lid', 'run', 'barcode', 'external_id']       
    #column_size = {'ID': '10%', 'lid': '30%', 'run': '30%', 'barcode': '10%', 'external_id': '20%'}
    
    
    
    header = ['ID', 'library', 'run', 'barcode', 'external_id']       
    column_size = {'ID': '10%', 'library': '30%', 'run': '30%', 'barcode': '10%', 'external_id': '20%'}
    
    
    
    
    Text.append(generate_table(FPR_info, header, column_size))            
    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
                
    #header = ['filename', 'md5sum', 'file_swid', 'ID', 'lid', 'run', 'barcode',
    #          'external_id', 'group_id', 'group_desc', 'tube_id', 'instrument',
    #          'read_count', 'coverage', 'on_target']       
 
    # add QC metrics table
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 2. QC metrics</p>')
    #header = ['ID', 'lid', 'run', 'read_count', 'coverage', 'on_target']       
    #column_size = {'ID': '10%', 'lid': '30%', 'run': '30%', 'read_count': '10%', 'coverage': '10%', 'on_target': '10%'}
    
    #header = ['ID', 'lid', 'run', 'read_count', 'coverage', 'on_target', 'percent_duplicate']       
    #column_size = {'ID': '10%', 'lid': '25%', 'run': '30%', 'read_count': '10%', 'coverage': '10%', 'on_target': '10%', 'percent_duplicate': '5%'}
    
    
    header = ['ID', 'library', 'run', 'reads', 'coverage', 'on_target', 'duplicate (%)']       
    column_size = {'ID': '10%', 'library': '24%', 'run': '29%', 'reads': '9%', 'coverage': '9%', 'on_target': '8%', 'duplicate (%)': '11%'}
    
    
    
    
    Text.append(generate_table(FPR_info, header, column_size))            
    # add page break between plots and tables
    Text.append('<div style="page-break-after: always;"></div>')
    
    # add md5sums
    Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 3. List of md5sums</p>')
    Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">A list of md5sums is also available in the accompanying file: <span style="color: black; font-style: italic">{0}</span></p>'.format(os.path.basename(md5_file)))
    header = ['filename', 'md5sum']       
    column_size = {'filename': '70%', 'md5sum': '30%'}
    Text.append(generate_table(FPR_info, header, column_size, 'md5sum'))            

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
#    for i in figure_files:
#        for j in figure_files[i]:
#            os.remove(j)

##################################################


#mark duplicates_ESTIMATED_LIBRARY_SIZE  qi
#mark duplicates_LIBRARY qs
#mark duplicates_PERCENT_DUPLICATION     qf
#mark duplicates_READ_PAIR_DUPLICATES    qi
#mark duplicates_READ_PAIR_OPTICAL_DUPLICATES    qi
#mark duplicates_READ_PAIRS_EXAMINED     qi
#mark duplicates_UNMAPPED_READS  qi
#mark duplicates_UNPAIRED_READ_DUPLICATES        qi
#mark duplicates_UNPAIRED_READS_EXAMINED qi

      
    
    
    
    
    
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
    files_release, files_withhold, md5sums = extract_files(project, runs, args.workflow, args.nomiseq, libraries, files_names, exclude)
    
    for run in files_release:
        # make a list of items to write to file
        # make a list of records items
        L, recorded = [], []
        # grab all the FPR records for that run
        records = subprocess.check_output('zcat /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz | grep {0}'.format(run), shell=True).decode('utf-8').rstrip().split('\n')
        # map files in run to external IDs
        mapped_files = map_file_ids(records)
        for file in files_release[run]:
            if file in mapped_files:
                R = mapped_files[file]
                if R not in recorded:
                    L.append(R)
                recorded.append(R)
        
        write_map_file(args.projects_dir, args.project_name, run, L, args.suffix)
    

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


def mark_files_nabu(args):
    '''
    (str, str, str) -> None

    Mark released files with user name and PASS and withheld files with user name and FAIL in Nabu

    Parameters
    ----------    
    - directory (str): Directory with links organized by project and run in gsi space 
    - release (str): Mark files accordingly when released or withheld. Valid options:
                     - pass: files that are released
                     - fail: files that are withheld
    - user (str): User name to appear in Nabu for each released or whitheld file
    - comment (str): A comment to used to tag the file. For instance the Jira ticket 
    '''
    
    # make a list of files in directory
    # get directory name
    run = os.path.basename(args.directory)
    if '.withhold' in run:
        if args.release != 'fail':
            raise ValueError('this run should not be released. Expected release value is fail')
    else:
        if args.release == 'fail':
            raise ValueError('this run should be released. Expected release value is pass')

    # get run ID
    run = run[:run.index('.')]
    
    # get the real path of the links in directory
    files = [os.path.realpath(os.path.join(args.directory, i)) for i in os.listdir(args.directory) if os.path.isfile(os.path.join(args.directory, i))]
    
    # make a list of swids
    swids = []
    records = subprocess.check_output('zcat /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz | grep {0}'.format(run), shell=True).decode('utf-8').rstrip().split('\n')
    
    mapped_files = map_swid_file(records)
    
    # get the swid Ids
    swids = [mapped_files[file] for file in files if file in mapped_files]

    # mark files il nabu
    for i in swids:
        print(i)
        if args.comment:
            subprocess.call('perl /.mounts/labs/gsiprojects/gsi/Nabu_DataSignoff/nabu.pl --post --swid {0} --status {1} --user {2} --comment \'{3}\''.format(i, args.release.upper(), args.user, args.comment), shell=True)
        else:
            subprocess.call('perl /.mounts/labs/gsiprojects/gsi/Nabu_DataSignoff/nabu.pl --post --swid {0} --status {1} --user {2}'.format(i, args.release.upper(), args.user), shell=True)
    
    
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
    m_parser.set_defaults(func=map_external_ids)

    # mark files in nabu 
    n_parser = subparsers.add_parser('mark', help="Mark released or withheld files in Nabu")
    n_parser.add_argument('-u', '--user', dest='user', help='User name to appear in Nabu for each released or whitheld file', required=True)
    n_parser.add_argument('-rl', '--release', dest='release', choices = ['fail', 'pass'], help='Mark files accordingly when released or withheld', required = True)
    n_parser.add_argument('-d', '--directory', dest='directory', help='Directory with links organized by project and run in gsi space', required=True)
    n_parser.add_argument('-c', '--comment', dest='comment', help='Comment to be added to the released file')
    n_parser.set_defaults(func=mark_files_nabu)
    
    # write a report
    r_parser = subparsers.add_parser('report', help="Write a PDF report for released FASTQs")
    r_parser.add_argument('-p', '--project', dest='project', help='Project name as it appears in File Provenance Report', required=True)
    r_parser.add_argument('-d', '--directory', dest='working_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Project name as it appears in File Provenance Report')
    r_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space', required = True)
    r_parser.add_argument('-c', '--code', dest='project_code', help='Project code from MISO', required = True)
    r_parser.add_argument('-r', '--runs', dest='run_directories', nargs='*', help='List of directories with released fastqs', required=True)
    r_parser.add_argument('-q', '--qc_table', dest='bamqc_table', help='Path to the bamqc table', required=True)
    r_parser.add_argument('-ct', '--contact', dest='contact_name', help='Name of the contact personn releasing the data', required=True)
    r_parser.add_argument('-e', '--email', dest='contact_email', help='Email of the contact personn releasing the data', required=True)
    r_parser.set_defaults(func=write_report)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)