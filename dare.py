# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 12:27:07 2020

@author: rjovelin
"""


import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import time
import math
import requests
import gzip
import sys
import json
import pathlib
import sqlite3
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from weasyprint import CSS
import re



def get_libraries(library_file):
    '''
    (str) -> dict
    
    Returns a dictionary with library, run or libray, lane, run key, value pairs.
    Note: runs need to be specified for all libraries
    
    Parameters
    ----------
    - sample_file (str): Path to sample file. Sample file is a tab delimited file
                         that includes 2 or 3 columns columns. The first column is always 
                         the library alias, and the second is lane or run id
    '''
    D = {}
    
    infile = open(library_file)
    content = infile.read().strip().split('\n')
    for i in range(len(content)):
        content[i] = content[i].split('\t')
    infile.close()
    
    # check that all lines have the same number of columns
    if all(map(lambda x: len(x) == 2 or len(x) == 3 , content)) == False:
        raise ValueError('File must have 2 or 3 columns')
    if all(map(lambda x: '_' in x[1], content)) == False:
        raise ValueError('Run id must be the second column')
        
    for i in content:
        library = i[0]
        run = i[1]
        if len(i) == 2:
            lane = ''
        elif len(i) == 3:
            lane = int(i[-1])
        
        if library not in D:
            D[library] = {}
        if run not in D[library]:
            D[library][run] = []
        D[library][run].append(lane)
       
    return D


def get_FPR_records(project, provenance):
    '''
    (str, str) -> list
    
    Returns a list with all the records from the File Provenance Report for a given project.
    Each individual record in the list is a list of fields    
    
    Parameters
    ----------
    - project (str): Name of a project or run as it appears in File Provenance Report
    - provenance (str): Path to File Provenance Report.
    '''
        
    # get the records for a single project
    records = []
    # open provenance for reading. allow gzipped file or not
    if is_gzipped(provenance):
        infile = gzip.open(provenance, 'rt', errors='ignore')
    else:
        infile = open(provenance)
    for line in infile:
        if project in line:
            line = line.rstrip().split('\t')
            if project == line[1]:
                records.append(line)
    infile.close()
    return records



def parse_fpr_records(provenance, project, workflow, prefix=None):
    '''
    (str, str, list, str | None) -> dict
  
    Returns a dictionary with file info extracted from FPR for a given project 
    and a given workflow if workflow is speccified. 
            
    Parameters
    ----------
    - provenance (str): Path to File Provenance Report
    - project (str): Project name as it appears in File Provenance Report. 
    - workflow (list): List of workflows used to generate the output files.
    - prefix (str | None): Prefix used to recover file full paths when File Provevance contains relative paths.
    '''
    
    # create a dict {file_swid: {file info}}
    D  = {}
    
    # get all the records for a single project
    records = get_FPR_records(project, provenance)
    
    # parse the records and get all the files for a given project
    for i in records:
        # keep records for project
        if project == i[1]:
            pipeline_workflow = i[30]
            # check workflow
            if len(workflow) == 1:
                if workflow[0].lower() == 'bcl2fastq':
                    # skip if not fastq-related workflows    
                    if pipeline_workflow.lower() not in ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']:
                        continue
                else:
                    # skip if not provided workflow
                    if workflow[0].lower() != pipeline_workflow.lower():
                        continue
            else:
                if pipeline_workflow.lower() not in list(map(lambda x: x.lower(), workflow)):
                    continue
            
            # get file path
            if prefix:
                file_path = os.path.join(prefix, i[46])
            else:
                file_path = i[46]
            # get file name
            file_name = os.path.basename(file_path)
            # get sample name
            sample_name = i[7]
            # get parent sample name
            parent_sample = i[9].split(':')[0]
            # get time stamp and convert to epoch
            creation_date = i[0]
            # remove milliseconds
            creation_date = creation_date.split('.')[0]
            pattern = '%Y-%m-%d %H:%M:%S'
            creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
            # record platform
            platform = i[22]
            
            # get md5sum
            md5 = i[47]
            # get workdlow swid
            workflow_run_id = i[36]
            # get workflow version
            workflow_version = i[31]
            # get file swid
            file_swid = i[44]
  
            # for merging workflows there will be multiple values for the variables below
            # get library aliquot
            aliquot = i[56].split('_')[-1]
            # get library aliases
            library = i[13]
            # get lims key
            limskey = i[56]
            # get run id
            run_id = i[18]
            # get run
            run = i[23]   
            # get lane
            lane = i[24]
            # get barcode
            barcode = i[27]  
          
            geo = i[12]
            if geo:
                geo = {k.split('=')[0]:k.split('=')[1] for k in geo.split(';')}
            else:
                geo = {}
            for j in ['geo_external_name', 'geo_group_id', 'geo_group_id_description',
                      'geo_targeted_resequencing', 'geo_library_source_template_type',
                      'geo_tissue_type', 'geo_tissue_origin']:
                if j not in geo:
                    geo[j] = 'NA'
                if j == 'geo_group_id':
                    # removes misannotations
                    geo[j] = geo[j].replace('&2011-04-19', '').replace('2011-04-19&', '')
       
            read_count = i[45]
            if read_count:
                read_count = {k.split('=')[0]:k.split('=')[1] for k in i[45].split(';')}
            if 'read_count' in read_count:
                read_count = int(read_count['read_count'])
            else:
                read_count = -1
       
            sample_id = sample_name + '_' + geo['geo_tissue_origin']+ '_' + geo['geo_tissue_type'] + '_' + geo['geo_library_source_template_type'] + '_' + geo['geo_group_id']
         
            d = {'workflow': pipeline_workflow, 'file_path': file_path, 'file_name': file_name,
                 'sample_name': sample_name, 'creation_date': creation_date, 'platform': platform,
                 'md5': md5, 'workflow_run_id': workflow_run_id, 'workflow_version': workflow_version,
                 'file_swid': file_swid, 'external_name': geo['geo_external_name'],
                 'panel': geo['geo_targeted_resequencing'], 'library_source': [geo['geo_library_source_template_type']],
                 'parent_sample': [parent_sample], 'run_id': [run_id], 'run': [run],
                 'limskey': [limskey], 'aliquot': [aliquot], 'library': [library],
                 'barcode': [barcode], 'tissue_type': [geo['geo_tissue_type']],
                 'tissue_origin': [geo['geo_tissue_origin']], 'groupdesc': [geo['geo_group_id_description']],
                 'groupid': [geo['geo_group_id']], 'read_count': read_count, 'sample_id': [sample_id], 'lane': [lane]}
            
            if file_swid not in D:
                D[file_swid] = d
            else:
                assert D[file_swid]['file_path'] == file_path
                assert D[file_swid]['external_name'] == geo['geo_external_name']
                assert D[file_swid]['read_count'] == read_count
                D[file_swid]['sample_id'].append(sample_id)
                D[file_swid]['parent_sample'].append(parent_sample)
                D[file_swid]['run_id'].append(run_id)
                D[file_swid]['run'].append(run)
                D[file_swid]['limskey'].append(limskey)
                D[file_swid]['aliquot'].append(aliquot)
                D[file_swid]['library'].append(library)
                D[file_swid]['barcode'].append(barcode)
                D[file_swid]['tissue_type'].append(geo['geo_tissue_type'])
                D[file_swid]['tissue_origin'].append(geo['geo_tissue_origin'])
                D[file_swid]['library_source'].append(geo['geo_library_source_template_type'])
                D[file_swid]['groupdesc'].append(geo['geo_group_id_description'])
                D[file_swid]['groupid'].append(geo['geo_group_id'])
                D[file_swid]['lane'].append(lane)
    
    return D    
        
            

def select_most_recent_workflow(files):
    '''
    (dict) -> dict
    
    Returns a new dictionary keeping files from the most recent workflow iteration if duplicate files exist
        
    Parameters
    ----------
    - files (dict): Dictionary with file records obtained from parsing FPR
    '''

    # find files with the same file name
    file_names = {}
    for file_swid in files:
        # get the file name
        file_name = files[file_swid]['file_name']
        creation_time = files[file_swid]['creation_date']
        if file_name in file_names:
            # compare creation times
            if creation_time >= file_names[file_name][0]:
                # keep the most recent file
                file_names[file_name] = [creation_time, file_swid]
        else:
            file_names[file_name] = [creation_time, file_swid]
      
    # select the most files
    D = {}
    # make a list of file swids to keep
    to_keep = list(map(lambda x: x[1], list(file_names.values())))
    for file_swid in files:
        if file_swid in to_keep:
            D[file_swid] = files[file_swid]
    return D        
    
    

def exclude_miseq_secords(files):
    '''
    (dict) -> dict
    
    Returns a new dictionary removing file records in files if sequencing was performed on a MiSeq platform.
        
    Parameters
    ----------
    - files (dict): Dictionary with file records obtained from parsing FPR
    '''
    
    D = {}
    exclude = [file_swid for file_swid in files if 'miseq' in files[file_swid]['platform'].lower()]
    for file_swid in files:
        if file_swid not in exclude:
            D[file_swid] = files[file_swid]        
    return D

      
def exclude_non_specified_runs(files, runs):
    '''
    (dict, list) -> dict
    
    Returns a new dictionary removing file records in files if sequencing was performed during a run not specified in runs,
    
    Parameters
    ----------
    - files (dict) : Dictionary with file records obtained from parsing FPR
    - runs (list): List of run ids to keep
    '''    
    
    D = {}
    exclude = [file_swid for file_swid in files if set(files[file_swid]['run_id']).intersection(set(runs)) != set(files[file_swid]['run_id'])]
    
    for file_swid in files:
        if file_swid not in exclude:
            D[file_swid] = files[file_swid]        
    return D
    

def exclude_non_specified_libraries(files, valid_libraries):
    '''
    (dict, dict) -> dict

    Returns a new dictionary removing file records corresponding to libraries not specified in libraries 

    Parameters
    ----------
    - files (dict) : Dictionary with file records obtained from parsing FPR
    - valid_libraries (dict): Dictionary with libraries tagged for release
    '''    
    
    D = {}
    exclude = []

    for file_swid in files:
        libraries = files[file_swid]['library']
        runs = files[file_swid]['run_id']
        lanes = files[file_swid]['lane']
        if len(libraries) != 1 and len(runs) != 1 and len(lanes) != 1:
            sys.exit('Use option -a to link merging-workflows output files')
        library, run, lane = libraries[0], runs[0], int(lanes[0])
        # exclude file if libraries is not included in the input sheet
        if library not in valid_libraries:
            exclude.append(file_swid)
        # exclude file if run is not included in the input sheet
        elif run not in valid_libraries[library]:
            exclude.append(file_swid)
        # exclude file if lane is specified but not included
        elif '' not in valid_libraries[library][run] and lane not in valid_libraries[library][run]:
            exclude.append(file_swid)
    
    exclude = list(set(exclude))    
    
    for file_swid in files:
        if file_swid not in exclude:
            D[file_swid] = files[file_swid]        
    
    return D
        
    
    
def get_libraries_for_non_release(files, exclude):
    '''
    (dict, dict) -> dict

    Returns a dictionary with file records corresponding to libraries tagged for non-release

    Parameters
    ----------
    - files (dict) : Dictionary with file records obtained from parsing FPR
    - exclude (dict): Dictionary with libraries tagged for non-release
    '''    

    # make a list of files to exclude
    L = []
    
    D = {}
    for file_swid in files:
        libraries = files[file_swid]['library']
        runs = files[file_swid]['run_id']
        if len(libraries) != 1 and len(runs) != 1:
            sys.exit('Use option -a to link merging-workflows output files')
        library, run = libraries[0], runs[0]
        if library in exclude and run in exclude[library]:
            L.append(file_swid)
        if file_swid in L:
            D[file_swid] = files[file_swid]
    
    return D
    
    
    
def create_working_dir(project, project_dir, project_name=None):
    '''    
    (str, str, str | None) -> str
    
    Creates and returns path to a sub-directory in project_dir named iether after project or after project_name if defined
        
    Parameters
    ----------
    - project (str): Project name as it appears in File Provenance Report.
    - projects_dir (str): Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/
    - project_name (str | None): Project name used to create the project directory in gsi space
    '''
    
    # use project as project name if not specified
    if project_name:
        name = project_name
    else:
        name = project

    working_dir = os.path.join(project_dir, name)
    
    os.makedirs(working_dir, exist_ok=True)
    return working_dir
    

def generate_links(files, release, project, working_dir, suffix):
    '''
    (dict, bool, str, str, str) -> None
    
    Link fastq files to run directories under the project dir
        
    Parameters
    ----------
    - files_release (dict): Dictionary with file information
    - release (bool): True if files were tagged for release and False otherwise
    - project (str): Name of the project as it appears in FPR
    - working_dir (str): Path to the project sub-directory in GSI space  
    - suffix (str): Indicates fastqs or datafiles
    '''
    
    for file_swid in files:
        if len(files[file_swid]['run_id']) != 1:
            sys.exit('Use parameter -a to specify how files should be linked')
        assert len(files[file_swid]['run_id']) == 1
        run = files[file_swid]['run_id'][0]
        run_name = run + '.{0}.{1}'.format(project, suffix)
        if release == False:
            run_name += '.withold'
        
        run_dir = os.path.join(working_dir, run_name)
        
        os.makedirs(run_dir, exist_ok=True)
        filename = files[file_swid]['file_name']
        link = os.path.join(run_dir, filename)
        file = files[file_swid]['file_path']
        if os.path.isfile(link) == False:
            os.symlink(file, link)



def link_pipeline_data(pipeline_data, working_dir):
    '''
    (dict, str) -> None
    
    Link pipeline files according to hierarchy encoded in pipeline_data structure under working_dir
        
    Parameters
    ----------
    - pipeline_data (dict): Dictionary with files for each sample and workflow
    - working_dir (str): Path to the project sub-directory in GSI space  
    '''
    
    for donor in pipeline_data:
        for sample_name in pipeline_data[donor]:
            # create sample dir
            donor_dir = os.path.join(working_dir, donor)
            sample_dir = os.path.join(donor_dir, sample_name)
            os.makedirs(donor_dir, exist_ok=True)
            os.makedirs(sample_dir, exist_ok=True)
            for workflow_name in pipeline_data[donor][sample_name]:
                # create workflow dir
                workflow_dir = os.path.join(sample_dir, workflow_name)
                os.makedirs(workflow_dir, exist_ok=True)
                # create worfkflow run id dir
                for wf_id in pipeline_data[donor][sample_name][workflow_name]:
                    wfrun_dir = os.path.join(workflow_dir, wf_id)
                    os.makedirs(wfrun_dir, exist_ok=True)
                    # loop over files within each workflow
                    for i in pipeline_data[donor][sample_name][workflow_name][wf_id]:
                        file = i['file_path']
                        filename = i['file_name']
                        link = os.path.join(wfrun_dir, filename)
                        if os.path.isfile(link) == False:
                            os.symlink(file, link)



def exclude_non_specified_files(files, file_names):
    '''
    (dict, list) -> dict
    
    Returns a new dictionary with file records keeping only the files listed in file_names
    
    Parameters
    ----------
    - files (dict) : Dictionary with file records obtained from parsing FPR
    - file_names (list): List of valid file names for release 
    '''
    
    D = {}
    for file_swid in files:
        file_name = files[file_swid]['file_name']
        if file_name in file_names:
            D[file_swid] = files[file_swid]
    return D


def collect_files_for_release(files, release_files, nomiseq, runs, libraries, exclude):
    '''
    (dict, str | None, bool, list | None, str | None, str | None) -> (dict, dict)
    
    Returns dictionaries with file records for files tagged for release and files that should not be released, if any, or an empty dictionary.
        
    Parameters
    ----------
    - files (str): Dictionary with file records for a entire project or a specific workflow for a given project
                        Default is '/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz'
    - release_files (str | None): Path to file with file names to be released 
    - nomiseq (bool): Exclude MiSeq runs if True
    - runs (list | None): List of run IDs. Include one or more run Id separated by white space.
                          Other runs are ignored if provided
    - libraries (str | None): Path to 2 or 3 columns tab-delimited file with library IDs.
                              The first column is always the library alias.
                              The second column is the lane number or the run identifier.
                              The last column is always the run identifier
                              Only the samples with these library aliases are used if provided'
    - exclude (str | None): File with sample name or libraries to exclude from the release
    '''
    
    # keep most recent workflows
    files = select_most_recent_workflow(files)
       
    # check if a list of valid file names is provided
    if release_files:
        infile = open(release_files)
        file_names = infile.read().rstrip().split('\n')
        file_names = list(map(lambda x: os.path.basename(x), file_names))
        infile.close()
    else:
        file_names = []
    if file_names:
        files = exclude_non_specified_files(files, file_names)
    
    # remove files sequenced on miseq if miseq runs are excluded
    if nomiseq:
        files = exclude_miseq_secords(files)
    
    # keep only files from specified runs
    if runs:
        files = exclude_non_specified_runs(files, runs)
    
    # keep only files for specified libraries
    # parse library file if exists 
    valid_libraries = get_libraries(libraries) if libraries else {}
    if valid_libraries:
        files = exclude_non_specified_libraries(files, valid_libraries)
        
    # find files corresponding to libraries tagged for non-release
    excluded_libraries = get_libraries(exclude) if exclude else {}
    if excluded_libraries:
        files_non_release = get_libraries_for_non_release(files, excluded_libraries)
    else:
        files_non_release = {}
    
    return files, files_non_release


def get_pipeline_data(data_structure, files):
    '''
    (dict, dict) -> dict   
    
    Returns a dictionary with the files for each sample and workflow specified in the data_structure
    
    Parameters
    ----------
    - data_structure (str): Dictionary with samples, workflows and workflow_run_id hierarchical structure
    - files (dict): Dictionary with all file records for a given project extracted from FPR
    '''
    
    
    D = {}
     
    for donor in data_structure:
        for sample_id in data_structure[donor]:
            for file_swid in files:
                sample = files[file_swid]['sample_name']
                workflow = files[file_swid]['workflow']
                version = files[file_swid]['workflow_version']
                wf_id = files[file_swid]['workflow_run_id']
                if donor == sample and sample in sample_id and workflow in data_structure[donor][sample_id]:
                    for d in data_structure[donor][sample_id][workflow]:
                        if version == d['workflow_version'] and wf_id == d['workflow_id']:
                            if donor not in D:
                                D[donor] = {}
                            if sample_id not in D[donor]:
                                D[donor][sample_id] = {}
                            # check if workflow name needs to be replaced
                            if 'name' in d:
                                workflow_name = d['name']
                            else:
                                workflow_name = workflow
                            if workflow_name not in D[donor][sample_id]:
                                D[donor][sample_id][workflow_name] = {}
                            if wf_id not in D[donor][sample_id][workflow_name]:
                                D[donor][sample_id][workflow_name][wf_id] = []
                            # check which files are collected
                            if 'extension' in d:
                                # files are collected based on file extension
                                # get the file extension
                                extension = pathlib.Path(files[file_swid]['file_path']).suffix
                                if extension in d['extension']:
                                    D[donor][sample_id][workflow_name][wf_id].append({'file_path': files[file_swid]['file_path'],
                                                                               'md5': files[file_swid]['md5'],
                                                                               'file_name': os.path.basename(files[file_swid]['file_path'])})
                            elif 'files' in d:
                                # files are collected based on file name
                                if os.path.basename(files[file_swid]['file_path']) in list(map(lambda x: os.path.basename(x),  d['files'])):
                                    D[donor][sample_id][workflow_name][wf_id].append({'file_path': files[file_swid]['file_path'],
                                                                               'md5': files[file_swid]['md5'],
                                                                               'file_name': os.path.basename(files[file_swid]['file_path'])})
                            elif 'rename_files' in d:
                                # files are collected based on file name and renamed
                                # make a list of file names
                                file_paths, file_names = [], []
                                for i in d['rename_files']:
                                    file_paths.append(os.path.basename(i['file_path']))
                                    file_names.append(i['file_name'])
                                if os.path.basename(files[file_swid]['file_path']) in file_paths:
                                    # collect file path, md5 and new file name used to name the link
                                    j = file_paths.index(os.path.basename(files[file_swid]['file_path']))
                                    D[donor][sample_id][workflow_name][wf_id].append({'file_path': files[file_swid]['file_path'],
                                                                     'md5': files[file_swid]['md5'],
                                                                     'file_name': file_names[j]})
                            else:
                                D[donor][sample_id][workflow_name][wf_id].append({'file_path': files[file_swid]['file_path'],
                                                                           'md5': files[file_swid]['md5'],
                                                                           'file_name': os.path.basename(files[file_swid]['file_path'])})
    
    return D                        
                        
   
def write_md5sum(data, outputfile):
    '''
    (dict, str) -> None   
    
    Write a file in working_dir with md5sums for all files contained in data
    
    Parameters
    ----------
    - data (dict): Dictionary holding data from a single workflow
    - outpufile (str): Path to the outputfile with md5sums
    '''    
    
    newfile = open(outputfile, 'w')
    for file_swid in data:
        md5 = data[file_swid]['md5']
        file_path = data[file_swid]['file_path']
        newfile.write('\t'.join([file_path, md5]) + '\n')
    newfile.close()


def write_pipeline_md5sum(data, outputfile):
    '''
    (dict, str) -> None   
    
    Write a file in working_dir with md5sums for all files contained in data
    
    Parameters
    ----------
    - data (dict): Dictionary holding pipeline data
    - outpufile (str): Path to the outputfile with md5sums
    '''    
    
    newfile = open(outputfile, 'w')
    # pipeline data
    for donor in data:
        for sample in data[donor]:
            for workflow_name in data[donor][sample]:
                for workflow_run in data[donor][sample][workflow_name]:
                    for i in data[donor][sample][workflow_name][workflow_run]:
                        newfile.write('\t'.join([i['file_path'], i['md5']]) + '\n')
    newfile.close()

def link_files(args):
    '''
    (str, str, str, str | None, str | None, bool, list | None, str | None, str | None, str, str, str, str | None)
    
    Parameters
    ----------
    - provenance (str): Path to File Provenance Report.
                        Default is /.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz
    - project (str): Project name as it appears in File Provenance Report.
    - workflow (str): Worflow used to generate the output files
    - prefix (str | None): Use of prefix assumes that file paths in File Provenance Report are relative paths.
                           Prefix is added to the relative path in FPR to determine the full file path.
    - release_files (str | None): Path to file with file names to be released 
    - nomiseq (bool): Exclude MiSeq runs if True
    - runs (list | None): List of run IDs. Include one or more run Id separated by white space.
                          Other runs are ignored if provided
    - libraries (str | None): Path to 2 columns tab-delimited file with library IDs.
                              The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM).
                              The second column is the run ID. Only the samples with these library aliases are used if provided
    - exclude (str | None): File with sample name or libraries to exclude from the release
    - suffix (str): Indicates map for fastqs or datafiles in the output file name
    - project_name (str): Project name used to create the project directory in gsi space
    - projects_dir (str): Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/
    - analysis (str): Path to file with json structure of samples, workflows, workflow run ids hierarchy
    '''
    
    if args.runs and args.libraries:
        sys.exit('-r and -l are exclusive parameters')

    if args.analysis and args.workflow:
        sys.exit('-w and -a are exclusive parameters. Specify a single worklow or use a json structure to collect files')
    if not args.analysis:
        if not args.workflow:
            sys.exit('Use -w to specify the workflow')
    if not args.workflow:
        if not args.analysis:
            sys.exit('Use -a to indicate the pipeline workflows')
   
    # create a working directory to link files and save md5sums 
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)

    if args.workflow:
        if not args.suffix:
            sys.exit('Suffix -s is required')
        # parse FPR records
        files = parse_fpr_records(provenance, args.project, [args.workflow], args.prefix)
        print('Extracted files from File Provenance Report')
        # get file information for release and eventually for files that should not be released
        files, files_non_release = collect_files_for_release(files, args.release_files, args.nomiseq, args.runs, args.libraries, args.exclude)
                
        # link files to project dir
        if args.suffix == 'fastqs':
            assert args.workflow.lower() in ['bcl2fastq', 'casava', 'fileimport', 'fileimportforanalysis', 'import_fastq']
        generate_links(files, True, args.project, working_dir, args.suffix)
        # generate links for files to be witheld from release
        if files_non_release:
            generate_links(files_non_release, False, args.project, working_dir, args.suffix)
    
    elif args.analysis:
        infile = open(args.analysis)
        data_structure = json.load(infile)
        infile.close()
        
        # parse FPR records
        # make a list of workflows
        workflows = []
        for i in data_structure:
            for j in data_structure[i]:
                workflows.extend(list(data_structure[i][j].keys()))
        
        workflows = list(set(workflows))    
        
        print('workflows', workflows)
        
        
        files = parse_fpr_records(provenance, args.project, workflows, args.prefix)
        print('Extracted files from File Provenance Report')
        pipeline_data = get_pipeline_data(data_structure, files)
        link_pipeline_data(pipeline_data, working_dir)
    
    # write summary md5sums
    # create a dictionary {run: [md5sum, file_path]}
    current_time = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
    if args.workflow:
        outputfile = os.path.join(working_dir, '{0}.release.{1}.{2}.md5sums'.format(args.project, current_time, args.suffix))
        write_md5sum(files, outputfile)
    elif args.analysis:
        outputfile = os.path.join(working_dir, '{0}.release.{1}.pipeline.md5sums'.format(args.project, current_time))
        write_pipeline_md5sum(pipeline_data, outputfile)
    print('Files were extracted from FPR {0}'.format(provenance))


def group_sample_info_mapping(files, add_panel):
    '''
    (dict, bool) -> dict    
    
    Returns a dictionary of sample information organized by run
        
    Parameters
    ----------
    - files (dict) : Dictionary with file records obtained from parsing FPR
    - add_panel (bool): Add panel if True
    '''

    # group info by run id
    D = {}
    
    for file_swid in files:
        file_path = files[file_swid]['file_path']
        assert len(files[file_swid]['run']) == 1
        run = files[file_swid]['run'][0]
        assert len(files[file_swid]['run_id']) == 1
        run_id = files[file_swid]['run_id'][0]
        sample = files[file_swid]['sample_name']
        assert len(files[file_swid]['library']) == 1
        library = files[file_swid]['library'][0]
        assert len(files[file_swid]['library_source']) == 1
        library_source = files[file_swid]['library_source'][0]
        assert len(files[file_swid]['tissue_type']) == 1
        tissue_type = files[file_swid]['tissue_type'][0]
        assert len(files[file_swid]['tissue_origin']) == 1
        tissue_origin = files[file_swid]['tissue_origin'][0]
        assert len(files[file_swid]['barcode']) == 1
        barcode = files[file_swid]['barcode'][0]
        external_id = files[file_swid]['external_name']
        assert len(files[file_swid]['groupdesc']) == 1
        group_description = files[file_swid]['groupdesc'][0]
        assert len(files[file_swid]['groupid']) == 1
        group_id = files[file_swid]['groupid'][0]
        panel = files[file_swid]['panel']
        
        
        # group files by sample, library, run and lane
        key = '_'.join([sample, library, run, barcode])
        
        
        L = [sample, external_id, library, library_source, tissue_type, tissue_origin, run, barcode, group_id, group_description]
        if add_panel:
            L.append(panel)
            
        if key not in D:
            D[key] = {'info': L, 'files': [file_path]}
        else:
            assert D[key]['info'] == L
            D[key]['files'].append(file_path)
        
    return D            



def write_sample_map(project, files, working_dir, add_panel):
    '''
    (str, dict, str, bool) -> None
    
    Writes a sample map in the working directory about the released files 
    for the project of interest
    
    Parameters
    ----------
    - project (str): Project of interest
    - files (dict): Dictionary with information about the released files
    - working_dir (str): Directory where the sample map is written
    - add_panel (bool): Adds panel in the sample map if True
    '''
        
    # group sample information by run            
    sample_info = group_sample_info_mapping(files, add_panel)        
   
    # write sample maps
    current_time = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
    outputfile = os.path.join(working_dir, '{0}.release.{1}.{2}.map.tsv'.format(project, current_time, 'fastqs'))
    newfile = open(outputfile, 'w')
    header = ['sample', 'sample_id', 'library', 'library_source', 'tissue_type', 'tissue_origin', 'run', 'barcode', 'group_id', 'group_description', 'files']
    if add_panel:
        #header.append('panel')
        header.insert(-1, 'panel')
    newfile.write('\t'.join(header) + '\n')
    # make a list of sorted samples
    keys = sorted(list(sample_info.keys()))
    for k in keys:
        files = ';'.join(list(map(lambda x: os.path.basename(x), sample_info[k]['files'])))
        info = sample_info[k]['info']
        info.append(files)        
        newfile.write('\t'.join(list(map(lambda x: str(x), info))) + '\n')
    newfile.close()
    

def map_external_ids(args):
    '''
    (str, str, str | None, str | None, bool, list | None, str | None, str | None, str, str, bool) -> None
        
    Generate sample maps with sample and sequencing information
    
    Parameters
    ----------    
    - provenance (str): Path to File Provenance Report.
    - project (str): Project name as it appears in File Provenance Report.
    - prefix (str | None): Use of prefix assumes that file paths in File Provenance Report are relative paths.
                           Prefix is added to the relative path in FPR to determine the full file path.
    - release_files (str | None): Path to file with file names to be released 
    - nomiseq (bool): Exclude MiSeq runs if True
    - runs (list | None): List of run IDs. Include one or more run Id separated by white space.
                          Other runs are ignored if provided
    - libraries (str | None): Path to 1 or 2 columns tab-delimited file with library IDs.
                              The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM).
                              The second and optional column is the library aliquot ID (eg. LDI32439).
                              Only the samples with these library aliases are used if provided'
    - exclude (str | None): File with sample name or libraries to exclude from the release
    - project_name (str): Project name used to create the project directory in gsi space
    - projects_dir (str): Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/
    - add_panel (bool): Add panel column to sample map if True
    '''

    if args.runs and args.libraries:
        sys.exit('-r and -l are exclusive parameters')
    
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)

    # sample maps are generated using fastq-generating workflows
    workflow, suffix = 'bcl2fastq', 'fastqs'
    
    # create a working directory to link files and save md5sums 
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # parse FPR records
    files = parse_fpr_records(provenance, args.project, [workflow], args.prefix)
    print('Extracted files from File Provenance Report')
    
    # get raw sequence file info
    files, files_non_release = collect_files_for_release(files, args.release_files, args.nomiseq, args.runs, args.libraries, args.exclude)
    
    # write sample map
    write_sample_map(args.project, files, working_dir, args.add_panel)
    print('Generated sample maps in {0}'.format(working_dir))
    print('Information was extracted from FPR {0}'.format(provenance))


def list_files_release_folder(directories):
    '''
    (list) -> list
    
    Returns a list of full file paths for all fastqs located in each folder of the list of directories
        
    Parameters
    ----------
    - directories (list): A list of full paths to directories containining released fastqs 
    '''

    # make a list of files
    files = []
    for run in directories:
        # make a list of links for the directory
        links = [os.path.join(run, i) for i in os.listdir(run) if os.path.isfile(os.path.join(run, i)) and '.fastq.gz' in i]
        files.extend(links)
        # expect only paired fastqs
        if len(links) % 2 != 0:
            sys.exit('Odd number of fastqs in {0}: {1}'.format(run, len(links)))
    return files


def resolve_links(filenames):
    '''
    (list) -> list

    Returns a list of file paths with links resolved

    Parameters
    ----------
    - filenames(list): List of file paths
    '''
    
    files = [os.path.realpath(i) for i in filenames]
    return files


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

                                           
            
def count_released_fastqs_by_instrument(FPR_info, reads):
    '''
    (dict, str) -> dict
    
    Returns the count of released fastqs for each run and instrument
    
    Parameters
    ----------
    - FPR_info (dict): Information about the released fastqs collected from File Provenance Report
    - reads (str): Count all files (read= all) or fastq pairs (read= read1)
    '''
        
    # count released fastqs by instrument and run
    D = {}
    for file in FPR_info:
        instrument = FPR_info[file]['platform']
        assert len(FPR_info[file]['run_id']) == 1
        run = FPR_info[file]['run_id'][0]
        if instrument not in D:
            D[instrument] = {}
        if reads == 'all_reads':
            if run not in D[instrument]:
                D[instrument][run] = 1
            else:
                D[instrument][run] += 1
        elif reads == 'read1':
            if 'R1' in FPR_info[file]['file_path']:
                if run not in D[instrument]:
                    D[instrument][run] = 1
                else:
                    D[instrument][run] += 1
    return D
    


def count_released_fastqs_by_library_type_instrument(FPR_info):
    '''
    (dict, str) -> dict
    
    Returns the count of released fastqs for each library type, run and instrument
    Precondition: Fastqs are paired, return the count of R1
    
    Parameters
    ----------
    - FPR_info (dict): Information about the released fastqs collected from File Provenance Report
    '''
        
    # count released fastqs by instrument and run
    D = {}
    for file in FPR_info:
        instrument = FPR_info[file]['platform']
        assert len(FPR_info[file]['run_id']) == 1
        run = FPR_info[file]['run_id'][0]
        library_type = FPR_info[file]['library_source'][0]
        if library_type not in D:
            D[library_type] = {}
        if instrument not in D[library_type]:
            D[library_type][instrument] = {}
        if 'R1' in FPR_info[file]['file_path']:
            if run not in D[library_type][instrument]:
                D[library_type][instrument][run] = 1
            else:
                D[library_type][instrument][run] += 1
    return D



def list_released_fastqs_project(api, project):
    '''
    (str, str) -> list
    
    Returns a list of fastqs for a given project that were previously released by interrogating QC status in Nabu
    Pre-condition: Released files need to be marked in Nabu
        
    Parameters
    ----------
    - api (str): URL of the nabu API
    - project (str): file_swid (str): File unique identifier
    '''
    
    # get end-point
    api += 'get-fileqcs' if api[-1] == '/' else '/get-fileqcs'
    
    R = []
 
    # check each fastq-generating workflow
    headers = {'accept': 'application/json','Content-Type': 'application/json'}
    json_data = {"project": "{0}".format(project)}
    response = requests.post(api, headers=headers, json=json_data)
    # check response code
    if response.status_code == 200:
        L = response.json()['fileqcs']
        if L:
            for i in L:
                #file_swid = i['fileid']
                qc_status = i['qcstatus']
                file_path = i['filepath']
                if qc_status.upper() == 'PASS':
                    R.append(file_path)
    return R    


def get_tissue_types():
    '''
    (None) -> dict
    
    Returns a dictionary mapping Tissue Type codes to their definitions
    Pre-condition: These definitions are obtained from the Configuration tab in MISO
    '''
    
    D = {'X': 'Xenograft derived from some tumour. Note: may not necessarily be a mouse xenograft',
         'U': 'Unspecified',
         'T': 'Unclassifed tumour',
         'S': 'Serum from blood where clotting proteins have been removed',
         'R': 'Reference or non-tumour, non-diseased tissue sample. Typically used as a donor-specific comparison to a diseased tissue, usually a cancer',
         'P': 'Primary tumour',
         'O': 'Organoid',
         'n': 'Unknown',
         'M': 'Metastatic tumour',
         'F': 'Fibroblast cells',
         'E': 'Endothelial cells',
         'C': 'Cell line derived from a tumour',
         'B': 'Benign tumour',
         'A': 'Cells taken from Ascites fluid'}
    
    return D

def get_tissue_origin():
    '''
    (None) -> dict
    
    Returns a dictionary mapping Tissue Origin codes to their definitions
    Pre-condition: These definitions are obtained from the Configuration tab in MISO
    '''
    
    D = {'Ab': 'Abdomen', 'Ad': 'Adipose', 'Ae': 'Adnexa', 'Ag': 'Adrenal', 'An': 'Anus',
         'Ao': 'Anorectal', 'Ap': 'Appendix', 'As': 'Ascites', 'At': 'Astrocytoma', 'Av': 'Ampulla',
         'Ax': 'Axillary', 'Ba': 'Back', 'Bd': 'Bile', 'Bi': 'Biliary', 'Bl': 'Bladder',
         'Bm': 'Bone', 'Bn': 'Brain', 'Bo': 'Bone', 'Br': 'Breast', 'Bu': 'Buccal',
         'Bw': 'Bowel', 'Cb': 'Cord', 'Cc': 'Cecum', 'Ce': 'Cervix', 'Cf': 'Cell-Free', 'Ch': 'Chest',
         'Cj': 'Conjunctiva', 'Ck': 'Cheek', 'Cn': 'Central', 'Co': 'Colon', 'Cr': 'Colorectal',
         'Cs': 'Cul-de-sac', 'Ct': 'Circulating', 'Di': 'Diaphragm', 'Du': 'Duodenum',
         'En': 'Endometrial', 'Ep': 'Epidural', 'Es': 'Esophagus', 'Ey': 'Eye', 'Fa': 'Fallopian',
         'Fb': 'Fibroid', 'Fs': 'Foreskin', 'Ft': 'Foot', 'Ga': 'Gastric', 'Gb': 'Gallbladder',
         'Ge': 'Gastroesophageal', 'Gi': 'Gastrointestinal', 'Gj': 'Gastrojejunal', 'Gn': 'Gingiva',
         'Gt': 'Genital', 'Hp': 'Hypopharynx', 'Hr': 'Heart', 'Ic': 'ileocecum', 'Il': 'Ileum',
         'Ki': 'Kidney', 'La': 'Lacrimal', 'Lb': 'Limb', 'Le': 'Leukocyte', 'Lg': 'Leg',
         'Li': 'Large', 'Ln': 'Lymph', 'Lp': 'Lymphoblast', 'Lu': 'Lung', 'Lv': 'Liver', 
         'Lx': 'Larynx', 'Ly': 'Lymphocyte', 'Md': 'Mediastinum', 'Me': 'Mesenchyme', 'Mn': 'Mandible',
         'Mo': 'Mouth', 'Ms': 'Mesentary', 'Mu': 'Muscle', 'Mx': 'Maxilla', 'Nk': 'Neck',
         'nn': 'Unknown', 'No': 'Nose', 'Np': 'Nasopharynx', 'Oc': 'Oral', 'Om': 'Omentum',
         'Or': 'Orbit', 'Ov': 'Ovary', 'Pa': 'Pancreas', 'Pb': 'Peripheral', 'Pc': 'Pancreatobiliary',
         'Pd': 'Parathyroid', 'Pe': 'Pelvic', 'Pg': 'Parotid', 'Ph': 'Paratracheal', 'Pi': 'Penis',
         'Pl': 'Plasma', 'Pm': 'Peritoneum', 'Pn': 'Peripheral', 'Po': 'Peri-aorta', 'Pr': 'Prostate',
         'Pt': 'Palate', 'Pu': 'Pleura', 'Py': 'periampullary', 'Ra': 'Right', 'Rc': 'Rectosigmoid',
         'Re': 'Rectum', 'Ri': 'Rib', 'Rp': 'Retroperitoneum', 'Sa': 'Saliva', 'Sb': 'Small',
         'Sc': 'Scalp', 'Se': 'Serum', 'Sg': 'Salivary', 'Si': 'Small', 'Sk': 'Skin', 'Sm': 'Skeletal',
         'Sn': 'Spine', 'So': 'Soft', 'Sp': 'Spleen', 'Sr': 'Serosa', 'Ss': 'Sinus', 'St': 'Stomach',
         'Su': 'Sternum', 'Ta': 'Tail', 'Te': 'Testes', 'Tg': 'Thymic', 'Th': 'Thymus',
         'Tn': 'Tonsil', 'To': 'Throat', 'Tr': 'Trachea', 'Tu': 'Tongue', 'Ty': 'Thyroid',
         'Uc': 'Urachus', 'Ue': 'Ureter', 'Um': 'Umbilical', 'Up': 'Urine', 'Ur': 'Urethra',
         'Us': 'Urine', 'Ut': 'Uterus', 'Uw': 'Urine', 'Vg': 'Vagina', 'Vu': 'Vulva', 'Wm': 'Worm'}

    return D


def get_library_design():
    '''
    (None) -> dict
    
    Returns a dictionary mapping Library Design codes to their definitions
    Pre-condition: These definitions are obtained from the Configuration tab in MISO
    '''

    D = {'WT': 'Whole Transcriptome', 'WG': 'Whole Genome', 'TS': 'Targeted Sequencing',
         'TR': 'Total RNA', 'SW': 'Shallow Whole Genome', 'SM': 'smRNA', 'SC': 'Single Cell',
         'NN': 'Unknown', 'MR': 'mRNA', 'EX': 'Exome', 'CT': 'ctDNA', 'CM': 'cfMEDIP',
         'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq',
         'PG': 'Plasma Whole Genome', 'MC': 'Methylation Detection (C to T) ctDNA',
         'MG': 'Methylation Detection (C to T) Whole Genome'}
       
    return D




  
 
def extract_bamqc_data(bamqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the bamqc table of qc-etl
        
    Parameters
    ----------
    - bamqc_db (str): Path to the bamqc SQLite database generated by qc-etl
    '''
            
    conn = sqlite3.connect(bamqc_db)
    cur = conn.cursor()
    
    # get table name
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    table_name = [i[0] for i in cur][0]
    # get all data from table
    conn.row_factory = sqlite3.Row
    data = conn.execute('select * from {0}'.format(table_name)).fetchall()
        
    columns = ['sample', 'Pinery Lims ID', 'Run Alias', 'instrument', 'library', 'Barcodes', 
               'coverage', 'coverage deduplicated', 'mark duplicates_PERCENT_DUPLICATION',
               'mapped reads', 'total bases on target', 'total reads', 'bases mapped', 'Lane Number']
    
    D = {}
    for i in data:
        d = {}
        for j in columns:
            if j in ['coverage', 'coverage deduplicated', 'mark duplicates_PERCENT_DUPLICATION']:
                d[j] = float(i[j])        
            elif j in ['mapped reads', 'total bases on target', 'total reads', 'bases mapped', 'Lane Number']:
                d[j] = int(i[j])  
            else:
                d[j] = i[j]
                
        # compute on_target rate, not available through qc-etl
        d['on_target'] = compute_on_target_rate(d['bases mapped'], d['total bases on target']) 
                        
        run = i['Run Alias']
        sample = i['Pinery Lims ID']
        
        if run not in D:
            D[run] = {}
        if sample not in D[run]:
            D[run][sample] = []
        D[run][sample].append(d)
    return D
     

def add_bamqc_metrics(FPR_info, file_swid, bamqc_info):
    '''
    (dict, str, dict) -> None
    
    Update the information obtained from File Provenance Report in place with QC information
    collected from bamqc for a given file determined by file_swid if appropriate library 
       
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - file_swid (str): Unique file identifier
    - bamqc_info (dict): QC information for each paired fastq from the bamqc table
    '''
    
    qc_found = False
    limskey = FPR_info[file_swid]['limskey'][0]
    run_alias = FPR_info[file_swid]['run_id'][0]
    sample_id = FPR_info[file_swid]['sample_id'][0]
    barcode = FPR_info[file_swid]['barcode'][0]
    lane = FPR_info[file_swid]['lane'][0]
    library_source = FPR_info[file_swid]['library_source'][0]
    instrument = FPR_info[file_swid]['platform'].replace('_', ' ')
    library = FPR_info[file_swid]['library'][0]
    groupid = FPR_info[file_swid]['groupid'][0]
    
    excluded_libraries = ['CM', 'WT', 'MC', 'MG']
    
    # check that run in recorded in bamqc
    if library_source not in excluded_libraries and run_alias in bamqc_info:
        if limskey in bamqc_info[run_alias]:
            # map file info with bamqc info
            for d in bamqc_info[run_alias][limskey]:
                if d['Pinery Lims ID'] == limskey:
                    assert '_'.join(sample_id.split('_')[:2]) == '_'.join(d['sample'].split('_')[:2])
                    if groupid == 'NA':
                        assert ('_'.join(sample_id.split('_')[:-1]) == d['sample'] or library == d['sample'])
                    assert int(lane) == int(d['Lane Number'])
                    assert barcode == d['Barcodes']
                    assert instrument  == d['instrument']
                    assert library == d['library']    
                    qc_found = True
                    FPR_info[file_swid]['coverage'] = round(d['coverage'], 2)
                    FPR_info[file_swid]['coverage_dedup'] = round(d['coverage deduplicated'], 2)
                    FPR_info[file_swid]['on_target'] = round(d['on_target'], 2)                
                    FPR_info[file_swid]['percent_duplicate'] = round(d['mark duplicates_PERCENT_DUPLICATION'], 2)
    if library_source not in excluded_libraries and qc_found == False:
        FPR_info[file_swid]['coverage'] = 'NA'
        FPR_info[file_swid]['coverage_dedup'] = 'NA'
        FPR_info[file_swid]['on_target'] = 'NA'                
        FPR_info[file_swid]['percent_duplicate'] = 'NA'



def add_cfmedipqc_metrics(FPR_info, file_swid, cfmedipqc_info):
    '''
    (dict, str, dict) -> None
    
    Update the information obtained from File Provenance Report in place with QC information
    collected from cfmedipqc for a given file determined by file_swid if library source is CM
    
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - file_swid (str): Unique file identifier
    - cfmedipqc_info (dict): QC information for each paired fastq from the cfmedip QC db
    '''
    
    qc_found = False
    run_alias = FPR_info[file_swid]['run_id'][0]
    limskey = FPR_info[file_swid]['limskey'][0]
    barcode = FPR_info[file_swid]['barcode'][0]
    lane = FPR_info[file_swid]['lane'][0]
    library_source = FPR_info[file_swid]['library_source'][0]
    
    
    # # check that run in recorded in cfmedipqc
    if library_source == 'CM' and run_alias in cfmedipqc_info:
        if limskey in cfmedipqc_info[run_alias]:
            for d in cfmedipqc_info[run_alias][limskey]:
                if d['Pinery Lims ID'] == limskey:
                    assert int(d['Lane Number']) == int(lane)
                    assert d['Barcodes'] == barcode
                    assert d['Run Alias'] == run_alias
                    qc_found = True
                    FPR_info[file_swid]['AT_dropout'] = d['AT Dropout']
                    FPR_info[file_swid]['methylation_beta'] = d['Methylation beta']
                    FPR_info[file_swid]['duplication'] = d['Percent Duplication']
                    FPR_info[file_swid]['CpG_enrichment'] = d['Relative CpG Frequency in Regions']
    if library_source == 'CM' and qc_found == False:
        FPR_info[file_swid]['AT_dropout'] = 'NA'
        FPR_info[file_swid]['methylation_beta'] = 'NA'
        FPR_info[file_swid]['duplication'] = 'NA'
        FPR_info[file_swid]['CpG_enrichment'] = 'NA'
 


def add_rnaseqqc_metrics(FPR_info, file_swid, rnaseqqc_info):
    '''
    (dict, str, dict) -> None
    
    Update the information obtained from File Provenance Report in place with QC information
    collected from cfmedipqc for a given file determined by file_swid if library source is WT
    
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - file_swid (str): Unique file identifier
    - rnaseqqc_info (dict): QC information for each paired fastq from the rnaseqqc QC db
    '''
    
    qc_found = False
    run_alias = FPR_info[file_swid]['run_id'][0]
    limskey = FPR_info[file_swid]['limskey'][0]
    barcode = FPR_info[file_swid]['barcode'][0]
    lane = FPR_info[file_swid]['lane'][0]
    library_source = FPR_info[file_swid]['library_source'][0]
    
    # check that run in recorded in rnaseqqc_db
    if library_source == 'WT' and run_alias in rnaseqqc_info:
        if limskey in rnaseqqc_info[run_alias]:
            for d in rnaseqqc_info[run_alias][limskey]:
                if d['Pinery Lims ID'] == limskey:
                    assert int(d['Lane Number']) == int(lane)
                    assert d['Barcodes'] == barcode
                    assert d['Run Alias'] == run_alias
                    qc_found = True
                    FPR_info[file_swid]["5'-3' bias"] = d['MEDIAN_5PRIME_TO_3PRIME_BIAS']
                    FPR_info[file_swid]['rRNA contamination'] = round((d['rrna contamination properly paired'] / d['rrna contamination in total (QC-passed reads + QC-failed reads)'] * 100), 3)
                    FPR_info[file_swid]['Coding (%)'] = d['PCT_CODING_BASES']
                    FPR_info[file_swid]['Correct strand reads (%)'] = d['PCT_CORRECT_STRAND_READS']
    if library_source == 'WT' and qc_found == False:
        FPR_info[file_swid]["5'-3' bias"] = 'NA'
        FPR_info[file_swid]['rRNA contamination'] = 'NA'
        FPR_info[file_swid]['Coding (%)'] = 'NA'
        FPR_info[file_swid]['Correct strand reads (%)'] = 'NA'

    



def add_emseqqc_metrics(FPR_info, file_swid, emseqqc_info):
    '''
    (dict, str, dict) -> None
      
    Update the information obtained from File Provenance Report in place with QC information
    collected from emseqqc for a given file determined by file_swid if library source is MC or MG
        
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - file_swid (str): Unique file identifier
    - emseqqc_info (dict): QC information for each paired fastq from the emseq QC db
    '''
        
    qc_found = False
    run_alias = FPR_info[file_swid]['run_id'][0]
    limskey = FPR_info[file_swid]['limskey'][0]
    barcode = FPR_info[file_swid]['barcode'][0]
    lane = FPR_info[file_swid]['lane'][0]
    library_source = FPR_info[file_swid]['library_source'][0]

    # check that run in recorded in emseqc_db
    if library_source in ['MC', 'MG'] and run_alias in emseqqc_info:
        if limskey in emseqqc_info[run_alias]:
            for d in emseqqc_info[run_alias][limskey]:
                if d['Pinery Lims ID'] == limskey:
                    assert int(d['Lane Number']) == int(lane)
                    assert d['Barcodes'] == barcode
                    assert d['Run Alias'] == run_alias
                    qc_found = True
                    FPR_info[file_swid]['Lambda_methylation'] = d['Lambda']
                    FPR_info[file_swid]['pUC19_methylation'] = d['pUC19']
                    FPR_info[file_swid]['percent_duplication'] = d['mark duplicates_PERCENT_DUPLICATION']
                    
    if library_source in ['MC', 'MG'] and qc_found == False:
        FPR_info[file_swid]['Lambda_methylation'] = 'NA'
        FPR_info[file_swid]['pUC19_methylation'] = 'NA'
        FPR_info[file_swid]['percent_duplication'] = 'NA'
    

def map_QC_metrics_to_fpr(FPR_info, bamqc_info, cfmedipqc_info, rnaseqqc_info, emseqqc_info):
    '''
    (dict, dict, dict, dict, dict) -> None
    
    Update the information obtained from File Provenance Report in place with information
    collected from the appropriate library source-specific QC db
    
    Parameters
    ----------
    - FPR_info (dict): Information for each released fastq collected from File Provenance Report
    - bamqc_info (dict): QC information for each paired fastq from the bamqc db
    - cfmedipqc_info (dict): QC information for each paired fastq from the cfmedipqc db
    -rnaseqqc_info (dict) QC information for each paired fastqs from the rnaseqqc db
    -emseqqc_info (dict) QC information for each paired fastqs from the emseqqc db
    '''
    
    for file_swid in FPR_info:
        # check library source
        library_source = FPR_info[file_swid]['library_source'][0]
        if library_source == 'CM':
            add_cfmedipqc_metrics(FPR_info, file_swid, cfmedipqc_info)
        elif library_source == 'WT':
            add_rnaseqqc_metrics(FPR_info, file_swid, rnaseqqc_info)
        elif library_source in ['WG', 'EX', 'TS', 'PG']:
            add_bamqc_metrics(FPR_info, file_swid, bamqc_info)
        elif library_source in ['MC', 'MG']:
            add_emseqqc_metrics(FPR_info, file_swid, emseqqc_info)              


def extract_cfmedipqc_data(cfmedipqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the cfmedip table of qc-etl
        
    Parameters
    ----------
    - cfmedipqc_db (str): Path to the cfmedip SQLIte database generated by qc-etl
    '''
    
    #conn = sqlite3.connect('prov_report_test.db')
    conn = sqlite3.connect(cfmedipqc_db)
    conn.row_factory = sqlite3.Row
    
    # get all data
    data = conn.execute('select * from cfmedipqc_cfmedipqc_4').fetchall()

    conn.close()
    
    # get clumns of interest
    columns = ['AT Dropout',
               'Methylation beta',
               'Percent Duplication',          
               'Relative CpG Frequency in Regions',
               'Barcodes',
               'Lane Number',
               'Pinery Lims ID',
               'Run Alias']
    
    D = {}

    for i in data:
        i = dict(i)
        run = i['Run Alias']
        if run not in D:
            D[run] = {}
        sample = i['Pinery Lims ID']
        if sample not in D[run]:
            D[run][sample] = []
        d = {}
        for j in columns:
            try:
                float(i[j])
            except:
                d[j] = i[j]
            else:
                d[j] = round(float(i[j]), 3)
        D[run][sample].append(d)
           
    return D



def extract_rnaseqqc_data(rnaseqqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the rnaseqqc table of qc-etl
        
    Parameters
    ----------
    - rnaseqqc_db (str): Path to the rnaseq SQLIte database generated by qc-etl
    '''
    
    #conn = sqlite3.connect('prov_report_test.db')
    conn = sqlite3.connect(rnaseqqc_db)
    cur = conn.cursor()
    
    # get table name
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    table_name = [i[0] for i in cur][0]
    
    # get all data from table
    conn.row_factory = sqlite3.Row
    data = conn.execute('select * from {0}'.format(table_name)).fetchall()
    conn.close()
    
    # get clumns of interest
    columns = ['Barcodes',
               'Lane Number',
               'Pinery Lims ID',
               'Run Alias',
               'PCT_CODING_BASES',
               'rrna contamination in total (QC-passed reads + QC-failed reads)',
               'rrna contamination properly paired',
               'PCT_CORRECT_STRAND_READS',
               'MEDIAN_5PRIME_TO_3PRIME_BIAS']
        
        
    D = {}

    for i in data:
        i = dict(i)
        run = i['Run Alias']
        if run not in D:
            D[run] = {}
        sample = i['Pinery Lims ID']
        if sample not in D[run]:
            D[run][sample] = []
        d = {}
        for j in columns:
            try:
                float(i[j])
            except:
                d[j] = i[j]
            else:
                d[j] = round(float(i[j]), 3)
        D[run][sample].append(d)
           
    return D


def extract_emseqqc_data(emseqqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the rnaseqqc table of qc-etl
        
    Parameters
    ----------
    - rnaseqqc_db (str): Path to the rnaseq SQLIte database generated by qc-etl
    '''
    
    #conn = sqlite3.connect('prov_report_test.db')
    conn = sqlite3.connect(emseqqc_db)
    conn.row_factory = sqlite3.Row
    # get all methylation data
    data = conn.execute('select * from emseqqc_methylation_2').fetchall()
    conn.close()
    
    # get clumns of interest
    columns = ['Barcodes',
               'Lane Number',
               'Pinery Lims ID',
               'Run Alias',
               'Lambda',
               'pUC19']
               
    D = {}

    for i in data:
        i = dict(i)
        run = i['Run Alias']
        if run not in D:
            D[run] = {}
        sample = i['Pinery Lims ID']
        if sample not in D[run]:
            D[run][sample] = []
        d = {}
        for j in columns:
            try:
                float(i[j])
            except:
                d[j] = i[j]
            else:
                d[j] = round(float(i[j]), 3)
        D[run][sample].append(d)
           
    # add duplication rate 
    conn = sqlite3.connect(emseqqc_db)
    conn.row_factory = sqlite3.Row
    data2 = conn.execute('select * from emseqqc_bamqc_2').fetchall()
    conn.close()
    
    for i in data2:
        i = dict(i)
        run = i['Run Alias']
        sample = i['Pinery Lims ID']
        barcodes = i['Barcodes']
        lane = i['Lane Number']
        duplicate = i['mark duplicates_PERCENT_DUPLICATION']
        
        # find dict in D
        assert run in D
        assert sample in D[run]
        for j in range(len(D[run][sample])):
            if D[run][sample][j]['Barcodes'] == barcodes and D[run][sample][j]['Lane Number'] == lane and \
                D[run][sample][j]['Pinery Lims ID'] == sample and D[run][sample][j]['Run Alias'] == run:
                D[run][sample][j]['mark duplicates_PERCENT_DUPLICATION'] = duplicate    

    return D


def get_file_prefix(file):
    '''
    (str) -> str
    
    Returns the file prefix of a fastq file.
    File prefix is the file name without the file extension and read number
    
    Parameters
    ----------
    - file (str): File path
    '''


    filename = os.path.basename(file)
    x = re.search("[\._]{1}R[1-2][\._]+.*fastq.gz", filename)
    
    assert x
    prefix = filename[:x.span()[0]]
    
    return prefix
     
    
def add_file_prefix(L):
    '''
    (list) -> None

    Adds file prefix to each dictionary of file info.
    L is a list of inner lists with information for paired fastqs
    
    Parameters
    ----------
    - L (list): List of lists with dictionaries of file information
    '''

    for i in L:
        for j in i:
            file = j['file_path']
            prefix = get_file_prefix(file)
            j['prefix'] = prefix


def find_fastq_pairs(files, platform):
    '''
    (dict, str) -> list
    
    Returns a list of 2-item lists with dictionary about file info of paired fastqs.
    Pre-condition: All reads are paired-reads and it exists 2 fastqs for read 1 and read 2
    
    Parameters
    ----------
    - files (dict): Dictionary with file information extracted from FPR 
    '''
     
    # make a list with dictionaries of file info
    file_info = [files[i] for i in files if files[i]['platform'] == platform]
    if len(file_info) % 2 != 0:
        sys.exit('Expecting paired fastqs for platform {0}, but got {1} files'.format(platform, len(file_info)))
    
    # sort according to file name
    file_info.sort(key = lambda x: x['file_path'])
    
    L = []
    for i in range(0, len(file_info), 2):
        assert file_info[i]['run'][0] == file_info[i+1]['run'][0] 
        assert file_info[i]['platform'] == file_info[i+1]['platform'] == platform 
        assert file_info[i]['barcode'][0] == file_info[i+1]['barcode'][0]
        assert file_info[i]['library'][0] == file_info[i+1]['library'][0]
        assert file_info[i]['sample_name'] == file_info[i+1]['sample_name']
        assert file_info[i]['external_name'] == file_info[i+1]['external_name']
        assert file_info[i]['library_source'][0] == file_info[i+1]['library_source'][0]
        assert file_info[i]['lane'][0] == file_info[i+1]['lane'][0] 
        assert file_info[i]['tissue_origin'][0] == file_info[i+1]['tissue_origin'][0] 
        assert file_info[i]['tissue_type'][0] == file_info[i+1]['tissue_type'][0] 
        assert file_info[i]['sample_id'][0] == file_info[i+1]['sample_id'][0] 
        assert file_info[i]['limskey'][0] == file_info[i+1]['limskey'][0]
        assert file_info[i]['read_count'] == file_info[i+1]['read_count']
        L.append([file_info[i], file_info[i+1]])
        
    assert len(L) == len(file_info) / 2
    for i in L:
        assert i[0]['file_path'] != i[1]['file_path']
        if 'R1' in i[0]['file_name']:
            assert 'R2' in i[1]['file_name']
            r1, r2 = 0, 1
        elif 'R2' in i[0]['file_name']:
            assert 'R1' in i[1]['file_name']
            r1, r2 = 1, 0
        assert i[r1]['file_name'][:i[r1]['file_name'].rindex('R1')] == i[r2]['file_name'][:i[r2]['file_name'].rindex('R2')]
    
    return L



def get_run_level_metrics(files, platform, library_source, metrics):
    '''
    (dict, str, str) -> (list, list, list, list, list)
    
    Returns a tuple with parallel lists of run-level metrics for a given instrument.
    Pre-condition: All reads are paired-reads and it exists 2 fastqs for read 1 and read 2
    
    Parameters
    ----------
    - files (dict): Dictionary with file information extracted from FPR 
    - platform (str): Sequencing platform
    - library_source (str): Type of library (eg: CM, WG, WT)
    - metrics (list): List of metrics of interest
    '''
     
    # find fastq pairs
    L = find_fastq_pairs(files, platform)
    # add file prefix
    add_file_prefix(L)
    
    QC_metrics = [[] for i in range(len(metrics))]
    
    for i in L:
        if i[0]['library_source'][0] == library_source:
            for j in range(len(metrics)):
                QC_metrics[j].append(i[0][metrics[j]])
       
    return QC_metrics                 


def clean_up_metrics(metrics):
    '''
    (list) -> list
    
    Returns a list of lists with metrics without missing values NA, keeping the original order of the lists
    
    Parameters
    ----------
    - metrics (list): List of metrics
    '''
    
    while any(list(map(lambda x: 'NA' in x, metrics))):
        for i in metrics:
            if 'NA' in i:
                pos = i.index('NA')
                break
        for i in range(len(metrics)):
            if metrics[i]:
                del metrics[i][pos]
    return metrics 


def sort_metrics(QC_metrics):
    '''
    (list) -> list
    
    Returns a list of lists of metrics preserving the order among lists and 
    and sorted according to the order of reads (ie, 1st metric)
    Pre-condition: There is no missing values and all lists have the same length
        
    Parameters
    ----------
    - QC_metrics (list): List list of metrics. The 1st list are the read counts
    '''
    
    # make a list with inner lists containing each the ith value of each metric
    a = list(zip(*QC_metrics))
    # sort the inner lists according to the read count, it the first value of each inner list
    a.sort(key=lambda x : x[0])
    # unpack the sorted list to get back the list of each metric, order-preserved and read-counts sorted
    if a:
        QC_metrics = list(zip(*a))
    return QC_metrics 


def get_x_axis_labels(data):
    '''
    (list) -> list
    
    Returns a list of x axis labels.
    Only the last defined metrics in data is labeled with Samples
        
    Parameters
    ----------
    - data (list): List of lists with metrics
    '''
    
    
    # make a list with x labels to find out wich subplot should show Samples label
    labels = ['Samples' if data[i] else None for i in range(len(data))]
    max_index = 0
    for i in range(len(labels)):
        if labels[i]:
            max_index = i
    for i in range(len(labels)):
        if i < max_index:
            labels[i] = None
    return labels


def count_subplots(data):
    '''
    (list) -> int
    
    Returns the number of expected subplots in figure.
        
    Parameters
    ----------
    - data (list): List of lists with metrics
    '''
    
    # determine how many subplots are expected with defined metrics
    return sum([1 if i else 0 for i in data])
     


def get_subplot_position(data):
    '''
    (list) -> list
    
    Returns the position of each subplot for metrics in data.
        
    Parameters
    ----------
    - data (list): List of lists with metrics
    '''
    
    subplot_pos = [1 if i else 0 for i in data]
    for i in range(1,len(subplot_pos)):
        subplot_pos[i] = subplot_pos[i-1] + subplot_pos[i]
    return subplot_pos


def create_ax(row, col, pos, figure, Data, YLabel, color, title = None, XLabel = None):
    '''
    (int, int, int, matplotlib.figure.Figure, list, str, str, str | None, str | None)
    
    Parameters
    ----------
    
    - row (int): Row position of the plot in figure
    - col (int): Column position of the plot in figure
    - figure (matplotlib.figure.Figure): Matplotlib figure
    - Data (list): List of metrics to plot in graph
    - YLabel (str): Label of the Y axis
    - color (str): Color of markers for Data
    - title (str | None): Title of the graph
    - XLabel (str | None): Label of the X axis    
    '''
    
    # create ax in figure to plot data
    ax = figure.add_subplot(row, col, pos)
    
    # plot data and median  
    xcoord = [i/10 for i in range(len(Data))]
        
    ax.plot(xcoord, Data, clip_on=False, linestyle='', marker= 'o', markerfacecolor = color, markeredgecolor = color, markeredgewidth = 1, markersize = 10, alpha=0.5)
    # compute median the data
    median = np.median(Data)
    # plot median and mean. use zorder to bring line to background
    ax.axhline(y=median, color=color, linestyle='-', linewidth=1.5, alpha=0.5, zorder=1)
    
    # start y axis at 0
    ax.set_ylim(ymin=0)
       
    # write axis labels
    if XLabel is not None:
        ax.set_xlabel(XLabel, color='black', size=18, ha='center', weight= 'normal')
    ax.set_ylabel(YLabel, color='black', size=18, ha='center', weight='normal')
    
    # add title 
    if title is not None:
        ax.set_title(title, weight='bold', pad =20, fontdict={'fontsize':20})

    # set xticks
    # get all the ticks and set labels to empty str
    plt.xticks(xcoord, ['' for i in range(len(xcoord))], ha='center', fontsize=12, rotation=0)
    # set every N ticks
    N = 3
    xticks_pos = ax.get_xticks()
    xticks_labels = ax.get_xticklabels()
    myticks = [j for i,j in enumerate(xticks_pos) if not i % N]  # index of selected ticks
    newlabels = [label for i,label in enumerate(xticks_labels) if not i % N]
    plt.xticks(myticks, newlabels, ha='center', fontsize=12, rotation=0)
    
    # add splace bewteen axis and tick labels
    ax.yaxis.labelpad = 17
    
    # do not show frame lines  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)    
        
    # offset the x axis
    for loc, spine in ax.spines.items():
        spine.set_position(('outward', 5))
        #spine.set_smart_bounds(True)
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    # disable scientific notation
    ax.ticklabel_format(style='plain', axis='y')
    
    return ax




def generate_figures(files, project, library_source, platform, metrics, Y_axis, colors, working_dir, height=16, width=13):
    '''
    (dict, str, str, str, str, int, int) -> str
    
    Generate a figure with metrics from FPR and QC-etl for a given library type and sequencing platform
    and returns the path to the figure file
        
    Parameters
    ----------
    - files (dict): Dictionary with file info extracted from FPR and with QC info extracted from qc-etl
    - project (str): Name of project
    - library_source (str): Type of library
    - platform (str): Sequencing platform
    - working_dir (str): Path to the folder where figure files are written
    - height (int): Height of the figure
    - width (int): Width of the figure
    '''
    
    # make lists with metrics for each instrument 
    QC_metrics = get_run_level_metrics(files, platform, library_source, metrics)
    # remove undefined metric values
    QC_metrics = clean_up_metrics(QC_metrics)
    # sort metrics according to read counts (read counts is always first metric)
    QC_metrics = sort_metrics(QC_metrics)
        
    # get the outputfile
    current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    outputfile = os.path.join(working_dir, '{0}.{1}.{2}.{3}.QC_plots.png'.format(project, platform, library_source, current_time))
    
    if QC_metrics[0]:
        figure = plt.figure()
        figure.set_size_inches(width, height)
        # make a list of with X axis labels to determine which subplot should display the Samples label
        x_labels = get_x_axis_labels(QC_metrics)
        # determine how many subplots are expected
        subplots = count_subplots(QC_metrics) 
        # determine the position of each subplot
        subplot_pos = get_subplot_position(QC_metrics)
            
        for i in range(len(QC_metrics)):
            # determine title
            title = platform + ' {0} libraries'.format(library_source) if i == 0 else None
            # plot data
            create_ax(subplots, 1, subplot_pos[i], figure, QC_metrics[i], Y_axis[i], colors[i], title = title, XLabel = x_labels[i])
                
        # make sure axes do not overlap
        plt.tight_layout(pad = 2.5)
        # write figure to file  
        figure.savefig(outputfile, bbox_inches = 'tight')
        plt.close()
                
        return outputfile
    else:
        return ''



def get_library_metrics(library_type):
    '''
    (str) -> list
    
    Returns a list of metrics of interest for library_type if library_type is assigned specific metrics
    and returns a list with read count only otherwise
    
    Parameters
    ----------
    - library_type (str): The Library code defined in MISO
    '''
    
    metrics = {'CM': ['read_count', 'methylation_beta', 'CpG_enrichment'],
               'WT': ['read_count', 'rRNA contamination', 'Coding (%)'],
               'WG': ['read_count', 'coverage_dedup'],
               'PG': ['read_count', 'coverage_dedup'],
               'TS': ['read_count', 'coverage_dedup', 'on_target'],
               'EX': ['read_count', 'coverage_dedup', 'on_target'],
               'MC': ['read_count', 'Lambda_methylation', 'pUC19_methylation', 'percent_duplication'],
               'MG': ['read_count', 'Lambda_methylation', 'pUC19_methylation', 'percent_duplication']}
               
    if library_type in metrics:
        return metrics[library_type]
    else:
        return ['read_count']
    


def count_samples_with_missing_values(files):
    '''
    (dict) -> dict
    
    Returns a dictionary with the number of samples with missing metric values for
    each library type and instrument
    
    Parameters
    ----------
    - files (dict): Dictionary with file info extracted from FPR and QC metrics extracted from qc-etl
    '''
    
    # count the number of samples with missing values for each instrument and library type
    D = {}
    for file_swid in files:
        platform = files[file_swid]['platform']
        sample = files[file_swid]['external_name']
        library_source = files[file_swid]['library_source'][0]
        metrics = get_library_metrics(library_source)
                       
        for i in metrics:
            if i in files[file_swid]:
                if files[file_swid][i] == 'NA':
                    if library_source not in D:
                        D[library_source] = {}
                    if platform not in D[library_source]:
                        D[library_source][platform] = set()
                    D[library_source][platform].add(sample)
    
    for library_source in D:
        for platform in D[library_source]:
            D[library_source][platform] = len(list(D[library_source][platform]))
    
    for library_source in D:
        to_remove = [platform for platform in D[library_source] if D[library_source][platform] == 0]
        for platform in to_remove:
            del D[library_source][platform]
    to_remove = [library_source for library_source in D if len(D[library_source]) == 0]
    for i in to_remove:
        del D[i]
    
    return D   


def fit_into_column(text):
    '''
    (str)- > str
    
    
    Returns text reformatted to fit the table column
    
    Parameters
    ----------
    - text (str): Value in column 
    '''
    
    text = '{0} {1}'.format(text[:len(text)//2], text[len(text)//2:])

    return text
       

def group_sample_metrics(files, table, metrics = None):
    '''
    (dict, str, dict | None) -> dict 
    
    Group sample identifiers or metrics 
    
    Parameters
    ----------
    - files (dict): Information about fastqs extracted from File Provenance Report
    - table (str): Table of interest in the report
    - metrics (dict or None): Dictionary with metrics to report or None
    '''
    
    
    # make a list of instruments
    instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
    
    # record sample metrics for each instrument
    D = {}
    
    # find pairs of fastqs
    for platform in instruments:
        pairs = find_fastq_pairs(files, platform)
        add_file_prefix(pairs)
        for i in pairs:
            # add comma to reads for readability
            case = i[0]['sample_name']
            external_name = i[0]['external_name']
            sample = i[0]['sample_id'][0]
            library_source = i[0]['library_source'][0]
            library = i[0]['library'][0]
            prefix = i[0]['prefix']
            groupid = i[0]['groupid'][0]
            group_description = i[0]['groupdesc'][0]
            run = i[0]['run_id']
            
            max_length = 20
            
            
            # reformat prefix to fit the table column
       
            if len(prefix) >= max_length:
                prefix = fit_into_column(prefix)
            if len(library) >= max_length:
                library = fit_into_column(library)
            if len(sample) >= max_length:
                sample = fit_into_column(sample)
            if len(groupid) >= max_length:
                groupid = fit_into_column(groupid)
            if len(group_description) >= max_length:
                group_description = fit_into_column(group_description)
                        
            
            tissue_origin = i[0]['tissue_origin'][0]
            tissue_type = i[0]['tissue_type'][0]
            sequencing_run = '{0} lane_{1}_{2}'.format(i[0]['run_id'][0], i[0]['lane'][0], i[0]['barcode'][0])
                        
            if table == 'sample_identifiers':
                L = [library, case, external_name, groupid, group_description, library_source, tissue_origin, tissue_type]
            elif table == 'qc_metrics':
                assert metrics
                QC_metrics = []
                for metric in metrics[library_source]:
                    if metric == 'read_count':
                        QC_metrics.append('{:,}'.format(i[0]['read_count']))
                    else:
                        QC_metrics.append(i[0][metric])
                L = [library, prefix]
                L.extend(QC_metrics)
                
            if library_source not in D:
                D[library_source] = {}
            if platform not in D[library_source]:
                D[library_source][platform] = []
            if L not in D[library_source][platform]:
                D[library_source][platform].append(L)
    # sort sample info on sample name
    for i in D:
        for j in D[i]:
            D[i][j].sort(key = lambda x: x[0])        
        
    return D            





def group_cumulative_samples(files):
    '''
    (dict) -> list
    
    Returns a list of sample identifiers sorted on donor id 
    
    Parameters
    ----------
    - files (dict): Information about fastqs extracted from File Provenance Report
    '''
        
    # make a list of instruments
    instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
    
    # record sample identifiers across library type and platform
    SL = []
    
    # find pairs of fastqs
    for platform in instruments:
        pairs = find_fastq_pairs(files, platform)
        add_file_prefix(pairs)
        for i in pairs:
            # add comma to reads for readability
            case = i[0]['sample_name']
            external_name = i[0]['external_name']
            sample = i[0]['sample_id'][0]
            library_source = i[0]['library_source'][0]
            library = i[0]['library'][0]
            prefix = i[0]['prefix']
            groupid = i[0]['groupid'][0]
            group_description = i[0]['groupdesc'][0]
            
            
            max_length = 20
            
            # reformat prefix to fit the table column
            if len(prefix) >= max_length:
                prefix = fit_into_column(prefix)
            if len(library) >= max_length:
                library = fit_into_column(library)
            if len(sample) >= max_length:
                sample = fit_into_column(sample)
            if len(groupid) >= max_length:
                groupid = fit_into_column(groupid)
            if len(group_description) >= max_length:
                group_description = fit_into_column(group_description)
                        
            tissue_origin = i[0]['tissue_origin'][0]
            tissue_type = i[0]['tissue_type'][0]
            L = [case, external_name, sample, group_description, library_source, tissue_origin, tissue_type]

            if L not in SL:
                SL.append(L)
    SL.sort(key = lambda x: x[0])
    
    return SL            






def get_library_tissue_types(files):
    '''
    (dict) -> dict
    
    Returns a dictionary with tissue origin, type and library source definitions from MISO
    for each releaed library
    
    Parameters
    ----------
    - files (dict): Dictionary with file information from released libraries extracted from FPR
    '''

    D = {'Library Type': [], 'Tissue Type': [], 'Tissue Origin': []}

    # definitions are from the Configuration tab in MISO
    tissue_types = get_tissue_types()
    tissue_origin = get_tissue_origin()
    library_design = get_library_design()
    
    for file_swid in files:
        D['Library Type'].extend(files[file_swid]['library_source'])
        D['Tissue Type'].extend(files[file_swid]['tissue_type'])
        D['Tissue Origin'].extend(files[file_swid]['tissue_origin'])
    
    for i in D:
        D[i] = sorted(list(set(D[i])))                
       
    for i in range(len(D['Library Type'])):
        D['Library Type'][i] = '{0}: {1}'.format(D['Library Type'][i], library_design[D['Library Type'][i]])
    for i in range(len(D['Tissue Type'])):
        D['Tissue Type'][i] = '{0}: {1}'.format(D['Tissue Type'][i], tissue_types[D['Tissue Type'][i]])
    for i in range(len(D['Tissue Origin'])):
        D['Tissue Origin'][i] = '{0}: {1}'.format(D['Tissue Origin'][i], tissue_origin[D['Tissue Origin'][i]])
    
    for i in D:
        D[i] = ', '.join(D[i])                
    
    return D


def get_identifiers_appendix(files, report):
    '''
    (dict, str) -> list
    
    Returns a list with definitions of columns in the sample identifier table
    
    Parameters
    ----------
    - files (dict): Dictionary with file information from released libraries extracted from FPR
    - report (str): cumulative or batch report
    '''

    # get the library type, tissue type and tissue origin 
    D = get_library_tissue_types(files)
    
    if report  == 'batch':
        L = ['Library Id: OICR-generated library identifier',
             'OICR Case Id: OICR-generated case identifier',
             'Donor Id: user supplied donor identifier',
             'Sample Id: user supplied sample, this distinguishes distinct samples of the same type from the same donor. If only one sample per donor is submitted the value may match the donor Id',
             'Sample Description: a description of the Sample Id',
             'Library Type (LT): {0}'.format(D['Library Type']),
             'Tissue Origin (TO): {0}'.format(D['Tissue Origin']),
             'Tissue Type (TT): {0}'.format(D['Tissue Type'])]         
    elif report == 'cumulative':
        L = ['OICR Case Id: OICR-generated case identifier',
             'Donor Id: user supplied donor identifier',
             'OICR Sample Id: The OICR generated sample identifier. The sample Id is formed from the following: 1. Case Id, 2. Tissue Origin, 3. Tissue Type, 4. Library Type and 5. User supplied Sample Id',
             'Sample Id: user supplied sample, this distinguishes distinct samples of the same type from the same donor. If only one sample per donor is submitted the value may match the donor Id',
             'Sample Description: a description of the Sample Id',
             'Library Type (LT): {0}'.format(D['Library Type']),
             'Tissue Origin (TO): {0}'.format(D['Tissue Origin']),
             'Tissue Type (TT): {0}'.format(D['Tissue Type'])]
       
    return L




def get_qc_metrics_table_names(library_sources, start=2):
    '''
    (list, int) -> dict
    
    Returns a dictionry with Names of QC metrics tables for each library type
    
    Parameters
    ----------
    - library_sources (list): Sorted list of library types 
    - start (int): Starting number to count the QC tables 
    '''

    # get the QC metrics sub-tables titles
    qc_subtables = {}
    for library_type in library_sources:
        qc_subtables[library_type] = 'Table {0} QC metrics for {1} libraries'.format(start, library_type)
    
    return qc_subtables


def metrics_definitions():
    '''
    (None) -> dict
    
    Returns a dictionary with definitions of metrics and identifiers
    '''
    
    definitions = {'CpG frequency': 'The frequency of CpGs within the captured regions.',
                   'Methylation {0}'.format(chr(946)): 'Proportion of the methylated signal over the total signal.',
                   'Coding (%)': 'Percentage of bases in the coding regions of the genome.',
                   'Coverage': 'Mean depth of coverage corrected for duplication rate.',
                   'On target': 'Percentage of mapped bases within the target region.',
                   'rRNA contamination': 'Percentage of reads aligning to ribosomal sequences.',
                   'Library Id': 'OICR generated library identifier.',
                   'File prefix': 'The common prefix, followed by the sequencing Read (R1, R2) and the file suffix .fastq.gz. The file prefix is formed from the following: 1. Library Id, 2. Run date, 3. Instrument Id, 4. Sequencing Instrument Run, 5. Flow cell identifier, 6. Lane number, 7. Demultiplex barcodes',
                   'Read pairs': 'Number of read pairs. The number of reads is twice the number of read pairs.',
                   '{0} methylation'.format(chr(955)): 'Ratio of methylated CpG over total CpG in the negative control',
                   'pUC19 methylation': 'Ratio of methylated CpG over total CpG in the positive control',
                   'Duplication rate': 'Percent duplication of Picard marked duplicates'
                   }
    
    return definitions


def get_metrics_appendix(library_sources):
    '''
    (list) -> dict
    
    Returns a dictionry with definitions of columns in the sample identifier table
    
    Parameters
    ----------
    - library_sources (list): Sorted list of library types 
    '''

    # get the QC metrics sub-tables and appendices
    

    definitions = metrics_definitions()    
    columns = [': '.join([i, definitions[i]]) for i in ['Library Id', 'File prefix', 'Read pairs']]
    
    counter = 1
    qc_appendix = {'tables': {}, 'metrics': {}}
    for library_type in library_sources:
        qc_appendix['tables'][library_type] = 'Appendix Table 2.{0}'.format(counter)
        if library_type == 'CM':
            qc_appendix['metrics'][library_type] = columns + [': '.join([i, definitions[i]]) for i in ['Methylation {0}'.format(chr(946)), 'CpG frequency']]
        elif library_type in ['EX', 'TS']:
            qc_appendix['metrics'][library_type] = columns + [': '.join([i, definitions[i]]) for i in ['Coverage', 'On target']]
        elif library_type in ['WG', 'PG']:
            qc_appendix['metrics'][library_type] = columns + ['{0}: {1}'.format('Coverage', definitions['Coverage'])]
        elif library_type == 'WT':
            qc_appendix['metrics'][library_type] = columns + [': '.join([i, definitions[i]]) for i in ['rRNA contamination', 'Coding (%)']]
        elif library_type in ['MC', 'MG']:
            qc_appendix['metrics'][library_type] = columns + [': '.join([i, definitions[i]]) for i in ['{0} methylation'.format(chr(955)), 'pUC19 methylation', 'Duplication rate']]
        else:
            qc_appendix['metrics'][library_type] = columns
        counter += 1
        
    return qc_appendix




def get_cumulative_metrics_appendix(library_sources):
    '''
    (list) -> dict
    
    Returns a dictionry with definitions of columns in the QC metrics table of the cumulative report
    
    Parameters
    ----------
    - library_sources (list): Sorted list of library types 
    '''

    # get the QC metrics sub-tables and appendices
    definitions = metrics_definitions()    
      
    counter = 1
    qc_appendix = {'tables': {}, 'metrics': {}}
    for library_type in library_sources:
        qc_appendix['tables'][library_type] = 'Appendix Table 2.{0}'.format(counter)
        
        L = ['OICR Case Id: OICR-generated case identifier',
             'OICR Sample Id: The OICR generated sample identifier. The sample Id is formed from the following: 1. Case Id, 2. Tissue Origin, 3. Tissue Type, 4. Library Type and 5. User supplied Sample Id',
             'Lane Count: Number of sequencing lane',
             'Total Read Pairs: {0}'.format(definitions['Read pairs'])]
        
        if library_type == 'CM':
            L.extend([': '.join([i, definitions[i]]) for i in ['Methylation {0}'.format(chr(946)), 'CpG frequency']])
        elif library_type in ['EX', 'TS']:
            L.extend([': '.join([i, definitions[i]]) for i in ['Coverage', 'On target']])
        elif library_type in ['WG', 'PG']:
            L.extend(['{0}: {1}'.format('Coverage', definitions['Coverage'])])
        elif library_type == 'WT':
            L.extend([': '.join([i, definitions[i]]) for i in ['rRNA contamination', 'Coding (%)']])
        elif library_type in ['MC', 'MG']:
            L.extend([': '.join([i, definitions[i]]) for i in ['{0} methylation'.format(chr(955)), 'pUC19 methylation', 'Duplication rate']])
        qc_appendix['metrics'][library_type] = L
        counter += 1
        
    return qc_appendix



def makepdf(html, outputfile):
    """
    (str) -> None
    
    Generates a PDF file from a string of HTML
   
    Parameters
    ----------
    - html (str) String of formated HTML
    - outputfile (str): Name of the output PDF file
    """
    
    #htmldoc = HTML(string=html, base_url=__file__)
    #htmldoc.write_pdf(outputfile, stylesheets=[CSS('./static/css/style.css')], presentational_hints=True)

    css_file = os.path.join(os.path.dirname(__file__), './static/css/style.css')
    htmldoc = HTML(string=html, base_url=__file__)
    htmldoc.write_pdf(outputfile, stylesheets=[CSS(css_file)], presentational_hints=True)


def merge_bamqc_dnaseqc(bamqc_info, dnaseqqc_info):
    '''
    (dict, dict) -> None
    
    Merge the QC information extracted from the bamqc and dnaseqqc databases.
    Note that records in dnseqqc have precedence over reords in bamqc in case of duplicates
        
    Parameters
    ----------
    - bamqc_info (dict): QC information for each paired fastq from the bamqc db
    - dnaseqqc_info (dict): QC information for each paired fastq from the dnaseqqc db
    '''
    
    D = {}
    
    # dnaseq qc has precedence over bam qc in case of duplicate records
    
    for run in bamqc_info:
        if run not in D:
            D[run] = {}
        for sample in bamqc_info[run]:
            D[run][sample] = bamqc_info[run][sample]
    
    for run in dnaseqqc_info:
        if run not in D:
            D[run] = {}
        for sample in dnaseqqc_info[run]:
            D[run][sample] = dnaseqqc_info[run][sample]
            
    return D        



def get_Y_axis_labels(library_source):
    '''
    (str) -> list
    
    Returns a list o Y axis labels for the library type
    
    Parameters
    ----------
    - library_source (str): Specific library type
    '''
    
    L = []
    if library_source == 'CM':
        L = ['Read pairs', 'Methylation {0}'.format(chr(946)), 'CpG frequency']
    elif library_source == 'WT':
        L = ['Read pairs', 'rRNA contamination', 'Coding (%)']
    elif library_source in ['WG', 'PG']:
        L = ['Read pairs', 'Coverage']
    elif library_source in ['TS', 'EX']:
        L = ['Read pairs', 'Coverage', 'On target']
    elif library_source in ['MC', 'MG']:
        L = ['Read pairs', '{0} methylation'.format(chr(955)), 'pUC19 methylation', 'Duplication rate']
    else:
        L = ['Read pairs']
    return L


def write_batch_report(args):
    '''
    (str, str, str, str, str, str, list, str, str | None, bool) -> None
    
    Write a PDF report with QC metrics and released fastqs for a given project

    - project (str): Project name as it appears in File Provenance Report
    - working-dir (str): Path to the directory with project directories and links to fastqs 
    - project_name (str): Project name used to create the project directory in gsi space
    - project_code (str): Project code from MISO
    - bamqc_db (str): Path to the bamqc db
    - cfmedipqc_db (str): Path to the cfmedipqc db 
    - run_directories (list): List of directories with links to fastqs
    - provenance (str): Path to File Provenance Report.
    - prefix (str | None): Use of prefix assumes that file paths in File Provenance Report are relative paths.
                           Prefix is added to the relative path in FPR to determine the full file path.
    - keep_html (bool): Writes html report to file if True
    '''
    
    if args.runs and args.libraries:
        sys.exit('-r and -l are exclusive parameters')
    if args.run_directories and (args.runs or args.libraries or args.release_files or args.exclude or args.nomiseq):
        sys.exit('-rn cannot be used with options -r, -l, -f, -e or --exclude_miseq')
    
    if all(map(lambda x: x is None, [args.runs, args.libraries, args.release_files, args.exclude, args.nomiseq]))  and args.run_directories is None:
        sys.exit('Please provide a list of run folders')
    # get the project directory with release run folders
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # get the records for the project of interest
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
    print('Information was extracted from FPR {0}'.format(provenance))
    # collect relevant information from File Provenance Report about fastqs for project 
    files = parse_fpr_records(provenance, args.project, ['bcl2fastq'], args.prefix)
    
    # keep only info about released fastqs
    if args.run_directories:
        released_files = list_files_release_folder(args.run_directories)
    else:
        released_files, _  = collect_files_for_release(files, args.release_files, args.nomiseq, args.runs, args.libraries, args.exclude)
        released_files = [released_files[i]['file_path'] for i in released_files]
    # resolve links
    released_files = resolve_links(released_files)
    to_remove = [file_swid for file_swid in files if os.path.realpath(files[file_swid]['file_path']) not in released_files]
    for file_swid in to_remove:
        del files[file_swid]
    
    # write sample map
    write_sample_map(args.project, files, working_dir, args.add_panel)
    print('Generated sample maps in {0}'.format(working_dir))
    print('Information was extracted from FPR {0}'.format(provenance))

    # count the number of released fastq pairs for each run and instrument
    fastq_counts = count_released_fastqs_by_instrument(files, 'read1')
    all_released_files = sum([fastq_counts[instrument][run] for instrument in fastq_counts for run in fastq_counts[instrument]])
    
    # collect information from bamqc table
    bamqc = extract_bamqc_data(args.bamqc_db)
    # collect information from dnaseqqc table
    dnaseqqc = extract_bamqc_data(args.dnaseqqc_db)
    # merge bamqc and dnaseqc
    bamqc_info = merge_bamqc_dnaseqc(bamqc, dnaseqqc)    
    
    # collect information from cfmedip table
    cfmedipqc_info = extract_cfmedipqc_data(args.cfmedipqc_db)
    # collect information from rnaseq table
    rnaseqqc_info = extract_rnaseqqc_data(args.rnaseqqc_db)
    # collect information from emseq cache
    emseqqc_info = extract_emseqqc_data(args.emseqqc_db)
    
    # update FPR info with QC metrics
    map_QC_metrics_to_fpr(files, bamqc_info, cfmedipqc_info, rnaseqqc_info, emseqqc_info)    
    
    # make a list of library types
    library_sources = sorted(list(set([files[i]['library_source'][0] for i in files])))
    
    # list all platforms for each library source
    libraries = {}
    for library_source in library_sources:
        instruments = sorted(list(set([files[file_swid]['platform'] for file_swid in files if files[file_swid]['library_source'][0] == library_source])))
        libraries[library_source] = instruments        
        
    # generate plots for each instrument and library source and keep track of figure files
    figure_files = {}
    metrics, Y_axis = {}, {}
    for library_source in libraries:
        metrics[library_source] = get_library_metrics(library_source)
        Y_axis[library_source] = get_Y_axis_labels(library_source)
        colors = ['#00CD6C', '#AF58BA', '#FFC61E', '#009ADE']
        for platform in libraries[library_source]:
            figure = generate_figures(files, args.project, library_source, platform, metrics[library_source], Y_axis[library_source], colors, working_dir)
            if library_source not in figure_files:
                figure_files[library_source] = {}
            figure_files[library_source][platform] = figure
             
    # count the number of samples with missing metric values
    samples_missing_metrics = count_samples_with_missing_values(files)
    
    # issue warning if samples with missing QC metrics
    if samples_missing_metrics:
        print('========')
        print('WARNING!!')
        print('Some samples have missing QC information. Please review')
        for i in samples_missing_metrics:
            for j in samples_missing_metrics[i]:
                print('Library type: {0} - Platform: {1} - {2} Samples'.format(i, j, samples_missing_metrics[i][j]))
        print('========')        
    
    # write md5sums to separate file
    current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    md5sum_file = os.path.join(working_dir, '{0}.batch.release.{1}.md5'.format(args.project, current_time))
    write_md5sum(files, md5sum_file)

    # write report
    # get the report template
    #environment = Environment(loader=FileSystemLoader("./templates/"))
    #template = environment.get_template("batch_report_template.html")
    
    template_dir = os.path.join(os.path.dirname(__file__), './templates')
    environment = Environment(loader = FileSystemLoader(template_dir), autoescape = True)
    template = environment.get_template("batch_report_template.html")
    
    # make a dict with project information
    projects = {'acronym': args.project, 'name': args.project_full_name, 'date': time.strftime('%Y-%m-%d', time.localtime(time.time()))}

    # group metrics by pairs of files
    header_identifiers = ['Library Id', 'Case Id', 'Donor Id', 'Sample Id', 'Sample Description', 'LT', 'TO', 'TT']
    
    sample_identifiers = group_sample_metrics(files, 'sample_identifiers')
    appendix_identifiers = get_identifiers_appendix(files, 'batch')
    
    qc_metrics = group_sample_metrics(files, 'qc_metrics', metrics)
    header_metrics = {}
    for i in library_sources:
        header_metrics[i] = ['Library Id', 'File prefix'] + Y_axis[i]
    
    
    # get the qc metrics subtables
    qc_subtables = get_qc_metrics_table_names(library_sources)
    # get the metrics appendix
    qc_appendices = get_metrics_appendix(library_sources)
    
    # fill in template
    context = {'projects' : projects,
               'file_count': all_released_files,
               'fastq_counts': fastq_counts,
               'figure_files': figure_files,
               'samples_missing_metrics': samples_missing_metrics,
               'header_identifiers': header_identifiers,
               'sample_identifiers': sample_identifiers,
               'appendix_identifiers': appendix_identifiers,
               'header_metrics': header_metrics,
               'qc_metrics': qc_metrics,
               'qc_subtables': qc_subtables,
               'qc_appendices': qc_appendices,
               'library_sources': library_sources,
               'libraries': libraries, 
               'user': args.user,
               'ticket': os.path.basename(args.ticket),
               'md5sum': os.path.basename(md5sum_file)}
       
    # render template html 
    content = template.render(context)

    # save html file to disk
    if args.keep_html:
        html_file = os.path.join(working_dir, '{0}_run_level_data_release_report.{1}.html'.format(args.project, current_time))
        newfile = open(html_file, 'w')
        newfile.write(content)
        newfile.close()

    # convert html to PDF
    report_file = os.path.join(working_dir,  '{0}_run_level_data_release_report.{1}.pdf'.format(args.project, current_time))
    makepdf(content, report_file)

    # remove figure files from disk    
    if args.keep_html == False:
        for i in figure_files:
            for j in figure_files[i]:
                if os.path.isfile(figure_files[i][j]):
                    os.remove(figure_files[i][j])



def extract_merged_bamqc_data(merged_bamqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the bamqc table of qc-etl
        
    Parameters
    ----------
    - bamqc_db (str): Path to the bamqc SQLite database generated by qc-etl
    '''
            
    conn = sqlite3.connect(merged_bamqc_db)
    cur = conn.cursor()
    
    # get table name
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    table_name = [i[0] for i in cur][0]
    # get all data from table
    conn.row_factory = sqlite3.Row
    data = conn.execute('select * from {0}'.format(table_name)).fetchall()
        
    columns = ['library', 'Library Design', 'Donor', 'Project', 'Tissue Origin', 
               'Tissue Type', 'Group ID', 'sample', 'Merged Pinery Lims ID', 'instrument',
               'sample', 'bases mapped', 'coverage', 'coverage deduplicated',
               'mapped reads', 'mark duplicates_PERCENT_DUPLICATION', 'total bases on target',
               'total reads']
    
    D = {}
    
    for i in data:
        d = {}
        for j in columns:
            if j in ['coverage', 'coverage deduplicated', 'mark duplicates_PERCENT_DUPLICATION']:
                d[j] = float(i[j])
            elif j in ['bases mapped', 'mapped reads', 'total bases on target', 'total reads']:
                d[j] = int(i[j])
            elif j == 'library':
                d[j] = i[j].replace('\"', '').split(',')
            elif j == 'instrument':
                d[j] = i[j].replace('\"', '')
            elif j == 'Merged Pinery Lims ID':
                d[j] = list(map(lambda x: x.strip(), i[j].replace('[', '').replace(']', '').replace('\"', '').split(',')))
            else:
                d[j] = i[j]
        # compute on_target rate, not available through qc-etl
        d['on_target'] = compute_on_target_rate(d['bases mapped'], d['total bases on target']) 
        
        merged_samples = i['Merged Pinery Lims ID']
        merged_samples = list(map(lambda x: x.strip(), merged_samples.replace('[', '').replace(']', '').replace('\"', '').split(',')))
        
        merged_samples = ';'.join(merged_samples)
        D[merged_samples] = d
        
    return D


def extract_merged_rnaseqqc_data(merged_rnaseqqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the rnaseqqc table of qc-etl
        
    Parameters
    ----------
    - rnaseqqc_db (str): Path to the rnaseq SQLIte database generated by qc-etl
    '''
    
    #conn = sqlite3.connect('prov_report_test.db')
    conn = sqlite3.connect(merged_rnaseqqc_db)
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    table_name = [i[0] for i in cur][0]
    conn.row_factory = sqlite3.Row
    
    # get all data
    data = conn.execute('select * from {0}'.format(table_name)).fetchall()
    conn.close()

    D = {}

    columns = ['Donor',
               'Project',
               'Merged Pinery Lims ID',
               'Group ID',
               'Library Design',       
               'Tissue Origin',       
               'Tissue Type',       
               'PCT_CODING_BASES', 
               'rrna contamination properly paired',
               'rrna contamination in total (QC-passed reads + QC-failed reads)',
               'PCT_CORRECT_STRAND_READS',
               'MEDIAN_5PRIME_TO_3PRIME_BIAS']
        
    
    for i in data:
        d = {}
        for j in columns:
            assert j in dict(i).keys()
            if j == 'Merged Pinery Lims ID':
                d[j] = list(map(lambda x: x.strip(), i[j].replace('[', '').replace(']', '').replace('\"', '').split(',')))
            else:
                d[j] = i[j]
        
        merged_samples = i['Merged Pinery Lims ID']
        merged_samples = list(map(lambda x: x.strip(), merged_samples.replace('[', '').replace(']', '').replace('\"', '').split(',')))
        merged_samples = ';'.join(merged_samples)
        
        D[merged_samples] = d
        
    return D
    

def group_fastq_pairs(files):
    '''
    (dict) -> dict
    
    Returns a dictionary with lists of fastq pairs for each instrument and library type
    Pre-condition: fastqs are paired with R1 and R2
    
    Parameters
    ----------
    - files (dict): Dictionary with fastq information extracted from FPR
    '''
  
    # find pairs of files for each library type and instrument
    pairs = {}
    instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
    for platform in instruments:
        P = find_fastq_pairs(files, platform)
        library_sources = list(set([i[0]['library_source'][0] for i in P]))
        pairs[platform] = {}
        for j in library_sources:
            L = [i for i in P if i[0]['library_source'][0] == j]
            add_file_prefix(L)
            pairs[platform][j] = L
    return pairs


def collect_bamqc_metrics(D, pairs, bamqc, platform, library_source):
    '''
    (dict, dict, dict, str, str) -> None
    
    Update dictionary D with QC metrics from bamqc for fastq pairs sorted by instrument and library type
    
    Parameters
    ----------
    - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
    - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
    - bamqc (dict): Dictionary with QC metrics from the merged bamqc cache 
    - platform (str): Sequencing platform
    - library_source (str): 2-letters code of library type
    '''
    
    for i in bamqc:
        merged_limskeys = sorted(bamqc[i]['Merged Pinery Lims ID'])            
        # collect limkskeys from FPR records
        L = []
        readcount = []
        prefix = []
        runs = []
        lanes = []
        
        # collect read count for all fastqs
        # track limskeys and file prefix for fastqs
        for j in pairs[platform][library_source]:
            for k in j:
                limskey = k['limskey'][0]
                if library_source == bamqc[i]['Library Design'] \
                and platform == '_'.join(bamqc[i]['instrument'].split()) \
                and limskey in merged_limskeys \
                and k['groupid'][0] == bamqc[i]['Group ID']:
                    assert k['sample_id'][0] == bamqc[i]['sample']
                    # get the number of read for each file
                    readcount.append(k['read_count'])
                    L.append(limskey)
                    prefix.append(k['prefix'])
                    runs.extend(k['run_id'])
                    lanes.extend(k['run'])
        # check that all limskeys have been recorded
        if sorted(list(set(L))) == sorted(merged_limskeys):
            assert sum(readcount) % 2 == 0
            readcount = sum(readcount)
            prefix = list(set(prefix))
            runs = sorted(list(set(runs)))
            lanes = list(set(lanes))
            libraries = bamqc[i]['library']
            donor = bamqc[i]['Donor']
            coverage = round(bamqc[i]['coverage'], 2)
            coverage_dedup = round(bamqc[i]['coverage deduplicated'], 2)
            on_target = round(bamqc[i]['on_target'], 2)
            percent_duplicate = round(bamqc[i]['mark duplicates_PERCENT_DUPLICATION'], 2)
            sample = bamqc[i]['sample']

            d = {'read_count': readcount,
                 'prefix': prefix,
                 'run': runs,
                 'lane': len(lanes),
                 'libraries': libraries,
                 'donor': donor,
                 'coverage': coverage,
                 'coverage_dedup': coverage_dedup,
                 'on_target': on_target,
                 'percent_duplicate': percent_duplicate,
                 'sample': sample}
                                   
            if platform not in D:
                D[platform] = {}
            if library_source not in D[platform]:
                D[platform][library_source] = {}
            assert i not in D[platform][library_source]
            D[platform][library_source][i] = d
    

def collect_rnaseqc_metrics(D, pairs, rnaseqqc, platform, library_source):
    '''
    (dict, dict, dict, str, str) -> None
    
    Update dictionary D with QC metrics from rnaseqc for fastq pairs sorted by instrument and library type
    
    Parameters
    ----------
    - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
    - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
    - rnaseqqc (dict): Dictionary with QC metrics from the merged rnaseqqc cache 
    - platform (str): Sequencing platform
    - library_source (str): 2-letters code of library type
    '''
    
    for i in rnaseqqc:
        merged_limskeys = sorted(rnaseqqc[i]['Merged Pinery Lims ID'])            
        # collect limkskeys from FPR records
        L = []
        readcount = []
        prefix = []
        # track library ids because not present in rnaseq merged cache
        libraries = []
        runs = []
        lanes = []
        
        # collect read count for all fastqs
        # track limskeys and file prefix for fastqs
        for j in pairs[platform][library_source]:
            for k in j:
                limskey = k['limskey'][0]
                if library_source == rnaseqqc[i]['Library Design'] \
                and limskey in merged_limskeys \
                and k['groupid'][0] == rnaseqqc[i]['Group ID']:
                    # get the number of read for each file
                    readcount.append(k['read_count'])
                    L.append(limskey)
                    prefix.append(k['prefix'])
                    libraries.extend(k['library'])
                    runs.extend(k['run_id'])
                    lanes.extend(k['run'])
        # check that all limskeys have been recorded
        if sorted(list(set(L))) == sorted(merged_limskeys):
            assert sum(readcount) % 2 == 0
            readcount = sum(readcount)
            prefix = list(set(prefix))
            runs = sorted(list(set(runs)))
            lanes = list(set(lanes))
            libraries = list(set(libraries))
            donor = rnaseqqc[i]['Donor']
            bias = rnaseqqc[i]['MEDIAN_5PRIME_TO_3PRIME_BIAS']
            contamination = round((rnaseqqc[i]['rrna contamination properly paired'] / rnaseqqc[i]['rrna contamination in total (QC-passed reads + QC-failed reads)'] * 100), 3)
            coding = round(rnaseqqc[i]['PCT_CODING_BASES'], 3)
            strand = round(rnaseqqc[i]['PCT_CORRECT_STRAND_READS'], 3)
            sample = '_'.join([donor, rnaseqqc[i]['Tissue Origin'], rnaseqqc[i]['Tissue Type'], rnaseqqc[i]['Library Design'], rnaseqqc[i]['Group ID']])

            d = {'read_count': readcount,
                 'prefix': prefix,
                 'run': runs,
                 'lane': len(lanes),
                 'libraries': libraries,
                 'donor': donor,
                 "5'-3' bias": bias,
                 'rRNA contamination': contamination,
                 'Coding (%)': coding,
                 'Correct strand reads (%)': strand,
                 'sample': sample}
                     
            if platform not in D:
                D[platform] = {}
            if library_source not in D[platform]:
                D[platform][library_source] = {}
            assert i not in D[platform][library_source]
            D[platform][library_source][i] = d


def collect_cfmedipqc_metrics(D, pairs, cfmedipqc, platform, library_source):
    '''
    (dict, dict, dict, str, str) -> None
    
    Update dictionary D with QC metrics from cfmedipqc for fastq pairs sorted by instrument and library type
    
    Parameters
    ----------
    - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
    - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
    - cfmedipqc (dict): Dictionary with QC metrics from the cfmedipqc cache 
    - platform (str): Sequencing platform
    - library_source (str): 2-letters code of library type
    '''
    
    for i in pairs[platform][library_source]:
        for j in i:
            limskey = j['limskey'][0]
            run = j['run_id'][0]
            prefix = j['prefix']
            barcode = j['barcode'][0]
            lane = j['lane'][0]
            library = j['library'][0]
            donor = j['sample_name']
            # get the read count for each file
            readcount = j['read_count']
            sample = j['sample_id'][0]
            lanes = j['run']
                        
            if run in cfmedipqc and limskey in cfmedipqc[run]:
                assert len(cfmedipqc[run][limskey]) == 1
                d = cfmedipqc[run][limskey][0]
                assert d['Pinery Lims ID'] == limskey
                assert d['Run Alias'] == run
                assert int(d['Lane Number']) == int(lane)
                assert d['Barcodes'] == barcode
                AT_dropout = round(d['AT Dropout'], 3)
                methylation_beta = round(d['Methylation beta'], 3)
                duplication = round(d['Percent Duplication'], 3)
                CpG_enrichment = round(d['Relative CpG Frequency in Regions'], 3)
                    
                if platform not in D:
                    D[platform] = {}
                if library_source not in D[platform]:
                    D[platform][library_source] = {}
                if limskey in D[platform][library_source]:
                    assert D[platform][library_source][limskey]['prefix'][0] == prefix
                    D[platform][library_source][limskey]['read_count'] += readcount
                    assert D[platform][library_source][limskey]['read_count'] % 2 == 0
                    assert D[platform][library_source][limskey]['sample'] == sample
                else:
                    D[platform][library_source][limskey] = {
                        'read_count': readcount,
                        'prefix': [prefix],
                        'run': run,
                        'lane': len(lanes),
                        'libraries': [library],
                        'donor': donor,
                        'AT_dropout': AT_dropout,
                        'methylation_beta': methylation_beta,
                        'duplication': duplication,
                        'CpG_enrichment': CpG_enrichment,
                        'sample': sample
                        }
                            

def collect_emseqqc_metrics(D, pairs, emseqqc, platform, library_source):
    '''
    (dict, dict, dict, str, str) -> None
    
    Update dictionary D with QC metrics from emseqqc for fastq pairs sorted by instrument and library type
    
    Parameters
    ----------
    - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
    - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
    - emseqqc (dict): Dictionary with QC metrics from the emseqqc cache 
    - platform (str): Sequencing platform
    - library_source (str): 2-letters code of library type
    '''
    
    
    for i in pairs[platform][library_source]:
        for j in i:
            limskey = j['limskey'][0]
            run = j['run_id'][0]
            prefix = j['prefix']
            barcode = j['barcode'][0]
            lane = j['lane'][0]
            lanes = j['run']
            library = j['library'][0]
            donor = j['sample_name']
            sample = j['sample_id'][0]
            # get the number of reads for each file
            readcount = j['read_count']
            if run in emseqqc and limskey in emseqqc[run]:
                assert len(emseqqc[run][limskey]) == 1
                d = emseqqc[run][limskey][0]
                assert d['Pinery Lims ID'] == limskey
                assert d['Run Alias'] == run
                assert int(d['Lane Number']) == int(lane)
                assert d['Barcodes'] == barcode
                Lambda_methylation = round(d['Lambda'], 3)
                pUC19_methylation = round(d['pUC19'], 3)
                percent_duplication = round(d['mark duplicates_PERCENT_DUPLICATION'], 3)

                if platform not in D:
                    D[platform] = {}
                if library_source not in D[platform]:
                    D[platform][library_source] = {}
                if limskey in D[platform][library_source]:
                    assert D[platform][library_source][limskey]['prefix'][0] == prefix
                    D[platform][library_source][limskey]['read_count'] += readcount
                    assert D[platform][library_source][limskey]['read_count'] % 2 == 0
                    assert D[platform][library_source][limskey]['sample'] == sample
                else:
                    D[platform][library_source][limskey] = {
                        'read_count': readcount,
                        'prefix': [prefix],
                        'run': run,
                        'lane': len(lanes),
                        'libraries': [library],
                        'donor': donor,
                        'Lambda_methylation': Lambda_methylation,
                        'pUC19_methylation': pUC19_methylation,
                        'percent_duplication': percent_duplication,
                        'sample': sample
                        }
                
                
def convert_read_count_to_read_pairs(D):
    '''
    (dict) -> None
    
    Convert the read counts to read pair counts in place in dictionary D
    
    Parameters
    ----------
    - D (dict): Dictionary with metrics of interest organized by sequencing plarform and library source
    '''
    
    for platform in D:
        for library_source in D[platform]:
            for limskeys in D[platform][library_source]:
                # get the number of read pairs
                assert D[platform][library_source][limskeys]['read_count'] % 2 == 0
                D[platform][library_source][limskeys]['read_count'] = int(D[platform][library_source][limskeys]['read_count'] / 2)  

    
def get_metrics_cumulative_report(files, bamqc, rnaseqqc, cfmedipqc, emseqqc):
    '''
    (dict, dict, dict, dict, dict) -> dict
    
    Returns a dictionary with identifiers and QC metrics extracted from the QC-etl caches
    
    Parameters
    ----------
    - files (dict): Fastq information extracted from FPR
    - bamqc (dict): Dictionary with bamqc QC metrics
    - rnaseqqc (dict): Dictionary with rnaseqqc QC metrics
    - cfmedipqc (dict): Dictionary with cfmedipqc QC metrics
    - emseqqc (dict): Dictionary with emseqqc QC metrics
    '''
    
    # organize metrics for cumulative report
    
    D = {}
    
    # find pairs of files for each library type and instrument
    pairs = group_fastq_pairs(files)
        
    # Group metrics by platform and library type
    for platform in pairs:
        for library_source in pairs[platform]:
            if library_source not in ['CM', 'WT', 'MC', 'MG']:
                # use the bamqc merged cache
                collect_bamqc_metrics(D, pairs, bamqc, platform, library_source)
            elif library_source == 'WT':
                # use the rnaseqc merged cache
                collect_rnaseqc_metrics(D, pairs, rnaseqqc, platform, library_source)
            elif library_source in ['MC', 'MG']:
                # use the lane level emseqqc cache
                collect_emseqqc_metrics(D, pairs, emseqqc, platform, library_source)
            elif library_source == 'CM':
                # use the lane level cfmedip cache
                collect_cfmedipqc_metrics(D, pairs, cfmedipqc, platform, library_source)
    
    return D                

def format_qc_metrics(qc_metrics):
    '''
    (dict) -> dict
    
    Returns a dictionary with relevant QC metrics and identifiers
    for each instrument and library source
    
    Parameters
    ----------
    - qc_metrics (dict): Dictionary with QC metrics collected from the QC-etl caches
    '''
    
    D = {}
    for instrument in qc_metrics:
        for library_source in qc_metrics[instrument]:
            metrics = get_library_metrics(library_source)
            if instrument not in D:
                D[instrument] = {}
            if library_source not in D[instrument]:
                D[instrument][library_source] = []
            # collect relevant information
            for limskey in qc_metrics[instrument][library_source]:
                donor = qc_metrics[instrument][library_source][limskey]['donor']
                libraries = list(qc_metrics[instrument][library_source][limskey]['libraries'])
                prefix = list(qc_metrics[instrument][library_source][limskey]['prefix'])
                run = list(qc_metrics[instrument][library_source][limskey]['run'])
                libraries = [fit_into_column(i) if len(i) >=40 else i for i in libraries]
                prefix = [fit_into_column(i) if len(i) >=40 else i for i in prefix]
                run = [fit_into_column(i) if len(i) >=40 else i for i in run]
                libraries.sort()
                prefix.sort()
                run.sort()
                sample = qc_metrics[instrument][library_source][limskey]['sample']
                lane_count = qc_metrics[instrument][library_source][limskey]['lane']
                
                L = [donor, sample, lane_count]
                
                for i in metrics:
                    if i == 'read_count':
                        L.append('{:,}'.format(qc_metrics[instrument][library_source][limskey][i]))                     
                    else:
                        L.append(qc_metrics[instrument][library_source][limskey][i])
                D[instrument][library_source].append(L)         
                # sort according to donor
                D[instrument][library_source].sort(key=lambda x: x[0])

    return D


def format_sequencing_table(files):
    '''
    (dict) -> list
    
    Returns a list with sequencing information
    
    Parameters
    ----------
    - files (dict): Dictionary with file information extracted from FPR
    '''

    pairs = group_fastq_pairs(files)

    SL = []
    
    for instrument in pairs:
        for library_source in pairs[instrument]:
            for i in pairs[instrument][library_source]:
                donor = i[0]['sample_name']
                sample = i[0]['sample_id'][0]
                library = i[0]['library'][0]
                library_type = i[0]['library_source'][0]
                prefix = i[0]['prefix']
                read_pairs = i[0]['read_count']
                L = [donor, sample, library, prefix, read_pairs]
                if L not in SL:
                    SL.append(L)
                # sort according to donor, sample, snf library
                SL.sort(key = lambda x: (x[0], x[1], x[2]))
    # adjust long names to fit in column
    for i in range(len(SL)):
        L = [fit_into_column(j) if len(j) >=30 else j for j in SL[i][:-1]]
        L.append('{:,}'.format(SL[i][-1]))
        SL[i] = L
    
    return SL
   

def generate_cumulative_figures(working_dir, project, qc_metrics, platform, library_type, Y_axis, colors, width=13, height=16):
    '''
    (str, str, dict, str, str, list, list, int, int) -> str
        
    Generate a figure with metrics from FPR and QC-etl for a specific sequencing platform
    and library type and returns the path to the figure file
        
    Parameters
    ----------
    - working_dir (str): Path to the folder where figure files are written
    - project (str): Name of project
    - qc_metrics (dict): Dictionary with QC metrics extracted from QC etl caches
    - platform (str): Sequencing platform
    - library_type (str): Type of library
    - Y_axis (list): List of Y_axis labels   
    - colors (list): List of colors
    - width (int): Width of the figure
    - height (int): Height of the figure
    '''
    
    metrics_interest = get_library_metrics(library_type)

    # make lists of metrics, sort according to read count
    D = {}
    for i in qc_metrics[platform][library_type]:
        for j in metrics_interest:
            if j not in D:
                D[j] = [qc_metrics[platform][library_type][i][j]]
            else:
                D[j].append(qc_metrics[platform][library_type][i][j])
    
    L = [D[i] for i in metrics_interest]
    L = sort_metrics(L)
    
    # get the outputfile
    current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    outputfile = os.path.join(working_dir, '{0}.{1}.{2}.{3}.QC_plots.png'.format(project, platform, library_type, current_time))
    
    if L[0]:
        figure = plt.figure()
        figure.set_size_inches(width, height)
        # make a list of with X axis labels to determine which subplot should display the Samples label
        x_labels = get_x_axis_labels(L)
        # determine how many subplots are expected
        subplots = count_subplots(L) 
        # determine the position of each subplot
        subplot_pos = get_subplot_position(L)
        
        for i in range(len(L)):
            # determine title
            title = platform + ' {0} libraries'.format(library_type) if i == 0 else None
            # plot data
            create_ax(subplots, 1, subplot_pos[i], figure, L[i], Y_axis[i], colors[i], title = title, XLabel = x_labels[i])
                
        # make sure axes do not overlap
        plt.tight_layout(pad = 2.5)
        # write figure to file  
        figure.savefig(outputfile, bbox_inches = 'tight')
        plt.close()
                
        return outputfile
    else:
        return ''

    

def write_cumulative_report(args):
    '''
    (str, str, str, str, str | None, str, str, str, str, str, str, str) -> None
        
    Write a cumulative PDF report with QC metrics for all released fastqs for a given project
  
    - project (str): Name of project of interest
    - projects_dir (str): Parent directory containing the project subdirectories with file links
    - project_name (str): Project name used to create the project directory in gsi space
    - provenance (str): Path to File Provenance Report
    - prefix (str | None): Use of prefix assumes that FPR containes relative paths.
                           Prefix is added to the relative paths in FPR to determine the full file paths
    - api (str): URL of the Nabu API
    - merged_bamqc_db (str): Path to the merged bamqc SQLite database
    - merged_rnaseqqc_db (str): Path to the merged rnaseq SQLite database
    - cfmedipqc_db (str): Path to the cfmedip SQLite database
    - emseqqc_db (str): Path to the emseq SQLite database
    - project_full_name (str): Full name of the project
    - ticket (str): Jira data release ticket code
    - user (str): Name of the GSI personnel generating the report
    '''
    
    # get current time
    current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
        
    # get the project directory with release run folders
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # get the records for the project of interest
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
    print('Information was extracted from FPR {0}'.format(provenance))
    
    # make a dict with project information
    projects = {'acronym': args.project, 'name': args.project_full_name, 'date': time.strftime('%Y-%m-%d', time.localtime(time.time()))}

    # get information about the released fastqs
    # collect relevant information from File Provenance Report about fastqs for project 
    files = parse_fpr_records(provenance, args.project, ['bcl2fastq'], args.prefix)
    
    # get the released files at the project level from nabu
    released_files = list_released_fastqs_project(args.api, args.project)
    
    # resolve links
    released_files = resolve_links(released_files)
    # remove files not released
    to_remove = [file_swid for file_swid in files if os.path.realpath(files[file_swid]['file_path']) not in released_files]
    for file_swid in to_remove:
        del files[file_swid]
       
    # count the number of released fastq pairs for each run and instrument
    fastq_counts = count_released_fastqs_by_library_type_instrument(files)
    all_released_files = 0 
    for i in fastq_counts:
        all_released_files += sum([sum(list(j.values())) for j in fastq_counts[i].values()])
    
    # sort library types, sequencing platforms, sequencing runs
    lt = sorted(list(fastq_counts.keys()))
    sp = sp = {i:sorted(fastq_counts[i].keys()) for i in lt}
    rn = {}
    for i in lt:
        rn[i] = {}
        for j in sp[i]:
            rn[i][j] = rn[i][j] = sorted(list(fastq_counts[i][j].keys()))
        
    # get the identifiers of all released files
    sample_identifiers = group_cumulative_samples(files)
    
    appendix_identifiers = get_identifiers_appendix(files, 'cumulative')
    # define identifier table header
    header_identifiers = ['OICR Case Id', 'Donor Id', 'OICR Sample Id', 'Sample Description', 'LT', 'TO', 'TT']
    
    # get information for lane level sequencing
    sequencing = format_sequencing_table(files)
        
    # collect information from merged bamqc table
    bamqc_info = extract_merged_bamqc_data(args.merged_bamqc_db)
    # collect information from rnaseq table
    rnaseqqc_info = extract_merged_rnaseqqc_data(args.merged_rnaseqqc_db)
    # collect information from cfmedip table
    cfmedipqc_info = extract_cfmedipqc_data(args.cfmedipqc_db)
    # collect information from emseq cache
    emseqqc_info = extract_emseqqc_data(args.emseqqc_db)
    print('Extracted QC metrics from QC-etl caches')
    
    # group identifiers, metrics by instrument and library type
    library_metrics = get_metrics_cumulative_report(files, bamqc_info, rnaseqqc_info, cfmedipqc_info, emseqqc_info)
    # convert read counts to counts of read pairs
    convert_read_count_to_read_pairs(library_metrics)
    # format qc metrics for template
    qc_metrics = format_qc_metrics(library_metrics)
    
    # list all platforms for each library source
    platforms = {}
    for i in qc_metrics:
        for j in qc_metrics[i]:
            if j in platforms:
                platforms[j].append(i)
                platforms[j].sort()
            else:
                platforms[j] = [i]

    # make a list of library types from data with QC metrics
    library_sources = sorted(list(platforms.keys()))
       
    # get the qc metrics subtables
    qc_subtables = get_qc_metrics_table_names(library_sources, 3)
    # get the metrics appendix
    qc_appendices = get_cumulative_metrics_appendix(library_sources)
    
      
    # generate plots for each instrument and library source and keep track of figure files
    figure_files = {}
    Y_axis = {}
    colors = ['#00CD6C', '#AF58BA', '#FFC61E', '#009ADE']
    for library_source in platforms:
        Y_axis[library_source] = get_Y_axis_labels(library_source)
        for i in range(len(Y_axis[library_source])):
            if Y_axis[library_source][i] == 'Read pairs':
                Y_axis[library_source][i] = 'Total Read Pairs'
            elif Y_axis[library_source][i] == 'Coverage':
                Y_axis[library_source][i] = 'Final Merged Coverage'
        for instrument in platforms[library_source]:
            figure = generate_cumulative_figures(working_dir, args.project, library_metrics, instrument, library_source, Y_axis[library_source], colors)
            if library_source not in figure_files:
                figure_files[library_source] = {}
            figure_files[library_source][instrument] = figure
    print('Generated QC plots')
        
    header_metrics = {}
    for i in library_sources:
        header_metrics[i] = ['OICR Case Id', 'OICR Sample Id', 'Lane Count'] + Y_axis[i]
         
    # write report
    # get the report template
    template_dir = os.path.join(os.path.dirname(__file__), './templates')
    environment = Environment(loader = FileSystemLoader(template_dir), autoescape = True)
    template = environment.get_template("cumulative_report_template.html")
       
    
    # fill in template
    context = {'projects' : projects,
               'file_count': all_released_files,
               'fastq_counts': fastq_counts,
               'header_identifiers': header_identifiers,
               'sample_identifiers': sample_identifiers,
               'appendix_identifiers': appendix_identifiers,
               'sequencing': sequencing,
               'user': args.user,
               'ticket': os.path.basename(args.ticket),
               'library_sources': library_sources,
               'qc_subtables': qc_subtables,
               'qc_appendices': qc_appendices,
               'header_metrics': header_metrics,
               'qc_metrics': qc_metrics,
               'figure_files': figure_files,
               'platforms': platforms,
               'lt': lt,
               'sp': sp,
               'rn': rn
               }
    
    # render template html 
    content = template.render(context)

    # convert html to PDF
    report_file = os.path.join(working_dir,  '{0}_cumulative_data_release_report.{1}.pdf'.format(args.project, current_time))
    makepdf(content, report_file)

    # remove figure files
    for i in figure_files:
        for j in figure_files[i]:
            if os.path.isfile(figure_files[i][j]):
                os.remove(figure_files[i][j])



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

    for j in L:
        swid = j[44]
        # get file path and resolve any links
        file = os.path.realpath(j[46])
        D[file] = swid         
    return D



def change_nabu_status(api, file_swids, qc_status, user_name, comment=None):
    '''
    (str, list, str, str, str | None) -> None
    
    Modify the record of file in file_swids in Nabu with QC status updated to qc_status,
    username updated to user_name and comment updated to comment if used 
    
    Parameters
    ----------
    - api (str): URL of the nabu API
    - file_swids (list): List of file unique identifiers
    - qc_status (str): File QC status: PASS, PENDING or FAIL
    - comment (str): Jira ticket of the release
    '''
    
    # get end-point
    api += 'add-fileqcs' if api[-1] == '/' else '/add-fileqcs'
        
    if qc_status not in ['PASS', 'PENDING', 'FAIL']:
        raise ValueError('QC status is PASS, FAIL or PENDING')
    
    headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
    json_data = {'fileqcs': []}
    for file_swid in file_swids:
        d = {'fileid': file_swid, 'qcstatus': qc_status, 'username': user_name}
        if comment:
            d['comment'] = comment    
        json_data['fileqcs'].append(d)
         
    response = requests.post(api, headers=headers, json=json_data)
    
    # check response code
    if response.status_code == 201:
        for file_swid in file_swids:
            # record created
            print('Successfully updated {0} status to {1}'.format(file_swid, qc_status))
    else:
        for file_swid in file_swids:
            print('Could not update {0} status. Nabu response code: {1}'.format(file_swid, response.status_code))





def get_pipeline_swids(data_structure, files):
    '''
    (dict, dict) -> list   
    
    Returns a list of file swids for each sample and workflow specified in the data_structure
    
    Parameters
    ----------
    - data_structure (str): Dictionary with samples, workflows and workflow_run_id hierarchical structure
    - files (dict): Dictionary with all file records for a given project extracted from FPR
    '''
    
    L = []
         
    for donor in data_structure:
        for sample_id in data_structure[donor]:
            for file_swid in files:
                sample = files[file_swid]['sample_name']
                workflow = files[file_swid]['workflow']
                version = files[file_swid]['workflow_version']
                wf_id = files[file_swid]['workflow_run_id']
                if donor == sample and sample in sample_id and workflow in data_structure[donor][sample_id]:
                    for d in data_structure[donor][sample_id][workflow]:
                        if version == d['workflow_version'] and wf_id == d['workflow_id']:
                            # check which files are collected
                            if 'extension' in d:
                                # files are collected based on file extension
                                # get the file extension
                                extension = pathlib.Path(files[file_swid]['file_path']).suffix
                                if extension in d['extension']:
                                    L.append(file_swid)
                            elif 'files' in d:
                                # files are collected based on file name
                                if os.path.basename(files[file_swid]['file_path']) in list(map(lambda x: os.path.basename(x),  d['files'])):
                                    L.append(file_swid)
                            elif 'rename_files' in d:
                                # files are collected based on file name and renamed
                                file_paths = []
                                for i in d['rename_files']:
                                    file_paths.append(os.path.basename(i['file_path']))
                                if os.path.basename(files[file_swid]['file_path']) in file_paths:
                                    L.append(file_swid)
                            else:
                                L.append(file_swid)
    
    L = list(set(L))
    
    return L                        
                        


def mark_files_nabu(args):
    '''
    (str, str, str) -> None

    Mark released files with user name and PASS and withheld files with user name and FAIL in Nabu

    Parameters
    ----------    
    - project (str): Name of Project in FPR
    - directory (str): Directory with links organized by project and run in gsi space 
    - status (str): Mark files accordingly when released or withheld. Valid options:
                     - pass: files that are released
                     - fail: files that are withheld
    - user (str): User name to appear in Nabu for each released or whitheld file
    - comment (str): A comment to used to tag the file. For instance the Jira ticket 
    - provenance (str): Path to File Provenance Report
    - run_directory (str): Directory with links organized by project and run in gsi space
    - runs (list): List of run IDs
    - libraries (str): Path to file with libraries tagged for release
    - workflow (str): Worflow used to generate the output files
    - nomiseq (bool): Exclude MiSeq runs if True
    - prefix (str): Use of prefix assumes that FPR containes relative paths.
                    Prefix is added to the relative paths in FPR to determine the full file paths
    - exclude (str): Path to file with libraries tagged for non-release
    - release_files (str): File with file names to be released
    - nabu (str):  URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca
    - provenance (str): Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz
    - analysis (str): Path to the file with hierarchical structure storing sample and workflow ids
    '''
    
    # check valid combinations of parameters
    if args.runs and args.libraries:
        sys.exit('-r and -l are exclusive parameters')
    if args.run_directory and (args.workflow or args.runs or args.libraries or args.release_files or args.exclude or args.prefix or args.nomiseq):
        sys.exit('-rn cannot be used with options -w, -r, -l, -f, -e, -px or --exclude_miseq')
    if args.analysis is None:
        if all(map(lambda x: x is None, [args.workflow, args.runs, args.libraries, args.release_files, args.exclude, args.prefix, args.nomiseq]))  and args.run_directory is None:
            sys.exit('Please provide the path to a run folders')
        if args.workflow is None and args.run_directory is None:
            sys.exit('Please provide the path to a run folder or a workflow')
    else:
        if args.runs or args.libraries or args.run_directory or args.workflow or args.release_files or args.exclude or args.prefix or args.nomiseq:
            sys.exit('''--analysis cannot be used with options -r, -l, -rn, -w, -f, -e, -px or --exclude_miseq
                     You are attempting to mark analysis pipeline files. Provide the same json structure used to link the files''')
       
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
           
    if args.run_directory:
        # get the swids of the files linked in run directory
        # check directory
        if os.path.isdir(args.run_directory) == False:
            sys.exit('{0} is not a valid directory'.format(args.run_directory))
        # make a list of files in directory
        # get the real path of the links in directory
        files = [os.path.realpath(os.path.join(args.run_directory, i)) for i in os.listdir(args.run_directory) if os.path.isfile(os.path.join(args.run_directory, i))]
        # make a list of swids
        swids = []
        records = get_FPR_records(args.project, provenance)
        print('Information was extracted from FPR {0}'.format(provenance))
        mapped_files = map_swid_file(records)
        # get the swid Ids
        swids = [mapped_files[file] for file in files if file in mapped_files]
        # list the released files without any swid
        no_swids = [file for file in files if file not in mapped_files]
        if no_swids:
            for file in no_swids:
                print('File {0} in directory {1} does not have a swid'.format(os.path.basename(file), args.run_directory))
    elif args.analysis:
        # get the file swids from the json structure
        infile = open(args.analysis)
        data_structure = json.load(infile)
        infile.close()
        
        # parse FPR records
        # make a list of workflows
        workflows = []
        for i in data_structure:
            for j in data_structure[i]:
                workflows.extend(list(data_structure[i][j].keys()))
        workflows = list(set(workflows))    
        print('workflows', workflows)
        files = parse_fpr_records(provenance, args.project, workflows, args.prefix)
        print('Extracted files from File Provenance Report')
        swids = get_pipeline_swids(data_structure, files)
    else:
        # collect relevant information from File Provenance Report about fastqs for project 
        files = parse_fpr_records(provenance, args.project, [args.workflow], args.prefix)
        print('Information was extracted from FPR {0}'.format(provenance))
        released_files, _  = collect_files_for_release(files, args.release_files, args.nomiseq, args.runs, args.libraries, args.exclude)
        swids = list(released_files.keys())
    
    # mark files in nabu
    change_nabu_status(args.nabu, swids, args.status.upper(), args.user, comment=args.comment)
        
    
if __name__ == '__main__':

    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'dare.py', description='A tool to manage data release')
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')
    
   	# link files in gsi space 
    l_parser = subparsers.add_parser('link', help="Link files extracted from FPR to gsi space")
    l_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The optional second column is the run id')
    l_parser.add_argument('-w', '--workflow', dest='workflow', help='Worflow used to generate the output files')
    l_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    l_parser.add_argument('-p', '--parent', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    l_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report. Used to parse the FPR by project. Files are further filtered by run is runs parameter if provided, or all files for the project and workflow are used', required=True)
    l_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run IDs. Include one or more run Id separated by white space. Other runs are ignored if provided')
    l_parser.add_argument('--exclude_miseq', dest='nomiseq', action='store_true', help='Exclude MiSeq runs if activated')
    l_parser.add_argument('-e', '--exclude', dest='exclude', help='File with libraries tagged for non-release. The first column is always the library. The optional second column is the run id')
    l_parser.add_argument('-s', '--suffix', dest='suffix', help='Indicates if fastqs or datafiles are released by adding suffix to the directory name. Use fastqs or workflow name')
    l_parser.add_argument('-f', '--files', dest='release_files', help='File with file names of the files to be released')
    l_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    l_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    l_parser.add_argument('-a', '--analysis', dest='analysis', help='Path to the file with hierarchical structure storing sample and workflow ids')
    l_parser.set_defaults(func=link_files)
    
   	# map external IDs 
    m_parser = subparsers.add_parser('map', help="Map files to external IDs")
    m_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The optional second column is the run id')
    m_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    m_parser.add_argument('-p', '--parent', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    m_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report. Used to parse the FPR by project. Files are further filtered by run is runs parameter if provided, or all files for the project and workflow are used', required=True)
    m_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run IDs. Include one or more run Id separated by white space. Other runs are ignored if provided')
    m_parser.add_argument('--exclude_miseq', dest='nomiseq', action='store_true', help='Exclude MiSeq runs if activated')
    m_parser.add_argument('-e', '--exclude', dest='exclude', help='File with libraries tagged for non-release. The first column is always the library. The optional second column is the run id')
    m_parser.add_argument('-f', '--files', dest='release_files', help='File with file names to be released')
    m_parser.add_argument('--panel', dest='add_panel', action='store_true', help='Add panel to sample if option is used. By default, panel is not added.')
    m_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    m_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    m_parser.set_defaults(func=map_external_ids)

    # mark files in nabu 
    n_parser = subparsers.add_parser('mark', help="Mark released or withheld files in Nabu")
    n_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report. Used to parse the FPR by project', required = True)
    n_parser.add_argument('-u', '--user', dest='user', help='User name to appear in Nabu for each released or whitheld file', required=True)
    n_parser.add_argument('-st', '--status', dest='status', choices = ['fail', 'pass'], help='Mark files accordingly when released or withheld', required = True)
    n_parser.add_argument('-rn', '--rundir', dest='run_directory', help='Directory with links organized by project and run in gsi space')
    n_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run IDs. Include one or more run Id separated by white space. Other runs are ignored if provided')
    n_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The optional second column is the run id')
    n_parser.add_argument('-w', '--workflow', dest='workflow', help='Worflow used to generate the output files')
    n_parser.add_argument('--exclude_miseq', dest='nomiseq', action='store_true', help='Exclude MiSeq runs if activated')
    n_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    n_parser.add_argument('-e', '--exclude', dest='exclude', help='File with libraries tagged for non-release. The first column is always the library. The optional second column is the run id')
    n_parser.add_argument('-f', '--files', dest='release_files', help='File with file names to be released')
    n_parser.add_argument('-c', '--comment', dest='comment', help='Comment to be added to the released file')
    n_parser.add_argument('-n', '--nabu', dest='nabu', default='https://nabu-prod.gsi.oicr.on.ca', help='URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca')
    n_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    n_parser.add_argument('-a', '--analysis', dest='analysis', help='Path to the file with hierarchical structure storing sample and workflow ids')
    n_parser.set_defaults(func=mark_files_nabu)
        
    # write a report
    r_parser = subparsers.add_parser('report', help="Write a PDF report for released FASTQs")
    r_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report', required=True)
    r_parser.add_argument('-p', '--parents', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    r_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    r_parser.add_argument('-fn', '--full_name', dest='project_full_name', help='Full name of the project', required = True)
    r_parser.add_argument('-rn', '--rundirs', dest='run_directories', nargs='*', help='List of directories with released fastqs')
    r_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run IDs. Include one or more run Id separated by white space. Other runs are ignored if provided')
    r_parser.add_argument('--exclude_miseq', dest='nomiseq', action='store_true', help='Exclude MiSeq runs if activated')
    r_parser.add_argument('-e', '--exclude', dest='exclude', help='File with libraries tagged for non-release. The first column is always the library. The optional second column is the run id')
    r_parser.add_argument('-f', '--files', dest='release_files', help='File with file names to be released')
    r_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The optional second column is the run id')
    r_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    r_parser.add_argument('-a', '--api', dest='api', default='https://nabu-prod.gsi.oicr.on.ca', help='URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca')
    r_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    r_parser.add_argument('-bq', '--bamqc', dest='bamqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest', help='Path to the bamqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest')
    r_parser.add_argument('-dq', '--dnaseqqc', dest='dnaseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/dnaseqqc/latest', help='Path to the dnaseqqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/dnaseqqc/latest')
    r_parser.add_argument('-cq', '--cfmedipqc', dest='cfmedipqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest', help='Path to the cfmedip SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest')
    r_parser.add_argument('-rq', '--rnaseqqc', dest='rnaseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2/latest', help='Path to the rnaseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2/latest')
    r_parser.add_argument('-eq', '--emseqqc', dest='emseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest', help='Path to the emseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest')
    r_parser.add_argument('-u', '--user', dest='user', help='Name of the GSI personnel generating the report', required = True)
    r_parser.add_argument('-t', '--ticket', dest='ticket', help='Jira data release ticket code', required = True)
    r_parser.add_argument('--keep_html', dest='keep_html', action='store_true', help='Write html report if activated.')
    r_parser.add_argument('--panel', dest='add_panel', action='store_true', help='Add panel to sample if option is used. By default, panel is not added.')
    r_parser.set_defaults(func=write_batch_report)
    
    
    # write a report
    c_parser = subparsers.add_parser('cumulative_report', help="Write a cumulative project PDF report for released FASTQs")
    c_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report', required=True)
    c_parser.add_argument('-p', '--parents', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    c_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    c_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    c_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    c_parser.add_argument('-a', '--api', dest='api', default='https://nabu-prod.gsi.oicr.on.ca', help='URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca')
    c_parser.add_argument('-bq', '--bamqc', dest='merged_bamqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/bamqc4merged/latest', help='Path to the merged bamqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/bamqc4merged/latest')
    c_parser.add_argument('-rq', '--rnaseqqc', dest='merged_rnaseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2merged/latest', help='Path to the merged rnaseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2merged/latest')
    c_parser.add_argument('-cq', '--cfmedipqc', dest='cfmedipqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest', help='Path to the cfmedip SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest')
    c_parser.add_argument('-eq', '--emseqqc', dest='emseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest', help='Path to the emseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest')
    c_parser.add_argument('-fn', '--full_name', dest='project_full_name', help='Full name of the project', required = True)
    c_parser.add_argument('-t', '--ticket', dest='ticket', help='Jira data release ticket code', required = True)
    c_parser.add_argument('-u', '--user', dest='user', help='Name of the GSI personnel generating the report', required = True)
    c_parser.set_defaults(func=write_cumulative_report)

    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)