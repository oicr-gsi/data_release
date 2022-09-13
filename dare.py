# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 12:27:07 2020

@author: rjovelin
"""


import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import time
#import xhtml2pdf.pisa as pisa
import mistune
import base64
from PIL import Image
import math
import requests
import gzip
import io
import sys
import json
import pathlib
import sqlite3
from jinja2 import Environment, FileSystemLoader



def get_libraries(library_file):
    '''
    (str) -> dict
    
    Returns a dictionary with library, run key, value pairs, or a dictionary of library, empty string 
    key, value pairs if runs are not specified. 
    Note: runs need to be specified or missing for all libraries
    
    Parameters
    ----------
    - sample_file (str): Path to sample file. Sample file is a tab delimited file
                         that includes 1 or 2 columns. The first column is always 
                         the library alias, and the second optional column is run id
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

    L = list(D.values())
    # check that run id is specified or missing for all libraries
    if not (all(map(lambda x: x != '', L)) or all(map(lambda x: x == '', L))):
        raise ValueError('Run id must be specified or absent for all libraries')    
    
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
                if workflow[0] == 'bcl2fastq':
                    # skip if not fastq-related workflows    
                    if pipeline_workflow.lower() not in ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport']:
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
                 'panel': geo['geo_targeted_resequencing'], 'library_source': geo['geo_library_source_template_type'],
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

    # make a list of valid runs if runs are specified in valid_libraries
    valid_runs = list(valid_libraries.values())
    
    for file_swid in files:
        libraries = files[file_swid]['library']
        if set(libraries).intersection(valid_libraries.keys()) != set(libraries):
            exclude.append(file_swid)
        # check if runs are defined
        if all(valid_runs):
            runs = files[file_swid]['run_id']
            if set(runs).intersection(set(valid_runs)) != set(runs):
                exclude.append(file_swid)
        
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

    # make a list of runs to be excluded
    exclude_runs = list(exclude.values())
    
    # make a list of files to exclude
    L = []
    
    D = {}
    for file_swid in files:
        library = files[file_swid]['library']
        if set(library).intersection(exclude.keys()):
            L.append(file_swid)
        # check if runs are defined
        if all(exclude_runs):
            runs = files[file_swid]['run_id']
            if set(runs).intersection(set(exclude_runs)):
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
      
    for sample_name in pipeline_data:
        sample_dir = os.path.join(working_dir, sample_name)
        os.makedirs(sample_dir, exist_ok=True)
        
        # create sample dir
        for workflow_name in pipeline_data[sample_name]:
            # create workflow dir
            workflow_dir = os.path.join(sample_dir, workflow_name)
            os.makedirs(workflow_dir, exist_ok=True)
            for i in pipeline_data[sample_name][workflow_name]:
                file = i['file_path']
                if 'file_name' in i:
                    # use new file name to name link
                    filename = i['file_name']
                else:
                    filename = os.path.basename(file)
                link = os.path.join(workflow_dir, filename)
                if os.path.isfile(link) == False:
                    os.symlink(file, link)



def collect_md5sums(files):
    '''
    (dict) -> dict
    
    Returns a dictionary with lists of md5sums, file_path for each run in dictionary files
        
    Parameters
    ----------
    - files (dict) : Dictionary with file records obtained from parsing FPR
    '''
    
    # create a dictionary {run: [md5sum, file_path]} 
    D = {}
    for file_swid in files:
        run_name = ';'.join(files[file_swid]['run_id'])
        md5 = files[file_swid]['md5']
        file_path = files[file_swid]['file_path']
        if run_name in D:
            D[run_name].append([md5, file_path])
        else:
            D[run_name] = [[md5, file_path]]
    return D


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


def collect_files_for_release(files, release_files, nomiseq, runs, libraries, exclude, suffix):
    '''
    (dict, str | None, bool, list | None, str | None, str | None, str) -> (dict, dict)
    
    Returns dictionaries with file records for files tagged for release and files that should not be released, if any, or an empty dictionary.
        
    Parameters
    ----------
    - files (str): Dictionary with file records for a entire project or a specific workflow for a given project
                        Default is '/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz'
    - release_files (str | None): Path to file with file names to be released 
    - nomiseq (bool): Exclude MiSeq runs if True
    - runs (list | None): List of run IDs. Include one or more run Id separated by white space.
                          Other runs are ignored if provided
    - libraries (str | None): Path to 1 or 2 columns tab-delimited file with library IDs.
                              The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM).
                              The second and optional column is the library aliquot ID (eg. LDI32439).
                              Only the samples with these library aliases are used if provided'
    - exclude (str | None): File with sample name or libraries to exclude from the release
    - suffix (str): Indicates map for fastqs or datafiles in the output file name
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
    (dict, dict)    
    
    Returns a dictionary with the files for each sample and workflow specified in the data_structure
    
    Parameters
    ----------
    - data_structure (str): Dictionary with samples, workflows and workflow_run_id hierarchical structure
    - files (dict): Dictionary with all file records for a given project extracted from FPR
    '''

    
    D = {}
    
    for file_swid in files:
        sample = files[file_swid]['sample_name']
        workflow = files[file_swid]['workflow']
        if sample in data_structure:
            if workflow in data_structure[sample]:
                if files[file_swid]['workflow_version'] == data_structure[sample][workflow]['workflow_version']:
                    # get workflow run id
                    if files[file_swid]['workflow_run_id'] == data_structure[sample][workflow]['workflow_id']:
                        if sample not in D:
                            D[sample] = {}
                        if 'name' in data_structure[sample][workflow]:
                            workflow_name = data_structure[sample][workflow]['name']
                        else:
                            workflow_name = workflow
                        if workflow_name not in D[sample]:
                            D[sample][workflow_name] = []
                        
                        # check which files are collected
                        if 'extension' in data_structure[sample][workflow]:
                            # files are collected based on file extension
                            # get the file extension
                            extension = pathlib.Path(files[file_swid]['file_path']).suffix
                            if extension in data_structure[sample][workflow]['extension']:
                                D[sample][workflow_name].append({'file_path': files[file_swid]['file_path'], 'md5': files[file_swid]['md5']})
                        elif 'files' in data_structure[sample][workflow]:
                            # files are collected based on file name
                            if os.path.basename(files[file_swid]['file_path']) in list(map(lambda x: os.path.basename(x),  data_structure[sample][workflow]['files'])):
                                D[sample][workflow_name].append({'file_path': files[file_swid]['file_path'], 'md5': files[file_swid]['md5']})
                        elif 'rename_files' in data_structure[sample][workflow]:
                            # files are collected based on file name and renamed
                            # make a list of file names
                            file_paths = list(map(lambda x: os.path.basename(x), [i['file_path'] for i in data_structure[sample][workflow]['rename_files']]))
                            file_names = [i['file_name'] for i in data_structure[sample][workflow]['rename_files']]
                            if os.path.basename(files[file_swid]['file_path']) in file_paths:
                                # collect file path, md5 and new file name used to name the link
                                j = file_paths.index(os.path.basename(files[file_swid]['file_path']))
                                D[sample][workflow_name].append({'file_path': files[file_swid]['file_path'],
                                                                 'md5': files[file_swid]['md5'],
                                                                 'file_name': file_names[j]})
                        else:
                            D[sample][workflow_name].append({'file_path': files[file_swid]['file_path'], 'md5': files[file_swid]['md5']})
    
    return D                        
                        
    
    
    
    
def write_md5sum(data, outputfile):
    '''
    (dict, str, str, str | None) -> None   
    
    Write a file in working_dir with md5sums for all files contained in data
    
    Parameters
    ----------
    - data (dict): Dictionary holding pipeline data or data from a single workflow
    - outpufile (str): Path to the outputfile with md5sums
    '''    
    
    newfile = open(outputfile, 'w')
    # check the data structure to determine if it contains pipeline analysis data or a single workflow data 
    if 'md5' in data[list(data.keys())[0]]:
        # single workflow data
        for file_swid in data:
            md5 = data[file_swid]['md5']
            file_path = data[file_swid]['file_path']
            newfile.write('\t'.join([file_path, md5]) + '\n')
    else:
        # pipeline data
        for sample in data:
            for workflow_name in data[sample]:
                for i in data[sample][workflow_name]:
                    newfile.write('\t'.join([i['file_path'], i['md5']]) + '\n')
    newfile.close()

    
    
def link_files(args):
    '''
    (str, str, str, str | None, str | None, bool, list | None, str | None, str | None, str, str, str, str | None)
    
    Parameters
    ----------
    - provenance (str): Path to File Provenance Report.
                        Default is '/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz'
    - project (str): Project name as it appears in File Provenance Report.
    - workflow (str): Worflow used to generate the output files
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
        files, files_non_release = collect_files_for_release(files, args.release_files, args.nomiseq, args.runs, args.libraries, args.exclude, args.suffix)
        # link files to project dir
        if args.suffix == 'fastqs':
            assert args.workflow.lower() in ['bcl2fastq', 'casava', 'fileimport', 'fileimportforanalysis']
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
            workflows.extend(list(data_structure[i].keys()))
        workflows = list(set(workflows))    
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
        write_md5sum(pipeline_data, outputfile)
    print('Files were extracted from FPR {0}'.format(provenance))


def group_sample_info_mapping(files, add_time_points, add_panel):
    '''
    (dict, dict, bool, bool) -> dict    
    
    Returns a dictionary of sample information organized by run
        
    Parameters
    ----------
    - files (dict) : Dictionary with file records obtained from parsing FPR
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
        library_source = files[file_swid]['library_source']
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
        timepoint = files[file_swid]['time_point']
        
        # group files by sample, library, run and lane
        key = '_'.join([sample, library, run, barcode])
        
        
        L = [sample, external_id, library, library_source, tissue_type, tissue_origin, run, barcode, group_id, group_description]
        if add_time_points:
            L.append(timepoint)
        if add_panel:
            L.append(panel)
            
        if key not in D:
            D[key] = {'info': L, 'files': [file_path]}
        else:
            assert D[key]['info'] == L
            D[key]['files'].append(file_path)
        
    return D            


def extract_sample_info(sample_provenance):
    '''
    (str) -> list
    
    Returns a list of dictionary with sample information pulled down from the
    sample_provenance Pinary API

    Parameters
    ----------
    - sample_provenance (str): Pinery API, http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance
    '''
    
    response = requests.get(sample_provenance)
    if response.ok:
        L = response.json()
    else:
        L = []
    return L


def get_time_points(sample_information):
    '''
    (list) -> dict
    
    Returns a dictionary of time points for each library
    
    Parameters
    ----------
    - sample_information (list): List of dictionary with sample information pulled
                                 down from the inery API 
    '''

    D = {}

    for i in sample_information:
        sample = i['sampleName']
        if 'timepoint' in i['sampleAttributes']:
            time_point = i['sampleAttributes']['timepoint']
        else:
            time_point = []
        if sample in D:
            D[sample].extend(time_point)
        else:
            D[sample] = time_point
        
    for i in D:
        D[i] = ';'.join(sorted(list(set(D[i]))))
        if not D[i]:
            D[i] = 'NA'
    return D



def add_time_points(sample_provenance, files):
    '''
    (str, dict) -> None
    
    Add time points retrieved from sample_provenance to each sample in sample_metrics.
    Updates files in place
    
    Parameters
    ----------
    - sample_provenance (str): Pinery API, http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance
    - files (dict): Dictionary with file information extracted from FPR
    '''

    
    time_points = get_time_points(extract_sample_info(sample_provenance))
    # add time points
    for file_swid in files:
        assert len(files[file_swid]['library']) == 1
        library = files[file_swid]['library'][0]
        if library in time_points:
            files[file_swid]['time_point'] = time_points[library]
        else:
            files[file_swid]['time_point'] = 'NA'




def map_external_ids(args):
    '''
    (str, str, str, str | None, str | None, bool, list | None, str | None, str | None, str, str, str, bool, bool, str) -> None
    
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
    - timepoints (bool): Add time points column to sample map if True
    - add_panel (bool): Add panel column to sample map if True
    - sample_provenance (str): Pinery API, http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance
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
    files, files_non_release = collect_files_for_release(files, args.release_files, args.nomiseq, args.runs, args.libraries, args.exclude, suffix)
    
    # add time points
    add_time_points(args.sample_provenance, files)
    
    # group sample information by run            
    sample_info = group_sample_info_mapping(files, args.timepoints, args.add_panel)        
    
    # write sample maps
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    current_time = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
    outputfile = os.path.join(working_dir, '{0}.release.{1}.{2}.map.tsv'.format(args.project, current_time, suffix))
    newfile = open(outputfile, 'w')
    header = ['sample', 'sample_id', 'library', 'library_source', 'tissue_type', 'tissue_origin', 'run', 'barcode', 'group_id', 'group_description', 'files']
    if args.timepoints:
       #header.append('time_points')
       header.insert(-1, 'time_points')  
    if args.add_panel:
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
    print('Generated sample maps in {0}'.format(working_dir))
    print('Information was extracted from FPR {0}'.format(provenance))


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


# def parse_bamqc(bamqc_table, project):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with project-level relevant information from the bamqc table of qc-etl
        
#     Parameters
#     ----------
#     - bamqc_table (str): Comma-separated table from bamqc generated by qc-etl
#     - project (str): Specific project of interest
#     '''
            
#     # get versions for bamqc4
#     #version = subprocess.check_output('/.mounts/labs/gsi/modulator/sw/Ubuntu18.04/gsi-qc-etl-0.53.0/bin/gsiqcetl versions bamqc4', shell=True).decode('utf-8').rstrip()
    
#     # get tables
#     #tables = subprocess.check_output('/.mounts/labs/gsi/modulator/sw/Ubuntu18.04/gsi-qc-etl-0.53.0/bin/gsiqcetl tables bamqc4 {0}'.format(version), shell=True).decode('utf-8').rstrip()
    
#     # get bamqc records
#     #bamqc_records =  subprocess.check_output('/.mounts/labs/gsi/modulator/sw/Ubuntu18.04/gsi-qc-etl-0.53.0/bin/gsiqcetl dump -f /scratch2/groups/gsi/production/qcetl/bamqc4/latest bamqc4 5 bamqc4', shell=True).decode('utf-8').rstrip().split('\n')
    
#     infile = open(bamqc_table)
#     bamqc_records = infile.read().rstrip().split('\n')
#     header = bamqc_records.pop(0).split(',')
#     infile.close()
    
#     # initiate dictionary, collect information for each file at the run level  
#     B = {}
    
#     # loop over bamqc records, collect relevant information
#     for j in bamqc_records:
#         if project in j:
#             i = j.split(',')
#             # some fields have internal commas. need to stitch back entry for that field
#             if len(i) != len(header):
#                 i =  stich_back_etl_record(i)
#             # check that fields have been stiched back correctly
#             assert len(i) == len(header)
#             # collect relevant information
#             run_alias = i[header.index('Run Alias')]
#             instrument = i[header.index('instrument')]
#             bases_mapped = int(i[header.index('bases mapped')])
#             coverage = float(i[header.index('coverage')])
#             coverage_dedup = float(i[header.index('coverage deduplicated')])
#             library = i[header.index('library')]
#             mapped_reads = int(i[header.index('mapped reads')])
#             barcodes = i[header.index('Barcodes')]
#             lane = i[header.index('Lane Number')]
#             sample = i[header.index('sample')]
#             percent_duplicate = float(i[header.index('mark duplicates_PERCENT_DUPLICATION')])
            
#             total_bases_on_target = int(i[header.index('total bases on target')])
#             total_reads = int(i[header.index('total reads')])
            
#             on_target = compute_on_target_rate(bases_mapped, total_bases_on_target)    
            
#             # map file information
#             d = {'run_alias': run_alias, 'instrument': instrument,
#                  'bases_mapped': bases_mapped, 'coverage': coverage,
#                  'coverage_dedup': coverage_dedup, 'library': library,
#                  'mapped_reads': mapped_reads, 'sample': sample,
#                  'total_bases_on_target': total_bases_on_target,
#                  'total_reads': total_reads, 'on_target': on_target,
#                  'barcodes':barcodes, 'lane': lane, 'percent_duplicate': percent_duplicate}
#             # collect information at the run level
#             if run_alias in B:
#                 B[run_alias].append(d)
#             else:
#                 B[run_alias] = [d]
#     return B


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
            case = i[header.index('Donor')]
            library_design = i[header.index('Library Design')]
            tissue_origin = i[header.index('Tissue Origin')]
            tissue_type = i[header.index('Tissue Type')]
            sample_name = case + '_' + tissue_origin + '_' + tissue_type + '_' + library_design
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
        run_alias = FPR_info[file]['run_alias']
        run = FPR_info[file]['run']
        library = FPR_info[file]['lid']
        instrument = FPR_info[file]['instrument']
        barcode = FPR_info[file]['barcode']
        ext_id = FPR_info[file]['external_id']
        case = FPR_info[file]['ID']
        library_source = FPR_info[file]['library_source']
        tissue_origin = FPR_info[file]['tissue_origin']
        tissue_type = FPR_info[file]['tissue_type']
        read_count = FPR_info[file]['read_count']
        coverage = FPR_info[file]['coverage']
        coverage_dedup = FPR_info[file]['coverage_dedup']
        on_target = FPR_info[file]['on_target']
        duplicate = FPR_info[file]['percent_duplicate']
        group_id = FPR_info[file]['group_id']


        if instrument not in D:
            D[instrument] = {}
        
        if sample not in D[instrument]:
            D[instrument][sample] = [{'sample': sample, 'lane': lane, 'run': run,
                                      'run_alias': run_alias, 'library': library,
                                      'instrument': instrument, 'barcode': barcode, 'ext_id': ext_id,
                                      'case': case, 'library_source': library_source,
                                      'tissue_type': tissue_type, 'tissue_origin': tissue_origin,
                                      'reads': read_count, 'coverage': coverage,
                                      'coverage_dedup': coverage_dedup, 'on_target': on_target,
                                      'duplicate (%)': duplicate, 'files': [file], 'group_id': group_id}]
        else:
            # find paired fastq
            for d in D[instrument][sample]:
                if run == d['run'] and lane == d['lane'] and library == d['library'] \
                     and instrument == d['instrument'] and barcode == d['barcode'] \
                     and ext_id == d['ext_id'] and case == d['case'] and library_source == d['library_source']:
                         d['reads'] += read_count
                         assert group_id == d['group_id']
                         assert coverage == d['coverage']
                         assert coverage_dedup == d['coverage_dedup']
                         assert duplicate == d['duplicate (%)']
                         assert on_target == d['on_target']
                         assert tissue_origin == d['tissue_origin']
                         assert tissue_type == d['tissue_type']
                         d['files'].append(file)
                         # update boolean. paired end fastq found
                         found = True
            if found == False:
                # record file info
                D[instrument][sample].append({'sample': sample, 'lane': lane, 'run': run,
                             'run_alias': run_alias, 'library': library,
                             'instrument': instrument, 'barcode': barcode, 'ext_id': ext_id,
                             'case': case, 'library_source': library_source,
                             'tissue_type': tissue_type, 'tissue_origin': tissue_origin,
                             'reads': read_count, 'coverage': coverage,
                             'coverage_dedup': coverage_dedup, 'on_target': on_target,
                             'duplicate (%)': duplicate, 'files': [file], 'group_id': group_id})
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
        run_alias = FPR_info[file]['run_alias']
        run = FPR_info[file]['run']
        library = FPR_info[file]['lid']
        instrument = FPR_info[file]['instrument']
        barcode = FPR_info[file]['barcode']
        ext_id = FPR_info[file]['external_id']
        case = FPR_info[file]['ID']
        library_source = FPR_info[file]['library_source']
        tissue_origin = FPR_info[file]['tissue_origin']
        tissue_type = FPR_info[file]['tissue_type']
        read_count = FPR_info[file]['read_count']
        coverage = FPR_info[file]['coverage']
        coverage_dedup = FPR_info[file]['coverage_dedup']
        on_target = FPR_info[file]['on_target']
        duplicate = FPR_info[file]['percent_duplicate']
        
        group_id = FPR_info[file]['group_id']
        
        if instrument not in D:
            D[instrument] = {}
               
        if sample not in D[instrument]:
            D[instrument][sample] = {'sample': sample, 'lane': [lane], 'run': [run],
                              'run_alias': [run_alias], 'library': [library],
                              'instrument': instrument, 'barcode': [barcode], 'ext_id': ext_id,
                              'case': case, 'library_source': library_source,
                              'tissue_origin': tissue_origin, 'tissue_type': tissue_type,
                              'reads': read_count, 'coverage': coverage,
                              'coverage_dedup': coverage_dedup, 'on_target': on_target,
                              'duplicate (%)': duplicate, 'files': [file], 'group_id': group_id}
        else:
            assert ext_id == D[instrument][sample]['ext_id']
            assert case == D[instrument][sample]['case']
            assert library_source == D[instrument][sample]['library_source']
            assert tissue_type == D[instrument][sample]['tissue_type']
            assert tissue_origin == D[instrument][sample]['tissue_origin']
            D[instrument][sample]['library'].append(library)  
            D[instrument][sample]['reads'] += read_count
            D[instrument][sample]['files'].append(file)
            D[instrument][sample]['run'].append(run)  
            D[instrument][sample]['lane'].append(lane)  
            D[instrument][sample]['run_alias'].append(run_alias)  
            
            D[instrument][sample]['barcode'].append(barcode)
            
            assert coverage == D[instrument][sample]['coverage']
            assert coverage_dedup == D[instrument][sample]['coverage_dedup']
            assert duplicate == D[instrument][sample]['duplicate (%)']
            assert on_target == D[instrument][sample]['on_target']
    
            assert group_id == D[instrument][sample]['group_id']
    
    
    
    # collapse lanes, runs and libraries
    for instrument in D:
        for sample in D[instrument]:
            D[instrument][sample]['run'] = ';'.join(list(set(D[instrument][sample]['run'])))
            D[instrument][sample]['lane'] = ';'.join(list(set(D[instrument][sample]['lane'])))
            D[instrument][sample]['run_alias'] = ';'.join(list(set(D[instrument][sample]['run_alias'])))
            D[instrument][sample]['library'] = ';'.join(list(set(D[instrument][sample]['library'])))
            D[instrument][sample]['barcode'] = ';'.join(list(set(D[instrument][sample]['barcode'])))
            
    return D                         

                                           
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
                    D[instrument] = [[sample_metrics[instrument][sample][metric], sample]]
    return D



def sort_metric2_according_to_metric1_order(sample_metrics, instrument, metric2, samples1, level):
    '''
    (dict, str, str, list, str) -> list
    
    Returns a list of metric2 sorted according to order of metric1 
    
    Parameters
    ----------
    - sample_metrics (dict): Run-level or cumulative QC metrics for all samples 
    - instrument (str): Sequencing intrument
    - metric2 (str): Name of QC metric 2
    - samples1 (str): List of samples sorted according to the order of metric 1
    - level (str): Report level: single or cumulative
    '''

    # group metric 2 by instrument
    D = group_qc_metric_by_instrument(sample_metrics, metric2, level)
    # make a list of metric2 sorted according to the order of metric1
    Q = []
    for i in samples1:
        for j in D[instrument]:
            if j[1] == i:
                Q.append(j[0])
    return Q

    
    
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
    - sample_metrics (dict): Dictionary with run-level sample metrics
    - header (list): List of column names
    - column_size (dict): Dictionary with columns width for each column in header
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
                    if i == 'Reads':
                        # add commas to read count
                        j = '{:,}'.format(d['reads'])
                    elif i == 'On target':
                        j = str(d['on_target'])
                    elif i == 'External ID':
                        j = str(d['ext_id'])
                    elif i == 'Group ID':
                        j = str(d['group_id'])
                    elif i == 'Library ID':
                        j = str(d['library'])
                    elif i == 'Library ID (time point)':
                        j = str(d['library']) + '\n' + '({0})'.format(str(d['timepoint']))
                    elif i == 'Run':
                        j = str(d[i.lower()])
                        if ';' in j:
                            j = j.replace(';', ';\n')
                    elif i == 'Indices':
                        j = str(d['barcode'])
                        if '-' in j:
                            j = j.replace('-', '-\n')
                    elif i == 'Sequencing Run':
                        j = str(d['run']) + '\n' + str(d['barcode']) 
                    elif i == 'O':
                        j = str(d['tissue_origin'])
                    elif i == 'S':
                        j = str(d['library_source'])
                    elif i == 'T':
                        j = str(d['tissue_type'])
                    else:
                        j = str(d[i.lower()])
                    if counter + 1 == cells:
                        table.append('<td style="border-bottom: 1px solid #000000; padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
                    else:
                        table.append('<td style="padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
                table.append('</tr>')
                # update counter
                counter += 1
    table.append('</table>')
    return ''.join(table)




def generate_cumulative_table(sample_metrics, header, column_size, table_type=None):
    '''
    (dict, list, dict, str, str) -> str
    
    Returns a html string representing a table
    
    Parameters
    ----------
    - sample_metrics (dict): Dictionary with cumulative sample metrics
    - header (list): List of column names
    - column_size (dict): Dictionary with columns width for each column in header
    - table_type (str): Report library and lane counts rather than Ids if set to metrics    
    '''
    
    # count the expected number of cells (excluding header) in tables
    cells = 0
    for instrument in sample_metrics:
        cells += len(sample_metrics[instrument].keys())
      
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
            if counter % 2 == 0:
                table.append('<tr style="background-color: #eee">')
            else:
                table.append('<tr style="background-color: #fff">')
            for i in header:
                if i == 'Reads':
                    j = '{:,}'.format(sample_metrics[instrument][sample][i.lower()])    
                elif i == 'Sample':
                    j = sample
                    k = len(j.split('_')) // 4
                    j = '_'.join(j.split('_')[0:k]) + '_\n' + '_'.join(j.split('_')[k:2*k]) + '_\n' + '_'.join(j.split('_')[2*k:3*k]) + '_\n' + '_'.join(j.split('_')[3*k:])
                elif i == 'External ID':
                    j = str(sample_metrics[instrument][sample]['ext_id'])
                elif i == 'S':
                    j = str(sample_metrics[instrument][sample]['library_source'])
                elif i == 'T':
                    j = str(sample_metrics[instrument][sample]['tissue_type'])
                elif i == 'O':
                    j = str(sample_metrics[instrument][sample]['tissue_origin'])
                elif i == 'Group ID':
                    j = str(sample_metrics[instrument][sample]['group_id'])
                    if table_type != 'metrics':
                        k = len(j.split('_')) // 2
                        j = '_'.join(j.split('_')[0:k]) + ' \n' + '_'.join(j.split('_')[k:]) 
                elif i == 'Library ID' or i == 'Libraries':
                    j = str(sample_metrics[instrument][sample]['library'])
                elif i == 'Runs':
                    j = str(sample_metrics[instrument][sample]['run'])
                elif i == 'Indices':
                    j = str(sample_metrics[instrument][sample]['barcode'])
                else:
                    j = str(sample_metrics[instrument][sample][i.lower()])
                    
                if i in ['Indices', 'Run', 'Library ID']:
                    if '-' in j:
                        j = j.replace('-', '-\n')
                    if ';' in j:
                        j = j.replace(';', ';\n')
                                        
                if table_type == 'metrics':
                    if i == 'Runs':
                        j = len(list(set(str(sample_metrics[instrument][sample]['run']).split(';'))))
                    elif i == 'Libraries':
                        j = len(list(set(str(sample_metrics[instrument][sample]['library']).split(';'))))
                    #j = len(list(set(j.split(';'))))
                                                      
                if counter + 1 == cells:
                    table.append('<td style="border-bottom: 1px solid #000000; padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
                else:
                    table.append('<td style="padding: {0}; font-size: 10px; text-align: left;">{1}</td>'.format(padding, j))
            table.append('</tr>')
            # update counter
            counter += 1
    table.append('</table>')
    return ''.join(table)




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


def generate_project_table(project_name, project_full_name, current_date):
    '''
    (str, str, str, str, str) -> str

    Returns a html string representing a table with project information

    Parameters
    ----------
    - project_name (str): Acronym of the project
    - project_full_name (str): Full name of the project 
    - current_date (str): Date of the release (Y-M-D)
    '''

    content = [project_name, project_full_name, current_date]
    column_width = [30, 50, 20]

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('<tr>')
    for i in ['Project Acronym', 'Project Name', 'Date']:
        table.append('<th style="padding: 3px; border: 0.75px solid grey; border-collapse:collapse; text-align: center; font-size:10px">{0}</th>'.format(i))
    table.append('</tr>')
    # add lines in table
    table.append('<tr>')
    for i in range(len(content)):
        table.append('<td style="width: {0}%; padding: 3px; border: 0.75px solid grey;  border-collapse:collapse; text-align: center; font-size:10px">{1}</td>'.format(column_width[i], content[i]))
    table.append('</tr>')
    table.append('</table>')
    return ''.join(table)



def generate_figure_table(file1, factor, file2=None):
   
    '''
    (str, str | None) -> str

    Returns a html string representing a table with figure files

    Parameters
    ----------
    - file1 (str): Path to the figure file with metrics 1 and 2 
    - file2 (str | None): Path to the figure file with metrics 3 and 4 if exists
    - factor (float): Resizing factor
    '''

    table = []
    # add table style
    table.append('<table style="width:100%; font-family: Arial, Helvetica, sans-serif">')
    # add header
    table.append('</tr>')
    table.append('<tr>')
    # resize figure 1
    height1, width1 = resize_image(file1, factor)
    if file2:
        width = '50%'
    else:
        width = '100%'
    table.append('<td style="width: {0}; padding: 3px; text-align: left"><img src="{1}" alt="{2}" title="{2}" style="padding-right: 0px; padding-left:0px; width:{3}; height:{4}"></td>'.format(width, file1, os.path.basename(file1), width1, height1))
    if file2:
        # resize figure 
        height2, width2 = resize_image(file2, factor)
        table.append('<td style="width: {0}; padding: 3px; text-align: left"><img src="{1}" alt="{2}" title="{2}" style="padding-right: 0px; padding-left:0px; width:{3}; height:{4}"></td>'.format(width, file2, os.path.basename(file2), width2, height2))
    table.append('</tr>')
    table.append('</table>')
    
    return ''.join(table)



def count_all_files(fastq_counts):
    '''
    (dict) -> int
    
    Returns the total number of fast files (or number of paired fastq files) across instruments and runs
    
    Parameters
    ----------
    - fastq_counts (dict): Counts of released fastqs for each run and instrument
    '''

    total = 0
    for instrument in fastq_counts:
        for run in fastq_counts[instrument]:
            total += fastq_counts[instrument][run]
    return total



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


def write_md5sums(md5_file, md5sums):
    '''
    (str, str, str, dict, bool) -> str
    
    Writes a table file with file names and md5sums and returns thefile path
        
    Parameters
    ----------
    - md5_file (str): Path to the outputfile with md5sums
    - md5sums (dict): Dictionary of file path, md5sum key, value pairs
    '''

    newfile = open(md5_file, 'w')
    for i in md5sums:
        newfile.write('\t'.join([i[0], i[1]]) + '\n')
    newfile.close()       


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
    
    D = {'Ab': 'Abdomen', 'Ad': 'Adipose', 'Ae': 'Adnexa', 'Ag': 'Adrenal',
         'An': 'Anus', 'Ao': 'Anorectal', 'Ap': 'Appendix', 'As': 'Ascites',
         'At': 'Astrocytoma', 'Av': 'Ampulla', 'Ax': 'Axillary', 'Ba': 'Back',
         'Bd': 'Bile', 'Bi': 'Biliary', 'Bl': 'Bladder', 'Bm': 'Bone', 'Bn': 'Brain',
         'Bo': 'Bone', 'Br': 'Breast', 'Bu': 'Buccal', 'Bw': 'Bowel', 'Cb': 'Cord',
         'Cc': 'Cecum', 'Ce': 'Cervix', 'Cf': 'Cell-Free', 'Ch': 'Chest', 'Cj': 'Conjunctiva',
         'Ck': 'Cheek', 'Cn': 'Central', 'Co': 'Colon', 'Cr': 'Colorectal', 'Cs': 'Cul-de-sac',
         'Ct': 'Circulating', 'Di': 'Diaphragm', 'Du': 'Duodenum', 'En': 'Endometrial',
         'Ep': 'Epidural', 'Es': 'Esophagus', 'Ey': 'Eye', 'Fa': 'Fallopian',
         'Fb': 'Fibroid', 'Fs': 'Foreskin', 'Ft': 'Foot', 'Ga': 'Gastric',
         'Gb': 'Gallbladder', 'Ge': 'Gastroesophageal', 'Gi': 'Gastrointestinal',
         'Gj': 'Gastrojejunal', 'Gn': 'Gingiva', 'Gt': 'Genital', 'Hp': 'Hypopharynx',
         'Hr': 'Heart', 'Ic': 'ileocecum', 'Il': 'Ileum', 'Ki': 'Kidney', 'La': 'Lacrimal',
         'Lb': 'Limb', 'Le': 'Leukocyte', 'Lg': 'Leg', 'Li': 'Large', 'Ln': 'Lymph',
         'Lp': 'Lymphoblast', 'Lu': 'Lung', 'Lv': 'Liver', 'Lx': 'Larynx',
         'Ly': 'Lymphocyte', 'Md': 'Mediastinum', 'Me': 'Mesenchyme', 'Mn': 'Mandible',
         'Mo': 'Mouth', 'Ms': 'Mesentary', 'Mu': 'Muscle', 'Mx': 'Maxilla',
         'Nk': 'Neck', 'nn': 'Unknown', 'No': 'Nose', 'Np': 'Nasopharynx',
         'Oc': 'Oral', 'Om': 'Omentum', 'Or': 'Orbit', 'Ov': 'Ovary',
         'Pa': 'Pancreas', 'Pb': 'Peripheral', 'Pc': 'Pancreatobiliary',
         'Pd': 'Parathyroid', 'Pe': 'Pelvic', 'Pg': 'Parotid', 'Ph': 'Paratracheal',
         'Pi': 'Penis', 'Pl': 'Plasma', 'Pm': 'Peritoneum', 'Pn': 'Peripheral',
         'Po': 'Peri-aorta', 'Pr': 'Prostate', 'Pt': 'Palate', 'Pu': 'Pleura',
         'Py': 'periampullary', 'Ra': 'Right', 'Rc': 'Rectosigmoid', 'Re': 'Rectum',
         'Ri': 'Rib', 'Rp': 'Retroperitoneum', 'Sa': 'Saliva', 'Sb': 'Small',
         'Sc': 'Scalp', 'Se': 'Serum', 'Sg': 'Salivary', 'Si': 'Small', 'Sk': 'Skin',
         'Sm': 'Skeletal', 'Sn': 'Spine', 'So': 'Soft', 'Sp': 'Spleen', 'Sr': 'Serosa',
         'Ss': 'Sinus', 'St': 'Stomach', 'Su': 'Sternum', 'Ta': 'Tail', 'Te': 'Testes',
         'Tg': 'Thymic', 'Th': 'Thymus', 'Tn': 'Tonsil', 'To': 'Throat', 'Tr': 'Trachea',
         'Tu': 'Tongue', 'Ty': 'Thyroid', 'Uc': 'Urachus', 'Ue': 'Ureter', 'Um': 'Umbilical',
         'Up': 'Urine', 'Ur': 'Urethra', 'Us': 'Urine', 'Ut': 'Uterus', 'Uw': 'Urine',
         'Vg': 'Vagina', 'Vu': 'Vulva', 'Wm': 'Worm'}

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
         'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq'}
    
    return D




  
 
def extract_bamqc_data(bamqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the bamqc table of qc-etl
        
    Parameters
    ----------
    - bamqc_table (str): Comma-separated table from bamqc generated by qc-etl
    - project (str): Specific project of interest
    '''
            
    conn = sqlite3.connect(bamqc_db)
    cur = conn.cursor()
    
    # get table name
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    table_name = [i[0] for i in cur][0]
    cur.execute('select * from {0}'.format(table_name))
    header = [i[0] for i in cur.description]
    
    data = cur.fetchall()
    conn.close()
    
    columns = ['sample', 'Run Alias', 'instrument', 'library', 'Barcodes', 
               'coverage', 'coverage deduplicated', 'mark duplicates_PERCENT_DUPLICATION',
               'mapped reads', 'total bases on target', 'total reads', 'bases mapped', 'Lane Number']
    
    D = {}
    for i in data:
        d = {}
        for j in range(len(columns)):
            key = columns[j]
            if j < 5:
                d[key] = i[header.index(key)]
            elif 5 <= j <= 7:
                d[key] = float(i[header.index(key)])
            elif j > 7:
                d[key] = int(i[header.index(key)])    
        
        # compute on_target rate, not available through qc-etl
        d['on_target'] = compute_on_target_rate(d['bases mapped'], d['total bases on target']) 
                        
        run = i[header.index('Run Alias')]
        sample = i[header.index('sample')]
        
        if run not in D:
            D[run] = {}
        if sample not in D[run]:
            D[run][sample] = []
        D[run][sample].append(d)
    return D
     


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
    
    for file_swid in FPR_info:
        qc_found = False
        run_alias = FPR_info[file_swid]['run_id'][0]
        sample_id = FPR_info[file_swid]['sample_id'][0]
        # check that run in recorded in bamqc
        if run_alias in bamqc_info:
            if sample_id in bamqc_info[run_alias]:
                # map file info with bamqc info
                for d in bamqc_info[run_alias][sample_id]:
                    if FPR_info[file_swid]['lane'] == d['lane'] \
                        and FPR_info[file_swid]['sample_name'][0] == d['sample'] \
                        and FPR_info[file_swid]['barcode'][0] == d['barcodes'] \
                        and FPR_info[file_swid]['platform'].replace('_', ' ') == d['instrument']:
                            assert FPR_info[file_swid]['library'][0] == d['library']    
                            qc_found = True
                            FPR_info[file_swid]['coverage'] = round(d['coverage'], 2)
                            FPR_info[file_swid]['coverage_dedup'] = round(d['coverage_dedup'], 2)
                            FPR_info[file_swid]['on_target'] = round(d['on_target'], 2)                
                            FPR_info[file_swid]['percent_duplicate'] = round(d['percent_duplicate'], 2)

        if qc_found == False:
            FPR_info[file_swid]['coverage'] = 'NA'
            FPR_info[file_swid]['coverage_dedup'] = 'NA'
            FPR_info[file_swid]['on_target'] = 'NA'                
            FPR_info[file_swid]['percent_duplicate'] = 'NA'




def find_fastq_pairs(files, platform):
    '''
    (dict, str) -> list
    
    Returns a list of 2-item lists with dictionary about file info of paired fastqs.
    Pre-condition: All reads are paired-reads and it exists 2 fastqs for read 1 and read 2
    
    Parameters
    ----------
    - files (dict): Dictionary with file information extracted from FPR 
    '''
     
    # make a sorted list of file_swids
    file_swids = sorted(list(files.keys()))
    
    L = []
    
    for i in range(len(file_swids) -1):
        for j in range(i+1, len(file_swids)):
            if files[file_swids[i]]['run'][0] == files[file_swids[j]]['run'][0] \
                and files[file_swids[i]]['platform'] == files[file_swids[j]]['platform'] == platform \
                and files[file_swids[i]]['barcode'][0] == files[file_swids[j]]['barcode'][0] \
                and files[file_swids[i]]['library'][0] == files[file_swids[j]]['library'][0] \
                and files[file_swids[i]]['sample_name'] == files[file_swids[j]]['sample_name'] \
                and files[file_swids[i]]['external_name'] == files[file_swids[j]]['external_name'] \
                and files[file_swids[i]]['library_source'][0] == files[file_swids[j]]['library_source'][0] \
                and files[file_swids[i]]['lane'][0] == files[file_swids[j]]['lane'][0] \
                and files[file_swids[i]]['tissue_origin'][0] == files[file_swids[j]]['tissue_origin'][0] \
                and files[file_swids[i]]['tissue_type'][0] == files[file_swids[j]]['tissue_type'][0] \
                and files[file_swids[i]]['sample_id'][0] == files[file_swids[j]]['sample_id'][0] \
                and files[file_swids[i]]['limskey'][0] == files[file_swids[j]]['limskey'][0]:
                    assert files[file_swids[i]]['read_count'] == files[file_swids[j]]['read_count']
                    assert files[file_swids[i]]['coverage'] == files[file_swids[j]]['coverage']
                    assert files[file_swids[i]]['coverage_dedup'] == files[file_swids[j]]['coverage_dedup']
                    assert files[file_swids[i]]['on_target'] == files[file_swids[j]]['on_target']
                    assert files[file_swids[i]]['percent_duplicate'] == files[file_swids[j]]['percent_duplicate']
                    L.append([files[file_swids[i]], files[file_swids[j]]])
    assert len(L) == len(files) / 2
    for i in L:
        assert i[0]['file_path'] != i[1]['file_path']
        if 'R1' in i[0]['file_name']:
            assert 'R2' in i[1]['file_name']
            r1, r2 = 0, 1
            
        elif 'R2' in i[0]['file_name']:
            assert 'R1' in i[1]['file_name']
            r1, r2 = 1, 0
        assert i[r1]['file_name'][:i[r1]['file_name'].index('R1')] == i[r2]['file_name'][:i[r2]['file_name'].index('R2')]
        
    return L



def get_run_level_metrics(files, platform):
    '''
    (dict, str) -> (list, list, list, list, list)
    
    Returns a tuple with parallel lists of run-level metrics for a given instrument.
    Pre-condition: All reads are paired-reads and it exists 2 fastqs for read 1 and read 2
    
    Parameters
    ----------
    - files (dict): Dictionary with file information extracted from FPR 
    '''
     
    # find fastq pairs
    L = find_fastq_pairs(files, platform)
    
    reads, coverage, coverage_dedup, on_target, percent_duplicate = [], [], [], [], []    
    
    for i in L:
        reads.append(i[0]['read_count'])
        coverage.append(i[0]['coverage'])
        coverage_dedup.append(i[0]['coverage_dedup'])                
        on_target.append(i[0]['on_target'])
        percent_duplicate.append(i[0]['percent_duplicate'])
    
    return reads, coverage, coverage_dedup, on_target, percent_duplicate                
           

def clean_up_metrics(reads, coverage, on_target, percent_duplicate):
    '''
    (list, list, list, list) -> (list, list, list, list)
    
    Returns a tuple of lists with metrics without missing values NA, keeping the original order of the lists
    
    Parameters
    ----------
    - reads (list): List of read counts
    - coverage (list): List of coverage
    - on_target (list): List of on_target rate 
    - percent_duplicate (list): List of percent_duplicate 
    '''
    
    a = [reads, coverage, on_target, percent_duplicate]
    while any(list(map(lambda x: 'NA' in x, a))):
        if 'NA' in reads:
            pos = reads.index('NA')
        elif 'NA' in coverage:
            pos = coverage.index('NA')
        elif 'NA' in on_target:
            pos = on_target.index('NA')
        elif 'NA' in percent_duplicate:
            pos = percent_duplicate.index('NA')
        for i in range(len(a)):
            if a[i]:
                del a[i][pos]
        a = [reads, coverage, on_target, percent_duplicate]
    
    return reads, coverage, on_target, percent_duplicate
    

    
def sort_metrics(reads, coverage, on_target, percent_duplicate):
    '''
    (list, list, list, list) -> (list, list, list, list)
    
    Returns a tuple of lists of metrics preserving the order among lists and 
    and sorted according to the order of reads
    Pre-condition: There is no missing values and all lists have the same length
        
    Parameters
    ----------
    - reads (list): List of read counts
    - coverage (list): List of coverage
    - on_target (list): List of on_target rate 
    - percent_duplicate (list): List of percent_duplicate 
    '''
    
    # make a list with inner lists containing each the ith value of each metric
    a = list(zip(reads, coverage, on_target, percent_duplicate))
    # sort the inner lists according to the read count, it the first value of each inner list
    a.sort(key=lambda x : x[0])
    # unpack the sorted list to get back the list of each metric, order-preserved and read-counts sorted
    if a:
        reads, coverage, on_target, percent_duplicate = list(zip(*a))
    return reads, coverage, on_target, percent_duplicate

    

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



def generate_figures(files, project, working_dir, height=16, width=13):
    '''
    (dict, str, str, int, int)
    
    Generate figures for each instrument with metrics from FPR and QC-etl and returns
    a dictionary with the figure path for each instrument
        
    Parameters
    ----------
    - files (dict): Dictionary with file info extracted from FPR and with QC info extracted from qc-etl
    - project (str): Name of project
    - working_dir (str): Path to the folder where figure files are written
    - height (int): Height of the figure
    - width (int): Width of the figure
    '''
    
    # make a list of instruments
    instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
    # track the figure file names for each instrument
    figure_files = {}
    for platform in instruments:
        # make lists with metrics for each instrument 
        reads, coverage, coverage_dedup, on_target, percent_duplicate = get_run_level_metrics(files, platform)
        # remove undefined metric values
        reads, coverage, on_target, percent_duplicate = clean_up_metrics(reads, coverage, on_target, percent_duplicate)
        # sort metrics according to read counts
        reads, coverage, on_target, percent_duplicate = sort_metrics(reads, coverage, on_target, percent_duplicate)
        
        # get the outputfile
        current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
        outputfile = os.path.join(working_dir, '{0}.{1}.{2}.QC_plots.png'.format(project, platform, current_time))
    
        data = [reads, coverage, on_target, percent_duplicate]
        
        if reads and coverage:
            figure = plt.figure()
            figure.set_size_inches(width, height)
            # make a list of with X axis labels to determine which subplot should display the Samples label
            x_labels = get_x_axis_labels(data)
            # determine how many subplots are expected
            subplots = count_subplots(data) 
            # determine the position of each subplot
            subplot_pos = get_subplot_position(data)
                        
            # plot read counts and coverage 
            ax1 = create_ax(subplots, 1, subplot_pos[0], figure, reads, 'Read counts', '#00CD6C', title = platform, XLabel = x_labels[0])
            ax2 = create_ax(subplots, 1, subplot_pos[1], figure, coverage, 'Coverage', '#AF58BA', title = None, XLabel = x_labels[1])
            # check if other metrics defined
            if on_target:
                ax3 = create_ax(subplots, 1, subplot_pos[2], figure, on_target, 'On target', '#FFC61E', title = None, XLabel = x_labels[2])
            if percent_duplicate:
                ax4 = create_ax(subplots, 1, subplot_pos[3], figure, percent_duplicate, 'Percent duplicate', '#009ADE', title = None, XLabel = x_labels[3])
                
            # make sure axes do not overlap
            plt.tight_layout(pad = 2.5)
            # write figure to file  
            figure.savefig(outputfile, bbox_inches = 'tight')
            plt.close()
                
            assert platform not in figure_files
            figure_files[platform] = outputfile
    return figure_files





def count_samples_with_missing_values(files):
    '''
    (dict) -> dict
    
    Returns a dictionary with the number of samples with missing metric values for each instrument
    
    Paraneters
    ----------
    - files (dict): Dictionary with file info extracted from FPR and QC metrics extracted from qc-etl
    '''
    
    # count the number of samples with missing values for each instrument
    D = {}
    # make a list of instruments
    instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
    for platform in instruments:
        # make lists with metrics for each instrument 
        reads, coverage, coverage_dedup, on_target, percent_duplicate = get_run_level_metrics(files, platform)
        
        # make a list of indices for whoch the value is NA
        # ie. the indices correspond to different samples, but the same index among lists is the same sample
        I = []
        for i in [reads, coverage, on_target, percent_duplicate]:
            for j in range(len(i)):
                if i[j] == 'NA':
                    I.append(j)
        if I:
            D[platform] = len(list(set(I)))
    return D   


def group_sample_metrics(files, table, add_time_points=None):
    '''
    (dict, str, bool) -> dict 
    
    Returns the header 
    
    
    
    
    
    '''
    
    
    # make a list of instruments
    instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
    
    # record sample metrics for each instrument
    D = {}
    
    # find pairs of fastqs
    for platform in instruments:
        pairs = find_fastq_pairs(files, platform)
        for i in pairs:
            # add comma to reads for readability
            case = i[0]['sample_name']
            external_name = i[0]['external_name']
            sample = i[0]['sample_id'][0]
            library = i[0]['library'][0]
            library_name = '{0} ({1})'.format(i[0]['library'][0], i[0]['time_point'])  
            library_source = i[0]['library_source'][0]
            tissue_origin = i[0]['tissue_origin'][0]
            tissue_type = i[0]['tissue_type'][0]
            reads = '{:,}'.format(i[0]['read_count'])
            on_target = i[0]['on_target']
            coverage = i[0]['coverage']
            duplicate = i[0]['percent_duplicate']
            sequencing_run = '{0} {1}'.format(i[0]['run_id'][0], i[0]['barcode'])
            
            if table == 'sample_identifiers':
                L = [external_name, case, sample, library, library_source, tissue_origin, tissue_type]
                # add time point if selected
                if add_time_points:
                    L[3] = library_name
            elif table == 'qc_metrics':
                L = [case, library, sequencing_run, reads, coverage, on_target, duplicate]
                           
            if platform not in D:
                D[platform] = []
            if L not in D[platform]:
                D[platform].append(L)
            
    return D            



def get_appendix_identifiers(files):
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
    tissue_types = {'X': 'Xenograft derived from some tumour. Note: may not necessarily be a mouse xenograft',
                    'U': 'Unspecified', 'T': 'Unclassifed tumour', 'S': 'Serum from blood where clotting proteins have been removed',
                    'R': 'Reference or non-tumour, non-diseased tissue sample. Typically used as a donor-specific comparison to a diseased tissue, usually a cancer',
                    'P': 'Primary tumour', 'O': 'Organoid', 'n': 'Unknown', 'M': 'Metastatic tumour',
                    'F': 'Fibroblast cells', 'E': 'Endothelial cells', 'C': 'Cell line derived from a tumour',
                    'B': 'Benign tumour', 'A': 'Cells taken from Ascites fluid'}
    tissue_origin = {'Ab': 'Abdomen', 'Ad': 'Adipose', 'Ae': 'Adnexa', 'Ag': 'Adrenal', 'An': 'Anus',
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
    library_design = {'WT': 'Whole Transcriptome', 'WG': 'Whole Genome', 'TS': 'Targeted Sequencing',
                      'TR': 'Total RNA', 'SW': 'Shallow Whole Genome', 'SM': 'smRNA', 'SC': 'Single Cell',
                      'NN': 'Unknown', 'MR': 'mRNA', 'EX': 'Exome', 'CT': 'ctDNA', 'CM': 'cfMEDIP',
                      'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq'}
        
    for file_swid in files:
        D['Library Type'].append(files[file_swid]['library_source'][0])
        D['Tissue Type'].append(files[file_swid]['tissue_type'][0])
        D['Tissue Origin'].append(files[file_swid]['tissue_origin'][0])
    
    for i in ['Library Type', 'Tissue Type', 'Tissue Origin']:
        D[i] = sorted(list(set(D[i])))                
    for i in range(len(D['Library Type'])):
        D['Library Type'][i] = '{0}: {1}'.format(D[i], library_design[D[i]])
    for i in range(len(D['Tissue Type'])):
        D['Tissue Type'][i] = '{0}: {1}'.format(D[i], tissue_types[D[i]])
    for i in range(len(D['Tissue Origin'])):
        D['Tissue Origin'][i] = '{0}: {1}'.format(D[i], tissue_origin[D[i]])
        
    return D

       


def write_batch_report(args):
    '''
    (str, str, str, str, str, str, str, list, str | None)

    Write a PDF report with QC metrics and released fastqs for a given project

    - project (str): Project name as it appears in File Provenance Report
    - working-dir (str): Path to the directory with project directories and links to fastqs 
    - project_name (str): Project name used to create the project directory in gsi space
    - project_code (str): Project code from MISO
    - bamqc_table (str): Path to the bamqc table of qc-etl
    - run_directories (list): List of directories with links to fastqs
    - provenance (str): Path to File Provenance Report.
    - level (str): Simgle release or cumulative project level report. Values: single or cumulative 
    - prefix (str | None): Use of prefix assumes that file paths in File Provenance Report are relative paths.
                           Prefix is added to the relative path in FPR to determine the full file path.
    '''
    
    # check that runs are specified for single data release report
    if args.run_directories is None:
        sys.exit('Please provide a list of run folders')
    # get the project directory with release run folders
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # get the records for the project of interest
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
    records = get_FPR_records(args.project, provenance)
    print('Information was extracted from FPR {0}'.format(provenance))
    # collect relevant information from File Provenance Report about fastqs for project 
    files = parse_fpr_records(provenance, args.project, ['bcl2fastq'], args.prefix)
    
    # keep only info about released fastqs
    # make a list of full paths to the released fastqs resolving the links in the run directories
    released_files = list_files_release_folder(args.run_directories)
    to_remove = [file_swid for file_swid in files if os.path.realpath(files[file_swid]['file_path']) not in released_files]
    for file_swid in to_remove:
        del files[file_swid]
    
    # add time points    
    add_time_points(args.sample_provenance, files)
        
    # count the number of released fastq pairs for each run and instrument
    fastq_counts = count_released_fastqs_by_instrument(files, 'read1')
    all_released_files = sum([fastq_counts[instrument][run] for instrument in fastq_counts for run in fastq_counts[instrument]])
    
    


    
    # collect information from bamqc table
    bamqc_info = extract_bamqc_data(args.bamqc_db)
    # update FPR info with QC info from bamqc table
    map_bamqc_info_to_fpr(files, bamqc_info)
    
    # generate plots for each instrument and keep track of figure files
    figure_files = generate_figures(files, args.project, working_dir)
    
    # write md5sums to separate file
    current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    md5sum_file = os.path.join(working_dir, '{0}.batch.release.{1}.md5'.format(args.project, current_time))
    write_md5sum(files, md5sum_file)

    # write report
    
    # get the report template
    environment = Environment(loader=FileSystemLoader("templates/"))
    template = environment.get_template("batch_report_template.html")
    
    # make a dict with project information
    projects = [{'acronym': args.project_name, 'name': args.project_full_name, 'date': time.strftime('%Y-%m-%d', time.localtime(time.time()))}]

    # count the number of samples with missing metric values
    samples_missing_metrics = count_samples_with_missing_values(files)

    # group metrics by pairs of files
    header_identifiers =  ['Sample', 'Case', 'Group ID', 'Library ID', 'S', 'O', 'T']
    if args.timepoints:
        header_identifiers[3] = 'Library ID (time point)'
    
    sample_identifiers = group_sample_metrics(files, 'sample_identifiers', add_time_points=args.timepoints)
    appendix_identifiers = get_appendix_identifiers(files)
    
    
    
    
    header_metrics = ['Case', 'Library', 'Sequencing Run', 'Reads', 'Coverage', 'On target', 'Duplicate (%)']            
    
    
    
    # fill in template
    context = {'projects' : projects,
               'file_count': all_released_files,
               'fastq_counts': fastq_counts,
               'figure_files': figure_files,
               'samples_missing_metrics': samples_missing_metrics,
               'header_identifiers': header_identifiers,
               'sample_identifiers': sample_identifiers,
               'appendix_identifiers': appendix_identifiers}
    
    
    ### check function for figure formatting
    ## generate_figure_table(figure_files[instrument], factor)
    
    ### generate_table for sample and metric tables
    
    # render template html 
    content = template.render(context)

    filename = "data_release_report_test.html"
    with open(filename, mode="w", encoding="utf-8") as message:
        message.write(content)
        print(f"... wrote {filename}")

    
        
    
        
    
    
    
    
    
    
    # # add appendix with library design, tissue origin and type
    # library_design, tissue_type, tissue_origin = list_library_tissue_codes(sample_metrics, args.level)
    # L = ['{0}: {1}'.format(i, get_library_design()[i]) for i in library_design]
    # T = ['{0}: {1}'.format(i, get_tissue_types()[i]) for i in tissue_type]
    # O = ['{0}: {1}'.format(i, get_tissue_origin()[i]) for i in tissue_origin]
    
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Appendix Table 1</p>')
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Library type:</span> {0}.<li/></ul>'.format(', '.join(L)))
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Tissue Type:</span> {0}.<li/></ul>'.format(', '.join(T)))
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Tissue Origin:</span> {0}.<li/></ul>'.format(', '.join(O)))
    # # add page break between plots and tables
    # Text.append('<div style="page-break-after: always;"></div>')
                
    # # add QC metrics table
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 2. QC metrics</p>')
    # Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Table 2 provides QC metrics about the raw sequences of each sample.</span></p>')
        
    # if args.level == 'single':
    #     header = ['Case', 'Library', 'Sequencing Run', 'Reads', 'Coverage', 'On target', 'Duplicate (%)']       
    #     column_size = {'Case': '10%', 'Library': '22%', 'Sequencing Run': '31%', 'Reads': '9%', 'Coverage': '9%', 'On target': '8%', 'Duplicate (%)': '11%'}
    #     Text.append(generate_table(sample_metrics, header, column_size))
    # # add space
    # Text.append('<br />')    
    
    # # add QC metrics appendix
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Appendix Table 2</p>')
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style= "font-weight: bold">Raw Coverage:</span> an estimate of the mean depth of coverage in the target space = total bases on target / size of the target space<li/></ul>')
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style= "font-weight: bold">On Target Rate:</span> percentage of reads that overlap the target space by at least one base = reads on target/total reads.<li/></ul>')
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Percent duplicate:</span> Percent of duplicate reads estimated by Picard MarkDuplicates.<li/></ul>')
    
    # # add md5sums
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 3. List of md5sums</p>')
    # Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">A list of md5sums is available in the accompanying file: <span style="color: black; font-style: italic">{0}</span></p>'.format(os.path.basename(md5_file)))
        
    # # convert to html
    # renderer = mistune.Markdown()
    # Text = '\n'.join(Text)
    # html_str = renderer(Text)
    
    # # convert html to pdf    
    # report_name = '{0}_run_level_data_release_report.{1}.pdf' if args.level == 'single' else '{0}_cumulative_data_release_report.{1}.pdf'
    # report_file = os.path.join(project_dir, report_name.format(args.project_name, current_date))
    # newfile = open(report_file, "wb")
    # pisa.CreatePDF(html_str, newfile)
    # newfile.close()

    # # remove figure files from disk
    # for i in figure_files:
    #     os.remove(figure_files[i])





        
def write_report(args):
    '''
    (str, str, str, str, str, str, str, list, str | None)

    Write a PDF report with QC metrics and released fastqs for a given project

    - project (str): Project name as it appears in File Provenance Report
    - working-dir (str): Path to the directory with project directories and links to fastqs 
    - project_name (str): Project name used to create the project directory in gsi space
    - project_code (str): Project code from MISO
    - bamqc_table (str): Path to the bamqc table of qc-etl
    - run_directories (list): List of directories with links to fastqs
    - provenance (str): Path to File Provenance Report.
    - level (str): Simgle release or cumulative project level report. Values: single or cumulative 
    - prefix (str | None): Use of prefix assumes that file paths in File Provenance Report are relative paths.
                           Prefix is added to the relative path in FPR to determine the full file path.
    '''
    
    # check that runs are specified for single data release report
    if args.level == 'single' and args.run_directories is None:
        sys.exit('Please provide a list of run folders')
    # emit warning if time points are used with cumulative report
    if args.level == 'cumulative' and args.timepoints:
        print('Option timepoint has no effect. Time points are only added to batch level reports')
        
    # get the project directory with release run folders
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # get the records for the project of interest
    # dereference link to FPR
    provenance = os.path.realpath(args.provenance)
    records = get_FPR_records(args.project, provenance)
    print('Information was extracted from FPR {0}'.format(provenance))
    # collect relevant information from File Provenance Report about fastqs for project 
    files = parse_fpr_records(provenance, args.project, ['bcl2fastq'], args.prefix)
    
    # keep only info about released fastqs
    if args.level == 'single':
        # make a list of full paths to the released fastqs resolving the links in the run directories
        released_files = list_files_release_folder(args.run_directories)
    #elif args.level == 'cumulative':
    #     released_files = list_released_fastqs_project(args.api, FPR_info)
    to_remove = [file_swid for file_swid in files if os.path.realpath(files[file_swid]['file_path']) not in released_files]
    for file_swid in to_remove:
        del files[file_swid]
        
    # count the number of released fastq pairs for each run and instrument
    fastq_counts = count_released_fastqs_by_instrument(files, 'read1')
    
    # collect information from bamqc table
    if args.level == 'single':
       bamqc_info = extract_bamqc_data(args.bamqc_db)
       # update FPR info with QC info from bamqc table
       map_bamqc_info_to_fpr(files, bamqc_info)
       
           
    
    
    
    #     # re-organize metrics per sample and instrument
    #     sample_metrics = get_run_level_sample_metrics(FPR_info)
    #     # add time points
    #     if args.timepoints:
    #         # update sample metrics by adding time points to each sample
    #         add_time_points(args.sample_provenance, sample_metrics)
            
    # elif args.level == 'cumulative':
    #     bamqc_info = parse_merged_bamqc(args.bamqc_table, args.project)   
    #     # update FPR info with QC info from bamqc merged table
    #     map_merged_bamqc_info_to_fpr(FPR_info, bamqc_info)
    #     # re-organize metrics per sample and instrument
    #     sample_metrics = get_cumulative_level_sample_metrics(FPR_info)
            
    # # generate figure files
    # if args.level == 'single':
    #     # plot read count, coverage, on target and duplicate rate
    #     figure_files = generate_figures(project_dir, args.level, args.project_name, sample_metrics, 'reads', 'coverage', 'Read counts', 'Coverage', '#00CD6C', '#AF58BA', 'Samples', 13, 16, 'duplicate (%)', 'on_target', '#009ADE', '#FFC61E', 'Percent duplicate', 'On target', keep_on_target=False)
    # elif args.level == 'cumulative':
    #     # plot read count and coverage
    #     figure_files = generate_figures(project_dir, args.level, args.project_name, sample_metrics, 'reads', 'coverage', 'Read counts', 'Coverage', '#00CD6C', '#AF58BA', 'Samples', 13, 8)
        
    # # get current date (year-month-day)
    # current_date = datetime.today().strftime('%Y-%m-%d')
    
    # # write md5sums to separate text file
    # if args.level == 'single':
    #     md5_file = os.path.join(project_dir, '{0}_fastqs_release_{1}.md5'.format(args.project_name, current_date))
    #     # make a list of file paths, md5sum inner lists
    #     write_md5sums(md5_file, [[FPR_info[file]['md5sum'], FPR_info[file]['filename']] for file in FPR_info])
            
    # # make a list to store report
    # Text = []
    # # add title and logo
    # logo = str(get_logo())
    # height, width = resize_image(logo, 0.085)
    # Text.append(generate_header_table(logo, width, height, args.level))
    # Text.append('<br />' * 3)
    # # add information about project and contact personn
    # Text.append(generate_project_table(args.project_name, args.project_full_name, current_date))
    # Text.append('<br />')           
    
    # # add description
    # if args.level == 'single':
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">This report provides detailed sample information and QC metrics about newly released raw sequences.</p>')
    # elif args.level == 'cumulative':
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">This report provides detailed sample information and QC metrics about all released raw sequences to date.</p>')
    # Text.append('<br />' * 2)
          
    # # list the file count            
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">1. File Count</p>')
    # if args.level == 'single':
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">This release includes {0} pairs of fastq files. File count is broken down by instrument and run as follow.</p>'.format(count_all_files(fastq_counts)))
    # elif args.level == 'cumulative':
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">{0} pairs of fastq files have been released. File count is broken down by instrument and run as follow.</p>'.format(count_all_files(fastq_counts)))
    # Text.append(generate_file_count_table(fastq_counts, ['Platform', 'Run', 'Paired fastq files'], {'Platform': '25%', 'Run': '30%', 'Paired fastq files': '25%'}))
    # Text.append('<br />')           
    
    # # add page break between plots and tables
    # Text.append('<div style="page-break-after: always;"></div>')
        
    # # add QC plots
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">2. QC plots</p>')
    # if args.level == 'single':
    #     # count samples with missing values
    #     discarded_samples = count_samples_with_missing_values(sample_metrics, 'reads', 'coverage', args.level, 'duplicate (%)', 'on_target')
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">QC plots are reported by instrument. Lines are the median of each metric. <span style="font-style: italic">Read count</span> is plotted by ascending order. Other metrics are plotted according to the order of <span style="font-style: italic">read counts</span></p>')
    # elif args.level == 'cumulative':
    #     # count samples with missing values
    #     discarded_samples = count_samples_with_missing_values(sample_metrics, 'reads', 'coverage', args.level)
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">QC plots are reported by instrument. Lines are the median of each metric. <span style="font-style: italic">Read count</span> is plotted by ascending order. <span style="font-style: italic">Coverage</span> is plotted according to the order of <span style="font-style: italic">read counts</span></p>')
    
    # if sum(discarded_samples.values()):
    #     S = ['{0}: {1}'.format(instrument, discarded_samples[instrument]) for instrument in sorted(list(discarded_samples.keys()))]
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Some samples could not be plotted because of missing QC values. Missing QC values appear as NA in the QC tables below. The number of discarded samples for each instrument is: {0}</span></p>'.format(', '.join(S)))
    # Text.append('<br />')
    
    
    # if args.level == 'single':
    #     #factor = 0.3
    #     factor = 1
        
    # elif args.level == 'cumulative':
    #     factor = 1.3
    
    # for instrument in sorted(list(sample_metrics.keys())):
    #     # check that figures exist for instrument
    #     if instrument in figure_files:
    #         Text.append(generate_figure_table(figure_files[instrument], factor))
    #         Text.append('<br />')
                
    # # add page break between plots and tables
    # Text.append('<div style="page-break-after: always;"></div>')
        
    # # add table with sample Ids
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 1. Sample identifiers</p>')
    # Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Table 1 provides information about the sequenced samples.</span></p>')
    # Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">S: Library type, T: Tissue type, O: Tissue origin.</span></p>')
       
    # if args.level == 'single':
    #     if args.timepoints:
    #         header = ['External ID', 'Case', 'Group ID', 'Library ID (time point)', 'S', 'O', 'T']
    #         column_size = {'Case': '10%', 'Group ID': '36%', 'Library ID (time point)': '22%', 'S': '5%', 'O': '5%', 'T': '5%', 'External ID': '17%'}
    #     else:
    #         header = ['External ID', 'Case', 'Group ID', 'Library ID', 'S', 'O', 'T']
    #         column_size = {'Case': '10%', 'Group ID': '36%', 'Library ID': '22%', 'S': '5%', 'O': '5%', 'T': '5%', 'External ID': '17%'}
    #     Text.append(generate_table(sample_metrics, header, column_size))            
    # elif args.level == 'cumulative':
    #     header = ['External ID', 'Case', 'Sample', 'Library ID', 'S', 'O', 'T', 'Run']
    #     column_size = {'Case': '9%', 'Library ID': '21%', 'S': '3%', 'O': '3%', 'T': '3%', 'Run': '32%', 'Sample': '18%', 'External ID': '11%'}
    #     Text.append(generate_cumulative_table(sample_metrics, header, column_size))
        
    # Text.append('<br />')
    
    # # add appendix with library design, tissue origin and type
    # library_design, tissue_type, tissue_origin = list_library_tissue_codes(sample_metrics, args.level)
    # L = ['{0}: {1}'.format(i, get_library_design()[i]) for i in library_design]
    # T = ['{0}: {1}'.format(i, get_tissue_types()[i]) for i in tissue_type]
    # O = ['{0}: {1}'.format(i, get_tissue_origin()[i]) for i in tissue_origin]
    
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Appendix Table 1</p>')
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Library type:</span> {0}.<li/></ul>'.format(', '.join(L)))
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Tissue Type:</span> {0}.<li/></ul>'.format(', '.join(T)))
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Tissue Origin:</span> {0}.<li/></ul>'.format(', '.join(O)))
    # # add page break between plots and tables
    # Text.append('<div style="page-break-after: always;"></div>')
                
    # # add QC metrics table
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 2. QC metrics</p>')
    # Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Table 2 provides QC metrics about the raw sequences of each sample.</span></p>')
        
    # if args.level == 'single':
    #     header = ['Case', 'Library', 'Sequencing Run', 'Reads', 'Coverage', 'On target', 'Duplicate (%)']       
    #     column_size = {'Case': '10%', 'Library': '22%', 'Sequencing Run': '31%', 'Reads': '9%', 'Coverage': '9%', 'On target': '8%', 'Duplicate (%)': '11%'}
    #     Text.append(generate_table(sample_metrics, header, column_size))
    # elif args.level == 'cumulative':
    #     header = ['Case', 'Sample', 'Libraries', 'Runs', 'Reads', 'Coverage', 'Coverage_dedup']
    #     column_size = {'Case': '15%', 'Sample': '36%', 'Libraries': '7%', 'Runs': '7%', 'Reads': '10%', 'Coverage': '10%', 'Coverage_dedup': '15%'}
    #     Text.append(generate_cumulative_table(sample_metrics, header, column_size, table_type='metrics'))        
    # # add space
    # Text.append('<br />')    
    
    # # add QC metrics appendix
    # Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:normal">Appendix Table 2</p>')
    # Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style= "font-weight: bold">Raw Coverage:</span> an estimate of the mean depth of coverage in the target space = total bases on target / size of the target space<li/></ul>')
    # if args.level == 'single':
    #     Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style= "font-weight: bold">On Target Rate:</span> percentage of reads that overlap the target space by at least one base = reads on target/total reads.<li/></ul>')
    #     Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Percent duplicate:</span> Percent of duplicate reads estimated by Picard MarkDuplicates.<li/></ul>')
    # elif args.level == 'cumulative':
    #     Text.append('<ul style="list-style-type: circle; text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-style:normal; font-weight:normal"><li><span style="font-weight: bold">Coverage Deduplicated:</span> an estimate of the mean depth of coverage after removal of marked pcr duplicates. = raw coverage / (1 – percent_duplicates).<li/></ul>')
    
    # # add md5sums
    # if args.level == 'single':
    #     Text.append('<p style="text-align: left; color: black; font-size:14px; font-family: Arial, Verdana, sans-serif; font-weight:bold">Table 3. List of md5sums</p>')
    #     Text.append('<p style="text-align: left; color: black; font-size:12px; font-family: Arial, Verdana, sans-serif; font-weight:normal">A list of md5sums is available in the accompanying file: <span style="color: black; font-style: italic">{0}</span></p>'.format(os.path.basename(md5_file)))
        
    # # convert to html
    # renderer = mistune.Markdown()
    # Text = '\n'.join(Text)
    # html_str = renderer(Text)
    
    # # convert html to pdf    
    # report_name = '{0}_run_level_data_release_report.{1}.pdf' if args.level == 'single' else '{0}_cumulative_data_release_report.{1}.pdf'
    # report_file = os.path.join(project_dir, report_name.format(args.project_name, current_date))
    # newfile = open(report_file, "wb")
    # pisa.CreatePDF(html_str, newfile)
    # newfile.close()

    # # remove figure files from disk
    # for i in figure_files:
    #     os.remove(figure_files[i])










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
    - project (str): Name of Project in FPR
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
        records = get_FPR_records(args.project, provenance)
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
    m_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report. Used to parse the FPR by project. Files are further filtered by run is runs parameter if provided, or all files for the project and workflow are used')
    m_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run IDs. Include one or more run Id separated by white space. Other runs are ignored if provided')
    m_parser.add_argument('--exclude_miseq', dest='nomiseq', action='store_true', help='Exclude MiSeq runs if activated')
    m_parser.add_argument('-e', '--exclude', dest='exclude', help='File with libraries tagged for non-release. The first column is always the library. The optional second column is the run id')
    m_parser.add_argument('-f', '--files', dest='release_files', help='File with file names to be released')
    m_parser.add_argument('--time_points', dest='timepoints', action='store_true', help='Add time points to sample map if option is used. By default, time points are not added.')
    m_parser.add_argument('--panel', dest='add_panel', action='store_true', help='Add panel to sample if option is used. By default, panel is not added.')
    m_parser.add_argument('-spr', '--sample_provenance', dest='sample_provenance', default='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance', help='Path to File Provenance Report. Default is http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance')
    m_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    m_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    m_parser.set_defaults(func=map_external_ids)

    # mark files in nabu 
    n_parser = subparsers.add_parser('mark', help="Mark released or withheld files in Nabu")
    n_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report. Used to parse the FPR by project', required = True)
    n_parser.add_argument('-u', '--user', dest='user', help='User name to appear in Nabu for each released or whitheld file', required=True)
    n_parser.add_argument('-s', '--status', dest='status', choices = ['fail', 'pass'], help='Mark files accordingly when released or withheld', required = True)
    n_parser.add_argument('-d', '--directory', dest='directory', help='Directory with links organized by project and run in gsi space', required=True)
    n_parser.add_argument('-c', '--comment', dest='comment', help='Comment to be added to the released file')
    n_parser.add_argument('-a', '--api', dest='api', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
    n_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    n_parser.set_defaults(func=mark_files_nabu)
    
    
    # write a report
    r_parser = subparsers.add_parser('report', help="Write a PDF report for released FASTQs")
    r_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report', required=True)
    r_parser.add_argument('-p', '--parents', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    r_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space', required = True)
    r_parser.add_argument('-fn', '--full_name', dest='project_full_name', help='Full name of the project', required = True)
    r_parser.add_argument('-r', '--runs', dest='run_directories', nargs='*', help='List of directories with released fastqs')
    r_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    r_parser.add_argument('-a', '--api', dest='api', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
    r_parser.add_argument('--time_points', dest='timepoints', action='store_true', help='Add time points to Identifiers Table if option is used. By default, time points are not added.')
    r_parser.add_argument('-spr', '--sample_provenance', dest='sample_provenance', default='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance', help='Path to File Provenance Report. Default is http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance')
    r_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    r_parser.add_argument('-bq', '--bamqc', dest='bamqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest', help='Path to the bamqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest')
    r_parser.set_defaults(func=write_batch_report)
    
    # # write a report
    # r_parser = subparsers.add_parser('report', help="Write a PDF report for released FASTQs")
    # r_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report', required=True)
    # r_parser.add_argument('-p', '--parents', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    # r_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space', required = True)
    # r_parser.add_argument('-fn', '--full_name', dest='project_full_name', help='Full name of the project', required = True)
    # r_parser.add_argument('-r', '--runs', dest='run_directories', nargs='*', help='List of directories with released fastqs')
    # r_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    # r_parser.add_argument('-a', '--api', dest='api', default='http://gsi-dcc.oicr.on.ca:3000', help='URL of the Nabu API. Default is http://gsi-dcc.oicr.on.ca:3000')
    # r_parser.add_argument('-l', '--level', dest='level', choices=['single', 'cumulative'], help='Generates a single release report or a cumulative project report', required = True)
    # r_parser.add_argument('--time_points', dest='timepoints', action='store_true', help='Add time points to Identifiers Table if option is used. By default, time points are not added.')
    # r_parser.add_argument('-spr', '--sample_provenance', dest='sample_provenance', default='http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance', help='Path to File Provenance Report. Default is http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance')
    # r_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    # r_parser.add_argument('-bq', '--bamqc', dest='bamqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest', help='Path to the bamqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest')
    # r_parser.set_defaults(func=write_batch_report)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)