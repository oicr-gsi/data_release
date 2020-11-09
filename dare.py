# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 12:27:07 2020

@author: rjovelin
"""


import argparse
import subprocess
import os


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
    


def extract_files(project, runs, workflow, nomiseq, library_aliases, exclude):
    '''
    (str, list, str, bool, dict, list) -> (dict, dict)
  
    Returns a tuple with dictionaries with files extracted from FPR and their corresponding run
    respectively from release and withheld from release
    
        
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
    - exclude (list): List of samples or libraries to exclude from release
    '''
    
    # create a dict {run: [files]}
    D, K = {}, {}
    
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
            if workflow != i[30]:
                continue
        # check if library aliases are provided
        if library_aliases:
            # skip libraries not included in file
            if library not in library_aliases:
                continue
            else: 
                # skip libraries if aliquot ID is provided and incorrect
                if library_aliases[library] and library_aliases[library] != aliquot:
                    continue
        # check if sample or library is excluded 
        if sample_name in exclude or library in exclude:
            if run_id not in K:
                K[run_id] = []
            K[run_id].append(file_path)
        else:
            # record files per run
            if run_id not in D:
                D[run_id] = []
            D[run_id].append(file_path)

    return D, K



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
    (str | None, list | None, str | None, str, bool, str, str, dict) -> None
    
    Parameters
    ----------
    - libraries (str | None): Path to 1 or 2 columns tab-delimited file with library IDs.
                              The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM).
                              The second and optional column is the library aliquot ID (eg. LDI32439).
                              Only the samples with these library aliases are used if provided'
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
    - run_name (str): Specifies the run folder name. Run Id or run.withhold as run folder name of not specified. 
    - keywords (str): Optional run name parameter. Valid option: run_name. Replaces run ID
    - exclude (str | list): File with sample name or libraries to exclude from the release,
                            or a list of sample name or libraries
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

    # extract files from FPR
    runs = args.runs if args.runs else []
    project = args.project if args.project else ''
    files_release, files_withhold = extract_files(project, runs, args.workflow, args.nomiseq, libraries, exclude)
    
    # link files to project dir
    if args.suffix == 'fastqs':
        assert args.workflow == 'bcl2fastq' or args.workflow.lower() == 'casava'
    generate_links(files_release, files_withhold, args.project_name, args.projects_dir, args.suffix, run_name = args.run_name)
        
    
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

    
    
def map_external_ids(args):
    '''
    (str | None, list | None, str | None, str, bool, str, str, str) -> None

    Parameters
    ----------    
    - libraries (str | None): Path to 1 or 2 columns tab-delimited file with library IDs.
                              The first column is always the library alias (TGL17_0009_Ct_T_PE_307_CM).
                              The second and optional column is the library aliquot ID (eg. LDI32439).
                              Only the samples with these library aliases are used if provided'
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

    # extract files from FPR
    runs = args.runs if args.runs else []
    project = args.project if args.project else ''
    files_release, files_withhold = extract_files(project, runs, args.workflow, args.nomiseq, libraries, exclude)
    
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
    l_parser.add_argument('-s', '--suffix', dest='suffix', choices=['fastqs', 'datafiles'], help='Indicates if fastqs or datafiles are released by adding suffix fastqs or datafiles in the directory name', required=True)
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
    m_parser.add_argument('-s', '--suffix', dest='suffix', choices=['fastqs', 'datafiles'], help='Indicates if fastqs or datafiles are released by adding suffix fastqs or datafiles in the directory name', required=True)
    m_parser.set_defaults(func=map_external_ids)

    # mark files in nabu 
    n_parser = subparsers.add_parser('mark', help="Mark released or withheld files in Nabu")
    n_parser.add_argument('-u', '--user', dest='user', help='User name to appear in Nabu for each released or whitheld file', required=True)
    n_parser.add_argument('-rl', '--release', dest='release', choices = ['fail', 'pass'], help='Mark files accordingly when released or withheld', required = True)
    n_parser.add_argument('-d', '--directory', dest='directory', help='Directory with links organized by project and run in gsi space', required=True)
    n_parser.set_defaults(func=mark_files_nabu)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)