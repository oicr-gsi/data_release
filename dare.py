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



#### deprecated functions for cumulative report

# def write_cumulative_report(args):
#     '''
#     (str, str, str, str, str | None, str, str, str, str, str, str, str, str) -> None
        
#     Write a cumulative PDF report with QC metrics for all released fastqs for a given project
  
#     - project (str): Name of project of interest
#     - projects_dir (str): Parent directory containing the project subdirectories with file links
#     - project_name (str): Project name used to create the project directory in gsi space
#     - provenance (str): Path to File Provenance Report
#     - prefix (str | None): Use of prefix assumes that FPR containes relative paths.
#                            Prefix is added to the relative paths in FPR to determine the full file paths
#     - api (str): URL of the Nabu API
#     - merged_bamqc_db (str): Path to the merged bamqc SQLite database
#     - archived_bamqc (str): Path to the archived merged bamqc json file
#     - merged_rnaseqqc_db (str): Path to the merged rnaseq SQLite database
#     - cfmedipqc_db (str): Path to the cfmedip SQLite database
#     - emseqqc_db (str): Path to the emseq SQLite database
#     - project_full_name (str): Full name of the project
#     - ticket (str): Jira data release ticket code
#     - user (str): Name of the GSI personnel generating the report
#     '''
    
#     # get current time
#     current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
        
#     # get the project directory with release run folders
#     working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
       
#     # get the records for the project of interest
#     # dereference link to FPR
#     provenance = os.path.realpath(args.provenance)
#     print('Information was extracted from FPR {0}'.format(provenance))
    
#     # make a dict with project information
#     projects = {'acronym': args.project, 'name': args.project_full_name, 'date': time.strftime('%Y-%m-%d', time.localtime(time.time()))}

#     # get information about the released fastqs
#     # collect relevant information from File Provenance Report about fastqs for project 
#     files = parse_fpr_records(provenance, args.project, ['bcl2fastq'], args.prefix)
    
#     print('files', len(files))
    
#     # get the released files at the project level from nabu
#     released_files = list_released_fastqs_project(args.api, args.project)
    
#     print('released files', len(released_files))
    
#     # resolve links
#     released_files = resolve_links(released_files)
#     # remove files not released
#     to_remove = [file_swid for file_swid in files if os.path.realpath(files[file_swid]['file_path']) not in released_files]
#     for file_swid in to_remove:
#         del files[file_swid]
    
#     print('files', len(files))    
        
#     # count the number of released fastq pairs for each run and instrument
#     fastq_counts = count_released_fastqs_by_library_type_instrument(files)
#     all_released_files = 0 
#     for i in fastq_counts:
#         all_released_files += sum([sum(list(j.values())) for j in fastq_counts[i].values()])
    
#     # sort library types, sequencing platforms, sequencing runs
#     lt = sorted(list(fastq_counts.keys()))
#     sp = sp = {i:sorted(fastq_counts[i].keys()) for i in lt}
#     rn = {}
#     for i in lt:
#         rn[i] = {}
#         for j in sp[i]:
#             rn[i][j] = rn[i][j] = sorted(list(fastq_counts[i][j].keys()))
        
#     # get the identifiers of all released files
#     sample_identifiers = group_cumulative_samples(files)
    
#     appendix_identifiers = get_identifiers_appendix(files, 'cumulative')
#     # define identifier table header
#     header_identifiers = ['OICR Case Id', 'Donor Id', 'OICR Sample Id', 'Sample Description', 'LT', 'TO', 'TT']
    
#     # get information for lane level sequencing
#     sequencing = format_sequencing_table(files)
        
#     # collect information from merged bamqc table
#     bamqc_info = extract_merged_bamqc_data(args.merged_bamqc_db)
#     print('bamqc', len(bamqc_info))
    
    
#     ## NEED TO IDENTIFY COVERAGE AND COVERAGE DEDUPLICATED IN ARCHIVED BAMQC
    
#     # collect information from archived bamqc
#     archived_bamqc_info = extract_archived_merged_bamqc_data(args.archived_bamqc, args.project)
#     # merged bamqc info
#     bamqc_info.update(archived_bamqc_info)
#     print('bamqc + archive', len(bamqc_info))    
    
    
#     # collect information from rnaseq table
#     rnaseqqc_info = extract_merged_rnaseqqc_data(args.merged_rnaseqqc_db)
#     print('rnaseqqc', len(rnaseqqc_info))
#     archived_rnaseq_info = extract_archived_merged_rnaseqqc_data(args.archived_rnaseqqc, args.project)
#     # merge rnaseqqc info
#     rnaseqqc_info.update(archived_rnaseq_info)
#     print('rnaseqqc + archive', len(rnaseqqc_info))    
    
#     # collect information from cfmedip table
#     cfmedipqc_info = extract_cfmedipqc_data(args.cfmedipqc_db)
#     print('cfmediqc', len(cfmedipqc_info))    
    
#     #### NEED TO IDENTIFY RELATIVE CPG FROM ARCHIVED CFMEDIP
    
#     archived_cfmedipqc = extract_archived_cfmedipqc_data(args.archived_cfmedipqc)
#     # merge cfmedipqc
#     cfmedipqc_info = merge_qc_info(cfmedipqc_info, archived_cfmedipqc)
#     print('cfmediqc + archive', len(cfmedipqc_info))    

#     # collect information from emseq cache
#     emseqqc_info = extract_emseqqc_data(args.emseqqc_db)
#     print('emseqqc', len(emseqqc_info))
        
#     ## NEED TO IDENTIFY LAMBDA AND PUC19 FROM ARCHIVED EMSEQQC
    
#     archived_emseqqc = extract_archived_emseqqc_data(args.archived_emseqqc)
#     # merge emseeqc info
#     emseqqc_info = merge_qc_info(emseqqc_info, archived_emseqqc)
#     print('emseqqc + archive', len(emseqqc_info))
    
#     print('Extracted QC metrics from QC-etl caches')
    
#     # group identifiers, metrics by instrument and library type
#     library_metrics = get_metrics_cumulative_report(files, bamqc_info, rnaseqqc_info, cfmedipqc_info, emseqqc_info)
#     # convert read counts to counts of read pairs
#     convert_read_count_to_read_pairs(library_metrics)
#     # format qc metrics for template
#     qc_metrics = format_qc_metrics(library_metrics)
    
#     # list all platforms for each library source
#     platforms = {}
#     for i in qc_metrics:
#         for j in qc_metrics[i]:
#             if j in platforms:
#                 platforms[j].append(i)
#                 platforms[j].sort()
#             else:
#                 platforms[j] = [i]

#     # make a list of library types from data with QC metrics
#     library_sources = sorted(list(platforms.keys()))
       
#     # get the qc metrics subtables
#     qc_subtables = get_qc_metrics_table_names(library_sources, 3)
#     # get the metrics appendix
#     qc_appendices = get_cumulative_metrics_appendix(library_sources)
    
      
#     # generate plots for each instrument and library source and keep track of figure files
#     figure_files = {}
#     Y_axis = {}
#     colors = ['#00CD6C', '#AF58BA', '#FFC61E', '#009ADE']
#     for library_source in platforms:
#         Y_axis[library_source] = get_Y_axis_labels(library_source)
#         for i in range(len(Y_axis[library_source])):
#             if Y_axis[library_source][i] == 'Read pairs':
#                 Y_axis[library_source][i] = 'Total Read Pairs'
#             elif Y_axis[library_source][i] == 'Coverage':
#                 Y_axis[library_source][i] = 'Final Merged Coverage'
#         for instrument in platforms[library_source]:
#             figure = generate_cumulative_figures(working_dir, args.project, library_metrics, instrument, library_source, Y_axis[library_source], colors)
#             if library_source not in figure_files:
#                 figure_files[library_source] = {}
#             figure_files[library_source][instrument] = figure
#     print('Generated QC plots')
        
#     header_metrics = {}
#     for i in library_sources:
#         header_metrics[i] = ['OICR Case Id', 'OICR Sample Id', 'Lane Count'] + Y_axis[i]
         
#     # write report
#     # get the report template
#     template_dir = os.path.join(os.path.dirname(__file__), './templates')
#     environment = Environment(loader = FileSystemLoader(template_dir), autoescape = True)
#     template = environment.get_template("cumulative_report_template.html")
       
    
#     # fill in template
#     context = {'projects' : projects,
#                'file_count': all_released_files,
#                'fastq_counts': fastq_counts,
#                'header_identifiers': header_identifiers,
#                'sample_identifiers': sample_identifiers,
#                'appendix_identifiers': appendix_identifiers,
#                'sequencing': sequencing,
#                'user': args.user,
#                'ticket': os.path.basename(args.ticket),
#                'library_sources': library_sources,
#                'qc_subtables': qc_subtables,
#                'qc_appendices': qc_appendices,
#                'header_metrics': header_metrics,
#                'qc_metrics': qc_metrics,
#                'figure_files': figure_files,
#                'platforms': platforms,
#                'lt': lt,
#                'sp': sp,
#                'rn': rn
#                }
    
#     # render template html 
#     content = template.render(context)

#     # convert html to PDF
#     report_file = os.path.join(working_dir,  '{0}_cumulative_data_release_report.{1}.pdf'.format(args.project, current_time))
#     makepdf(content, report_file)

#     # remove figure files
#     for i in figure_files:
#         for j in figure_files[i]:
#             if os.path.isfile(figure_files[i][j]):
#                 os.remove(figure_files[i][j])



# def generate_cumulative_figures(working_dir, project, qc_metrics, platform, library_type, Y_axis, colors, width=13, height=16):
#     '''
#     (str, str, dict, str, str, list, list, int, int) -> str
        
#     Generate a figure with metrics from FPR and QC-etl for a specific sequencing platform
#     and library type and returns the path to the figure file
        
#     Parameters
#     ----------
#     - working_dir (str): Path to the folder where figure files are written
#     - project (str): Name of project
#     - qc_metrics (dict): Dictionary with QC metrics extracted from QC etl caches
#     - platform (str): Sequencing platform
#     - library_type (str): Type of library
#     - Y_axis (list): List of Y_axis labels   
#     - colors (list): List of colors
#     - width (int): Width of the figure
#     - height (int): Height of the figure
#     '''
    
#     metrics_interest = get_library_metrics(library_type)

#     # make lists of metrics, sort according to read count
#     D = {}
#     for i in qc_metrics[platform][library_type]:
#         for j in metrics_interest:
#             if j not in D:
#                 D[j] = [qc_metrics[platform][library_type][i][j]]
#             else:
#                 D[j].append(qc_metrics[platform][library_type][i][j])
    
#     L = [D[i] for i in metrics_interest]
#     L = sort_metrics(L)
    
#     # get the outputfile
#     current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
#     outputfile = os.path.join(working_dir, '{0}.{1}.{2}.{3}.QC_plots.png'.format(project, platform, library_type, current_time))
    
#     if L[0]:
#         figure = plt.figure()
#         figure.set_size_inches(width, height)
#         # make a list of with X axis labels to determine which subplot should display the Samples label
#         x_labels = get_x_axis_labels(L)
#         # determine how many subplots are expected
#         subplots = count_subplots(L) 
#         # determine the position of each subplot
#         subplot_pos = get_subplot_position(L)
        
#         for i in range(len(L)):
#             # determine title
#             title = platform + ' {0} libraries'.format(library_type) if i == 0 else None
#             # plot data
#             create_ax(subplots, 1, subplot_pos[i], figure, L[i], Y_axis[i], colors[i], title = title, XLabel = x_labels[i])
                
#         # make sure axes do not overlap
#         plt.tight_layout(pad = 2.5)
#         # write figure to file  
#         figure.savefig(outputfile, bbox_inches = 'tight')
#         plt.close()
                
#         return outputfile
#     else:
#         return ''


# def get_cumulative_metrics_appendix(library_sources):
#     '''
#     (list) -> dict
    
#     Returns a dictionry with definitions of columns in the QC metrics table of the cumulative report
    
#     Parameters
#     ----------
#     - library_sources (list): Sorted list of library types 
#     '''

#     # get the QC metrics sub-tables and appendices
#     definitions = metrics_definitions()    
      
#     counter = 1
#     qc_appendix = {'tables': {}, 'metrics': {}}
#     for library_type in library_sources:
#         qc_appendix['tables'][library_type] = 'Appendix Table 2.{0}'.format(counter)
        
#         L = ['OICR Case Id: OICR-generated case identifier',
#              'OICR Sample Id: The OICR generated sample identifier. The sample Id is formed from the following: 1. Case Id, 2. Tissue Origin, 3. Tissue Type, 4. Library Type and 5. User supplied Sample Id',
#              'Lane Count: Number of sequencing lane',
#              'Total Read Pairs: {0}'.format(definitions['Read pairs'])]
        
#         if library_type == 'CM':
#             L.extend([': '.join([i, definitions[i]]) for i in ['Methylation {0}'.format(chr(946)), 'CpG frequency']])
#         elif library_type in ['EX', 'TS']:
#             L.extend([': '.join([i, definitions[i]]) for i in ['Coverage', 'On target']])
#         elif library_type in ['WG', 'PG']:
#             L.extend(['{0}: {1}'.format('Coverage', definitions['Coverage'])])
#         elif library_type == 'WT':
#             L.extend([': '.join([i, definitions[i]]) for i in ['rRNA contamination', 'Coding (%)']])
#         elif library_type in ['MC', 'MG']:
#             L.extend([': '.join([i, definitions[i]]) for i in ['{0} methylation'.format(chr(955)), 'pUC19 methylation', 'Duplication rate']])
#         qc_appendix['metrics'][library_type] = L
#         counter += 1
        
#     return qc_appendix




# def format_qc_metrics(qc_metrics):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with relevant QC metrics and identifiers
#     for each instrument and library source
    
#     Parameters
#     ----------
#     - qc_metrics (dict): Dictionary with QC metrics collected from the QC-etl caches
#     '''
    
#     D = {}
#     for instrument in qc_metrics:
#         for library_source in qc_metrics[instrument]:
#             metrics = get_library_metrics(library_source)
#             if instrument not in D:
#                 D[instrument] = {}
#             if library_source not in D[instrument]:
#                 D[instrument][library_source] = []
#             # collect relevant information
#             for limskey in qc_metrics[instrument][library_source]:
#                 donor = qc_metrics[instrument][library_source][limskey]['donor']
#                 libraries = list(qc_metrics[instrument][library_source][limskey]['libraries'])
#                 prefix = list(qc_metrics[instrument][library_source][limskey]['prefix'])
#                 run = list(qc_metrics[instrument][library_source][limskey]['run'])
#                 libraries = [fit_into_column(i) if len(i) >=40 else i for i in libraries]
#                 prefix = [fit_into_column(i) if len(i) >=40 else i for i in prefix]
#                 run = [fit_into_column(i) if len(i) >=40 else i for i in run]
#                 libraries.sort()
#                 prefix.sort()
#                 run.sort()
#                 sample = qc_metrics[instrument][library_source][limskey]['sample']
#                 lane_count = qc_metrics[instrument][library_source][limskey]['lane']
                
#                 L = [donor, sample, lane_count]
                
#                 for i in metrics:
#                     if i == 'read_count':
#                         L.append('{:,}'.format(qc_metrics[instrument][library_source][limskey][i]))                     
#                     else:
#                         L.append(qc_metrics[instrument][library_source][limskey][i])
#                 D[instrument][library_source].append(L)         
#                 # sort according to donor
#                 D[instrument][library_source].sort(key=lambda x: x[0])

#     return D




# def convert_read_count_to_read_pairs(D):
#     '''
#     (dict) -> None
    
#     Convert the read counts to read pair counts in place in dictionary D
    
#     Parameters
#     ----------
#     - D (dict): Dictionary with metrics of interest organized by sequencing plarform and library source
#     '''
    
#     for platform in D:
#         for library_source in D[platform]:
#             for limskeys in D[platform][library_source]:
#                 # get the number of read pairs
#                 assert D[platform][library_source][limskeys]['read_count'] % 2 == 0
#                 D[platform][library_source][limskeys]['read_count'] = int(D[platform][library_source][limskeys]['read_count'] / 2)  

    
# def get_metrics_cumulative_report(files, bamqc, rnaseqqc, cfmedipqc, emseqqc):
#     '''
#     (dict, dict, dict, dict, dict) -> dict
    
#     Returns a dictionary with identifiers and QC metrics extracted from the QC-etl caches
    
#     Parameters
#     ----------
#     - files (dict): Fastq information extracted from FPR
#     - bamqc (dict): Dictionary with bamqc QC metrics
#     - rnaseqqc (dict): Dictionary with rnaseqqc QC metrics
#     - cfmedipqc (dict): Dictionary with cfmedipqc QC metrics
#     - emseqqc (dict): Dictionary with emseqqc QC metrics
#     '''
    
#     # organize metrics for cumulative report
    
#     D = {}
    
#     # find pairs of files for each library type and instrument
#     pairs = group_fastq_pairs(files)
        
#     # Group metrics by platform and library type
#     for platform in pairs:
#         for library_source in pairs[platform]:
#             if library_source not in ['CM', 'WT', 'MC', 'MG']:
#                 # use the bamqc merged cache
#                 collect_bamqc_metrics(D, pairs, bamqc, platform, library_source)
#             elif library_source == 'WT':
#                 # use the rnaseqc merged cache
#                 collect_rnaseqc_metrics(D, pairs, rnaseqqc, platform, library_source)
#             elif library_source in ['MC', 'MG']:
#                 # use the lane level emseqqc cache
#                 collect_emseqqc_metrics(D, pairs, emseqqc, platform, library_source)
#             elif library_source == 'CM':
#                 # use the lane level cfmedip cache
#                 collect_cfmedipqc_metrics(D, pairs, cfmedipqc, platform, library_source)
    
#     return D                


# def extract_archived_emseqqc_data(archived_emseqqc):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with emseqqc from the archived qc-etl
        
#     Parameters
#     ----------
#     - archived_emseqqc (str): Path to the archived qc-etl emseqqc json 
#     '''
    
#     infile = open(archived_emseqqc)
#     data = json.load(infile)
#     infile.close()
    
#     D = {}

#     while len(data) != 0:
#         i = data.pop(0)
#         run = i['run']
#         barcode = i['barcode']
#         lane = i['lane']       
#         sample = i['pinery_lims_id']    
               
#         if run not in D:
#             D[run] = {}
#         if sample not in D[run]:
#             D[run][sample] = []
        
#         ### NEED TO IDENTIFY LAMBDA and pUC19
               
#         d = {'Barcodes': barcode,
#              'Lane Number': lane,
#              'Pinery Lims ID': sample,
#              'Run Alias': run,
#              'Lambda': 'NA',
#              'pUC19': 'NA'}
        
#         # ectract bamqc metrics
#         infile = open(i['file_bamqc'])
#         bamqc = json.load(infile)
#         infile.close()
#         d['mark duplicates_PERCENT_DUPLICATION'] = float(bamqc['mark duplicates']['PERCENT_DUPLICATION'])    

#         D[run][sample].append(d)
        
#     return D


# def extract_emseqqc_data(emseqqc_db):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with emseqqc from the qc-etl 
        
#     Parameters
#     ----------
#     - emseqqc_db (str): Path to the emseqqc SQLIte database generated by qc-etl
#     '''
    
#     conn = sqlite3.connect(emseqqc_db)
#     conn.row_factory = sqlite3.Row
#     # get all methylation data
#     data = conn.execute('select * from emseqqc_methylation_2').fetchall()
#     conn.close()
    
#     # get clumns of interest
#     columns = ['Barcodes',
#                'Lane Number',
#                'Pinery Lims ID',
#                'Run Alias',
#                'Lambda',
#                'pUC19']
               
#     D = {}

#     for i in data:
#         i = dict(i)
#         run = i['Run Alias']
#         if run not in D:
#             D[run] = {}
#         sample = i['Pinery Lims ID']
#         if sample not in D[run]:
#             D[run][sample] = []
#         d = {}
#         for j in columns:
#             try:
#                 float(i[j])
#             except:
#                 d[j] = i[j]
#             else:
#                 d[j] = round(float(i[j]), 3)
#         D[run][sample].append(d)
           
#     # add duplication rate 
#     conn = sqlite3.connect(emseqqc_db)
#     conn.row_factory = sqlite3.Row
#     data2 = conn.execute('select * from emseqqc_bamqc_2').fetchall()
#     conn.close()
    
#     for i in data2:
#         i = dict(i)
#         run = i['Run Alias']
#         sample = i['Pinery Lims ID']
#         barcodes = i['Barcodes']
#         lane = i['Lane Number']
#         duplicate = i['mark duplicates_PERCENT_DUPLICATION']
        
#         # find dict in D
#         assert run in D
#         assert sample in D[run]
#         for j in range(len(D[run][sample])):
#             if D[run][sample][j]['Barcodes'] == barcodes and D[run][sample][j]['Lane Number'] == lane and \
#                 D[run][sample][j]['Pinery Lims ID'] == sample and D[run][sample][j]['Run Alias'] == run:
#                 D[run][sample][j]['mark duplicates_PERCENT_DUPLICATION'] = duplicate    

#     return D


# def merge_qc_info(qc_info, archived_qc):
#     '''
#     (dict, dict) -> dict
    
#     Update the qc_info with data archived qc data
    
#     Parameters
#     ----------
#     - qc_info (dict): Dictionary with qc info 
#     - archived_qc (dict): Dictionary with archived qc data
#     '''
    
#     for i in archived_qc:
#         if i in qc_info:
#             for j in archived_qc[i]:
#                 if j in qc_info[i]:
#                     for d in archived_qc[i][j]:
#                         if d not in qc_info[i][j]:
#                             qc_info[i][j].append(d)
#                 else:
#                     qc_info[i][j] = archived_qc[i][j]
#         else:
#             qc_info[i] = archived_qc[i]
    
#     return qc_info




# def extract_archived_cfmedipqc_data(archived_cfmedipqc):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with archived qc-etl cfmedipqc 
        
#     Parameters
#     ----------
#     - archived_cfmedipqc (str): Path to the archived qc-etl cfmedipqc json
#     '''
    
#     infile = open(archived_cfmedipqc)
#     data = json.load(infile)
#     infile.close()
    
#     D = {}

#     while len(data) != 0:
#         i = data.pop(0)
#         run = i['run']
#         lane = i['lane']
#         barcode = i['barcode']
#         sample = i['pinery_lims_id'] 
#         metrics = collect_archived_cfmediqc_metrics(i['file'])

#         if run not in D:
#             D[run] = {}
#         if sample not in D[run]:
#             D[run][sample] = []
            
#         d = {'Barcodes': barcode, 
#              'Lane Number': lane,
#              'Pinery Lims ID': sample,
#              'Run Alias': run}
#         for j in metrics:
#             d[j] = metrics[j]
  
#         D[run][sample].append(d)
           
#     return D



# def extract_cfmedipqc_data(cfmedipqc_db):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with project-level relevant information from the cfmedip table of qc-etl
        
#     Parameters
#     ----------
#     - cfmedipqc_db (str): Path to the cfmedip SQLIte database generated by qc-etl
#     '''
    
#     #conn = sqlite3.connect('prov_report_test.db')
#     conn = sqlite3.connect(cfmedipqc_db)
#     conn.row_factory = sqlite3.Row
    
#     # get all data
#     data = conn.execute('select * from cfmedipqc_cfmedipqc_4').fetchall()

#     conn.close()
    
#     # get clumns of interest
#     columns = ['AT Dropout',
#                'Methylation beta',
#                'Percent Duplication',          
#                'Relative CpG Frequency in Regions',
#                'Barcodes',
#                'Lane Number',
#                'Pinery Lims ID',
#                'Run Alias']
    
#     D = {}

#     for i in data:
#         i = dict(i)
#         run = i['Run Alias']
#         if run not in D:
#             D[run] = {}
#         sample = i['Pinery Lims ID']
#         if sample not in D[run]:
#             D[run][sample] = []
#         d = {}
#         for j in columns:
#             try:
#                 float(i[j])
#             except:
#                 d[j] = i[j]
#             else:
#                 d[j] = round(float(i[j]), 3)
#         D[run][sample].append(d)
           
#     return D


# def extract_archived_merged_rnaseqqc_data(archived_merged_rnaseqqc, project_name):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with project-level relevant information from the archivied qc-etl rnaseqqc 
        
#     Parameters
#     ----------
#     - archived_merged_rnaseqqc (str): Path to the archived rnaseqqc json 
#     - project_name (str): Project of interest
#     '''
    
#     infile = open(archived_merged_rnaseqqc)
#     data = json.load(infile)
#     infile.close()
    
    
#     D = {}

#     while len(data) != 0:
#         i = data.pop(0)
#         d = {}
    
#         donor = i['donor']
#         group_id = i['group_id']
#         library_design = i['library_design']
#         project = i['project']
#         lims_id = i['pinery_lims_ids']
#         tissue_origin = i['tissue_origin']
#         tissue_type = i['tissue_type']
#         metrics_file = i['path'] 
        
#         d = {'Donor': donor,
#              'Project': project,
#              'Group ID': group_id,
#              'Library Design': library_design,
#              'Tissue Origin': tissue_origin,
#              'Tissue Type': tissue_type,
#              'Merged Pinery Lims ID': lims_id}
             
#         if project == project_name:
#             metrics = collect_archived_rnaseq_metrics(metrics_file)
#             for i in metrics:
#                 d[i] = metrics[i]
            
#         merged_samples = ';'.join(lims_id)
#         D[merged_samples] = d
        
#     return D




# def extract_merged_rnaseqqc_data(merged_rnaseqqc_db):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with project-level relevant information from the rnaseqqc table of qc-etl
        
#     Parameters
#     ----------
#     - rnaseqqc_db (str): Path to the rnaseq SQLIte database generated by qc-etl
#     '''
    
#     #conn = sqlite3.connect('prov_report_test.db')
#     conn = sqlite3.connect(merged_rnaseqqc_db)
#     cur = conn.cursor()
#     cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
#     table_name = [i[0] for i in cur][0]
#     conn.row_factory = sqlite3.Row
    
#     # get all data
#     data = conn.execute('select * from {0}'.format(table_name)).fetchall()
#     conn.close()

#     D = {}

#     columns = ['Donor',
#                'Project',
#                'Merged Pinery Lims ID',
#                'Group ID',
#                'Library Design',       
#                'Tissue Origin',       
#                'Tissue Type',       
#                'PCT_CODING_BASES', 
#                'rrna contamination properly paired',
#                'rrna contamination in total (QC-passed reads + QC-failed reads)',
#                'PCT_CORRECT_STRAND_READS',
#                'MEDIAN_5PRIME_TO_3PRIME_BIAS']
        
    
#     for i in data:
#         d = {}
#         for j in columns:
#             assert j in dict(i).keys()
#             if j == 'Merged Pinery Lims ID':
#                 d[j] = list(map(lambda x: x.strip(), i[j].replace('[', '').replace(']', '').replace('\"', '').split(',')))
#             else:
#                 d[j] = i[j]
        
#         merged_samples = i['Merged Pinery Lims ID']
#         merged_samples = list(map(lambda x: x.strip(), merged_samples.replace('[', '').replace(']', '').replace('\"', '').split(',')))
#         merged_samples = ';'.join(merged_samples)
        
#         D[merged_samples] = d
        
#     return D
    

# def extract_archived_merged_bamqc_data(archived_merged_bamqc, project_name):
#     '''
#     (str, str) -> dict
    
#     Returns a dictionary with project-level relevant information from the archived qc-etl bamqc
        
#     Parameters
#     ----------
#     - archived_merged_bamqc (str): Path to the archived bamqc json 
#     - project_name (str): Project of interest
#     '''
            
#     infile = open(archived_merged_bamqc)
#     data = json.load(infile)
#     infile.close()
    
#     D = {}

#     while len(data) != 0:
#         i = data.pop(0)
#         d = {}
#         donor = i['donor']
#         project = i['project']
#         tissue_origin = i[ 'tissue_origin']
#         library_design = i['library_design']
#         tissue_type = i['tissue_type']
#         group_id = i['group_id']
#         pinery_lims_ids = i['pinery_lims_ids']
        
#         # only collect information for the relevant project
#         if project == project_name:
#             metrics = collect_archived_bamqc_metrics(i['path'])
#             bases_mapped = int(metrics['bases mapped'])
#             instrument = metrics['instrument']
#             library = metrics['library'].split(',') 
#             library = list(map(lambda x: x.strip(), library))
#             sample = metrics['sample']
#             total_reads = int(metrics['total reads'])
#             total_bases_on_target = int(metrics['total bases on target'])
#             mark_duplicates = float(metrics['mark duplicates']['PERCENT_DUPLICATION'])
#             # compute on_target rate, not available through qc-etl
#             on_target = compute_on_target_rate(bases_mapped, total_bases_on_target) 
        
#             # 'coverage histogram',
#             ### need to replace coverage and coverage deduplicated
        
#             d = {'Project': project, 'Donor': donor, 'Tissue Origin': tissue_origin,
#                  'Library Design': library_design, 'Tissue Type': tissue_type,
#                  'Group ID': group_id, 'Merged Pinery Lims ID': pinery_lims_ids,
#                  'library': library, 'sample': sample, 'instrument': instrument,
#                  'total reads': total_reads, 'mapped reads': bases_mapped,
#                  'total bases on target': total_bases_on_target, 
#                  'mark duplicates_PERCENT_DUPLICATION': mark_duplicates,
#                  'on_target': on_target,
#                  'coverage': 0, 'coverage deduplicated': 0}
               
#             merged_samples = ';'.join(pinery_lims_ids)
#             D[merged_samples] = d
        
#     return D


# def extract_merged_bamqc_data(merged_bamqc_db):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with project-level relevant information from the bamqc table of qc-etl
        
#     Parameters
#     ----------
#     - bamqc_db (str): Path to the bamqc SQLite database generated by qc-etl
#     '''
            
#     conn = sqlite3.connect(merged_bamqc_db)
#     cur = conn.cursor()
    
#     # get table name
#     cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
#     table_name = [i[0] for i in cur][0]
#     # get all data from table
#     conn.row_factory = sqlite3.Row
#     data = conn.execute('select * from {0}'.format(table_name)).fetchall()
        
#     columns = ['library', 'Library Design', 'Donor', 'Project', 'Tissue Origin', 
#                'Tissue Type', 'Group ID', 'sample', 'Merged Pinery Lims ID', 'instrument',
#                'sample', 'bases mapped', 'coverage', 'coverage deduplicated',
#                'mapped reads', 'mark duplicates_PERCENT_DUPLICATION', 'total bases on target',
#                'total reads']
    
#     D = {}
    
#     for i in data:
#         d = {}
#         for j in columns:
#             if j in ['coverage', 'coverage deduplicated', 'mark duplicates_PERCENT_DUPLICATION']:
#                 d[j] = float(i[j])
#             elif j in ['bases mapped', 'mapped reads', 'total bases on target', 'total reads']:
#                 d[j] = int(i[j])
#             elif j == 'library':
#                 d[j] = i[j].replace('\"', '').split(',')
#             elif j == 'instrument':
#                 d[j] = i[j].replace('\"', '')
#             elif j == 'Merged Pinery Lims ID':
#                 d[j] = list(map(lambda x: x.strip(), i[j].replace('[', '').replace(']', '').replace('\"', '').split(',')))
#             else:
#                 d[j] = i[j]
#         # compute on_target rate, not available through qc-etl
#         d['on_target'] = compute_on_target_rate(d['bases mapped'], d['total bases on target']) 
        
#         merged_samples = i['Merged Pinery Lims ID']
#         merged_samples = list(map(lambda x: x.strip(), merged_samples.replace('[', '').replace(']', '').replace('\"', '').split(',')))
        
#         merged_samples = ';'.join(merged_samples)
#         D[merged_samples] = d
        
#     return D


# def format_sequencing_table(files):
#     '''
#     (dict) -> list
    
#     Returns a list with sequencing information
    
#     Parameters
#     ----------
#     - files (dict): Dictionary with file information extracted from FPR
#     '''

#     pairs = group_fastq_pairs(files)

#     SL = []
    
#     for instrument in pairs:
#         for library_source in pairs[instrument]:
#             for i in pairs[instrument][library_source]:
#                 donor = i[0]['sample_name']
#                 sample = i[0]['sample_id'][0]
#                 library = i[0]['library'][0]
#                 library_type = i[0]['library_source'][0]
#                 prefix = i[0]['prefix']
#                 read_pairs = i[0]['read_count']
#                 L = [donor, sample, library, prefix, read_pairs]
#                 if L not in SL:
#                     SL.append(L)
#                 # sort according to donor, sample, snf library
#                 SL.sort(key = lambda x: (x[0], x[1], x[2]))
#     # adjust long names to fit in column
#     for i in range(len(SL)):
#         L = [fit_into_column(j) if len(j) >=30 else j for j in SL[i][:-1]]
#         L.append('{:,}'.format(SL[i][-1]))
#         SL[i] = L
    
#     return SL



# def get_identifiers_appendix(files, report):
#     '''
#     (dict, str) -> list
    
#     Returns a list with definitions of columns in the sample identifier table
    
#     Parameters
#     ----------
#     - files (dict): Dictionary with file information from released libraries extracted from FPR
#     - report (str): cumulative or batch report
#     '''

#     # get the library type, tissue type and tissue origin 
#     D = get_library_tissue_types(files)
    
#     if report  == 'batch':
#         L = ['Library Id: OICR-generated library identifier',
#              'OICR Case Id: OICR-generated case identifier',
#              'Donor Id: user supplied donor identifier',
#              'Sample Id: user supplied sample, this distinguishes distinct samples of the same type from the same donor. If only one sample per donor is submitted the value may match the donor Id',
#              'Sample Description: a description of the Sample Id',
#              'Library Type (LT): {0}'.format(D['Library Type']),
#              'Tissue Origin (TO): {0}'.format(D['Tissue Origin']),
#              'Tissue Type (TT): {0}'.format(D['Tissue Type'])]         
#     elif report == 'cumulative':
#         L = ['OICR Case Id: OICR-generated case identifier',
#              'Donor Id: user supplied donor identifier',
#              'OICR Sample Id: The OICR generated sample identifier. The sample Id is formed from the following: 1. Case Id, 2. Tissue Origin, 3. Tissue Type, 4. Library Type and 5. User supplied Sample Id',
#              'Sample Id: user supplied sample, this distinguishes distinct samples of the same type from the same donor. If only one sample per donor is submitted the value may match the donor Id',
#              'Sample Description: a description of the Sample Id',
#              'Library Type (LT): {0}'.format(D['Library Type']),
#              'Tissue Origin (TO): {0}'.format(D['Tissue Origin']),
#              'Tissue Type (TT): {0}'.format(D['Tissue Type'])]
       
#     return L



# def group_cumulative_samples(files):
#     '''
#     (dict) -> list
    
#     Returns a list of sample identifiers sorted on donor id 
    
#     Parameters
#     ----------
#     - files (dict): Information about fastqs extracted from File Provenance Report
#     '''
        
#     # make a list of instruments
#     instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
    
#     # record sample identifiers across library type and platform
#     SL = []
    
#     # find pairs of fastqs
#     for platform in instruments:
#         pairs = find_fastq_pairs(files, platform)
#         add_file_prefix(pairs)
#         for i in pairs:
#             # add comma to reads for readability
#             case = i[0]['sample_name']
#             external_name = i[0]['external_name']
#             sample = i[0]['sample_id'][0]
#             library_source = i[0]['library_source'][0]
#             library = i[0]['library'][0]
#             prefix = i[0]['prefix']
#             groupid = i[0]['groupid'][0]
#             group_description = i[0]['groupdesc'][0]
            
            
#             max_length = 20
            
#             # reformat prefix to fit the table column
#             if len(prefix) >= max_length:
#                 prefix = fit_into_column(prefix)
#             if len(library) >= max_length:
#                 library = fit_into_column(library)
#             if len(sample) >= max_length:
#                 sample = fit_into_column(sample)
#             if len(groupid) >= max_length:
#                 groupid = fit_into_column(groupid)
#             if len(group_description) >= max_length:
#                 group_description = fit_into_column(group_description)
                        
#             tissue_origin = i[0]['tissue_origin'][0]
#             tissue_type = i[0]['tissue_type'][0]
#             L = [case, external_name, sample, group_description, library_source, tissue_origin, tissue_type]

#             if L not in SL:
#                 SL.append(L)
#     SL.sort(key = lambda x: x[0])
    
#     return SL            


# def count_released_fastqs_by_library_type_instrument(FPR_info):
#     '''
#     (dict, str) -> dict
    
#     Returns the count of released fastqs for each library type, run and instrument
#     Precondition: Fastqs are paired, return the count of R1
    
#     Parameters
#     ----------
#     - FPR_info (dict): Information about the released fastqs collected from File Provenance Report
#     '''
        
#     # count released fastqs by instrument and run
#     D = {}
#     for file in FPR_info:
#         instrument = FPR_info[file]['platform']
#         assert len(FPR_info[file]['run_id']) == 1
#         run = FPR_info[file]['run_id'][0]
#         library_type = FPR_info[file]['library_source'][0]
#         if library_type not in D:
#             D[library_type] = {}
#         if instrument not in D[library_type]:
#             D[library_type][instrument] = {}
#         if 'R1' in FPR_info[file]['file_path']:
#             if run not in D[library_type][instrument]:
#                 D[library_type][instrument][run] = 1
#             else:
#                 D[library_type][instrument][run] += 1
#     return D


# def resolve_links(filenames):
#     '''
#     (list) -> list

#     Returns a list of file paths with links resolved

#     Parameters
#     ----------
#     - filenames(list): List of file paths
#     '''
    
#     files = [os.path.realpath(i) for i in filenames]
#     return files


# def list_released_fastqs_project(api, project):
#     '''
#     (str, str) -> list
    
#     Returns a list of fastqs for a given project that were previously released by interrogating QC status in Nabu
#     Pre-condition: Released files need to be marked in Nabu
        
#     Parameters
#     ----------
#     - api (str): URL of the nabu API
#     - project (str): file_swid (str): File unique identifier
#     '''
    
#     # get end-point
#     api += 'get-fileqcs' if api[-1] == '/' else '/get-fileqcs'
    
#     R = []
 
#     # check each fastq-generating workflow
#     headers = {'accept': 'application/json','Content-Type': 'application/json'}
#     json_data = {"project": "{0}".format(project)}
#     response = requests.post(api, headers=headers, json=json_data)
#     # check response code
#     if response.status_code == 200:
#         L = response.json()['fileqcs']
#         if L:
#             for i in L:
#                 #file_swid = i['fileid']
#                 qc_status = i['qcstatus']
#                 file_path = i['filepath']
#                 if qc_status.upper() == 'PASS':
#                     R.append(file_path)
#     return R    


# def parse_fpr_records(provenance, project, workflow, prefix=None):
#     '''
#     (str, str, list, str | None) -> dict
  
#     Returns a dictionary with file info extracted from FPR for a given project 
#     and a given workflow if workflow is speccified. 
            
#     Parameters
#     ----------
#     - provenance (str): Path to File Provenance Report
#     - project (str): Project name as it appears in File Provenance Report. 
#     - workflow (list): List of workflows used to generate the output files.
#     - prefix (str | None): Prefix used to recover file full paths when File Provevance contains relative paths.
#     '''
    
#     # create a dict {file_swid: {file info}}
#     D  = {}
    
#     # get all the records for a single project
#     records = get_FPR_records(project, provenance)
    
#     # parse the records and get all the files for a given project
#     for i in records:
#         # keep records for project
#         if project == i[1]:
#             pipeline_workflow = i[30]
#             # check workflow
#             if len(workflow) == 1:
#                 if workflow[0].lower() == 'bcl2fastq':
#                     # skip if not fastq-related workflows    
#                     if pipeline_workflow.lower() not in ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']:
#                         continue
#                 else:
#                     # skip if not provided workflow
#                     if workflow[0].lower() != pipeline_workflow.lower():
#                         continue
#             else:
#                 if pipeline_workflow.lower() not in list(map(lambda x: x.lower(), workflow)):
#                     continue
            
#             # get file path
#             if prefix:
#                 file_path = os.path.join(prefix, i[46])
#             else:
#                 file_path = i[46]
#             # get file name
#             file_name = os.path.basename(file_path)
#             # get sample name
#             sample_name = i[7]
#             # get parent sample name
#             parent_sample = i[9].split(':')[0]
#             # get time stamp and convert to epoch
#             creation_date = i[0]
#             # remove milliseconds
#             creation_date = creation_date.split('.')[0]
#             pattern = '%Y-%m-%d %H:%M:%S'
#             creation_date = int(time.mktime(time.strptime(creation_date, pattern)))
#             # record platform
#             platform = i[22]
            
#             # get md5sum
#             md5 = i[47]
#             # get workdlow swid
#             workflow_run_id = i[36]
#             # get workflow version
#             workflow_version = i[31]
#             # get file swid
#             file_swid = i[44]
  
#             # for merging workflows there will be multiple values for the variables below
#             # get library aliquot
#             aliquot = i[56].split('_')[-1]
#             # get library aliases
#             library = i[13]
#             # get lims key
#             limskey = i[56]
#             # get run id
#             run_id = i[18]
#             # get run
#             run = i[23]   
#             # get lane
#             lane = i[24]
#             # get barcode
#             barcode = i[27]  
          
#             geo = i[12]
#             if geo:
#                 geo = {k.split('=')[0]:k.split('=')[1] for k in geo.split(';')}
#             else:
#                 geo = {}
#             for j in ['geo_external_name', 'geo_group_id', 'geo_group_id_description',
#                       'geo_targeted_resequencing', 'geo_library_source_template_type',
#                       'geo_tissue_type', 'geo_tissue_origin']:
#                 if j not in geo:
#                     geo[j] = 'NA'
#                 if j == 'geo_group_id':
#                     # removes misannotations
#                     geo[j] = geo[j].replace('&2011-04-19', '').replace('2011-04-19&', '')
       
#             read_count = i[45]
#             if read_count:
#                 read_count = {k.split('=')[0]:k.split('=')[1] for k in i[45].split(';')}
#             if 'read_count' in read_count:
#                 read_count = int(float(read_count['read_count']))
#             else:
#                 read_count = -1
       
#             sample_id = sample_name + '_' + geo['geo_tissue_origin']+ '_' + geo['geo_tissue_type'] + '_' + geo['geo_library_source_template_type'] + '_' + geo['geo_group_id']
         
#             d = {'workflow': pipeline_workflow, 'file_path': file_path, 'file_name': file_name,
#                  'sample_name': sample_name, 'creation_date': creation_date, 'platform': platform,
#                  'md5': md5, 'workflow_run_id': workflow_run_id, 'workflow_version': workflow_version,
#                  'file_swid': file_swid, 'external_name': geo['geo_external_name'],
#                  'panel': geo['geo_targeted_resequencing'], 'library_source': [geo['geo_library_source_template_type']],
#                  'parent_sample': [parent_sample], 'run_id': [run_id], 'run': [run],
#                  'limskey': [limskey], 'aliquot': [aliquot], 'library': [library],
#                  'barcode': [barcode], 'tissue_type': [geo['geo_tissue_type']],
#                  'tissue_origin': [geo['geo_tissue_origin']], 'groupdesc': [geo['geo_group_id_description']],
#                  'groupid': [geo['geo_group_id']], 'read_count': read_count, 'sample_id': [sample_id], 'lane': [lane]}
            
#             if file_swid not in D:
#                 D[file_swid] = d
#             else:
#                 assert D[file_swid]['file_path'] == file_path
#                 assert D[file_swid]['external_name'] == geo['geo_external_name']
#                 assert D[file_swid]['read_count'] == read_count
#                 D[file_swid]['sample_id'].append(sample_id)
#                 D[file_swid]['parent_sample'].append(parent_sample)
#                 D[file_swid]['run_id'].append(run_id)
#                 D[file_swid]['run'].append(run)
#                 D[file_swid]['limskey'].append(limskey)
#                 D[file_swid]['aliquot'].append(aliquot)
#                 D[file_swid]['library'].append(library)
#                 D[file_swid]['barcode'].append(barcode)
#                 D[file_swid]['tissue_type'].append(geo['geo_tissue_type'])
#                 D[file_swid]['tissue_origin'].append(geo['geo_tissue_origin'])
#                 D[file_swid]['library_source'].append(geo['geo_library_source_template_type'])
#                 D[file_swid]['groupdesc'].append(geo['geo_group_id_description'])
#                 D[file_swid]['groupid'].append(geo['geo_group_id'])
#                 D[file_swid]['lane'].append(lane)
    
#     return D    



# def collect_archived_cfmediqc_metrics(metrics_file):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with cfmedipqc metrics 
    
#     Parameters
#     ----------
#     - metrics_file (str): Path to the json file with metrics
#     '''
    
#     infile = open(metrics_file)
#     metrics = json.load(infile)
#     infile.close()
    
#     ## need to identify 'Relative CpG Frequency in Regions'
    
#     D ={"AT Dropout": float(metrics["AT_DROPOUT"]),
#         "Percent Duplication": float(metrics["PERCENT_DUPLICATION"]),
#         'Methylation beta': float(metrics["THALIANA_BETA"]),
#         'Relative CpG Frequency in Regions': 'NA'}
              
#     return D    
      


# def collect_archived_bamqc_metrics(metrics_file):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with bamqc metrics 
    
#     Parameters
#     ----------
#     - metrics_file (str): Path to the json file with metrics
#     '''
    
    
#     infile = open(metrics_file)
#     metrics = json.load(infile)
#     infile.close()
    
#     L = ['bases mapped', 'instrument', 'library', 'group id', 'donor',
#          'sample', 'tissue origin', 'tissue type', 'total reads', 
#          'total bases on target', 'bases mapped', 'mark duplicates']
        
#     D = {i: metrics[i]  for i in L}
        
#     return D    




# def collect_archived_rnaseq_metrics(metrics_file):
#     '''
#     (str) -> dict
    
#     Returns a dictionary with rnaseqqc metrics 
    
#     Parameters
#     ----------
#     - metrics_file (str): Path to the json file with metrics
#     '''
    
    
#     infile = open(metrics_file)
#     metrics = json.load(infile)
#     infile.close()
    
#     rrna_total = int(metrics['rrna_contamination']['in total (QC-passed reads + QC-failed reads)'])
#     rrna_properly_paired = int(metrics['rrna_contamination']['properly paired'])                 
#     rrna_contamination = round(rrna_properly_paired / rrna_total * 100, 3)
#     pct_coding_bases = float(metrics['picard']['metrics']['PCT_CODING_BASES'])
#     bias = float(metrics['picard']['metrics']['MEDIAN_5PRIME_TO_3PRIME_BIAS'])
#     pct_correct_strand = float(metrics['picard']['metrics']['PCT_CORRECT_STRAND_READS'])
       
#     D = {'rrna_contamination': rrna_contamination,
#          'rrna contamination properly paired': rrna_properly_paired,
#          'rrna contamination in total (QC-passed reads + QC-failed reads)': rrna_total,
#          'PCT_CODING_BASES': pct_coding_bases,
#          'PCT_CORRECT_STRAND_READS': pct_correct_strand,
#          'MEDIAN_5PRIME_TO_3PRIME_BIAS': bias}
        
#     return D


# def group_fastq_pairs(files):
#     '''
#     (dict) -> dict
    
#     Returns a dictionary with lists of fastq pairs for each instrument and library type
#     Pre-condition: fastqs are paired with R1 and R2
    
#     Parameters
#     ----------
#     - files (dict): Dictionary with fastq information extracted from FPR
#     '''
  
#     # find pairs of files for each library type and instrument
#     pairs = {}
#     instruments = list(set([files[file_swid]['platform'] for file_swid in files]))
#     for platform in instruments:
#         P = find_fastq_pairs(files, platform)
#         library_sources = list(set([i[0]['library_source'][0] for i in P]))
#         pairs[platform] = {}
#         for j in library_sources:
#             L = [i for i in P if i[0]['library_source'][0] == j]
#             add_file_prefix(L)
#             pairs[platform][j] = L
#     return pairs


# def collect_bamqc_metrics(D, pairs, bamqc, platform, library_source):
#     '''
#     (dict, dict, dict, str, str) -> None
    
#     Update dictionary D with QC metrics from bamqc for fastq pairs sorted by instrument and library type
    
#     Parameters
#     ----------
#     - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
#     - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
#     - bamqc (dict): Dictionary with QC metrics from the merged bamqc cache 
#     - platform (str): Sequencing platform
#     - library_source (str): 2-letters code of library type
#     '''
    
#     for i in bamqc:
#         merged_limskeys = sorted(bamqc[i]['Merged Pinery Lims ID'])            
#         # collect limkskeys from FPR records
#         L = []
#         readcount = []
#         prefix = []
#         runs = []
#         lanes = []
        
#         # collect read count for all fastqs
#         # track limskeys and file prefix for fastqs
#         for j in pairs[platform][library_source]:
#             for k in j:
#                 limskey = k['limskey'][0]
#                 if library_source == bamqc[i]['Library Design'] \
#                 and platform == '_'.join(bamqc[i]['instrument'].split()) \
#                 and limskey in merged_limskeys \
#                 and k['groupid'][0] == bamqc[i]['Group ID']:
#                     assert k['sample_id'][0] == bamqc[i]['sample']
#                     # get the number of read for each file
#                     readcount.append(k['read_count'])
#                     L.append(limskey)
#                     prefix.append(k['prefix'])
#                     runs.extend(k['run_id'])
#                     lanes.extend(k['run'])
#         # check that all limskeys have been recorded
#         if sorted(list(set(L))) == sorted(merged_limskeys):
#             assert sum(readcount) % 2 == 0
#             readcount = sum(readcount)
#             prefix = list(set(prefix))
#             runs = sorted(list(set(runs)))
#             lanes = list(set(lanes))
#             libraries = bamqc[i]['library']
#             donor = bamqc[i]['Donor']
#             coverage = round(bamqc[i]['coverage'], 2)
#             coverage_dedup = round(bamqc[i]['coverage deduplicated'], 2)
#             on_target = round(bamqc[i]['on_target'], 2)
#             percent_duplicate = round(bamqc[i]['mark duplicates_PERCENT_DUPLICATION'], 2)
#             sample = bamqc[i]['sample']

#             d = {'read_count': readcount,
#                  'prefix': prefix,
#                  'run': runs,
#                  'lane': len(lanes),
#                  'libraries': libraries,
#                  'donor': donor,
#                  'coverage': coverage,
#                  'coverage_dedup': coverage_dedup,
#                  'on_target': on_target,
#                  'percent_duplicate': percent_duplicate,
#                  'sample': sample}
                                   
#             if platform not in D:
#                 D[platform] = {}
#             if library_source not in D[platform]:
#                 D[platform][library_source] = {}
#             assert i not in D[platform][library_source]
#             D[platform][library_source][i] = d

# def collect_rnaseqc_metrics(D, pairs, rnaseqqc, platform, library_source):
#     '''
#     (dict, dict, dict, str, str) -> None
    
#     Update dictionary D with QC metrics from rnaseqc for fastq pairs sorted by instrument and library type
    
#     Parameters
#     ----------
#     - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
#     - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
#     - rnaseqqc (dict): Dictionary with QC metrics from the merged rnaseqqc cache 
#     - platform (str): Sequencing platform
#     - library_source (str): 2-letters code of library type
#     '''
    
#     for i in rnaseqqc:
#         merged_limskeys = sorted(rnaseqqc[i]['Merged Pinery Lims ID'])            
#         # collect limkskeys from FPR records
#         L = []
#         readcount = []
#         prefix = []
#         # track library ids because not present in rnaseq merged cache
#         libraries = []
#         runs = []
#         lanes = []
        
#         # collect read count for all fastqs
#         # track limskeys and file prefix for fastqs
#         for j in pairs[platform][library_source]:
#             for k in j:
#                 limskey = k['limskey'][0]
#                 if library_source == rnaseqqc[i]['Library Design'] \
#                 and limskey in merged_limskeys \
#                 and k['groupid'][0] == rnaseqqc[i]['Group ID']:
#                     # get the number of read for each file
#                     readcount.append(k['read_count'])
#                     L.append(limskey)
#                     prefix.append(k['prefix'])
#                     libraries.extend(k['library'])
#                     runs.extend(k['run_id'])
#                     lanes.extend(k['run'])
#         # check that all limskeys have been recorded
#         if sorted(list(set(L))) == sorted(merged_limskeys):
#             assert sum(readcount) % 2 == 0
#             readcount = sum(readcount)
#             prefix = list(set(prefix))
#             runs = sorted(list(set(runs)))
#             lanes = list(set(lanes))
#             libraries = list(set(libraries))
#             donor = rnaseqqc[i]['Donor']
#             bias = rnaseqqc[i]['MEDIAN_5PRIME_TO_3PRIME_BIAS']
#             contamination = round((rnaseqqc[i]['rrna contamination properly paired'] / rnaseqqc[i]['rrna contamination in total (QC-passed reads + QC-failed reads)'] * 100), 3)
#             coding = round(rnaseqqc[i]['PCT_CODING_BASES'], 3)
#             strand = round(rnaseqqc[i]['PCT_CORRECT_STRAND_READS'], 3)
#             sample = '_'.join([donor, rnaseqqc[i]['Tissue Origin'], rnaseqqc[i]['Tissue Type'], rnaseqqc[i]['Library Design'], rnaseqqc[i]['Group ID']])

#             d = {'read_count': readcount,
#                  'prefix': prefix,
#                  'run': runs,
#                  'lane': len(lanes),
#                  'libraries': libraries,
#                  'donor': donor,
#                  "5'-3' bias": bias,
#                  'rRNA contamination': contamination,
#                  'Coding (%)': coding,
#                  'Correct strand reads (%)': strand,
#                  'sample': sample}
                     
#             if platform not in D:
#                 D[platform] = {}
#             if library_source not in D[platform]:
#                 D[platform][library_source] = {}
#             assert i not in D[platform][library_source]
#             D[platform][library_source][i] = d



# def collect_cfmedipqc_metrics(D, pairs, cfmedipqc, platform, library_source):
#     '''
#     (dict, dict, dict, str, str) -> None
    
#     Update dictionary D with QC metrics from cfmedipqc for fastq pairs sorted by instrument and library type
    
#     Parameters
#     ----------
#     - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
#     - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
#     - cfmedipqc (dict): Dictionary with QC metrics from the cfmedipqc cache 
#     - platform (str): Sequencing platform
#     - library_source (str): 2-letters code of library type
#     '''
    
#     for i in pairs[platform][library_source]:
#         for j in i:
#             limskey = j['limskey'][0]
#             run = j['run_id'][0]
#             prefix = j['prefix']
#             barcode = j['barcode'][0]
#             lane = j['lane'][0]
#             library = j['library'][0]
#             donor = j['sample_name']
#             # get the read count for each file
#             readcount = j['read_count']
#             sample = j['sample_id'][0]
#             lanes = j['run']
                        
#             if run in cfmedipqc and limskey in cfmedipqc[run]:
#                 assert len(cfmedipqc[run][limskey]) == 1
#                 d = cfmedipqc[run][limskey][0]
#                 assert d['Pinery Lims ID'] == limskey
#                 assert d['Run Alias'] == run
#                 assert int(d['Lane Number']) == int(lane)
#                 assert d['Barcodes'] == barcode
#                 AT_dropout = round(d['AT Dropout'], 3)
#                 methylation_beta = round(d['Methylation beta'], 3)
#                 duplication = round(d['Percent Duplication'], 3)
#                 CpG_enrichment = round(d['Relative CpG Frequency in Regions'], 3)
                    
#                 if platform not in D:
#                     D[platform] = {}
#                 if library_source not in D[platform]:
#                     D[platform][library_source] = {}
#                 if limskey in D[platform][library_source]:
#                     assert D[platform][library_source][limskey]['prefix'][0] == prefix
#                     D[platform][library_source][limskey]['read_count'] += readcount
#                     assert D[platform][library_source][limskey]['read_count'] % 2 == 0
#                     assert D[platform][library_source][limskey]['sample'] == sample
#                 else:
#                     D[platform][library_source][limskey] = {
#                         'read_count': readcount,
#                         'prefix': [prefix],
#                         'run': run,
#                         'lane': len(lanes),
#                         'libraries': [library],
#                         'donor': donor,
#                         'AT_dropout': AT_dropout,
#                         'methylation_beta': methylation_beta,
#                         'duplication': duplication,
#                         'CpG_enrichment': CpG_enrichment,
#                         'sample': sample
#                         }



# def collect_emseqqc_metrics(D, pairs, emseqqc, platform, library_source):
#     '''
#     (dict, dict, dict, str, str) -> None
    
#     Update dictionary D with QC metrics from emseqqc for fastq pairs sorted by instrument and library type
    
#     Parameters
#     ----------
#     - D (dict): Dictionary with QC metrics and identifiers for each sequencing platform and library type
#     - pairs (dict): Dictionary of lists of fastq pairs for each sequencing platform and library type
#     - emseqqc (dict): Dictionary with QC metrics from the emseqqc cache 
#     - platform (str): Sequencing platform
#     - library_source (str): 2-letters code of library type
#     '''
    
    
#     for i in pairs[platform][library_source]:
#         for j in i:
#             limskey = j['limskey'][0]
#             run = j['run_id'][0]
#             prefix = j['prefix']
#             barcode = j['barcode'][0]
#             lane = j['lane'][0]
#             lanes = j['run']
#             library = j['library'][0]
#             donor = j['sample_name']
#             sample = j['sample_id'][0]
#             # get the number of reads for each file
#             readcount = j['read_count']
#             if run in emseqqc and limskey in emseqqc[run]:
#                 assert len(emseqqc[run][limskey]) == 1
#                 d = emseqqc[run][limskey][0]
#                 assert d['Pinery Lims ID'] == limskey
#                 assert d['Run Alias'] == run
#                 assert int(d['Lane Number']) == int(lane)
#                 assert d['Barcodes'] == barcode
#                 Lambda_methylation = round(d['Lambda'], 3)
#                 pUC19_methylation = round(d['pUC19'], 3)
#                 percent_duplication = round(d['mark duplicates_PERCENT_DUPLICATION'], 3)

#                 if platform not in D:
#                     D[platform] = {}
#                 if library_source not in D[platform]:
#                     D[platform][library_source] = {}
#                 if limskey in D[platform][library_source]:
#                     assert D[platform][library_source][limskey]['prefix'][0] == prefix
#                     D[platform][library_source][limskey]['read_count'] += readcount
#                     assert D[platform][library_source][limskey]['read_count'] % 2 == 0
#                     assert D[platform][library_source][limskey]['sample'] == sample
#                 else:
#                     D[platform][library_source][limskey] = {
#                         'read_count': readcount,
#                         'prefix': [prefix],
#                         'run': run,
#                         'lane': len(lanes),
#                         'libraries': [library],
#                         'donor': donor,
#                         'Lambda_methylation': Lambda_methylation,
#                         'pUC19_methylation': pUC19_methylation,
#                         'percent_duplication': percent_duplication,
#                         'sample': sample
#                         }



###### end of functions for cumulative report 




def load_data(provenance_data_file):
    '''
    (str) -> list
    
    Returns the list of data contained in the provenance_data_file
    
    Parameters
    ----------
    - provenance_data_file (str): Path to the file with production data extracted from Shesmu
    '''

    infile = open(provenance_data_file, encoding='utf-8')
    provenance_data = json.load(infile)
    infile.close()
    
    return provenance_data



def get_donor_name(case_data):
    '''
    (str) -> str
    
    Returns the name of the donor in case_data
    
    Parameters
    ----------
    - case_data (dict): Dictionary with a single case data 
    '''

    donor = list(set([i['donor'] for i in case_data['sample_info']]))
    assert len(donor) == 1
    donor = donor[0]

    return donor


def get_file_attributes(attributes):
    '''
    (str) -> dict
    
    Returns a dictionary of file attributes
    
    Parameters
    ----------
    - attributes (str): String representation of a list of file attributes
    '''
    
    if attributes:
        attributes = json.loads(attributes)
        file_attributes = {}
        for k in attributes:
            file_attributes[k] = attributes[k][0]
        if len(file_attributes) == 0:
            file_attributes = {}
        else:
            file_attributes = json.dumps(file_attributes)
    else:
        file_attributes = {}
    
    return file_attributes


def extract_sample_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with sample information for each lims_id of a case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case production data
    '''
        
    D = {}
    
    for d in case_data['sample_info']:
        lims_id = d['limsId']
        barcode = d['barcode']
        donor = d['donor']
        external_id = d['externalId']
        if d['groupId']:
            group_id = d['groupId']
        else:
            group_id = 'NA'
        if d['groupDesc']:
            group_description = d['groupDesc']
        else:
            group_description = 'NA'
        lane = d['lane']
        library = d['library']
        library_design = d['libraryDesign']
        run = d['run']
        sample_id = d['sampleId']
        tissue_origin = d['tissueOrigin']
        tissue_type = d['tissueType']
        instrument = d['instrument']
        
        assert lims_id not in D
        D[lims_id] = {'lims_id': lims_id, 'barcode': barcode, 'donor': donor,
                      'external_id': external_id, 'group_id': group_id,
                      'group_description': group_description, 'lane': lane,
                      'library': library, 'library_design': library_design,
                      'run': run, 'sample_id': sample_id, 'tissue_origin': tissue_origin,
                      'tissue_type': tissue_type, 'instrument': instrument}
            
    return D            
        
        
        
def extract_file_info(case_data):
    '''
    (dict) -> dict
    
    Returns a dictionary with file information for each file of a case
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case production data
    '''

    D = {}
    
    projects = [{i['project']: i['deliverables']} for i in case_data['project_info']]
    
    for d in case_data['workflow_runs']:
        files = json.loads(d['files'])
        lims = d['limsIds'].split(',')
        workflow = d['wf']
        wfrun_id = d['wfrunid']
        version = d['wfv']
        for k in files:
            file = k['path']
            md5sum = k['md5']
            accession = k['accession']
            file_attributes = json.loads(k['file_attributes'])
            
            assert file not in D
            D[file] = {'lims': lims, 'workflow': workflow, 'wfrun_id': wfrun_id,
                       'version': version, 'md5sum': md5sum, 'accession': accession,
                       'attributes': file_attributes, 'case_id': case_data['case'],
                       'project': projects}
            
    return D                   
            
        
def add_sample_info(file_info, sample_info):
    '''
    (dict, dict) -> dict
    
    Returns a dictionary with file information including the corresponding sample
    information for all files in a case
    
    Parameters
    ----------
    - file_info (dict): Dictionary with information about all files in case
    - sample_info (dict): Dictionarty with information about all samples in case
    '''
    
    for file in file_info:
        for limsid in file_info[file]['lims']:
            assert limsid in sample_info
            if 'samples'  not in file_info[file]:
                file_info[file]['samples'] = [sample_info[limsid]] 
            else:
                if sample_info[limsid] not in file_info[file]['samples']:
                    file_info[file]['samples'].append(sample_info[limsid])
            
    return file_info            
 
    
 
    
def is_correct_project(projects, expected_project):
    '''
    (dict, str) -> bool     
    
    Returns True if expected_project is member of the projects of a case
    
    Parameters
    ----------
    - projects (list): List of dictionaries project: deliverables associated with a file
    - expected_project (str): Project name 
    '''
    
    L = []
    for i in projects:
        L.extend(list(i.keys()))
        
    return expected_project in L
            
 
    
def is_correct_workflow(workflow, workflows):
    '''
    (str, list) -> bool
    
    Returns True if worfklow is a member of workflows
        
    Parameters
    ----------
    - workflow (str): Name of workflow
    - workflows (list): List of workflows generating the data to release
    '''
    
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
    
    # select on workflow
    L = list(map(lambda x: x.lower(), workflows))
    if 'bcl2fastq' not in L:
        return workflow.lower() in L
    else:
        w = L + sequencing_workflows
        return workflow.lower() in w
    
    
def is_correct_case(case, cases):
    '''
    (str, list) -> bool
       
    Returns True if case is a member of cases
        
    Parameters
    ----------
    - case (str): Case identifier
    - cases (list): List of cases     
    '''
    
    return case in cases
        
    
    
def is_correct_run(samples, runs):
    '''
    (dict, list) -> bool
    
    Returns True if all the runs in the samples dictionary are included
    in the list of runs to release
    
    Parameters
    ----------
    - samples (dict): Dictionary with sample information for a given file
    - runs (list): List of runs to release
    '''    
    
    L = [d['run'] for d in samples]
    
    return set(L).intersection(set(runs)) == set(L)
    
 
    
def map_libraries_to_file(samples):
    '''
    (dict) -> dict
    
    Returns a dictionary of libraries and runs mapped to a given file
    
    Parameters
    ----------
    - samples (dict): Dictionary with sample information for a given file
    '''
    
    D = {}
    
    for d in samples:
        library = d['library']
        run = d['run']
        lane = d['lane']
        if library not in D:
            D[library] = {}
        if run not in D[library]:
            D[library][run] = []
        D[library][run].append(lane)
        D[library][run].sort()
    
    return D
    
 
    
def is_correct_library(samples, valid_libraries):
    '''
    (dict, dict) -> bool
    
    Returns True if the libraries, runs and lanes (if specified) of a file are flagged to release 
    
    Parameters
    ----------
    - samples (dict): Dictionary with sample information for a given file
    - valid_libraries (dict): Dictionary with libraries tagged for release
    '''
    
    keep = [] 
     
    libraries = map_libraries_to_file(samples)
    for library in libraries:
        if library not in valid_libraries:
            keep.append(False)
        else:
            for run in libraries[library]:
                if run not in valid_libraries[library]:
                    keep.append(False)
                else:
                    for lane in libraries[library][run]:
                        # check if lane is specified in valid libraries
                        if '' not in valid_libraries[library][run]:
                            if lane not in valid_libraries[library][run]:
                                keep.append(False)
    
    return all(keep)
    

def get_analysis_files(analysis_pipeline):
    '''
    (str) -> list

    Returns a list of files generated by the analysis workflows 
    Pre-condition: The data is organized by case, workflow and files
    
    Parameters
    ----------
    - analysis_pipeline (str): Path to the json file with pipeline analysis data
    '''
    
    infile = open(analysis_pipeline)
    data = json.load(infile)
    infile.close()
    
    L = []
    for case in data:
        for workflow in data[case]:
            L.extend(data[case][workflow])
    
    return L
    
    
def select_files(file_info, project, workflows=None, runs = None, cases = None, libraries=None, release_files=None):
    '''
    (dict, str, list, list | None, list | None, dict | None, list | None) -> dict
     
    Restrict the file information of the files to release based
     
    Parameters
    ----------
    - file_info (dict): Dictionary with file information of all files in a case
    - project (str): Project of interest
    - workflows (list | None): List of workflows generating the files to release
    - runs (list | None): List of sequencing runs of generating the underlying sequences of the data 
    - cases (List | None): List of case identifiers
    - libraries (dict | None): Dictionary with libraries and runs for which data is to be released
    - release_files (list | None): List of file names or file paths to release 
    '''

    # selecting on files takes precedence, but files need to belong to project
    if release_files:
        # keep only files in the list of release files
        # check if the files are file names or full paths
        if all(map(lambda x: os.path.dirname(x) == '', release_files)):
            # evaluating file names
            to_remove = [i for i in file_info if os.path.basename(i) not in release_files or is_correct_project(file_info[i]['project'], project) == False]
        elif all(map(lambda x: os.path.dirname(x), release_files)):
            # dereference links
            release_files = list(map(lambda x: os.path.realpath(x), release_files))
            # evaluating full paths
            to_remove = [i for i in file_info if i not in release_files or is_correct_project(file_info[i]['project'], project) == False]    
    else:
        # select based on other options
        to_remove = []
        for file in file_info:
            # select the correct project
            L = [is_correct_project(file_info[file]['project'], project)]
            # select the correct workflow
            if workflows:
                L.append(is_correct_workflow(file_info[file]['workflow'], workflows))
            # select based on the list of cases
            if cases:
                L.append(is_correct_case(file_info[file]['case_id'], cases))
            # select based on the list of runs
            if runs:
                L.append(is_correct_run(file_info[file]['samples'], runs))
            # select based on libraries and runs and optionally on lanes
            if libraries:
                L.append(is_correct_library(file_info[file]['samples'], libraries))
            if all(L) == False: 
                to_remove.append(file)
                
    if to_remove:
        for i in to_remove:
            del file_info[i]

    return file_info


    
def update_file_info(D, file_info):
    '''
    (dict, dict) -> dict
    
    Returns an updated dictionary with the file information of a given case
    
    Parameters
    ----------
    - D (dict):
    - file_info (dict):
    '''
    
    for file in file_info:
        assert file not in D
        D[file] = file_info[file]
    
    return D



def extract_data(provenance_data, project, workflows=None, runs=None, cases=None, libraries=None, release_files=None):
    '''
    (list, str, list, list | None, list | None, dict | None, list | None) -> dict
    
    Returns a dictionary of file information for files selected for release
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries with case information
    - project (str): Project of interest
    - workflows (list | None): List of workflows generating the data to release
    - runs (list | None): List of sequencing runs
    - cases (list | None): List of case identifiers for which data need to be released
    - libraries (dict | None): Dictionary with libraries, runs and optional lanes 
    - release_files (list | None): List of files to release
    '''
    
    D = {}
    
    for case_data in provenance_data:
        case_id = case_data['case'] 
        # extract file info, sample info and map samples to files
        file_info = extract_file_info(case_data)
        sample_info = extract_sample_info(case_data)
        file_info = add_sample_info(file_info, sample_info)
        # select files
        file_info = select_files(file_info, project, workflows, runs, cases, libraries, release_files)        
        # update dict 
        if file_info:
            assert case_id not in D
            D[case_id] = file_info
            
    return D                  
                    

      
def create_working_dir(project, project_dir, project_name=None):
    '''    
    (str, str, str | None) -> str
    
    Creates and returns path to a sub-directory in project_dir named either after project or after project_name if defined
        
    Parameters
    ----------
    - project (str): Project name as it appears in File Provenance Report.
    - project_dir (str): Absolute or relative path to the project directory.
    - project_name (str | None): Project name used to create the project directory in gsi space
    '''

    # check if project dir is abolsute or relative path
    if os.path.isabs(project_dir):
        projectdir = project_dir
    else:
        projectdir = os.path.join(os.getcwd(), project_dir)
        projectdir = os.path.realpath(projectdir)

    # use project as project name if not specified
    if project_name:
        name = project_name
    else:
        name = project

    working_dir = os.path.join(projectdir, name)
    os.makedirs(working_dir, exist_ok=True)
    
    return working_dir
        

def is_case_info_incomplete(case_data):
    '''
    (dict) -> bool
    
    Returns True if the case information is complete
    
    Parameters
    ----------
    - case_data (dict): Dictionary with case information from production
    '''
    
    incomplete = [len(case_data[i]) == 0 for i in case_data]
    return any(incomplete)
    

def clean_up_provenance(provenance_data):
    '''
    (list) -> list, list
    
    Returns a list of dictionaries removing cases for which some information is not defined
    
    Parameters
    ----------
    - provenance_data (list): List of dictionaries with production data for cases
    '''    
    
    to_remove = [i for i in provenance_data if is_case_info_incomplete(i)]
    for i in to_remove:
        provenance_data.remove(i)
    
    return provenance_data, to_remove


def write_md5sum(file_info, outputfile):
    '''
    (dict, str) -> None   
    
    Writes the md5sum of each file in dictionary file info
        
    Parameters
    ----------
    - file_info (dict): Dictionary with information about the released files
    - outpufile (str): Path to the outputfile with md5sums
    '''    
    
    L = [] 
    for case_id in file_info:
        for file in file_info[case_id]:
            md5 = file_info[case_id][file]['md5sum']
            L.append([file, md5])
    L = list(map(lambda x: '\t'.join(x), L))
    L = sorted(list(set(L)))    
            
    newfile = open(outputfile, 'w')
    for i in L:
        newfile.write(i + '\n')
    newfile.close()


def group_sample_info_mapping(file_info):
    '''
    (dict) -> dict    
    
    Returns a dictionary of sample information organized by run
        
    Parameters
    ----------
    - file_info (dict) : Dictionary with file information for each file in case
    '''

    # group info by run id
    D = {}
    
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
    
    for case_id in file_info:
        for file in file_info[case_id]:
            workflow = file_info[case_id][file]['workflow']
            assert case_id == file_info[case_id][file]['case_id']
            if workflow.lower() in sequencing_workflows:
                for d in file_info[case_id][file]['samples']:
                    run = d['run']
                    sample_id = d['sample_id']
                    library = d['library']
                    barcode = d['barcode']                
                    external_id = d['external_id']    
                    donor= d['donor']
                    group_id = d['group_id']
                    library_design = d['library_design']
                    tissue_origin = d['tissue_origin']
                    tissue_type = d['tissue_type']
                    lane = d['lane']  
                    run = run + '_' + lane

                    # group files by sample, library, run and lane
                    key = '_'.join([case_id, donor, sample_id, library, run, barcode])
            
                    L = [case_id, donor, external_id, library, library_design, tissue_type, tissue_origin, run, barcode, group_id]
        
                    if key not in D:
                        D[key] = {'info': L, 'files': [file]}
                    else:
                        assert D[key]['info'] == L
                        D[key]['files'].append(file)
        
    return D            


def write_sample_map(project, file_info, working_dir):
    '''
    (str, dict, str) -> str
    
    Writes a sample map in the working directory about the released files 
    for the project of interest and returns the path to the sample map
    
    Parameters
    ----------
    - project (str): Project of interest
    - file_info (dict) : Dictionary with file information for each file in case
    - working_dir (str): Directory where the sample map is written
    '''

    # group sample information by run            
    sample_info = group_sample_info_mapping(file_info)        

    # write sample maps
    current_time = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
    outputfile = os.path.join(working_dir, '{0}.release.{1}.{2}.map.tsv'.format(project, current_time, 'fastqs'))
    newfile = open(outputfile, 'w')
    header = ['oicr_case', 'oicr_donor', 'donor_id', 'library_id', 'library_type', 'tissue_type', 'tissue_origin', 'run', 'barcode', 'sample_id', 'files']
    
    newfile.write('\t'.join(header) + '\n')
    # make a list of sorted samples
    keys = sorted(list(sample_info.keys()))
    for k in keys:
        files = ';'.join(sorted(list(map(lambda x: os.path.basename(x), sample_info[k]['files']))))
        info = sample_info[k]['info']
        info.append(files)        
        newfile.write('\t'.join(list(map(lambda x: str(x), info))) + '\n')
    newfile.close()

    return outputfile


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
    
    if library_file:
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
                lane = i[-1]
        
            if library not in D:
                D[library] = {}
            if run not in D[library]:
                D[library][run] = []
            D[library][run].append(lane)
            D[library][run] = sorted(list(set(D[library][run])))       
        
    return D



def get_release_files(release_files):
    '''
    (str) -> list
    
    Returns a list of file names or file paths to be released
    Pre-condition: The list should contain only file paths or only file names
    
    Parameters
    ----------
    - release_files (str): Path to the file containing the names or paths of the files to release
    '''
    
    L = []
    
    if release_files:
        infile = open(release_files)
        L = infile.read().rstrip().split('\n')
        # check that files are all full paths or file names 
        if all(map(lambda x: os.path.dirname(x) == '', L)) or all(map(lambda x: os.path.dirname(x), L)) == False:
            raise ValueError('Expecting only full paths or only file names')
              
    return L


def generate_links(file_info, project_dir):
    '''
    (dict, str) -> None
    
    Link files in case and workflow subdirectories of the project directory
        
    Parameters
    ----------
    - file_info (dict): Dictionary with file information
    - project_dir (str): Path to the project directory in GSI space  
    '''
    
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
    
    for case_id in file_info:
        for file in file_info[case_id]:
            assert case_id == file_info[case_id][file]['case_id']
            workflow = file_info[case_id][file]['workflow']
            case_dir = os.path.join(project_dir, case_id)
            os.makedirs(case_dir, exist_ok=True)
            donor = list(set([d['donor'] for d in file_info[case_id][file]['samples']]))
            assert len(donor) == 1
            donor = donor[0]
            donor_dir = os.path.join(case_dir, donor)
            os.makedirs(donor_dir, exist_ok=True)
            if workflow in sequencing_workflows:
                workflow_dir = os.path.join(donor_dir, 'sequences')
                os.makedirs(workflow_dir, exist_ok=True)
                run = list(set([d['run'] for d in file_info[case_id][file]['samples']]))
                assert len(run) == 1
                run = run[0]
                run_dir = os.path.join(workflow_dir, run + '.fastqs')
                os.makedirs(run_dir, exist_ok=True)
                filename = os.path.basename(file)
                link = os.path.join(run_dir, filename)
            else:
                workflow_dir = os.path.join(donor_dir, workflow)
                os.makedirs(workflow_dir, exist_ok=True)
                filename = os.path.basename(file)
                link = os.path.join(workflow_dir, filename)
            if os.path.isfile(link) == False:
                os.symlink(file, link)




def change_nabu_status(api, file_swids, qc_status, user_name, ticket):
    '''
    (str, list, str, str, str) -> None
    
    Modifies the file qc status in Nabu to qc_status, the username to user_name
    and comment to ticket if used 
    
    Parameters
    ----------
    - api (str): URL of the nabu API
    - file_swids (list): List of file unique identifiers
    - qc_status (str): File QC status: PASS, PENDING or FAIL
    - ticket (str): Jira ticket of the release
    '''
    
    # get end-point
    api += 'add-fileqcs' if api[-1] == '/' else '/add-fileqcs'
        
    if qc_status not in ['PASS', 'PENDING', 'FAIL']:
        raise ValueError('QC status is PASS, FAIL or PENDING')
    
    headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
    json_data = {'fileqcs': []}
    for file_swid in file_swids:
        d = {'fileid': file_swid, 'qcstatus': qc_status, 'username': user_name, 'comment': ticket}    
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



def list_files(directory):
    '''
    (str) -> list
    
    Returns a list of files contained in directory and all its subdirectories.
    If the files are links, the source file is returned
    
    Parameters
    - directory (str): Path to a directory
    '''
    
    L = os.walk(directory)    
    
    files = []
    
    for i in L:
        parentdir = i[0]
        # check if files are in parent directory
        if i[-1]:
            # list all the files
            files.extend([os.path.join(parentdir, j) for j in i[-1]])
    
    # get the real path of the files
    files = list(map(lambda x: os.path.realpath(x), files))

    return files





def count_released_fastqs_by_instrument(file_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with fastq files organized by instrument and run
        
    Parameters
    ----------
    - file_info (dict): Dictionary with file information for a case
    '''
        
    sequencing_workflows = ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport', 'import_fastq']
        
    D = {}
    
    for case_id in file_info:
        for file in file_info[case_id]:
            assert file_info[case_id][file]['workflow'].lower() in sequencing_workflows
            for d in file_info[case_id][file]['samples']:
                run = d['run']
                instrument = d['instrument']
                if instrument not in D:
                    D[instrument] = {}
                if run not in D[instrument]:
                    D[instrument][run] = []
                D[instrument][run].append(file)
                D[instrument][run] = list(set(D[instrument][run]))
    
    return D
    




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
                        
        sample = i['Pinery Lims ID']
        assert sample not in D
        D[sample] = d
    return D





def merge_bamqc_dnaseqc(bamqc_info, dnaseqqc_info):
    '''
    (dict, dict) -> None
    
    Merge the QC information extracted from the bamqc and dnaseqqc databases.
    Note that records in dnseqqc have precedence over records in bamqc in case of duplicates
        
    Parameters
    ----------
    - bamqc_info (dict): QC information for each paired fastq from the bamqc db
    - dnaseqqc_info (dict): QC information for each paired fastq from the dnaseqqc db
    '''
    
    D = {}
    
    # dnaseq qc has precedence over bam qc in case of duplicate records
    for sample in dnaseqqc_info:
        D[sample] =  dnaseqqc_info[sample]
       
    for sample in bamqc_info:
        if sample not in D:
            D[sample] = bamqc_info[sample]
    
    return D        


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
        sample = i['Pinery Lims ID']
        d = {}
        for j in columns:
            try:
                float(i[j])
            except:
                d[j] = i[j]
            else:
                d[j] = round(float(i[j]), 3)
        assert sample not in D
        D[sample] = d
                  
    return D


def extract_rnaseqqc_data(rnaseqqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with project-level relevant information from the rnaseqqc table of qc-etl
        
    Parameters
    ----------
    - rnaseqqc_db (str): Path to the rnaseq SQLIte database generated by qc-etl
    '''
    
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
        sample = i['Pinery Lims ID']
        d = {}
        for j in columns:
            try:
                float(i[j])
            except:
                d[j] = i[j]
            else:
                d[j] = round(float(i[j]), 3)
        assert sample not in D
        D[sample] = d
           
    return D



def extract_emseqqc_data(emseqqc_db):
    '''
    (str) -> dict
    
    Returns a dictionary with emseqqc from the qc-etl 
        
    Parameters
    ----------
    - emseqqc_db (str): Path to the emseqqc SQLIte database generated by qc-etl
    '''
    
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
        sample = i['Pinery Lims ID']
        d = {}
        for j in columns:
            try:
                float(i[j])
            except:
                d[j] = i[j]
            else:
                d[j] = round(float(i[j]), 3)
        assert sample not in D
        D[sample] = d           
    
    
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
        assert sample in D
        assert D[sample]['Barcodes'] == barcodes and D[sample]['Run Alias'] == run
        assert D[sample]['Lane Number'] == lane
        D[sample]['mark duplicates_PERCENT_DUPLICATION'] = duplicate    

    return D


def add_cfmedipqc_metrics(file_info, case_id, file, cfmedipqc_info):
    '''
    (dict, str, str, dict) -> None
    
    Update the file information in place with the QC information
    collected from cfmedipqc for a given file if library source is CM
    
    Parameters
    ----------
    - file_info (dict): Dictionary with information about the released fastqs
    - case_id (str): Case identifier
    - file (str): File path
    - cfmedipqc_info (dict): QC information for each paired fastq from the cfmedip QC db
    '''
    
    qc_found = False
    assert len(file_info[case_id][file]['samples']) == 1
    run = file_info[case_id][file]['samples'][0]['run']
    limskey = file_info[case_id][file]['samples'][0]['lims_id']
    barcode = file_info[case_id][file]['samples'][0]['barcode']     
    lane = file_info[case_id][file]['samples'][0]['lane']
    library_design = file_info[case_id][file]['samples'][0]['library_design']
    read_count = int(file_info[case_id][file]['attributes']['read_count'][0])
        
    if library_design ==' CM' and limskey in cfmedipqc_info:
        assert limskey == cfmedipqc_info[limskey]['Pinery Lims ID']
        assert lane == cfmedipqc_info[limskey]['Lane Number']
        assert barcode == cfmedipqc_info[limskey]['Barcodes']
        assert run == cfmedipqc_info[limskey]['Run Alias']
        qc_found = True
        assert 'metrics' not in file_info[case_id][file]
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      'AT_dropout': cfmedipqc_info[limskey]['AT Dropout'],
                                      'methylation_beta': cfmedipqc_info[limskey]['Methylation beta'],
                                      'duplication': cfmedipqc_info[limskey]['Percent Duplication'],
                                      'CpG_enrichment': cfmedipqc_info[limskey]['Relative CpG Frequency in Regions']}
    
    if library_design == 'CM' and qc_found == False:
        assert 'metrics' not in file_info[case_id][file]
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      'AT_dropout': 'NA',
                                      'methylation_beta': 'NA',
                                      'duplication': 'NA',
                                      'CpG_enrichment': 'NA'}
    
    


def add_rnaseqqc_metrics(file_info, case_id, file, rnaseqqc_info):
    '''
    (dict, str, str, dict) -> None
    
    Update the file information in place with the QC information
    collected from rnaseqqc for a given file if library source is WT
    
    Parameters
    ----------
    - file_info (dict): Dictionary with information about the released fastqs
    - case_id (str): Case identifier
    - file (str): File path
    - rnaseqqc_info (dict): QC information for each paired fastq from the rnaseqqc QC db
    '''
    
    
    qc_found = False
    assert len(file_info[case_id][file]['samples']) == 1
    run = file_info[case_id][file]['samples'][0]['run']
    limskey = file_info[case_id][file]['samples'][0]['lims_id']
    barcode = file_info[case_id][file]['samples'][0]['barcode']     
    lane = file_info[case_id][file]['samples'][0]['lane']
    library_design = file_info[case_id][file]['samples'][0]['library_design']
    read_count = int(file_info[case_id][file]['attributes']['read_count'][0])
        
    # check that limskey in recorded in rnaseqqc_db
    if library_design == 'WT' and limskey in rnaseqqc_info:
        assert limskey == rnaseqqc_info[limskey]['Pinery Lims ID']
        assert lane == rnaseqqc_info[limskey]['Lane Number']
        assert barcode == rnaseqqc_info[limskey]['Barcodes']
        assert run == rnaseqqc_info[limskey]['Run Alias']
        qc_found = True
        assert 'metrics' not in file_info[case_id][file]
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      "5'-3' bias": rnaseqqc_info[limskey]['MEDIAN_5PRIME_TO_3PRIME_BIAS'],
                                      'rRNA contamination': round((rnaseqqc_info[limskey]['rrna contamination properly paired'] / rnaseqqc_info[limskey]['rrna contamination in total (QC-passed reads + QC-failed reads)'] * 100), 3),
                                      'Coding (%)': rnaseqqc_info[limskey]['PCT_CODING_BASES'],
                                      'Correct strand reads (%)': rnaseqqc_info[limskey]['PCT_CORRECT_STRAND_READS']}
    
    if library_design == 'WT' and qc_found == False:
        assert 'metrics' not in file_info[case_id][file]    
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      "5'-3' bias": 'NA',
                                      'rRNA contamination': 'NA',
                                      'Coding (%)': 'NA',
                                      'Correct strand reads (%)': 'NA'}
    
        

def add_bamqc_metrics(file_info, case_id, file, bamqc_info):
    '''
    (dict, str, str, dict) -> None
    
    Update the file information with QC information collected from bamqc for a given file 
           
    Parameters
    ----------
    - file_info (dict): Dictionary with information about the released fastqs
    - case_id (str): Case identifier
    - file (str): File path
    - bamqc_info (dict): QC information for each paired fastq from the bamqc table
    '''
    
    qc_found = False
    assert len(file_info[case_id][file]['samples']) == 1
    run = file_info[case_id][file]['samples'][0]['run']
    limskey = file_info[case_id][file]['samples'][0]['lims_id']
    barcode = file_info[case_id][file]['samples'][0]['barcode']     
    lane = file_info[case_id][file]['samples'][0]['lane']
    library_design = file_info[case_id][file]['samples'][0]['library_design']
    sample_id = file_info[case_id][file]['samples'][0]['sample_id']
    group_id = file_info[case_id][file]['samples'][0]['group_id']
    library = file_info[case_id][file]['samples'][0]['library']
    instrument = file_info[case_id][file]['samples'][0]['instrument']
    read_count = int(file_info[case_id][file]['attributes']['read_count'][0])
        
    excluded_libraries = ['CM', 'WT', 'MC', 'MG']
    
    # check that limskey in recorded in bamqc
    if library_design not in excluded_libraries and limskey in bamqc_info:
        assert limskey == bamqc_info[limskey]['Pinery Lims ID']
        assert lane == bamqc_info[limskey]['Lane Number']
        assert barcode == bamqc_info[limskey]['Barcodes']
        assert run == bamqc_info[limskey]['Run Alias']
        assert '_'.join(sample_id.split('_')[:2]) == '_'.join(bamqc_info[limskey]['sample'].split('_')[:2])
        if group_id == 'NA':
            assert ('_'.join(sample_id.split('_')[:-1]) == bamqc_info[limskey]['sample'] or library == bamqc_info[limskey]['sample'])
        assert instrument  == bamqc_info[limskey]['instrument']
        qc_found = True
        assert 'metrics' not in file_info[case_id][file]
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      'coverage': round(bamqc_info[limskey]['coverage'], 2),
                                      'coverage_dedup': round(bamqc_info[limskey]['coverage deduplicated'], 2),
                                      'on_target': round(bamqc_info[limskey]['on_target'], 2),
                                      'percent_duplicate': round(bamqc_info[limskey]['mark duplicates_PERCENT_DUPLICATION'], 2)}
        
    if library_design not in excluded_libraries and qc_found == False:
        assert 'metrics' not in file_info[case_id][file]
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      'coverage': 'NA',
                                      'coverage_dedup': 'NA',
                                      'on_target': 'NA',
                                      'percent_duplicate': 'NA'}



def add_emseqqc_metrics(file_info, case_id, file, emseqqc_info):
    '''
    (dict, str, str, dict) -> None
      
    Update the file information with QC information collected from emseqqc for a given file 
           
    Parameters
    ----------
    - file_info (dict): Dictionary with information about the released fastqs
    - case_id (str): Case identifier
    - file (str): File path
    - emseqqc_info (dict): QC information for each paired fastq from the emseq QC db
    '''
        
    qc_found = False
    assert len(file_info[case_id][file]['samples']) == 1
    run = file_info[case_id][file]['samples'][0]['run']
    limskey = file_info[case_id][file]['samples'][0]['lims_id']
    barcode = file_info[case_id][file]['samples'][0]['barcode']     
    lane = file_info[case_id][file]['samples'][0]['lane']
    library_design = file_info[case_id][file]['samples'][0]['library_design']
    read_count = int(file_info[case_id][file]['attributes']['read_count'][0])
    
    # check that limskey has qc
    if library_design in ['MC', 'MG'] and limskey in emseqqc_info:
        assert limskey == emseqqc_info[limskey]['Pinery Lims ID']
        assert lane ==  emseqqc_info[limskey]['Lane Number']
        assert barcode == emseqqc_info[limskey]['Barcodes']
        assert run == emseqqc_info[limskey]['Run Alias']
        qc_found = True
        assert 'metrics' not in file_info[case_id][file]
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      'Lambda_methylation': emseqqc_info[limskey]['Lambda'],
                                      'pUC19_methylation': emseqqc_info[limskey]['pUC19'],
                                      'percent_duplication': emseqqc_info[limskey]['mark duplicates_PERCENT_DUPLICATION']}
        
    if library_design in ['MC', 'MG'] and qc_found == False:
        assert 'metrics' not in file_info[case_id][file]
        file_info[case_id][file]['metrics'] = {'read_count': read_count,
                                      'Lambda_methylation': 'NA',
                                      'pUC19_methylation': 'NA',
                                      'percent_duplication': 'NA'}
          
        

def add_read_count(file_info, case_id, file):
    '''
    (dict, str, str) -> None
      
    Records read count metrics  
           
    Parameters
    ----------
    - file_info (dict): Dictionary with information about the released fastqs
    - case_id (str): Case identifier
    - file (str): File path
    '''
        
    read_count = int(file_info[case_id][file]['attributes']['read_count'][0]) 
    
    if 'metrics' not in file_info[case_id][file]:
        file_info[case_id][file]['metrics'] = {}
    file_info[case_id][file]['metrics']['read_count'] = read_count
    


def add_QC_metrics(file_info, bamqc_info, cfmedipqc_info, rnaseqqc_info, emseqqc_info):
    '''
    (dict, dict, dict, dict, dict) -> None
    
    Update the file information with information collected from the appropriate library source-specific QC db
    
    Parameters
    ----------
    - file_info (dict): File information for fastqs 
    - bamqc_info (dict): QC information for each paired fastq from the bamqc db
    - cfmedipqc_info (dict): QC information for each paired fastq from the cfmedipqc db
    - rnaseqqc_info (dict) QC information for each paired fastqs from the rnaseqqc db
    - emseqqc_info (dict) QC information for each paired fastqs from the emseqqc db
    '''
    
    for case_id in file_info:
        for file in file_info[case_id]:
            assert len(file_info[case_id][file]['samples']) == 1
            library_design = file_info[case_id][file]['samples'][0]['library_design']
            if library_design == 'CM':
                add_cfmedipqc_metrics(file_info, case_id, file, cfmedipqc_info)
            elif library_design == 'WT':
                add_rnaseqqc_metrics(file_info, case_id, file, rnaseqqc_info)
            elif library_design in ['WG', 'EX', 'TS', 'PG']:
                add_bamqc_metrics(file_info, case_id, file, bamqc_info)
            elif library_design in ['MC', 'MG']:
                add_emseqqc_metrics(file_info, case_id, file, emseqqc_info)              
            else:
                add_read_count(file_info, case_id, file)
            


def map_library_design_to_instrument(file_info):
    '''
    (dict) -> dict
    
    Returns a dictionary of library design and its sequencing instruments
        
    Parameters
    ----------
    - file_info (dict): Dictionary with information about the released fastqs
    '''
    
    # list all instruments for each library design
    library_designs = {}
    for case_id in file_info:
        for file in file_info[case_id]:
            for d in file_info[case_id][file]['samples']:
                instrument = d['instrument']
                library_design = d['library_design']
                if library_design not in library_designs:
                    library_designs[library_design] = []
                library_designs[library_design].append(instrument)
    
    for i in library_designs:
        library_designs[i] = sorted(list(set(library_designs[i])))
                                                          
    return library_designs            
        



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
    
    # injected data have a different file format
    if x is None:
        x = re.search("[\._]{1}[1-2][\._]+.*fastq.gz" , filename)
        
    assert x
    prefix = filename[:x.span()[0]]
    
    return prefix
     

def add_file_prefix(file_info):
    '''
    (dict) -> None

    Adds file prefix to each fastq file in the file info dictionary in place.
    Pre-condition: The dictionary only stores information for FASTQ files   
        
    Parameters
    ----------
    - file_info (dict): Dictrionary of fastq file information
    '''

    for case_id in file_info:
        for file in file_info[case_id]:
            prefix = get_file_prefix(file)
            file_info[case_id][file]['prefix'] = prefix


def find_fastq_pairs(file_info, platform):
    '''
    (dict, str) -> list
    
    Returns a list of 2-item lists with dictionary about file info of paired fastqs
    for a specific sequencing instrument.
    Pre-condition: All reads are paired-reads and it exists 2 fastqs for read 1 and read 2
    
    Parameters
    ----------
    - file_info (dict): Dictionary with fastq file information 
    = platform (str): Sequencing instrument
    '''
     
    found = []
    
    # make a list with dictionaries of file info
    L = []
    
    for case_id in file_info:
        for file1 in file_info[case_id]:
            for file2 in file_info[case_id]:
                if file1 not in found and file2 not in found:
                    if file1 != file2 and file_info[case_id][file1]['wfrun_id'] == file_info[case_id][file2]['wfrun_id']:
                        assert file_info[case_id][file1]['metrics'] == file_info[case_id][file2]['metrics']
                        assert file_info[case_id][file1]['prefix'] == file_info[case_id][file2]['prefix']
                        read_num = [[file1, int(file_info[case_id][file1]['attributes']['read_number'][0])],
                                    [file2, int(file_info[case_id][file2]['attributes']['read_number'][0])]]
                        read_num.sort(key=lambda x: x[1])
                        assert file_info[case_id][file1]['samples'][0]['instrument'] == file_info[case_id][file2]['samples'][0]['instrument']
                        if file_info[case_id][file1]['samples'][0]['instrument'] == platform:
                            L.append([{read_num[0][0]: file_info[case_id][read_num[0][0]]},
                                      {read_num[1][0]: file_info[case_id][read_num[1][0]]}])
                            found.extend([file1, file2])
                            break
    return L


def get_run_level_metrics(file_info, platform, library_source):
    '''
    (dict, str, str) -> (list)
    
    Returns a tuple with parallel lists of run-level metrics for a given instrument.
    Pre-condition: All reads are paired-reads and it exists 2 fastqs for read 1 and read 2
    
    Parameters
    ----------
    - file_info (dict): Dictionary with fastq file information  
    - platform (str): Sequencing platform
    - library_source (str): Type of library (eg: CM, WG, WT)
    '''
    
    # get the list of metrics of interest 
    metrics = get_library_metrics(library_source)
     
    # find fastq pairs
    L = find_fastq_pairs(file_info, platform)
    
    # make parallel lists for each metrics
    QC_metrics = [[] for i in range(len(metrics))]
    
    for i in L:
        # record only metrics from read 1 file
        file = list(i[0].keys())[0]
        if i[0][file]['samples'][0]['library_design'] == library_source:
            for j in range(len(metrics)):
                QC_metrics[j].append(i[0][file]['metrics'][metrics[j]])            
    
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


def generate_library_platform_figure(file_info, project, library_source, platform, colors, working_dir, height=16, width=13):
    '''
    (dict, str, str, str, str, int, int) -> str
    
    Generate a figure with metrics QC-etl for a given library type and sequencing platform
    and returns the path to the figure file
        
    Parameters
    ----------
    - file_info (dict): Dictionary with fastq file information with metrics extracted from qc-etl
    - project (str): Name of project
    - library_source (str): Type of library
    - platform (str): Sequencing platform
    - working_dir (str): Path to the folder where figure files are written
    - height (int): Height of the figure
    - width (int): Width of the figure
    '''
    
    # make lists with metrics for each instrument 
    QC_metrics = get_run_level_metrics(file_info, platform, library_source)
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
            
        # get the Y axis labels
        Y_axis = get_Y_axis_labels(library_source)
                
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


def generate_report_plots(file_info, project, library_designs, working_dir):
    '''
    (dict, str, dict, str) -> dict
    
    Returns a dictionary with figure paths for each library design and sequencing platorm

    Parameters
    ----------
    - file_info (dict): Dictionary with fastq file information with metrics extracted from qc-etl
    - project (str): Name of project
    - library_designs (dict): Dictionary instruments for each library design
    - working_dir (str): Path to the folder where figure files are written
    '''

    # generate a figure for each instrument and library source and keep track of figure files
    figure_files = {}
    for library_source in library_designs:
        colors = ['#00CD6C', '#AF58BA', '#FFC61E', '#009ADE']
        for platform in library_designs[library_source]:
            figure = generate_library_platform_figure(file_info, project, library_source, platform, colors, working_dir)
            if library_source not in figure_files:
                figure_files[library_source] = {}
            figure_files[library_source][platform] = figure

    return figure_files



def count_samples_with_missing_values(file_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with the number of samples with missing metric values for
    each library type and instrument
    
    Parameters
    ----------
    - file_info (dict): Dictionary with file information and QC metrics extracted from qc-etl
    '''
    
    # count the number of samples with missing values for each instrument and library type
    D = {}
    
    for case_id in file_info:
        for file in file_info[case_id]:
            platform = file_info[case_id][file]['samples'][0]['instrument']
            sample = file_info[case_id][file]['samples'][0]['external_id']
            library_design = file_info[case_id][file]['samples'][0]['library_design']
            metrics = get_library_metrics(library_design)
                       
            for i in metrics:
                if i in file_info[case_id][file]['metrics'] and file_info[case_id][file]['metrics'][i] == 'NA':
                    if library_design not in D:
                        D[library_design] = {}
                    if platform not in D[library_design]:
                        D[library_design][platform] = set()
                    D[library_design][platform].add(sample)
    
    for library_design in D:
        for platform in D[library_design]:
            D[library_design][platform] = len(list(D[library_design][platform]))
    
    for library_design in D:
        to_remove = [platform for platform in D[library_design] if D[library_design][platform] == 0]
        for platform in to_remove:
            del D[library_design][platform]
    to_remove = [library_design for library_design in D if len(D[library_design]) == 0]
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


def group_sample_metrics(file_info, table):
    '''
    (dict, str) -> dict 
    
    Group sample identifiers or metrics 
    
    Parameters
    ----------
    - file_info (dict): Information about fastq files with metrics extracted from QC etl
    - table (str): Table of interest in the report
    '''
    
    # make a list of instruments
    instruments = []
    for case_id in file_info:
        for file in file_info[case_id]:
            for d in file_info[case_id][file]['samples']:
                instruments.append(d['instrument'])
    instruments = list(set(instruments))
    
    # record sample metrics for each instrument
    D = {}
    
    # find pairs of fastqs
    for platform in instruments:
        pairs = find_fastq_pairs(file_info, platform)
        for i in pairs:
            file = list(i[0].keys())[0]
            
            library = i[0][file]['samples'][0]['library']
            library_design = i[0][file]['samples'][0]['library_design']
            prefix = i[0][file]['prefix']
            group_id = i[0][file]['samples'][0]['group_id']
            run =  i[0][file]['samples'][0]['run']
            case_id = i[0][file]['case_id']
            external_name = i[0][file]['samples'][0]['external_id']      
            sample = i[0][file]['samples'][0]['sample_id']     
            donor = i[0][file]['samples'][0]['donor']
            tissue_origin = i[0][file]['samples'][0]['tissue_origin']
            tissue_type = i[0][file]['samples'][0]['tissue_type']
            lane = i[0][file]['samples'][0]['lane']
            barcode = i[0][file]['samples'][0]['barcode']
            sequencing_run = '{0} lane_{1}_{2}'.format(run, lane, barcode)
            
            max_length = 20
                        
            # reformat prefix to fit the table column
            if len(prefix) >= max_length:
                prefix = fit_into_column(prefix)
            if len(library) >= max_length:
                library = fit_into_column(library)
            if len(sample) >= max_length:
                sample = fit_into_column(sample)
            if len(group_id) >= max_length:
                group_id = fit_into_column(group_id)
                                    
            if table == 'sample_identifiers':
                L = [case_id, library, donor, external_name, group_id, library_design, tissue_origin, tissue_type]
            elif table == 'qc_metrics':
                metrics = get_library_metrics(library_design)
                QC_metrics = []
                for metric in metrics:
                    if metric == 'read_count':
                        QC_metrics.append('{:,}'.format(i[0][file]['metrics'][metric]))
                    else:
                        QC_metrics.append(i[0][file]['metrics'][metric])
                L = [library, prefix]
                L.extend(QC_metrics)
                
            if library_design not in D:
                D[library_design] = {}
            if platform not in D[library_design]:
                D[library_design][platform] = []
            if L not in D[library_design][platform]:
                D[library_design][platform].append(L)
    # sort sample info on sample name
    for i in D:
        for j in D[i]:
            D[i][j].sort(key = lambda x: x[0])        
        
    return D            


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
         'MG': 'Methylation Detection (C to T) Whole Genome', '6B': 'Biomodal 6 base genomes',
         '5B': 'Biomodal 5 base genomes'}
       
    return D



def get_library_tissue_types(file_info):
    '''
    (dict) -> dict
    
    Returns a dictionary with tissue origin, type and library source definitions from MISO
    for each released library
    
    Parameters
    ----------
    - file_info (dict): Dictionary with fastq file information 
    '''

    D = {'Library Type': [], 'Tissue Type': [], 'Tissue Origin': []}

    # definitions are from the Configuration tab in MISO
    tissue_types = get_tissue_types()
    tissue_origin = get_tissue_origin()
    library_design = get_library_design()

    for case_id in file_info:
        for file in file_info[case_id]:
            for d in file_info[case_id][file]['samples']:
                D['Library Type'].append(d['library_design'])
                D['Tissue Type'].append(d['tissue_type'])
                D['Tissue Origin'].append(d['tissue_origin'])
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




def get_identifiers_appendix(file_info, report):
    '''
    (dict, str) -> list
    
    Returns a list with definitions of columns in the sample identifier table
    
    Parameters
    ----------
    - file_info (dict): Dictionary with fastq file information
    - report (str): cumulative or batch report
    '''

    # get the library type, tissue type and tissue origin 
    D = get_library_tissue_types(file_info)
    
    if report  == 'batch':
        L = ['OICR Case: OICR-generated case identifier',
            'Library Id: OICR-generated library identifier',
             'Case Id: OICR-generated donor identifier',
             'Donor Id: user supplied donor identifier',
             'Sample Id: user supplied sample, this distinguishes distinct samples of the same type from the same donor. If only one sample per donor is submitted the value may match the donor Id',
             'Library Type (LT): {0}'.format(D['Library Type']),
             'Tissue Origin (TO): {0}'.format(D['Tissue Origin']),
             'Tissue Type (TT): {0}'.format(D['Tissue Type'])]         
    elif report == 'cumulative':
        L = ['OICR Case Id: OICR-generated case identifier',
             'Donor Id: user supplied donor identifier',
             'OICR Sample Id: The OICR generated sample identifier. The sample Id is formed from the following: 1. Case Id, 2. Tissue Origin, 3. Tissue Type, 4. Library Type and 5. User supplied Sample Id',
             'Sample Id: user supplied sample, this distinguishes distinct samples of the same type from the same donor. If only one sample per donor is submitted the value may match the donor Id',
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



def create_nabu_signoff(cases, nabu_key_file, user_name, ticket, signoff_step_name, deliverable, deliverable_type, nabu='https://nabu.gsi.oicr.on.ca/case/sign-off'):
    '''
    (list, str, str) -> dict
    
    Creates a signoff record  in Nabu at signoff_step_name for the deliverable and deliverable_type,
    identifying the user and the ticket, for all the case. Returns the response as a json       
        
    Parameters
    ----------
    - cases (list): List of case identifiers
    - nabu_key_file (str): File storing the nabu API key
    - user_name (str):
    - ticket (str): 
    - signoff_step_name (str): 
    - deliverable (str):
    - deliverable_type (str):
    - nabu (str) URL to access the signoffs in Nabu
    '''
    
    infile = open(nabu_key_file)
    nabu_key = infile.read().rstrip()
    infile.close()
    
    headers = {
        'accept': 'application/json',
        'X-API-KEY': nabu_key,
        'Content-Type': 'application/json',
    }

    json_data = {'caseIdentifiers': cases,
                 'qcPassed': True,
                 'username': user_name,
                 'signoffStepName': 'RELEASE',
                 'deliverableType': 'Data Release',
                 'deliverable': deliverable,
                 'comment': ticket,
                 'release': True}


    response = requests.post(nabu, headers=headers, json=json_data)

    return response



def get_deliverables(file_info):
    '''
    (dict) -> dict
    
    Returns a list of deliverables for each case in file_info
    
    Parameters
    ----------
    - file_info (dict): Dictionary with file information
    '''
    
    deliverables = {}
    
    for case_id in file_info:
        for file in file_info[case_id]:
            assert case_id == file_info[file]['case_id']
            if case_id not in deliverables:
                deliverables[case_id] = []
            for d in file_info[file]['project']:
                v = list(d.values())
                for i in v:
                    deliverables[case_id].extend(i.split(','))
            deliverables[case_id] = list(set(deliverables[case_id]))

    return deliverables



def keep_files_for_deliverable(file_info, deliverable):
    '''
    (dict, str) -> dict
    
    Returns the dictionary of file information keeping only the files
    corresponding to the deliverable
    
    Parameters
    ----------
    - file_info (dict): Dictionary with file information
    - deliverable (str): Deliverable associated with the case (FastQ or Full Pipeline)
    '''
    
    # keep opnly files corresponding to deliverable
    fastqs, analysis_files = [], []
    for case_id in file_info:
        for file in file_info[case_id]:
            if 'fastq.gz' in file:
                fastqs.append(file)
            else:
                analysis_files.append(file)
             
    if deliverable == 'FastQ':
        print('{0}/{1} are fastq files'.format(len(fastqs), len(file_info)))
        if analysis_files:
            print('removing {0} analysis files'.format(len(analysis_files)))
            for i in analysis_files:
                for case_id in file_info:
                    if i in file_info[case_id]:
                        del file_info[case_id][i]
    elif deliverable == 'Full Pipeline':
        print('{0}/{1} are analysis files'.format(len(analysis_files), len(file_info)))
        if fastqs:
            print('removing {0} fastqs'.format(len(fastqs)))
            for i in fastqs:
                for case_id in file_info:
                    if i in file_info[case_id]:
                        del file_info[case_id][i]

    return file_info


def collect_qc_metrics(library_designs, bamqc_db, dnaseqqc_db, rnaseqqc_db, emseqqc_db, cfmedipqc_db):
    '''
    (list, str, str, str, str, str) -> list
    
    Returns a list of dictionaries with QC metrics extracted from the various caches
        
    Parameters
    ----------
    - bamqc_db (str): Path to the bamqc SQLite database
    - dnaseqqc_db (str): Path to the dnaseqqc SQLite database
    - rnaseqqc_db (str): Path to the rnaseq SQLite database
    - emseqqc_db (str): Path to the emseq SQLite database
    - cfmedipqc_db (str): Path to the cfmedip SQLite database
    '''
    
    
    if 'CM' in library_designs:
        # collect information from cfmedip table
        cfmedipqc_info = extract_cfmedipqc_data(cfmedipqc_db)
    else:
        cfmedipqc_info = {}


    if 'WT' in library_designs:
        # collect information from rnaseq table
        rnaseqqc_info = extract_rnaseqqc_data(rnaseqqc_db)
    else:
        rnaseqqc_info = {}

    if 'MC' in library_designs or 'MG' in library_designs:
        # collect information from emseq cache
        emseqqc_info = extract_emseqqc_data(emseqqc_db)    
    else:
        emseqqc_info = {}

    if 'WG' in library_designs or 'EX' in library_designs \
        or 'TS' in library_designs or 'PG' in library_designs:
            # collect information from bamqc table
            bamqc = extract_bamqc_data(bamqc_db)
            # collect information from dnaseqqc table
            dnaseqqc = extract_bamqc_data(dnaseqqc_db)
            # merge bamqc and dnaseqc
            bamqc_info = merge_bamqc_dnaseqc(bamqc, dnaseqqc)    
    else:
        bamqc_info = {}

    return bamqc_info, rnaseqqc_info, emseqqc_info, cfmedipqc_info



def link_files(args):
    '''
    (str, str, str, str, str | None, list | None, list | None, List | None, str | None, str | None ) -> None
    
    Links files to release in the project directory, with data organized by case and workflow name
        
    Parameters
    ----------
    - provenance (str): Path to json with production data. Default is
                        /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json
    - project (str): Project of interest
    - projects_dir (str): Parent directory containing the project subdirectories with file links.
                          Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/
    - project_name (str): Project name used to create the project directory in gsi space
    - libraries (str | None): File with libraries tagged for release.
                              The first column is always the library.
                              The second column is the run id.
                              The third optional column is the lane number.
    - workflows (list | None): List of workflows generating the data to release
    - runs (list | None): List of run Ids
    - cases (List | None): List of case Ids
    - release_files (str | None): Path to file with file names or paths to be released 
    - analyses (str | None): Path to the file with hierarchical structure storing sample and workflow ids
    '''
    
    # check options
    if args.release_files:
        a = [args.workflows, args.runs, args.cases, args.libraries, args.analyses]
        if any(a):
            c = ['-w', '-r', '-c', '-l', '-a']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-f cannot be used with options {0}'.format(err))
    else:
        if args.analyses:
            a = [args.release_files, args.workflows, args.runs, args.cases, args.libraries]
            if any(a):
                c = ['-f', '-w', '-r', '-c', '-l']
                err = ','.join([c[i] for i in range(len(c)) if a[i]])
                sys.exit('-a cannot be used with options {0}'.format(err))
        else:
            if not args.workflows:
                sys.exit('Use -w to indicate the pipeline workflows')
            if args.runs and args.libraries:
                sys.exit('-r and -l are exclusive parameters')    
        
    # create a working directory to link files and save md5sums 
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # load data
    provenance_data = load_data(args.provenance)
    print('loaded data')
    # clean up data
    provenance_data, deleted_cases = clean_up_provenance(provenance_data)
    print('removed {0} incomplete cases'.format(len(deleted_cases)))
    # extract data to release
    libraries = get_libraries(args.libraries)
    release_files = []
    if args.analyses:
        release_files = get_analysis_files(args.analyses)
    if args.release_files:
        release_files = get_release_files(args.release_files)
    file_info = extract_data(provenance_data, args.project, workflows=args.workflows, runs=args.runs, cases=args.cases, libraries=libraries, release_files=release_files)
    print('extracted data for {0} files'.format(len(file_info))) 

    # create sample map
    sample_map = write_sample_map(args.project, file_info, working_dir)
    print('wrote sample map {0}'.format(sample_map))
      
    # write summary md5sums
    # create a dictionary {run: [md5sum, file_path]}
    current_time = time.strftime('%Y-%m-%d_%H:%M', time.localtime(time.time()))
    outputfile = os.path.join(working_dir, '{0}.release.{1}.md5sums'.format(args.project, current_time))
    write_md5sum(file_info, outputfile)
    print('wrote md5sums: {0}'.format(outputfile))
    
    # link data
    generate_links(file_info, working_dir)
    print('linked data to {0}'.format(working_dir))
    
    
       
def map_external_ids(args):
    '''
    (str, str, str, str, str | None, list | None, list | None) -> None
        
    Generate sample maps with sample and sequencing information
    
    Parameters
    ----------
    - provenance (str): Path to json with production data. Default is
                        /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json
    - project (str): Project of interest
    - projects_dir (str): Parent directory containing the project subdirectories with file links.
                          Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/
    - project_name (str): Project name used to create the project directory in gsi space
    - libraries (str | None): File with libraries tagged for release.
                              The first column is always the library.
                              The second column is the run id.
                              The third optional column is the lane number.
    - runs (list | None): List of run Ids
    - cases (List | None): List of case Ids
    - release_files (str | None): File with file names or full paths of files to release
    - analyses (str | None): Path to the json file storing analysis data
    - directories (list | None): List of directories with links or files to mark in Nabu
    '''
    
    # check options
    if args.release_files:
        a = [args.runs, args.cases, args.libraries, args.analyses, args.directories]
        if any(a):
            c = ['-r', '-c', '-l', '-a', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-f cannot be used with options {0}'.format(err))
    elif args.analyses:
        a = [args.release_files, args.runs, args.cases, args.libraries, args.directories]
        if any(a):
            c = ['-f', '-r', '-c', '-l', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-a cannot be used with options {0}'.format(err))
    elif args.directories:
        # check that all directories are valid
        if all(list(map(lambda x: os.path.isdir(x), args.directories))) == False:
            sys.exit('Please provide valid directories with -d')        
        a = [args.release_files, args.runs, args.cases, args.libraries, args.analyses]
        if any(a):
            c = ['-f', '-r', '-c', '-l', '-a']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-d cannot be used with options {0}'.format(err))
    else:
        if args.runs and args.libraries:
           sys.exit('-r and -l are exclusive parameters')    
    
    # create a working directory to link files and save md5sums 
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # load data
    provenance_data = load_data(args.provenance)
    print('loaded data')
    # clean up data
    provenance_data, deleted_cases = clean_up_provenance(provenance_data)
    print('removed {0} incomplete cases'.format(len(deleted_cases)))
    
    
    if args.directories:
        # list of the linked files
        linked_files = []
        # list all files in directories and subdirectories
        # list path of target files if files are links
        for directory in args.directories:
            linked_files.extend(list_files(directory))
        # keep only fastq files
        if linked_files:
            linked_files = [i for i in linked_files if 'fastq.gz' in i]
            if linked_files:
                file_info = extract_data(provenance_data, args.project, release_files=linked_files)
            else:
                file_info = {}
    else:
        # extract data to release
        libraries = get_libraries(args.libraries)
        release_files = []
        if args.analyses:
            release_files = get_analysis_files(args.analyses)
        if args.release_files:
            release_files = get_release_files(args.release_files)
        # keep only the fastq files 
        if release_files:
            release_files = [i for i in release_files if 'fastq.gz' in i]
        file_info = extract_data(provenance_data, args.project, ['bcl2fastq'], runs=args.runs, cases=args.cases, libraries=libraries, release_files=release_files)
    print('extracted data for {0} files'.format(len(file_info))) 

    # create sample map
    sample_map = write_sample_map(args.project, file_info, working_dir)
    print('wrote sample map {0}'.format(sample_map))
      


def mark_files_nabu(args):
    '''
    (str, str, str, list|None, List|None,
     list|None, str|None, str|None, list|None,
     str|None, str, str, str) -> None
    
    Mark released files with user name and PASS and withheld files with user name and FAIL in Nabu

    Parameters
    ----------    
    - provenance (str): Path to json with production data. Default is
                        /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json
    - nabu (str): URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca
    - project (str): Project of interest
    - workflows (list | None): List of workflows generating the data to release
    - cases (List | None): List of case Ids
    - runs (list | None): List of run Ids
    - libraries (str | None): File with libraries tagged for release.
                              The first column is always the library.
                              The second column is the run id.
                              The third optional column is the lane number.
    - release_files (str | None): File with file names or full paths of files to release
    - directories (list | None): List of directories with links or files to mark in Nabu
    - analyses (str | None): Path to the json with pipeline data to release
    - status (str): New file QC status. Values are pass, fail or pending
    - user (str): User name to appear in Nabu for each released or whitheld file
    - ticket (str): Jira ticket 
    '''
    
    # check options
    if args.release_files:
        a = [args.workflows, args.runs, args.cases, args.libraries, args.analyses, args.directories]
        if any(a):
            c = ['-w', '-r', '-c', '-l', '-a', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-f cannot be used with options {0}'.format(err))
    elif args.analyses:
        a = [args.release_files, args.workflows, args.runs, args.cases, args.libraries, args.directories]
        if any(a):
            c = ['-f', '-w', '-r', '-c', '-l', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-a cannot be used with options {0}'.format(err))
    elif args.directories:
        # check that all directories are valid
        if all(list(map(lambda x: os.path.isdir(x), args.directories))) == False:
            sys.exit('Please provide valid directories with -d')        
        a = [args.release_files, args.workflows, args.runs, args.cases, args.libraries, args.analyses]
        if any(a):
            c = ['-f', '-w', '-r', '-c', '-l', '-a']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-d cannot be used with options {0}'.format(err))
    else:
        if not args.workflows:
            sys.exit('Use -w to indicate the pipeline workflows')
        if args.runs and args.libraries:
           sys.exit('-r and -l are exclusive parameters')    
    

    # get the files to mark
    # load data
    provenance_data = load_data(args.provenance)
    print('loaded data')
    # clean up data
    provenance_data, deleted_cases = clean_up_provenance(provenance_data)
    print('removed {0} incomplete cases'.format(len(deleted_cases)))
    
    if args.directories:
        # list of the linked files
        linked_files = []
        # list all files in directories and subdirectories
        # list path of target files if files are links
        for directory in args.directories:
            linked_files.extend(list_files(directory))
        if linked_files:
            file_info = extract_data(provenance_data, args.project, release_files=linked_files)
        else:
            file_info = {}
    else:
        # extract data to release
        libraries = get_libraries(args.libraries)
        release_files = []
        if args.analyses:
            release_files = get_analysis_files(args.analyses)
        if args.release_files:
            release_files = get_release_files(args.release_files)
        file_info = extract_data(provenance_data, args.project, args.workflows, runs=args.runs, cases=args.cases, libraries=libraries, release_files=release_files)
    print('extracted data for {0} files'.format(len(file_info))) 

    # make a list of swids
    if file_info:
        swids = []
        for case_id in file_info:
            for file in file_info[case_id]:
                swids.append(file_info[case_id][file]['accession'])
        swids = list(set(swids))
        # mark files in nabu
        change_nabu_status(args.nabu, swids, args.status.upper(), args.user, args.ticket)


def case_signoff(args):
    '''
    (str, str, str, list|None, List|None,
     list|None, str|None, str|None, list|None,
     str|None, str, str, str) -> None
    
    Signoff case record in Nabu for specfific deliverable
    
    Parameters
    ----------    
    - provenance (str): Path to json with production data. Default is
                        /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json
    - project (str): Project of interest
    - workflows (list | None): List of workflows generating the data to release
    - cases (List | None): List of case Ids
    - runs (list | None): List of run Ids
    - libraries (str | None): File with libraries tagged for release
    - release_files (str | None): File with file names or full paths of files to release
    - analyses (str | None): Path to the json file storing analysis data
    - directories (list | None): List of directories with links or files to mark in Nabu
    - user_name (str): Name of user
    - ticket (str): Jira ticket
    - nabu (str): URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca
    - nabu_key_file (str): Path to the file with the nabu key. Default is /.mounts/labs/gsi/secrets/nabu-prod_qc-gate-etl_api-key
    - deliverable (str): Deliverable. Choice is FastQ
    - signoff_step_name (str): Signoff step. Default is RELEASE
    - deliverable_type (str): Deliverable type. Default is Data Release
    '''
    
    # check options
    if args.release_files:
        a = [args.workflows, args.runs, args.cases, args.libraries, args.analyses, args.directories]
        if any(a):
            c = ['-w', '-r', '-c', '-l', '-a', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-f cannot be used with options {0}'.format(err))
    elif args.analyses:
        a = [args.release_files, args.workflows, args.runs, args.cases, args.libraries, args.directories]
        if any(a):
            c = ['-f', '-w', '-r', '-c', '-l', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-a cannot be used with options {0}'.format(err))
    elif args.directories:
        # check that all directories are valid
        if all(list(map(lambda x: os.path.isdir(x), args.directories))) == False:
            sys.exit('Please provide valid directories with -d')        
        a = [args.release_files, args.workflows, args.runs, args.cases, args.libraries, args.analyses]
        if any(a):
            c = ['-f', '-w', '-r', '-c', '-l', '-a']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-d cannot be used with options {0}'.format(err))
    else:
        if not args.workflows:
            sys.exit('Use -w to indicate the pipeline workflows')
        if args.runs and args.libraries:
           sys.exit('-r and -l are exclusive parameters')    
    

    # get the files to mark
    # load data
    provenance_data = load_data(args.provenance)
    print('loaded data')
    # clean up data
    provenance_data, deleted_cases = clean_up_provenance(provenance_data)
    print('removed {0} incomplete cases'.format(len(deleted_cases)))
    
    if args.directories:
        # list of the linked files
        linked_files = []
        # list all files in directories and subdirectories
        # list path of target files if files are links
        for directory in args.directories:
            linked_files.extend(list_files(directory))
        if linked_files:
            file_info = extract_data(provenance_data, args.project, release_files=linked_files)
        else:
            file_info = {}
    else:
        # extract data to release
        libraries = get_libraries(args.libraries)
        release_files = []
        if args.analyses:
            release_files = get_analysis_files(args.analyses)
        if args.release_files:
            release_files = get_release_files(args.release_files)
        file_info = extract_data(provenance_data, args.project, args.workflows, runs=args.runs, cases=args.cases, libraries=libraries, release_files=release_files)
        print('extracted data for {0} files'.format(len(file_info))) 

        # get end-point
        if args.nabu[-1] == '/':
            api = args.nabu + 'case/sign-off'
        else:
            api = args.nabu + '/case/sign-off'
    
    # list deliverables for each case
    cases = get_deliverables(file_info)
    print('identified {0} cases for {1} signoff'.format(len(cases), args.deliverable))
    # keep opnly files corresponding to deliverable
    file_info = keep_files_for_deliverable(file_info, args.deliverable)
    # list deliverables for each case
    cases = get_deliverables(file_info)
    print('keeping {0} cases for {1} signoff'.format(len(cases), args.deliverable))
    
    if cases:
        print('Signing off for {0} cases at step {1} for deliverable {2} of {3}'.format(len(cases), args.signoff_step_name, args.deliverable, args.deliverable_type))
        # list cases with expected deliverables
        signoff_cases = [i for i in cases if args.deliverable in cases[i]]
               
        # check that all expected cases have the correct deliverables
        if sorted(list(cases.keys())) == sorted(signoff_cases):
            # proceed to signoff
            response = create_nabu_signoff(signoff_cases, args.nabu_key_file, args.user_name, args.ticket, args.signoff_step_name, args.deliverable, args.deliverable_type, api)
            # check the response
            if response.ok:
                print('Signed off {0} for {1} deliverable at step {2} for cases: {3}'.format(args.deliverable_type,
                       args.deliverable, args.signoff_step_name, ','.join(cases)))
            else:
                print('Could not proceed with the signoff. Status code is {0}'.format(response.status_code))
        else:
            print('Only {0} cases out of {1} cases with released files have the expected deliverable'.format(len(signoff_cases), len(cases)))
            


def write_batch_report(args):
    '''
    (str, str, str, str, str, str, list, str, str | None, bool) -> None
    
    Write a PDF report with QC metrics and released fastqs for a given project

    Parameters
    ----------
    - project (str): Project name
    - project_name (str): Project name used to create the project directory in gsi space
    - projects_dir (str): Parent directory containing the project subdirectories with file links
    - project_full_name (str): Full name of the project
    - user (str): Name of the GSI personnel generating the report
    - ticket (str): Jira data release ticket code
    - provenance (str): Path to the json with production data
    - runs (list): List of run Ids
    - cases (list): List of case Ids
    - libraries (str): File with libraries tagged for release.
                       The first column is always the library.
                       The second column is the run id.
                       The third column optional is optional and is the lane
    - release_files (str): File with file names or full paths of files to release
    - analyses (str): Path to the json file storing analysis data
    - directories (list): List of directories with links or files
    - nabu (str): URL of the Nabu API
    - bamqc_db (str): Path to the bamqc SQLite database
    - dnaseqqc_db (str): Path to the dnaseqqc SQLite database
    - cfmedipqc_db (str): Path to the cfmedip SQLite database
    - rnaseqqc_db (str): Path to the rnaseq SQLite database
    - emseqqc_db (str): Path to the emseq SQLite database
    - keep_html (bool): Write html report if activated
    '''

    # check options
    if args.release_files:
        a = [args.runs, args.cases, args.libraries, args.analyses, args.directories]
        if any(a):
            c = ['-r', '-c', '-l', '-a', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-f cannot be used with options {0}'.format(err))
    elif args.analyses:
        a = [args.release_files, args.runs, args.cases, args.libraries, args.directories]
        if any(a):
            c = ['-f', '-r', '-c', '-l', '-d']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-a cannot be used with options {0}'.format(err))
    elif args.directories:
        # check that all directories are valid
        if all(list(map(lambda x: os.path.isdir(x), args.directories))) == False:
            sys.exit('Please provide valid directories with -d')        
        a = [args.release_files, args.runs, args.cases, args.libraries, args.analyses]
        if any(a):
            c = ['-f', '-r', '-c', '-l', '-a']
            err = ','.join([c[i] for i in range(len(c)) if a[i]])
            sys.exit('-d cannot be used with options {0}'.format(err))
    else:
        if args.runs and args.libraries:
           sys.exit('-r and -l are exclusive parameters')    
    
    # create working directory
    working_dir = create_working_dir(args.project, args.projects_dir, args.project_name)
    
    # get the files to mark
    # load data
    provenance_data = load_data(args.provenance)
    print('loaded data')
    # clean up data
    provenance_data, deleted_cases = clean_up_provenance(provenance_data)
    print('removed {0} incomplete cases'.format(len(deleted_cases)))
    
    if args.directories:
        # list of the linked files
        linked_files = []
        # list all files in directories and subdirectories
        # list path of target files if files are links
        for directory in args.directories:
            linked_files.extend(list_files(directory))
        # keep only fastq files
        if linked_files:
            linked_files = [i for i in linked_files if 'fastq.gz' in i]
            if linked_files:
                file_info = extract_data(provenance_data, args.project, release_files=linked_files)
            else:
                file_info = {}
                print('Expecting fastqs in directories, but found none')
        else:
            file_info = {}
            print('Expecting fastqs in directories, but found none')

    else:
        # extract data to release
        libraries = get_libraries(args.libraries)
        release_files = []
        if args.analyses:
            release_files = get_analysis_files(args.analyses)
        if args.release_files:
            release_files = get_release_files(args.release_files)
        # keep only the fastq files 
        if release_files:
            release_files = [i for i in release_files if 'fastq.gz' in i]
            if release_files:
                file_info = extract_data(provenance_data, args.project, ['bcl2fastq'], runs=args.runs, cases=args.cases, libraries=libraries, release_files=release_files)
            else:
                file_info = {}
                print('Expecting fastqs but found none')
        else:
            file_info = extract_data(provenance_data, args.project, ['bcl2fastq'], runs=args.runs, cases=args.cases, libraries=libraries, release_files=release_files)
        print('extracted data for {0} files'.format(len(file_info))) 

    # write sample map
    sample_map = write_sample_map(args.project, file_info, working_dir)
    print('wrote sample map {0}'.format(sample_map))
    
    # write summary md5sums
    # create a dictionary {run: [md5sum, file_path]}
    current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    md5sum_file = os.path.join(working_dir, '{0}.batch.release.{1}.md5'.format(args.project, current_time))
    write_md5sum(file_info, md5sum_file)
    print('wrote md5sums: {0}'.format(md5sum_file))
    # generate md5sum
    
    # count the number of fastq files by instrument and run
    fastq_counts = count_released_fastqs_by_instrument(file_info)
    all_released_files = sum([len(fastq_counts[instrument][run]) for instrument in fastq_counts for run in fastq_counts[instrument]])
    all_released_files = int(all_released_files / 2)

    #list all platforms for each library design
    library_designs = map_library_design_to_instrument(file_info)
    
    # collect QC metrics    
    bamqc_info, rnaseqqc_info, emseqqc_info, cfmedipqc_info = collect_qc_metrics(library_designs.keys(), args.bamqc_db, args.dnaseqqc_db, args.rnaseqqc_db, args.emseqqc_db, args.cfmedipqc_db)
    # update file info with QC metrics
    add_QC_metrics(file_info, bamqc_info, cfmedipqc_info, rnaseqqc_info, emseqqc_info)    
    print('collected QC metrics')
    
    # add file prefix to each fastq file
    add_file_prefix(file_info)
    
    # generate plots for each instrument and library source and keep track of figure files
    figure_files = generate_report_plots(file_info, args.project, library_designs, working_dir)
    print('generated figure files')
    
    # count the number of samples with missing metric values
    samples_missing_metrics = count_samples_with_missing_values(file_info)
    
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
    #header_identifiers = ['Library Id', 'Case Id', 'Donor Id', 'Sample Id', 'Sample Description', 'LT', 'TO', 'TT']
    
    header_identifiers = ['OICR case', 'Library Id', 'Case Id', 'Donor Id', 'Sample Id', 'LT', 'TO', 'TT']
        
    sample_identifiers = group_sample_metrics(file_info, 'sample_identifiers')
    appendix_identifiers = get_identifiers_appendix(file_info, 'batch')
    
    qc_metrics = group_sample_metrics(file_info, 'qc_metrics')
    library_sources = sorted(list(library_designs.keys()))
    header_metrics = {}
    for i in library_sources:
        header_metrics[i] = ['Library Id', 'File prefix'] + get_Y_axis_labels(i)    
    
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
               'library_designs': library_designs, 
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


    
if __name__ == '__main__':

    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'dare.py', description='A tool to manage data release')
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')
    
   	# link files in gsi space 
    l_parser = subparsers.add_parser('link', help="Link files to release")
    l_parser.add_argument('-p', '--parent', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    l_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    l_parser.add_argument('-pr', '--project', dest='project', help='Project name', required=True)
    l_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The second column is the run id and the third optional column is the lane number')
    l_parser.add_argument('-pv', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json', help='Path to the json with production data. Default is /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json')
    l_parser.add_argument('-w', '--workflows', dest='workflows', nargs='*', help='List of workflows')
    l_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run Ids')
    l_parser.add_argument('-f', '--files', dest='release_files', help='File with file names or full paths of files to release')
    l_parser.add_argument('-c', '--cases', dest='cases', nargs='*', help='List of case Ids')
    l_parser.add_argument('-a', '--analyses', dest='analyses', help='Path to the json file storing analysis data')
    l_parser.set_defaults(func=link_files)
    
   	# map external IDs 
    m_parser = subparsers.add_parser('map', help="Map files to external IDs")
    m_parser.add_argument('-p', '--parent', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    m_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    m_parser.add_argument('-pr', '--project', dest='project', help='Project name', required=True)
    m_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The second column is the run id and the third optional column is the lane number')
    m_parser.add_argument('-pv', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json', help='Path to the json with production data. Default is /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json')
    m_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run Ids')
    m_parser.add_argument('-c', '--cases', dest='cases', nargs='*', help='List of case Ids')
    m_parser.add_argument('-f', '--files', dest='release_files', help='File with file names or full paths of files to release')
    m_parser.add_argument('-a', '--analyses', dest='analyses', help='Path to the json file storing analysis data')
    m_parser.add_argument('-d', '--directories', dest='directories', nargs='*', help='List of directories with links or files to mark in Nabu')
    m_parser.set_defaults(func=map_external_ids)

    # mark files in nabu 
    qc_parser = subparsers.add_parser('qc', help="Updates FileQC status in Nabu")
    qc_parser.add_argument('-pv', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json', help='Path to the json with production data. Default is /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json')
    qc_parser.add_argument('-pr', '--project', dest='project', help='Project name', required=True)
    qc_parser.add_argument('-w', '--workflows', dest='workflows', nargs='*', help='List of workflows')
    qc_parser.add_argument('-c', '--cases', dest='cases', nargs='*', help='List of case Ids')
    qc_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run Ids')
    qc_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The second column is the run id and the third optional column is the lane number')
    qc_parser.add_argument('-f', '--files', dest='release_files', help='File with file names or full paths of files to release')
    qc_parser.add_argument('-a', '--analyses', dest='analyses', help='Path to the json file storing analysis data')
    qc_parser.add_argument('-d', '--directories', dest='directories', nargs='*', help='List of directories with links or files to mark in Nabu')
    qc_parser.add_argument('-nb', '--nabu', dest='nabu', default='https://nabu-prod.gsi.oicr.on.ca', help='URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca')
    qc_parser.add_argument('-st', '--status', dest='status', choices = ['fail', 'pass'], help='Mark files accordingly when released or withheld', required = True)
    qc_parser.add_argument('-u', '--user', dest='user', help='Name of user', required=True)
    qc_parser.add_argument('-t', '--ticket', dest='ticket', help='Ticket associated with the file QC change', required=True)
    qc_parser.set_defaults(func=mark_files_nabu)
    
    
    # signoff cases in nabu 
    s_parser = subparsers.add_parser('signoff', help="Create case signoff in Nabu")
    s_parser.add_argument('-pv', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json', help='Path to the json with production data. Default is /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json')
    s_parser.add_argument('-pr', '--project', dest='project', help='Project name', required=True)
    s_parser.add_argument('-w', '--workflows', dest='workflows', nargs='*', help='List of workflows')
    s_parser.add_argument('-c', '--cases', dest='cases', nargs='*', help='List of case Ids')
    s_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run Ids')
    s_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The second column is the run id and the third optional column is the lane number')
    s_parser.add_argument('-f', '--files', dest='release_files', help='File with file names or full paths of files to release')
    s_parser.add_argument('-a', '--analyses', dest='analyses', help='Path to the json file storing analysis data')
    s_parser.add_argument('-d', '--directories', dest='directories', nargs='*', help='List of directories with links or files to mark in Nabu')
    s_parser.add_argument('-u', '--user', dest='user_name', help='Name of user', required=True)
    s_parser.add_argument('-t', '--ticket', dest='ticket', help='Ticket associated with the file QC change', required=True)
    s_parser.add_argument('-nb', '--nabu', dest='nabu', default='https://nabu-prod.gsi.oicr.on.ca', help='URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca')
    s_parser.add_argument('-nk', '--nabu_key', dest='nabu_key_file', default='/.mounts/labs/gsi/secrets/nabu-prod_case-etl_api-key', help='Path to the file with the nabu key. Default is /.mounts/labs/gsi/secrets/nabu-prod_qc-gate-etl_api-key')
    s_parser.add_argument('-dv', '--deliverable', dest='deliverable', choices=['FastQ', 'Full Pipeline'], help='Deliverable', required=True)
    s_parser.add_argument('-s', '--signoff_step', dest='signoff_step_name', default='RELEASE', choices=['RELEASE'], help='Signoff step. Default is RELEASE')
    s_parser.add_argument('-dt', '--deliverable_type', dest='deliverable_type', default='Data Release', choices=['Data Release'], help='Deliverable type. Default is Data Release')
    s_parser.set_defaults(func=case_signoff)
    
    
    # write a report
    r_parser = subparsers.add_parser('report', help="Write a PDF report for released FASTQs")
    r_parser.add_argument('-pr', '--project', dest='project', help='Project name', required=True)
    r_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    r_parser.add_argument('-p', '--parents', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    r_parser.add_argument('-fn', '--full_name', dest='project_full_name', help='Full name of the project', required = True)
    r_parser.add_argument('-u', '--user', dest='user', help='Name of the GSI personnel generating the report', required = True)
    r_parser.add_argument('-t', '--ticket', dest='ticket', help='Jira data release ticket code', required = True)
    r_parser.add_argument('-pv', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json', help='Path to the json with production data. Default is /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json')
    r_parser.add_argument('-r', '--runs', dest='runs', nargs='*', help='List of run Ids')
    r_parser.add_argument('-c', '--cases', dest='cases', nargs='*', help='List of case Ids')
    r_parser.add_argument('-l', '--libraries', dest='libraries', help='File with libraries tagged for release. The first column is always the library. The second column is the run id and the third optional column is the lane number')
    r_parser.add_argument('-f', '--files', dest='release_files', help='File with file names or full paths of files to release')
    r_parser.add_argument('-a', '--analyses', dest='analyses', help='Path to the json file storing analysis data')
    r_parser.add_argument('-d', '--directories', dest='directories', nargs='*', help='List of directories with links or files to mark in Nabu')
    r_parser.add_argument('-bq', '--bamqc', dest='bamqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest', help='Path to the bamqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest')
    r_parser.add_argument('-dq', '--dnaseqqc', dest='dnaseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/dnaseqqc/latest', help='Path to the dnaseqqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/dnaseqqc/latest')
    r_parser.add_argument('-cq', '--cfmedipqc', dest='cfmedipqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest', help='Path to the cfmedip SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest')
    r_parser.add_argument('-rq', '--rnaseqqc', dest='rnaseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2/latest', help='Path to the rnaseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2/latest')
    r_parser.add_argument('-eq', '--emseqqc', dest='emseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest', help='Path to the emseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest')
    r_parser.add_argument('--keep_html', dest='keep_html', action='store_true', help='Write html report if activated.')
    r_parser.set_defaults(func=write_batch_report)
    
    
    # # write a report
    # c_parser = subparsers.add_parser('cumulative_report', help="Write a cumulative project PDF report for released FASTQs")
    # c_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report', required=True)
    # c_parser.add_argument('-p', '--parents', dest='projects_dir', default='/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/', help='Parent directory containing the project subdirectories with file links. Default is /.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS/')
    # c_parser.add_argument('-n', '--name', dest='project_name', help='Project name used to create the project directory in gsi space')
    # c_parser.add_argument('-fpr', '--provenance', dest='provenance', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    # c_parser.add_argument('-px', '--prefix', dest='prefix', help='Use of prefix assumes that FPR containes relative paths. Prefix is added to the relative paths in FPR to determine the full file paths')
    # c_parser.add_argument('-a', '--api', dest='api', default='https://nabu-prod.gsi.oicr.on.ca', help='URL of the Nabu API. Default is https://nabu-prod.gsi.oicr.on.ca')
    # c_parser.add_argument('-bq', '--bamqc', dest='merged_bamqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/bamqc4merged/latest', help='Path to the merged bamqc SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/bamqc4merged/latest')
    # c_parser.add_argument('-abq', '--archived_bamqc', dest='archived_bamqc', default = '/scratch2/groups/gsi/production/qcetl_archival/bamqc4merged/lastinput.json', help='Path to the archived merged bamqc json file. Default is /scratch2/groups/gsi/production/qcetl_archival/bamqc4merged/lastinput.json')
    # c_parser.add_argument('-rq', '--rnaseqqc', dest='merged_rnaseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2merged/latest', help='Path to the merged rnaseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2merged/latest')
    # c_parser.add_argument('-arq', '--archived_rnaseqqc', dest='archived_rnaseqqc', default = '/scratch2/groups/gsi/production/qcetl_archival/rnaseqqc2merged/lastinput.json', help='Path to the archived merged rnaseqqc json file. Default is /scratch2/groups/gsi/production/qcetl_archival/rnaseqqc2merged/lastinput.json') 
    # c_parser.add_argument('-cq', '--cfmedipqc', dest='cfmedipqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest', help='Path to the cfmedip SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest')
    # c_parser.add_argument('-acq', '--archived_cfmedipqc', dest='archived_cfmedipqc', default = '/scratch2/groups/gsi/production/qcetl_archival/cfmedipqc/lastinput.json', help='Path to the archived qc-etl json. Default is /scratch2/groups/gsi/production/qcetl_archival/cfmedipqc/lastinput.json')
    # c_parser.add_argument('-eq', '--emseqqc', dest='emseqqc_db', default = '/scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest', help='Path to the emseq SQLite database. Default is /scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest')
    # c_parser.add_argument('-aeq', '--archived_emseqqc', dest='archived_emseqqc', default = '/scratch2/groups/gsi/production/qcetl_archival/emseqqc/lastinput.json', help='Path to the archived qc-etl emseq json. Default is /scratch2/groups/gsi/production/qcetl_archival/emseqqc/lastinput.json')
    # c_parser.add_argument('-fn', '--full_name', dest='project_full_name', help='Full name of the project', required = True)
    # c_parser.add_argument('-t', '--ticket', dest='ticket', help='Jira data release ticket code', required = True)
    # c_parser.add_argument('-u', '--user', dest='user', help='Name of the GSI personnel generating the report', required = True)
    # c_parser.set_defaults(func=write_cumulative_report)

    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)




