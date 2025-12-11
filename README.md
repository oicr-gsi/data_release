# Data release #

Scripts in this repo are used for releasing data generated at OICR to stakeholders.
`dare.py` is used for data release of FASTQ and analysis files generated, organizing data by project and case. 

## Linking files ##

Links files in GSI space in directories under `/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS`
This option also create a list of md5sums for all the files that are being linked.

Data is linked in the project folder and organized as:

case_id
 |- donor_id
   |- sequences
     |- run_folder
       |- fastq1
       |- fastq2
       .
       .
   |- analysis_workflow
     |- file1
     |- file2
     .
     .

Example 1: Link all the raw sequences from 2 runs for project JYCF 

```dare link -w bcl2fastq -pr JYCF -r 251104_M06816_0279_000000000-DVVRD 251029_M00753_0927_000000000-DVWBH```

Example 2: Link only the raw sequences for the libraries, runs and lanes specified in file LIBRARIES for project JYCF

```dare link -w bcl2fastq -l /path/to/LIBRARIES -pr JYCF```

Example 3: Link fastq and bmpp bams for project COMBAT for 2 runs

```dare link -w bcl2fastq bammergepreprocessing -pr COMBAT -r 20251119_LH00130_0289_A237M22LT3 20251029_LH00130_0285_B23CNVJLT4```

Example 4: Link data contained in the analysis.json for project COMBAT

```dare link -pr COMBAT -a /path/to/analysis.json```

Example 5: Link all the sequences for specific cases for project JYCF from run 251126_M06816_0286_000000000-GR3B7:

```dare link -pr JYCF -r 251126_M06816_0286_000000000-GR3B7 -w bcl2fastq -c R4547_a140_JYCF_0109_Es_P R4547_a140_JYCF_0111_Es_P R4547_a140_JYCF_0113_Es_P R4547_a140_JYCF_0122_Es_P R4547_a140_JYCF_0123_Es_P R4547_a140_JYCF_0124_Es_P R4547_a140_JYCF_0125_Es_P R4547_a140_JYCF_0126_Es_P R4547_a140_JYCF_0127_Es_P```


Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project name   | required              |
| -pv | Json file with production data | required             | 
| -p | Project directory | required             |
| -n | Project name | optional             |
| -w | List of worfklows | optional             |
| -r | List of runs | optional             |
| -l | File with libraries tagged for release | optional             |
| -f | File with names or paths of files to release | optional             |
| -c | List of case identifiers | optional             |
| -a | Path to json file storing data to release | optional             |


Option `-f` cannot be used with options `-w`, `-r`, `-c`, `-l` and `-a`.

Option `-a` cannot be used with options `-w`, `-r`, `-c`, `-l` and `-f`.

Option `-w` is required if `-a` and `-f` are not used.

Options `-r` and `-l` are mutually exclusive

- project `-pr/ --project`:
Required parameter. The name of the project. 

- provenance_json `-pv/ --provenance`: 
Path to the json with production data.
Required with default value: /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json

- project directory `-p/ --parent`:
Required parameter with default value `/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS`
Parent directory containing project directories with links to released files.

- project name `-n/ --name`:
Rename the project folder when creating the case folders and the symlinks.

- workflows `-w/ --workflows`:
White space separated list of workflows. Use bcl2fastq for any FASTQ-generating workflow (bcl2fastq, CASAVA and fastq-importing workflows).
Workflow name is case-incensitive.

- runs `-r/ --runs`:
White space-separated run identifiers.

- libraries `-l/ --libraries`:
Tab-delimited file with 2 or 3 columns, without header. The First column and second columns are always respectively the library and run identifiers. The third and optional column is the lane number. The lane number is either omitted or specified for all libraries. If the lane number is omitted, all the lanes for a given library and run will be linked.

- release files  `-f/ --files`:
File with names of paths of the files to be released

- cases `-c/ --cases`:
White space-separated case identifiers.

- analyses `-a/ --analyses`:
Path to the json file storing data to release (eg. from waterzooi)


## Generating a sample map ##

Create a table file with internal and external identifiers mapped to fastq files.

Example 1: Create a sample map for 2 runs for project JYCF

```dare map -pr JYCF -r 251104_M06816_0279_000000000-DVVRD 251029_M00753_0927_000000000-DVWBH```

Example 2: Create a sample map using only the libraries, runs and lanes specified in file LIBRARIES for project JYCF

```dare map -l /path/to/LIBRARIES -pr JYCF```

Example 3: Create a sample map from the fastqs listed in file release_files for project JYCF. Note that only fastqs are used to extract identifiers. Any other files are ignored.

```dare map -pr JYCF -f /path/to/release_files```

Example 4: Create a sample map from fastqs listed in file analysis.json for project COMBAT. Note that only fastqs are used to extract identifiers. Any other files are ignored.

```dare map -pr COMBAT -a /path/to/analysis.json```

Example 5: Create a sample map for specific cases for project JYCF and run 251126_M06816_0286_000000000-GR3B7:

```dare map -pr JYCF -r 251126_M06816_0286_000000000-GR3B7 -c R4547_a140_JYCF_0109_Es_P R4547_a140_JYCF_0111_Es_P R4547_a140_JYCF_0113_Es_P R4547_a140_JYCF_0122_Es_P R4547_a140_JYCF_0123_Es_P R4547_a140_JYCF_0124_Es_P R4547_a140_JYCF_0125_Es_P R4547_a140_JYCF_0126_Es_P R4547_a140_JYCF_0127_Es_P```

Example 6: Create a sample map from fastqs linked in directories for project JYCF:

```dare map -pr JYCF -d R4547_a140_JYCF_0109_Es_P R4547_a140_JYCF_0111_Es_P R4547_a140_JYCF_0113_Es_P```



Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project name   | required              |
| -pv | Json file with production data | required             |
| -p | Project directory | required             |
| -n | Project name | optional             |
| -r | List of runs | optional             |
| -l | File with libraries tagged for release | optional             |
| -f | File with names or paths of files to release | optional             |
| -c | List of case identifiers | optional             |
| -a | Path to json file storing data to release | optional             |
| -d | List of directories with links or files| optional             |


Option `-f` cannot be used with options `-d`, `-r`, `-c`, `-l` and `-a`.

Option `-a` cannot be used with options `-f`, `-r`, `-c`, `-l` and `-d`.

Option `-d` cannot be used with options `-f`, `-r`, `-c`, `-l` and `-a`.

Options `-r` and `-l` are mutually exclusive


- project `-pr/ --project`:
Required parameter. The name of the project.

- provenance_json `-pv/ --provenance`:
Path to the json with production data.
Required with default value: /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json

- project directory `-p/ --parent`:
Required parameter with default value `/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS`
Parent directory containing project directories with links to released files.

- project name `-n/ --name`:
Rename the project folder when creating the case folders and the symlinks.

- runs `-r/ --runs`:
White space-separated run identifiers.

- libraries `-l/ --libraries`:
Tab-delimited file with 2 or 3 columns, without header. The First column and second columns are always respectively the library and run identifiers. The third and optional column is the lane number. The lane number is either omitted or specified for all libraries. If the lane number is omitted, all the lanes for a given library and run will be linked.

- release files  `-f/ --files`:
File with names of paths of the files to be released

- cases `-c/ --cases`:
White space-separated case identifiers.

- analyses `-a/ --analyses`:
Path to the json file storing data to release (eg. from waterzooi)

- direcories `-d/ --directories``:
White space-separated directories with linked files. Traverses all the subdirectories and extracts identifiers only from linked fastqs.


## Marking files in Nabu ##

Changes FileQC status in Nabu following release.  

Example 1: Mark fastq files from 2 runs for project JYCF

```dare qc -pr JYCF -w bcl2fastq -r 251104_M06816_0279_000000000-DVVRD 251029_M00753_0927_000000000-DVWBH -u rjovelin -t GDR-0001 -st pass```

Example 2: Mark fastqs and bmpp output files for project COMBAT for 2 runs

```dare qc -w bcl2fastq bammergepreprocessing -pr COMBAT -r 20251119_LH00130_0289_A237M22LT3 20251029_LH00130_0285_B23CNVJLT4 -u rjovelin -t GDR-0001 -st pass```

Example 3: Mark files listed in file release_files for project JYCF

```dare qc -pr JYCF -f /path/to/release_files -u rjovelin -t GDR-0001 -st pass```

Example 4: Mark all files from the file analysis.json for project COMBAT

```dare qc -pr COMBAT -a /path/to/analysis.json -u rjovelin -t GDR-0001 -st pass```

Example 5: Mark all fastq files for specific cases for run 251126_M06816_0286_000000000-GR3B7 and project JYCF

```dare qc -pr JYCF -r 251126_M06816_0286_000000000-GR3B7 -c R4547_a140_JYCF_0109_Es_P R4547_a140_JYCF_0111_Es_P -w bcl2fastq -u rjovelin -t GDR-0001 -st pass```

Example 6: Mark all files linked in case directories for project JYCF:

```dare qc -pr JYCF -d R4547_a140_JYCF_0109_Es_P R4547_a140_JYCF_0111_Es_P R4547_a140_JYCF_0113_Es_P -u rjovelin -t GDR-0001 -st pass```


Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project name   | required              |
| -pv | Json file with production data | required             |
| -nb | URL of the Nabu API | required             |
| -st | QC status | required             |
| -u | User name | required             |
| -t | Jira ticket | required             |
| -w | List of worfklows | optional             |
| -r | List of runs | optional             |
| -l | File with libraries tagged for release | optional             |
| -f | File with names or paths of files to release | optional             |
| -c | List of case identifiers | optional             |
| -a | Path to json file storing data to release | optional             |
| -d | List of directories with links or files| optional             |


Option `-f` cannot be used with options `-w`, `-r`, `-c`, `-l`, `-a` and `-d`.

Option `-a` cannot be used with options `-w`, `-r`, `-c`, `-l`, `-f` and `-d`.

Option `-d` cannot be used with options `-w`, `-r`, `-c`, `-l`, `-f` and `-a`.

Option `-w` is required if `-a`, `-f` `-d` are not used.

Options `-r` and `-l` are mutually exclusive


- project `-pr/ --project`:
Required parameter. The name of the project.

- provenance_json `-pv/ --provenance`:
Path to the json with production data.
Required with default value: /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json

- nabu `-nb/ --nabu`:
URL of the Nabu API. Required with default value.

- user name `-u/ --user`:
Name of the GSI user releasing the date

- ticket `-t/ --ticket`:
Jira ticket

- runs `-r/ --runs`:
White space-separated run identifiers.

- libraries `-l/ --libraries`:
Tab-delimited file with 2 or 3 columns, without header. The First column and second columns are always respectively the library and run identifiers. The third and optional column is the lane number. The lane number is either omitted or specified for all libraries. If the lane number is omitted, all the lanes for a given library and run will be linked.

- release files  `-f/ --files`:
File with names of paths of the files to be released

- cases `-c/ --cases`:
White space-separated case identifiers.

- analyses `-a/ --analyses`:
Path to the json file storing data to release (eg. from waterzooi)

- direcories `-d/ --directories``:
White space-separated directories with linked files. Traverses all the subdirectories and extracts identifiers only from linked fastqs.


## Generating a batch release report ##

Generates a PDF report of a batch release of sequence data.

Example 1: Generate a batch report for project JYCF and a single run

```dare report -pr JYCF -r 251126_M06816_0286_000000000-GR3B7 -fn "Jonathan Yeung cfDNA project trial" -u "Richard Jovelin" -t GDR-0001```

Example 2: Generate a batch report for a single run and specific cases for project JYCF  

```dare report -pr JYCF -r 251126_M06816_0286_000000000-GR3B7 -c R4547_a140_JYCF_0109_Es_P R4547_a140_JYCF_0111_Es_P -fn "Jonathan Yeung cfDNA project trial" -u "Richard Jovelin" -t GDR-0001```


Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project name   | required              |
| -pv | Json file with production data | required             |
| -u | User name | required             |
| -t | Jira ticket | required             |
| -p | Project directory | required             |
| -fn | Full name of the project | required             |
| -n | Project name | optional             |
| -bq | Bamqc cache | required             |
| -dq | DNASeqQC cache | required             |
| -cq | CfMedipQC cache | required             |
| -rq | RNASeqQC cache | required             |
| -eq | EMSeqQC cache | required             |
| -r | List of runs | optional             |
| -l | File with libraries tagged for release | optional             |
| -f | File with names or paths of files to release | optional             |
| -c | List of case identifiers | optional             |
| -a | Path to json file storing data to release | optional             |
| -d | List of directories with links or files| optional             |
| --keep_html | Keep hmtl after PDF conversion | optional             |



Option `-f` cannot be used with options `-r`, `-c`, `-l`, `-a` and `-d`.

Option `-a` cannot be used with options `-r`, `-c`, `-l`, `-f` and `-d`.

Option `-d` cannot be used with options `-r`, `-c`, `-l`, `-f` and `-a`.

Options `-r` and `-l` are mutually exclusive


- project `-pr/ --project`:
Required parameter. The name of the project.

- provenance_json `-pv/ --provenance`:
Path to the json with production data.
Required with default value: /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json

- user name `-u/ --user`:
Name of the GSI user releasing the date

- ticket `-t/ --ticket`:
Jira ticket

- project directory `-p/ --parent`:
Required parameter with default value `/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS`
Parent directory containing project directories with links to released files.

- full name `-fn/ --full_name`:
Full name of the project

- project name `-n/ --name`:
Rename the project folder when creating the case folders and the symlinks.

- bamqc cache `-bq/ --bamqc`:
Path to the bamnqc cache. Required with default value /scratch2/groups/gsi/production/qcetl_v1/bamqc4/latest

- dnaseqqc cache `-dq/ --dnaseqqc`:
Path to the dnaseqqc cqache. Required with default /scratch2/groups/gsi/production/qcetl_v1/dnaseqqc/latest

- cfmedipqc cache `-cq/ --cfmedipqc`:
Path to the cfmedipqc cache. Required with default value /scratch2/groups/gsi/production/qcetl_v1/cfmedipqc/latest

- bamqc cache `-rq/ --rnaseqqc`:
Path to the rnaseqqc cache. Required with default value /scratch2/groups/gsi/production/qcetl_v1/rnaseqqc2/latest

- emseqqc cache `-eq/ --emseqqc`:
Path to the emseqqc cache. Required with default value /scratch2/groups/gsi/production/qcetl_v1/emseqqc/latest

- runs `-r/ --runs`:
White space-separated run identifiers.

- libraries `-l/ --libraries`:
Tab-delimited file with 2 or 3 columns, without header. The First column and second columns are always respectively the library and run identifiers. The third and optional column is the lane number. The lane number is either omitted or specified for all libraries. If the lane number is omitted, all the lanes for a given library and run will be linked.

- release files  `-f/ --files`:
File with names of paths of the files to be released

- cases `-c/ --cases`:
White space-separated case identifiers.

- analyses `-a/ --analyses`:
Path to the json file storing data to release (eg. from waterzooi)

- direcories `-d/ --directories``:
White space-separated directories with linked files. Traverses all the subdirectories and extracts identifiers only from linked fastqs.






## Case signoff in Nabu ##

Case signoff for FastQ and Full pipeline deliverables in Nabu. 

Example 1: Signoff of Fastq deliverable for a single run for project JYCF

```dare signoff -pr JYCF -r 251126_M06816_0286_000000000-GR3B7 -dv FastQ -u "Richard Jovelin" -t GDR-0001```

Example 2: Signoff of pipeline data for specific cases for project COMBAT

```dare signoff -pr COMBAT -c R1391_a94_COMBAT_0001_Pl_T_tT0 R1391_a94_COMBAT_0002_Pl_T_tT0 R1391_a94_COMBAT_0003_Pl_T_tT0 R1391_a94_COMBAT_0004_Pl_P_tT0 R1391_a94_COMBAT_0005_Pl_P_tT0 -dv "Full Pipeline" -u "Richard Jovelin" -t GDR-0001```


Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project name   | required              |
| -pv | Json file with production data | required             |
| -u | User name | required             |
| -t | Jira ticket | required             |
| -nb | URL of the Nabu API | required             |
| -nk | File with nabu key | required             |
| -dv | Deliverable | required             |
| -dt | Deliverable type | required             |
| -s | Signoff step | required             |
| -r | List of runs | optional             |
| -l | File with libraries tagged for release | optional             |
| -f | File with names or paths of files to release | optional             |
| -c | List of case identifiers | optional             |
| -a | Path to json file storing data to release | optional             |
| -d | List of directories with links or files| optional             |


Option `-f` cannot be used with options `-w`, `-r`, `-c`, `-l`, `-a` and `-d`.

Option `-a` cannot be used with options `-w`, `-r`, `-c`, `-l`, `-f` and `-d`.

Option `-d` cannot be used with options `-w`, `-r`, `-c`, `-l`, `-f` and `-a`.

Option `-w` is required if `-a`, `-f` `-d` are not used.

Options `-r` and `-l` are mutually exclusive


- project `-pr/ --project`:
Required parameter. The name of the project.

- provenance_json `-pv/ --provenance`:
Path to the json with production data.
Required with default value: /scratch2/groups/gsi/production/pr_refill_v2/provenance_reporter.json

- user name `-u/ --user`:
Name of the GSI user releasing the date

- ticket `-t/ --ticket`:
Jira ticket

- nabu `-nb/ --nabu`:
URL of the Nabu API

- nabu key `-nk/ --nabu_key`:
File with nabu key required to record signoffs. 
Required with default value

- deliverable `-dv'/ --deliverable`:
Deliverable. Value is FastQ or Full Pipeline.

- sign off step `-s/ --signoff_step`:
Name of the signoff step. Required value is RELEASE

- deliveranle type `-dt/ --deliverable_type`
Name of the deliverable type. Required value is Data Release.

- workflows `-w/ --workflows`:
White space separated list of workflows. Use bcl2fastq for any FASTQ-generating workflow (bcl2fastq, CASAVA and fastq-importing workflows).
Workflow name is case-incensitive.

- runs `-r/ --runs`:
White space-separated run identifiers.

- libraries `-l/ --libraries`:
Tab-delimited file with 2 or 3 columns, without header. The First column and second columns are always respectively the library and run identifiers. The third and optional column is the lane number. The lane number is either omitted or specified for all libraries. If the lane number is omitted, all the lanes for a given library and run will be linked.

- release files  `-f/ --files`:
File with names of paths of the files to be released

- cases `-c/ --cases`:
White space-separated case identifiers.

- analyses `-a/ --analyses`:
Path to the json file storing data to release (eg. from waterzooi)

- direcories `-d/ --directories``:
White space-separated directories with linked files. Traverses all the subdirectories and extracts identifiers only from linked fastqs.


