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

Create a table file referencing internal IDs with external IDs privided by stakeholders.

usage: ```dare map -l LIBRARIES -w WORKFLOW -n PROJECT_NAME -p PROJECTS_DIR -pr PROJECT -r RUNS --exclude_miseq -rn RUN_NAME -e EXCLUDE -s SUFFIX -f FILES```

Parameters are the same as described above

Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project name indicated in FPR  | required              |
| -l | Path to file with libraries to be released   | optional                |
| -n | Project name | optional              |
| -p | Project directory   | required              |
| -r runs | Single or white space-separated run IDs  | optional              |
| --exclude_miseq | Excludes MiSeq runs | optional              |
| -e | Path to file with libraries to be excluded | optional              |
| -f | Path to file with file names/paths to be released | optional              |
| --time_points | Add column with time points | optional              |
| --panel | Add column with panel | optional              |
| -px | Prefix | optional              |
| -fpr | Path to the File Provenance Report | required              |
| -spr | Path to the Sample Provenance | required              |


Options as defined above. In addition to the `map` specific options:

- time points `--time_points`:
By default time points are not added. Add a column with time points retrieved from Pinery is used.

- time points `--time_points`:
By default time points are not added. Add a column with time points retrieved from Pinery is used.

- sample provenance `-spr/ -sample_provenance`:
Required parameter with default value: http://pinery.gsi.oicr.on.ca/provenance/v9/sample-provenance
Path to the sample provenance 


## Marking files in Nabu with dare.py ##

Tag released files in Nabu.

There are different ways to mark files in Nabu
For fastq files typically organized by run folder, one possibility is to point to the run directory containing the file links with the following command:

```dare mark -u USER -st pass -rn DIRECTORY -c TICKET -pr PROJECT```

For analysis pipelines, it is still possible to point to a directory containing links provided that links are not in sub-folders.
However, it is more convinient to mark files using the same json structure that was used to link the files (see baove section).

```dare mark -u USER -st pass -c TICKET -pr PROJECT -a JSON_FILE``` 



Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -u | Name of user handling the data release   | required              |
| -st | Mark files fail or pass | required              |
| -pr | Project name | required              |
| -rn | Directory with file links | optional              |
| -c | Jira ticket | optional              |
| -r | List of space-separated run Ids | optional              |
| -l | File with libraries tagged for release | optional              |
| -w | Workflow used to generate the output files | optional              |
| -px | Prefix to file paths if FPR contains relative paths | optional              |
| --exclude_miseq | Exclude miseq runs | optional              |
| -n | Nabu api | default              |
| -fpr | Path to FPR | default              |
| -e | Files with libraries tagged for non-release | optional              |
| -f | File with file names to be released | optional              |
| -a | Files with hierarchical structure storing sample and workflow ids | optional              |



## Generating a batch release report with dare.py ##

example usage: ```dare report -u USER -rn DIRECTORIES -t TICKET -pr PROJECT -fn PROJECT_FULL_NAME```

Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project acronym | required              |
| -fn | Project full name | required              |
| -u | Name of user handling the data release   | required              |
| -st | Mark files fail or pass | required              |
| -rn | List of space-sparated directories with file links | optional              |
| -t | Jira ticket | optional              |
| -r | List of space-separated run Ids | optional              |
| -l | File with libraries tagged for release | optional              |
| -w | Workflow used to generate the output files | optional              |
| -px | Prefix to file paths if FPR contains relative paths | optional              |
| --exclude_miseq | Exclude miseq runs | optional              |
| --time_points | Include time points | optional              |
| -a | Nabu api | default              |
| -spr | Sample provenance | default              |
| -fpr | Path to FPR | default              |
| -bq | BamQC cache | default              |
| -dq | DNASeqQC cache | default              |
| -cq | CfMedipQC cache | default              |
| -rq | RNASeqQC cache | default              |


## Transferring files through ociwire ##

Transfers files in `run_directory` through ociwire to UHN.

```bash transfer_out_ociwire.sh run_directory```

## Uploading files to the transfer server ##

Uploads files in directory `run_directory` to the transfer server.

```
tar -cvzhf NAMEORTAR.tar.gz run_directory
bash copy_to_transfer.sh NAMEOFTAR.tar.gz 
```

