﻿# Data release #

Scripts in this repo are used for releasing data generated at OICR to stakeholders.
`dare.py` is used for data release of FASTQ and analysis files generated by pipeline and workflows and recorded in the File Provenance report (FPR). 

## Linking out files with dare.py ##

Links files in GSI space in directories under `/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS`
This option also create a list of md5sums for all the files that are being linked.

Use `-w` or `-a` to link fastqs or pipeline data respectively.


Example 1: raw sequence data
```python dare.py link -l LIBRARIES -w WORKFLOW -pr PROJECT -s SUFFIX```

Example 2: WGS pipeline data
```python dare.py link -pr PROJECT -a path_to_json_file```


Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| -pr | Project name indicated in FPR  | required              |
| -w | Workflow used to generate the files   | optional              |
| -n | Project name | optional              |
| -p | Project directory   | optional              |
| -s | Suffix appended to the run folder name | optional              |
| -l | Path to file with libraries to be released   | optional                |
| -e | Path to file with libraries to be excluded | optional              |
| -r runs | Single or white space-separated run IDs  | optional              |
| --exclude_miseq | Excludes MiSeq runs | optional              |
| --fpr | Path to the File Provenance Report | optional              |
| -f | Path to file with file names/paths to be released | optional              |
| -px | Prefix | optional              |
| -a | Path to json file storing information about analysis pipeline data | optional              |

- project `-pr/ --project`:
Required parameter. The name of the project as it appears in FPR.
If options `-r` or `-l` are not used then all the files from that project and for the specified `workflow` will be linked.

- workflow `-w/ --workflow`:
Required parameter. Workflow name.
Use `workflow` for non-merging workflows. data from merging workflows or analysis workflows is best handled with the `analysis` option (see below). 
The accepted value is bcl2fastq for when releasing FASTQs. This will pull out FASTQs generated by bcl2fastq, CASAVA and fastq-importing workflows from FPR.
Workflow name is case-incensitive.

- project name `-n/ --name`:
Optional parameter. The name of the project directory containing the run folders with symlinks in GSI space.
The project name specified with `-pr` is used to name the project directory by default. Use this option to change the name of the project folder.

- project directory `-p/ --parent`: 
Required parameter with default value `/.mounts/labs/gsiprojects/gsi/Data_Transfer/Release/PROJECTS`
Parent directory containing project directories with links to released files.

- suffix `-s/ --suffix`:
Optional parameter, but required when `-w` is used. 
The recommended value is `fastqs` for FASTQs release and a workflow name or data type for analysis files from non-merging workflows.
This parameter appends the suffix to the run folder name containing the file links. 

- libraries `-l/ --libraries` :
Optional parameter. Libraries approved for release.
Tab-delimited file with 2 columns with library and run Ids. 

- exclude `-e/ --exclude`:
Optional parameter. Libraries excluded from release.
Tab-delimited file with 2 columns with library and run Ids.

- runs `-r/ --runs`:
Optional parameter. Accepts one or more white space-separated run IDs.
Using the `-r` parameter will retrieve all the files corresponding to these runs for a given project and workflow.

- exclude MiSeq runs `--exclude_miseq`:
Use this flag if MiSeq runs are excluded from the release. 
If the flag is not activated then files sequenced on the MiSeq corresponding to a given project and library list will be included.

- fpr  `-fpr/ --fpr`:
Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz

- release files  `-f/ --files`:
File with file names of file paths of the files to be released

- prefix `-px/ --prefix`:
Use of prefix assumes that FPR containes relative paths.
Prefix is added to the relative paths in FPR to determine the full file paths.

- pipeline data `-a/ --analysis`:
Path to the file with hierarchical structure storing sample and workflow ids


### Linking pipeline data with the `-a` option ###

Option `-a` takes a json file with a hierarchical structure that describes the workflows and pipeline data
to be released. Different behaviors can be specified using keywords in the json structure.


#### Basic json structure ####

The basic structure is organized with sample ids, and workflow names and ids.

```
{donor_id:
   {sample_id:
      {workflow_1:
         [{"workflow_id":"19168526", "workflow_version":"2.0.2"}],
       workflow_2:
         [{"workflow_id":"16962244", "workflow_version":"2.0.2"}]
      }
   }
}   
```

This will create the following directory structure under the directory specified by the `-p` parameter:

```
donor_id
 |- sample_id
   |- workflow_1
     |- workflow_id
       |- file 1
         |- file 2
         |- ...
   |- workflow_2
     |- workflow_id
       |- file 1
       |- file 2
       |- ...
```

#### Renaming workflow sub-directories with the name keyword ####

Alternatively, the json can include an optional "name" key to replace the name of the sub-directories from workflows to names:

```
{donor_id:
   {sample_id:
     {workflow_1:
        [{"workflow_id": "19168526", "workflow_version":"2.0.2", "name": "calls.copynumber" }],
      workflow_2:
        [{"workflow_id":"16962244", "workflow_version":"2.0.2", "name": "calls.mutations"}]
     }
   }
}   
```

This json will create the following directory hierarchy

```
donor_id
  |- sample_id
    |- calls.copynumber
      |- workflow_id 
        |- file 1
        |- file 2
        |- ...
    |- calls.mutations
      |- workflow_id
        |- file 1
        |- file 2
        |- ...
``` 

#### Selecting output files with the extension keyword ####

By default, all the files generated by a workflow are collected and linked out according to the hierarchy specified by the json.Accepted
Files can be filtered and selected using optional keys "extension" or "files". Note that "extension" and "files" are mutually exclusive.

```
{donor_id:
  {sample_id:
    {workflow_1:
       [{"workflow_id": "19168526", "workflow_version":"2.0.2", "name": "calls.copynumber", "extension": [".gz"] }],
     workflow_2:
       [{"workflow_id":"16962244", "workflow_version":"2.0.2", "name": "calls.mutations", "extension": [".bai", ".bam"]}]
    }
  }  
}
```

This will only link `gz` files for worklow_1 and `bam` and `bai` files for workflow_2:

```
donor_id
  |- sample_id
    |- calls.copynumber
      |- workflow_id
        |- file 1.gz
        |- file 2.gz
        |- ....gz
    |- calls.mutations
      |- workflow_id
        |- file 1.bam
        |- file 2.bai
        | - ...
``` 


#### Selecting output files with the files keyword ####

Alternatively, outputfiles can be selected by specifying lists of files with the `files` keyword for each worklow.

```
{donor_id:
  {sample_id:
    {workflow_1:
      [{"workflow_id": "19168526", "workflow_version":"2.0.2", "name": "calls.copynumber", "files": [file1, file2, file3]}],
     workflow_2:
      [{"workflow_id":"16962244", "workflow_version":"2.0.2", "name": "calls.mutations", "files": [file1, file2]}]
    }
  }  
}
```
 
This will only link files specified with the `files` keyword.


```
donor_id:
  |- sample_id
    |- calls.copynumber
      |- workflow_id
        |- file 1
        |- file 2
        |- file 3
    |- calls.mutations
      |- workflow_id
        |- file 1
        |- file 2
``` 

#### Renaming output files with the rename_files keyword ####

Another way to select the workflow output files and to optionally rename them while doing so is to use the `rename_files` keyword.
Note that keywords `extension`, `files` and `rename_files` are mutually exclusive for a given workflow but can be mixed within the json.

The value of key `rename_files` is a list of 2-item dictionaries containing `file_path` and `file_name` key, value pairs.
This option allows changing the name of the link. By default, the name of the link is the file's basename.
Here the name of the link is specified by key `file_name`.

```
{donor_id:
  {sample_id:
    {workflow_1:
      [{"workflow_id": "19168526", "workflow_version":"2.0.2", "name": "calls.copynumber",
       "rename_files": [
         {"file_path": "file1", "file_name": "myfile1"}, {"file_path": "file2", "file_name": "this_is_a_file"}
        ]
      }],
    workflow_2:
      [{"workflow_id":"16962244", "workflow_version":"2.0.2", "name": "calls.mutations", "files": [file1, file2]}]
    }
  }  
}
```

This will result in selecting output files from workflows 1 and 2 and renaming the links of workflow 1.

```
donor_id:
  |- sample_id
    |- calls.copynumber
      |- workflow_id
        |- myfile1
        |- this_is_a_file
    |- calls.mutations
      |- workflow_id
        |- file 1
        |- file 2
``` 


## Mapping files to external IDs with dare.py ##

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

example usage: ```dare mark -u USER -st pass -rn DIRECTORY -c COMMENT -pr PROJECT```

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
| -a | Nabu api | default              |
| -fpr | Path to FPR | default              |


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

