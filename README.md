# hpv-fgbio
## FGBIO Pipeline for processing HPV Datasets
### This pipeline is intended to run on 02 (HMS HPC).

## Directory Setup
1.	[02] Create directory for the current sample. I create folders in ```/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files```.
2.	[DropBox/02] Transfer fastq files from DropBox folder to current sample directory created in step 1. I use ```rclone``` to transfer files from DropBox to 02. Install ```rclone``` on your 02 account and sync with your DropBox account. Then run: ```rclone copy dropbox:"dropbox_folder" /n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers```.


## Directory Structure 
For each sample directory, subdirectory should look follwo this structure:
```
-sample
  -fastq_files
  -logs
  -scripts
  -UMI_extracted_bam
    -option1
    -option2
    -option3
    -option4
    -filtered_bams
      -option1
      -option2
      -option3
      -option4
      -merged_bams
        -count_matrices
        -qc_plots
        -simplex
        -simplex_filtered
        -spectra_plots
```
Here option1/2/3/4 refer to folders for different read orientation options. The UMI read orientation for the 4 options are as follows:
1.	Option1: 4M+T 6M+T
2.	Option2: 6M+T 6M+T
3.	Option3: 4M+T 4M+T
4.	Option4: 6M+T 4M+T

## Processing
Now we are ready to begin processing. Refer to ```fgbio_hpv.sh``` script for code. 
Optional if you are running on 02: You can setup a conda env and install necessary libs by executing Step 0.
For running QC and spectra plotting scripts, navigate to dir: ```/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/qc_scripts``` & ```/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/motif_spectra_scripts``` respectively.
