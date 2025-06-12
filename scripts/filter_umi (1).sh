#!/bin/bash
#SBATCH -J filter_umi
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -p medium
#SBATCH -o /n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/filter_umi.%J.out
#SBATCH -e /n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/filter_umi.%J.err
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=asingh46@mgh.harvard.edu

source /n/app/miniconda3/23.1.0/etc/profile.d/conda.sh 
conda activate /n/data1/hms/dbmi/gulhan/lab/ankit/conda_envs/SNVCurate

# Input BAM directory is provided as an argument when submitting the job
INPUT_BAM_DIR=$1
OUTPUT_BAM_DIR=$2

python /n/data1/hms/dbmi/gulhan/lab/ankit/scripts/fgbio/filter_umi.py "$INPUT_BAM_DIR" "$OUTPUT_BAM_DIR"