#!/bin/bash
#SBATCH -J merge_bams
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -p short
#SBATCH -o /n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/merge_bams.%J.out
#SBATCH -e /n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/merge_bams.%J.err
#SBATCH --mem-per-cpu=4G

module load gcc samtools

# Define folders
folder1=$1
folder2=$2
folder3=$3
folder4=$4
output_folder=$5
# csv="/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/Diagnostic/scripts/samples.csv"

# Create output directory if it doesn't exist
mkdir -p "$output_folder"

# Loop through BAM files in folder1 (assuming all folders have the same filenames)
# count=0
for bam_file in "$folder1"/*.bam; do
    # ((count++))
    # if [[ $count -gt 40 ]]; then
    #     break
    # fi
    filename=$(basename "$bam_file")  # Extract BAM filename
    echo "$filename"
    output_bam="$output_folder/$filename"
    
    # Skip if BAM is not listed in the CSV
    # if ! grep -q "$filename" "$csv"; then
    #     continue
    # fi

    # Check if BAM exists in all folders
    if [[ -f "$folder2/$filename" && -f "$folder3/$filename" && -f "$folder4/$filename" ]]; then
        echo "Merging $filename..."
        
        # Merge BAM files
        samtools merge -o "$output_bam" "$folder1/$filename" "$folder2/$filename" "$folder3/$filename" "$folder4/$filename"

        # Index the merged BAM
        # samtools index "$output_bam"
    else
        echo "Skipping $filename - Not found in all folders."
    fi
done

echo "BAM merging completed!"