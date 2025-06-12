#!/bin/bash
#SBATCH -J read_extraction
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -p short
#SBATCH -o /n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/Controls/logs/read_extraction.%J.out
#SBATCH -e /n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/Controls/logs/read_extraction.%J.err
#SBATCH --mem-per-cpu=4G

# Load required modules
module load gcc samtools

# Check that input_dir is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch $0 <input_dir>"
    exit 1
fi

input_dir="$1"

# Paths
# input_dir="/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/UMI_extracted_new_bam/filtered_bams/merged_bams/simplex"
output_dir_human="${input_dir}/human"
output_dir_hpv="${input_dir}/hpv"
output_dir_chimeric="${input_dir}/chimeric"

region_file_human="/n/data1/hms/dbmi/gulhan/lab/ankit/misc_files/HPV/human_refs.txt"
region_file_hpv="/n/data1/hms/dbmi/gulhan/lab/ankit/misc_files/HPV/hpv_refs.txt"

# Ensure output directories exist
mkdir -p "${output_dir_human}"
mkdir -p "${output_dir_hpv}"
mkdir -p "${output_dir_chimeric}"

# Loop through each .bam file in the input directory
for bam_file in "${input_dir}"/*.bam; do
    bam_basename=$(basename "${bam_file}" .bam)
    output_bam_human="${output_dir_human}/${bam_basename}_human.bam"
    output_bam_hpv="${output_dir_hpv}/${bam_basename}_hpv.bam"
    output_bam_chimeric="${output_dir_chimeric}/${bam_basename}_chimeric.bam"

    # Index BAM if needed
    if [ ! -f "${bam_file}.bai" ]; then
        echo "Indexing ${bam_file}..."
        samtools index "${bam_file}"
    fi

    # Extract human reads
    if [ ! -f "${output_bam_human}" ]; then
        echo "Extracting human reads from ${bam_file} -> ${output_bam_human}"
        samtools view -b "${bam_file}" $(cat "${region_file_human}") -o "${output_bam_human}"
    else
        echo "Human BAM already exists: ${output_bam_human}, skipping."
    fi

    # Extract HPV reads
    if [ ! -f "${output_bam_hpv}" ]; then
        echo "Extracting HPV reads from ${bam_file} -> ${output_bam_hpv}"
        samtools view -b "${bam_file}" $(cat "${region_file_hpv}") -o "${output_bam_hpv}"
    else
        echo "HPV BAM already exists: ${output_bam_hpv}, skipping."
    fi

    # Extract chimeric reads
    if [ ! -f "${output_bam_chimeric}" ]; then
        echo "Extracting chimeric reads from ${bam_file} -> ${output_bam_chimeric}"
        
        # Find reads mapped to human and HPV
        samtools view "${bam_file}" | \
            awk -v human_refs="$(paste -sd'|' "${region_file_human}")" \
                -v hpv_refs="$(paste -sd'|' "${region_file_hpv}")" \
                '
                BEGIN {
                    split(human_refs, human_arr, "|")
                    for (i in human_arr) human[human_arr[i]] = 1
                    split(hpv_refs, hpv_arr, "|")
                    for (i in hpv_arr) hpv[hpv_arr[i]] = 1
                }
                {
                    if ($3 in human || $3 in hpv) {
                        print $1 "\t" $3
                    }
                }
                ' | sort -k1,1 > "${output_dir_chimeric}/${bam_basename}_read_ref_map.txt"

        # Get reads appearing in both human and HPV
        awk '{print $1}' "${output_dir_chimeric}/${bam_basename}_read_ref_map.txt" | \
            uniq -d > "${output_dir_chimeric}/${bam_basename}_chimeric_read_ids.txt"

        # Extract chimeric reads
        samtools view -N "${output_dir_chimeric}/${bam_basename}_chimeric_read_ids.txt" -b "${bam_file}" -o "${output_bam_chimeric}"

        # Clean up intermediate files
        rm "${output_dir_chimeric}/${bam_basename}_read_ref_map.txt"
        rm "${output_dir_chimeric}/${bam_basename}_chimeric_read_ids.txt"
    else
        echo "Chimeric BAM already exists: ${output_bam_chimeric}, skipping."
    fi
done