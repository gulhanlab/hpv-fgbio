#!/bin/bash

# Check if a step argument is provided; if not, exit with usage message
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <step>"
  exit 1
fi

# Step to run
step=$1  

# Email for job notifications
email=

# Define paths for input files
# Path for the folder containing fastq files
fastq_path=/path/to/fastq
# Path for csv file for samples to proces in different steps
csv=/path/to/csv
# Only for step 2
out_bams=/path/to/unaligned_bams    # usually stored in ${fastq_path}/.FastqToSam as output from step1

# Define output directory t ostore UMI Extracted BAMs
out_dir_ext=/path/to/UMI_extracted_bam

# Alignment input files (Not needed until step 7 and 7b)
alignment_wdl=/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/alignment_scripts/processing-for-variant-discovery-gatk4.short.wdl
alignment_json=/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/alignment_scripts/processing-for-variant-discovery-gatk4.hg38.wgs.cgap.inputs_modified.json

# Paths to required software
fgbio_path=/n/data1/hms/dbmi/gulhan/lab/software/fgbio
picard_path=/n/data1/hms/dbmi/gulhan/lab/software/picard

# Reference genome paths
# ref=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/b37/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta    # b37
# ref=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/hg38/cgap_matches/Homo_sapiens_assembly38.fa                          # hg38
ref=/n/data1/hms/dbmi/gulhan/lab/REFERENCE/hg38/HPV_hg38/HPV_hg38.fasta                          # HPV

# Define output directories for each step

out1=${fastq_path}/.FastqToSam
out2=${out_dir_ext}/option1
out3=${out_dir_ext}/filtered_bams/option1
out4=${out3}/.PreProcessing
out5=${out3}/step5_out
out6=${out3}/step6_out
out7=${out3}/step7_out
out8=${out3}/step8_out
out9=${out3}/step9_out
out10=${out9}/.PreProcessing
out11=${out3}/step11_out
out12=${out3}/step12_out
out13=${out3}/step13_out
out14=${out3}/merged_bams/simplex
out13_option1=${out3_option1}/step13_out    #Step 13 output for option1
out13_option2=${out3_option2}/step13_out    #Step 13 output for option2
out13_option3=${out3_option3}/step13_out    #Step 13 output for option3
out13_option4=${out3_option4}/step13_out    #Step 13 output for option4


# Define actions for each step
case $step in
  0)
    # ** Step 0: Set up Conda environment **
    echo "Set up conda env (only need to do it once). If already set up, use it."
    echo "conda env create -f /n/data1/hms/dbmi/gulhan/lab/ankit/scripts/mutation_calling/fgbio/general3.yml"
    echo "To activate the environment:"
    echo "conda activate general3; module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa"
  ;;
  
  1)
    ## ** Step 1: converting fastq.gz files -> unaligned bam files. Creates .FastqToSam folder in the output **
    echo "converting fastq.gz files -> unaligned bam files. Creates .FastqToSam folder in the output"
    echo "change the -stall param as needed and the -r1 and -r2 as needed"
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa 
    source /n/app/miniconda3/23.1.0/etc/profile.d/conda.sh 
    conda activate /n/data1/hms/dbmi/gulhan/lab/ankit/conda_envs/SNVCurate
    python3 /n/data1/hms/dbmi/gulhan/lab/ankit/scripts/alignment_scripts/FastqToSam2.py \
        -in_dir ${fastq_path} \
        -out ${fastq_path} \
        -csv ${csv} \
        --mem_per_cpu 10G -t 12:00:00 -stall 1 -r1 0 -r2 220
    ;;
    

  2)
    # ** Step 2: Extract UMI from unmapped bams **
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out2"

    for bam in "$out1"/*.bam; do
        name=$(basename "$bam" .bam)
        out_file="$out2/${name}.bam"

        # # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        sbatch -J "step11_${name}" -p short -t 60:00 --mem=5000 \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step2_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step2_%j.err" \
            --wrap="echo $out_file; \
            java -jar ${fgbio_path}/fgbio-2.1.0.jar ExtractUmisFromBam \
            -i ${bam} \
            -o ${out_file} \
            -r 4M+T 6M+T \
            -t RX \
            -a true"
    done
  ;;
  

  3)
    # ** Step 3: Filter UMIs **
    mkdir -p "$out3"
    sbatch filter_umi.sh ${out2} ${out3}
  ;;


  4)
    # ** Step 4: align to HPV_hg38 **
    echo "unaligned bam files -> aligned bam files. Creates a .PreProcessing folder"
    echo "change the -t param if running smaller files, and can change the processing "
    echo "change the -stall param as needed and the -r1 and -r2 as needed"
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa
    source /n/app/miniconda3/23.1.0/etc/profile.d/conda.sh
    conda activate /home/ans4371/.conda/envs/general3
    python3 /n/data1/hms/dbmi/gulhan/lab/ankit/scripts/alignment_scripts/PreProcessing.py \
        -in_dir ${out3} \
        -out ${out3} \
        -input_json ${alignment_json} \
        -gatk_wdl ${alignment_wdl} \
        -n 2 -t 12:00:00 \
        -p short \
        --mem_per_cpu 3G \
        -stall 100 -r1 0 -r2 80
  ;;
  

  5)
    # ** Step 5: Fix Mate Information (Picard: FixMateInformation) **
    module load gcc java/jdk-21.0.2 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out5" 

    for i in ${out4}/*.bam; do
        bam_name=$(basename "$i")
        out_file="$out5/$bam_name"

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        sbatch -J "step5_${bam_name}" -p short -t 01:00:00 --mem=4G \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step5_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step5_%j.err" \
            --wrap="echo $out_file; \
            java -jar ${picard_path}/picard.jar FixMateInformation \
            I=$i \
            O=$out_file \
            SORT_ORDER=queryname \
            TMP_DIR=${out_dir}/tmp_fgbio"
    done
  ;;


  6) 
    # ** Step 6: Group Reads by UMI (fgbio: GroupReadsByUmi) **
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out6"  

    for i in "$out5"/*.bam; do
        bam_name=$(basename "$i")
        echo "$bam_name"

        # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        out_file="$out6/$bam_name"

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        sbatch -J "step6_${bam_name}" -p short -t 45:00 --mem=3G \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step6_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step6_%j.err" \
            --wrap="echo $out_file; \
            export _JAVA_OPTIONS=-Djava.io.tmpdir=${out_dir}/tmp_fgbio; \
            java -jar ${fgbio_path}/fgbio-2.1.0.jar GroupReadsByUmi -s adjacency \
            -i $i \
            -o $out_file \
            -f ${out_file}.family_histogram.txt \
            -e 1"
    done
  ;;
  

  7)
    # **Step 7: Generate Molecular Consensus Reads (fgbio: CallMolecularConsensusReads) **
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out7" 

    for i in "$out6"/*.bam; do
        bam_name=$(basename "$i")
        out_file="$out7/$bam_name"

        # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        sbatch -J "step7_${bam_name}" -p short -t 45:00 --mem=3000 \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/fgbio_step7_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/fgbio_step7_%j.err" \
            --wrap="echo $out_file; \
            java -jar ${fgbio_path}/fgbio-2.1.0.jar CallMolecularConsensusReads \
            -i $i \
            -o $out_file \
            -r ${out_file}.reject.bam \
            -M 1 \
            --max-reads 50 \
            --min-input-base-quality 20"
    done
  ;;
  

  8)
    # ** Step 8: Filter Consensus Reads (fgbio: FilterConsensusReads) **
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out8"

    for i in "$out7"/*.bam; do
        [[ "$i" == *.reject.bam ]] && continue
        bam_name=$(basename "$i")
        out_file="$out8/$bam_name"

        # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

         sbatch -J "step8_${bam_name}" -p short -t 45:00 --mem=6000 \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step8_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step8_%j.err" \
            --wrap="echo $out_file; \
            java -jar ${fgbio_path}/fgbio-2.1.0.jar FilterConsensusReads \
            -i $i \
            -o $out_file \
            -r $ref \
            -M 1 -N 3"
    done
  ;;

  9)
    # ** Step 9: Sort BAMs for simplex consensus reads (Picard: SortSam) **
    module load gcc java/jdk-21.0.2 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out9" 
    
    for i in "$out8"/*.bam; do
        bam_name=$(basename "$i")
        out_file="$out9/$bam_name"

        # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        sbatch -J "step9_${bam_name}" -p short -t 30:00 --mem=6000 \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step9_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step9_%j.err" \
            --wrap="echo $out_file; \
            java -Xmx12G -jar ${picard_path}/picard.jar SortSam \
            I=$i \
            O=$out_file \
            SORT_ORDER=coordinate \
            TMP_DIR=${out_dir}/tmp_fgbio"
    done
  ;;

  10)
    # ** Step 10: Align BAM files to reference genome (HPV_hg38) post simplex-collapse **
    source /n/app/miniconda3/23.1.0/etc/profile.d/conda.sh
    conda activate /home/ans4371/.conda/envs/general3
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa
    python3 /n/data1/hms/dbmi/gulhan/lab/ankit/scripts/alignment_scripts/PreProcessing.py \
        -in_dir ${out9} \
        -out ${out9} \
        -gatk_wdl ${alignment_wdl} \
        -input_json ${alignment_json} \
        -n 1 -t 12:00:00 -p short --mem_per_cpu 4G -stall 120 -r1 0 -r2 80
  ;;

  11) 
    # ** Step 11: Sort BAMs by query name (Picard: SortSam) **
    module load gcc java/jdk-21.0.2 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p $out11

    for i in ${out10}/*bam; do
        bam_name=$(basename $i)
        out_file=$out11/$bam_name

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        sbatch -J "step11_${bam_name}" -p short -t 30:00 --mem=4000 --mail-user=$email --mail-type=FAIL -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step11_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step11_%j.err" --wrap="echo $out_file; \
          java -Xmx12G -jar ${picard_path}/picard.jar SortSam \
          I=$i \
          O=$out_file \
          SORT_ORDER=queryname \
          TMP_DIR=${out_dir}/tmp_fgbio"
    done
  ;;

  12)
    # ** Step 12: Filter Consensus Reads (fgbio: FilterConsensusReads) **
    module load gcc java/jdk-1.8u112 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out12" 

    for i in "$out11"/*.bam; do
        bam_name=$(basename "$i")
        out_file="$out12/$bam_name"

        # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        echo "Running ${bam_name}..."

        sbatch -J "step12_${bam_name}" -p short -t 20:00 --mem=6500 \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step12_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step12_%j.err" \
            --wrap="echo $out_file; \
            java -jar ${fgbio_path}/fgbio-2.1.0.jar FilterConsensusReads \
            -i $i \
            -o $out_file \
            -r $ref \
            -M 1 -N 3"
    done
  
  ;;

  13)
    # ** Step 13: Final simplex BAM sorting by coordinate (Picard: SortSam) **
    module load gcc java/jdk-21.0.2 python/3.6.0 samtools htslib bcftools bwa
    mkdir -p "$out13"

    for i in "$out12"/*.bam; do
        bam_name=$(basename "$i")
        out_file="$out13/$bam_name"

        # Skip if BAM is not listed in the CSV
        # if ! grep -q "$bam_name" "$csv"; then
        #     continue
        # fi

        # Skip if output file already exists
        if [ -e "$out_file" ]; then
            echo "$out_file exists, skipping."
            continue
        fi

        sbatch -J "step13_${bam_name}" -p short -t 20:00 --mem=5000 \
            --mail-user="$email" --mail-type=FAIL \
            -o "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step13_%j.out" -e "/n/data1/hms/dbmi/gulhan/lab/DATA/HPV-ctDNA/FASTQ_files/PreCancers/logs/fgbio_step13_%j.err" \
            --wrap="echo $out_file; \
            java -Xmx12G -jar ${picard_path}/picard.jar SortSam \
            I=$i \
            O=$out_file \
            SORT_ORDER=coordinate \
            TMP_DIR=${out_dir}/tmp_fgbio"
    done
  ;;

  14)
    # ** Step 14: Merge BAMs across the 4 options **
    mkdir -p "$out14"
    sbatch merge_bams.sh ${out13_option1} ${out13_option2} ${out13_option3} ${out13_option4} ${out14}
  ;;

  15)
    # ** Step 15: Extract human, HPV & chimeric reads **
    sbatch extract_reads.sh ${out14}
  ;;

  16)
    # ** Step 16: Filter BAMs **
    # ** Modify config file appropriately for qc script. Demonstarted is an example for QC filtering chimeric read BAM **
    # ** NOTE: Remove the length filter for chimeric reads. Modify script filter_bam_v4_hpv.py to comment out ```abs(read.template_length) <= max_length``` (line 59) **
    ./run_qc_and_filtering.sh ./config_Controls_chimeric.sh 3 false
  ;;

  17)
    # ** Step 17: Intersect human reads with PIK3CA and GAPDH **
    ./extract_PIK3CA.sh ${out14}/simplex_filtered/human ${out14}/simplex_filtered/human/intersected_bams
  ;;

  18)
    # ** Step 18: Compute count matrices: Step 0 **
    # ** Modify config file appropriately for motif spectra script. Demonstarted is an example for step 0 for chimeric read BAM. This step creates bed file files per sample with counts metadata **
    # ** Stores output in 'counts_dir' defined in config file **
    ./run_extract_shift_count_collect.sh ./config_Followup_chimeric.sh 0
  ;;

  19)
    # ** Step 19: Compute count matrices: Step 1 **
    # ** Modify config file appropriately for motif spectra script. Demonstarted is an example for step 1 for chimeric read BAM. This step creates count matrices file for out5p, in5p and length motifs each. **
    # ** Stores output in .... **
    ./run_extract_shift_count_collect.sh ./config_Followup_chimeric.sh 1
  ;;

  20)
    # ** Step 20: Directory setup for motif spectra plotting **
    # ** 'collected_counts_dir' parameter should be adjusted in the config file motif spectra script **
    cp ${collected_counts_dir}/final_counts/*length.csv ${collected_counts_dir}/final_counts_length/
    cp ${collected_counts_dir}/final_counts/*in5p.csv ${collected_counts_dir}/final_counts_motif/
    cp ${collected_counts_dir}/final_counts/*out5p.csv ${collected_counts_dir}/final_counts_motif/
    mv ${collected_counts_dir}/final_counts/*.csv ${out14}/count_matrices/chimeric
  ;;

  21)
    # ** Step 21: Plotting out5p & in5p motif spectra plots **
    ./run_extract_shift_count_collect.sh ./config_Followup_chimeric.sh 2
  ;;

  22)
    # ** Step 22: Plotting length motif spectra plots **
    ./run_extract_shift_count_collect.sh ./config_Followup_chimeric.sh 3
  ;;

  23)
    # ** Step 23: QC Plots step 0 **
    # ** Modify config file appropriately for qc script. Demonstarted is an example for step 0 for chimeric read BAM. This step extracts information for QC plot. **
    # ** Stores output in 'extracted_features_dir' defined in the config file **
    ./run_qc_and_filtering.sh ./config_Followup_chimeric.sh 0 false
  ;;

  24)
    # ** Step 23: QC Plots step 1 **
    # ** Modify config file appropriately for qc script. Demonstarted is an example for step 0 for chimeric read BAM. This step creates the QC plots **
    ./run_qc_and_filtering.sh ./config_Followup_chimeric.sh 1 false
  ;;
  
  *)
    # ** Invalid step handling **
    echo "Invalid step: $step"
    exit 1
  ;;

esac
