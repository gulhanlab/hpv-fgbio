import pysam
import os
import sys
import Levenshtein

# Paths
umi_list_file = "/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/fgbio/umi_list.txt"
# bam_files_to_process = ['05_14_2024_DD-029_S29_L001.bam', '11_01_2024_21_S21.bam']
# bam_files_to_process = ['03_01_2024_9_S9_L001.bam', '03_01_2024_40_S40_L001.bam', '03_01_2024_39_S39_L001.bam', '03_01_2024_38_S38_L001.bam', '03_01_2024_37_S37_L001.bam', '03_01_2024_10_S10_L001.bam', '03_01_2024_27_S27_L001.bam', '03_29_2023_6_S6_L001.bam', '03_29_2023_5_S5_L001.bam', '03_29_2023_4_S4_L001.bam', '03_29_2023_18_S18_L001.bam', '03_29_2023_55_S55_L001.bam', '03_29_2023_19_S19_L001.bam', '03_29_2023_8_S8_L001.bam', '03_29_2023_9_S9_L001.bam', '03_29_2023_7_S7_L001.bam']

# Read input arguments from SLURM
if len(sys.argv) != 3:
    print("Usage: python3 filter_umi.py <input_bam_dir> <output_bam_dir>")
    sys.exit(1)

input_bam_dir = sys.argv[1]
output_bam_dir = sys.argv[2]

# Read valid UMIs into a set for fast lookup
with open(umi_list_file, "r") as f:
    valid_umis = set(line.strip() for line in f if line.strip())
    
# Function to find close matches allowing one letter mismatch
def find_near_matches(target1, target2, umi_list, max_distance=1):
    """Finds UMIs in umi_list that match target1 or target2 with at most one character mismatch."""
    matches1 = [umi for umi in umi_list if Levenshtein.distance(target1, umi) <= max_distance]
    matches2 = [umi for umi in umi_list if Levenshtein.distance(target2, umi) <= max_distance]
    return matches1, matches2

# Function to filter BAM based on UMIs
def filter_bam_by_umi(input_bam, output_bam):
    """Filters reads in BAM file based on valid UMI pairs."""
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as bam_in, \
         pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:

        for read in bam_in:
            # Extract UMI (RX tag)
            if read.has_tag("RX"):
                umi = read.get_tag("RX")  # Get the UMI sequence
                
                # Ensure the UMI contains a hyphen and split into parts
                if "-" in umi:
                    umi_part1, umi_part2 = umi.split("-")

                    # Find near matches for both parts
                    matches1, matches2 = find_near_matches(umi_part1, umi_part2, valid_umis)

                    # Ensure at least one close match for each part while keeping order intact
                    if matches1 and matches2:
                        bam_out.write(read)  # Write the read if both UMI parts have a valid near match

# Process each BAM file
for bam_file in os.listdir(input_bam_dir):
# for bam_file in bam_files_to_process:
    if bam_file.endswith(".bam") and not bam_file.startswith("filtered_"):
        input_bam_path = os.path.join(input_bam_dir, bam_file)
        output_bam_path = os.path.join(output_bam_dir, f"filtered_{bam_file}")
        
        # Skip if output file already exists
        if os.path.exists(output_bam_path):
            print(f"Skipping {bam_file} — output already exists.")
            continue

        print(f"Processing: {bam_file} → {output_bam_path}")
        filter_bam_by_umi(input_bam_path, output_bam_path)

print("UMI filtering completed!")