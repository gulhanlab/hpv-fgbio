#!/usr/bin/env bash
# set -euo pipefail

module load gcc sambamba samtools

if [ $# -lt 1 ]; then
  echo "Usage: $0 <INPUT_DIR> [THREADS]"
  exit 1
fi

INPUT_DIR=$1
# Where sliced BAMs will go:
OUTPUT_DIR=$2
mkdir -p "$OUTPUT_DIR"

THREADS=${2:-4}

BED_FILE="./regions.bed"

echo "Using BED: $BED_FILE"
echo "Writing outputs to: $OUTPUT_DIR"
echo "Using $THREADS threads per job"

for bam in "${INPUT_DIR}"/*.bam; do
  sample=$(basename "$bam" .bam)
  out_bam="${OUTPUT_DIR}/${sample}_PIK3CA.bam"

  echo "[$(date)] Processing $sample …"

  samtools index $bam
  sambamba view \
    -t "$THREADS" \
    -f bam \
    -L "$BED_FILE" \
    "$bam" | \
  sambamba sort \
    -t "$THREADS" \
    -o "$out_bam" \
    /dev/stdin

  sambamba index \
    -t "$THREADS" \
    "$out_bam"

done

echo "All done. Sliced BAMs are in $OUTPUT_DIR"