#!/usr/bin/env bash
# Input : $WORK/<run>/work/fastq/<run>_{1,2}.fastq (SRA Toolkit 3.4.1: fasterq-dump --split-files)
# Output: $WORK/<run>/work/markdup.bam (coordinate-sorted, duplicates marked), $WORK/<run>/duplicate_metrics.txt
set -euo pipefail
source "$(dirname "$0")/config.sh"
for s in $RUNS; do
  S=$WORK/$s; W=$S/work; mkdir -p "$W/trimmed" "$W/tmp"
  "$NGMERGE" -a -q 33 -n 4 -1 "$W/fastq/${s}_1.fastq" -2 "$W/fastq/${s}_2.fastq" -o "$W/trimmed/$s"
  "$BOWTIE2" --phred33 --sensitive --no-discordant --no-mixed -k 2 -X 2500 -p 4 -x "$BT2_INDEX" \
    -1 "$W/trimmed/${s}_1.fastq" -2 "$W/trimmed/${s}_2.fastq" | "$SAMTOOLS" view -b -o "$W/aligned.unsorted.bam" -
  "$SAMTOOLS" sort -@ 2 -m 1G -T "$W/tmp/sort" -o "$W/coordinate_sorted.bam" "$W/aligned.unsorted.bam"
  "$SAMTOOLS" collate -@ 2 -T "$W/tmp/collate" -o "$W/name_collated.bam" "$W/coordinate_sorted.bam"
  "$SAMTOOLS" fixmate -@ 2 -m "$W/name_collated.bam" "$W/mate_fixed.bam"
  "$SAMTOOLS" sort -@ 2 -m 1G -O BAM -T "$W/tmp/fixsort" -o "$W/mate_fixed.coordinate_sorted.bam" "$W/mate_fixed.bam"
  "$SAMTOOLS" markdup -@ 2 -m t -S -T "$W/tmp/markdup" -f "$S/duplicate_metrics.txt" "$W/mate_fixed.coordinate_sorted.bam" "$W/markdup.bam"
done
