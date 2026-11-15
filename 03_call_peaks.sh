#!/usr/bin/env bash
# Input : $WORK/<run>/filtered/<run>.primary.filtered.bam (02)
# Output: $WORK/<run>/peaks/<run>_peaks.broadPeak (MACS2 2.2.7.1, broad, paired-end fragments)
set -euo pipefail
source "$(dirname "$0")/config.sh"
for s in $RUNS; do
  "$MACS2" callpeak -t "$WORK/$s/filtered/$s.primary.filtered.bam" -f BAMPE -g mm --keep-dup all --broad -n "$s" --outdir "$WORK/$s/peaks"
done
