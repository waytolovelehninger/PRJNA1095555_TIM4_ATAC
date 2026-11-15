#!/usr/bin/env bash
# Input : filtered male BAMs (02), $RESULTS/05_deseq2/timd4_normalized_counts.tsv (size factors), $TSS_BED
# Output: $RESULTS/07_fig1B/bigwig/*.sfnorm.bw and matrix/TSS_2kb_50bp.male.mat.gz (input of 07_fig1B_tss_heatmap.R)
# - deepTools 3.5.5 in its own venv (created on first run)
# - one fragment per pair (read1 only, extended to the mate: --samFlagInclude 64 --extendReads)
# - --maxFragmentLength 2500: without it deepTools truncates fragments > 4x median length (~720-816 bp)
#   to the median length; filtered BAMs contain no fragment > 2,500 bp, so this keeps every fragment's true span
# - scaled by 1 / DESeq2 size factor of the male common-peak fit (same normalization as the figures)
# - 10-bp bins, no Tn5 shift, no further normalization
set -euo pipefail
source "$(dirname "$0")/config.sh"
HERE=$RESULTS/07_fig1B
VENV=$DEEPTOOLS_VENV
mkdir -p "$HERE/bigwig" "$HERE/matrix" "$HERE/logs"
LOG="$HERE/logs/make_bigwigs.log"
exec > >(tee -a "$LOG") 2>&1
echo "== start $(date -Is)"

if [ ! -x "$VENV/bin/bamCoverage" ]; then
  PY=${PYTHON_FOR_VENV:-python3}                    # any Python 3 with venv/ensurepip
  rm -rf "$VENV"
  "$PY" -m venv "$VENV"
  "$VENV/bin/pip" install --quiet "deeptools==3.5.5"
fi
"$VENV/bin/bamCoverage" --version

bam_of() { echo "$WORK/$1/filtered/$1.primary.filtered.bam"; }

# sample, age group, size factor (male common-peak fit) from the male count table
SF="$RESULTS/05_deseq2/timd4_normalized_counts.tsv"
awk -F'\t' 'NR>1 && $1=="common_peak" && $2=="common_0010622" {print $3"\t"$5"\t"$8}' "$SF" | sort -k2,2r -k1,1 > "$HERE/samples.tsv"
[ "$(wc -l < "$HERE/samples.tsv")" -eq 9 ] || { echo "expected 9 male samples"; exit 1; }
cat "$HERE/samples.tsv"

while IFS=$'\t' read -r s grp sf; do
  bw="$HERE/bigwig/${s}.${grp}.sfnorm.bw"
  if [ -s "$bw" ]; then echo "exists: $bw"; continue; fi
  scale=$(awk -v x="$sf" 'BEGIN{printf "%.10f", 1/x}')
  echo "== bamCoverage $s ($grp) scale=$scale $(date -Is)"
  "$VENV/bin/bamCoverage" -b "$(bam_of "$s")" -o "$bw.tmp" --outFileFormat bigwig \
    --binSize 10 --extendReads --maxFragmentLength 2500 --samFlagInclude 64 --scaleFactor "$scale" --normalizeUsing None \
    -p "$THREADS"
  mv "$bw.tmp" "$bw"
done < "$HERE/samples.tsv"

# genome-wide TSS matrix: one TSS per protein-coding gene (GENCODE M25), +-2 kb, 50-bp bins
TSS="$TSS_BED"
BWS=$(awk -F'\t' -v d="$HERE/bigwig" '{printf "%s/%s.%s.sfnorm.bw ", d, $1, $2}' "$HERE/samples.tsv")
LABS=$(awk -F'\t' '{printf "%s_%s ", $1, $2}' "$HERE/samples.tsv")
if [ ! -s "$HERE/matrix/TSS_2kb_50bp.male.mat.gz" ]; then
  echo "== computeMatrix $(date -Is)"
  "$VENV/bin/computeMatrix" reference-point --referencePoint TSS -R "$TSS" -S $BWS --samplesLabel $LABS \
    -b 2000 -a 2000 --binSize 50 --missingDataAsZero -p "$THREADS" \
    -o "$HERE/matrix/TSS_2kb_50bp.male.mat.gz" --outFileSortedRegions "$HERE/matrix/TSS_2kb_50bp.regions.bed"
fi

# cross-check: dump 10-bp bigWig values over the Timd4 window for comparison with the R coverage
"$VENV/bin/python" - "$HERE" <<'PY'
import sys, csv, pyBigWig, glob, os
here = sys.argv[1]
start, end = 46800000, 46852000          # same window as the track figures (1-based start)
with open(os.path.join(here, "logs", "timd4_window_bigwig_values.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["sample", "pos_center", "bw_value"])
    for bw in sorted(glob.glob(os.path.join(here, "bigwig", "*.sfnorm.bw"))):
        s = os.path.basename(bw).split(".")[0]
        b = pyBigWig.open(bw)
        for x in range(start - 1, end - 9, 10):           # 0-based 10-bp bins
            v = b.stats("chr11", x, x + 10, type="mean")[0]
            w.writerow([s, x + 1 + 4.5, 0.0 if v is None else v])
        b.close()
print("timd4 window dumped")
PY
echo "== done $(date -Is)"
touch "$HERE/SUCCESS"
