# Input : none (edit paths/tools below or override them as environment variables)
# Output: exported variables used by 01-04 and 07_fig1B_bigwig.sh
HERE_CFG=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export SAMPLES=${SAMPLES:-$HERE_CFG/data/samples.tsv}   # run_accession, ..., sex_verified, age_group
export WORK=${WORK:-$HERE_CFG/work}                      # per run: $WORK/<run>/{work,filtered,peaks}
export RESULTS=${RESULTS:-$HERE_CFG/results}
export NGMERGE=${NGMERGE:-NGmerge}                       # 0.3
export BOWTIE2=${BOWTIE2:-bowtie2}                       # 2.5.1
export SAMTOOLS=${SAMTOOLS:-samtools}                    # 1.19.2
export MACS2=${MACS2:-macs2}                             # 2.2.7.1
export DEEPTOOLS_VENV=${DEEPTOOLS_VENV:-$HOME/deeptools_3.5.5}
export BT2_INDEX=${BT2_INDEX:-$HERE_CFG/ref/mm10}        # Bowtie2 index prefix (mm10)
export BLACKLIST=${BLACKLIST:-$HERE_CFG/ref/mm10-blacklist.v2.bed}
export TSS_BED=${TSS_BED:-$HERE_CFG/data/TSS.gencode_M25.protein_coding.unique.bed}
export THREADS=${THREADS:-10}
RUNS=${RUNS:-$(tail -n +2 "$SAMPLES" | cut -f1)}
