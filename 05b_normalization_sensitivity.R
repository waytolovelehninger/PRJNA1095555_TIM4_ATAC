#!/usr/bin/env Rscript
# Input : $RESULTS/05_deseq2/common_peak_Male_dds.rds and timd4_normalized_counts.tsv (05_deseq2.R),
#         $RESULTS/04_consensus/common_peaks.bed, TSS_BED
# Output: $RESULTS/05b_normalization_sensitivity/normalization_sensitivity.tsv, size_factors.tsv, sessionInfo.txt
# Repeats the male DESeq2 test with two alternative size factors: (1) median-of-ratios on peaks within
# +-500 bp of a protein-coding TSS only, (2) total counts in peaks. Sensitivity check, not the main analysis.
.d <- { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)   # Rscript, RStudio Source, or console
  if (length(a)) dirname(sub("^--file=", "", a)) else if (!is.null(sys.frame(1)$ofile)) dirname(sys.frame(1)$ofile) else getwd() }
source(file.path(.d, "config.R"))
suppressPackageStartupMessages({ library(DESeq2); library(GenomicRanges) })
out <- file.path(RESULTS, "05b_normalization_sensitivity"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(file.path(RESULTS, "05_deseq2/common_peak_Male_dds.rds"))
timd4 <- unique(read.delim(file.path(RESULTS, "05_deseq2/timd4_normalized_counts.tsv"))$feature_id); stopifnot(length(timd4) == 1L)
bed <- read.delim(file.path(RESULTS, "04_consensus/common_peaks.bed"), header = FALSE, col.names = c("chr", "start0", "end0", "name"))
stopifnot(identical(bed$name, rownames(dds)))
tss <- read.delim(TSS_BED, header = FALSE, col.names = c("chr", "start0", "end0", "gene_id", "score", "strand"))

ctr <- c("age_group", "Old", "Young")
res0 <- results(dds, contrast = ctr, alpha = 0.05)

# (1) size factors from TSS-proximal peaks only (peak overlaps TSS +-500 bp)
pk <- GRanges(bed$chr, IRanges(bed$start0 + 1, bed$end0))
prox <- overlapsAny(pk, GRanges(tss$chr, IRanges(tss$end0 - 500, tss$end0 + 500)))
d1 <- DESeq(estimateSizeFactors(dds, controlGenes = prox), quiet = TRUE)
res1 <- results(d1, contrast = ctr, alpha = 0.05)

# (2) size factors from total counts in peaks (scaled to geometric mean 1)
d2 <- dds; tot <- colSums(counts(dds)); sizeFactors(d2) <- tot / exp(mean(log(tot)))
d2 <- nbinomWaldTest(estimateDispersions(d2, quiet = TRUE), quiet = TRUE)
res2 <- results(d2, contrast = ctr, alpha = 0.05)

row <- function(name, r, n_sf) data.frame(size_factors = name, peaks_used_for_size_factors = n_sf,
  Timd4_peak = timd4, Timd4_log2FC = r[timd4, "log2FoldChange"], Timd4_padj = r[timd4, "padj"],
  n_sig = sum(r$padj < 0.05, na.rm = TRUE), n_down = sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE),
  n_up = sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE))
summ <- rbind(row("median_of_ratios_all_peaks (main)", res0, nrow(dds)),
              row("median_of_ratios_TSS_proximal_peaks", res1, sum(prox)),
              row("total_counts_in_peaks", res2, nrow(dds)))
write.table(summ, file.path(out, "normalization_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = colnames(dds), main = sizeFactors(dds), TSS_proximal = sizeFactors(d1), total_counts = sizeFactors(d2)),
            file.path(out, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
