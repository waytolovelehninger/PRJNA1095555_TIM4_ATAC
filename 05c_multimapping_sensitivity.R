#!/usr/bin/env Rscript
# Input : $RESULTS/05c_multimapping/counts_{noXS,ASminusXS_gt5}.tsv (05c .py), $RESULTS/05_deseq2/ results (05)
# Output: $RESULTS/05c_multimapping/multimapping_sensitivity.tsv (Timd4 peak and genome-wide comparison with 05)
.d <- { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(sub("^--file=", "", a)) else if (!is.null(sys.frame(1)$ofile)) dirname(sys.frame(1)$ofile) else getwd() }
source(file.path(.d, "config.R"))
suppressPackageStartupMessages(library(DESeq2))
out <- file.path(RESULTS, "05c_multimapping")

md <- read.delim(SAMPLES)
md <- md[md$sex_verified == "Male", ]; rownames(md) <- md$run_accession
md$age_group <- factor(md$age_group, levels = c("Young", "Old"))
ref <- read.delim(file.path(RESULTS, "05_deseq2/common_peak_Male_full_results.tsv"), row.names = 1)
timd4 <- read.delim(file.path(RESULTS, "05_deseq2/timd4_normalized_counts.tsv"))$feature_id[1]
sig <- function(r) rownames(r)[!is.na(r$padj) & r$padj < 0.05]

summ <- function(name, res, nc) {
  s <- sig(res); s0 <- sig(ref); ok <- !is.na(res$log2FoldChange) & !is.na(ref$log2FoldChange)
  data.frame(counts = name, timd4_log2FC = res[timd4, "log2FoldChange"], timd4_padj = res[timd4, "padj"],
             timd4_norm_young = paste(round(range(nc[md$age_group == "Young"]), 1), collapse = "-"),
             timd4_norm_old = paste(round(range(nc[md$age_group == "Old"]), 1), collapse = "-"),
             fdr05 = length(s), fdr05_down = sum(res[s, "log2FoldChange"] < 0), fdr05_up = sum(res[s, "log2FoldChange"] > 0),
             retained_of_main = length(intersect(s, s0)), new_vs_main = length(setdiff(s, s0)),
             pearson_log2FC_vs_main = cor(res$log2FoldChange[ok], ref$log2FoldChange[ok]))
}
rows <- lapply(c("noXS", "ASminusXS_gt5"), function(k) {
  tab <- read.delim(file.path(out, paste0("counts_", k, ".tsv")), check.names = FALSE)
  mat <- t(as.matrix(tab[, -1])); colnames(mat) <- tab[[1]]; storage.mode(mat) <- "integer"
  dds <- DESeq(DESeqDataSetFromMatrix(mat[, md$run_accession], md, design = ~ age_group))
  res <- results(dds, contrast = c("age_group", "Old", "Young"), alpha = 0.05)[rownames(ref), ]
  summ(k, res, counts(dds, normalized = TRUE)[timd4, ])
})
main <- read.delim(file.path(RESULTS, "05_deseq2/timd4_normalized_counts.tsv"))
tab <- rbind(summ("main", ref, setNames(main$normalized_count, main$sample)[md$run_accession]), do.call(rbind, rows))
write.table(tab, file.path(out, "multimapping_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
