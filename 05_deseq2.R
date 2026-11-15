#!/usr/bin/env Rscript
# Input : $RESULTS/04_consensus/sample_by_common_peak.raw_fragment_counts.tsv, data/samples.tsv
# Output: $RESULTS/05_deseq2/common_peak_Male_full_results.tsv, timd4_normalized_counts.tsv, common_peak_Male_dds.rds
.d <- { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)   # Rscript, RStudio Source, or console
  if (length(a)) dirname(sub("^--file=", "", a)) else if (!is.null(sys.frame(1)$ofile)) dirname(sys.frame(1)$ofile) else getwd() }
source(file.path(.d, "config.R"))
suppressPackageStartupMessages(library(DESeq2))
out <- file.path(RESULTS, "05_deseq2"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

md <- read.delim(SAMPLES)
md <- md[md$sex_verified == "Male", ]; rownames(md) <- md$run_accession
md$age_group <- factor(md$age_group, levels = c("Young", "Old"))
tab <- read.delim(file.path(RESULTS, "04_consensus/sample_by_common_peak.raw_fragment_counts.tsv"), check.names = FALSE)
mat <- t(as.matrix(tab[, -1])); colnames(mat) <- tab[[1]]; storage.mode(mat) <- "integer"

dds <- DESeq(DESeqDataSetFromMatrix(mat[, md$run_accession], md, design = ~ age_group))
res <- results(dds, contrast = c("age_group", "Old", "Young"), alpha = 0.05)
write.table(as.data.frame(res), file.path(out, "common_peak_Male_full_results.tsv"), sep = "\t", quote = FALSE, col.names = NA)
saveRDS(dds, file.path(out, "common_peak_Male_dds.rds"))

timd4 <- "common_0010622"   # Timd4 TSS peak, chr11:46,810,221-46,812,849
write.table(data.frame(metric = "common_peak", feature_id = timd4, sample = md$run_accession, sex = "Male",
                       age_group = md$age_group, raw_count = counts(dds)[timd4, ],
                       normalized_count = counts(dds, normalized = TRUE)[timd4, ], size_factor = sizeFactors(dds)),
            file.path(out, "timd4_normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res[timd4, ])
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
