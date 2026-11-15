#!/usr/bin/env Rscript
# Input : $RESULTS/04_consensus/sample_by_common_peak.raw_fragment_counts.tsv and common_peaks.bed, data/samples.tsv, GTF
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

# Timd4 TSS peak = the consensus peak containing the Timd4 TSS (GENCODE M25)
gtf <- rtracklayer::import(GTF, which = GRanges("chr11", IRanges(46700000, 46900000)))
tx <- gtf[gtf$type == "transcript" & gtf$gene_name == "Timd4"]; stopifnot(length(tx) == 1L)
tss <- if (as.character(strand(tx)) == "+") start(tx) else end(tx)
bed <- read.delim(file.path(RESULTS, "04_consensus/common_peaks.bed"), header = FALSE, col.names = c("chr", "start0", "end0", "name"))
timd4 <- bed$name[bed$chr == as.character(seqnames(tx)) & bed$start0 < tss & bed$end0 >= tss]
if (length(timd4) != 1L) stop("expected exactly one consensus peak containing the Timd4 TSS, found ", length(timd4))
write.table(data.frame(metric = "common_peak", feature_id = timd4, sample = md$run_accession, sex = "Male",
                       age_group = md$age_group, raw_count = counts(dds)[timd4, ],
                       normalized_count = counts(dds, normalized = TRUE)[timd4, ], size_factor = sizeFactors(dds)),
            file.path(out, "timd4_normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res[timd4, ])
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
