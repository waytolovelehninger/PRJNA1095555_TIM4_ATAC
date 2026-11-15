#!/usr/bin/env Rscript
# Input : filtered BAMs (02), $RESULTS/04_consensus/common_peaks.bed, $RESULTS/05_deseq2/*, GENCODE M25 GTF, Lu et al. Table S2
# Output: $RESULTS/06_fig1A/figures/timd4_locus_track_*.png|pdf and results/timd4_coverage_10bp_*.tsv (Fig. 1A)
options(stringsAsFactors = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "config.R"))
out <- file.path(RESULTS, "06_fig1A")
suppressPackageStartupMessages({
  library(GenomicAlignments); library(rtracklayer); library(ggplot2); library(patchwork)
})
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), showWarnings = FALSE)

# ---- samples, BAMs, size factors -------------------------------------------------------
meta <- read.delim(SAMPLES); meta <- meta[meta$sex_verified == "Male", ]
bam <- setNames(bam_path(meta$run_accession), meta$run_accession)
stopifnot(all(file.exists(bam)), all(file.exists(paste0(bam, ".bai"))))

nc <- read.delim(file.path(RESULTS, "05_deseq2/timd4_normalized_counts.tsv"))
nc <- nc[nc$metric == "common_peak" & nc$feature_id == "common_0010622", ]
md <- data.frame(sample = nc$sample, age_group = nc$age_group, size_factor = nc$size_factor,
                 raw_peak_count = nc$raw_count)
md$title <- sub("_PeriMac", "", sub("_ATACseq$", "", meta$sample_title[match(md$sample, meta$run_accession)]))
md <- md[order(factor(md$age_group, c("Young", "Old")), md$title), ]
stopifnot(nrow(md) == 9L, setequal(md$sample, names(bam)),
          identical(as.vector(table(factor(md$age_group, c("Young", "Old")))), c(5L, 4L)),
          all(meta$sex_verified == "Male"))

# Cross-check: size factors equal those stored in the male DESeq2 object
dds_path <- file.path(RESULTS, "05_deseq2/common_peak_Male_dds.rds")
if (file.exists(dds_path)) {
  sf <- DESeq2::sizeFactors(readRDS(dds_path))
  stopifnot(max(abs(sf[md$sample] - md$size_factor)) < 1e-9)
}

# ---- Timd4 annotation (GENCODE M25) and peaks ------------------------------------------
gtf <- import(GTF,
              which = GRanges("chr11", IRanges(46700000, 46900000)))
tx <- gtf[gtf$type == "transcript" & gtf$gene_name == "Timd4"]
ex <- gtf[gtf$type == "exon" & gtf$gene_name == "Timd4"]
cds <- gtf[gtf$type == "CDS" & gtf$gene_name == "Timd4"]
stopifnot(length(tx) == 1L, as.character(strand(tx)) == "+")
TSS <- start(tx)  # 46810800

peaks <- import(file.path(RESULTS, "04_consensus/common_peaks.bed"),
                which = GRanges("chr11", IRanges(46790000, 46860000)))
res <- read.delim(file.path(RESULTS, "05_deseq2/common_peak_Male_full_results.tsv"))
names(res)[1] <- "feature"
pk <- data.frame(name = peaks$name, start = start(peaks), end = end(peaks))  # 1-based closed
pk <- merge(pk, res[, c("feature", "log2FoldChange", "padj")], by.x = "name", by.y = "feature", all.x = TRUE)
# Published ATAC peaks (Lu et al., Table S2 Sheet B) in the window; Excel Start/End are 1-based closed
sb <- as.data.frame(readxl::read_excel(PUBLISHED_XLSX, sheet = "B"))
pub <- sb[sb$Chr == "chr11" & sb$End >= 46790000 & sb$Start <= 46860000, c("PeakName", "Start", "End", "padj_M")]
names(pub) <- c("name", "start", "end", "padj_M")
pub$tss <- pub$name == "chr11:46810157:46812849"
stopifnot(sum(pub$tss) == 1L)

# ---- fragment coverage ------------------------------------------------------------------
REG <- GRanges("chr11", IRanges(46800000, 46852000))      # region read from the BAMs
query <- resize(REG, width(REG) + 5000, fix = "center")   # max retained fragment <= 2,500 bp
BIN <- 10L
bins <- seq(start(REG), end(REG) - BIN + 1L, by = BIN)
peak_tgt <- GRanges("chr11", IRanges(46810221, 46812849))  # common_0010622, 1-based closed

cov_tab <- list(); check <- list()
for (s in md$sample) {
  pr <- readGAlignmentPairs(bam[[s]], param = ScanBamParam(which = query))
  fr <- granges(pr, on.discordant.seqnames = "drop")       # outer fragment span, unshifted
  cv <- coverage(fr)$chr11
  v <- as.numeric(cv[start(REG):end(REG)])
  binned <- colMeans(matrix(v[seq_len(length(bins) * BIN)], nrow = BIN))
  sfac <- md$size_factor[md$sample == s]
  cov_tab[[s]] <- data.frame(sample = s, pos = bins + BIN / 2, raw = binned, norm = binned / sfac)
  check[[s]] <- data.frame(sample = s, fragments_overlapping_common_0010622 = sum(overlapsAny(fr, peak_tgt)),
                           stored_raw_count = md$raw_peak_count[md$sample == s])
}
check <- do.call(rbind, check)
check$match <- check$fragments_overlapping_common_0010622 == check$stored_raw_count
write.table(check, file.path(out, "results/peak_count_crosscheck.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(check)
if (!all(check$match)) warning("Fragment counts at common_0010622 differ from the stored count matrix; see results/peak_count_crosscheck.tsv")

cv <- do.call(rbind, cov_tab)
cv$age_group <- md$age_group[match(cv$sample, md$sample)]
write.table(cv, file.path(out, "results/timd4_coverage_10bp_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
grp <- aggregate(norm ~ pos + age_group, cv, function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
grp <- data.frame(pos = grp$pos, age_group = grp$age_group, mean = grp$norm[, "mean"], se = grp$norm[, "se"])
write.table(grp, file.path(out, "results/timd4_coverage_10bp_group_mean.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- plotting ---------------------------------------------------------------------------
COL <- c(Young = "#3B6FB6", Old = "#D1495B")
male_stat <- res[res$feature == "common_0010622", ]
stat_lab <- sprintf("common_0010622  log2FC %s,  padj %s", sub("^-", "−", sprintf("%.2f", male_stat$log2FoldChange)),
                    formatC(male_stat$padj, format = "g", digits = 2))
base_theme <- theme_classic(base_size = 12) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), plot.margin = margin(2, 8, 2, 8),
        strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0))

make_track <- function(xlim, file_stem, title, step_kb, show_pub = TRUE) {
  g <- grp[grp$pos >= xlim[1] & grp$pos <= xlim[2], ]
  ind <- cv[cv$pos >= xlim[1] & cv$pos <= xlim[2], ]
  ymax <- max(g$mean + g$se) * 1.08
  hl <- data.frame(xmin = 46810221, xmax = 46812849)
  lab_n <- c(Young = sprintf("Young male (4 mo, n = %d)", sum(md$age_group == "Young")),
             Old = sprintf("Old male (20 mo, n = %d)", sum(md$age_group == "Old")))
  g$facet <- factor(lab_n[g$age_group], lab_n); ind$facet <- factor(lab_n[ind$age_group], lab_n)

  p_cov <- ggplot(g, aes(pos)) +
    geom_rect(data = hl, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
              fill = "#F2E6B8", alpha = 0.55) +
    geom_line(data = ind, aes(y = norm, group = sample, colour = age_group), linewidth = 0.25, alpha = 0.35) +
    geom_ribbon(aes(ymin = 0, ymax = mean, fill = age_group), alpha = 0.85) +
    geom_ribbon(aes(ymin = pmax(mean - se, 0), ymax = mean + se), fill = "grey10", alpha = 0.18) +
    facet_wrap(~ facet, ncol = 1, strip.position = "top") +
    scale_fill_manual(values = COL, guide = "none") + scale_colour_manual(values = COL, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, ymax), oob = scales::squish) +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    labs(y = "Normalized ATAC signal\n(fragment coverage ÷ size factor)", title = title,
         subtitle = paste0("PRJNA1095555 ATAC-seq reanalysis, male peritoneal macrophages
", stat_lab, "  (Old vs Young)")) +
    base_theme + theme(plot.title = element_text(face = "bold", size = 15),
                       plot.subtitle = element_text(size = 10.5, colour = "grey25"))

  # overlay of the two means
  p_ov <- ggplot(g, aes(pos, mean, colour = age_group)) +
    geom_rect(data = hl, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
              fill = "#F2E6B8", alpha = 0.55) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = COL, name = NULL, labels = c(Young = "Young mean", Old = "Old mean")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, ymax), oob = scales::squish) +
    coord_cartesian(xlim = xlim, expand = FALSE) + labs(y = "Normalized ATAC\nsignal (overlay)") +
    base_theme + theme(legend.position = c(0.99, 0.95), legend.justification = c(1, 1),
                       legend.background = element_blank(), legend.key.height = unit(0.35, "cm"))

  # peaks: ours (coloured by male result) + published TSS peak
  pkv <- pk[pk$end >= xlim[1] & pk$start <= xlim[2], ]
  pkv$cls <- ifelse(pkv$name == "common_0010622", "Timd4 TSS peak (common_0010622)",
                    ifelse(!is.na(pkv$padj) & pkv$padj < 0.05, "other common peak, padj < 0.05", "other common peak, n.s."))
  pkv$y <- 2
  pubv <- transform(pub[pub$end >= xlim[1] & pub$start <= xlim[2], ], y = 1)
  PCOL <- c("Timd4 TSS peak (common_0010622)" = "#B8860B", "other common peak, padj < 0.05" = "grey35",
            "other common peak, n.s." = "grey75")
  p_pk <- ggplot() +
    geom_rect(data = pkv, aes(xmin = start, xmax = end, ymin = y - 0.32, ymax = y + 0.32, fill = cls)) +
    scale_fill_manual(values = PCOL, guide = "none") +
    scale_y_continuous(breaks = if (show_pub) c(1, 2) else 2,
                       labels = if (show_pub) c("Published peaks", "Reanalysis peaks") else "Reanalysis peaks",
                       limits = if (show_pub) c(0.5, 2.5) else c(1.5, 2.5)) +
    coord_cartesian(xlim = xlim, expand = FALSE) + labs(y = NULL) +
    base_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                       axis.text.y = element_text(size = 10))

  if (show_pub) p_pk <- p_pk +
    geom_rect(data = pubv, aes(xmin = start, xmax = end, ymin = y - 0.32, ymax = y + 0.32),
              fill = ifelse(pubv$tss, "grey35", "grey75"))

  # Completeness checks against independently re-read sources: every feature overlapping
  # the window must be drawn, and the window must lie inside the loaded data ranges.
  stopifnot(xlim[1] >= start(REG), xlim[2] <= end(REG), xlim[1] >= 46790000, xlim[2] <= 46860000)
  n_bed <- sum(all_bed$chr == "chr11" & all_bed$end0 >= xlim[1] & all_bed$start0 + 1 <= xlim[2])
  n_sb <- sum(sb$Chr == "chr11" & sb$End >= xlim[1] & sb$Start <= xlim[2])
  if (n_bed != nrow(pkv)) stop(sprintf("%s: %d common peaks overlap the window but %d drawn", file_stem, n_bed, nrow(pkv)))
  if (show_pub && n_sb != nrow(pubv)) stop(sprintf("%s: %d Sheet B peaks overlap the window but %d drawn", file_stem, n_sb, nrow(pubv)))
  el <- rbind(
    data.frame(figure = file_stem, track = "Reanalysis peaks", name = pkv$name, start_1based = pkv$start,
               end = pkv$end, male_log2FC = pkv$log2FoldChange, male_padj = pkv$padj),
    if (show_pub) data.frame(figure = file_stem, track = "Published peaks (Sheet B)", name = pubv$name,
                             start_1based = pubv$start, end = pubv$end, male_log2FC = NA, male_padj = pubv$padj_M))
  el$male_log2FC[el$track != "Reanalysis peaks"] <- NA  # Sheet B log2FC is per month; padj kept as published
  elements[[file_stem]] <<- el

  # gene model
  exv <- as.data.frame(ex); cdsv <- as.data.frame(cds)
  arrows_x <- seq(max(start(tx), xlim[1]) + 600, min(end(tx), xlim[2]) - 300, by = 1500)
  p_gene <- ggplot() +
    geom_segment(aes(x = max(start(tx), xlim[1]), xend = min(end(tx), xlim[2]), y = 0, yend = 0), linewidth = 0.5) +
    geom_segment(data = data.frame(x = arrows_x), aes(x = x, xend = x + 1, y = 0, yend = 0),
                 arrow = arrow(length = unit(0.12, "cm"), type = "open"), linewidth = 0.4, colour = "grey20") +
    geom_rect(data = exv, aes(xmin = start, xmax = end, ymin = -0.18, ymax = 0.18), fill = "grey20") +
    geom_rect(data = cdsv, aes(xmin = start, xmax = end, ymin = -0.32, ymax = 0.32), fill = "grey20") +
    geom_segment(aes(x = TSS, xend = TSS, y = 0.32, yend = 0.75), linewidth = 0.5) +
    geom_segment(aes(x = TSS, xend = TSS + diff(xlim) * 0.025, y = 0.75, yend = 0.75), linewidth = 0.5,
                 arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
    annotate("text", x = TSS + diff(xlim) * 0.03, y = 0.75, label = "Timd4 (Tim4)  TSS", hjust = 0,
             fontface = "italic", size = 3.6) +
    scale_y_continuous(limits = c(-0.5, 1.05)) +
    scale_x_continuous(breaks = TSS + seq(-50, 50, by = step_kb) * 1000,
                       labels = function(x) sub("^-", "−", sprintf("%g", (x - TSS) / 1000))) +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    labs(x = sprintf("Distance from Timd4 TSS (kb)   [chr11:%s–%s, mm10]",
                     format(xlim[1], big.mark = ","), format(xlim[2], big.mark = ",")), y = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          plot.margin = margin(2, 8, 4, 8))

  cap <- paste(
    "Coverage = unshifted paired-end fragments from final filtered BAMs (duplicates, chrM, blacklist removed), 10-bp bins,",
    "divided by DESeq2 size factors of the male common-peak fit (y = number of fragments covering each base, depth-normalized).",
    "Filled area = group mean; dark band = ±SE; thin lines = individual mice.",
    if (show_pub) paste0("Peak tracks: gold = common_0010622 (Timd4 TSS peak), light grey = other reanalysis peaks (all padj > 0.05 or NA in males);\n",
                         "published peaks (Lu et al., Table S2 Sheet B): dark grey = TSS peak, light grey = others.")
    else "Peak track: gold = common_0010622 (Timd4 TSS peak), light grey = other reanalysis peaks (all padj > 0.05 or NA in males).",
    "Yellow band = common_0010622 (chr11:46,810,221–46,812,849). Statistics: DESeq2 ~ age_group, male only (5 Young vs 4 Old), BH over 77,823 peaks.",
    sep = "\n")
  p <- p_cov / p_ov / p_pk / p_gene + plot_layout(heights = c(4.2, 1.6, 0.9, 0.9)) +
    plot_annotation(caption = cap, theme = theme(plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30")))
  ggsave(file.path(out, "figures", paste0(file_stem, ".png")), p, device = ragg::agg_png,
         width = 9, height = 8.2, dpi = 300, bg = "white")
  ggsave(file.path(out, "figures", paste0(file_stem, ".pdf")), p, device = cairo_pdf, width = 9, height = 8.2, bg = "white")
  cat("written", file_stem, "\n")
}

# Full peak BED re-read without region restriction, for the completeness checks
all_bed <- read.delim(file.path(RESULTS, "04_consensus/common_peaks.bed"),
                      header = FALSE, col.names = c("chr", "start0", "end0", "name"))
elements <- list()
for (sp in c(TRUE, FALSE)) {
  sfx <- if (sp) "" else "_no_published"
  make_track(c(TSS - 4000, TSS + 6000), paste0("timd4_locus_track_TSS_zoom", sfx),
             "Tim4 (Timd4) promoter accessibility, male mice", step_kb = 1, show_pub = sp)
  make_track(c(TSS - 6000, TSS + 38000), paste0("timd4_locus_track_whole_gene", sfx),
             "Tim4 (Timd4) locus accessibility, male mice", step_kb = 5, show_pub = sp)
}

write.table(do.call(rbind, elements), file.path(out, "results/figure_elements.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# IGV files: all Sheet B peaks genome-wide (0-based BED, name = PeakName, score unused)
dir.create(file.path(out, "igv"), showWarnings = FALSE)
sbb <- sb[order(sb$Chr, sb$Start), ]
writeLines(c('track name="Lu2025_SheetB_ATAC_peaks" description="Published ATAC peaks (Table S2 Sheet B)" color=90,90,90',
             sprintf("%s\t%d\t%d\t%s_maleFDR=%s", sbb$Chr, as.integer(sbb$Start) - 1L, as.integer(sbb$End), sbb$PeakName,
                     formatC(sbb$padj_M, format = "g", digits = 2))),
           file.path(out, "igv/published_peaks_SheetB.bed"))
cat("igv BED rows:", nrow(sbb), "\n")

writeLines(capture.output(sessionInfo()), file.path(out, "results/sessionInfo.txt"))
