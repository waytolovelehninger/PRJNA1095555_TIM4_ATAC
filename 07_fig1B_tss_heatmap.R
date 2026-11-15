#!/usr/bin/env Rscript
# Input : $RESULTS/07_fig1B/matrix/TSS_2kb_50bp.male.mat.gz and logs (07_fig1B_bigwig.sh), 06_fig1A per-sample coverage
# Output: $RESULTS/07_fig1B/figures/TSS_genomewide_heatmap_male_v3.png|pdf, results/genomewide_TSS_summary.tsv (Fig. 1B)
options(stringsAsFactors = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "config.R"))
out <- file.path(RESULTS, "07_fig1B")
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid); library(ggplot2); library(patchwork) })
dir.create(file.path(out, "results"), showWarnings = FALSE)
dir.create(file.path(out, "figures"), showWarnings = FALSE)
bwdir <- out
stopifnot(file.exists(file.path(bwdir, "SUCCESS")))

# ---- 1. bigWig cross-check against the R fragment coverage of the Timd4 window ---------------
bwv <- read.delim(file.path(bwdir, "logs/timd4_window_bigwig_values.tsv"))
rcv <- read.delim(file.path(RESULTS, "06_fig1A/results/timd4_coverage_10bp_per_sample.tsv"))
bwv$bin_start <- round(bwv$pos_center - 4.5); rcv$bin_start <- round(rcv$pos - 5)
cc <- merge(bwv, rcv[, c("sample", "bin_start", "norm")], by = c("sample", "bin_start"))
chk <- do.call(rbind, lapply(split(cc, cc$sample), function(x)
  data.frame(sample = x$sample[1], bins = nrow(x), pearson = cor(x$bw_value, x$norm),
             mean_bigwig = mean(x$bw_value), mean_R = mean(x$norm), max_abs_diff = max(abs(x$bw_value - x$norm)))))
write.table(chk, file.path(out, "results/bigwig_vs_R_coverage_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(chk)
# deepTools counts a fragment fully in every 10-bp bin it touches, the R coverage averages per base,
# so bigWig values run slightly higher; require near-identical shape and <5% mean difference.
stopifnot(nrow(chk) == 9L, all(chk$bins == 5200L), all(chk$pearson > 0.995),
          all(abs(chk$mean_bigwig / chk$mean_R - 1) < 0.05))

# ---- 2. read matrix -------------------------------------------------------------------------
mf <- file.path(bwdir, "matrix/TSS_2kb_50bp.male.mat.gz")
hdr <- jsonlite::fromJSON(sub("^@", "", readLines(gzfile(mf), n = 1)))
x <- read.delim(gzfile(mf), header = FALSE, skip = 1)
reg <- x[, 1:6]; names(reg) <- c("chr", "start", "end", "gene_id", "score", "strand")
v <- as.matrix(x[, -(1:6)]); v[is.na(v)] <- 0
labs <- hdr$sample_labels; bnd <- hdr$sample_boundaries
nb <- diff(bnd)[1]
stopifnot(length(labs) == 9L, all(diff(bnd) == nb), ncol(v) == 9L * nb, nb == 80L, nrow(v) > 20000)
grp <- sub("^.*_", "", labs); stopifnot(sum(grp == "Young") == 5L, sum(grp == "Old") == 4L)
sl <- lapply(seq_along(labs), function(i) v[, (bnd[i] + 1):bnd[i + 1]])
Ym <- Reduce(`+`, sl[grp == "Young"]) / 5
Om <- Reduce(`+`, sl[grp == "Old"]) / 4
pos <- seq(-2000 + 25, 2000 - 25, by = 50)

# per-gene TSS +-500 bp signal and Old/Young log2 ratio (descriptive, not a test)
core <- which(abs(pos) <= 500)
gs <- data.frame(gene_id = reg$gene_id, young = rowMeans(Ym[, core]), old = rowMeans(Om[, core]))
gs$log2_old_young <- log2((gs$old + 1) / (gs$young + 1))
expressed <- (gs$young + gs$old) / 2 >= 5          # TSSs with appreciable signal
tim <- which(gs$gene_id == "ENSMUSG00000055546.6"); stopifnot(length(tim) == 1L)
rank_pct <- mean(gs$log2_old_young[expressed] <= gs$log2_old_young[tim]) * 100
summ <- data.frame(
  n_TSS = nrow(gs), n_TSS_signal_ge5 = sum(expressed),
  median_log2_old_young_signal_ge5 = median(gs$log2_old_young[expressed]),
  Timd4_young = gs$young[tim], Timd4_old = gs$old[tim], Timd4_log2_old_young = gs$log2_old_young[tim],
  Timd4_percentile_among_signal_ge5 = rank_pct)
write.table(summ, file.path(out, "results/genomewide_TSS_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gs, file.path(out, "results/per_gene_TSS500_signal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

# ---- 3. drop low-signal TSSs ----------------------------------------------------------------
# TSSs with TSS +-500 bp mean signal < 5 cannot be judged (ratio dominated by noise). They used to be
# greyed out in the ratio panel; now they are removed from the profile and from every heatmap.
keep <- which(expressed)
n_drop <- nrow(gs) - length(keep)
Yk <- Ym[keep, ]; Ok <- Om[keep, ]
tim_k <- match(tim, keep); stopifnot(!is.na(tim_k))
fmt <- function(n) format(n, big.mark = ",")

# ---- 4. log2 (Old / Young) heatmap (rows sorted by mean signal) ------------------------------
ht_opt$message <- FALSE
COL <- c(Young = "#3B6FB6", Old = "#D1495B")
ord <- order(-rowMeans((Yk + Ok) / 2))
col_lr <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
LR <- log2((Ok + 1) / (Yk + 1))
tim_row <- match(tim_k, ord)
HW_CM <- 7
TICK <- -2:2; TLAB <- c("−2", "−1", "TSS", "+1", "+2")
ht <- Heatmap(LR[ord, ], name = "lr", col = col_lr, cluster_rows = FALSE, cluster_columns = FALSE,
              use_raster = TRUE, raster_quality = 4, show_row_names = FALSE, show_column_names = FALSE,
              column_title = "log2 (Old / Young)", column_title_gp = gpar(fontsize = 10.5, fontface = "bold"),
              heatmap_legend_param = list(title = "log2\n(Old+1)/(Young+1)", title_gp = gpar(fontsize = 8.5), labels_gp = gpar(fontsize = 8)),
              bottom_annotation = HeatmapAnnotation(ax_lr = anno_empty(border = FALSE, height = unit(5, "mm")), show_annotation_name = FALSE),
              right_annotation = rowAnnotation(mark = anno_mark(at = tim_row, labels = "Tim4 (Timd4)",
                                                                labels_gp = gpar(fontsize = 9.5, fontface = "bold"))),
              width = unit(HW_CM, "cm"), height = unit(11, "cm"))

# x axis in kb, drawn in a viewport whose native x scale is the +-2 kb span of the heatmap columns
xaxis <- function(fs = 8) {
  grid.segments(unit(TICK, "native"), unit(1, "npc"), unit(TICK, "native"), unit(1, "npc") - unit(1.2, "mm"))
  grid.text(TLAB, unit(TICK, "native"), unit(1, "npc") - unit(1.8, "mm"), just = "top", gp = gpar(fontsize = fs))
}

# ---- 5. profile (same TSS set) and figure ---------------------------------------------------
yp <- colMeans(Yk); op <- colMeans(Ok)
yt <- pretty(c(0, max(yp, op))); yr <- range(yt)

W <- 7; H <- 10.3                                         # inches
CAP_CHARS <- 128                                          # caption wrap width at 7.8 pt for a 7-in page
capt <- c(
  paste("Male peritoneal macrophages, PRJNA1095555 ATAC-seq reanalysis. One TSS per protein-coding gene (GENCODE M25), ±2 kb, 50-bp bins.",
        "Signal = fragment coverage (deepTools bamCoverage, read1 extended to fragment, 10-bp bins) divided by DESeq2 size factors of the male",
        "common-peak fit; group means (Young n = 5, Old n = 4). Heatmap = log2((Old+1)/(Young+1)) per bin, rows sorted by mean signal;",
        "descriptive, not a statistical test."),
  paste(sprintf("Of %s protein-coding TSSs, %s with TSS ±500 bp mean signal < 5 were removed (too weak to judge; their ratio is dominated by noise).",
                fmt(nrow(gs)), fmt(n_drop)),
        sprintf("The profile and the heatmap show the remaining %s TSSs; the profile x-axis has the same ±2 kb span and width as the heatmap.",
                fmt(length(keep)))),
  sprintf("TSS ±500 bp: median log2(Old/Young) = %.2f across these %s TSSs; Tim4 (Timd4) = %.2f (%.1f%% of TSSs change as much or more negatively).",
          summ$median_log2_old_young_signal_ge5, fmt(summ$n_TSS_signal_ge5), summ$Timd4_log2_old_young, rank_pct))
capt <- paste(unlist(lapply(capt, strwrap, width = CAP_CHARS)), collapse = "\n")
title <- sprintf("Genome-wide TSS accessibility is unchanged with age;\nTim4 is in the most decreased %.1f%%", rank_pct)

fig <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text(title, x = unit(0.17, "in"), y = unit(H - 0.2, "in"), just = c("left", "top"),
            gp = gpar(fontsize = 13, fontface = "bold", lineheight = 1.1))
  # heatmap
  pushViewport(viewport(x = 0, y = unit(1.55, "in"), width = 1, height = unit(5.4, "in"), just = c("left", "bottom")))
  draw(ht, newpage = FALSE, heatmap_legend_side = "right",
       column_title = "Distance from TSS (kb); rows = TSSs sorted by mean signal",
       column_title_side = "bottom", column_title_gp = gpar(fontsize = 9.5))
  decorate_annotation("ax_lr", { pushViewport(viewport(xscale = c(-2, 2))); xaxis(); popViewport() })
  seekViewport("lr_heatmap_body_1_1")
  lr <- deviceLoc(unit(c(0, 1), "npc"), unit(c(0, 1), "npc"), valueOnly = TRUE)
  upViewport(0)
  stopifnot(abs(diff(lr$x) * 2.54 - HW_CM) < 1e-3)        # heatmap body is exactly HW_CM wide
  # profile directly above the heatmap: same left/right edges and the same +-2 kb scale
  PH <- 1.85; pb <- lr$y[2] + 0.75
  stopifnot(pb + PH + 0.3 < H - 0.75)
  pushViewport(viewport(x = unit(lr$x[1], "in"), y = unit(pb, "in"), width = unit(diff(lr$x), "in"), height = unit(PH, "in"),
                        just = c("left", "bottom"), xscale = c(-2, 2), yscale = yr, clip = "off"))
  grid.lines(unit(pos / 1000, "native"), unit(yp, "native"), gp = gpar(col = COL[["Young"]], lwd = 2))
  grid.lines(unit(pos / 1000, "native"), unit(op, "native"), gp = gpar(col = COL[["Old"]], lwd = 2))
  grid.lines(unit(c(0, 0, 1), "npc"), unit(c(1, 0, 0), "npc"))
  pushViewport(viewport(y = 0, height = 0, just = "bottom", xscale = c(-2, 2))); xaxis(); popViewport()
  grid.text("Distance from TSS (kb)", y = unit(-5.6, "mm"), just = "top", gp = gpar(fontsize = 8.5))
  grid.segments(unit(0, "npc"), unit(yt, "native"), unit(-1.2, "mm"), unit(yt, "native"))
  grid.text(yt, unit(-1.8, "mm"), unit(yt, "native"), just = "right", gp = gpar(fontsize = 8))
  grid.text("Mean normalized\nATAC signal", x = unit(-9, "mm"), rot = 90, just = "bottom", gp = gpar(fontsize = 9))
  grid.text(sprintf("Mean TSS profile (%s TSSs with signal ≥ 5)", fmt(length(keep))), x = 0, y = unit(1, "npc") + unit(2.5, "mm"),
            just = c("left", "bottom"), gp = gpar(fontsize = 10.5, fontface = "bold"))
  # legend to the right of the profile (outside the panel, clear of the peak)
  leg <- c(Young = "Young mean (n = 5)", Old = "Old mean (n = 4)")
  for (i in seq_along(leg)) {
    yy <- unit(1, "npc") - unit(2 + (i - 1) * 4.5, "mm")
    grid.segments(unit(1, "npc") + unit(5, "mm"), yy, unit(1, "npc") + unit(11, "mm"), yy, gp = gpar(col = COL[[names(leg)[i]]], lwd = 2.4))
    grid.text(leg[[i]], unit(1, "npc") + unit(13, "mm"), yy, just = "left", gp = gpar(fontsize = 8.5))
  }
  popViewport()
  # caption
  grid.text(capt, x = unit(0.17, "in"), y = unit(1.32, "in"), just = c("left", "top"),
            gp = gpar(fontsize = 7.8, col = "grey30", lineheight = 1.1))
}
fb <- file.path(out, "figures/TSS_genomewide_heatmap_male_v3")
ragg::agg_png(paste0(fb, ".png"), width = W, height = H, units = "in", res = 300, background = "white"); fig(); invisible(dev.off())
cairo_pdf(paste0(fb, ".pdf"), width = W, height = H, bg = "white"); fig(); invisible(dev.off())
writeLines(capture.output(sessionInfo()), file.path(out, "results/sessionInfo.txt"))
cat("done\n")
