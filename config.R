# Input : none (edit paths below or set them as environment variables)
# Output: path variables used by 05-07 R scripts
env <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) v else default }
HERE    <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), winslash = "/"))
SAMPLES <- env("SAMPLES", file.path(HERE, "data/samples.tsv"))
WORK    <- env("WORK", file.path(HERE, "work"))          # BAMs: WORK/<run>/filtered/<run>.primary.filtered.bam
RESULTS <- env("RESULTS", file.path(HERE, "results"))
GTF     <- env("GTF", file.path(HERE, "ref/gencode.vM25.annotation.gtf.gz"))
PUBLISHED_XLSX <- env("PUBLISHED_XLSX", file.path(HERE, "ref/media-2.xlsx"))  # Lu et al. Table S2 (peak track)
bam_path <- function(s) file.path(WORK, s, "filtered", paste0(s, ".primary.filtered.bam"))
