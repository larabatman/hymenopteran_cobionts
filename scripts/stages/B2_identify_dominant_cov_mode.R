#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
INPUT <- args[2] # gc_cov.tsv
OUTDIR <- args[3]

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# gc_cov contains columns conting len gc mean_cov
# Each row is one contig with len, gc and mean_cov
gc_cov <- read_tsv(INPUT, show_col_types = FALSE)

EPS <- 1e-6 # small delta in case dividing by 0
NORM <- 1.4826 # ensure consistency of MAD, assuming a normal distribution (default parameter)

#------------------------------------------------------
# 1) Coverage modeling weighted by contig length
# Log-transform mean_cov and add small delta
# Cov distirbutions can be skewed, and log(0) is undefined
gc_cov <- gc_cov %>%
    filter(!is.na(mean_cov), mean_cov >= 0, !is.na(len), len > 0) %>%
    mutate(logcov = log10(mean_cov + EPS))

# Store as a vector
logcov_vec <- gc_cov$logcov

# Length-weighted statistics: long contgis influence the estimate more (more reliable than short contigs frrom assembly noise, repeats and fragments)
# Define weights as contig length
weights <- gc_cov$len

# Compute host median per-contig mean coverage, weighted by contig length
# Estimate hos coverage center: largest fraction of the assembled DNA is host, and the median coverage across contigs is assumed host
# Median is more robust than mean
host_median <- weightedMedian(logcov_vec, w = weights, na.rm = TRUE)

# Compute host MAD distribution with a normal constant, weighted by contig length
# Measure spread: Median Absolute Deviation, an outlier-robust version of the standard deviation
abs_dev <- abs(logcov_vec - host_median)

# Using weighted median and multiplying by NORM as correction to approximate sd
host_mad <- weightedMedian(abs_dev, w = weights, na.rm = TRUE) * NORM

# Check for MAD == 0
if (host_mad == 0) {
  stop("MAD is zero: coverage distribution too tight or invalid.")
}

#------------------------------------------------------
# 2) Coverage classification
# Compute z-score to measure how far a contig's coverage is from the host coverage
# Compare logcov of all contigs against the MAD distribution from log contigs
# Label cases: 
# - if the difference is smaller than 2 MADs, then it is host like
# - if the difference is smaller than 4 MADs, then it is ambiguous
# - if the difference is greater than 4 MADs, then it is a coverage outlier (cobiont-like)
gc_cov <- gc_cov %>%
    mutate(
        z_cov = (logcov - host_median) / host_mad,
        abs_z = abs(z_cov),
        coverage_class = case_when(
            abs_z < 2 ~ "host_like",
            abs_z < 4 ~ "ambiguous",
            TRUE ~ "coverage_outlier"
        )
    )

#------------------------------------------------------
# 3) Define backbone 
# The host-backbone is within 2 abs_z from the host_mad
host_backbone <- gc_cov %>%
    filter(coverage_class == "host_like")

#------------------------------------------------------
# 4) Summary table
# The backbone percentage is given by the sum of the host_backbone length divided by the total length of the contigs times 100
summary_tbl <- tibble(
    species = SPECIES,
    n_total_contigs = nrow(gc_cov),
    host_median_logcov = host_median,
    host_mad_logcov = host_mad,
    n_backbone = nrow(host_backbone),
    percent_backbone_bp = sum(host_backbone$len) / sum(gc_cov$len) * 100
)

print(summary_tbl)

#------------------------------------------------------
# 5) Writing files
write_tsv(gc_cov, file.path(OUTDIR, "coverage_classification.tsv"))
write_tsv(host_backbone, file.path(OUTDIR, "host_backbone.tsv"))
write_tsv(summary_tbl, file.path(OUTDIR, "coverage_backbone_summary.tsv"))

message("Coverage backbone definition finished.")
