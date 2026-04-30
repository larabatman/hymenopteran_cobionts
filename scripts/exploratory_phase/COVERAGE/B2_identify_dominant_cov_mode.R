#!/usr/bin/env Rscript

# =============================================================================
# Stage B2: Coverage modeling and host backbone definition
#
# Models the coverage distribution using length-weighted MAD to define the
# host backbone. Produces coverage classification, backbone table, summary
# statistics, and a base blobplot (GC vs coverage, colored by coverage class).
#
# Usage:
#   Rscript B1b_identify_dominant_cov_mode.R <SPECIES> <GC_COV_TSV> <OUTDIR>
# After assembly: mixed contigs from host genome and potential cobionts. 
# Host contigs will cluster together in coverage space: they are all replicated once per cell, and should have the same sequencing depth.
# Cobionts, at higher or lower abundance, appear at different coverage. 
# Modelling that host coverage cluster staitsitcally to draw boundary and label host vs non-host 

# Does not assume coverage distribution to be normal, but uses median and AD for mixture of distributions
# Normality is assumed when scaling the MAD to make z-score thresholds interpretable

library(dplyr)
library(readr)
library(matrixStats) # provides weightedMedian() needed for length-weighted statistics
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
INPUT   <- args[2]   # gc_cov.tsv
OUTDIR  <- args[3]

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

gc_cov <- read_tsv(INPUT, show_col_types = FALSE)

EPS  <- 1e-6    # small delta for log(0)
NORM <- 1.4826  # MAD consistency constant for normal distribution


# 1) Coverage modeling weighted by contig length
# Filter out >= 0 mean_cov, NA lengths and negative or 0 lengths
# Log10-transform the raw coverage as it is right-skewed and spans orders of magnitudes 
# On a linear scale, the host cluster gets copressed and outliers dominate
# Avoid small error to avoid log10(0) = - Inf for any 0-cov contigs that survied the filter
gc_cov <- gc_cov %>%
  filter(!is.na(mean_cov), mean_cov >= 0, !is.na(len), len > 0) %>%
  mutate(logcov = log10(mean_cov + EPS))

logcov_vec <- gc_cov$logcov
weights <- gc_cov$len

# Length-weighted median: host coverage center
# Assembly has many contigs of different lengths; length.weighted median finds the coverage calue such that half of the total assembly bp are below it and half are above it
host_median <- weightedMedian(logcov_vec, w = weights, na.rm = TRUE)

# Length-weighted MAD with normal consistency constant
# Median Absolute Deviation measures the spread as standard deviation, but using medians instead of means, such that it is robust to outliers
# Distirbution is not Gaussian, but tailed 
abs_dev  <- abs(logcov_vec - host_median)
# NORM constnat rescales MAD to be a consistent estimator of the standard deviation under a normal distirbution, meaning that if data was truly normal, host_mad = sd()
# Makes the z-scores interpretable 
host_mad <- weightedMedian(abs_dev, w = weights, na.rm = TRUE) * NORM 

if (host_mad == 0) {
  stop("MAD is zero: coverage distribution too tight or invalid.")
}


# 2) Coverage classification
# Compute z-score: how many MADs away from the host median is the contig?
# log-scale: z=2 means 2 MADs above the log-median
# |z| < 2: within core host coverage cluster, captures 95% or a normal distirbution
# 2 ≤ |z| < 4: outside the host core, but not clearly foreign
# |z| ≥ 4: different coverage, outlier
# Adding an is_extreme flag for contigs |z| > 6, highly abundant cobionts or assembly artifacts for downstream filtering
# cov_direction: which way the outlier relative to the log median. High-coverage outliers are potential cobionts, while low-coverage outliers are more likely sequencing arrtifacts
gc_cov <- gc_cov %>%
  mutate(
    z_cov = (logcov - host_median) / host_mad,
    abs_z = abs(z_cov),
    coverage_class = case_when(
      abs_z < 2 ~ "host_like",
      abs_z < 4 ~ "ambiguous",
      TRUE      ~ "coverage_outlier"
    ),
    is_extreme = abs_z > 6,
    cov_direction = case_when(
      coverage_class == "host_like" ~ "neutral",
      z_cov > 0 ~ "high",
      z_cov < 0 ~ "low",
      TRUE      ~ "neutral"
    )
  )

# 3) Define backbone

host_backbone <- gc_cov %>% filter(coverage_class == "host_like")

# 4) Summary table

summary_tbl <- tibble(
  species = SPECIES,
  n_total_contigs = nrow(gc_cov),
  host_median_logcov = host_median,
  host_mad_logcov = host_mad,
  n_backbone = nrow(host_backbone),
  n_ambiguous = sum(gc_cov$coverage_class == "ambiguous"),
  n_coverage_outlier = sum(gc_cov$coverage_class == "coverage_outlier"),
  n_extreme = sum(gc_cov$is_extreme),
  percent_backbone_bp = sum(host_backbone$len) / sum(gc_cov$len) * 100
)

print(summary_tbl)


# 5) Base blobplot: GC vs coverage, colored by coverage class
# Filtering out GC < 13% & GC > 80% as artifacts biologically implausible
df_plot <- gc_cov %>% filter(mean_cov > 0, len > 0, gc >= 13, gc <= 80)
# Each point is a contig
# Three variables:
# x: GC content
# y: mean coverage
# Point size: log10(contig length)
# Fill color: coverage class
# shape = 21: filled circle with a separate border color 
# pmax(len, 1) guards against log10(0) for any zero-length edge cases before size mapping
# scale_y_log10() for y on log scale
p <- ggplot(df_plot, aes(gc, mean_cov)) +
  geom_point(
    aes(size = log10(pmax(len, 1)), fill = coverage_class),
    shape = 21, color = "black", alpha = 0.4
  ) +
  scale_y_log10() +
  scale_fill_manual(
    values = c(
      host_like        = "steelblue",
      coverage_outlier = "firebrick",
      ambiguous        = "grey70"
    ),
    na.value = "grey85"
  ) +
  theme_bw() +
  labs(
    title = paste("GC vs Coverage:", SPECIES),
    x     = "GC (%)",
    y     = "Mean coverage",
    fill  = "Coverage class",
    size  = "log10(contig length)"
  )

ggsave(file.path(OUTDIR, "blobplot_coverage_class.png"), p, width = 7, height = 5, dpi = 300)
message("[INFO] Written: ", file.path(OUTDIR, "blobplot_coverage_class.png"))

# Coverage histogram
# Add dashed vertical line at host-median to visually confirm that the dashed line lands at the peak of the blue, host_like distirbution
p_hist <- ggplot(df_plot, aes(x = log10(mean_cov), fill = coverage_class)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      host_like        = "steelblue",
      coverage_outlier = "firebrick",
      ambiguous        = "grey70"
    )
  ) +
  geom_vline(xintercept = host_median, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    title = paste("Coverage distribution:", SPECIES),
    subtitle = sprintf("Host median = %.2f log10(cov), MAD = %.3f", host_median, host_mad),
    x    = "log10(mean coverage)",
    y    = "Count",
    fill = "Coverage class"
  )

ggsave(file.path(OUTDIR, "coverage_histogram.png"), p_hist, width = 7, height = 5, dpi = 300)
message("[INFO] Written: ", file.path(OUTDIR, "coverage_histogram.png"))

# 6) Write files

write_tsv(gc_cov, file.path(OUTDIR, "coverage_classification.tsv"))
write_tsv(host_backbone, file.path(OUTDIR, "host_backbone.tsv"))
write_tsv(summary_tbl, file.path(OUTDIR, "coverage_backbone_summary.tsv"))

message("Coverage backbone definition finished.")