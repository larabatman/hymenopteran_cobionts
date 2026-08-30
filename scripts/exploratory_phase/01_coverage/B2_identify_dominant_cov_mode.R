#!/usr/bin/env Rscript

# Create a model of coverage distribution: length-weighted mean absolute deviation defines the host backbone
# Classification into host_like, ambiguous and coverage_outlier
# Usage:
# Rscript B1b_identify_dominant_cov_mode.R species gc_cov.tsv outdir

library(dplyr)
library(readr)
# weightedMedian() from matrixStats package
library(matrixStats)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
INPUT <- args[2]
OUTDIR <- args[3]

# Read gc_cov.tsv table 
gc_cov <- read_tsv(INPUT, show_col_types = FALSE)
# Define a small error in case of log(0)
EPS <- 1e-6
# Define the constant for normal distribution and z interpretability
NORM <- 1.4826

# Modelling coverage weighted by contig length
n_raw <- nrow(gc_cov) # count rows
gc_cov <- gc_cov %>%
  filter(!is.na(mean_cov), is.finite(mean_cov), mean_cov >= 0, !is.na(len), len > 0) %>% # Filter out mean_cov that are smaller than 0, lengths that are NA or smaller than 0
  mutate(logcov = log10(mean_cov + EPS)) %>% # log10 transform raw coverage, add the small error for log(0)
  filter(is.finite(logcov)) # filter out infinite coverage
# Trace dropped contigs
n_dropped <- n_raw - nrow(gc_cov)
if (n_dropped > 0)
  message(sprintf("Dropped contigs: ", n_dropped))

# Select log coverage and their lengths, which will be their weights
logcov_vec <- gc_cov$logcov
weights <- gc_cov$len

# Host coverage center as the median: weightedMedian finds the coverage value at which half of the total assembly bp are below it, and half are above it
host_median <- weightedMedian(logcov_vec, w = weights, na.rm = TRUE)

# Median Absolute Deviation: spread as standard deviation using medians instead of means 
abs_dev  <- abs(logcov_vec - host_median)
# Scaling by the normality constant for interpretability using the estimated of the standard deviation under a normal distribution
host_mad <- weightedMedian(abs_dev, w = weights, na.rm = TRUE) * NORM

# Classification: computing z scores and appreciating the number of MADs away from the host median
# Since in log scale: z=2 means 2 MADs above the log median 
# |z| < 2: host_like coverage, interpreted as 95% of spread within normal distribution 
# 2 ≤ |z| < 4: in between host_like and coverage_outlier, ambiguous
# |z| ≥ 4: coverage_outlier
# |z| > 6: is_extreme
# Recording the direction of deviation with cov_direction: above or below log median 
gc_cov <- gc_cov %>%
  mutate(
    z_cov = (logcov - host_median) / host_mad,
    abs_z = abs(z_cov),
    coverage_class = case_when(abs_z < 2 ~ "host_like", abs_z < 4 ~ "ambiguous", TRUE ~ "coverage_outlier"),
    is_extreme = abs_z > 6,
    cov_direction = case_when(coverage_class == "host_like" ~ "neutral", z_cov > 0 ~ "high", z_cov < 0 ~ "low"))

# The host backbone are the contigs that are host_like
host_backbone <- gc_cov %>% 
  filter(coverage_class == "host_like")

# Write summary table
summary_table <- tibble(
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
print(summary_table)

# Blobplot: GC against coverage coloured by coverage class
# Filter out GC that are below 13 or above 80%
df_plot <- gc_cov %>% filter(mean_cov > 0, len > 0, gc >= 13, gc <= 80)
# Each point is a contig, with three variables: x as GC content, y as mean coverage, point size as log10 contig length 
# shape = 21 to fill the circles with a different border color
p_blobplot <- ggplot(df_plot, aes(gc, mean_cov)) +
  geom_point(aes(size = log10(pmax(len, 1)), fill = coverage_class), shape = 21, color = "black", alpha = 0.4) +
  scale_y_log10() +
  scale_fill_manual(values = c(host_like = "steelblue", coverage_outlier = "firebrick",ambiguous = "grey70"), na.value = "grey85") +
  theme_bw() +
  labs(title = paste("GC vs Coverage:", SPECIES), x = "GC (%)", y = "Mean coverage", fill = "Coverage class", size = "log10(contig length)")

ggsave(file.path(OUTDIR, "blobplot_coverage_class.pdf"), p_blobplot, width = 7, height = 5, dpi = 300)

# Coverage histogram with dashed vertical lines at host median 
p_hist <- ggplot(df_plot, aes(x = log10(mean_cov), fill = coverage_class)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  scale_fill_manual(values = c(host_like = "steelblue", coverage_outlier = "firebrick", ambiguous = "grey70")) +
  geom_vline(xintercept = host_median, linetype = "dashed", color = "black") +
  theme_bw()

ggsave(file.path(OUTDIR, "coverage_histogram.pdf"), p_hist, width = 7, height = 5, dpi = 300)

# Tables
write_tsv(gc_cov, file.path(OUTDIR, "coverage_classification.tsv"))
write_tsv(host_backbone, file.path(OUTDIR, "host_backbone.tsv"))
write_tsv(summary_tbl, file.path(OUTDIR, "coverage_backbone_summary.tsv"))