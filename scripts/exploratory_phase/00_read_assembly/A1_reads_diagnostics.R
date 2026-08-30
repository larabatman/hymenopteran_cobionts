#!/usr/bin/env Rscript

# This R script producesbefore and after read length histograms 
# It is launched by A1a_hifi_Qc.sh

# Rscript A1_reads_diagnostics.R species hifi_raw_read_lengths.txt hifi_filtered_read_lengths.txt species/hifi_qc

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)


SPECIES <- args[1]
RAW_FILE <- args[2]
FIL_FILE <- args[3]
OUTDIR <- args[4]

raw <- scan(RAW_FILE, quiet=TRUE)
file <- scan(FIL_FILE, quiet=TRUE)

# Create dataframe with both raw and filtered read lengths
df <- rbind(data.frame(length=raw, set="raw"), data.frame(length=file, set="filtered")) # define the raw and filtered sets

# Each bin is the same size 0.5kbp for interpretability
p <- ggplot(df, aes(x=length/1000, fill=set)) +
  geom_histogram(binwidth=0.5, alpha=0.6, position="identity") +
  facet_wrap(~set, ncol=1, scales="free_y") + # split plot into separate panels, feaceting formula with one value for each set column to have raw and filtered panels. Stacked into a single column with ncol = 1 and let each panel pick the y axis range with free_y
  theme_bw() +
  labs(title=paste("HiFi read length distribution:", SPECIES), x="Read length (kbp)", y="Count")
outfile <- file.path(OUTDIR, "hifi_length_distribution_before_after.pdf")
ggsave(outfile, p, width=10, height=8, dpi=150)
