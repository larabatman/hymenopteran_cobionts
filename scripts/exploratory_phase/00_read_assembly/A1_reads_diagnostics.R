#!/usr/bin/env Rscript

# This R script produces before and after read length histograms and generic stats
# It is launched through A1_reads_QC.sh

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)


SPECIES <- args[1]
RAW_FILE <- args[2]
FIL_FILE <- args[3]
OUTDIR <- args[4]

raw <- scan(RAW_FILE, quiet=TRUE)
file <- scan(FIL_FILE, quiet=TRUE)

# Create dataframe with both raw and filtered read lengths
df <- rbind(
  data.frame(length=raw, set="raw"),
  data.frame(length=file, set="filtered")
)
# Each bin is the same size 0.5kbp for interpretability
p <- ggplot(df, aes(x=length/1000, fill=set)) +
  geom_histogram(binwidth=0.5, alpha=0.6, position="identity") +
  facet_wrap(~set, ncol=1, scales="free_y") +
  theme_bw() +
  labs(
    title=paste("HiFi read length distribution:", SPECIES), x="Read length (kb)", y="Count")
outfile <- file.path(OUTDIR, "hifi_length_distribution_before_after.png")
ggsave(outfile, p, width=10, height=8, dpi=150)
