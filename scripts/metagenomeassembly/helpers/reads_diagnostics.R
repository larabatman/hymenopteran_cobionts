#!/usr/bin/env Rscript

# This R script produces before and after read length histograms and generic stats

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

SPECIES <- args[1]
RAW_FILE <- args[2]
FIL_FILE <- args[3]
OUTDIR <- args[4]

raw <- scan(RAW_FILE, quiet=TRUE)
fil <- scan(FIL_FILE, quiet=TRUE)

# Create dataframe with both raw and filtered read lengths
df <- rbind(
  data.frame(length=raw, set="raw"),
  data.frame(length=fil, set="filtered")
)

# Changed histogram binning to binwicdth 0.5 kb, each bin is thus 0.5 kb to be interpretable among species 
p <- ggplot(df, aes(x=length/1000, fill=set)) +
  geom_histogram(binwidth=0.5, alpha=0.6, position="identity") +
  facet_wrap(~set, ncol=1, scales="free_y") +
  theme_bw() +
  labs(
    title=paste("HiFi read length distribution:", SPECIES),
    x="Read length (kb)",
    y="Count"
  )
outfile <- file.path(OUTDIR, "hifi_length_distribution_before_after.png")

ggsave(outfile, p, width=10, height=8, dpi=150)
