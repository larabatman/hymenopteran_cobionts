#!/usr/bin/env Rscript

# Adds coverage cross-tabulations to the BUSCO contigs
# Usage:
# Rscript C2c_busco_anchor.R species busco_per_contig.tsv outdir

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]
BUSCO <- args[2]
OUTDIR <- args[3]

busco_table <- read_tsv(BUSCO, show_col_types = FALSE)

# Host evidence: which lineages placed BUSCOs on a contig:
# strong_hymenoptera: hymenotpera BUSCOs
# supportive_arthropoda: arthropoda but no hymenoptera BUSCOs
# none
anchor_eval <- busco_table %>%
    mutate(anchor_type = case_when(
        has_hym ~ "strong_hymenoptera",
        has_arth ~ "supportive_arthropoda", 
        TRUE ~ "none"),
        bacterial_signal = has_bact
    )

# Cross tabulation against the coverage classes: do host anchored contigs land in host_like coverage
distribution_host <- anchor_eval %>%
    count(anchor_type, coverage_class) %>%
    group_by(anchor_type) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup()
distribution_bact <- anchor_eval %>%
    count(bacterial_signal, coverage_class) %>%
    group_by(bacterial_signal) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup()

# Tables:
write_tsv(anchor_eval, file.path(OUTDIR, "busco_anchor_full_eval.tsv"))
write_tsv(distribution_host, file.path(OUTDIR, "busco_anchor_coverage_distribution.tsv"))
write_tsv(distribution_bact, file.path(OUTDIR, "busco_bacteria_coverage_distribution.tsv"))