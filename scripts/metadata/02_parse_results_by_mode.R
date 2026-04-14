#!/usr/bin/env Rscript
# parse_inventory.R
# Reads the raw data inventory, deduplicates to one row per species,
# cross-references biosamples between HiFi and Hi-C,
# classifies by assembly mode, and flags species needing manual decisions.
#
# Usage: Rscript parse_inventory.R <inventory.tsv> <output_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript parse_inventory.R <inventory.tsv> <output_dir>")
}

infile  <- args[1]
outdir  <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Read and deduplicate to one row per species
# ------------------------------------------------------------------

raw <- read_tsv(infile, show_col_types = FALSE)

cat("\n========== RAW INVENTORY ==========\n")
cat("Total rows:", nrow(raw), "\n")
cat("Unique species:", n_distinct(raw$species), "\n")

# Runs and biosamples are identical across duplicate rows for the same species,
# so we just take the first row per species.
samples <- raw %>%
  group_by(species) %>%
  summarise(
    assembly_accession        = first(assembly_accession),
    taxid                     = first(taxid),
    has_pacbio                = any(has_pacbio == TRUE | has_pacbio == "True"),
    has_hic                   = any(has_hic == TRUE | has_hic == "True"),
    has_rnaseq                = any(has_rnaseq == TRUE | has_rnaseq == "True"),
    pacbio_runs               = first(pacbio_runs),
    hic_runs                  = first(hic_runs),
    pacbio_biosamples         = first(pacbio_biosamples),
    hic_biosamples            = first(hic_biosamples),
    shared_pacbio_hic_biosample = first(shared_pacbio_hic_biosample),
    n_biosamples              = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------------
# 2. Count runs, find shared biosamples
# ------------------------------------------------------------------

count_items <- function(s) {
  ifelse(is.na(s) | s == "", 0L, str_count(s, ",") + 1L)
}

split_items <- function(s) {
  if (is.na(s) || s == "") return(character(0))
  str_split(s, ",")[[1]] %>% str_trim()
}

# Find the intersection of PacBio and Hi-C biosamples
find_shared_biosample <- function(pb_bs, hic_bs) {
  pb  <- split_items(pb_bs)
  hic <- split_items(hic_bs)
  shared <- intersect(pb, hic)
  if (length(shared) == 0) return(NA_character_)
  paste(shared, collapse = ",")
}

samples <- samples %>%
  mutate(
    n_pacbio_runs = count_items(pacbio_runs),
    n_hic_runs    = count_items(hic_runs),
    # Identify which biosample(s) are shared between HiFi and Hi-C
    matched_biosample = mapply(
      find_shared_biosample, pacbio_biosamples, hic_biosamples
    ),
    biosample_match = !is.na(matched_biosample)
  )

# ------------------------------------------------------------------
# 3. Classify assembly mode and flag issues
# ------------------------------------------------------------------

samples <- samples %>%
  mutate(
    # Assembly mode
    asm_mode = case_when(
      !has_pacbio                        ~ "skip",
      has_pacbio & has_hic & biosample_match  ~ "hic",
      has_pacbio & has_hic & !biosample_match ~ "hic_unmatched",
      has_pacbio & !has_hic              ~ "nohic",
      TRUE                               ~ "skip"
    ),
    # Flags (multiple can apply, separated by ;)
    flag = {
      flags <- character(nrow(.))
      for (i in seq_len(nrow(.))) {
        f <- c()
        if (!.$has_pacbio[i])                   f <- c(f, "NO_PACBIO")
        if (.$n_pacbio_runs[i] > 1)             f <- c(f, "MULTI_HIFI")
        if (.$has_pacbio[i] & .$has_hic[i] &
            !.$biosample_match[i])              f <- c(f, "BIOSAMPLE_MISMATCH")
        if (length(f) == 0) f <- "OK"
        flags[i] <- paste(f, collapse = ";")
      }
      flags
    }
  )

# ------------------------------------------------------------------
# 4. Summary
# ------------------------------------------------------------------

cat("\n========== SUMMARY ==========\n")
cat("\nSpecies by assembly mode:\n")
print(table(samples$asm_mode))

cat("\nBiosample matching (species with both HiFi + Hi-C):\n")
both <- samples %>% filter(has_pacbio & has_hic)
cat("  Matched biosample:  ", sum(both$biosample_match), "\n")
cat("  Unmatched biosample:", sum(!both$biosample_match), "\n")

cat("\nFlag breakdown:\n")
# Expand flags for counting
all_flags <- unlist(str_split(samples$flag, ";"))
print(sort(table(all_flags), decreasing = TRUE))

cat("\nSpecies with multiple HiFi runs:\n")
multi <- samples %>% filter(n_pacbio_runs > 1)
cat("  Count:", nrow(multi), "\n")
cat("  Run counts:", paste(sort(unique(multi$n_pacbio_runs)), collapse = ", "), "\n")

# ------------------------------------------------------------------
# 5. Write outputs
# ------------------------------------------------------------------

cols_full <- c("species", "assembly_accession", "taxid", "asm_mode",
               "has_pacbio", "has_hic", "biosample_match",
               "n_pacbio_runs", "n_hic_runs",
               "pacbio_runs", "hic_runs",
               "pacbio_biosamples", "hic_biosamples", "matched_biosample",
               "flag")

# Full sample sheet
out_full <- file.path(outdir, "sample_sheet_full.tsv")
samples %>%
  select(all_of(cols_full)) %>%
  arrange(asm_mode, species) %>%
  write_tsv(out_full)
cat("\nFull sample sheet:", out_full, "\n")

# Ready to go: matched biosample, single HiFi run
out_ready <- file.path(outdir, "sample_sheet_ready.tsv")
samples %>%
  filter(asm_mode %in% c("hic", "nohic"), flag == "OK") %>%
  select(species, assembly_accession, taxid, asm_mode,
         pacbio_runs, hic_runs, matched_biosample) %>%
  arrange(asm_mode, species) %>%
  write_tsv(out_ready)
cat("Ready to run (no flags):", out_ready, "\n")

# Multi-HiFi: user needs to pick one run
out_multi_hifi <- file.path(outdir, "sample_sheet_multi_hifi.tsv")
samples %>%
  filter(str_detect(flag, "MULTI_HIFI")) %>%
  select(species, asm_mode, n_pacbio_runs, pacbio_runs,
         pacbio_biosamples, hic_biosamples,
         biosample_match, matched_biosample, flag) %>%
  arrange(desc(n_pacbio_runs), species) %>%
  write_tsv(out_multi_hifi)
cat("Multi-HiFi (choose one run):", out_multi_hifi, "\n")

# Biosample mismatch: HiFi and Hi-C from different individuals
out_mismatch <- file.path(outdir, "sample_sheet_biosample_mismatch.tsv")
samples %>%
  filter(str_detect(flag, "BIOSAMPLE_MISMATCH")) %>%
  select(species, asm_mode, pacbio_runs, hic_runs,
         pacbio_biosamples, hic_biosamples, flag) %>%
  arrange(species) %>%
  write_tsv(out_mismatch)
cat("Biosample mismatch:", out_mismatch, "\n")

# Skipped (no PacBio)
out_skip <- file.path(outdir, "sample_sheet_skip.tsv")
samples %>%
  filter(asm_mode == "skip") %>%
  select(species, has_hic, hic_runs, hic_biosamples, flag) %>%
  arrange(species) %>%
  write_tsv(out_skip)
cat("Skipped (no PacBio):", out_skip, "\n")

# Species lists for batch scripting
out_sp_ready <- file.path(outdir, "species_list_ready.txt")
samples %>%
  filter(asm_mode %in% c("hic", "nohic"), flag == "OK") %>%
  pull(species) %>% sort() %>%
  writeLines(out_sp_ready)
cat("Ready species list:", out_sp_ready, "\n")

cat("\n========== DONE ==========\n")