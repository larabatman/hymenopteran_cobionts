#!/usr/bin/env Rscript

# This scripts aims at summarizing the data produced with ncbi.GCA_query.py
# We want to know which species has PacBio, Hi-C, or both, and if they also have RNA-seq data. 

library(readr)
library(dplyr)
library(stringr)

# Load the file
inventory <- read_tsv("data/raw_data_inventory.tsv", show_col_types = FALSE)

# Handling column type
inventory <- inventory %>%
    mutate(
        has_pacbio = as.logical(has_pacbio),
        has_hic = as.logical(has_hic),
        has_rnaseq = as.logical(has_rnaseq),
        shared_pacbio_hic_biosample = as.logical(shared_pacbio_hic_biosample)
    )

# Core counts
counts <- tibble(
    n_species = nrow(inventory),
    n_pacbio = sum(inventory$has_pacbio, na.rm = TRUE),
    n_hic = sum(inventory$has_hic, na.rm = TRUE),
    n_rnaseq = sum(inventory$has_rnaseq, na.rm = TRUE),

    n_pacbio_and_hic = sum(inventory$has_pacbio & inventory$has_hic, na.rm = TRUE),
    n_hic_only = sum(!inventory$has_pacbio & inventory$has_hic, na.rm = TRUE),
    n_pacbio_only = sum(inventory$has_pacbio & !inventory$has_hic, na.rm = TRUE), 

    # Among species with both PacBio and HiC, how many share at least one BioSample
    n_pacbio_hic_shared_biosample = sum(inventory$has_pacbio & inventory$has_hic & inventory$shared_pacbio_hic_biosample, na.rm = TRUE),
    n_pacbio_hic_not_shared_biosample = sum(inventory$has_pacbio & inventory$has_hic & !inventory$shared_pacbio_hic_biosample, na.rm = TRUE)
)

print(counts)

# Split RNA-Seq among the PacBio + HiC group:
pacbio_hic_rnaseq <- inventory %>%
    filter(has_pacbio & has_hic) %>%
    summarise(
        n = n(),
        with_rnaseq = sum(has_rnaseq, na.rm = TRUE),
        without_rnaseq = sum(!has_rnaseq & !is.na(has_rnaseq)),
        shared_biosample = sum(shared_pacbio_hic_biosample, na.rm = TRUE),
        shared_biosample_and_rnaseq = sum(shared_pacbio_hic_biosample & has_rnaseq, na.rm = TRUE)
        )

print(pacbio_hic_rnaseq)

# Species lists of interest
species_pacbio <- inventory %>%
    filter(has_pacbio) %>%
    arrange(species) %>%
    select(species, assembly_accession, taxid)

species_hic <- inventory %>%
    filter(has_hic) %>%
    arrange(species) %>%
    select(species, assembly_accession, taxid)

species_pacbio_hic <- inventory %>%
    filter(has_pacbio & has_hic) %>%
    arrange(desc(shared_pacbio_hic_biosample), species) %>%
    select(species, assembly_accession, taxid, has_rnaseq, shared_pacbio_hic_biosample, pacbio_runs, hic_runs, pacbio_biosamples, hic_biosamples)

# Write down shortlists
write_tsv(species_pacbio, "data/species_with_pacbio.tsv")
write_tsv(species_hic, "data/species_with_hic.tsv")
write_tsv(species_pacbio_hic, "data/species_with_pacbio_and_hic.tsv")

species_interest <- inventory %>%
    filter(has_pacbio & has_hic & shared_pacbio_hic_biosample) %>%
    arrange(species) %>%
    select(species, assembly_accession, taxid, has_rnaseq, pacbio_runs, hic_runs, pacbio_biosamples, hic_biosamples)

write_tsv(species_interest, "data/species_interest_shared_biosamples.tsv")

cat("Wrote: \n",
    "- data/sepcies_with_pacbio.tsv\n",
    "- data/sepcies_with_hic.tsv\n",
    "- data/sepcies_with_pacbio_and_hic.tsv\n",
    "- data/species_interest_shared_biosamples.tsv\n"
    )
