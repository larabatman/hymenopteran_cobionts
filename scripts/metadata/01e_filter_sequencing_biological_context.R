#!/usr/bin/env Rscript
#------
# Build a biologically restricted BioSample dataset
# It selects suitable BioSamples for downstream analysis by enforcing: 
# - PacBio and Hi-C evidence shating at least one BioSample
# - female samples only
# - adult or unkown developmental stage
# - whole body or abdomen tissue
# - excludes head and thorax tissues
# Input: data/annotation/data_inventory_biosample_annotation.tsv
# Output: data/annotation/restricted_dataset.tsv and printed summary statistics
#-----

library(readr)
library(dplyr)
library(stringr)

# Load annotated BioSample data
input_file <- "data/annotation/data_inventory_biosample_annotation.tsv"
output_file <- "data/annotation/restricted_dataset.tsv"

annot <- read_tsv(input_file, show_col_types = FALSE)

# Normalize fields: 
annot <- annot %>%
    mutate(
        sex = str_to_lower(str_trim(sex)),
        tissue = str_to_lower(str_trim(tissue)),
        dev_stage = str_to_lower(str_trim(dev_stage))
    )

# Biological and technical filtering
restricted_df <- annot %>%
    filter(
        # Technical coherence: same biological source for PacBio and Hi-C
        shared_pacbio_hic_biosample == TRUE,
        # Biological constraints
        sex == "female",
        dev_stage == "adult" | is.na(dev_stage) | dev_stage == "",
        str_detect(tissue, "whole|abdomen"),
        !str_detect(tissue, "head|thorax")
    )


# Write restricted dataset
write_tsv(restricted_df, output_file)

# Summary statistics
summary <- tibble(
    n_total_biosamples = nrow(annot),
    n_shared_pacbio_hic = sum(annot$shared_pacbio_hic_biosample == TRUE, na.rm = TRUE),
    n_female = sum(annot$sex == "female", na.rm = TRUE),
    n_adult = sum(annot$dev_stage == "adult" | is.na(annot$dev_stage) | annot$dev_stage == "", na.rm = TRUE),
    n_restricted_biosamples = nrow(restricted_df),
    n_species_restricted = n_distinct(restricted_df$species)
)

cat("\n== Restricted dataset summary ==\n")
print(summary)

cat("\nWrote proper dataset to:\n", output_file, "\n")


