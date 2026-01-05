#!/usr/bin/env Rscript

#----- 
# Summary of sequencing availablility across species
# PacBio and Hi-C availability is the primary filter. The evidence strength refines the interpretation, and BioSample sharing is used as biological coherence signal.

#Input: data/data.inventory.tsv
#Output: a printed summary tample and opotional TSV exports for downstream inspection
#-----

# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Load inventory file
dir.create("data/inventory_analysis", showWarnings = FALSE)
inventory <- read_tsv("data/data_inventory.tsv", show_col_types = FALSE)

# Type normalization: columns are handled as logical signals rather than strings
inventory <- inventory %>%
    mutate(
        has_pacbio = as.logical(has_pacbio),
        has_hic = as.logical(has_hic),
        has_rnaseq = as.logical(has_rnaseq),
        shared_pacbio_hic_biosample = as.logical(shared_pacbio_hic_biosample),
        run_mapping_level = factor(
            run_mapping_level, 
            levels = c("BIOPROJECT", "BIOSAMPLE", "TAXID_FALLBACK", "NONE"))
    )

# Summary table: 
summary <- inventory %>%
    group_by(run_mapping_level)%>%
    summarise(
        n_species = n(),
        # Primary signals
        n_pacbio = sum(has_pacbio, na.rm = TRUE),
        n_hic = sum(has_hic, na.rm = TRUE),
        n_pacbio_hic = sum(has_pacbio & has_hic, na.rm = TRUE),
        # Biological coherence
        n_pacbio_hic_shared_biosample = sum(has_pacbio & has_hic & shared_pacbio_hic_biosample, na.rm = TRUE),
        # RNA-seq context
        n_rnaseq = sum(has_pacbio & has_hic & has_rnaseq, na.rm = TRUE),
        # Negative space
        n_no_pacbio_no_hic = sum(!has_pacbio & !has_hic, na.rm = TRUE),
        # Remove the grouping for downstream operations 
        .groups = "drop"
    ) %>%
    mutate(
    pct_pacbio = n_pacbio / n_species * 100,
    pct_hic = n_hic / n_species * 100,
    pct_pacbio_hic = n_pacbio_hic / n_species * 100,
    pct_pacbio_hic_shared_biosample =
        ifelse(
            n_pacbio_hic > 0,
            n_pacbio_hic_shared_biosample / n_pacbio_hic * 100,
            NA_real_
        ),
    pct_no_pacbio_no_hic = n_no_pacbio_no_hic / n_species * 100
)


cat("\n== SUMMARY == \n")
print(summary %>% mutate(across(starts_with("pct_"), ~ round(.x, 1))))

# Targeted summaries

# PacBio + Hi-C cases:
pacbio_hic_cases <- inventory %>%
    filter(has_pacbio & has_hic) %>%
    arrange(run_mapping_level, desc(shared_pacbio_hic_biosample), species)

# PacBio + HiC + shared BioSample:
pacbio_hic_shared_biosample <- inventory %>%
    filter(has_pacbio & has_hic & shared_pacbio_hic_biosample) %>%
    arrange(run_mapping_level, species)


# Species without PacBio nor Hi-C, but other strategies:
other_context <- inventory %>% 
    filter(!has_pacbio & !has_hic)%>% 
    separate_rows(other_sequencing_methods, sep = ",")%>% 
    filter(!is.na(other_sequencing_methods), 
    other_sequencing_methods != "") %>% 
    count(other_sequencing_methods, sort = TRUE)

#other_context <- inventory %>%
#    filter(!has_pacbio & !has_hic)%>%
#    arrange(run_mapping_level, species)

cat("\n== Other sequencing strategies ==")

# Exports
write_tsv(summary, "data/inventory_analysis/summary.tsv")
write_tsv(pacbio_hic_cases, "data/inventory_analysis/pacbio_hic.tsv")
write_tsv(pacbio_hic_shared_biosample, "data/inventory_analysis/pacbio_hic_shared_biosample.tsv")
write_tsv(other_context, "data/inventory_analysis/other_sequencing_strategies.tsv")

cat("Wrote: \n",
    "- inventory_analysis/summary.tsv\n",
    "- inventory_analysis/pacbio_hic.tsv\n",
    "- inventory_analysis/other_sequencing_strategies.tsv\n"
    )
