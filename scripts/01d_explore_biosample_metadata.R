#!/usr/bin/env Rscript
# This script performs explotaroty analysis of BioSample annotations. Thus results are descriptive and do not necessarily reflect assembly inputs. 

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# Paths
input_file <- "data/annotation/data_inventory_biosample_annotation.tsv"
out_dir <- "data/annotation/exploration"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Data loading
annot <- read_tsv(input_file, show_col_types = FALSE)

# Normalize raw fields 
annot <- annot %>%
    mutate(
        sex = str_to_lower(str_trim(sex)),
        tissue = str_to_lower(str_trim(tissue)),
        dev_stage = str_to_lower(str_trim(dev_stage))
    )

# Explicit classification
classified <- annot %>%
    mutate(
        sex_class = case_when(
            sex == "female" ~ "female",
            sex == "male" ~ "male",
            TRUE ~ "unknown_sex"
        ),
        tissue_class = case_when(
            str_detect(tissue, "whole") ~ "whole_body",
            str_detect(tissue, "abdomen") ~ "abdomen",
            str_detect(tissue, "mid body") ~ "mid_body",
            str_detect(tissue, "head|thorax") ~"head_thorax",
            TRUE ~ "unknown_tissue"
        ),
        stage_class = case_when(
            dev_stage == "adult" ~ "adult",
            dev_stage == "larva" ~ "larva", 
            dev_stage == "pupa" ~ "pupa", 
            TRUE ~ "unknown_stage"
        )
    )

# Glboal combination summary 
combo_summary <- classified %>%
    count(sex_class, tissue_class, stage_class, name = "n_samples") %>%
    arrange(desc(n_samples))

write_tsv(combo_summary, file.path(out_dir, "combination_summary.tsv"))

# One file per combination
classified %>%
    group_by(sex_class, tissue_class, stage_class) %>%
    group_walk(~{
        fname <- paste0(
            .y$sex_class, "__",
            .y$tissue_class, "__", 
            .y$stage_class,
            ".tsv"
        )
        write_tsv(.x, file.path(out_dir, fname))
    })

# Dimension-wise summaries
write_tsv(
    classified %>% count(sex_class),
    file.path(out_dir, "sex_summary.tsv")
)

write_tsv(
    classified %>% count(tissue_class),
    file.path(out_dir, "tissue_summary.tsv")
)

write_tsv(
    classified %>% count(stage_class),
    file.path(out_dir, "stage_summary.tsv")
)

cat("Exploration written to:", out_dir, "\n")