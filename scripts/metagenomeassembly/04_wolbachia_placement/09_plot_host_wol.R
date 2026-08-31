#!/usr/bin/env Rscript

# Draws a host cladogram with branches coloured by host lineages and tips marker by Wolbachia status 
# Usage: 
# module load R-bundle-Bioconductor/3.15-foss-2021a-R-4.2.1
# Rscript 09_plot_host_wol.R

library(ape)
library(ggplot2)
library(ggtree)

HOST_CLAD <- "/data/users/lland/cobionts/hym_species_tree.nwk"
METADATA_FILE <- "/data/projects/p2025-0083_mining_cobionts/sample_sheets/Hymenopteran_genomes_metadata.tsv"
MAG_PLACEMENT <- "/data/users/lland/placement_out/fig2_wolb_tree_grade_mag_supergroup.tsv"
PROCESSED_SPECIES <- "/data/users/lland/processed_species.txt"
OUT <- "/data/users/lland/placement_out/fig2b_host_grade_wolbachia"

# Define the host that carried Wolbachia from two different supergroups
coinfected_hosts <- c("Ametastegia_equiseti", "Seladonia_tumulorum", "Stilbops_ruficornis")
# Host lineages color 
superfamily_grade <- c(Tenthredinoidea = "Symphyta", Cephoidea = "Symphyta", Xiphydrioidea = "Symphyta", Ichneumonoidea  = "Parasitica", Chalcidoidea = "Parasitica", Platygastroidea = "Parasitica", Diaprioidea = "Parasitica", Cynipoidea = "Parasitica", Evanioidea = "Parasitica", Apoidea = "Aculeata", Vespoidea = "Aculeata", Pompiloidea = "Aculeata", Tiphioidea = "Aculeata", Formicoidea = "Aculeata", Chrysidoidea = "Aculeata")
grade_levels <- c("Symphyta", "Parasitica", "Aculeata")
grade_pal <- c(Symphyta = "#c16328", Parasitica = "#a42f94", Aculeata = "#33c16e")  
# Figure configuration
fan_size <- 22 
label_size <- 1.9
point_size <- 2.5 # the status dots 
star_size <- 5.0 # a star for coinfections
branch_width <- 0.4
label_space <- 0.35

# Read host cladogram tree
host_tree <- read.tree(HOST_CLAD)
tips <- host_tree$tip.label
n_tips <- length(tips)
n_all <- n_tips + host_tree$Nnode # tips + internal nodes

# Read metadata
meta<- read.delim(METADATA_FILE, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
family_host <- meta$family[match(tips, meta$species)]
superfamily_host <- meta$superfamily[match(tips, meta$species)]
grade <- unname(superfamily_grade[superfamily_host])

# Define the processed species: 
screened <- trimws(readLines(PROCESSED_SPECIES))
screened <- screened[screened != ""]

# Read placement file, one row for each species
place_mag <- read.delim(MAG_PLACEMENT, stringsAsFactors = FALSE)
names(place_mag)[names(place_mag) == "host_species"] <- "species"
# Make A, B and A+B
supergroup_by_species <- tapply(place_mag$supergroup, place_mag$species, function(x) paste(sort(unique(x)), collapse = "+"))
n_mags_by_species <- table(place_mag$species)

# dataframe with the tip statuses
dat <- data.frame(
  species = tips,
  family = family_host,
  superfamily = superfamily_host,
  grade = grade,
  screened = tips %in% screened,
  supergroups = unname(supergroup_by_species[tips]),
  stringsAsFactors = FALSE)
dat$n_mags <- as.integer(n_mags_by_species[tips])
dat$n_mags[is.na(dat$n_mags)] <- 0L

# Write the different statuses
dat$status <- ifelse(!dat$screened, "not screened", ifelse(is.na(dat$supergroups), "no Wolbachia MAG", dat$supergroups))
dat$coinfected <- dat$species %in% coinfected_hosts
supergroup_levels <- sort(setdiff(unique(dat$status), c("not screened", "no Wolbachia MAG")))
status_levels <- c("not screened", "no Wolbachia MAG", supergroup_levels)
dat$status <- factor(dat$status, levels = status_levels)

# Status color
supergroup_cols <- c("#D55E00", "#dc31d6", "#0072B2")
pal <- c("grey85", "grey55", rep_len(supergroup_cols, length(supergroup_levels)))
names(pal) <- status_levels

# Branch colours: 
tip_grade <- grade
node_grade <- rep(NA_character_, n_all) # vector holding the lineage for every node in the tree
node_grade[seq_len(n_tips)] <- tip_grade # the node IDs 1 to n_tips in ape are the tips, so we can fill them while internal nodes stay NA
# Two column matrix of the tree, sorting by postorder: every edge appears after all edges in the child subtree
E <- reorder(host_tree, "postorder")$edge
for (i in seq_len(nrow(E))) { 
  parent <- E[i, 1] # walk edges
  child <- E[i, 2]
  child_grade <- node_grade[child]
  parent_grade <- node_grade[parent]
  node_grade[parent] <- if (is.na(parent_grade)) child_grade else if (identical(parent_grade, child_grade)) parent_grade else "mixed" # parent still empty takes child grade, parent has already grade and they agree with child so keep it, and if they disagree, mixed
}
node_grade[node_grade == "mixed"] <- NA  
mixed_col <- "grey75"

# Node table for ggtree 
idx <- match(tips, dat$species)
node_df <- data.frame(node = seq_len(n_all), branch_grade = factor(node_grade, levels = grade_levels), stringsAsFactors = FALSE)
node_df$status <- factor(NA, levels = status_levels)
node_df$status[seq_len(n_tips)] <- dat$status[idx]
node_df$coinfected <- NA
node_df$coinfected[seq_len(n_tips)] <- dat$coinfected[idx]
node_df$n_mags <- NA_integer_
node_df$n_mags[seq_len(n_tips)] <- dat$n_mags[idx]

write.table(dat, paste0(OUT, "_tip_status.tsv"),  sep = "\t", quote = FALSE, row.names = FALSE)
# Draw the plot
p <- ggtree(host_tree, layout = "fan", size = branch_width) %<+% node_df +
  aes(colour = branch_grade) +
  geom_tippoint(aes(fill = status), shape = 21, colour = "grey25", stroke = 0.2, size = point_size, na.rm = TRUE) +
  geom_tippoint(aes(subset = !is.na(coinfected) & coinfected), shape = 8, size = star_size, colour = "black", stroke = 0.5, na.rm = TRUE) +
  geom_tiplab(size = label_size, offset = 0.3, colour = "black", na.rm = TRUE) +
  scale_colour_manual(values = grade_pal, name = "Host grade", breaks = grade_levels, na.value = mixed_col, na.translate = FALSE, guide = guide_legend(order = 1)) +
  scale_fill_manual(values = pal, name = "Wolbachia status", na.translate = FALSE, guide = guide_legend(order = 2)) +
  theme(legend.position = "right") +
  hexpand(label_space)

ggsave(paste0(OUT, "_fan.pdf"), p, width = fan_size, height = fan_size, limitsize = FALSE)