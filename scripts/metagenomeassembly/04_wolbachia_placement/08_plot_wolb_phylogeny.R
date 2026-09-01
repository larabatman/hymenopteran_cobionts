#!/usr/bin/env Rscript

# Plot the Wolbachia phylogram as a circular fan where tips are filled by supergroup colours, with the reference strains in circle and the recovered MAGs as triangles
# The MAG text is also coloured by host lineages Symphyta / Parasitica / Aculeata
# Computes the distance matrix in patristic distances to place the Wolbachia MAGs into their nearest reference supergroup
# Usage:
# Rscript 08_plot_wolb_phylogeny.R

library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(readr)

# .treefile from the IQTREE inferrence
TREE_FILE <- "placement_out/wolb_mag_tree_clean.treefile"
# File with the list of reference strians and their corresponding supergroups
ACC2SUPERGROUP <- "/data/projects/p2025-0083_mining_cobionts/nuwt_scan/wolbachia_database/accession_to_supergroup.tsv"
# Metadata file for the host lineages
METADATA_FILE <- "/data/projects/p2025-0083_mining_cobionts/sample_sheets/Hymenopteran_genomes_metadata.tsv"
# Prefix for the figure
OUT <- "placement_out/fig2_wolb_tree_grade"

ASSEMBLER_TAGS <- "metamdbg|myloasm|metaflye"
OUTGROUP_TAGS <- c("outgroup")

# Colouring by host lineages
superfamily_grade <- c(Tenthredinoidea = "Symphyta", Cephoidea = "Symphyta", Xiphydrioidea = "Symphyta",Ichneumonoidea  = "Parasitica", Chalcidoidea = "Parasitica", Platygastroidea = "Parasitica", Diaprioidea = "Parasitica", Cynipoidea = "Parasitica", Evanioidea = "Parasitica", Apoidea = "Aculeata", Vespoidea = "Aculeata", Pompiloidea = "Aculeata", Tiphioidea = "Aculeata", Formicoidea = "Aculeata", Chrysidoidea = "Aculeata")
grade_levels <- c("Symphyta", "Parasitica", "Aculeata")
grade_pal <- c(Symphyta   = "#c16328", Parasitica = "#a42f94", Aculeata   = "#33c16e")

# Figure configuration
fig_size <- 30 # size of the canvas
branch_width <- 0.35 # thickness of the branches
open_angle <- 5 # the fan opening
cladogram <- FALSE # phylogram
mag_label_size <- 2.5 # size of the mags, they will be labelled with a W prefix
ref_label_size <- 1.5
mag_point_size <- 3.5 # triangle size
ref_point_size <- 3.0 # circle size

# Read the tree and the accession to supergroups map
phylogram <- read.tree(TREE_FILE)
map <- read_tsv(ACC2SUPERGROUP, col_names = c("tip", "supergroup"), show_col_types = FALSE)
ref_supergroup <- setNames(map$supergroup, map$tip)
# Define the outgroups and root the tree
ref_outgroup <- intersect(map$tip[map$supergroup %in% OUTGROUP_TAGS], phylogram$tip.label)
ref_supergroup[ref_outgroup] <- "OG"
phylogram <- root(phylogram, outgroup = ref_outgroup, resolve.root = TRUE)
# Define the reference strains 
is_ref <- phylogram$tip.label %in% map$tip
# Define the mags: they are not in the accession to supergroup map
mags <- phylogram$tip.label[!is_ref]
refs <- setdiff(phylogram$tip.label[is_ref], ref_outgroup)

# Each MAG takes the supergroup of its nearest reference
# cophenetic() gives patristic distances which are the branch lengths summed along the path between two tips.
# D <- cophenetic(phylogram) is a 255x255 table of distances between every pair of tips, where rows and columns are both tip names
#                     Ace     GCA_000167475   GCA_008033215 ...
# Ace                 0.0     1.523           1.521
# GCA_000167475       1.523   0.0             0.009
# Acromyrmex_..._1    1.498   0.021           0.019
# Each number is the distance along the tree
# D[m, refs] for instance for m = "Acromyrmex_octospinosus_metamdbg_magscot_1" it takes that one row and keeps only the 156 reference columns
# The named vector:
# GCA_000167475 GCA_008033215 GCA_902648465 ... 156 entries
# 0.021         0.019         0.034
# Do which.min(...) to return the position of the smallest value, and names(...) converts that position to a name and the reference this MAG sits closest on the tree
# vapply to repeat that for all the MAGs
# Resulting in: 
# mag                       nearest_ref
# Acromyrmex_..._1          GCA_008033215
# Andrena_barbilaris_..._1  GCA_000167475
# Then look up which supergroups the references were from 
# we have ref_supergroup that comes from the accession_supergroup.tsv: ref_supergroup["GCA_008033215"] returns "A"
D <- cophenetic(phylogram)
place <- tibble(
  mag = mags,
  nearest_ref = vapply(mags, function(m) names(which.min(D[m, refs])), character(1)))
place$supergroup <- unname(ref_supergroup[place$nearest_ref])

# Trim the MAG ids
place$host_species <- sub(paste0("_(", ASSEMBLER_TAGS, ")_.*$"), "", place$mag) # Acromyrmex_octospinosus_metamdbg_magscot_1 becomes Acromyrmex_octospinosus
# Read the metadata table 
meta <- read.delim(METADATA_FILE, stringsAsFactors = FALSE, na.strings = c("NA", "NaN", ""))
# To make sure that everything matches the place tibble, use match(x, table) that returns, for each eleemnt of x, its position in table or NA if absent
# if Acromyrmex_octospinosus is row 47 in the metadata, match gives 47 and its superfamily column gives Formicoidea
place$superfamily <- meta$superfamily[match(place$host_species, meta$species)]
# Now match that superfamily to its larger lineage for plotting 
place$grade <- unname(superfamily_grade[place$superfamily]) # unname takes out the top row of the named vector
write_tsv(place, paste0(OUT, "_mag_supergroup.tsv"))

# Build the tip annotation table
tip <- phylogram$tip.label # tip names in the phylogram order
kind <- ifelse(tip %in% mags, "MAG", ifelse(tip %in% ref_outgroup, "OG", "ref")) # three categories: MAG, OG if it is an outgroup, and ref if it is anything else
# If it is a mag, return its row number in place else NA
mag_supergroup <- place$supergroup[match(tip, place$mag)]
supergroup <- ifelse(kind == "MAG", mag_supergroup, unname(ref_supergroup[tip])) # for MAG tips, take their supergroups and for everything else, look up the supergroup in ref_supergroup
# Text colour for MAGs and labels: the host lineages
grade <- factor(place$grade[match(tip, place$mag)], levels = grade_levels)
# MAG labels: create two columns, lab_mag and lab_ref. if MAG, paste W and its host, otherwise NA
host <- sub(paste0("_(", ASSEMBLER_TAGS, ")_.*$"), "", tip)
lab_mag <- ifelse(kind == "MAG", paste0("W", host), NA_character_)
lab_ref <- ifelse(kind != "MAG", tip, NA_character_)
# Tibble for the plot
tipdat <- tibble(tip, kind, supergroup, grade, lab_mag, lab_ref)

# Color the tips for each supergroup
# Supergroups in the tree, without the outgroups:
supergroups <- sort(unique(na.omit(tipdat$supergroup[tipdat$supergroup != "OG"])))
# Supergroups A and B get specific colours to distinguish them
supA_supB <- c(A = "#D55E00", B = "#0072B2")
supA_supB <- supA_supB[names(supA_supB) %in% supergroups]
rest_supergroups <- setdiff(supergroups, names(supA_supB))
supergroup_palette <- c(supA_supB, setNames(hcl.colors(length(rest_supergroups), "Dark 3"), rest_supergroups))
supergroup_palette["OG"] <- "grey60"
# Figure
p <- ggtree(phylogram, layout = "fan", open.angle = open_angle, size = branch_width, branch.length = "branch.length") %<+% tipdat +
  geom_tippoint(aes(fill = supergroup, shape = kind, size = kind), colour = "grey20", stroke = 0.25, alpha = 0.9) +   # one point layer: fill = supergroup, shape & size = kind, so a single MAG/reference shape legend appears alongside the supergroup legend.
  geom_tiplab(aes(label = lab_mag, colour = grade), size = mag_label_size, offset = 0.012, fontface = "bold", na.rm = TRUE) +
  geom_tiplab(aes(label = lab_ref), size = ref_label_size, offset = 0.008, colour = "grey35", na.rm = TRUE) + # reference accession labels: small, grey
  scale_fill_manual(values = supergroup_palette, na.value = "grey70", name = "Wolbachia supergroup", na.translate = FALSE, guide = guide_legend(order = 1, override.aes = list(shape = 21, size = 3.5))) +
  scale_shape_manual(values = c(MAG = 24, ref = 21, OG = 21), breaks = c("MAG", "ref"),labels = c(MAG = "recovered MAG", ref = "reference strain"), name = NULL,guide = guide_legend(order = 2, override.aes = list(fill = "grey50", size = 3.5))) +
  scale_size_manual(values = c(MAG = mag_point_size, ref = ref_point_size, OG = ref_point_size), guide = "none") +
  scale_colour_manual(values = grade_pal, na.value = "grey50", guide = "none") +
  ggtitle("Wolbachia MAG phylogeny") +
  theme(legend.position = "right")

ggsave(paste0(OUT, ".pdf"), p, width = fig_size, height = fig_size, limitsize = FALSE)