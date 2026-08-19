#!/usr/bin/env Rscript

# Annotation and filtering step between the placement and figures. This tajes the per-fragment supergroup calls and attaches their annotations to them
# The annotations are there further analyzed as mito OG, ankyrinn repeats, and residual mitochondrial OG 
# The tow outputs, with and without the mitochomdrial OGs, are kept for completenedd
#
# residual_mito.txt lists orthogroups with the same characteristics but could not be 100% sure
# removed_TEs.txt lists orthogroups excluded before placement, that should not be here
#
# Input:
# _per_nuwt.tsv from aggregate_supergroups.R: one row per assigned fragment as name, species, supergroup, og, length
# og_lookup_table.tsv: og_class is the four-valued SwissProt classification, sp_description the product name.
# files dir with removed_mito.txt, residual_mito.txt and removed_TEs.txt, one OG id per line.
#
# Output
# <out_prefix>_annotated.tsv every assigned fragment, plus og_class, product, og_type and excluded. Feeds the master table and the supergroup composition figure.
# <out_prefix>_kept.tsv the same, minus excluded fragments. Feeds merge_nuwt_regions.R and therefore Figure 2.
#
# USAGE
# Rscript annotate_nuwt_fragments.R <per_nuwt_tsv> <og_lookup_tsv> <files_dir> <out_prefix>

args <- commandArgs(trailingOnly = TRUE)

per_nuwt <- args[1]
lookup_f <- args[2]
files_dir <- args[3]
out_prefix <- args[4]

read_tsv0 <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

# Read they excluded/ weird OG text files
# one OG id per line but still trimmed with blanks dropped, so a missing trailing newline or stray whitespace cannot silently break the %in% tests below
read_ids <- function(f) {
  x <- trimws(readLines(file.path(files_dir, f), warn = FALSE))
  x[nzchar(x)]
}
removed_mito  <- read_ids("removed_mito.txt")
residual_mito <- read_ids("residual_mito.txt")
removed_TEs   <- read_ids("removed_TEs.txt")

# Attach the orthogroup annotation
d <- read_tsv0(per_nuwt)
lk <- read_tsv0(lookup_f)[, c("og_id", "og_class", "sp_description")]
names(lk) <- c("og", "og_class", "product")

d <- merge(d, lk, by = "og", all.x = TRUE) # all.x to keep every fragment
d$og_class[is.na(d$og_class)] <- "no_swissprot_hit"
d$product[is.na(d$product)] <- "no_hit"

# Classify the orthogroups: First match wins. 
# The curated lists takes priority over an automatic annotation, and no annotation are above unnanotated ones so they are keep as such rather than other wolbachia genes
UNINFORMATIVE <- "^no_hit$|hypothetical|uncharacteri[sz]ed|unknown function|^$"
has_product <- !grepl(UNINFORMATIVE, d$product, ignore.case = TRUE)
# First match wins rule
d$og_type <- ifelse(
  d$og %in% removed_mito, "mitochondrial (removed)",
  ifelse(d$og %in% residual_mito, "mitochondrial (residual)",
  ifelse(grepl("ankyrin", d$product, ignore.case = TRUE), "ankyrin", # check for ankyrin in the description
  ifelse(d$og_class == "bacterial_IS", "bacterial IS",
  ifelse(d$og_class == "no_swissprot_hit" & !has_product,"no annotation", "other Wolbachia gene")))))

# The removed mitochondrial OGs
d$excluded <- d$og %in% removed_mito

# Write tables before and after removal
cols <- c("name", "species", "supergroup", "og", "length", "og_class", "product", "og_type", "excluded")
d <- d[order(d$species, d$supergroup, d$og), cols]

write.table(d, paste0(out_prefix, "_annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(d[!d$excluded, ], paste0(out_prefix, "_kept.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)