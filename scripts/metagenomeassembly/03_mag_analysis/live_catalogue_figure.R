#!/usr/bin/env Rscript

# Make plots for live catalogue
# This input is the MAGs_complete_list.tsv that has been manually curated to retain bins that were > 90% complete and <5% contaminated
# The file contains one row for each MAG
# Usage:
# Rscript live_catalogue_figure.R  mag_file outdir

args <- commandArgs(trailingOnly = TRUE)
MAG_FILE <- args[1]
OUTDIR <- args[2]

mag_list <- read.delim(MAG_FILE, stringsAsFactors = FALSE)

# Take the genus names: retain Genus column
# Two Genus names are relabelled and the MAGs with no genus level assignment were taken to their highest resolved ranks
# Order alphabetically 
ORDER <- "alpha"
genus <- mag_list$Cobiont.genus
genus[genus == "G964019805"] <- "Ca. Lariskella" # taken from NCBI classification
genus[genus == "Blochmanniella"] <- "Ca. Blochmanniella"
mag_list$genus <- genus

# Counts for each genus
# table counts rows for each genus and returns a table in 1D with c(genus). Since each row is one MAG, this counts the MAGs in each genus name 
n_mags <- c(table(mag_list$genus))

# plotting order
# barplot(horiz=TRUE) draws the first element at the bottom of the panel and works upward so if we want order, we need to give it backwards 
ord <- if (ORDER == "count") order(n_mags, decreasing = FALSE) else order(names(n_mags), decreasing = TRUE) # From Z-A, from smallest to largest
n_mags <- n_mags[ord]

# Catalogue: one bar for each genus
pdf(file.path(OUTDIR, "Figure_catalogue.pdf"), width = 6.5, height = 7.5)
par(mar = c(4, 11, 1, 3))
bar <- barplot(n_mags, horiz = TRUE, las = 1, names.arg = "", xlim = c(0, max(n_mags) * 1.1), col = "#4a7ebb", border = NA, xlab = "MAGs", space = 0.35)

# genus names italic on the left, counts printed at the end of each bar
axis(2, at = bar, labels = names(n_mags), las = 1, tick = FALSE, font = 3, cex.axis = 0.75, line = -0.5)
text(n_mags, bar, labels = n_mags, pos = 4, cex = 0.65, offset = 0.3)
dev.off()

# Co-occurrence and multiple MAGs
# Hosts carrying more than one MAG, as a heatmap
count_host <- table(mag_list$Host.species) # MAGs in each host species
multi_mags <- names(count_host)[count_host > 1] # select the host names with more than one MAG
genus_multi_mags <- mag_list[mag_list$Host.species %in% multi_mags, ] # go back in dataframe and keep every row that has multi to have the genus information as well
mat <- table(genus_multi_mags$Host.species, genus_multi_mags$genus) # rows = hosts, cols = genera
mat <- as.matrix(mat)

# choose ordering: most cobionts first then alphabetical
mat <- mat[order(rowSums(mat), rownames(mat), decreasing = c(TRUE, FALSE), method = "radix"), , drop = FALSE]
first_host <- apply(mat > 0, 2, function(x) which(x)[1])
mat <- mat[, order(-colSums(mat), first_host), drop = FALSE]
mat <- mat[nrow(mat):1, , drop = FALSE]

# Colour palette
col_pal <- c("grey95", "#c6dbef", "#6baed6", "#08519c")
breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)

pdf(file.path(OUTDIR, "Figure_cooccurrence.pdf"), width = 8, height = 7)
par(mar = c(9, 12, 2, 2))
# Heatmap that gives a matrix of numbers and paints one coloured rectangle per cell, 
# transpose: image(x, y, z) reads z[i, j] as position i along x, position j along y such that z row = horizontal axis and z columns = vertical one and since we want genera along the bottom and hosts up, t(m) because m is hosts x genera
image(x = 1:ncol(mat), y = 1:nrow(mat), z = t(mat), col = col_pal, breaks = breaks, axes = FALSE, xlab = "", ylab = "")
# thin white grid so the tiles read as separate cells
abline(v = seq(0.5, ncol(mat) + 0.5), col = "white", lwd = 1.5)
abline(h = seq(0.5, nrow(mat) + 0.5), col = "white", lwd = 1.5)

axis(1, at = 1:ncol(mat), labels = colnames(mat), las = 2, tick = FALSE, font = 3, cex.axis = 0.7)
axis(2, at = 1:nrow(mat), labels = rownames(mat), las = 1, tick = FALSE,  font = 3, cex.axis = 0.7)

legend("topright", inset = c(0, -0.06), xpd = TRUE, horiz = TRUE, bty = "n", cex = 0.7, legend = c("1", "2", "3 MAGs"), fill = col_pal[2:4], border = NA)
dev.off()