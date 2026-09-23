#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: 09_plot_tree.R TREEFILE METADATA.csv OUTPUT.pdf")

suppressPackageStartupMessages(library(ape))
tree <- read.tree(args[1])
metadata <- read.csv(args[2], check.names = FALSE, fileEncoding = "UTF-8-BOM")
metadata <- metadata[match(tree$tip.label, metadata$sampleID), ]
if (any(is.na(metadata$sampleID))) stop("Some tree tips are absent from metadata")

outgroups <- metadata$sampleID[is.na(metadata$s380) | metadata$s380 == ""]
outgroups <- intersect(outgroups, tree$tip.label)
if (length(outgroups)) {
  tree <- tryCatch(root(tree, outgroup = outgroups, resolve.root = TRUE), error = function(e) tree)
  metadata <- metadata[match(tree$tip.label, metadata$sampleID), ]
}

groups <- sort(unique(as.character(metadata$spID)))
colors <- setNames(hcl.colors(length(groups), "Dark 3"), groups)
tip_colors <- colors[as.character(metadata$spID)]

pdf(args[3], width = 14, height = 24, useDingbats = FALSE)
par(mar = c(1, 1, 2, 1), xpd = NA)
plot(
  tree, type = "phylogram", direction = "rightwards",
  show.tip.label = TRUE, tip.color = tip_colors,
  cex = 0.28, label.offset = 0.0001,
  main = "Pan-Dryas maximum-likelihood tree"
)
add.scale.bar(cex = 0.7)
legend("topleft", legend = groups, col = colors, pch = 19, cex = 0.65, bty = "n", ncol = 2)
dev.off()

