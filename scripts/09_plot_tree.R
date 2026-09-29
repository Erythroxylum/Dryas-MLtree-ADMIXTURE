#!/usr/bin/env Rscript

# to run:
#TIP_LABEL_COLUMN=phyloID \
#S170_COLUMN=s170_BPP \
#S47_COLUMN=s47-p9 \
#TIP_LABEL_CEX=0.35 \
#NODE_LABEL_CEX=0.4 \
#MEMBERSHIP_DOT_CEX=0.65 \
#LEGEND_CEX=0.65 \
#MIN_SUPPORT=70 \
#Rscript scripts/09_plot_tree.R \
#  results/ml_tree/dryas_cleaned_asc.treefile \
#  data/Dryas_sampledata.csv \
#  figures/09_dryas_tree.pdf

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: 09_plot_tree.R TREEFILE METADATA.csv OUTPUT.pdf")

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phytools))
tree <- read.tree(args[1])
metadata <- read.csv(args[2], check.names = FALSE, fileEncoding = "UTF-8-BOM")

tip_label_cex <- as.numeric(Sys.getenv("TIP_LABEL_CEX", "0.40"))
node_label_cex <- as.numeric(Sys.getenv("NODE_LABEL_CEX", "0.75"))
membership_dot_cex <- as.numeric(Sys.getenv("MEMBERSHIP_DOT_CEX", "0.65"))
legend_cex <- as.numeric(Sys.getenv("LEGEND_CEX", "0.65"))
minimum_support <- as.numeric(Sys.getenv("MIN_SUPPORT", "50"))
tip_label_column <- Sys.getenv("TIP_LABEL_COLUMN", "phyloID")
s170_column <- Sys.getenv("S170_COLUMN", "s170_BPP")
s47_column <- Sys.getenv("S47_COLUMN", "s47-p9")

is_member <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  z <- tolower(trimws(as.character(x)))
  !is.na(x) & !(z %in% c("", "false", "f", "0", "na", "nan", "none"))
}
extract_support <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(NA_real_)
  z <- suppressWarnings(as.numeric(strsplit(x, "/", fixed = TRUE)[[1]]))
  z <- z[is.finite(z)]
  if (length(z)) tail(z, 1) else NA_real_
}

for (z in c("sampleID", "s380", "spID", s170_column, s47_column)) {
  if (!(z %in% names(metadata))) stop("Metadata column not found: ", z)
}
if (!(tip_label_column %in% names(metadata))) tip_label_column <- "sampleID"
if (anyDuplicated(tree$tip.label)) stop("Tree has duplicated tip labels")
if (is.null(tree$edge.length)) stop("Tree lacks branch lengths")

outgroups <- metadata$sampleID[is.na(metadata$s380) | metadata$s380 == ""]
outgroups <- intersect(outgroups, tree$tip.label)
if (length(outgroups) != 4) stop("Expected four outgroup tips; found ", length(outgroups))
if (is.rooted(tree)) tree <- unroot(tree)
ingroup <- setdiff(tree$tip.label, outgroups)
descendant_tips <- function(phy, node) {
  if (node <= Ntip(phy)) phy$tip.label[node] else extract.clade(phy, node)$tip.label
}
edges <- which(vapply(seq_len(nrow(tree$edge)), function(i) {
  d <- descendant_tips(tree, tree$edge[i, 2])
  setequal(d, outgroups) || setequal(d, ingroup)
}, logical(1)))
if (length(edges) != 1) stop("Could not find one edge separating the four outgroups")
e <- edges[1]
tree <- phytools::reroot(tree, tree$edge[e, 2], tree$edge.length[e] / 2)
tree <- ladderize(tree, right = TRUE)

metadata <- metadata[match(tree$tip.label, metadata$sampleID), , drop = FALSE]
if (any(is.na(metadata$sampleID))) stop("Some tree tips are absent from metadata")

# Recover the exact bottom-to-top plotted order and export both full-tree and
# ADMIXTURE (s380) orders for scripts 10--12.
pdf(NULL)
plot(tree, show.tip.label = FALSE, plot = FALSE)
probe <- get("last_plot.phylo", envir = .PlotPhyloEnv)
dev.off()
plot_order <- tree$tip.label[order(probe$yy[seq_len(Ntip(tree))])]
prefix <- sub("\\.pdf$", "", args[3], ignore.case = TRUE)
writeLines(plot_order, paste0(prefix, "_tree_tip_order.txt"))
admix_ids <- metadata$sampleID[!is.na(metadata$s380) & metadata$s380 != ""]
writeLines(plot_order[plot_order %in% admix_ids],
           paste0(prefix, "_admixture_sample_order.txt"))

groups <- sort(unique(as.character(metadata$spID)))
colors <- setNames(hcl.colors(length(groups), "Dark 3"), groups)
tip_colors <- colors[as.character(metadata$spID)]
display_labels <- as.character(metadata[[tip_label_column]])
display_labels[is.na(display_labels) | display_labels == ""] <-
  metadata$sampleID[is.na(display_labels) | display_labels == ""]
display_tree <- tree
display_tree$tip.label <- display_labels
s170 <- is_member(metadata[[s170_column]])
s47 <- is_member(metadata[[s47_column]])
max_depth <- max(node.depth.edgelength(tree))

pdf(args[3], width = 10, height = 24, useDingbats = FALSE)
par(mar = c(1, 1, 2, 1), xpd = NA)
plot(
  display_tree, type = "phylogram", direction = "rightwards",
  show.tip.label = TRUE, tip.color = tip_colors,
  cex = tip_label_cex, label.offset = max_depth * 0.006,
  align.tip.label = TRUE, x.lim = c(0, max_depth * 1.35),
  main = "Pan-Dryas maximum-likelihood tree"
)
add.scale.bar(cex = 0.7)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- lp$yy[seq_len(Ntip(tree))]
if (!is.null(tree$node.label)) {
  support <- vapply(tree$node.label, extract_support, numeric(1))
  keep <- which(is.finite(support) & support >= minimum_support)
  if (length(keep)) nodelabels(round(support[keep]), Ntip(tree) + keep,
                                frame = "none", cex = node_label_cex,
                                adj = c(1.05, -0.15))
}
track_x <- c(max_depth * 1.19, max_depth * 1.23)
points(rep(track_x[1], sum(s170)), tip_y[s170], pch = 16,
       cex = membership_dot_cex, col = "#2166ac")
points(rep(track_x[2], sum(s47)), tip_y[s47], pch = 16,
       cex = membership_dot_cex, col = "#b2182b")
text(track_x, max(tip_y) + 2, c("s170", "s47"), cex = 0.70, font = 2)
legend("left", legend = groups, col = colors, pch = 19,
       cex = legend_cex, bty = "n", ncol = 2, title = "Tip groups")
#legend("topright", legend = c("s170", "s47"), col = c("#2166ac", "#b2182b"),
#       pch = 16, pt.cex = membership_dot_cex, cex = legend_cex, bty = "n",
#       title = "Analysis membership")
dev.off()

message("Tree: ", args[3])
message("Tree order: ", paste0("09_", prefix, "_tree_tip_order.txt"))
message("ADMIXTURE order: ", paste0("09_", prefix, "_admixture_sample_order.txt"))
