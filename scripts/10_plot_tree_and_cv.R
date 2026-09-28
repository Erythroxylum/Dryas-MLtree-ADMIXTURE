#!/usr/bin/env Rscript

# Produce the rooted Pan-Dryas ML tree, full ADMIXTURE cross-validation
# summary, and active-tree sample order used by the ancestry plotting script.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage: 10_plot_tree_admixture.R TREEFILE ADMIXTURE_DIR",
    "METADATA.csv OUTPUT_PREFIX_OR_PDF [K_MIN] [K_MAX]"
  ))
}

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phytools))

tree_file <- args[1]
admixture_dir <- args[2]
metadata_file <- args[3]
output_argument <- args[4]
k_min <- if (length(args) >= 5) as.integer(args[5]) else -Inf
k_max <- if (length(args) >= 6) as.integer(args[6]) else Inf
s170_column <- Sys.getenv("S170_COLUMN", "s170_BPP")
s47_column <- Sys.getenv("S47_COLUMN", "s47-p9")
tip_label_column <- Sys.getenv("TIP_LABEL_COLUMN", "phyloID")
minimum_support <- as.numeric(Sys.getenv("MIN_SUPPORT", "70"))

# Retain compatibility with the previous command, which supplied a .pdf name.
output_prefix <- sub("\\.pdf$", "", output_argument, ignore.case = TRUE)
tree_output <- paste0(output_prefix, "_tree.pdf")
cv_output <- paste0(output_prefix, "_cv.pdf")
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

read_metadata <- function(path) {
  read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM",
           stringsAsFactors = FALSE)
}

is_member <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  value <- tolower(trimws(as.character(x)))
  !is.na(x) & !(value %in% c("", "false", "f", "0", "na", "nan", "none"))
}

extract_support <- function(label) {
  if (is.na(label) || !nzchar(trimws(as.character(label)))) return(NA_real_)
  pieces <- strsplit(as.character(label), "/", fixed = TRUE)[[1]]
  values <- suppressWarnings(as.numeric(pieces))
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  tail(values, 1)
}

metadata <- read_metadata(metadata_file)
required_metadata <- c("sampleID", "s384", "s380", "spID")
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) {
  stop("Metadata columns missing: ", paste(missing_metadata, collapse = ", "))
}
if (!(s170_column %in% names(metadata))) {
  stop("s170 membership column not found: ", s170_column)
}
if (!(s47_column %in% names(metadata))) {
  stop("s47 membership column not found: ", s47_column)
}
if (!(tip_label_column %in% names(metadata))) tip_label_column <- "sampleID"

tree <- read.tree(tree_file)
if (is.null(tree) || !length(tree$tip.label)) stop("Tree could not be read")
if (anyDuplicated(tree$tip.label)) stop("Tree contains duplicated tip labels")

# Remove the arbitrary Newick root, locate the edge whose bipartition separates
# all four non-Dryas samples from the Dryas ingroup, and place the root at the
# midpoint of that edge. This does not assume which outgroup genus is closest
# to Dryas and does not alter relationships within either side of the split.
outgroups <- metadata$sampleID[
  !is.na(metadata$s384) & metadata$s384 != "" &
    (is.na(metadata$s380) | metadata$s380 == "")
]
outgroups <- intersect(outgroups, tree$tip.label)
if (length(outgroups) != 4) {
  stop("Expected four non-Dryas outgroups but found: ",
       paste(outgroups, collapse = ", "))
}
if (is.null(tree$edge.length)) {
  stop("Tree has no branch lengths; cannot calculate the midpoint root")
}
if (is.rooted(tree)) tree <- unroot(tree)
ingroup <- setdiff(tree$tip.label, outgroups)

descendant_tips <- function(phy, node) {
  if (node <= Ntip(phy)) return(phy$tip.label[node])
  extract.clade(phy, node = node)$tip.label
}

separating_edges <- which(vapply(seq_len(nrow(tree$edge)), function(i) {
  descendants <- descendant_tips(tree, tree$edge[i, 2])
  setequal(descendants, outgroups) || setequal(descendants, ingroup)
}, logical(1)))

if (length(separating_edges) != 1) {
  stop(
    "Could not identify one branch separating the four outgroups from Dryas; ",
    "candidate edge count = ", length(separating_edges),
    ". Inspect the inferred outgroup topology before rooting."
  )
}
root_edge <- separating_edges[1]
root_node <- tree$edge[root_edge, 2]
root_position <- tree$edge.length[root_edge] / 2
tree <- phytools::reroot(
  tree, node.number = root_node, position = root_position
)
tree <- ladderize(tree, right = TRUE)
tip_ids <- tree$tip.label

# ladderize() changes the edge traversal used for plotting but does not simply
# rewrite tree$tip.label into vertical plot order. Probe the exact coordinates
# that ape will use and explicitly recover the bottom-to-top order of the tips.
grDevices::pdf(NULL)
plot(
  tree, type = "phylogram", direction = "rightwards",
  show.tip.label = TRUE, align.tip.label = TRUE, plot = FALSE
)
tree_order_probe <- get("last_plot.phylo", envir = .PlotPhyloEnv)
grDevices::dev.off()
plotted_tip_ids <- tip_ids[
  order(tree_order_probe$yy[seq_len(Ntip(tree))])
]
writeLines(plotted_tip_ids, paste0(output_prefix, "_tree_tip_order.txt"))

tree_metadata <- metadata[match(tip_ids, metadata$sampleID), , drop = FALSE]
if (any(is.na(tree_metadata$sampleID))) {
  missing <- tip_ids[is.na(tree_metadata$sampleID)]
  stop("Tree tips absent from metadata: ", paste(missing, collapse = ", "))
}

display_labels <- tree_metadata[[tip_label_column]]
display_labels[is.na(display_labels) | display_labels == ""] <-
  tip_ids[is.na(display_labels) | display_labels == ""]

groups <- sort(unique(as.character(tree_metadata$spID)))
tip_palette <- setNames(hcl.colors(length(groups), "Dark 3"), groups)
tip_colors <- unname(tip_palette[as.character(tree_metadata$spID)])
s170_membership <- is_member(tree_metadata[[s170_column]])
s47_membership <- is_member(tree_metadata[[s47_column]])

sample_order <- readLines(file.path(admixture_dir, "sample_order.txt"))
if (anyDuplicated(sample_order)) stop("sample_order.txt contains duplicated IDs")
cv <- read.delim(file.path(admixture_dir, "cross_validation.tsv"),
                 stringsAsFactors = FALSE)
cv <- cv[cv$K >= k_min & cv$K <= k_max, , drop = FALSE]
if (!nrow(cv)) stop("No ADMIXTURE K values remain after filtering")
if (any(!is.finite(cv$CV_error))) stop("Non-numeric or missing CV errors found")

# Validate that all requested K values contain the same replicate set.
runs <- cv[order(cv$K, cv$replicate), c("K", "replicate", "CV_error"),
           drop = FALSE]
if (anyDuplicated(runs[c("K", "replicate")])) {
  stop("cross_validation.tsv has duplicated K/replicate combinations")
}
replicates_by_k <- split(runs$replicate, runs$K)
all_replicates <- sort(unique(runs$replicate))
bad_k <- names(replicates_by_k)[!vapply(
  replicates_by_k, function(x) setequal(x, all_replicates), logical(1)
)]
if (length(bad_k)) {
  stop(
    "The same replicate set was not found for every K. Incomplete K = ",
    paste(bad_k, collapse = ", ")
  )
}

# Save the active-tree order for the standalone ancestry plotting script.
ingroup_tip_ids <- plotted_tip_ids[plotted_tip_ids %in% sample_order]
if (length(ingroup_tip_ids) != length(sample_order)) {
  missing <- setdiff(sample_order, tip_ids)
  stop("ADMIXTURE samples absent from tree: ", paste(missing, collapse = ", "))
}
writeLines(ingroup_tip_ids, paste0(output_prefix, "_admixture_sample_order.txt"))

# Figure 1: rooted ML tree with analysis-membership annotations.
display_tree <- tree
display_tree$tip.label <- display_labels
max_depth <- max(node.depth.edgelength(tree))

pdf(tree_output, width = 13, height = 30, useDingbats = FALSE)
par(mar = c(1.5, 1, 2.2, 1), xpd = NA)
plot(
  display_tree, type = "phylogram", direction = "rightwards",
  show.tip.label = TRUE, tip.color = tip_colors,
  cex = 0.25, align.tip.label = TRUE,
  label.offset = max_depth * 0.006,
  x.lim = c(0, max_depth * 1.45), no.margin = FALSE
)
title("Rooted and ladderized maximum-likelihood tree", cex.main = 1.0)
add.scale.bar(cex = 0.6, lwd = 0.8)
last_tree_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- last_tree_plot$yy[seq_len(Ntip(tree))]

if (!is.null(tree$node.label)) {
  support <- vapply(tree$node.label, extract_support, numeric(1))
  show_support <- which(is.finite(support) & support >= minimum_support)
  if (length(show_support)) {
    nodelabels(
      text = round(support[show_support]),
      node = Ntip(tree) + show_support,
      frame = "none", cex = 0.35, adj = c(1.05, -0.15)
    )
  }
}

track_x <- c(max_depth * 1.34, max_depth * 1.39)
points(rep(track_x[1], sum(s170_membership)), tip_y[s170_membership],
       pch = 15, cex = 0.35, col = "#2166ac")
points(rep(track_x[2], sum(s47_membership)), tip_y[s47_membership],
       pch = 15, cex = 0.35, col = "#b2182b")
text(track_x, max(tip_y) + 2, labels = c("s170", "s47"),
     cex = 0.65, font = 2)
legend(
  "topleft", legend = groups, col = tip_palette[groups], pch = 15,
  ncol = min(4, ceiling(length(groups) / 3)), bty = "n", cex = 0.55,
  title = "Tip-label groups"
)
dev.off()

# Figure 2: all CV curves plus their mean and among-run standard deviation.
replicate_ids <- sort(unique(cv$replicate))
replicate_colors <- setNames(
  viridisLite::turbo(length(replicate_ids)), replicate_ids
)
summary_k <- sort(unique(cv$K))
cv_mean <- vapply(summary_k, function(k) mean(cv$CV_error[cv$K == k]), numeric(1))
cv_sd <- vapply(summary_k, function(k) sd(cv$CV_error[cv$K == k]), numeric(1))
cv_sd[!is.finite(cv_sd)] <- 0
cv_n <- vapply(summary_k, function(k) sum(cv$K == k), integer(1))
cv_summary <- data.frame(
  K = summary_k, mean_CV_error = cv_mean,
  SD_CV_error = cv_sd, n_replicates = cv_n
)
cv_summary_output <- paste0(output_prefix, "_cv_summary.tsv")
write.table(cv_summary, cv_summary_output, sep = "\t", row.names = FALSE,
            quote = FALSE)
cv_range <- range(c(cv$CV_error, cv_mean - cv_sd, cv_mean + cv_sd))
padding <- max(diff(cv_range) * 0.08, 0.0001)

pdf(cv_output, width = 9.5, height = 5.8, useDingbats = FALSE)
par(mar = c(4.2, 4.5, 2.8, 7.0), mgp = c(2.5, 0.75, 0),
    tcl = -0.25, xpd = NA)
plot(
  NA, xlim = range(cv$K), ylim = cv_range + c(-padding, padding),
  xlab = "K", ylab = "10-fold CV error", xaxt = "n",
  main = "ADMIXTURE cross-validation"
)
axis(1, at = sort(unique(cv$K)))
abline(h = pretty(cv_range), col = "#e3e3e3", lwd = 0.7)
polygon(
  c(summary_k, rev(summary_k)),
  c(cv_mean - cv_sd, rev(cv_mean + cv_sd)),
  col = adjustcolor("#737373", alpha.f = 0.22), border = NA
)
for (replicate in replicate_ids) {
  values <- cv[cv$replicate == replicate, , drop = FALSE]
  values <- values[order(values$K), ]
  lines(values$K, values$CV_error, type = "b", pch = 16, cex = 0.48,
        lwd = 0.8,
        col = adjustcolor(replicate_colors[as.character(replicate)], 0.68))
}
lines(summary_k, cv_mean, type = "b", pch = 16, cex = 0.8,
      lwd = 2.2, col = "#111111")
legend(
  "topright", inset = c(-0.34, 0),
  legend = c(paste0("R", replicate_ids), "Mean", "Mean +/- 1 SD"),
  col = c(replicate_colors, "#111111", adjustcolor("#737373", 0.35)),
  pch = c(rep(16, length(replicate_ids)), 16, 15),
  lty = c(rep(1, length(replicate_ids)), 1, NA),
  lwd = c(rep(0.8, length(replicate_ids)), 2.2, NA),
  ncol = 2, bty = "n", cex = 0.72, title = "Replicate"
)
dev.off()

message("Tree figure:      ", tree_output)
message("Tree tip order:   ", paste0(output_prefix, "_tree_tip_order.txt"))
message("CV figure:        ", cv_output)
message("CV summary:       ", cv_summary_output)
message("ADMIXTURE order:  ", paste0(output_prefix, "_admixture_sample_order.txt"))
