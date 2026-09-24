#!/usr/bin/env Rscript

# Produce three separate Pan-Dryas figures:
#   1. rooted and ladderized ML tree with s170/s47 annotations;
#   2. ADMIXTURE cross-validation curves;
#   3. labeled ADMIXTURE proportions ordered by the ML tree.

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
plot_replicates_text <- Sys.getenv("PLOT_REPLICATES", "1,2")

# Retain compatibility with the previous command, which supplied a .pdf name.
output_prefix <- sub("\\.pdf$", "", output_argument, ignore.case = TRUE)
tree_output <- paste0(output_prefix, "_tree.pdf")
cv_output <- paste0(output_prefix, "_cv.pdf")
admixture_output <- paste0(output_prefix, "_admixture.pdf")
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

safe_correlation <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(-Inf)
  cor(x, y, use = "pairwise.complete.obs")
}

# Greedily match components at adjacent K values by their sample profiles.
align_to_previous <- function(previous, current) {
  n_previous <- ncol(previous)
  n_current <- ncol(current)
  scores <- matrix(-Inf, n_previous, n_current)
  for (i in seq_len(n_previous)) {
    for (j in seq_len(n_current)) {
      scores[i, j] <- safe_correlation(previous[, i], current[, j])
    }
  }
  previous_used <- rep(FALSE, n_previous)
  current_used <- rep(FALSE, n_current)
  current_for_previous <- rep(NA_integer_, n_previous)
  for (step in seq_len(min(n_previous, n_current))) {
    available <- scores
    available[previous_used, ] <- -Inf
    available[, current_used] <- -Inf
    chosen <- which(available == max(available), arr.ind = TRUE)[1, ]
    current_for_previous[chosen[1]] <- chosen[2]
    previous_used[chosen[1]] <- TRUE
    current_used[chosen[2]] <- TRUE
  }
  remaining <- which(!current_used)
  current[, c(current_for_previous, remaining), drop = FALSE]
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

# Plot every run in K-major, replicate-minor order.
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
runs$key <- paste0("K", runs$K, "_R", runs$replicate)
k_values <- sort(unique(runs$K))

plot_replicates <- suppressWarnings(as.integer(trimws(strsplit(
  plot_replicates_text, ",", fixed = TRUE
)[[1]])))
if (any(!is.finite(plot_replicates)) || !length(plot_replicates)) {
  stop("PLOT_REPLICATES must be a comma-separated list such as 1,2")
}
missing_plot_replicates <- setdiff(plot_replicates, all_replicates)
if (length(missing_plot_replicates)) {
  stop("Requested PLOT_REPLICATES not found: ",
       paste(missing_plot_replicates, collapse = ", "))
}
plot_runs <- runs[runs$replicate %in% plot_replicates, , drop = FALSE]
plot_runs <- plot_runs[order(plot_runs$K, plot_runs$replicate), , drop = FALSE]
message("ADMIXTURE panel order: ", paste(plot_runs$key, collapse = ", "))

q_matrices <- list()
for (i in seq_len(nrow(plot_runs))) {
  k <- plot_runs$K[i]
  replicate <- plot_runs$replicate[i]
  q_file <- file.path(
    admixture_dir, paste0("rep", replicate),
    paste0("K", k, ".rep", replicate, ".Q")
  )
  q <- as.matrix(read.table(q_file, header = FALSE))
  if (nrow(q) != length(sample_order)) {
    stop("Sample count does not match sample_order.txt in ", q_file)
  }
  rownames(q) <- sample_order
  q_matrices[[plot_runs$key[i]]] <- q
}

# Establish component order at the first K from the tree-ordered samples, then
# preserve comparable colors across adjacent K values.
ingroup_tip_ids <- tip_ids[tip_ids %in% sample_order]
if (length(ingroup_tip_ids) != length(sample_order)) {
  missing <- setdiff(sample_order, tip_ids)
  stop("ADMIXTURE samples absent from tree: ", paste(missing, collapse = ", "))
}
first_key <- plot_runs$key[1]
first_q <- q_matrices[[first_key]]
tree_ordered_q <- first_q[match(ingroup_tip_ids, rownames(first_q)), , drop = FALSE]
first_order <- order(apply(tree_ordered_q, 2, which.max))
q_matrices[[first_key]] <- first_q[, first_order, drop = FALSE]
if (nrow(plot_runs) > 1) {
  for (i in 2:nrow(plot_runs)) {
    previous <- q_matrices[[plot_runs$key[i - 1]]]
    current <- q_matrices[[plot_runs$key[i]]]
    q_matrices[[plot_runs$key[i]]] <- align_to_previous(previous, current)
  }
}

# Use the vivid palette from the earlier ADMIXTURE plotting workflow.
if (requireNamespace("viridisLite", quietly = TRUE)) {
  ancestry_colors <- viridisLite::turbo(max(k_values))
} else {
  warning("viridisLite not installed; using the base-R Turbo palette")
  ancestry_colors <- hcl.colors(max(k_values), "Turbo")
}

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
      frame = "none", cex = 0.18, adj = c(1.05, -0.15)
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

pdf(cv_output, width = 7.5, height = 5.5, useDingbats = FALSE)
par(mar = c(4.2, 4.5, 2.8, 1), mgp = c(2.5, 0.75, 0), tcl = -0.25)
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
  "topright", legend = c("Individual runs", "Mean", "Mean ± 1 SD"),
  col = c("#636363", "#111111", adjustcolor("#737373", 0.35)),
  pch = c(16, 16, 15), lty = c(1, 1, NA), lwd = c(0.8, 2.2, NA),
  bty = "n", cex = 0.8
)
dev.off()

# Figure 3: tree-ordered ADMIXTURE panels with sample labels on the left.
ingroup_metadata <- metadata[match(ingroup_tip_ids, metadata$sampleID), , drop = FALSE]
ingroup_labels <- ingroup_metadata[[tip_label_column]]
ingroup_labels[is.na(ingroup_labels) | ingroup_labels == ""] <-
  ingroup_tip_ids[is.na(ingroup_labels) | ingroup_labels == ""]
ingroup_colors <- unname(tip_palette[as.character(ingroup_metadata$spID)])
n_samples <- length(ingroup_tip_ids)
n_panels <- nrow(plot_runs)
sample_y <- seq_len(n_samples)

pdf(admixture_output, width = max(24, 7.2 + 0.82 * n_panels), height = 32,
    useDingbats = FALSE)
layout(matrix(seq_len(n_panels + 1), nrow = 1),
       widths = c(7.2, rep(0.82, n_panels)))

par(mar = c(0.6, 0.2, 1.6, 0.1), xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.5, n_samples + 0.5))
text(0.99, sample_y, labels = ingroup_labels, adj = c(1, 0.5),
     cex = 0.30, col = ingroup_colors)
title("Samples in ML-tree order", cex.main = 0.85, line = 0.25)

for (i in seq_len(nrow(plot_runs))) {
  k <- plot_runs$K[i]
  replicate <- plot_runs$replicate[i]
  q <- q_matrices[[plot_runs$key[i]]]
  par(mar = c(0.6, 0.03, 1.6, 0.03), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.5, n_samples + 0.5))
  rect(0, 0.5, 1, n_samples + 0.5, col = "#eeeeee", border = NA)
  for (sample_index in seq_along(ingroup_tip_ids)) {
    q_row <- match(ingroup_tip_ids[sample_index], rownames(q))
    left <- 0
    for (component in seq_len(ncol(q))) {
      right <- left + q[q_row, component]
      if (right > left) {
        rect(left, sample_y[sample_index] - 0.49,
             right, sample_y[sample_index] + 0.49,
             col = ancestry_colors[component], border = NA)
      }
      left <- right
    }
  }
  title(paste0("K=", k, "\nR", replicate), cex.main = 0.78, line = 0.05)
  box(col = "white")
}
dev.off()

message("Tree figure:      ", tree_output)
message("CV figure:        ", cv_output)
message("CV summary:       ", cv_summary_output)
message("ADMIXTURE figure: ", admixture_output)
