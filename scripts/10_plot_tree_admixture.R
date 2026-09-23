#!/usr/bin/env Rscript

# Composite Pan-Dryas figure: rooted ML tree, analysis-membership tracks,
# ADMIXTURE cross-validation, and ancestry proportions aligned to tree tips.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage: 10_plot_tree_admixture.R TREEFILE ADMIXTURE_DIR",
    "METADATA.csv OUTPUT.pdf [K_MIN] [K_MAX]"
  ))
}

suppressPackageStartupMessages(library(ape))

tree_file <- args[1]
admixture_dir <- args[2]
metadata_file <- args[3]
output_file <- args[4]
k_min <- if (length(args) >= 5) as.integer(args[5]) else -Inf
k_max <- if (length(args) >= 6) as.integer(args[6]) else Inf
s170_column <- Sys.getenv("S170_COLUMN", "s170_BPP")
s47_column <- Sys.getenv("S47_COLUMN", "s47-p9")
tip_label_column <- Sys.getenv("TIP_LABEL_COLUMN", "phyloID")
minimum_support <- as.numeric(Sys.getenv("MIN_SUPPORT", "70"))

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

# Greedily match current components to the aligned components at the preceding K.
# Since consecutive K values differ by one component, this provides stable colors
# without introducing a separate CLUMPP dependency.
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

root_override <- trimws(Sys.getenv("OUTGROUPS", ""))
if (nzchar(root_override)) {
  outgroups <- trimws(strsplit(root_override, ",", fixed = TRUE)[[1]])
} else {
  outgroups <- metadata$sampleID[metadata$s384 != "" & !is.na(metadata$s384) &
                                  (is.na(metadata$s380) | metadata$s380 == "")]
}
outgroups <- intersect(outgroups, tree$tip.label)
if (!length(outgroups)) {
  stop("No outgroups found. Set OUTGROUPS to comma-separated tree-tip IDs.")
}
tree <- tryCatch(
  root(tree, outgroup = outgroups, resolve.root = TRUE),
  error = function(e) {
    stop(
      "The selected outgroups are not a rootable clade: ", conditionMessage(e),
      ". Set OUTGROUPS to the appropriate comma-separated tip IDs."
    )
  }
)
tree <- ladderize(tree, right = TRUE)
tip_ids <- tree$tip.label

tree_metadata <- metadata[match(tip_ids, metadata$sampleID), , drop = FALSE]
if (any(is.na(tree_metadata$sampleID))) {
  missing <- tip_ids[is.na(tree_metadata$sampleID)]
  stop("Tree tips absent from metadata: ", paste(missing, collapse = ", "))
}

sample_order <- readLines(file.path(admixture_dir, "sample_order.txt"))
cv <- read.delim(file.path(admixture_dir, "cross_validation.tsv"),
                 stringsAsFactors = FALSE)
cv <- cv[cv$K >= k_min & cv$K <= k_max, , drop = FALSE]
if (!nrow(cv)) stop("No ADMIXTURE K values remain after filtering")
if (any(!is.finite(cv$CV_error))) stop("Non-numeric or missing CV errors found")
best <- do.call(rbind, lapply(split(cv, cv$K), function(x) {
  x[which.min(x$CV_error), , drop = FALSE]
}))
best <- best[order(best$K), , drop = FALSE]
k_values <- best$K

q_matrices <- list()
for (i in seq_len(nrow(best))) {
  k <- best$K[i]
  replicate <- best$replicate[i]
  q_file <- file.path(
    admixture_dir, paste0("rep", replicate),
    paste0("K", k, ".rep", replicate, ".Q")
  )
  q <- as.matrix(read.table(q_file, header = FALSE))
  if (nrow(q) != length(sample_order)) {
    stop("Sample count does not match sample_order.txt in ", q_file)
  }
  rownames(q) <- sample_order
  q_matrices[[as.character(k)]] <- q
}

# Establish the first component order by the vertical position of its maximum,
# then inherit component colors across successive K values.
ingroup_tip_order <- tip_ids[tip_ids %in% sample_order]
first_key <- as.character(k_values[1])
first_q <- q_matrices[[first_key]]
tree_ordered_q <- first_q[match(ingroup_tip_order, rownames(first_q)), , drop = FALSE]
first_order <- order(apply(tree_ordered_q, 2, which.max))
q_matrices[[first_key]] <- first_q[, first_order, drop = FALSE]
if (length(k_values) > 1) {
  for (i in 2:length(k_values)) {
    previous <- q_matrices[[as.character(k_values[i - 1])]]
    current <- q_matrices[[as.character(k_values[i])]]
    q_matrices[[as.character(k_values[i])]] <- align_to_previous(previous, current)
  }
}

groups <- sort(unique(as.character(tree_metadata$spID)))
tip_palette <- setNames(hcl.colors(length(groups), "Dark 3"), groups)
tip_colors <- unname(tip_palette[as.character(tree_metadata$spID)])
ancestry_colors <- hcl.colors(max(k_values), "Dynamic")
s170_membership <- is_member(tree_metadata[[s170_column]])
s47_membership <- is_member(tree_metadata[[s47_column]])

display_tree <- tree
display_labels <- tree_metadata[[tip_label_column]]
display_labels[is.na(display_labels) | display_labels == ""] <-
  tip_ids[is.na(display_labels) | display_labels == ""]
display_tree$tip.label <- display_labels

n_k <- length(k_values)
layout_matrix <- rbind(
  c(1, rep(2, 2 + n_k)),
  c(3, 4, 5, seq.int(6, 5 + n_k))
)
figure_width <- max(15, 8.0 + 0.82 * n_k)
figure_height <- 30

pdf(output_file, width = figure_width, height = figure_height,
    useDingbats = FALSE, onefile = TRUE)
on.exit(dev.off(), add = TRUE)
layout(
  layout_matrix,
  widths = c(7.2, 0.32, 0.32, rep(0.82, n_k)),
  heights = c(4.2, 25.8)
)

# Panel 1: CV error for both replicates and the selected minimum-CV run.
par(mar = c(3.2, 3.6, 2.0, 0.6), mgp = c(2.0, 0.65, 0), tcl = -0.25)
cv_range <- range(cv$CV_error)
padding <- max(diff(cv_range) * 0.08, 0.0001)
plot(
  NA, xlim = range(cv$K), ylim = cv_range + c(-padding, padding),
  xlab = "K", ylab = "10-fold CV error", xaxt = "n",
  main = "ADMIXTURE cross-validation", cex.main = 0.9
)
axis(1, at = sort(unique(cv$K)), cex.axis = 0.75)
abline(h = pretty(cv_range), col = "#e6e6e6", lwd = 0.6)
replicate_colors <- setNames(
  hcl.colors(length(unique(cv$replicate)), "Blues 3"),
  sort(unique(cv$replicate))
)
for (replicate in sort(unique(cv$replicate))) {
  values <- cv[cv$replicate == replicate, , drop = FALSE]
  values <- values[order(values$K), ]
  lines(values$K, values$CV_error, type = "b", pch = 1, cex = 0.65,
        lwd = 0.8, col = replicate_colors[as.character(replicate)])
}
lines(best$K, best$CV_error, type = "b", pch = 16, cex = 0.7,
      lwd = 1.4, col = "#111111")
legend(
  "topright",
  legend = c(paste0("Replicate ", sort(unique(cv$replicate))), "Selected"),
  col = c(replicate_colors, "#111111"),
  pch = c(rep(1, length(replicate_colors)), 16),
  lty = 1, bty = "n", cex = 0.65
)

# Panel 2: compact legends and selected replicate summary.
par(mar = c(0.2, 0.5, 0.2, 0.2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
legend(
  "topleft", legend = groups, col = tip_palette[groups], pch = 15,
  ncol = min(5, ceiling(length(groups) / 2)), bty = "n",
  cex = 0.62, title = "Tip-label groups"
)
legend(
  "bottomleft",
  legend = c("s170", "s47"), col = c("#2166ac", "#b2182b"),
  pch = 15, bty = "n", horiz = TRUE, cex = 0.72,
  title = "Analysis membership"
)
replicate_summary <- paste0("K", best$K, ":R", best$replicate, collapse = "   ")
text(
  0.99, 0.12,
  labels = paste0("Lowest-CV runs:  ", replicate_summary),
  adj = c(1, 0), cex = 0.62
)

# Panel 3: rooted, ladderized tree. Tip coordinates define all later panels.
par(mar = c(0.4, 0.4, 1.5, 0.2), xpd = NA)
plot(
  display_tree, type = "phylogram", direction = "rightwards",
  show.tip.label = TRUE, tip.color = tip_colors,
  cex = 0.23, align.tip.label = TRUE,
  label.offset = max(node.depth.edgelength(tree)) * 0.006,
  no.margin = TRUE
)
title("Rooted and ladderized maximum-likelihood tree", cex.main = 0.9, line = 0.2)
add.scale.bar(cex = 0.55, lwd = 0.8)
last_tree_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- last_tree_plot$yy[seq_len(Ntip(tree))]
tree_ylim <- range(last_tree_plot$yy) + c(-0.5, 0.5)

if (!is.null(tree$node.label)) {
  extract_support <- function(label) {
    pieces <- strsplit(as.character(label), "/", fixed = TRUE)[[1]]
    suppressWarnings(as.numeric(tail(pieces, 1)))
  }
  support <- vapply(tree$node.label, extract_support, numeric(1))
  show_support <- which(is.finite(support) & support >= minimum_support)
  if (length(show_support)) {
    nodelabels(
      text = round(support[show_support]),
      node = Ntip(tree) + show_support,
      frame = "none", cex = 0.16, adj = c(1.05, -0.15)
    )
  }
}

plot_membership_track <- function(values, heading, color) {
  par(mar = c(0.4, 0, 1.5, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = tree_ylim)
  rect(0, tree_ylim[1], 1, tree_ylim[2], col = "#f3f3f3", border = NA)
  points(rep(0.5, sum(values)), tip_y[values], pch = 15, cex = 0.34, col = color)
  title(heading, cex.main = 0.7, line = 0.15)
  box(col = "white")
}

# Panels 4-5: s170 and s47 analysis membership.
plot_membership_track(s170_membership, "s170", "#2166ac")
plot_membership_track(s47_membership, "s47", "#b2182b")

# Remaining panels: ancestry proportions in the tree-tip order.
for (k in k_values) {
  q <- q_matrices[[as.character(k)]]
  par(mar = c(0.4, 0.03, 1.5, 0.03), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = tree_ylim)
  rect(0, tree_ylim[1], 1, tree_ylim[2], col = "#eeeeee", border = NA)
  for (tip_index in seq_along(tip_ids)) {
    q_row <- match(tip_ids[tip_index], rownames(q))
    if (is.na(q_row)) next
    left <- 0
    for (component in seq_len(ncol(q))) {
      right <- left + q[q_row, component]
      if (right > left) {
        rect(
          left, tip_y[tip_index] - 0.49,
          right, tip_y[tip_index] + 0.49,
          col = ancestry_colors[component], border = NA
        )
      }
      left <- right
    }
  }
  title(paste0("K=", k), cex.main = 0.72, line = 0.15)
  box(col = "white")
}

dev.off()
on.exit(NULL, add = FALSE)
message("Composite figure written to: ", output_file)
