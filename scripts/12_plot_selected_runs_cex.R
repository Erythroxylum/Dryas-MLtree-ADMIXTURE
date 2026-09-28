#!/usr/bin/env Rscript

# Plot the seven selected ADMIXTURE solutions in ascending K order.
# No tree or CV-error panel is produced.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 5) {
  stop(paste(
    "Usage: 12_plot_selected_runs_cex.R CROSS_VALIDATION.tsv METADATA.csv",
    "ORDER_FILE OUTPUT.pdf [SAMPLE_LABEL_CEX]"
  ))
}

cv_file <- args[1]
metadata_file <- args[2]
order_file <- args[3]
output_file <- args[4]
tip_label_column <- Sys.getenv("TIP_LABEL_COLUMN", "phyloID")
sample_label_cex <- if (length(args) == 5) {
  as.numeric(args[5])
} else {
  as.numeric(Sys.getenv("SAMPLE_LABEL_CEX", "0.40"))
}
if (!is.finite(sample_label_cex) || sample_label_cex <= 0) {
  stop("SAMPLE_LABEL_CEX must be a positive number")
}
admixture_dir <- dirname(normalizePath(cv_file))

selected <- data.frame(
  K = 16:22,
  replicate = c(4L, 1L, 5L, 10L, 7L, 4L, 3L)
)
selected$key <- paste0("K", selected$K, "_R", selected$replicate)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

safe_correlation <- function(x, y) {
  sx <- sd(x, na.rm = TRUE)
  sy <- sd(y, na.rm = TRUE)
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(-Inf)
  value <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
  if (!is.finite(value)) return(-Inf)
  value
}

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
  matched_previous <- integer(0)
  matched_current <- integer(0)

  for (step in seq_len(min(n_previous, n_current))) {
    available <- scores
    available[previous_used, ] <- -Inf
    available[, current_used] <- -Inf
    chosen <- which(available == max(available), arr.ind = TRUE)[1, ]
    matched_previous <- c(matched_previous, chosen[1])
    matched_current <- c(matched_current, chosen[2])
    previous_used[chosen[1]] <- TRUE
    current_used[chosen[2]] <- TRUE
  }

  matched_current <- matched_current[order(matched_previous)]
  remaining <- which(!current_used)
  new_order <- c(matched_current, remaining)
  if (length(new_order) != n_current || anyDuplicated(new_order)) {
    stop("Internal error while matching ancestry components")
  }
  current[, new_order, drop = FALSE]
}

find_q_file <- function(k, replicate) {
  candidates <- c(
    file.path(
      admixture_dir, paste0("rep", replicate),
      paste0("K", k, ".rep", replicate, ".Q")
    ),
    file.path(admixture_dir, paste0("K", k, ".rep", replicate, ".Q"))
  )
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("Q file not found. Checked: ", paste(candidates, collapse = ", "))
  }
  found[1]
}

cv <- read.delim(cv_file, stringsAsFactors = FALSE)
required_cv <- c("K", "replicate", "CV_error")
missing_cv <- setdiff(required_cv, names(cv))
if (length(missing_cv)) {
  stop("CV columns missing: ", paste(missing_cv, collapse = ", "))
}

selected$CV_error <- NA_real_
for (i in seq_len(nrow(selected))) {
  hit <- cv$K == selected$K[i] & cv$replicate == selected$replicate[i]
  if (sum(hit) != 1) {
    stop(
      "Expected exactly one CV entry for K=", selected$K[i],
      ", replicate=", selected$replicate[i]
    )
  }
  selected$CV_error[i] <- cv$CV_error[hit]
}

message("Selected runs:")
for (i in seq_len(nrow(selected))) {
  message(sprintf(
    "  K=%d  replicate=%d  CV_error=%.5f",
    selected$K[i], selected$replicate[i], selected$CV_error[i]
  ))
}

q_order_file <- file.path(admixture_dir, "sample_order.txt")
if (!file.exists(q_order_file)) stop("Missing Q-row order file: ", q_order_file)
q_sample_order <- readLines(q_order_file)
plot_sample_order <- readLines(order_file)

if (anyDuplicated(q_sample_order)) stop("sample_order.txt contains duplicate IDs")
if (anyDuplicated(plot_sample_order)) stop("ORDER_FILE contains duplicate IDs")
if (!setequal(plot_sample_order, q_sample_order)) {
  stop("ORDER_FILE and sample_order.txt do not contain exactly the same samples")
}

metadata <- read.csv(
  metadata_file, check.names = FALSE, fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE
)
required_metadata <- c("sampleID", "spID")
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) {
  stop("Metadata columns missing: ", paste(missing_metadata, collapse = ", "))
}
if (!(tip_label_column %in% names(metadata))) tip_label_column <- "sampleID"

plot_metadata <- metadata[match(plot_sample_order, metadata$sampleID), , drop = FALSE]
if (any(is.na(plot_metadata$sampleID))) {
  stop("Some ORDER_FILE samples are absent from the metadata")
}

plot_labels <- plot_metadata[[tip_label_column]]
missing_labels <- is.na(plot_labels) | plot_labels == ""
plot_labels[missing_labels] <- plot_sample_order[missing_labels]

groups <- sort(unique(as.character(plot_metadata$spID)))
tip_palette <- setNames(hcl.colors(length(groups), "Dark 3"), groups)
label_colors <- unname(tip_palette[as.character(plot_metadata$spID)])

q_matrices <- list()
for (i in seq_len(nrow(selected))) {
  q_file <- find_q_file(selected$K[i], selected$replicate[i])
  q <- as.matrix(read.table(q_file, header = FALSE))
  if (nrow(q) != length(q_sample_order)) {
    stop("Sample count does not match sample_order.txt in ", q_file)
  }
  rownames(q) <- q_sample_order
  q_matrices[[selected$key[i]]] <- q
}

# Align ancestry-component columns in the displayed K16-to-K22 sequence.
first_key <- selected$key[1]
first_q <- q_matrices[[first_key]]
ordered_first_q <- first_q[match(plot_sample_order, rownames(first_q)), , drop = FALSE]
first_order <- order(apply(ordered_first_q, 2, which.max))
q_matrices[[first_key]] <- first_q[, first_order, drop = FALSE]

for (i in 2:nrow(selected)) {
  previous <- q_matrices[[selected$key[i - 1]]]
  current <- q_matrices[[selected$key[i]]]
  q_matrices[[selected$key[i]]] <- align_to_previous(previous, current)
}

max_k <- max(selected$K)
if (requireNamespace("viridisLite", quietly = TRUE)) {
  ancestry_colors <- viridisLite::turbo(max_k)
} else {
  ancestry_colors <- hcl.colors(max_k, "Turbo")
}

n_samples <- length(plot_sample_order)
n_panels <- nrow(selected)
sample_y <- seq_len(n_samples)

pdf(
  output_file,
  width = max(17, 7.2 + 1.15 * n_panels),
  height = 32,
  useDingbats = FALSE
)
layout(
  matrix(seq_len(n_panels + 1), nrow = 1),
  widths = c(7.2, rep(1.15, n_panels))
)

par(mar = c(0.6, 0.2, 2.3, 0.1), xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.5, n_samples + 0.5))
text(
  0.99, sample_y, labels = plot_labels, adj = c(1, 0.5),
  cex = sample_label_cex, col = label_colors
)
title("Samples in active-tree order", cex.main = 0.85, line = 0.25)

for (i in seq_len(nrow(selected))) {
  k <- selected$K[i]
  replicate <- selected$replicate[i]
  q <- q_matrices[[selected$key[i]]]

  par(mar = c(0.6, 0.03, 2.3, 0.03), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.5, n_samples + 0.5))
  rect(0, 0.5, 1, n_samples + 0.5, col = "#eeeeee", border = NA)

  for (sample_index in seq_along(plot_sample_order)) {
    q_row <- match(plot_sample_order[sample_index], rownames(q))
    left <- 0
    for (component in seq_len(ncol(q))) {
      right <- left + q[q_row, component]
      if (right > left) {
        rect(
          left, sample_y[sample_index] - 0.49,
          right, sample_y[sample_index] + 0.49,
          col = ancestry_colors[component], border = NA
        )
      }
      left <- right
    }
  }

  title(
    sprintf("K=%d\nR%d\nCV=%.5f", k, replicate, selected$CV_error[i]),
    cex.main = 0.72, line = 0.05
  )
  box(col = "black", lwd = 0.8)
}

dev.off()
message("ADMIXTURE bar chart written to: ", output_file)
