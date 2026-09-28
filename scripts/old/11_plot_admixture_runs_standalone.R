#!/usr/bin/env Rscript

# Plot ADMIXTURE runs without reading or plotting a tree.
# The supplied order file must contain one sample ID per line in the desired
# bottom-to-top plot order, typically *_admixture_sample_order.txt generated
# previously from the active rooted and ladderized tree.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste(
    "Usage: 11_plot_admixture_runs.R MODE ADMIXTURE_DIR METADATA.csv",
    "ORDER_FILE OUTPUT.pdf K_MIN K_MAX",
    "where MODE is 'all' or 'lowest_cv'"
  ))
}

mode <- tolower(trimws(args[1]))
admixture_dir <- args[2]
metadata_file <- args[3]
order_file <- args[4]
output_file <- args[5]
k_min <- as.integer(args[6])
k_max <- as.integer(args[7])
tip_label_column <- Sys.getenv("TIP_LABEL_COLUMN", "phyloID")

if (!(mode %in% c("all", "lowest_cv"))) {
  stop("MODE must be either 'all' or 'lowest_cv'")
}
if (!is.finite(k_min) || !is.finite(k_max) || k_min > k_max) {
  stop("K_MIN and K_MAX must define a valid integer range")
}

output_prefix <- sub("\\.pdf$", "", output_file, ignore.case = TRUE)
cv_output <- paste0(output_prefix, "_CV.pdf")
cv_summary_output <- paste0(output_prefix, "_CV_summary.tsv")
selected_runs_output <- paste0(output_prefix, "_plotted_runs.tsv")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

safe_correlation <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(-Inf)
  cor(x, y, use = "pairwise.complete.obs")
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

q_sample_order <- readLines(file.path(admixture_dir, "sample_order.txt"))
plot_sample_order <- readLines(order_file)
if (anyDuplicated(q_sample_order)) stop("sample_order.txt contains duplicate IDs")
if (anyDuplicated(plot_sample_order)) stop("ORDER_FILE contains duplicate IDs")
if (!setequal(plot_sample_order, q_sample_order)) {
  stop("ORDER_FILE and sample_order.txt do not contain exactly the same samples")
}

plot_metadata <- metadata[match(plot_sample_order, metadata$sampleID), , drop = FALSE]
if (any(is.na(plot_metadata$sampleID))) {
  stop("Some ORDER_FILE samples are absent from the metadata")
}
plot_labels <- plot_metadata[[tip_label_column]]
plot_labels[is.na(plot_labels) | plot_labels == ""] <-
  plot_sample_order[is.na(plot_labels) | plot_labels == ""]
groups <- sort(unique(as.character(plot_metadata$spID)))
tip_palette <- setNames(hcl.colors(length(groups), "Dark 3"), groups)
label_colors <- unname(tip_palette[as.character(plot_metadata$spID)])

cv <- read.delim(
  file.path(admixture_dir, "cross_validation.tsv"),
  stringsAsFactors = FALSE
)
cv <- cv[cv$K >= k_min & cv$K <= k_max, , drop = FALSE]
if (!nrow(cv)) stop("No CV results remain in the requested K range")
if (any(!is.finite(cv$CV_error))) stop("CV_error contains missing values")
runs <- cv[order(cv$K, cv$replicate), c("K", "replicate", "CV_error"),
           drop = FALSE]
if (anyDuplicated(runs[c("K", "replicate")])) {
  stop("cross_validation.tsv contains duplicate K/replicate combinations")
}

replicates <- sort(unique(runs$replicate))
replicates_by_k <- split(runs$replicate, runs$K)
incomplete_k <- names(replicates_by_k)[!vapply(
  replicates_by_k, function(x) setequal(x, replicates), logical(1)
)]
if (length(incomplete_k)) {
  stop("Incomplete replicate set at K = ", paste(incomplete_k, collapse = ", "))
}
runs$key <- paste0("K", runs$K, "_R", runs$replicate)

lowest_cv_runs <- do.call(rbind, lapply(split(runs, runs$K), function(x) {
  x[which.min(x$CV_error), , drop = FALSE]
}))
lowest_cv_runs <- lowest_cv_runs[order(lowest_cv_runs$K), , drop = FALSE]

if (mode == "all") {
  plot_runs <- runs
} else {
  plot_runs <- lowest_cv_runs
  if (nrow(plot_runs) != length(unique(runs$K)) || anyDuplicated(plot_runs$K)) {
    stop("Internal error: lowest_cv mode did not select exactly one run per K")
  }
}
write.table(
  plot_runs[c("K", "replicate", "CV_error")], selected_runs_output,
  sep = "\t", row.names = FALSE, quote = FALSE
)
message("Plot mode: ", mode)
message("Panel order: ", paste(plot_runs$key, collapse = ", "))

q_matrices <- list()
for (i in seq_len(nrow(plot_runs))) {
  k <- plot_runs$K[i]
  replicate <- plot_runs$replicate[i]
  q_file <- file.path(
    admixture_dir, paste0("rep", replicate),
    paste0("K", k, ".rep", replicate, ".Q")
  )
  if (!file.exists(q_file)) stop("Q file not found: ", q_file)
  q <- as.matrix(read.table(q_file, header = FALSE))
  if (nrow(q) != length(q_sample_order)) {
    stop("Sample count does not match sample_order.txt in ", q_file)
  }
  rownames(q) <- q_sample_order
  q_matrices[[plot_runs$key[i]]] <- q
}

# Align component labels sequentially across the displayed run order.
first_key <- plot_runs$key[1]
first_q <- q_matrices[[first_key]]
ordered_first_q <- first_q[match(plot_sample_order, rownames(first_q)), , drop = FALSE]
first_order <- order(apply(ordered_first_q, 2, which.max))
q_matrices[[first_key]] <- first_q[, first_order, drop = FALSE]
if (nrow(plot_runs) > 1) {
  for (i in 2:nrow(plot_runs)) {
    previous <- q_matrices[[plot_runs$key[i - 1]]]
    current <- q_matrices[[plot_runs$key[i]]]
    q_matrices[[plot_runs$key[i]]] <- align_to_previous(previous, current)
  }
}

max_k <- max(runs$K)
if (requireNamespace("viridisLite", quietly = TRUE)) {
  ancestry_colors <- viridisLite::turbo(max_k)
  replicate_colors <- viridisLite::turbo(length(replicates))
} else {
  ancestry_colors <- hcl.colors(max_k, "Turbo")
  replicate_colors <- hcl.colors(length(replicates), "Dark 3")
}
replicate_colors <- setNames(replicate_colors, replicates)

# CV plot: every replicate has an explicit color entry in the legend.
summary_k <- sort(unique(cv$K))
cv_mean <- vapply(summary_k, function(k) mean(cv$CV_error[cv$K == k]), numeric(1))
cv_sd <- vapply(summary_k, function(k) sd(cv$CV_error[cv$K == k]), numeric(1))
cv_sd[!is.finite(cv_sd)] <- 0
cv_n <- vapply(summary_k, function(k) sum(cv$K == k), integer(1))
cv_summary <- data.frame(
  K = summary_k, mean_CV_error = cv_mean,
  SD_CV_error = cv_sd, n_replicates = cv_n
)
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
axis(1, at = summary_k)
abline(h = pretty(cv_range), col = "#e3e3e3", lwd = 0.7)
polygon(
  c(summary_k, rev(summary_k)),
  c(cv_mean - cv_sd, rev(cv_mean + cv_sd)),
  col = adjustcolor("#737373", alpha.f = 0.20), border = NA
)
for (replicate in replicates) {
  values <- cv[cv$replicate == replicate, , drop = FALSE]
  values <- values[order(values$K), ]
  lines(values$K, values$CV_error, type = "b", pch = 16, cex = 0.55,
        lwd = 1.0, col = replicate_colors[as.character(replicate)])
}
lines(summary_k, cv_mean, type = "b", pch = 16, cex = 0.78,
      lwd = 2.3, col = "black")
if (mode == "lowest_cv") {
  points(
    lowest_cv_runs$K, lowest_cv_runs$CV_error,
    pch = 21, cex = 1.15, lwd = 1.2, col = "black",
    bg = replicate_colors[as.character(lowest_cv_runs$replicate)]
  )
}
legend(
  "topright", inset = c(-0.34, 0),
  legend = c(paste0("R", replicates), "Mean", "Mean +/- 1 SD"),
  col = c(replicate_colors, "black", adjustcolor("#737373", 0.35)),
  pch = c(rep(16, length(replicates)), 16, 15),
  lty = c(rep(1, length(replicates)), 1, NA),
  lwd = c(rep(1, length(replicates)), 2.3, NA),
  ncol = 2, bty = "n", cex = 0.72, title = "Replicate"
)
dev.off()

# Ancestry panels in the supplied bottom-to-top sample order.
n_samples <- length(plot_sample_order)
n_panels <- nrow(plot_runs)
sample_y <- seq_len(n_samples)
panel_width <- if (mode == "all") 0.72 else 1.0
pdf(
  output_file,
  width = max(17, 7.2 + panel_width * n_panels), height = 32,
  useDingbats = FALSE
)
layout(
  matrix(seq_len(n_panels + 1), nrow = 1),
  widths = c(7.2, rep(panel_width, n_panels))
)

par(mar = c(0.6, 0.2, 1.7, 0.1), xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0.5, n_samples + 0.5))
text(0.99, sample_y, labels = plot_labels, adj = c(1, 0.5),
     cex = 0.40, col = label_colors)
title("Samples in active-tree order", cex.main = 0.85, line = 0.25)

for (i in seq_len(nrow(plot_runs))) {
  k <- plot_runs$K[i]
  replicate <- plot_runs$replicate[i]
  q <- q_matrices[[plot_runs$key[i]]]
  par(mar = c(0.6, 0.02, 1.7, 0.02), xaxs = "i", yaxs = "i")
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
  title(paste0("K=", k, "\nR", replicate), cex.main = 0.78, line = 0.05)
  box(col = "white")
  end_of_k <- i == nrow(plot_runs) || plot_runs$K[i + 1] != k
  if (end_of_k && i < nrow(plot_runs)) {
    segments(1, 0.5, 1, n_samples + 0.5,
             col = "black", lwd = 3, xpd = NA)
  }
}
dev.off()

message("Ancestry figure:  ", output_file)
message("CV figure:        ", cv_output)
message("CV summary:       ", cv_summary_output)
message("Plotted runs:     ", selected_runs_output)
