#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 06_plot_admixture.R ADMIXTURE_DIR METADATA.csv OUTPUT.pdf [K_MIN] [K_MAX]")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

admixture_dir <- args[1]
metadata_file <- args[2]
output_file <- args[3]
k_min <- if (length(args) >= 4) as.integer(args[4]) else -Inf
k_max <- if (length(args) >= 5) as.integer(args[5]) else Inf

cv <- read.delim(file.path(admixture_dir, "cross_validation.tsv"))
cv <- cv[cv$K >= k_min & cv$K <= k_max, ]
best <- do.call(rbind, lapply(split(cv, cv$K), function(x) x[which.min(x$CV_error), ]))
sample_order <- readLines(file.path(admixture_dir, "sample_order.txt"))
metadata <- read.csv(metadata_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
metadata <- metadata[match(sample_order, metadata$sampleID), ]
if (any(is.na(metadata$sampleID))) stop("Some ADMIXTURE samples are absent from metadata")

long <- list()
index <- 1
for (i in seq_len(nrow(best))) {
  k <- best$K[i]
  rep <- best$replicate[i]
  q_file <- file.path(admixture_dir, paste0("rep", rep), paste0("K", k, ".rep", rep, ".Q"))
  q <- as.matrix(read.table(q_file, header = FALSE))
  if (nrow(q) != length(sample_order)) stop(paste("Sample count mismatch in", q_file))
  for (component in seq_len(ncol(q))) {
    long[[index]] <- data.frame(
      sample_index = seq_along(sample_order),
      sampleID = sample_order,
      K = factor(paste0("K = ", k), levels = paste0("K = ", sort(unique(best$K)))),
      component = factor(component),
      ancestry = q[, component]
    )
    index <- index + 1
  }
}
plot_data <- do.call(rbind, long)
palette <- hcl.colors(max(best$K), "Dynamic")

p <- ggplot(plot_data, aes(sample_index, ancestry, fill = component)) +
  geom_col(width = 1) +
  facet_grid(K ~ ., scales = "free_y") +
  scale_fill_manual(values = palette, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Samples in VCF order", y = "Ancestry proportion") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.spacing.y = grid::unit(0.08, "lines"),
    strip.background = element_blank(), strip.text.y = element_text(angle = 0)
  )

height <- max(5, 1.15 * length(unique(plot_data$K)))
ggsave(output_file, p, width = 14, height = height, limitsize = FALSE)

cv_plot <- ggplot(cv, aes(K, CV_error, group = factor(replicate), color = factor(replicate))) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = sort(unique(cv$K))) +
  labs(x = "K", y = "10-fold CV error", color = "Replicate") +
  theme_classic()
cv_file <- sub("\\.pdf$", "_CV.pdf", output_file, ignore.case = TRUE)
ggsave(cv_file, cv_plot, width = 7, height = 4.5)

