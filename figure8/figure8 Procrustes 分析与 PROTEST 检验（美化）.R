## Figure 8: Procrustes（代谢组 PCA vs 微生物 PCoA, Bray–Curtis）
## R 4.4.3
## Methods: vegan::procrustes + vegan::protest(permutations = 10000)
## Inputs (repo root /data):
## - 完整数据-分泌物.csv（代谢矩阵；样本列：YS0_F-1 ... YS100_F-3）
## - 完整数据-微生物.csv（OTU/ASV counts + taxonomy；样本列：RDYS_1-1 ... RDYS_5-3）
## Outputs (repo root /figure8):
## - Fig8_main.pdf + Fig8_main.tiff
## - (optional) FigS8_procrustes_residuals.pdf

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE)
}

find_project_root <- function(start_dir) {
  candidate <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in 0:10) {
    exu <- file.path(candidate, "data", "完整数据-分泌物.csv")
    mic <- file.path(candidate, "data", "完整数据-微生物.csv")
    if (file.exists(exu) && file.exists(mic)) return(candidate)

    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  stop(
    "Cannot find project root containing data/完整数据-分泌物.csv and data/完整数据-微生物.csv. ",
    "Start: ", start_dir
  )
}

script_dir <- get_script_dir()
project_root <- find_project_root(script_dir)
data_dir <- file.path(project_root, "data")
out_dir <- file.path(project_root, "figure8")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

exudate_path <- file.path(data_dir, "完整数据-分泌物.csv")
microbe_path <- file.path(data_dir, "完整数据-微生物.csv")

extract_metabolome_pca <- function(path) {
  exu <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  sample_cols <- grep("^YS\\d+_F-\\d+$", colnames(exu), value = TRUE)
  if (length(sample_cols) == 0) stop("No metabolome sample columns like 'YS0_F-1' in: ", path)

  mat <- t(as.matrix(exu[, sample_cols, drop = FALSE]))
  storage.mode(mat) <- "double"
  mat <- log2(mat + 1)

  sds <- apply(mat, 2, sd)
  keep <- is.finite(sds) & (sds > 0)
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 2) stop("Metabolome matrix has <2 varying features after filtering.")

  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
  scores <- pca$x[, 1:2, drop = FALSE]
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  list(scores = scores, var_explained = var_explained[1:2])
}

map_rdys_to_ys <- function(rdys_sample) {
  m <- regexec("^RDYS_(\\d+)-(\\d+)$", rdys_sample)
  parts <- regmatches(rdys_sample, m)[[1]]
  if (length(parts) == 0) return(NA_character_)

  group_idx <- as.integer(parts[2])
  rep_idx <- as.integer(parts[3])
  ratios <- c(0, 25, 50, 75, 100)

  if (is.na(group_idx) || group_idx < 1 || group_idx > length(ratios)) return(NA_character_)
  if (is.na(rep_idx) || rep_idx < 1) return(NA_character_)
  sprintf("YS%d_F-%d", ratios[group_idx], rep_idx)
}

extract_microbiome_pcoa <- function(path) {
  mic <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  sample_cols <- grep("^RDYS_\\d+-\\d+$", colnames(mic), value = TRUE)
  if (length(sample_cols) == 0) stop("No microbiome sample columns like 'RDYS_1-1' in: ", path)

  counts <- as.matrix(mic[, sample_cols, drop = FALSE])
  storage.mode(counts) <- "double"
  counts[is.na(counts)] <- 0

  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  if (nrow(counts) == 0) stop("Microbiome count table has no non-zero features after filtering.")

  counts_t <- t(counts) # samples x features
  sample_totals <- rowSums(counts_t)
  rel <- sweep(counts_t, 1, sample_totals, "/")
  rel[!is.finite(rel)] <- 0

  bray <- vegdist(rel, method = "bray")
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  if (any(!is.finite(pcoa$points))) {
    pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)
  }
  scores <- pcoa$points

  mapped <- vapply(rownames(scores), map_rdys_to_ys, character(1))
  if (any(is.na(mapped))) {
    bad <- rownames(scores)[is.na(mapped)]
    stop("Failed to map some RDYS sample names to YS-style: ", paste(bad, collapse = ", "))
  }
  if (anyDuplicated(mapped)) stop("Mapped microbiome sample names are not unique.")
  rownames(scores) <- mapped

  eigvals <- pcoa$eig
  pos <- eigvals > 0
  var_explained <- eigvals / sum(eigvals[pos])
  list(scores = scores, var_explained = var_explained[1:2])
}

parse_ratio_from_ys <- function(ys_sample) {
  as.integer(sub("^YS(\\d+)_.*$", "\\1", ys_sample))
}

met <- extract_metabolome_pca(exudate_path)
mic <- extract_microbiome_pcoa(microbe_path)

pca_scores <- met$scores
pcoa_scores <- mic$scores

common_samples <- intersect(rownames(pca_scores), rownames(pcoa_scores))
if (length(common_samples) < 3) stop("Too few matched samples between metabolome and microbiome: ", length(common_samples))
pca_scores <- pca_scores[common_samples, , drop = FALSE]
pcoa_scores <- pcoa_scores[common_samples, , drop = FALSE]

set.seed(123)
proc <- procrustes(pca_scores, pcoa_scores, symmetric = TRUE)
prot <- protest(pca_scores, pcoa_scores, permutations = 10000, symmetric = TRUE)

m2 <- proc$ss
p_text <- format.pval(prot$signif, digits = 3, eps = 1e-4)
subtitle <- sprintf("m^2 = %.2f, P = %s (PROTEST, 10,000 permutations)", m2, p_text)

plot_df <- data.frame(
  Sample = rownames(proc$X),
  X1 = proc$X[, 1],
  X2 = proc$X[, 2],
  Y1 = proc$Yrot[, 1],
  Y2 = proc$Yrot[, 2],
  stringsAsFactors = FALSE
)
plot_df$ArrowLength <- sqrt((plot_df$X1 - plot_df$Y1)^2 + (plot_df$X2 - plot_df$Y2)^2)
plot_df$Ratio <- vapply(plot_df$Sample, parse_ratio_from_ys, integer(1))
plot_df$Group <- factor(
  plot_df$Ratio,
  levels = c(0, 25, 50, 75, 100),
  labels = c("0%", "25%", "50%", "75%", "100%")
)

p_main <- ggplot(plot_df) +
  geom_segment(
    aes(x = X1, y = X2, xend = Y1, yend = Y2, color = Group),
    arrow = grid::arrow(length = grid::unit(0.12, "cm")),
    linewidth = 0.5,
    alpha = 0.85
  ) +
  geom_point(aes(x = X1, y = X2, fill = Group), shape = 21, size = 3.2, color = "black", stroke = 0.3) +
  geom_point(aes(x = Y1, y = Y2, fill = Group), shape = 24, size = 3.2, color = "black", stroke = 0.3) +
  coord_equal() +
  labs(
    title = "Procrustes: Metabolome PCA vs Microbiome PCoA (Bray–Curtis)",
    subtitle = subtitle,
    x = "Dimension 1",
    y = "Dimension 2",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

ggsave(file.path(out_dir, "Fig8_main.pdf"), plot = p_main, width = 7, height = 6, units = "in")
ggsave(
  file.path(out_dir, "Fig8_main.tiff"),
  plot = p_main,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

make_supplement <- TRUE
if (isTRUE(make_supplement)) {
  p_resid <- ggplot(plot_df, aes(x = Group, y = ArrowLength, color = Group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.6, alpha = 0.12) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
    labs(
      title = "Procrustes residuals (arrow length)",
      subtitle = subtitle,
      x = "Treatment",
      y = "Arrow length"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")

  ggsave(file.path(out_dir, "FigS8_procrustes_residuals.pdf"), plot = p_resid, width = 7, height = 4, units = "in")
}

message("Done. Outputs written to: ", out_dir)
