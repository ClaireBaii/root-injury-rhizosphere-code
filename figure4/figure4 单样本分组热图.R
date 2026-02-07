# Figure 4: Phylum-level heatmap (ComplexHeatmap)
# R version: 4.4.3
# Data: ../data/完整数据-微生物.csv
# Outputs:
#   - Fig4_main.pdf + Fig4_main.tiff (main text, simplified)
#   - (optional) FigS4_full_phylum_heatmap.pdf (supplement, full heatmap)

# -------------------- config --------------------
top_n_main <- 27
make_full_supp <- TRUE

# RDYS_1..5 ↔ treatment labels (RS = Rhizosphere Soil %)
treatment_labels <- c("RS0", "RS25", "RS50", "RS75", "RS100")
treatment_colors <- c(
  RS0 = "#1B9E77",
  RS25 = "#D95F02",
  RS50 = "#7570B3",
  RS75 = "#E7298A",
  RS100 = "#66A61E"
)

module_colors <- c("Module 1" = "#E41A1C", "Module 2" = "#377EB8")

font_sizes <- list(
  row_names_main = 12,
  col_names_main = 13,
  row_names_full = 10,
  col_names_full = 10,
  legend_title = 12,
  legend_labels = 11
)

# -------------------- helpers --------------------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Missing package '%s'. Please install it first.", pkg), call. = FALSE)
  }
}

extract_phylum <- function(taxonomy) {
  if (is.na(taxonomy) || !nzchar(taxonomy)) return("Unclassified")
  m <- regmatches(taxonomy, regexpr("(?:^|;\\s*)p__[^;]+", taxonomy, perl = TRUE))
  if (length(m) == 0 || !nzchar(m)) return("Unclassified")
  phylum <- sub("^.*p__", "", m)
  phylum <- trimws(phylum)
  if (!nzchar(phylum)) return("Unclassified")
  phylum
}

parse_treatment_num <- function(sample_id) {
  m <- regmatches(sample_id, regexpr("^RDYS_\\d+", sample_id, perl = TRUE))
  if (length(m) == 0 || !nzchar(m)) return(NA_integer_)
  suppressWarnings(as.integer(sub("^RDYS_", "", m)))
}

parse_replicate <- function(sample_id) {
  m <- regmatches(sample_id, regexpr("-\\d+$", sample_id, perl = TRUE))
  if (length(m) == 0 || !nzchar(m)) return(NA_integer_)
  suppressWarnings(as.integer(sub("^-", "", m)))
}

cluster_to_modules <- function(mat, k = 2) {
  # cluster based on row-wise Z-score (pattern across columns)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  hc <- hclust(dist(z), method = "ward.D2")
  module_id <- cutree(hc, k = k)
  module <- factor(
    module_id,
    levels = sort(unique(module_id)),
    labels = paste0("Module ", seq_along(sort(unique(module_id))))
  )
  list(hc = hc, module = module)
}

# -------------------- packages --------------------
require_pkg("ComplexHeatmap")
require_pkg("circlize")
require_pkg("grid")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# -------------------- IO paths --------------------
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."))
data_path <- file.path(project_root, "data", "完整数据-微生物.csv")
out_dir <- script_dir

if (!file.exists(data_path)) {
  stop(sprintf("Data file not found: %s", data_path), call. = FALSE)
}

# -------------------- read + aggregate to phylum --------------------
if (requireNamespace("data.table", quietly = TRUE)) {
  df <- data.table::fread(data_path, data.table = FALSE)
} else {
  df <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
}

if (!("taxonomy" %in% names(df))) {
  stop("Missing column 'taxonomy' in the input file.", call. = FALSE)
}

sample_cols <- grep("^RDYS_", names(df), value = TRUE)
if (length(sample_cols) == 0) {
  stop("No sample columns detected (expected names starting with 'RDYS_').", call. = FALSE)
}

df[sample_cols] <- lapply(df[sample_cols], as.numeric)
phylum <- vapply(df$taxonomy, extract_phylum, character(1))

counts_mat <- rowsum(as.matrix(df[sample_cols]), group = phylum, na.rm = TRUE, reorder = FALSE)
counts_mat <- counts_mat[rownames(counts_mat) != "", , drop = FALSE]

col_totals <- colSums(counts_mat)
if (any(col_totals == 0)) {
  stop("At least one sample has zero total counts after aggregation.", call. = FALSE)
}

rel_mat <- sweep(counts_mat, 2, col_totals, "/")
rel_mat[is.na(rel_mat)] <- 0

# -------------------- sample metadata + ordering --------------------
sample_meta <- data.frame(
  Sample = colnames(rel_mat),
  Treatment_num = vapply(colnames(rel_mat), parse_treatment_num, integer(1)),
  Replicate = vapply(colnames(rel_mat), parse_replicate, integer(1)),
  stringsAsFactors = FALSE
)

if (any(is.na(sample_meta$Treatment_num))) {
  stop("Failed to parse treatment index from sample names (expected 'RDYS_1-1' style).", call. = FALSE)
}

sample_meta$Treatment <- factor(
  sample_meta$Treatment_num,
  levels = seq_along(treatment_labels),
  labels = treatment_labels
)

ord <- order(sample_meta$Treatment_num, sample_meta$Replicate)
rel_mat <- rel_mat[, ord, drop = FALSE]
sample_meta <- sample_meta[ord, , drop = FALSE]

# treatment-level mean matrix (simplified, n=3 per treatment)
treat_levels <- levels(sample_meta$Treatment)
rel_mean_mat <- sapply(treat_levels, function(tr) {
  cols <- sample_meta$Treatment == tr
  rowMeans(rel_mat[, cols, drop = FALSE])
})
colnames(rel_mean_mat) <- treat_levels

# -------------------- pick key phyla (for main text) --------------------
phylum_mean <- rowMeans(rel_mat)
phylum_mean <- phylum_mean[phylum_mean > 0]
phylum_rank <- sort(phylum_mean, decreasing = TRUE)
top_n_main <- min(top_n_main, length(phylum_rank))
key_phyla <- names(phylum_rank)[seq_len(top_n_main)]

key_tbl <- data.frame(
  Phylum = key_phyla,
  Mean_rel_abundance = unname(phylum_rank[key_phyla]),
  Mean_percent = unname(phylum_rank[key_phyla]) * 100,
  stringsAsFactors = FALSE
)
write.csv(key_tbl, file.path(out_dir, "Fig4_key_phyla.csv"), row.names = FALSE)
writeLines(paste(key_phyla, collapse = ", "), file.path(out_dir, "Fig4_key_phyla.txt"))

# -------------------- Figure 4 main: key phyla × treatment means (Z-score) ----------
mat_rel <- rel_mean_mat[key_phyla, , drop = FALSE]
# Row-wise Z-score standardization
mat_main <- t(scale(t(mat_rel)))
mat_main[is.na(mat_main)] <- 0

cl_main <- cluster_to_modules(mat_rel, k = 2)  # cluster on original values
mat_main <- mat_main[cl_main$hc$order, , drop = FALSE]
module_main <- cl_main$module[cl_main$hc$order]

# Symmetric diverging color scale: blue-white-red
col_fun_main <- colorRamp2(
  c(-1.5, 0, 1.5),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

ha_col_main <- HeatmapAnnotation(
  Treatment = factor(colnames(mat_main), levels = treatment_labels),
  col = list(Treatment = treatment_colors),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = font_sizes$legend_title),
    labels_gp = gpar(fontsize = font_sizes$legend_labels)
  )
)

ha_row_main <- rowAnnotation(
  Module = module_main,
  col = list(Module = module_colors),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = font_sizes$legend_title),
    labels_gp = gpar(fontsize = font_sizes$legend_labels)
  )
)

ht_main <- Heatmap(
  mat_main,
  name = "Z-score",
  col = col_fun_main,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = module_main,
  row_gap = unit(2, "mm"),
  column_title = "Treatment (n=3)",
  column_title_gp = gpar(fontsize = font_sizes$col_names_main + 1, fontface = "bold"),
  top_annotation = ha_col_main,
  left_annotation = ha_row_main,
  row_names_gp = gpar(fontsize = font_sizes$row_names_main),
  column_names_gp = gpar(fontsize = font_sizes$col_names_main),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = font_sizes$legend_title),
    labels_gp = gpar(fontsize = font_sizes$legend_labels)
  )
)

pdf(file.path(out_dir, "Fig4_main.pdf"), width = 7, height = 9, useDingbats = FALSE)
draw(ht_main, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

tiff(
  file.path(out_dir, "Fig4_main.tiff"),
  width = 7,
  height = 9,
  units = "in",
  res = 600,
  compression = "lzw"
)
draw(ht_main, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# -------------------- Figure S4 (optional): full phylum heatmap (all samples) --------------------
if (isTRUE(make_full_supp)) {
  nonzero_phyla <- names(phylum_mean)
  mat_rel_full <- rel_mat[nonzero_phyla, , drop = FALSE]
  mean_rel_full <- rel_mean_mat[nonzero_phyla, , drop = FALSE]

  # Row-wise Z-score for display matrix
  mat_full <- t(scale(t(mat_rel_full)))
  mat_full[is.na(mat_full)] <- 0

  cl_full <- cluster_to_modules(mean_rel_full, k = 2)  # cluster on original
  mat_full <- mat_full[cl_full$hc$order, , drop = FALSE]
  module_full <- cl_full$module[cl_full$hc$order]

  # Symmetric diverging color scale: blue-white-red
  col_fun_full <- colorRamp2(
    c(-1.5, 0, 1.5),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )

  ha_col_full <- HeatmapAnnotation(
    Treatment = sample_meta$Treatment,
    col = list(Treatment = treatment_colors),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = font_sizes$legend_title),
      labels_gp = gpar(fontsize = font_sizes$legend_labels)
    )
  )

  ha_row_full <- rowAnnotation(
    Module = module_full,
    col = list(Module = module_colors),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = font_sizes$legend_title),
      labels_gp = gpar(fontsize = font_sizes$legend_labels)
    )
  )

  ht_full <- Heatmap(
    mat_full,
    name = "Z-score",
    col = col_fun_full,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_split = module_full,
    row_gap = unit(1.5, "mm"),
    column_split = sample_meta$Treatment,
    column_gap = unit(1.5, "mm"),
    top_annotation = ha_col_full,
    left_annotation = ha_row_full,
    row_names_gp = gpar(fontsize = font_sizes$row_names_full),
    column_names_gp = gpar(fontsize = font_sizes$col_names_full),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = font_sizes$legend_title),
      labels_gp = gpar(fontsize = font_sizes$legend_labels)
    )
  )

  pdf(file.path(out_dir, "FigS4_full_phylum_heatmap.pdf"), width = 11, height = 10, useDingbats = FALSE)
  draw(ht_full, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}
