# fill_module_column.R
# 目的：复现 figure4 热图的聚类逻辑，填补 Supplementary_Table_S3 中缺失的 Module 列
# 依赖：无额外包（仅使用 base R）

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

# -------------------- paths --------------------
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."))
data_path <- file.path(project_root, "data", "完整数据-微生物.csv")
supp_csv <- file.path(script_dir, "Supplementary_Table_S3_Top27_Phyla_Fig4.csv")

stopifnot(file.exists(data_path))
stopifnot(file.exists(supp_csv))

# -------------------- 复现 figure4 数据处理流程 --------------------
top_n_main <- 27
treatment_labels <- c("RS0", "RS25", "RS50", "RS75", "RS100")

if (requireNamespace("data.table", quietly = TRUE)) {
  df <- data.table::fread(data_path, data.table = FALSE)
} else {
  df <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
}

sample_cols <- grep("^RDYS_", names(df), value = TRUE)
df[sample_cols] <- lapply(df[sample_cols], as.numeric)
phylum <- vapply(df$taxonomy, extract_phylum, character(1))

counts_mat <- rowsum(as.matrix(df[sample_cols]), group = phylum, na.rm = TRUE, reorder = FALSE)
counts_mat <- counts_mat[rownames(counts_mat) != "", , drop = FALSE]
col_totals <- colSums(counts_mat)
rel_mat <- sweep(counts_mat, 2, col_totals, "/")
rel_mat[is.na(rel_mat)] <- 0

# sample metadata + ordering
sample_meta <- data.frame(
  Sample = colnames(rel_mat),
  Treatment_num = vapply(colnames(rel_mat), parse_treatment_num, integer(1)),
  Replicate = vapply(colnames(rel_mat), parse_replicate, integer(1)),
  stringsAsFactors = FALSE
)
sample_meta$Treatment <- factor(
  sample_meta$Treatment_num,
  levels = seq_along(treatment_labels),
  labels = treatment_labels
)
ord <- order(sample_meta$Treatment_num, sample_meta$Replicate)
rel_mat <- rel_mat[, ord, drop = FALSE]
sample_meta <- sample_meta[ord, , drop = FALSE]

# treatment-level mean
treat_levels <- levels(sample_meta$Treatment)
rel_mean_mat <- sapply(treat_levels, function(tr) {
  cols <- sample_meta$Treatment == tr
  rowMeans(rel_mat[, cols, drop = FALSE])
})
colnames(rel_mean_mat) <- treat_levels

# top 27 phyla
phylum_mean <- rowMeans(rel_mat)
phylum_mean <- phylum_mean[phylum_mean > 0]
phylum_rank <- sort(phylum_mean, decreasing = TRUE)
top_n_main <- min(top_n_main, length(phylum_rank))
key_phyla <- names(phylum_rank)[seq_len(top_n_main)]

# -------------------- 聚类得到 Module 分配 --------------------
mat_rel <- rel_mean_mat[key_phyla, , drop = FALSE]
cl_main <- cluster_to_modules(mat_rel, k = 2)

# 构建 phylum -> module 映射
module_map <- setNames(as.character(cl_main$module), key_phyla)

cat("聚类结果：\n")
for (m in sort(unique(module_map))) {
  cat(sprintf("  %s: %s\n", m, paste(names(module_map[module_map == m]), collapse = ", ")))
}

# -------------------- 更新 CSV --------------------
supp_tbl <- read.csv(supp_csv, stringsAsFactors = FALSE, check.names = FALSE)

supp_tbl$Module_in_Fig4_annotation <- module_map[supp_tbl$Phylum]

# 检查是否所有 Phylum 都匹配到了 Module
missing <- is.na(supp_tbl$Module_in_Fig4_annotation)
if (any(missing)) {
  warning(sprintf("以下 Phylum 未匹配到 Module: %s",
                  paste(supp_tbl$Phylum[missing], collapse = ", ")))
}

write.csv(supp_tbl, supp_csv, row.names = FALSE)
cat(sprintf("\n已更新: %s\n", supp_csv))
