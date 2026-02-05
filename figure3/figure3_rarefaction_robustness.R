#!/usr/bin/env Rscript

# FigS｜rarefied vs non-rarefied 的 PCoA / PERMANOVA 稳健性对照
# R 4.4.3
#
# 输入：
# - data/完整数据-微生物.csv
#
# 输出（写入 figure3/）：
# - FigS_rarefaction_robustness.pdf + FigS_rarefaction_robustness.tiff
# - permanova_robustness.csv

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

get_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

require_cran <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message(sprintf("Installing missing package: %s", p))
      install.packages(p, repos = "https://cran.r-project.org")
    }
  }
}

infer_group <- function(sample_id) {
  m <- regexec("^RDYS_(\\d+)-", sample_id)
  mm <- regmatches(sample_id, m)[[1]]
  if (length(mm) >= 2) {
    grp_num <- as.integer(mm[2])
    return(switch(
      as.character(grp_num),
      "1" = "Control",
      "2" = "Treat1",
      "3" = "Treat2",
      "4" = "Treat3",
      "5" = "Treat4",
      paste0("Group", grp_num)
    ))
  }
  if (grepl("control", sample_id, ignore.case = TRUE)) return("Control")
  m2 <- regexec("treat\\s*([0-9]+)", sample_id, ignore.case = TRUE)
  mm2 <- regmatches(sample_id, m2)[[1]]
  if (length(mm2) >= 2) return(paste0("Treat", mm2[2]))
  "Unknown"
}

compute_pcoa <- function(dist_obj, meta, permutations, method_label) {
  pcoa <- cmdscale(dist_obj, k = 2, eig = TRUE)
  coords <- as.data.frame(pcoa$points)
  colnames(coords) <- c("PC1", "PC2")
  coords$Sample <- rownames(coords)
  coords <- coords %>% left_join(meta, by = "Sample")

  eig <- pcoa$eig
  eig_pos <- eig[eig > 0]
  var_expl <- eig_pos / sum(eig_pos) * 100
  pc1_pct <- ifelse(length(var_expl) >= 1, var_expl[1], NA_real_)
  pc2_pct <- ifelse(length(var_expl) >= 2, var_expl[2], NA_real_)

  adon <- vegan::adonis2(dist_obj ~ Group, data = meta, permutations = permutations)
  adon_df <- as.data.frame(adon)
  grp_row <- adon_df[rownames(adon_df) == "Group", , drop = FALSE]
  r2 <- as.numeric(grp_row$R2[1])
  p <- as.numeric(grp_row$`Pr(>F)`[1])

  list(
    coords = coords,
    pc1_pct = pc1_pct,
    pc2_pct = pc2_pct,
    permanova_r2 = r2,
    permanova_p = p,
    method_label = method_label
  )
}

repo_root <- get_repo_root()
in_csv <- file.path(repo_root, "data", "完整数据-微生物.csv")
out_dir <- file.path(repo_root, "figure3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

permutations <- 9999
rarefy_seed <- 1

require_cran(c("vegan", "ggplot2", "dplyr", "tidyr", "patchwork"))
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

if (!file.exists(in_csv)) stop(sprintf("Input not found: %s", in_csv))

raw <- read.csv(
  in_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

name_lower <- tolower(names(raw))
otu_id_col <- which(name_lower %in% c("otu id", "otu_id", "feature id", "feature_id", "asv id", "asv_id"))
otu_id_col <- if (length(otu_id_col) > 0) otu_id_col[1] else 1
tax_col <- which(name_lower == "taxonomy")
tax_col <- if (length(tax_col) > 0) tax_col[1] else ncol(raw)

raw <- raw[!is.na(raw[[otu_id_col]]) & raw[[otu_id_col]] != "", , drop = FALSE]
raw <- raw[!duplicated(raw[[otu_id_col]]), , drop = FALSE]

sample_cols <- setdiff(seq_len(ncol(raw)), c(otu_id_col, tax_col))
sample_ids <- names(raw)[sample_cols]

otu_mat <- as.matrix(raw[, sample_cols, drop = FALSE])
mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0
otu_mat <- round(otu_mat)
rownames(otu_mat) <- raw[[otu_id_col]]
otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]

meta <- data.frame(
  Sample = sample_ids,
  Group = vapply(sample_ids, infer_group, character(1)),
  stringsAsFactors = FALSE
)
group_order <- c("Control", "Treat1", "Treat2", "Treat3", "Treat4")
meta$Group <- factor(meta$Group, levels = unique(c(group_order, sort(unique(meta$Group)))))

group_colors <- c(
  Control = "#F8766D",
  Treat1 = "#7CAE00",
  Treat2 = "#00BFC4",
  Treat3 = "#00B0F6",
  Treat4 = "#C77CFF"
)

# -------- non-rarefied：相对丰度 Bray --------
otu_rel <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")
otu_rel[is.na(otu_rel)] <- 0
dist_nonref <- vegan::vegdist(t(otu_rel), method = "bray")
res_nonref <- compute_pcoa(dist_nonref, meta, permutations, method_label = "Non-rarefied (relative abundance)")

# -------- rarefied：等深度稀释后 Bray --------
min_depth <- min(rowSums(t(otu_mat)))
set.seed(rarefy_seed)
otu_raref <- vegan::rrarefy(t(otu_mat), sample = min_depth)
dist_raref <- vegan::vegdist(otu_raref, method = "bray")
res_raref <- compute_pcoa(dist_raref, meta, permutations, method_label = sprintf("Rarefied (depth=%d)", min_depth))

robust_tbl <- data.frame(
  method = c("non_rarefied", "rarefied"),
  label = c(res_nonref$method_label, res_raref$method_label),
  R2 = c(res_nonref$permanova_r2, res_raref$permanova_r2),
  p_value = c(res_nonref$permanova_p, res_raref$permanova_p),
  permutations = permutations,
  rarefaction_depth = c(NA_integer_, min_depth),
  distance = "Bray-Curtis",
  stringsAsFactors = FALSE
)
write.csv(robust_tbl, file.path(out_dir, "permanova_robustness.csv"), row.names = FALSE)

permanova_text <- function(r2, p, perm) sprintf("PERMANOVA: R² = %.3f, P = %.3g (perm=%d)", r2, p, perm)

make_plot <- function(res, xlim = NULL, ylim = NULL) {
  p <- ggplot(res$coords, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3.0, alpha = 0.95) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(
      title = res$method_label,
      subtitle = permanova_text(res$permanova_r2, res$permanova_p, permutations),
      x = sprintf("PC1 (%.1f%%)", res$pc1_pct),
      y = sprintf("PC2 (%.1f%%)", res$pc2_pct),
      color = "Group"
    ) +
    coord_equal() +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  if (!is.null(xlim) && !is.null(ylim)) {
    p <- p + xlim(xlim) + ylim(ylim)
  }
  p
}

# 统一坐标范围，方便对照
all_coords <- bind_rows(
  res_nonref$coords %>% transmute(PC1, PC2),
  res_raref$coords %>% transmute(PC1, PC2)
)
x_pad <- diff(range(all_coords$PC1)) * 0.06
y_pad <- diff(range(all_coords$PC2)) * 0.06
xlim <- range(all_coords$PC1) + c(-x_pad, x_pad)
ylim <- range(all_coords$PC2) + c(-y_pad, y_pad)

p1 <- make_plot(res_nonref, xlim = xlim, ylim = ylim) + theme(legend.position = "none")
p2 <- make_plot(res_raref, xlim = xlim, ylim = ylim)

fig <- p1 + p2 + patchwork::plot_layout(widths = c(1, 1))

pdf_path <- file.path(out_dir, "FigS_rarefaction_robustness.pdf")
tiff_path <- file.path(out_dir, "FigS_rarefaction_robustness.tiff")

ggsave(pdf_path, fig, width = 12, height = 5.2, units = "in")
ggsave(
  tiff_path,
  fig,
  width = 12,
  height = 5.2,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

message("Done. Outputs written to: ", out_dir)

