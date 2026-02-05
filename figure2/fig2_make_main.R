#!/usr/bin/env Rscript

## Figure 2（代谢物热图 + Top响应）
## 需求：TIC → log2(x+1) → Z-score；Heatmap 用 ComplexHeatmap；
##       Top50 按 max(|log2FC|)（各处理组 vs YS0）排序。
## 产出：Fig2_main.pdf（矢量）+ Fig2_main.tiff（600 dpi）+ Fig2_data_top50.csv（可选）

fig2_local_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

source(file.path(fig2_local_script_dir(), "fig2_helpers.R"))

fig2_make_main <- function(
  data_csv = NULL,
  out_dir = NULL,
  top_n = 50,
  tic_scale_factor = 1e6,
  control_treatment = "YS0"
) {
  fig2_require_pkgs(c("ComplexHeatmap", "circlize", "grid"))

  script_dir <- fig2_get_script_dir()
  project_root <- fig2_guess_project_root(script_dir)

  if (is.null(data_csv)) {
    data_csv <- file.path(project_root, "data", "完整数据-分泌物.csv")
  }
  if (is.null(out_dir)) {
    out_dir <- script_dir
  }

  out_pdf <- file.path(out_dir, "Fig2_main.pdf")
  out_tiff <- file.path(out_dir, "Fig2_main.tiff")
  out_top50 <- file.path(out_dir, "Fig2_data_top50.csv")

  d <- fig2_read_exudate_table(data_csv)
  feature_meta <- d$feature_meta
  mat_raw <- d$mat
  sample_meta <- d$sample_meta
  group_n <- table(sample_meta$Treatment)

  # --- 1) TIC 归一化 ---
  mat_tic <- fig2_tic_normalize(mat_raw, scale_factor = tic_scale_factor)

  # --- 2) 组均值（n=3 trees per treatment）---
  treatment_order <- fig2_default_treatment_order(sample_meta$Treatment)
  mat_mean <- fig2_group_means(mat_tic, sample_meta, treatment_order = treatment_order)

  if (!(control_treatment %in% colnames(mat_mean))) {
    stop("未找到对照组列：", control_treatment, call. = FALSE)
  }
  non_control <- setdiff(colnames(mat_mean), control_treatment)
  if (length(non_control) == 0) {
    stop("除对照外未找到其他处理组。", call. = FALSE)
  }

  # --- 3) log2FC（各处理 vs 对照；基于 TIC 后均值强度）---
  log2fc_mat <- log2((mat_mean[, non_control, drop = FALSE] + 1) / (mat_mean[, control_treatment] + 1))
  abs_log2fc_mat <- abs(log2fc_mat)

  max_abs_log2fc <- apply(abs_log2fc_mat, 1, max, na.rm = TRUE)
  which_max <- apply(abs_log2fc_mat, 1, which.max)
  max_treatment <- non_control[which_max]
  max_log2fc <- log2fc_mat[cbind(seq_len(nrow(log2fc_mat)), which_max)]

  names(max_abs_log2fc) <- rownames(log2fc_mat)
  names(max_treatment) <- rownames(log2fc_mat)
  names(max_log2fc) <- rownames(log2fc_mat)

  # --- 4) TopN 代谢物 ---
  valid <- is.finite(max_abs_log2fc) & !is.na(max_abs_log2fc)

  # --- 5) Heatmap 矩阵：TIC → log2(x+1) → Z-score（按代谢物行标准化）---
  mat_log <- log2(mat_mean + 1)
  mat_z <- t(scale(t(mat_log)))

  # 过滤标准化失败行（sd=0）
  keep <- rownames(mat_z)[apply(mat_z, 1, function(x) all(is.finite(x)))]
  mat_z <- mat_z[keep, , drop = FALSE]
  feature_meta <- feature_meta[keep, , drop = FALSE]

  valid2 <- valid
  valid2[!names(valid2) %in% keep] <- FALSE
  ord2 <- order(max_abs_log2fc[valid2], decreasing = TRUE)
  top_names <- names(max_abs_log2fc[valid2])[ord2][seq_len(min(top_n, sum(valid2)))]
  if (length(top_names) < 2) {
    stop("Top 代谢物数量不足（<2）。请检查数据/预处理。", call. = FALSE)
  }

  mat_z_top <- mat_z[top_names, , drop = FALSE]
  top_meta <- feature_meta[top_names, , drop = FALSE]

  # --- 输出 Top50 数据表（可用于图注/补充）---
  top_out <- cbind(
    top_meta,
    max_abs_log2FC = max_abs_log2fc[top_names],
    max_Treatment = max_treatment[top_names],
    max_log2FC = max_log2fc[top_names]
  )
  write.csv(top_out, out_top50, row.names = FALSE)

  # --- 6) 颜色与注释 ---
  set.seed(1)
  superclass_levels <- sort(unique(top_meta$Super_Class))
  superclass_cols <- setNames(
    circlize::rand_color(length(superclass_levels), luminosity = "bright"),
    superclass_levels
  )

  treatment_levels <- fig2_default_treatment_order(colnames(mat_z_top))
  base_treatment_cols <- c("#1B9E77", "#66A61E", "#E6AB02", "#D95F02", "#E7298A")
  if (length(treatment_levels) <= length(base_treatment_cols)) {
    treatment_cols <- setNames(base_treatment_cols[seq_along(treatment_levels)], treatment_levels)
  } else {
    extra_cols <- circlize::rand_color(length(treatment_levels) - length(base_treatment_cols), luminosity = "bright")
    treatment_cols <- setNames(c(base_treatment_cols, extra_cols), treatment_levels)
  }

  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
  bar_fill <- ifelse(max_log2fc[top_names] >= 0, "#B2182B", "#2166AC")

  left_anno <- ComplexHeatmap::rowAnnotation(
    Super_Class = top_meta$Super_Class,
    col = list(Super_Class = superclass_cols),
    annotation_name_gp = grid::gpar(fontsize = 9)
  )

  right_anno <- ComplexHeatmap::rowAnnotation(
    Max_Treatment = factor(max_treatment[top_names], levels = setdiff(treatment_levels, control_treatment)),
    `Max log2FC` = ComplexHeatmap::anno_barplot(
      max_log2fc[top_names],
      baseline = 0,
      gp = grid::gpar(fill = bar_fill, col = NA),
      border = FALSE,
      axis_param = list(side = "top", gp = grid::gpar(fontsize = 8))
    ),
    col = list(Max_Treatment = treatment_cols),
    annotation_name_gp = grid::gpar(fontsize = 9)
  )

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    Treatment = factor(colnames(mat_z_top), levels = treatment_levels),
    col = list(Treatment = treatment_cols),
    annotation_name_gp = grid::gpar(fontsize = 9)
  )

  ht <- ComplexHeatmap::Heatmap(
    mat_z_top,
    name = "Z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    column_dend_height = grid::unit(18, "mm"),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 9),
    top_annotation = top_anno,
    left_annotation = left_anno,
    right_annotation = right_anno,
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9)
    ),
    use_raster = FALSE
  )

  # --- 7) 输出 PDF / TIFF ---
  message("Figure 2 预处理：TIC → log2(x+1) → Z-score（行）")
  message("每组样本数：", paste(names(group_n), as.integer(group_n), sep = "=", collapse = ", "),
          "（按 Treatment 取均值绘图）")
  message("输出：", out_pdf)
  grDevices::pdf(out_pdf, width = 11, height = 8, useDingbats = FALSE)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()

  message("输出：", out_tiff)
  grDevices::tiff(out_tiff, width = 11, height = 8, units = "in", res = 600, compression = "lzw")
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()

  invisible(list(
    out_pdf = out_pdf,
    out_tiff = out_tiff,
    out_top50 = out_top50
  ))
}

if (sys.nframe() == 0) {
  fig2_make_main()
}
