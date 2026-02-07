#!/usr/bin/env Rscript

## Supplementary｜代谢组 PCA loadings/biplot
## 数据：data/完整数据-分泌物.csv
## 预处理：TIC → log2(x+1)；PCA：prcomp(center=TRUE, scale.=TRUE)
## 产出：FigS_metabolome_PCA_loadings.pdf + PCA_loadings_table.csv

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

fig2_make_pca_loadings <- function(
  data_csv = NULL,
  out_dir = NULL,
  tic_scale_factor = 1e6,
  rsd_xlsx = NULL,
  rsd_threshold = NULL,
  top_n_arrows = 5,
  top_n_bar = 30
) {
  fig2_require_pkgs(c("ggplot2", "grid", "ggrepel"))
  if (!is.null(rsd_threshold)) {
    fig2_require_pkgs("readxl")
  }

  script_dir <- fig2_get_script_dir()
  project_root <- fig2_guess_project_root(script_dir)

  if (is.null(data_csv)) {
    data_csv <- file.path(project_root, "data", "完整数据-分泌物.csv")
  }
  if (is.null(out_dir)) {
    out_dir <- script_dir
  }

  out_pdf <- file.path(out_dir, "FigS_metabolome_PCA_loadings.pdf")
  out_table <- file.path(out_dir, "PCA_loadings_table.csv")

  d <- fig2_read_exudate_table(data_csv)
  feature_meta <- d$feature_meta
  mat_raw <- d$mat
  sample_meta <- d$sample_meta
  group_n <- table(sample_meta$Treatment)

  # --- TIC → log2(x+1) ---
  mat_tic <- fig2_tic_normalize(mat_raw, scale_factor = tic_scale_factor)
  mat_log <- log2(mat_tic + 1)

  # --- 可选：用 xlsx 过滤高 RSD 特征（需要明确阈值）---
  if (!is.null(rsd_threshold)) {
    if (is.null(rsd_xlsx)) {
      rsd_xlsx <- file.path(project_root, "data", "代谢物定性定量结果表.xlsx")
    }
    if (file.exists(rsd_xlsx)) {
      rsd_df <- readxl::read_excel(rsd_xlsx)
      rsd_names <- names(rsd_df)

      name_col <- rsd_names[match(TRUE, tolower(rsd_names) %in% tolower(c("Name", "Compound", "Metabolite")))]
      rsd_col <- grep("RSD", rsd_names, ignore.case = TRUE, value = TRUE)[1]

      if (!is.na(name_col) && !is.na(rsd_col)) {
        rsd_vec <- suppressWarnings(as.numeric(rsd_df[[rsd_col]]))
        keep_names <- rsd_df[[name_col]][is.finite(rsd_vec) & (rsd_vec <= rsd_threshold)]
        keep_names <- intersect(rownames(mat_log), keep_names)
        if (length(keep_names) > 10) {
          mat_log <- mat_log[keep_names, , drop = FALSE]
          feature_meta <- feature_meta[keep_names, , drop = FALSE]
          message("已按 RSD≤", rsd_threshold, " 过滤特征：保留 ", length(keep_names), " 个。")
        } else {
          message("RSD 过滤后特征过少（≤10），跳过过滤。")
        }
      } else {
        message("xlsx 未找到可用的 Name/RSD 列，跳过 RSD 过滤。")
      }
    } else {
      message("未找到 xlsx：", rsd_xlsx, "（跳过 RSD 过滤）")
    }
  }

  # --- PCA：去掉零方差变量 ---
  sd_vec <- apply(mat_log, 1, sd, na.rm = TRUE)
  keep <- is.finite(sd_vec) & (sd_vec > 0)
  mat_log <- mat_log[keep, , drop = FALSE]
  feature_meta <- feature_meta[keep, , drop = FALSE]

  pca <- stats::prcomp(t(mat_log), center = TRUE, scale. = TRUE)
  var_expl <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$Sample <- rownames(scores)
  scores$Treatment <- sample_meta[scores$Sample, "Treatment"]
  scores$Rep <- sample_meta[scores$Sample, "Rep"]
  treatment_levels <- fig2_default_treatment_order(sample_meta$Treatment)
  scores$Treatment <- factor(scores$Treatment, levels = treatment_levels)

  loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE])
  loadings$Name <- rownames(loadings)
  loadings$Strength_PC1_PC2 <- sqrt(loadings$PC1 ^ 2 + loadings$PC2 ^ 2)
  loadings <- loadings[order(-loadings$Strength_PC1_PC2), , drop = FALSE]
  loadings$Rank_Strength <- seq_len(nrow(loadings))

  # 导出 loadings 表（含代谢物注释列）
  meta_no_name <- feature_meta[loadings$Name, setdiff(colnames(feature_meta), "Name"), drop = FALSE]
  loadings_out <- cbind(
    Name = loadings$Name,
    meta_no_name,
    loadings[, setdiff(colnames(loadings), "Name"), drop = FALSE]
  )
  write.csv(loadings_out, out_table, row.names = FALSE)

  # --- biplot（仅标注 top_n_arrows 个 loading）---
  top_n_arrows <- min(top_n_arrows, nrow(loadings))
  arrows_df <- loadings[seq_len(top_n_arrows), , drop = FALSE]

  # 缩放箭头到与 scores 同量级
  score_rng1 <- range(scores$PC1, finite = TRUE)
  score_rng2 <- range(scores$PC2, finite = TRUE)
  load_rng1 <- range(arrows_df$PC1, finite = TRUE)
  load_rng2 <- range(arrows_df$PC2, finite = TRUE)
  denom1 <- diff(load_rng1); denom2 <- diff(load_rng2)
  if (!is.finite(denom1) || denom1 == 0) denom1 <- 1
  if (!is.finite(denom2) || denom2 == 0) denom2 <- 1
  scale_factor <- 0.75 * min(diff(score_rng1) / denom1, diff(score_rng2) / denom2)
  arrows_df$PC1_scaled <- arrows_df$PC1 * scale_factor
  arrows_df$PC2_scaled <- arrows_df$PC2 * scale_factor

  base_treatment_cols <- c("#1B9E77", "#66A61E", "#E6AB02", "#D95F02", "#E7298A")
  if (length(treatment_levels) <= length(base_treatment_cols)) {
    treatment_cols <- setNames(base_treatment_cols[seq_along(treatment_levels)], treatment_levels)
  } else {
    fig2_require_pkgs(c("circlize"))
    extra_cols <- circlize::rand_color(length(treatment_levels) - length(base_treatment_cols), luminosity = "bright")
    treatment_cols <- setNames(c(base_treatment_cols, extra_cols), treatment_levels)
  }

  pc1_lab <- sprintf("PC1 (%.1f%%)", var_expl[1] * 100)
  pc2_lab <- sprintf("PC2 (%.1f%%)", var_expl[2] * 100)

  p_biplot <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = Treatment)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
    ggplot2::geom_segment(
      data = arrows_df,
      ggplot2::aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed", angle = 20),
      linewidth = 0.5,
      color = "grey20",
      alpha = 0.7
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = treatment_cols) +
    ggrepel::geom_text_repel(
      data = arrows_df,
      ggplot2::aes(x = PC1_scaled, y = PC2_scaled, label = Name),
      inherit.aes = FALSE,
      size = 2.3,
      color = "grey20",
      min.segment.length = 0,
      segment.color = "grey50",
      segment.size = 0.3,
      box.padding = 0.3,
      max.overlaps = Inf
    ) +
    ggplot2::labs(
      title = "Metabolome PCA biplot",
      subtitle = "TIC → log2(x+1); prcomp(center=TRUE, scale.=TRUE)",
      x = pc1_lab,
      y = pc2_lab
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.margin = grid::unit(c(1, 1, 1, 1), "cm") 
    )

  # Clean IUPAC names for display
  clean_iupac_names <- function(x) {
    cleaned <- sub("^\\([0-9a-zA-Z,]+\\)-", "", x)
    dupes <- cleaned[duplicated(cleaned) | duplicated(cleaned, fromLast = TRUE)]
    if (length(dupes) > 0) {
      idx <- which(cleaned %in% dupes)
      cleaned[idx] <- x[idx]
    }
    cleaned
  }

  # Wrap long names (force wrap on special chars or length)
  wrap_chemical_name <- function(x, width = 30) {
    x <- gsub("([-,])", "\\1 ", x)
    x <- vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
    return(x)
  }

  # --- Top loadings strength Lollipop Chart ---
  top_n_bar <- min(20, nrow(loadings)) # Force reduce to 20 for cleaner layout per user feedback
  bar_df <- loadings[seq_len(top_n_bar), c("Name", "Strength_PC1_PC2"), drop = FALSE]
  
  # 1. Clean names
  bar_df$Name_Display <- clean_iupac_names(bar_df$Name)
  
  # 2. Wrap names
  bar_df$Name_Wrapped <- wrap_chemical_name(bar_df$Name_Display, width = 40)
  bar_df$Name_Wrapped <- factor(bar_df$Name_Wrapped, levels = rev(bar_df$Name_Wrapped))

  # Dynamic baseline to "zoom in"
  min_val <- min(bar_df$Strength_PC1_PC2)
  max_val <- max(bar_df$Strength_PC1_PC2)
  # Start axis at 90% of the min value to exaggerate differences
  lower_limit <- min_val * 0.9 

  p_bar <- ggplot2::ggplot(bar_df, ggplot2::aes(x = Name_Wrapped, y = Strength_PC1_PC2, color = Strength_PC1_PC2)) +
    # Lollipop stem starts from lower_limit
    ggplot2::geom_segment(ggplot2::aes(x = Name_Wrapped, xend = Name_Wrapped, y = lower_limit, yend = Strength_PC1_PC2), 
                          linewidth = 0.8) +
    # Lollipop head
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_viridis_c(option = "C", direction = -1, guide = "none") + # Heatmap style color
    ggplot2::scale_y_continuous(limits = c(lower_limit, NA), oob = scales::squish, expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Top ", top_n_bar, " PCA loadings (PC1/PC2 strength)"),
      subtitle = "Zoomed in view (axis does not start at 0)",
      x = NULL,
      y = "sqrt(PC1^2 + PC2^2)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        # Add some vertical spacing between axis labels if possible? 
        # Reducing top_n to 20 is the main fix.
        plot.margin = grid::unit(c(1, 1, 1, 0), "cm")
    )


  # --- 输出 PDF（左右两栏）---
  message("每组样本数：", paste(names(group_n), as.integer(group_n), sep = "=", collapse = ", "))
  message("输出：", out_pdf)
  # Increase width slightly more to be safe? 16 is already huge.
  grDevices::pdf(out_pdf, width = 16, height = 9, useDingbats = FALSE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(0.4, 0.6), "npc"))))
  print(p_biplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_bar, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()

  invisible(list(out_pdf = out_pdf, out_table = out_table))
}

if (sys.nframe() == 0) {
  fig2_make_pca_loadings()
}
