fig2_require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "缺少 R 包：", paste(missing, collapse = ", "),
      "\n请先安装后再运行（示例）：",
      "\n  install.packages(c(\"ggplot2\", \"readxl\"))",
      "\n  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")",
      "\n  BiocManager::install(c(\"ComplexHeatmap\", \"circlize\"))",
      call. = FALSE
    )
  }
}

fig2_get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

fig2_guess_project_root <- function(script_dir) {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
}

fig2_read_exudate_table <- function(csv_path) {
  if (!file.exists(csv_path)) {
    stop("找不到数据文件：", csv_path, call. = FALSE)
  }

  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)

  required_cols <- c("Name", "Super_Class")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("数据缺少必要列：", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  sample_cols <- grep("^YS[0-9]+_F-[0-9]+$", names(df), value = TRUE)
  if (length(sample_cols) == 0) {
    stop("未找到样本列（期望格式如 YS0_F-1）。请检查：", csv_path, call. = FALSE)
  }

  # --- feature meta ---
  feature_meta <- df[, setdiff(names(df), sample_cols), drop = FALSE]
  rownames(feature_meta) <- feature_meta$Name
  if ("Super_Class" %in% names(feature_meta)) {
    feature_meta$Super_Class[is.na(feature_meta$Super_Class) | feature_meta$Super_Class == ""] <- "Unknown"
  }

  # --- sample matrix ---
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  mat <- apply(mat, 2, function(x) as.numeric(x))
  rownames(mat) <- df$Name
  colnames(mat) <- sample_cols
  mat[is.na(mat)] <- 0

  # --- sample meta ---
  treatment <- sub("(_F-[0-9]+)$", "", sample_cols)
  rep_id <- sub("^.*_F-", "F-", sample_cols)
  treatment_num <- suppressWarnings(as.numeric(sub("^YS", "", treatment)))

  sample_meta <- data.frame(
    Sample = sample_cols,
    Treatment = treatment,
    Treatment_num = treatment_num,
    Rep = rep_id,
    stringsAsFactors = FALSE,
    row.names = sample_cols
  )

  list(feature_meta = feature_meta, mat = mat, sample_meta = sample_meta)
}

fig2_default_treatment_order <- function(treatments) {
  uniq <- unique(treatments)
  nums <- suppressWarnings(as.numeric(sub("^YS", "", uniq)))
  if (all(is.finite(nums))) {
    uniq[order(nums)]
  } else {
    sort(uniq)
  }
}

fig2_tic_normalize <- function(mat, scale_factor = 1e6) {
  tic <- colSums(mat, na.rm = TRUE)
  if (any(!is.finite(tic)) || any(tic <= 0)) {
    stop("TIC 归一化失败：存在 TIC<=0 或非有限值。", call. = FALSE)
  }
  sweep(mat, 2, tic, FUN = "/") * scale_factor
}

fig2_group_means <- function(mat, sample_meta, treatment_order = NULL) {
  treatments <- sample_meta$Treatment
  names(treatments) <- rownames(sample_meta)

  if (is.null(treatment_order)) {
    treatment_order <- fig2_default_treatment_order(treatments)
  }

  means <- sapply(treatment_order, function(trt) {
    cols <- names(treatments)[treatments == trt]
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  colnames(means) <- treatment_order
  means
}
