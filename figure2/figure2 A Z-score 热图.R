#!/usr/bin/env Rscript

# Figure 2A/2B 已按 figure2/ToDO_fig2.md 重构（TIC → log2(x+1) → Z-score；ComplexHeatmap；Top50 按 max(|log2FC|)）。
# 请直接运行：Rscript figure2/fig2_make_main.R
# 若仍想从此脚本入口运行，也可直接执行本文件。

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

source(file.path(fig2_local_script_dir(), "fig2_make_main.R"))

if (sys.nframe() == 0) {
  fig2_make_main()
}
