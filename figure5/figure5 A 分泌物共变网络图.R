# Figure 5A：分泌物共变网络图
# 说明：已按 ToDO_fig5.md 将“数据路径 + BH-FDR + 拓扑参数/敏感性/置换检验”等整合到统一脚本。
# 建议直接运行：
#   Rscript figure5/fig5_make_all.R
#
# 如果只想生成 Figure 5A（阈值 0.7 的代谢物网络），可运行本脚本：
#   Rscript "figure5/figure5 A 分泌物共变网络图.R"

this_file <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
  }
  stop("Cannot determine the current script path. Please run via Rscript --file=... or set working directory to figure5/.", call. = FALSE)
}

source(file.path(dirname(this_file()), "fig5_make_all.R"))
run_fig5("exudate")
