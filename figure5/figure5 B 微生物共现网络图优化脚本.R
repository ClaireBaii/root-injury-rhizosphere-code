# Figure 5B：微生物共现网络图
# 说明：已按 ToDO_fig5.md 将“taxonomy 解析→phylum/genus 聚合 + BH-FDR + 阈值敏感性 + 拓扑参数”等整合到统一脚本。
# 建议直接运行：
#   Rscript figure5/fig5_make_all.R
#
# 如果只想生成 Figure 5B（阈值 0.7 的微生物网络），可运行本脚本：
#   Rscript "figure5/figure5 B 微生物共现网络图优化脚本.R"

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
run_fig5("microbe")
