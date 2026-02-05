# NOTE (2026-02): This script is superseded by `figure9_make_fig9.R`.
# It now acts as a thin wrapper so existing workflows still work.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
source(file.path(script_dir, "figure9_make_fig9.R"))
