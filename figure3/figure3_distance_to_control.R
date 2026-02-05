#!/usr/bin/env Rscript

# FigS（建议：distance-to-control）
# 目的：量化“处理组到 Control 的 Bray–Curtis distance-to-control”
# R 4.4.3
#
# 输入：
# - data/完整数据-微生物.csv
#
# 输出（写入 figure3/）：
# - FigS_distance_to_control.pdf + FigS_distance_to_control.tiff
# - distance_to_control_stats.csv（均值±SD + KW p + Dunn(BH) q，包含两种 distance-to-control 定义）
#
# distance-to-control 两种定义（都会计算并写入 stats）：
# 1) PairwiseMean：样本到所有 Control 样本的 Bray–Curtis 距离均值（Control 自身为 leave-one-out）
# 2) Centroid：样本到 Control 组“组成均值”（centroid）的 Bray–Curtis 距离

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

repo_root <- get_repo_root()
in_csv <- file.path(repo_root, "data", "完整数据-微生物.csv")
out_dir <- file.path(repo_root, "figure3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

control_group <- "Control"
permutations <- 9999 # 仅用于记录；distance-to-control 这里用 KW + Dunn
plot_definition <- "Centroid" # "PairwiseMean" 或 "Centroid"

require_cran(c("vegan", "ggplot2", "dplyr", "tidyr", "rstatix", "ggpubr"))
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rstatix)
  library(ggpubr)
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

control_samples <- meta$Sample[meta$Group == control_group]
if (length(control_samples) < 2) {
  stop("Need at least 2 Control samples for leave-one-out pairwise distance-to-control.")
}

# 使用相对丰度计算 Bray–Curtis
otu_rel <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")
otu_rel[is.na(otu_rel)] <- 0
otu_rel_t <- t(otu_rel) # samples x ASVs

bray <- vegan::vegdist(otu_rel_t, method = "bray")
bray_mat <- as.matrix(bray)

dist_pairwise <- vapply(sample_ids, function(s) {
  cs <- control_samples
  if (s %in% cs) cs <- setdiff(cs, s) # leave-one-out for Control
  if (length(cs) == 0) return(NA_real_)
  mean(bray_mat[s, cs], na.rm = TRUE)
}, numeric(1))

control_centroid <- colMeans(otu_rel_t[control_samples, , drop = FALSE])
dist_centroid <- vapply(sample_ids, function(s) {
  x <- otu_rel_t[s, ]
  # Bray–Curtis = sum(|x - y|) / sum(x + y); for relative abundances sum(x)=sum(y)=1 -> denom=2
  0.5 * sum(abs(x - control_centroid))
}, numeric(1))

dctl_df <- meta %>%
  mutate(
    PairwiseMean = as.numeric(dist_pairwise[Sample]),
    Centroid = as.numeric(dist_centroid[Sample])
  )

dctl_long <- dctl_df %>%
  pivot_longer(cols = c("PairwiseMean", "Centroid"), names_to = "Definition", values_to = "Distance") %>%
  mutate(
    Definition = factor(Definition, levels = c("Centroid", "PairwiseMean")),
    Group = factor(Group, levels = levels(meta$Group))
  )

kw <- dctl_long %>%
  group_by(Definition) %>%
  kruskal_test(Distance ~ Group) %>%
  ungroup() %>%
  select(Definition, kw_p = p)

dunn_all <- dctl_long %>%
  group_by(Definition) %>%
  dunn_test(Distance ~ Group, p.adjust.method = "BH") %>%
  ungroup()

dunn_vs_ctrl <- dunn_all %>%
  filter(group1 == control_group | group2 == control_group) %>%
  mutate(
    other = if_else(group1 == control_group, group2, group1),
    group1 = control_group,
    group2 = other
  ) %>%
  select(-other)

ranges <- dctl_long %>%
  group_by(Definition) %>%
  summarise(
    y_min = min(Distance, na.rm = TRUE),
    y_max = max(Distance, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    y_step = pmax(y_range, abs(y_max) * 0.05, 1e-6),
    y_base = y_max + y_step * 0.10
  )

dunn_annot <- dunn_vs_ctrl %>%
  left_join(ranges, by = "Definition") %>%
  group_by(Definition) %>%
  arrange(group2, .by_group = TRUE) %>%
  mutate(y.position = y_base + (row_number() - 1) * y_step * 0.18) %>%
  ungroup()

summary_tbl <- dctl_long %>%
  group_by(Definition, Group) %>%
  summarise(
    n = dplyr::n(),
    mean = mean(Distance, na.rm = TRUE),
    sd = sd(Distance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(kw, by = "Definition") %>%
  left_join(
    dunn_vs_ctrl %>%
      transmute(
        Definition,
        Group = group2,
        dunn_q_vs_control = p.adj,
        dunn_q_signif = p.adj.signif
      ),
    by = c("Definition", "Group")
  )

write.csv(summary_tbl, file.path(out_dir, "distance_to_control_stats.csv"), row.names = FALSE)

# ---------- 绘图（默认只画一种定义；另一定义在 stats 里） ----------
plot_df <- dctl_long %>% filter(Definition == plot_definition)
plot_ann <- dunn_annot %>% filter(Definition == plot_definition)
kw_p <- kw %>% filter(Definition == plot_definition) %>% pull(kw_p)
kw_txt <- if (length(kw_p) == 1) sprintf("Kruskal–Wallis: P = %.3g", kw_p) else NULL

p <- ggplot(plot_df, aes(x = Group, y = Distance)) +
  geom_jitter(width = 0.12, height = 0, size = 2.0, alpha = 0.90, color = "grey35") +
  stat_summary(fun = mean, geom = "point", size = 2.8, color = "black") +
  stat_summary(
    fun.data = ggplot2::mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    linewidth = 0.7,
    width = 0.20,
    color = "black"
  ) +
  labs(
    x = NULL,
    y = "Bray–Curtis distance-to-control",
    subtitle = kw_txt
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

if (nrow(plot_ann) > 0) {
  p <- p +
    ggpubr::stat_pvalue_manual(
      plot_ann,
      label = "p.adj.signif",
      y.position = "y.position",
      xmin = "group1",
      xmax = "group2",
      size = 3,
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

pdf_path <- file.path(out_dir, "FigS_distance_to_control.pdf")
tiff_path <- file.path(out_dir, "FigS_distance_to_control.tiff")

ggsave(pdf_path, p, width = 6.5, height = 4.2, units = "in")
ggsave(
  tiff_path,
  p,
  width = 6.5,
  height = 4.2,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

message("Done. Stats: ", file.path(out_dir, "distance_to_control_stats.csv"))

