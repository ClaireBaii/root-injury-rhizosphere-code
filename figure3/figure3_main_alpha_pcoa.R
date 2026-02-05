#!/usr/bin/env Rscript

# Figure 3（α多样性 + PCoA）
# R 4.4.3
#
# 输入数据：
# - data/完整数据-微生物.csv（ASV/OTU counts + taxonomy）
#
# 输出（写入 figure3/）：
# - Fig3_main.pdf + Fig3_main.tiff
# - Fig3_alpha_stats.csv（组均值±SD + KW p + Dunn(BH) q）
# - Fig3_pcoa_scores.csv（可选：PCoA坐标）
# - Fig3_permanova.csv（PERMANOVA 结果）
#
# 说明：
# - Shannon 指数使用 log2（与原图坐标轴一致：~11–12）。
# - Faith's PD 需要系统发育树；本脚本提供“可选读取外部 FaithPD 表”或“提供 tree.nwk 计算”的接口。

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

chao1_single <- function(x) {
  x <- x[!is.na(x) & x > 0]
  s_obs <- length(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if (f2 > 0) {
    s_obs + (f1 * f1) / (2 * f2)
  } else {
    # bias-corrected when f2 = 0
    s_obs + (f1 * (f1 - 1)) / 2
  }
}

repo_root <- get_repo_root()
in_csv <- file.path(repo_root, "data", "完整数据-微生物.csv")
out_dir <- file.path(repo_root, "figure3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 可配置项 ----------
permutations <- 9999
control_group <- "Control"

# 若你已经用 QIIME2/其他流程算好了 Faith's PD，可在这里填路径（示例：figure3/faith_pd.tsv）
# 支持两列：Sample, Faith_pd（列名大小写不敏感）
faith_pd_table <- NA_character_

# 若你有系统发育树（Newick），也可以在这里填路径，然后用 picante::pd 计算 Faith's PD
phylo_tree_file <- NA_character_
# ------------------------------

require_cran(c(
  "vegan", "ggplot2", "dplyr", "tidyr", "rstatix", "ggpubr",
  "patchwork", "readr"
))

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
  library(readr)
})

if (!file.exists(in_csv)) {
  stop(sprintf("Input not found: %s", in_csv))
}

message("Reading OTU/ASV table: ", in_csv)
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

# 去掉全零 ASV（加速后续计算）
otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]

meta <- data.frame(
  Sample = sample_ids,
  Group = vapply(sample_ids, infer_group, character(1)),
  stringsAsFactors = FALSE
)
group_order <- c("Control", "Treat1", "Treat2", "Treat3", "Treat4")
meta$Group <- factor(meta$Group, levels = unique(c(group_order, sort(unique(meta$Group)))))

# ---------- α 多样性 ----------
message("Computing alpha diversity metrics ...")
observed <- colSums(otu_mat > 0)
chao1 <- apply(otu_mat, 2, chao1_single)
shannon <- vegan::diversity(t(otu_mat), index = "shannon", base = 2)
simpson <- vegan::diversity(t(otu_mat), index = "simpson")

faith_pd <- rep(NA_real_, length(sample_ids))
names(faith_pd) <- sample_ids

read_faith_pd_table <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  x <- suppressWarnings(readr::read_delim(path, delim = "\t", show_col_types = FALSE))
  if (!is.data.frame(x) || nrow(x) == 0) return(NULL)
  names(x) <- tolower(names(x))
  if (!all(c("sample", "faith_pd") %in% names(x))) {
    return(NULL)
  }
  x <- x[, c("sample", "faith_pd")]
  x$sample <- as.character(x$sample)
  x$faith_pd <- as.numeric(x$faith_pd)
  x
}

faith_tbl <- read_faith_pd_table(faith_pd_table)
if (!is.null(faith_tbl)) {
  message("Using external Faith_pd table: ", faith_pd_table)
  idx <- match(sample_ids, faith_tbl$sample)
  faith_pd[!is.na(idx)] <- faith_tbl$faith_pd[idx[!is.na(idx)]]
}

# 如果提供了 phylo tree，则尝试用 picante::pd 计算 Faith_pd（需要安装 ape + picante）
if (all(is.na(faith_pd)) && !is.na(phylo_tree_file) && nzchar(phylo_tree_file) && file.exists(phylo_tree_file)) {
  if (!requireNamespace("ape", quietly = TRUE) || !requireNamespace("picante", quietly = TRUE)) {
    message("phylo_tree_file provided but missing packages: ape / picante. Skip Faith_pd.")
  } else {
    message("Computing Faith_pd from tree: ", phylo_tree_file)
    tree <- ape::read.tree(phylo_tree_file)
    # picante::pd expects sites x species
    pd_res <- picante::pd(t(otu_mat), tree, include.root = TRUE)
    faith_pd <- pd_res$PD
    names(faith_pd) <- rownames(pd_res)
    faith_pd <- faith_pd[sample_ids]
  }
}

alpha_df <- meta %>%
  mutate(
    Observed = as.numeric(observed[Sample]),
    Chao1 = as.numeric(chao1[Sample]),
    Shannon = as.numeric(shannon[Sample]),
    Simpson = as.numeric(simpson[Sample]),
    Faith_pd = as.numeric(faith_pd[Sample])
  )

metric_colors <- c(
  Chao1 = "#F8766D",
  Faith_pd = "#7CAE00",
  Shannon = "#00BFC4",
  Simpson = "#C77CFF"
)

alpha_metrics <- c("Chao1", "Faith_pd", "Shannon", "Simpson")
alpha_metrics <- alpha_metrics[alpha_metrics %in% names(alpha_df)]
if (all(is.na(alpha_df$Faith_pd))) {
  alpha_metrics <- setdiff(alpha_metrics, "Faith_pd")
  message("Faith_pd is NA for all samples. Alpha plot will omit Faith_pd (provide faith_pd_table or phylo_tree_file to enable).")
}

alpha_long <- alpha_df %>%
  select(Sample, Group, all_of(alpha_metrics)) %>%
  pivot_longer(cols = all_of(alpha_metrics), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = alpha_metrics),
    Group = factor(Group, levels = levels(meta$Group))
  )

message("Running Kruskal-Wallis + Dunn(BH) ...")
kw <- alpha_long %>%
  group_by(Metric) %>%
  kruskal_test(Value ~ Group) %>%
  ungroup() %>%
  select(Metric, kw_p = p)

dunn_all <- alpha_long %>%
  group_by(Metric) %>%
  dunn_test(Value ~ Group, p.adjust.method = "BH") %>%
  ungroup()

dunn_vs_ctrl <- dunn_all %>%
  filter(group1 == control_group | group2 == control_group) %>%
  mutate(
    other = if_else(group1 == control_group, group2, group1),
    group1 = control_group,
    group2 = other
  ) %>%
  select(-other)

metric_ranges <- alpha_long %>%
  group_by(Metric) %>%
  summarise(
    y_min = min(Value, na.rm = TRUE),
    y_max = max(Value, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  mutate(
    y_step = pmax(y_range, abs(y_max) * 0.05, 1e-6),
    y_base = y_max + y_step * 0.10
  )

dunn_annot <- dunn_vs_ctrl %>%
  left_join(metric_ranges, by = "Metric") %>%
  group_by(Metric) %>%
  arrange(factor(group2, levels = levels(meta$Group)), .by_group = TRUE) %>%
  mutate(y.position = y_base + (row_number() - 1) * y_step * 0.18) %>%
  ungroup() %>%
  mutate(
    group1 = as.character(group1),
    group2 = as.character(group2)
  )

alpha_summary <- alpha_long %>%
  group_by(Metric, Group) %>%
  summarise(
    n = dplyr::n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(kw, by = "Metric") %>%
  left_join(
    dunn_vs_ctrl %>%
      transmute(
        Metric,
        Group = factor(group2, levels = levels(meta$Group)),
        dunn_q_vs_control = p.adj,
        dunn_q_signif = p.adj.signif
      ),
    by = c("Metric", "Group")
  )

write.csv(alpha_summary, file.path(out_dir, "Fig3_alpha_stats.csv"), row.names = FALSE)
write.csv(dunn_all, file.path(out_dir, "Fig3_alpha_dunn_all_pairs.csv"), row.names = FALSE)

# ---------- α 多样性绘图 ----------
p_alpha <- ggplot(alpha_long, aes(x = Group, y = Value)) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.85, color = "grey35") +
  stat_summary(aes(color = Metric), fun = mean, geom = "point", size = 2.6) +
  stat_summary(
    aes(color = Metric),
    fun.data = ggplot2::mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    linewidth = 0.7,
    width = 0.20
  ) +
  scale_color_manual(values = metric_colors, guide = "none") +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Value") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey85", color = "grey45"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

if (nrow(dunn_annot) > 0) {
  p_alpha <- p_alpha +
    ggpubr::stat_pvalue_manual(
      dunn_annot,
      label = "p.adj.signif",
      y.position = "y.position",
      xmin = "group1",
      xmax = "group2",
      size = 3,
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

# ---------- PCoA（Bray-Curtis）+ PERMANOVA ----------
message("Computing Bray-Curtis PCoA + PERMANOVA ...")
otu_rel <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")
otu_rel[is.na(otu_rel)] <- 0
bray <- vegan::vegdist(t(otu_rel), method = "bray")

pcoa <- cmdscale(bray, k = 2, eig = TRUE)
coords <- as.data.frame(pcoa$points)
colnames(coords) <- c("PC1", "PC2")
coords$Sample <- rownames(coords)
coords <- coords %>% left_join(meta, by = "Sample")

eig <- pcoa$eig
eig_pos <- eig[eig > 0]
var_expl <- eig_pos / sum(eig_pos) * 100
pc1_pct <- ifelse(length(var_expl) >= 1, var_expl[1], NA_real_)
pc2_pct <- ifelse(length(var_expl) >= 2, var_expl[2], NA_real_)

adon <- vegan::adonis2(bray ~ Group, data = meta, permutations = permutations)
adon_df <- as.data.frame(adon)
grp_row <- adon_df[rownames(adon_df) == "Group", , drop = FALSE]
permanova_r2 <- as.numeric(grp_row$R2[1])
permanova_p <- as.numeric(grp_row$`Pr(>F)`[1])

permanova_tbl <- data.frame(
  term = "Group",
  R2 = permanova_r2,
  p_value = permanova_p,
  permutations = permutations,
  distance = "Bray-Curtis (relative abundance)",
  stringsAsFactors = FALSE
)
write.csv(permanova_tbl, file.path(out_dir, "Fig3_permanova.csv"), row.names = FALSE)

write.csv(
  coords %>%
    mutate(
      PC1_percent = pc1_pct,
      PC2_percent = pc2_pct
    ) %>%
    select(Sample, Group, PC1, PC2, PC1_percent, PC2_percent),
  file.path(out_dir, "Fig3_pcoa_scores.csv"),
  row.names = FALSE
)

permanova_text <- sprintf(
  "PERMANOVA: R² = %.3f, P = %.3g (perm=%d)",
  permanova_r2, permanova_p, permutations
)

group_colors <- c(
  Control = "#F8766D",
  Treat1 = "#7CAE00",
  Treat2 = "#00BFC4",
  Treat3 = "#00B0F6",
  Treat4 = "#C77CFF"
)

p_pcoa <- ggplot(coords, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3.2, alpha = 0.95) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pc1_pct),
    y = sprintf("PC2 (%.1f%%)", pc2_pct),
    subtitle = permanova_text,
    color = "Group"
  ) +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ---------- 合并输出 ----------
fig_main <- p_alpha + p_pcoa +
  patchwork::plot_layout(widths = c(1.2, 1)) +
  patchwork::plot_annotation(tag_levels = "A")

pdf_path <- file.path(out_dir, "Fig3_main.pdf")
tiff_path <- file.path(out_dir, "Fig3_main.tiff")

message("Saving: ", pdf_path)
ggsave(pdf_path, fig_main, width = 12, height = 5.6, units = "in")

message("Saving: ", tiff_path)
ggsave(
  tiff_path,
  fig_main,
  width = 12,
  height = 5.6,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)

message("Done.")
