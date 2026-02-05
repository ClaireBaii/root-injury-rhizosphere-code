#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_pkgs <- c("dplyr", "ggplot2", "tidyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing R packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))",
      call. = FALSE
    )
  }

  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 0) {
    return(NULL)
  }
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
  dirname(script_path)
}

script_dir <- get_script_dir()
find_project_root <- function(start_dir, max_up = 6) {
  dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_up + 1)) {
    if (file.exists(file.path(dir, "data", "完整数据-微生物.csv"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (parent == dir) {
      break
    }
    dir <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = FALSE)
}

base_dir <- if (!is.null(script_dir)) script_dir else getwd()
project_root <- find_project_root(base_dir)
out_dir <- file.path(project_root, "figure6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_candidates <- c(
  file.path(project_root, "data", "完整数据-微生物.csv"),
  file.path(getwd(), "data", "完整数据-微生物.csv"),
  "/mnt/data/完整数据-微生物.csv"
)
existing_inputs <- input_candidates[file.exists(input_candidates)]
if (length(existing_inputs) == 0) {
  stop(
    "Cannot find input CSV.\nTried:\n- ",
    paste(input_candidates, collapse = "\n- "),
    call. = FALSE
  )
}
input_path <- existing_inputs[[1]]

message("Reading: ", input_path)
otu <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)

id_col <- if ("OTU ID" %in% names(otu)) "OTU ID" else names(otu)[[1]]
if (!("taxonomy" %in% names(otu))) {
  stop("Input must contain a `taxonomy` column.", call. = FALSE)
}

sample_cols <- setdiff(names(otu), c(id_col, "taxonomy"))
if (length(sample_cols) == 0) {
  stop("No sample columns found (expected columns besides ID and taxonomy).", call. = FALSE)
}

extract_phylum <- function(taxonomy) {
  has_phylum <- grepl("p__", taxonomy, fixed = TRUE)
  phylum <- ifelse(has_phylum, sub(".*p__([^;]+).*", "\\1", taxonomy), NA_character_)
  phylum <- trimws(phylum)
  phylum[is.na(phylum) | phylum == "" | phylum %in% c("unclassified", "Unclassified")] <- "Unassigned"
  phylum
}

otu$Phylum <- extract_phylum(otu$taxonomy)

counts <- as.matrix(otu[, sample_cols, drop = FALSE])
storage.mode(counts) <- "numeric"
counts[is.na(counts)] <- 0

message("Aggregating to phylum ...")
phylum_counts <- rowsum(counts, group = otu$Phylum, reorder = TRUE, na.rm = TRUE)
sample_totals <- colSums(phylum_counts)
rel_abund <- sweep(phylum_counts, 2, sample_totals, FUN = "/")
rel_abund[, sample_totals == 0] <- NA_real_

rel_df <- as.data.frame(rel_abund, check.names = FALSE)
rel_df$Phylum <- rownames(rel_df)
rel_long <- rel_df %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "RelAbundance") %>%
  mutate(RelAbundance = as.numeric(RelAbundance))

parse_sample <- function(sample_name) {
  m <- regexec("^RDYS_([0-9]+)-([0-9]+)$", sample_name)
  parts <- regmatches(sample_name, m)[[1]]
  if (length(parts) == 3) {
    list(treat_id = as.integer(parts[[2]]), tree_id = as.integer(parts[[3]]))
  } else {
    list(treat_id = NA_integer_, tree_id = NA_integer_)
  }
}

sample_meta <- tibble(Sample = unique(rel_long$Sample)) %>%
  rowwise() %>%
  mutate(parsed = list(parse_sample(Sample))) %>%
  mutate(treat_id = parsed$treat_id, tree_id = parsed$tree_id) %>%
  ungroup() %>%
  select(-parsed)

if (any(is.na(sample_meta$treat_id))) {
  bad <- sample_meta$Sample[is.na(sample_meta$treat_id)]
  stop(
    "Unexpected sample names (cannot parse treatment/tree): ",
    paste(utils::head(bad, 10), collapse = ", "),
    if (length(bad) > 10) sprintf(" ... (+%d more)", length(bad) - 10) else "",
    call. = FALSE
  )
}

rs_map <- c(`1` = 0, `2` = 25, `3` = 50, `4` = 75, `5` = 100)
sample_meta <- sample_meta %>%
  mutate(
    RS = unname(rs_map[as.character(treat_id)]),
    RS = as.integer(RS)
  )

if (any(is.na(sample_meta$RS))) {
  stop("RS mapping failed: found treatment IDs outside 1..5.", call. = FALSE)
}

rs_levels <- sort(unique(sample_meta$RS))
rel_long <- rel_long %>%
  left_join(sample_meta, by = "Sample") %>%
  mutate(RS = factor(RS, levels = rs_levels))

summary_df <- rel_long %>%
  group_by(Phylum, RS) %>%
  summarise(
    n = sum(!is.na(RelAbundance)),
    mean = mean(RelAbundance, na.rm = TRUE),
    sd = sd(RelAbundance, na.rm = TRUE),
    .groups = "drop"
  )

overall_df <- rel_long %>%
  group_by(Phylum) %>%
  summarise(overall_mean = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

trend_df <- summary_df %>%
  group_by(Phylum) %>%
  summarise(
    rho = suppressWarnings(cor(as.numeric(as.character(RS)), mean, method = "spearman", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  left_join(overall_df, by = "Phylum")

TOP_N <- 10
RHO_THRESHOLD <- 0.7
MIN_OVERALL_MEAN <- 0.001

trend_selected <- trend_df %>%
  filter(!is.na(rho), overall_mean >= MIN_OVERALL_MEAN, abs(rho) >= RHO_THRESHOLD) %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = TOP_N)

if (nrow(trend_selected) > 0) {
  phyla_to_plot <- trend_selected$Phylum
  selection_note <- sprintf("Selected %d trend phyla (|rho| >= %.1f).", length(phyla_to_plot), RHO_THRESHOLD)
} else {
  phyla_to_plot <- overall_df %>%
    arrange(desc(overall_mean)) %>%
    slice_head(n = TOP_N) %>%
    pull(Phylum)
  selection_note <- sprintf("No phyla passed trend filter; using Top %d by overall mean.", TOP_N)
}
message(selection_note)

rel_plot_long <- rel_long %>% filter(Phylum %in% phyla_to_plot)
summary_plot <- summary_df %>%
  filter(Phylum %in% phyla_to_plot) %>%
  mutate(
    mean_pct = mean * 100,
    sd_pct = sd * 100
  )

safe_kw_p <- function(values, groups) {
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]
  groups <- droplevels(groups[ok])
  if (length(values) < 2 || nlevels(groups) < 2 || length(unique(values)) < 2) {
    return(1)
  }
  suppressWarnings(kruskal.test(values ~ groups)$p.value)
}

dunn_test <- function(values, groups, p_adjust_method = "BH") {
  empty <- tibble(group1 = character(), group2 = character(), z = numeric(), p_value = numeric(), p_adj = numeric())
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]
  groups <- droplevels(groups[ok])
  if (length(values) < 2 || nlevels(groups) < 2 || length(unique(values)) < 2) {
    return(empty)
  }

  ranks <- rank(values, ties.method = "average")
  n <- length(values)

  tie_counts <- table(values)
  tie_term <- sum(tie_counts^3 - tie_counts)
  tie_correction <- if (n > 1) 1 - tie_term / (n^3 - n) else 1
  if (!is.finite(tie_correction) || tie_correction <= 0) {
    tie_correction <- 1
  }

  mean_rank <- tapply(ranks, groups, mean)
  n_i <- tapply(ranks, groups, length)
  levels_g <- names(mean_rank)

  pairs <- combn(levels_g, 2, simplify = FALSE)
  res <- lapply(pairs, function(pair) {
    g1 <- pair[[1]]
    g2 <- pair[[2]]
    denom <- sqrt((n * (n + 1) / 12) * tie_correction * (1 / n_i[[g1]] + 1 / n_i[[g2]]))
    if (!is.finite(denom) || denom <= 0) {
      return(tibble(group1 = g1, group2 = g2, z = NA_real_, p_value = NA_real_))
    }
    z <- (mean_rank[[g1]] - mean_rank[[g2]]) / denom
    p <- 2 * pnorm(-abs(z))
    tibble(group1 = g1, group2 = g2, z = z, p_value = p)
  }) %>%
    bind_rows()

  res %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
}

kw_df <- rel_plot_long %>%
  group_by(Phylum) %>%
  summarise(kw_p = safe_kw_p(RelAbundance, RS), .groups = "drop") %>%
  mutate(kw_q = p.adjust(kw_p, method = "BH"))

dunn_df <- rel_plot_long %>%
  group_by(Phylum) %>%
  group_modify(~ dunn_test(.x$RelAbundance, .x$RS)) %>%
  ungroup()

stats_out <- bind_rows(
  kw_df %>%
    transmute(
      Phylum,
      test = "Kruskal-Wallis",
      group1 = NA_character_,
      group2 = NA_character_,
      z = NA_real_,
      p_value = kw_p,
      p_adj = kw_q,
      p_adj_method = "BH",
      p_adj_scope = "across_phyla"
    ),
  dunn_df %>%
    transmute(
      Phylum,
      test = "Dunn",
      group1 = as.character(group1),
      group2 = as.character(group2),
      z = as.numeric(z),
      p_value = as.numeric(p_value),
      p_adj = as.numeric(p_adj),
      p_adj_method = "BH",
      p_adj_scope = "within_phylum"
    )
) %>%
  arrange(Phylum, factor(test, levels = c("Kruskal-Wallis", "Dunn")), group1, group2)

stats_path <- file.path(out_dir, "Fig6_phylum_stats.csv")
write.csv(stats_out, stats_path, row.names = FALSE, fileEncoding = "UTF-8")
message("Wrote: ", stats_path)

sig_levels <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "NA",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

annot_df <- summary_plot %>%
  group_by(Phylum) %>%
  summarise(y = max(mean_pct + sd_pct, na.rm = TRUE), .groups = "drop") %>%
  left_join(kw_df, by = "Phylum") %>%
  mutate(
    y = ifelse(is.finite(y) & y > 0, y * 1.08, 0.1),
    RS = factor(max(rs_levels), levels = rs_levels),
    label = paste0("q=", signif(kw_q, 2), " ", sig_levels(kw_q))
  )

n_phyla <- length(phyla_to_plot)
ncol_facets <- 4
nrow_facets <- ceiling(n_phyla / ncol_facets)
fig_w <- 10
fig_h <- max(4.8, nrow_facets * 2.2)

caption_text <- paste(
  "Treatment-level mean across trees (n=3); error bars indicate ±SD.",
  "P-values: Kruskal–Wallis across treatments; BH-adjusted across plotted phyla."
)

p <- ggplot(summary_plot, aes(x = RS, y = mean_pct, group = 1)) +
  geom_line(linewidth = 0.4, color = "grey30") +
  geom_pointrange(
    aes(ymin = pmax(mean_pct - sd_pct, 0), ymax = mean_pct + sd_pct),
    linewidth = 0.35,
    fatten = 1.4
  ) +
  geom_text(
    data = annot_df,
    aes(x = RS, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 3.0
  ) +
  facet_wrap(~Phylum, scales = "free_y", ncol = ncol_facets) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  labs(
    x = "RS (%)",
    y = "Relative abundance (%)",
    title = "Figure 6. Phylum relative abundance trajectories",
    subtitle = selection_note,
    caption = caption_text
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 9)
  )

pdf_path <- file.path(out_dir, "Fig6_main.pdf")
tiff_path <- file.path(out_dir, "Fig6_main.tiff")

if (capabilities("cairo")) {
  cairo_pdf(filename = pdf_path, width = fig_w, height = fig_h)
} else {
  pdf(file = pdf_path, width = fig_w, height = fig_h)
}
print(p)
dev.off()
message("Wrote: ", pdf_path)

if (capabilities("tiff")) {
  tiff(filename = tiff_path, width = fig_w, height = fig_h, units = "in", res = 300, compression = "lzw")
} else if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_tiff(filename = tiff_path, width = fig_w, height = fig_h, units = "in", res = 300)
} else {
  stop(
    "TIFF device is not available in this R build.\n",
    "Try installing `ragg` (install.packages('ragg')) or use an R build with TIFF support.",
    call. = FALSE
  )
}
print(p)
dev.off()
message("Wrote: ", tiff_path)
