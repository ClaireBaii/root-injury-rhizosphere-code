# Figure 5（网络：代谢物网络 + 微生物网络）
# R >= 4.4.3
#
# 输入数据（来自 ../data/）：
#   - 完整数据-分泌物.csv
#   - 完整数据-微生物.csv
#
# 主要输出（默认写入本脚本所在目录 figure5/）：
#   - Fig5_main.pdf / Fig5_main.tiff
#   - Fig5_network_topology.csv
#   - FigS5_threshold_sensitivity.pdf / FigS5_threshold_sensitivity_table.csv
#   - FigS_posneg_permutation.pdf / posneg_test_results.csv

options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "igraph",
  "Hmisc",
  "ggraph",
  "ggplot2",
  "RColorBrewer",
  "scales"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
      "Please install them first, e.g.:\n",
      "  install.packages(c(", paste(sprintf('\"%s\"', missing_pkgs), collapse = ", "), "))\n"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(igraph)
  library(Hmisc)
  library(ggraph)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
})

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

script_path <- this_file()
script_dir <- dirname(script_path)

pick_existing <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }
  stop("File not found. Tried:\n  - ", paste(paths, collapse = "\n  - "), call. = FALSE)
}

data_exudate <- pick_existing(c(
  file.path(script_dir, "..", "data", "完整数据-分泌物.csv"),
  file.path(getwd(), "data", "完整数据-分泌物.csv"),
  file.path(getwd(), "..", "data", "完整数据-分泌物.csv")
))

data_microbe <- pick_existing(c(
  file.path(script_dir, "..", "data", "完整数据-微生物.csv"),
  file.path(getwd(), "data", "完整数据-微生物.csv"),
  file.path(getwd(), "..", "data", "完整数据-微生物.csv")
))

out_dir <- script_dir

config <- list(
  fdr_alpha = 0.05,
  thr_exudate = 0.80,          # 代谢物网络阈值（稍高以减少密度）
  thr_microbe = 0.60,          # 微生物网络阈值（降低以获得更多边）
  thr_sensitivity_high = 0.7,  # 敏感性分析高阈值
  thr_sensitivity_low = 0.5,   # 敏感性分析低阈值
  hub_top_n = 8,               # 标注的 top hub 数量
  seed = 42,
  # 为避免 1800+ 代谢物网络"毛球"过密，默认仅取 max(|log2FC|) Top N 代谢物构网；
  # 如需使用全部代谢物：将 exudate_filter = NULL（或把 exudate_top_n 调大）。
  exudate_filter = "max_abs_log2fc", # NULL / "max_abs_log2fc"
  exudate_top_n = 60,         # 适当增加节点数量以显示更丰富的网络结构
  microbe_rank = "genus",     # 改为 genus 以获得更多节点
  microbe_top_n = 40,         # 增加 taxa 数量
  microbe_min_prevalence = 0.15, # 降低 prevalence 阈值
  perm_n = 10000
)

extract_tax_rank <- function(taxonomy, rank = c("phylum", "genus")) {
  rank <- match.arg(rank)
  prefix <- switch(rank, phylum = "p__", genus = "g__")
  out <- rep(NA_character_, length(taxonomy))

  for (i in seq_along(taxonomy)) {
    tx <- taxonomy[i]
    if (is.na(tx) || !nzchar(tx)) next
    parts <- strsplit(tx, ";", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    hit <- parts[startsWith(parts, prefix)]
    if (length(hit) == 0) next
    val <- sub(paste0("^", prefix), "", hit[1])
    val <- trimws(val)
    if (!nzchar(val) || grepl("^unclassified|^uncultured", val, ignore.case = TRUE)) next
    out[i] <- val
  }

  out[is.na(out)] <- paste0("Unclassified_", rank)
  out
}

read_exudate_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  if (!("Name" %in% names(df))) stop("Exudate file missing column: Name", call. = FALSE)

  if (!is.null(config$exudate_filter) && identical(config$exudate_filter, "max_abs_log2fc")) {
    log2_cols <- grep("^log2FC_", names(df), value = TRUE)
    if (length(log2_cols) > 0 && is.finite(config$exudate_top_n) && config$exudate_top_n > 0) {
      log2_mat <- as.matrix(df[, log2_cols, drop = FALSE])
      storage.mode(log2_mat) <- "double"
      score <- apply(abs(log2_mat), 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
      ord <- order(score, decreasing = TRUE, na.last = NA)
      if (length(ord) == 0) {
        warning("All log2FC_ values are NA; using all exudate features.", call. = FALSE)
      } else {
        if (length(ord) > config$exudate_top_n) ord <- ord[seq_len(config$exudate_top_n)]
        df <- df[ord, , drop = FALSE]
      }
    } else {
      warning("No usable log2FC_ columns found; using all exudate features.", call. = FALSE)
    }
  }

  sample_cols <- grep("^YS", names(df), value = TRUE)
  if (length(sample_cols) < 3) stop("Exudate sample columns not found (expected columns starting with 'YS').", call. = FALSE)

  node_name <- df$Name
  node_name <- ifelse(is.na(node_name) | !nzchar(node_name), paste0("Metabolite_", seq_len(nrow(df))), node_name)
  node_name <- make.unique(node_name)

  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- node_name

  keep <- apply(mat, 1, function(x) isTRUE(sd(x, na.rm = TRUE) > 0))
  mat <- mat[keep, , drop = FALSE]

  meta_cols <- intersect(c("Super_Class", "Class", "Sub_Class"), names(df))
  meta <- df[keep, c("Name", meta_cols), drop = FALSE]
  rownames(meta) <- rownames(mat)

  list(mat = mat, meta = meta, sample_cols = sample_cols)
}

read_microbe_matrix <- function(path,
                               rank = c("phylum", "genus"),
                               top_n = 30,
                               min_prevalence = 0.2) {
  rank <- match.arg(rank)
  df <- read.csv(path, check.names = FALSE)

  if (!("taxonomy" %in% names(df))) stop("Microbe file missing column: taxonomy", call. = FALSE)
  if (!("OTU ID" %in% names(df))) stop("Microbe file missing column: OTU ID", call. = FALSE)

  sample_cols <- setdiff(names(df), c("OTU ID", "taxonomy"))
  if (length(sample_cols) < 3) stop("Microbe sample columns not found.", call. = FALSE)

  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"

  group <- extract_tax_rank(df$taxonomy, rank = rank)
  agg <- rowsum(mat, group = group, reorder = FALSE)

  # 将 NA 值替换为 0，避免 colSums 返回 NA
  agg[is.na(agg)] <- 0
  
  col_sums <- colSums(agg)
  if (any(is.na(col_sums)) || any(col_sums <= 0)) {
    warning("Some samples have zero total counts; filtering out these samples.", call. = FALSE)
    valid_cols <- !is.na(col_sums) & col_sums > 0
    agg <- agg[, valid_cols, drop = FALSE]
    col_sums <- col_sums[valid_cols]
    sample_cols <- sample_cols[valid_cols]
  }
  rel <- sweep(agg, 2, col_sums, "/")

  prevalence <- rowMeans(rel > 0)
  rel <- rel[prevalence >= min_prevalence, , drop = FALSE]

  keep <- apply(rel, 1, function(x) isTRUE(sd(x, na.rm = TRUE) > 0))
  rel <- rel[keep, , drop = FALSE]

  if (!is.null(top_n) && is.finite(top_n) && nrow(rel) > top_n) {
    m <- rowMeans(rel, na.rm = TRUE)
    top_taxa <- names(sort(m, decreasing = TRUE))[seq_len(top_n)]
    rel <- rel[top_taxa, , drop = FALSE]
  }

  list(mat = rel, sample_cols = sample_cols)
}

bh_adjust_matrix <- function(p_mat) {
  ut <- upper.tri(p_mat, diag = FALSE)
  p_vec <- p_mat[ut]
  p_vec[is.na(p_vec)] <- 1
  p_adj_vec <- p.adjust(p_vec, method = "BH")

  p_adj <- matrix(NA_real_, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))
  p_adj[ut] <- p_adj_vec
  p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]
  p_adj
}

build_correlation_network <- function(mat,
                                      threshold,
                                      fdr_alpha = 0.05,
                                      hub_top_n = 10,
                                      seed = 42) {
  if (nrow(mat) < 3) stop("Need at least 3 nodes to build a network.", call. = FALSE)
  if (ncol(mat) < 5) warning("Sample size < 5; correlation network may be unstable.", call. = FALSE)

  cor_res <- Hmisc::rcorr(t(mat), type = "spearman")
  r_mat <- cor_res$r
  p_mat <- cor_res$P
  p_adj <- bh_adjust_matrix(p_mat)

  ut <- upper.tri(r_mat, diag = FALSE)
  keep <- ut & (abs(r_mat) >= threshold) & (p_adj < fdr_alpha)

  idx <- which(keep, arr.ind = TRUE)
  edges <- data.frame(
    from = rownames(r_mat)[idx[, 1]],
    to = colnames(r_mat)[idx[, 2]],
    rho = r_mat[idx],
    abs_rho = abs(r_mat[idx]),
    p = p_mat[idx],
    p_adj = p_adj[idx],
    sign = ifelse(r_mat[idx] >= 0, "positive", "negative")
  )

  if (nrow(edges) == 0) {
    stop(
      sprintf("No edges passed the filter (|rho| >= %.2f and BH-FDR < %.3f).", threshold, fdr_alpha),
      call. = FALSE
    )
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  # 保留边属性，因为 simplify 会移除重复边
  # 首先创建映射以在 simplify 后恢复属性
  edge_attr_df <- edges[, c("from", "to", "rho", "abs_rho", "sign"), drop = FALSE]
  
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, 
                        edge.attr.comb = list(rho = "first", abs_rho = "first", 
                                              p = "first", p_adj = "first", 
                                              sign = "first"))
  
  # 确保边属性是数值型
  if (is.null(E(g)$rho) || any(is.na(E(g)$rho))) {
    # 如果属性丢失，从 edges 重新匹配
    el <- as_edgelist(g)
    for (i in seq_len(ecount(g))) {
      from_v <- el[i, 1]
      to_v <- el[i, 2]
      match_idx <- which((edges$from == from_v & edges$to == to_v) | 
                         (edges$from == to_v & edges$to == from_v))[1]
      if (!is.na(match_idx)) {
        E(g)[i]$rho <- edges$rho[match_idx]
        E(g)[i]$abs_rho <- edges$abs_rho[match_idx]
        E(g)[i]$sign <- edges$sign[match_idx]
      }
    }
  }
  
  E(g)$abs_rho <- as.numeric(E(g)$abs_rho)
  E(g)$rho <- as.numeric(E(g)$rho)
  E(g)$sign <- ifelse(E(g)$rho >= 0, "positive", "negative")

  cl <- igraph::cluster_louvain(g, weights = E(g)$abs_rho)
  V(g)$module <- as.factor(igraph::membership(cl))

  V(g)$degree <- igraph::degree(g, mode = "all")
  V(g)$strength <- igraph::strength(g, weights = E(g)$abs_rho)
  V(g)$betweenness <- igraph::betweenness(g, weights = 1 / E(g)$abs_rho, normalized = FALSE)

  deg <- V(g)$degree
  names(deg) <- V(g)$name
  hub_top_n <- min(hub_top_n, length(deg))
  hubs <- names(sort(deg, decreasing = TRUE))[seq_len(hub_top_n)]

  V(g)$is_hub <- V(g)$name %in% hubs
  V(g)$label <- ifelse(V(g)$is_hub, V(g)$name, "")

  list(graph = g, edges = edges, cluster = cl, hubs = hubs)
}

network_topology <- function(g, cluster) {
  comp <- igraph::components(g)
  giant_id <- which.max(comp$csize)
  giant_nodes <- sum(comp$membership == giant_id)

  g_giant <- igraph::induced_subgraph(g, vids = which(comp$membership == giant_id))
  mean_dist <- if (vcount(g_giant) > 1) igraph::mean_distance(g_giant, directed = FALSE) else NA_real_

  data.frame(
    nodes = igraph::vcount(g),
    edges = igraph::ecount(g),
    avg_degree = mean(igraph::degree(g)),
    density = igraph::edge_density(g, loops = FALSE),
    transitivity_global = igraph::transitivity(g, type = "global"),
    modularity_Q = igraph::modularity(cluster),
    n_modules = length(unique(igraph::membership(cluster))),
    pos_edges = sum(E(g)$rho > 0),
    neg_edges = sum(E(g)$rho < 0),
    pos_neg_ratio = ifelse(sum(E(g)$rho < 0) == 0, NA_real_, sum(E(g)$rho > 0) / sum(E(g)$rho < 0)),
    largest_component_nodes = giant_nodes,
    mean_distance_giant = mean_dist
  )
}

make_module_palette <- function(mod_levels) {
  n <- length(mod_levels)
  # Okabe-Ito colorblind-friendly palette (extended)
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#000000", "#88CCEE",
    "#44AA99", "#117733", "#332288", "#AA4499", "#882255",
    "#661100", "#6699CC", "#DDCC77"
  )
  cols <- if (n <= length(okabe_ito)) {
    okabe_ito[seq_len(n)]
  } else {
    grDevices::colorRampPalette(okabe_ito)(n)
  }
  names(cols) <- mod_levels
  cols
}


plot_network <- function(g,
                         title,
                         subtitle,
                         node_size_attr = c("degree", "betweenness"),
                         label_size = 2.8,
                         filter_isolates = TRUE,
                         layout_algo = c("stress", "fr", "graphopt"),
                         seed = 42) {
  node_size_attr <- match.arg(node_size_attr)
  layout_algo <- match.arg(layout_algo)
  
  # 过滤孤立节点（度=0）
  if (filter_isolates) {
    isolates <- which(igraph::degree(g) == 0)
    if (length(isolates) > 0) {
      g <- igraph::delete_vertices(g, isolates)
    }
  }
  
  # 如果过滤后没有节点，返回空图
  if (igraph::vcount(g) == 0) {
    return(ggplot() + theme_void() + 
           labs(title = title, subtitle = "No connected nodes after filtering"))
  }
  
  mod_levels <- levels(V(g)$module)
  if (is.null(mod_levels)) mod_levels <- as.character(unique(V(g)$module))
  module_cols <- make_module_palette(mod_levels)

  V(g)$size_plot <- if (node_size_attr == "degree") {
    scales::rescale(V(g)$degree, to = c(4, 14))
  } else {
    scales::rescale(log1p(V(g)$betweenness), to = c(4, 14))
  }

  set.seed(seed)
  
  # 选择布局算法
  layout_spec <- switch(layout_algo,
    "stress" = "stress",
    "fr" = "fr",
    "graphopt" = "graphopt"
  )
  
  p <- ggraph(g, layout = layout_spec) +
    geom_edge_link(
      aes(edge_colour = sign, edge_alpha = abs_rho, edge_width = abs_rho),
      show.legend = TRUE,
      lineend = "round"
    ) +
    geom_node_point(
      aes(fill = module, size = size_plot),
      shape = 21,
      color = "grey20",
      stroke = 0.4,
      alpha = 0.9,
      show.legend = TRUE
    )

  # 标签样式优化
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + geom_node_text(
      aes(label = label), 
      repel = TRUE, 
      size = label_size,
      fontface = "bold.italic",
      color = "grey10",
      max.overlaps = 20,
      point.padding = unit(0.3, "lines"),
      box.padding = unit(0.4, "lines"),
      segment.color = "grey40",
      segment.size = 0.25,
      segment.alpha = 0.6,
      min.segment.length = 0,
      force = 2,
      force_pull = 0.5
    )
  } else {
    p <- p + geom_node_text(aes(label = label), repel = FALSE, size = label_size, check_overlap = TRUE)
  }

  p +
    scale_fill_manual(values = module_cols, name = "Module") +
    scale_edge_colour_manual(
      values = c(positive = "#4393C3", negative = "#D6604D"),
      name = "Edge"
    ) +
    scale_edge_alpha(range = c(0.15, 0.75), name = expression("|"*rho*"|")) +
    scale_edge_width(range = c(0.15, 1.2), guide = "none") +
    scale_size_identity(guide = "none") +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.margin = margin(t = 8, b = 5),
      legend.spacing.x = unit(0.8, "cm"),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7.5),
      legend.key.size = unit(0.4, "cm"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, hjust = 0, color = "grey40", margin = margin(b = 8)),
      plot.margin = margin(15, 15, 10, 15),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(
      fill = guide_legend(nrow = 1, override.aes = list(size = 5, alpha = 1)),
      edge_colour = guide_legend(override.aes = list(edge_width = 2, edge_alpha = 0.8)),
      edge_alpha = guide_legend(nrow = 1, override.aes = list(edge_width = 2))
    ) +
    labs(title = title, subtitle = subtitle)
}



save_grid <- function(plots, nrow, ncol, path, width, height, device = c("pdf", "tiff", "png")) {
  device <- match.arg(device)

  if (device == "pdf") {
    grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
  } else if (device == "png") {
    grDevices::png(path, width = width, height = height, units = "in", res = 300, bg = "white")
  } else {
    grDevices::tiff(path, width = width, height = height, units = "in", res = 300, compression = "lzw")
  }

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
  for (i in seq_along(plots)) {
    r <- ((i - 1) %/% ncol) + 1
    c <- ((i - 1) %% ncol) + 1
    print(plots[[i]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
  }
  grDevices::dev.off()
}

posneg_permutation <- function(pos_edges, neg_edges, n_perm = 10000, seed = 42) {
  total <- pos_edges + neg_edges
  if (total <= 0) stop("No edges for permutation test.", call. = FALSE)

  obs_prop <- pos_edges / total
  set.seed(seed)
  null_prop <- stats::rbinom(n_perm, size = total, prob = 0.5) / total
  p_two_sided <- mean(abs(null_prop - 0.5) >= abs(obs_prop - 0.5))

  list(
    total_edges = total,
    pos_edges = pos_edges,
    neg_edges = neg_edges,
    obs_pos_prop = obs_prop,
    p_two_sided = p_two_sided,
    null_prop = null_prop
  )
}

run_fig5 <- function(mode = c("all", "exudate", "microbe", "sensitivity", "posneg")) {
  mode <- match.arg(mode)

  message("Reading inputs...")
  exu <- read_exudate_matrix(data_exudate)
  mic <- read_microbe_matrix(
    data_microbe,
    rank = config$microbe_rank,
    top_n = config$microbe_top_n,
    min_prevalence = config$microbe_min_prevalence
  )

  if (mode %in% c("all", "exudate")) {
    message("Building exudate network (threshold = ", config$thr_exudate, ")...")
    exu_net <- build_correlation_network(
      exu$mat,
      threshold = config$thr_exudate,
      fdr_alpha = config$fdr_alpha,
      hub_top_n = config$hub_top_n,
      seed = config$seed
    )
    exu_topo <- network_topology(exu_net$graph, exu_net$cluster)
    exu_topo$dataset <- "exudate"
    exu_topo$threshold <- config$thr_exudate
    exu_topo$fdr_alpha <- config$fdr_alpha

    write.csv(exu_net$edges, file.path(out_dir, sprintf("Fig5A_exudate_edges_thr%.1f.csv", config$thr_exudate)), row.names = FALSE)
    write.csv(
      data.frame(
        node = V(exu_net$graph)$name,
        module = V(exu_net$graph)$module,
        degree = V(exu_net$graph)$degree,
        betweenness = V(exu_net$graph)$betweenness,
        is_hub = V(exu_net$graph)$is_hub
      ),
      file.path(out_dir, sprintf("Fig5A_exudate_nodes_thr%.1f.csv", config$thr_exudate)),
      row.names = FALSE
    )

    p_exu <- plot_network(
      exu_net$graph,
      title = "A. Exudate co-variance network",
      subtitle = sprintf("Spearman |rho| >= %.2f; BH-FDR < %.2f; n = %d nodes", 
                         config$thr_exudate, config$fdr_alpha, igraph::vcount(exu_net$graph)),
      node_size_attr = "degree",
      label_size = 2.5,
      filter_isolates = TRUE,
      layout_algo = "stress",
      seed = config$seed
    )

    # 输出 PNG 格式（更小且适合论文）
    ggsave(file.path(out_dir, "Fig5A_exudate_network.png"), p_exu, width = 9, height = 8, dpi = 300, bg = "white")
    # 也输出 PDF 版本
    ggsave(file.path(out_dir, "Fig5A_exudate_network.pdf"), p_exu, width = 9, height = 8, device = cairo_pdf)
  }

  if (mode %in% c("all", "microbe")) {
    message("Building microbe network (threshold = ", config$thr_microbe, ")...")
    mic_net <- build_correlation_network(
      mic$mat,
      threshold = config$thr_microbe,
      fdr_alpha = config$fdr_alpha,
      hub_top_n = config$hub_top_n,
      seed = config$seed
    )
    mic_topo <- network_topology(mic_net$graph, mic_net$cluster)
    mic_topo$dataset <- sprintf("microbe_%s", config$microbe_rank)
    mic_topo$threshold <- config$thr_microbe
    mic_topo$fdr_alpha <- config$fdr_alpha
    mic_topo$top_n_taxa <- config$microbe_top_n
    mic_topo$min_prevalence <- config$microbe_min_prevalence

    write.csv(mic_net$edges, file.path(out_dir, sprintf("Fig5B_microbe_edges_thr%.1f.csv", config$thr_microbe)), row.names = FALSE)
    write.csv(
      data.frame(
        node = V(mic_net$graph)$name,
        module = V(mic_net$graph)$module,
        degree = V(mic_net$graph)$degree,
        betweenness = V(mic_net$graph)$betweenness,
        is_hub = V(mic_net$graph)$is_hub
      ),
      file.path(out_dir, sprintf("Fig5B_microbe_nodes_thr%.1f.csv", config$thr_microbe)),
      row.names = FALSE
    )

    p_mic <- plot_network(
      mic_net$graph,
      title = "B. Microbiome co-occurrence network",
      subtitle = sprintf("Spearman |rho| >= %.2f; BH-FDR < %.2f; %s level (Top %d)", 
                         config$thr_microbe, config$fdr_alpha, config$microbe_rank, config$microbe_top_n),
      node_size_attr = "degree",
      label_size = 3.2,
      filter_isolates = FALSE,   # 微生物网络保留所有节点
      layout_algo = "stress",
      seed = config$seed
    )

    # 输出 PNG 格式（更小且适合论文）
    ggsave(file.path(out_dir, "Fig5B_microbe_network.png"), p_mic, width = 9, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(out_dir, "Fig5B_microbe_network.pdf"), p_mic, width = 9, height = 8, device = cairo_pdf)
  }

  if (mode == "all") {
    message("Saving combined main figure...")
    save_grid(
      plots = list(p_exu, p_mic),
      nrow = 1,
      ncol = 2,
      path = file.path(out_dir, "Fig5_main.pdf"),
      width = 15,
      height = 7,
      device = "pdf"
    )
    save_grid(
      plots = list(p_exu, p_mic),
      nrow = 1,
      ncol = 2,
      path = file.path(out_dir, "Fig5_main.tiff"),
      width = 15,
      height = 7,
      device = "tiff"
    )
    # PNG 格式（更小且适合论文投稿）
    save_grid(
      plots = list(p_exu, p_mic),
      nrow = 1,
      ncol = 2,
      path = file.path(out_dir, "Fig5_main.png"),
      width = 15,
      height = 7,
      device = "png"
    )

    exu_topo$top_n_taxa <- NA_integer_
    exu_topo$min_prevalence <- NA_real_
    exu_topo$feature_filter <- sprintf("%s_top%d", config$exudate_filter, config$exudate_top_n)
    mic_topo$feature_filter <- sprintf("%s_top%d_prev%.1f", config$microbe_rank, config$microbe_top_n, config$microbe_min_prevalence)

    all_cols <- union(names(exu_topo), names(mic_topo))

    align_cols <- function(df, cols) {
      missing <- setdiff(cols, names(df))
      for (m in missing) df[[m]] <- NA
      df[, cols, drop = FALSE]
    }

    exu_topo <- align_cols(exu_topo, all_cols)
    mic_topo <- align_cols(mic_topo, all_cols)
    topo <- rbind(exu_topo, mic_topo)
    topo <- topo[, c("dataset", "threshold", "fdr_alpha", setdiff(names(topo), c("dataset", "threshold", "fdr_alpha"))), drop = FALSE]
    write.csv(topo, file.path(out_dir, "Fig5_network_topology.csv"), row.names = FALSE)
  }

  if (mode %in% c("all", "sensitivity")) {
    message("Running threshold sensitivity (", config$thr_sensitivity_high, " vs ", config$thr_sensitivity_low, ")...")
    thr_pair <- c(config$thr_sensitivity_high, config$thr_sensitivity_low)

    exu_thr <- lapply(thr_pair, function(th) build_correlation_network(exu$mat, threshold = th, fdr_alpha = config$fdr_alpha, hub_top_n = config$hub_top_n, seed = config$seed))
    mic_thr <- lapply(thr_pair, function(th) build_correlation_network(mic$mat, threshold = th, fdr_alpha = config$fdr_alpha, hub_top_n = config$hub_top_n, seed = config$seed))

    exu_metrics <- do.call(rbind, lapply(seq_along(thr_pair), function(i) {
      m <- network_topology(exu_thr[[i]]$graph, exu_thr[[i]]$cluster)
      m$dataset <- "exudate"
      m$threshold <- thr_pair[i]
      m
    }))

    mic_metrics <- do.call(rbind, lapply(seq_along(thr_pair), function(i) {
      m <- network_topology(mic_thr[[i]]$graph, mic_thr[[i]]$cluster)
      m$dataset <- sprintf("microbe_%s", config$microbe_rank)
      m$threshold <- thr_pair[i]
      m
    }))

    hubs_overlap <- function(h1, h2) {
      inter <- length(intersect(h1, h2))
      uni <- length(union(h1, h2))
      data.frame(
        hub_top_n = config$hub_top_n,
        hub_overlap_n = inter,
        hub_jaccard = ifelse(uni == 0, NA_real_, inter / uni)
      )
    }

    exu_hub <- hubs_overlap(exu_thr[[1]]$hubs, exu_thr[[2]]$hubs)
    mic_hub <- hubs_overlap(mic_thr[[1]]$hubs, mic_thr[[2]]$hubs)

    sens_tbl <- rbind(
      data.frame(
        dataset = "exudate",
        thr_high = thr_pair[1],
        thr_low = thr_pair[2],
        modularity_high = igraph::modularity(exu_thr[[1]]$cluster),
        modularity_low = igraph::modularity(exu_thr[[2]]$cluster),
        avg_degree_high = mean(igraph::degree(exu_thr[[1]]$graph)),
        avg_degree_low = mean(igraph::degree(exu_thr[[2]]$graph)),
        edges_high = igraph::ecount(exu_thr[[1]]$graph),
        edges_low = igraph::ecount(exu_thr[[2]]$graph),
        nodes_high = igraph::vcount(exu_thr[[1]]$graph),
        nodes_low = igraph::vcount(exu_thr[[2]]$graph)
      ),
      data.frame(
        dataset = sprintf("microbe_%s", config$microbe_rank),
        thr_high = thr_pair[1],
        thr_low = thr_pair[2],
        modularity_high = igraph::modularity(mic_thr[[1]]$cluster),
        modularity_low = igraph::modularity(mic_thr[[2]]$cluster),
        avg_degree_high = mean(igraph::degree(mic_thr[[1]]$graph)),
        avg_degree_low = mean(igraph::degree(mic_thr[[2]]$graph)),
        edges_high = igraph::ecount(mic_thr[[1]]$graph),
        edges_low = igraph::ecount(mic_thr[[2]]$graph),
        nodes_high = igraph::vcount(mic_thr[[1]]$graph),
        nodes_low = igraph::vcount(mic_thr[[2]]$graph)
      )
    )
    sens_tbl$modularity_delta <- sens_tbl$modularity_low - sens_tbl$modularity_high
    sens_tbl <- cbind(
      sens_tbl,
      rbind(exu_hub, mic_hub)
    )
    write.csv(sens_tbl, file.path(out_dir, "FigS5_threshold_sensitivity_table.csv"), row.names = FALSE)

    p_exu_high <- plot_network(
      exu_thr[[1]]$graph,
      title = sprintf("Exudate |ρ| ≥ %.1f", thr_pair[1]),
      subtitle = sprintf("BH-FDR < %.2f", config$fdr_alpha),
      seed = config$seed
    )
    p_exu_low <- plot_network(
      exu_thr[[2]]$graph,
      title = sprintf("Exudate |ρ| ≥ %.1f", thr_pair[2]),
      subtitle = sprintf("BH-FDR < %.2f", config$fdr_alpha),
      seed = config$seed
    )
    p_mic_high <- plot_network(
      mic_thr[[1]]$graph,
      title = sprintf("Microbe (%s) |ρ| ≥ %.1f", config$microbe_rank, thr_pair[1]),
      subtitle = sprintf("BH-FDR < %.2f", config$fdr_alpha),
      seed = config$seed
    )
    p_mic_low <- plot_network(
      mic_thr[[2]]$graph,
      title = sprintf("Microbe (%s) |ρ| ≥ %.1f", config$microbe_rank, thr_pair[2]),
      subtitle = sprintf("BH-FDR < %.2f", config$fdr_alpha),
      seed = config$seed
    )

    save_grid(
      plots = list(p_exu_high, p_exu_low, p_mic_high, p_mic_low),
      nrow = 2,
      ncol = 2,
      path = file.path(out_dir, "FigS5_threshold_sensitivity.pdf"),
      width = 14,
      height = 12,
      device = "pdf"
    )
  }

  if (mode %in% c("all", "posneg")) {
    message("Running pos/neg edge permutation test...")
    if (!exists("exu_net", inherits = FALSE)) {
      exu_net <- build_correlation_network(
        exu$mat,
        threshold = config$thr_exudate,
        fdr_alpha = config$fdr_alpha,
        hub_top_n = config$hub_top_n,
        seed = config$seed
      )
    }
    if (!exists("mic_net", inherits = FALSE)) {
      mic_net <- build_correlation_network(
        mic$mat,
        threshold = config$thr_microbe,
        fdr_alpha = config$fdr_alpha,
        hub_top_n = config$hub_top_n,
        seed = config$seed
      )
    }

    exu_perm <- posneg_permutation(
      pos_edges = sum(E(exu_net$graph)$rho > 0),
      neg_edges = sum(E(exu_net$graph)$rho < 0),
      n_perm = config$perm_n,
      seed = config$seed
    )
    mic_perm <- posneg_permutation(
      pos_edges = sum(E(mic_net$graph)$rho > 0),
      neg_edges = sum(E(mic_net$graph)$rho < 0),
      n_perm = config$perm_n,
      seed = config$seed
    )

    res_tbl <- rbind(
      data.frame(
        dataset = "exudate",
        total_edges = exu_perm$total_edges,
        pos_edges = exu_perm$pos_edges,
        neg_edges = exu_perm$neg_edges,
        pos_neg_ratio = ifelse(exu_perm$neg_edges == 0, NA_real_, exu_perm$pos_edges / exu_perm$neg_edges),
        obs_pos_prop = exu_perm$obs_pos_prop,
        p_two_sided = exu_perm$p_two_sided
      ),
      data.frame(
        dataset = sprintf("microbe_%s", config$microbe_rank),
        total_edges = mic_perm$total_edges,
        pos_edges = mic_perm$pos_edges,
        neg_edges = mic_perm$neg_edges,
        pos_neg_ratio = ifelse(mic_perm$neg_edges == 0, NA_real_, mic_perm$pos_edges / mic_perm$neg_edges),
        obs_pos_prop = mic_perm$obs_pos_prop,
        p_two_sided = mic_perm$p_two_sided
      )
    )
    write.csv(res_tbl, file.path(out_dir, "posneg_test_results.csv"), row.names = FALSE)

    perm_df <- rbind(
      data.frame(dataset = "exudate", pos_prop = exu_perm$null_prop),
      data.frame(dataset = sprintf("microbe_%s", config$microbe_rank), pos_prop = mic_perm$null_prop)
    )

    obs_df <- res_tbl
    p_perm <- ggplot(perm_df, aes(x = pos_prop)) +
      geom_histogram(bins = 40, fill = "grey80", color = "grey30") +
      geom_vline(data = obs_df, aes(xintercept = obs_pos_prop), color = "#e31a1c", linewidth = 1) +
      facet_wrap(~dataset, scales = "free_y") +
      theme_classic(base_size = 12) +
      labs(
        title = "Permutation test for positive-edge proportion (null: p=0.5)",
        x = "Positive-edge proportion",
        y = "Count"
      )

    ggsave(file.path(out_dir, "FigS_posneg_permutation.pdf"), p_perm, width = 10, height = 5)
  }

  message("Done. Outputs saved to: ", out_dir)
  invisible(TRUE)
}

if (sys.nframe() == 0) {
  run_fig5("all")
}
