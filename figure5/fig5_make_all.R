# Figure 5 Redesign: Strict Readability Rules
# 
# Outputs:
# - Fig5_main.pdf/png/vector.pdf (180x120mm)
# - FigS5A_exudate_network.pdf (Full network, no labels)
# - FigS5B_microbe_network.pdf (Full network, high res)
#
# Logic:
# - A: Exudate -> Louvain -> Aggregated Modules (Top 8 + Other) -> Top 15 Inter-module Edges
# - B: Microbe -> Top 40 Genera -> Top 2 Edges/Node -> Top 12 Labelled Hubs
# - C: Topology Table (Exudate thr0.8, Microbe thr0.6)

options(stringsAsFactors = FALSE)

# --- Packages ---
suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(gridExtra)
  library(stringr)
  library(scales)
  library(ggrepel)
  library(cowplot)
  library(grid) # Essential for textGrob
})

# --- Configuration ---
config <- list(
  # Canvas
  width_mm = 200,    # Increased from 180
  height_mm = 135,   # Increased from 120
  base_size = 11,
  
  # Inputs
  thr_exudate = 0.8,
  thr_microbe = 0.6,
  
  # Panel A Settings
  A_top_modules = 8,
  A_top_edges = 15,
  
  # Panel B Settings
  B_top_nodes = 40,
  B_sparsify_k = 2,
  B_label_hubs = 12,
  B_top_modules = 6,
  
  # Seeds
  seed = 123
)

# --- Helpers ---

load_data <- function() {
  # Construct paths
  f_ex_edge <- list.files(pattern = paste0("Fig5A_exudate_edges_thr", config$thr_exudate), full.names = TRUE)[1]
  f_ex_node <- list.files(pattern = paste0("Fig5A_exudate_nodes_thr", config$thr_exudate), full.names = TRUE)[1]
  f_mic_edge <- list.files(pattern = paste0("Fig5B_microbe_edges_thr", config$thr_microbe), full.names = TRUE)[1]
  f_mic_node <- list.files(pattern = paste0("Fig5B_microbe_nodes_thr", config$thr_microbe), full.names = TRUE)[1]
  f_topo <- "Fig5_network_topology.csv"
  
  if (any(is.na(c(f_ex_edge, f_ex_node, f_mic_edge, f_mic_node)))) stop("Missing input files used pattern matching.")
  
  list(
    ex_edges = read.csv(f_ex_edge),
    ex_nodes = read.csv(f_ex_node),
    mic_edges = read.csv(f_mic_edge),
    mic_nodes = read.csv(f_mic_node),
    topo = read.csv(f_topo)
  )
}

  build_igraph <- function(nodes, edges) {
  # Standardize 'name' column
  if (!"name" %in% names(nodes)) {
    candidates <- c("node", "id", "OTU", "ASV", "Gene", "name")
    match <- intersect(candidates, names(nodes))
    if (length(match) > 0) {
      nodes$name <- nodes[[match[1]]]
    } else {
      nodes$name <- as.character(nodes[[1]])
    }
  }
  
  # igraph uses the first column of 'vertices' df as the symbolic data info.
  # But specifically it looks for 'name' attribute.
  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
  # Force name attribute if not set (sometimes graph_from_data_frame is tricky)
  if (is.null(V(g)$name)) {
     V(g)$name <- as.character(nodes$name)
  }
  
  E(g)$weight <- abs(E(g)$rho)
  E(g)$sign <- ifelse(E(g)$rho > 0, "Pos", "Neg")
  V(g)$degree <- degree(g)
  return(g)
}

# --- Panel A Logic ---

process_panel_A <- function(g) {
  # 1. Louvain Clustering (if not already strictly module attribute, recalc to be safe or use existing)
  # User file has 'module'. Let's trust it or re-run if needed. 
  # Request says "用 igraph ... Louvain 分模块". Let's re-run to ensure consistency or use existing if compatible.
  # Better to use existing 'module' column from CSV if valid, to match previous analyses.
  # Checking if 'module' exists in V(g).
  if (is.null(V(g)$module)) {
    cl <- cluster_louvain(g)
    V(g)$module <- membership(cl)
  }
  
  # 2. Module Stats & Top K
  mod_counts <- table(V(g)$module)
  top_mods <- names(sort(mod_counts, decreasing = TRUE))[1:min(config$A_top_modules, length(mod_counts))]
  
  # Relabel
  # Use local variable to preserve factor levels (igraph may strip them)
  mod_char <- as.character(V(g)$module)
  mod_lbl <- ifelse(mod_char %in% top_mods, paste0("M", mod_char), "Other")
  
  # Ensure "Other" is last factor level for color consistency
  lvl <- c(paste0("M", top_mods), "Other")
  mod_lbl <- factor(mod_lbl, levels = lvl)
  
  # Store character in graph for plotting later (if needed)
  V(g)$mod_label <- as.character(mod_lbl)
  
  # Console stats
  cat("\n[Panel A] Top 8 Modules Sizes:\n")
  print(mod_counts[top_mods])
  
  # 3. Aggregate
  # We need a new graph where nodes are modules
  # Edges are sum of counts
  
  # Extract edge list with source/target modules
  el <- as_edgelist(g, names = FALSE)
  
  # Map node indices to module labels using local factor
  m1 <- mod_lbl[el[,1]]
  m2 <- mod_lbl[el[,2]]
  
  edge_df <- data.frame(from = pmin(as.character(m1), as.character(m2)),
                        to = pmax(as.character(m1), as.character(m2)),
                        rho = E(g)$rho) %>%
    filter(from != to) %>% # Remove intra-module edges for the summary map
    group_by(from, to) %>%
    summarise(
      edge_count = n(),
      mean_abs_rho = mean(abs(rho)),
      .groups = "drop"
    ) %>%
    arrange(desc(edge_count), desc(mean_abs_rho))
  
  # Keep top edges
  top_edges <- head(edge_df, config$A_top_edges)
  
  cat("\n[Panel A] Top 15 Inter-Module Edges:\n")
  print(top_edges)
  
  # Build Summary Graph
  # Build Summary Graph
  # Nodes: Module names + Sizes
  mod_sizes <- table(mod_lbl) # Use local factor
  node_df <- data.frame(name = levels(mod_lbl))
  node_df$size <- as.numeric(mod_sizes[node_df$name])
  
  g_sum <- graph_from_data_frame(top_edges, vertices = node_df, directed = FALSE)
  
  return(g_sum)
}

plot_panel_A <- function(g_sum) {
  set.seed(config$seed)
  # Layout adjustment: Circle and Spread
  layout <- create_layout(g_sum, layout = "circle")
  layout$x <- layout$x * 1.6
  layout$y <- layout$y * 1.6
  
  p <- ggraph(layout) +
    geom_edge_link(aes(width = edge_count), color = "grey60", alpha = 0.7) +
    scale_edge_width_continuous(range = c(0.5, 3), name = "Edge Count") +
    geom_node_point(aes(size = size, fill = name), shape = 21, color = "white", stroke = 1.5, alpha = 0.95) +
    scale_size_continuous(range = c(6, 22), name = "Module Size") +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set2"), "grey80"), name = "Module") +
    geom_node_label(aes(label = name), fill = "white", alpha = 0.8, 
                    size = 4, fontface = "bold", label.size = 0) + # label.size=0 removes border
    theme_void(base_size = config$base_size) +
    labs(tag = "A", 
         title = "Exudates Co-occurrence (Module Summary)",
         subtitle = paste0("Summarized by Louvain modules (thr=", config$thr_exudate, "); full network in Fig. S5A")) +
    theme(
      plot.margin = margin(8, 8, 8, 8),
      plot.tag = element_text(face = "bold", size = 16),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      legend.position = "right",
      legend.key.size = unit(4, "mm")
    ) +
    guides(size = guide_legend(order = 1), fill = guide_legend(order = 2, override.aes = list(size=4)))
  
  return(p)
}

# --- Panel B Logic ---

process_panel_B <- function(g) {
  # 1. Top 40 Genera
  if (vcount(g) > config$B_top_nodes) {
    keep <- names(sort(degree(g), decreasing = TRUE))[1:config$B_top_nodes]
    g <- induced_subgraph(g, keep)
  }
  
  # 2. Sparsify Edges (Top-k per node)
  # Helper to find mask
  el <- as_edgelist(g, names = FALSE)
  w <- E(g)$weight
  
  # For each node 1..N, find incident edges, keep top k
  keep_edges <- rep(FALSE, ecount(g))
  
  # Use incident_edges which returns edge IDs
  incs <- incident_edges(g, V(g))
  
  for (i in seq_along(incs)) {
    e_ids <- as.numeric(incs[[i]])
    if (length(e_ids) > 0) {
      e_w <- w[e_ids]
      # indices in e_ids to keep
      top_local <- order(e_w, decreasing = TRUE)[1:min(length(e_w), config$B_sparsify_k)]
      keep_edges[e_ids[top_local]] <- TRUE
    }
  }
  
  # Use subgraph_from_edges (newer igraph)
  # If using older igraph, subgraph.edges is fine but warns.
  # Safer to check or just use subgraph.edges and suppress warning? 
  # Let's use specific function if available or fallback.
  # Use subgraph.edges for better compatibility across versions
  g_sparse <- subgraph.edges(g, which(keep_edges), delete.vertices = FALSE)
  
  # Recalc degree on sparse graph for sizing

  V(g_sparse)$deg_sparse <- degree(g_sparse)
  
  # 3. Top Labels
  # Robust sort and name extraction
  vals <- V(g_sparse)$deg_sparse
  names(vals) <- V(g_sparse)$name # Ensure named vector
  
  top_hubs <- names(sort(vals, decreasing = TRUE))[1:min(length(vals), config$B_label_hubs)]
  
  cat("\n[Panel B] Debug: First 5 names: ", head(V(g_sparse)$name), "\n")
  cat("[Panel B] Top 12 Hubs:\n")
  print(top_hubs)
  
  # Create label attribute
  if (length(top_hubs) > 0) {
    V(g_sparse)$label <- ifelse(V(g_sparse)$name %in% top_hubs, V(g_sparse)$name, NA)
  } else {
    V(g_sparse)$label <- NA
  }
  
  # 4. Modules (Recalc on sparse or use original?)
  # Request: "点颜色=模块（Louvain），模块图例不超过 6 个"
  # Let's run Louvain on this sparse graph for cleaner visuals
  cl <- cluster_louvain(g_sparse)
  mem <- membership(cl)
  
  mod_counts <- table(mem)
  top_mods <- names(sort(mod_counts, decreasing = TRUE))[1:min(config$B_top_modules, length(mod_counts))]
  
  V(g_sparse)$mod_col <- ifelse(as.character(mem) %in% top_mods, as.character(mem), "Other")
  V(g_sparse)$mod_col <- factor(V(g_sparse)$mod_col, levels = c(top_mods, "Other"))
  
  return(g_sparse)
}

plot_panel_B <- function(g) {
  set.seed(config$seed)
  # Layout: FR + Spread
  layout <- create_layout(g, layout = "fr")
  layout$x <- layout$x * 1.8
  layout$y <- layout$y * 1.8
  
  p <- ggraph(layout) +
    geom_edge_link(aes(color = sign), alpha = 0.65, width = 0.5) +
    scale_edge_color_manual(values = c("Pos" = "red", "Neg" = "blue"), name = "Sign") +
    geom_node_point(aes(size = deg_sparse, fill = mod_col), shape = 21, color = "white", stroke = 0.5, alpha = 0.95) +
    
    # Labeling Hubs
    geom_text_repel(aes(x = x, y = y, label = label), 
                    size = 3.6,
                    bg.color = "white", bg.r = 0.15,
                    box.padding = 0.35, point.padding = 0.25,
                    max.overlaps = Inf, seed = 42) +
                    
    scale_fill_brewer(palette = "Paired", name = "Module") +
    scale_size_continuous(range = c(2.5, 6.5), name = "Degree") +
    theme_void(base_size = config$base_size) +
    labs(tag = "B",
         title = "Microbiome Co-occurrence",
         subtitle = stringr::str_wrap(paste0("Top ", config$B_top_nodes, " genera (thr=", config$thr_microbe, "); top-", config$B_sparsify_k, " edges per node; labels show top-", config$B_label_hubs, " hubs"), 50)) +
    coord_cartesian(clip = "off") +
    theme(
      plot.tag = element_text(face = "bold", size = 16),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      legend.position = "right",
      plot.margin = margin(8, 40, 8, 8) 
    )
  return(p)
}

# --- Panel C Logic ---

plot_panel_C <- function(topo_df) {
  # Filter Data
  df <- topo_df %>%
    filter((grepl("exudate", dataset, ignore.case = TRUE) & threshold == config$thr_exudate) |
             (grepl("microbe", dataset, ignore.case = TRUE) & threshold == config$thr_microbe)) %>%
    select(Dataset = dataset, Nodes = nodes, Edges = edges, AvgDegree = avg_degree, 
           ModularityQ = modularity_Q, PosNeg = pos_neg_ratio) %>%
    mutate(
      Dataset = ifelse(grepl("exudate", Dataset, ignore.case = TRUE), "Exudate", "Microbe"),
      AvgDegree = round(AvgDegree, 2),
      ModularityQ = round(ModularityQ, 3),
      PosNeg = round(PosNeg, 2)
    )
  
  # Format as Grob
  tt <- ttheme_default(
    core = list(fg_params = list(fontsize = 9, hjust=0, x=0.1),
                bg_params = list(fill = c("grey95", "white"))),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold", hjust=0, x=0.1),
                   bg_params = list(fill = "grey80")),
    padding = unit(c(4, 4), "mm")
  )
  
  tbl <- tableGrob(df, rows = NULL, theme = tt)
  
  # Wrap in ggplot to align with patchwork if needed, or return grob
  # Patchwork handles grobs via wrap_elements
  return(tbl)
}

# --- Supplementary Plots ---

plot_S5A_full <- function(g) {
  # Full Exudate: No labels, just structure
  cl <- cluster_louvain(g)
  V(g)$module <- as.factor(membership(cl))
  E(g)$color <- ifelse(E(g)$rho > 0, "red", "blue")
  
  set.seed(config$seed)
  ggraph(g, layout = "nicely") +
    geom_edge_link(aes(color = sign), alpha = 0.15, width = 0.3) +
    geom_node_point(aes(color = module), size = 1.5, alpha = 0.8) +
    scale_edge_color_manual(values = c("Pos"="red", "Neg"="blue")) +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = "Figure S5A: Full Exudate Network",
         subtitle = paste0("1686 Nodes, Thr=", config$thr_exudate, ". Labels omitted for readability."))
}

plot_S5B_full <- function(g) {
  # Full Microbe: Labels allowed on Hubs
  cl <- cluster_louvain(g)
  V(g)$module <- as.factor(membership(cl))
  V(g)$label <- ifelse(degree(g) > quantile(degree(g), 0.90), V(g)$name, NA)
  
  set.seed(config$seed)
  ggraph(g, layout = "nicely") +
    geom_edge_link(aes(color = sign), alpha = 0.2, width = 0.3) +
    geom_node_point(aes(color = module, size = degree), alpha = 0.9) +
    geom_text_repel(aes(x = x, y = y, label = label), size = 2, max.overlaps = 30) +
    scale_edge_color_manual(values = c("Pos"="red", "Neg"="blue")) +
    scale_size_continuous(range = c(1, 5)) +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = "Figure S5B: Full Microbiome Network",
         subtitle = paste0("Thr=", config$thr_microbe, ". Top 10% hubs labeled."))
}


# --- Main ---

main <- function() {
  cat("Loading data...\n")
  data <- load_data()
  
  # Build base graphs
  g_ex <- build_igraph(data$ex_nodes, data$ex_edges)
  g_mic <- build_igraph(data$mic_nodes, data$mic_edges)
  
  # --- Panel A ---
  cat("Processing Panel A...\n")
  g_ex_sum <- process_panel_A(g_ex)
  pA <- plot_panel_A(g_ex_sum)
  
  # --- Panel B ---
  cat("Processing Panel B...\n")
  g_mic_sparse <- process_panel_B(g_mic)
  pB <- plot_panel_B(g_mic_sparse)
  
  # --- Panel C ---
  cat("Processing Panel C...\n")
  pC <- plot_panel_C(data$topo)
  
  # --- Split Layout & Export ---
  cat("Creating Separated Panels (Clip Off)...\n")
  
  # 1. Prepare Individual Panels (Custom margins for standalone output)
  pA_only <- pA + 
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  pB_only <- pB + 
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(5, 35, 5, 5), # Extra right margin for labels/legend
      legend.position = "right"
    )
  
  # Ensure pC is a ggplot object for modifications if it was a Grob
  # pC passed here is a gtable/grob from tableGrob. 
  # We need to wrap it or just save it. tableGrob doesn't take theme(plot.margin).
  # But we can use grid.arrange or wrap_elements.
  # For consistent export, let's wrap it in wrap_elements (patchwork) which returns a ggplot.
  pC_only <- wrap_elements(pC) + 
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # 2. Export Separate Files
  files_to_check <- c()
  
  # Fig5A
  f_a_pdf <- "Fig5A_module_summary.pdf"
  f_a_png <- "Fig5A_module_summary.png"
  ggsave(f_a_pdf, pA_only, width = 100, height = 80, units = "mm", device = "pdf", limitsize = FALSE)
  ggsave(f_a_png, pA_only, width = 100, height = 80, units = "mm", dpi = 600, limitsize = FALSE)
  files_to_check <- c(files_to_check, f_a_pdf, f_a_png)
  
  # Fig5B
  f_b_pdf <- "Fig5B_microbiome_network.pdf"
  f_b_png <- "Fig5B_microbiome_network.png"
  ggsave(f_b_pdf, pB_only, width = 140, height = 80, units = "mm", device = "pdf", limitsize = FALSE)
  ggsave(f_b_png, pB_only, width = 140, height = 80, units = "mm", dpi = 600, limitsize = FALSE)
  files_to_check <- c(files_to_check, f_b_pdf, f_b_png)
  
  # Fig5C
  f_c_pdf <- "Fig5C_network_metrics.pdf"
  f_c_png <- "Fig5C_network_metrics.png"
  ggsave(f_c_pdf, pC_only, width = 120, height = 45, units = "mm", device = "pdf", limitsize = FALSE)
  ggsave(f_c_png, pC_only, width = 120, height = 45, units = "mm", dpi = 600, limitsize = FALSE)
  files_to_check <- c(files_to_check, f_c_pdf, f_c_png)
  
  # 3. Export Combined Big Version
  cat("Assembling Combined Big Figure...\n")
  # Use the _only versions but with layout
  # Need to ensure title C is included if we want it.
  
  # Re-make Title C for big plot
  title_C_grob <- textGrob("C", x = unit(0.01, "npc"), y = unit(0.9, "npc"), 
                      just = "left", gp = gpar(fontsize = 16, fontface = "bold"))
  pC_grob_big <- arrangeGrob(title_C_grob, pC, ncol = 1, heights = c(0.15, 0.85))
  
  final_big <- (pA_only + pB_only + plot_layout(widths = c(1, 1.45))) /
               wrap_elements(pC_grob_big) + 
               plot_layout(heights = c(1, 0.4))
               
  ggsave("Fig5_main.pdf", final_big, width = 220, height = 140, units = "mm", device = "pdf", limitsize = FALSE)
  ggsave("Fig5_main.vector.pdf", final_big, width = 220, height = 140, units = "mm", device = "pdf", limitsize = FALSE)
  ggsave("Fig5_main.png", final_big, width = 220, height = 140, units = "mm", dpi = 600, limitsize = FALSE)
  files_to_check <- c(files_to_check, "Fig5_main.pdf", "Fig5_main.vector.pdf", "Fig5_main.png")
  
  # --- Supplements ---
  cat("Generating Supplements...\n")
  
  pS5A <- plot_S5A_full(g_ex)
  ggsave("FigS5A_exudate_network.pdf", pS5A, width = 10, height = 10, device = "pdf")
  files_to_check <- c(files_to_check, "FigS5A_exudate_network.pdf")
  
  pS5B <- plot_S5B_full(g_mic)
  ggsave("FigS5B_microbe_network.pdf", pS5B, width = 12, height = 12, device = "pdf")
  files_to_check <- c(files_to_check, "FigS5B_microbe_network.pdf")
  
  # --- Verification Output ---
  cat("\n=== Verification ===\n")
  cat("1. Output Files:\n")
  for (f in files_to_check) {
    if (file.exists(f)) {
      size_mb <- file.size(f) / 1024 / 1024
      cat(sprintf("   - %-30s : %.2f MB\n", f, size_mb))
    } else {
      cat(sprintf("   - %-30s : [MISSING]\n", f))
    }
  }
  
  cat("\n2. Topology Metrics (Panel C Content):\n")
  # Filter again to show exactly what's in the table
  df_check <- data$topo %>%
    filter((grepl("exudate", dataset, ignore.case = TRUE) & threshold == config$thr_exudate) |
             (grepl("microbe", dataset, ignore.case = TRUE) & threshold == config$thr_microbe)) %>%
    select(Dataset = dataset, Nodes = nodes, Edges = edges, AvgDegree = avg_degree, 
           Modularity = modularity_Q, PosNeg=pos_neg_ratio)
  print(df_check)
  
  cat("All tasks completed.\n")
}

# Run
if (!interactive()) main() else main()
