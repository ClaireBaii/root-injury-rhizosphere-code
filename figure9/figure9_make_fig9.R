# Figure 9: sPLS-DA (scores + loadings) + chord (Spearman correlations)
# R version: 4.4.3
#
# Methods (per ToDO_fig9.md)
# - sPLS-DA: mixOmics::splsda + tune.splsda (keepX/ncomp) + perf (CV performance)
# - Correlation for chord: Spearman + BH-FDR; draw links with |rho| > 0.6 & q < 0.05
# - Chord plotting: circlize::chordDiagram (+ log2FC heat ring)
#
# Inputs (repo root /data)
# - 完整数据-分泌物.csv (metabolome; sample columns: YS0_F-1 ... YS100_F-3)
# - 完整数据-微生物.csv (microbiome; sample columns: RDYS_1-1 ... RDYS_5-3 + taxonomy)
#
# Outputs (repo root /figure9)
# - Fig9_main.pdf + Fig9_main.tiff
# - FigS9_splsda_CV.pdf
# - splsda_tuning_results.csv
# - FigS9_selected_features.csv (15 metabolites + 12 taxa)
# - FigS9_chord_links.csv (edges used in chord)

## -------------------- config --------------------
set.seed(123)

tax_rank <- "phylum" # "phylum" or "genus" (match the manuscript)

# sPLS-DA tuning / CV
ncomp_max <- 2
keepX_grid <- c(5, 10, 15, 20, 25, 30)
cv_folds <- 3
cv_repeats <- 50
cv_dist <- "centroids.dist"
cv_measure <- "BER"

# Loadings display (main Fig.9)
loadings_top_n <- 10
loadings_label_cex <- 1.0

# Feature list for Figure 9B / chord
n_select_metabolites <- 15
n_select_taxa <- 12

# Correlation filtering for chord (per manuscript)
rho_threshold <- 0.6
q_threshold <- 0.05

# Chord aesthetics
highlight_abs_rho <- 0.8
label_cex <- 0.85

## -------------------- helpers --------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE)
}

find_project_root <- function(start_dir) {
  candidate <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in 0:10) {
    exu <- file.path(candidate, "data", "完整数据-分泌物.csv")
    mic <- file.path(candidate, "data", "完整数据-微生物.csv")
    if (file.exists(exu) && file.exists(mic)) return(candidate)

    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  stop(
    "Cannot find project root containing data/完整数据-分泌物.csv and data/完整数据-微生物.csv. ",
    "Start: ", start_dir,
    call. = FALSE
  )
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Missing package '%s'. Please install it first.", pkg), call. = FALSE)
  }
}

read_csv_any <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE))
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

map_rdys_to_ys <- function(rdys_sample) {
  m <- regexec("^RDYS_(\\d+)-(\\d+)$", rdys_sample)
  parts <- regmatches(rdys_sample, m)[[1]]
  if (length(parts) == 0) return(NA_character_)

  group_idx <- as.integer(parts[2])
  rep_idx <- as.integer(parts[3])
  ratios <- c(0, 25, 50, 75, 100)

  if (is.na(group_idx) || group_idx < 1 || group_idx > length(ratios)) return(NA_character_)
  if (is.na(rep_idx) || rep_idx < 1) return(NA_character_)
  sprintf("YS%d_F-%d", ratios[group_idx], rep_idx)
}

extract_tax_rank <- function(taxonomy, rank = c("phylum", "genus")) {
  rank <- match.arg(rank)
  prefix <- if (rank == "phylum") "p__" else "g__"

  if (is.na(taxonomy) || !nzchar(taxonomy)) return("Unclassified")
  m <- regmatches(taxonomy, regexpr(paste0("(?:^|;\\s*)", prefix, "[^;]+"), taxonomy, perl = TRUE))
  if (length(m) == 0 || !nzchar(m)) return("Unclassified")

  val <- sub(paste0("^.*", prefix), "", m)
  val <- trimws(val)
  if (!nzchar(val)) return("Unclassified")
  val
}

parse_ratio_from_ys <- function(ys_sample) {
  suppressWarnings(as.integer(sub("^YS(\\d+)_.*$", "\\1", ys_sample)))
}

make_group_ratio_5 <- function(sample_ids) {
  ratio <- vapply(sample_ids, parse_ratio_from_ys, integer(1))
  if (any(is.na(ratio))) {
    bad <- sample_ids[is.na(ratio)]
    stop("Failed to parse YS ratio from sample ids: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (any(!ratio %in% c(0, 25, 50, 75, 100))) {
    bad <- sample_ids[!ratio %in% c(0, 25, 50, 75, 100)]
    stop("Unexpected YS ratio (expect 0/25/50/75/100) for samples: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  factor(
    ratio,
    levels = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%")
  )
}

read_metabolome <- function(path) {
  df <- read_csv_any(path)
  if (!("Name" %in% names(df))) stop("Metabolome file missing column 'Name': ", path, call. = FALSE)

  sample_cols <- grep("^YS\\d+_F-\\d+$", names(df), value = TRUE)
  if (length(sample_cols) == 0) stop("No metabolome sample columns like 'YS0_F-1' in: ", path, call. = FALSE)

  met_name <- as.character(df$Name)
  met_name[is.na(met_name) | !nzchar(met_name)] <- "Unknown"
  met_id <- make.unique(met_name)

  raw_mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(raw_mat) <- "double"
  raw_mat[!is.finite(raw_mat)] <- 0
  rownames(raw_mat) <- met_id

  log_mat <- log2(raw_mat + 1)

  sds <- apply(log_mat, 1, sd)
  keep <- is.finite(sds) & (sds > 0)
  log_mat <- log_mat[keep, , drop = FALSE]
  raw_mat <- raw_mat[keep, , drop = FALSE]

  x_log <- t(log_mat) # samples x metabolites
  x_raw <- t(raw_mat)
  rownames(x_log) <- sample_cols
  rownames(x_raw) <- sample_cols

  meta <- df[keep, c("Name", "Super_Class", "Class", "Sub_Class"), drop = FALSE]
  meta$Feature <- rownames(log_mat)
  rownames(meta) <- meta$Feature
  list(x_log = x_log, x_raw = x_raw, meta = meta)
}

read_microbiome <- function(path, rank = c("phylum", "genus"), min_mean_rel = 1e-4, max_features = 200) {
  rank <- match.arg(rank)
  df <- read_csv_any(path)
  if (!("taxonomy" %in% names(df))) stop("Microbiome file missing column 'taxonomy': ", path, call. = FALSE)

  sample_cols <- grep("^RDYS_\\d+-\\d+$", names(df), value = TRUE)
  if (length(sample_cols) == 0) stop("No microbiome sample columns like 'RDYS_1-1' in: ", path, call. = FALSE)

  counts <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(counts) <- "double"
  counts[!is.finite(counts)] <- 0

  tax <- vapply(df$taxonomy, extract_tax_rank, character(1), rank = rank)
  agg <- rowsum(counts, group = tax, na.rm = TRUE, reorder = FALSE)
  agg <- agg[rownames(agg) != "", , drop = FALSE]

  col_totals <- colSums(agg)
  if (any(col_totals == 0)) stop("At least one microbiome sample has zero total counts.", call. = FALSE)
  rel <- sweep(agg, 2, col_totals, "/")
  rel[!is.finite(rel)] <- 0

  mean_rel <- rowMeans(rel)
  keep <- is.finite(mean_rel) & (mean_rel > 0)
  rel <- rel[keep, , drop = FALSE]
  mean_rel <- mean_rel[keep]

  if (!is.null(min_mean_rel) && is.finite(min_mean_rel)) {
    keep2 <- mean_rel >= min_mean_rel
    rel <- rel[keep2, , drop = FALSE]
    mean_rel <- mean_rel[keep2]
  }
  if (!is.null(max_features) && is.finite(max_features) && nrow(rel) > max_features) {
    ord <- order(mean_rel, decreasing = TRUE)
    rel <- rel[ord[seq_len(max_features)], , drop = FALSE]
  }

  # map RDYS_* sample names to YS*_F-* (match metabolome)
  mapped <- vapply(colnames(rel), map_rdys_to_ys, character(1))
  if (any(is.na(mapped))) {
    bad <- colnames(rel)[is.na(mapped)]
    stop("Failed to map some RDYS sample names to YS-style: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(mapped)) stop("Mapped microbiome sample names are not unique.", call. = FALSE)
  colnames(rel) <- mapped

  x_rel <- t(rel) # samples x taxa
  x_log <- log2(x_rel + 1e-6)

  sds <- apply(x_log, 2, sd)
  keep_var <- is.finite(sds) & (sds > 0)
  x_log <- x_log[, keep_var, drop = FALSE]

  list(x_log = x_log)
}

select_top_by_loading <- function(model, met_features, tax_features, met_meta, n_met = 15, n_tax = 12) {
  load <- model$loadings$X
  if (is.null(load) || nrow(load) == 0) stop("Model loadings are empty; cannot select features.", call. = FALSE)

  abs_load <- abs(load)
  importance <- apply(abs_load, 1, max)
  comp_max <- apply(abs_load, 1, which.max)
  loading_at_max <- vapply(seq_along(importance), function(i) load[i, comp_max[i]], numeric(1))

  tbl <- data.frame(
    Feature = rownames(load),
    Importance = importance,
    Component = comp_max,
    Loading = loading_at_max,
    stringsAsFactors = FALSE
  )
  tbl$Type <- ifelse(tbl$Feature %in% met_features, "Metabolite",
    ifelse(tbl$Feature %in% tax_features, "Taxon", "Other")
  )

  met_tbl <- tbl[tbl$Type == "Metabolite", , drop = FALSE]
  tax_tbl <- tbl[tbl$Type == "Taxon", , drop = FALSE]
  met_tbl <- met_tbl[order(met_tbl$Importance, decreasing = TRUE), , drop = FALSE]
  tax_tbl <- tax_tbl[order(tax_tbl$Importance, decreasing = TRUE), , drop = FALSE]

  met_tbl <- met_tbl[seq_len(min(n_met, nrow(met_tbl))), , drop = FALSE]
  tax_tbl <- tax_tbl[seq_len(min(n_tax, nrow(tax_tbl))), , drop = FALSE]

  if (!is.null(met_meta) && nrow(met_meta) > 0) {
    met_tbl$Super_Class <- met_meta[met_tbl$Feature, "Super_Class"]
    met_tbl$Class <- met_meta[met_tbl$Feature, "Class"]
    met_tbl$Sub_Class <- met_meta[met_tbl$Feature, "Sub_Class"]
    
    tax_tbl$Super_Class <- NA_character_
    tax_tbl$Class <- NA_character_
    tax_tbl$Sub_Class <- NA_character_
  }
  out <- rbind(met_tbl, tax_tbl)
  out$Rank <- ave(out$Importance, out$Type, FUN = function(x) rank(-x, ties.method = "first"))
  out <- out[order(out$Type, out$Rank), , drop = FALSE]
  rownames(out) <- NULL
  out
}

make_loadings_long <- function(model, met_features, tax_features, met_meta = NULL, zero_tol = 1e-12) {
  load <- model$loadings$X
  if (is.null(load) || nrow(load) == 0) stop("Model loadings are empty; cannot export loadings table.", call. = FALSE)

  comps <- seq_len(ncol(load))
  out <- vector("list", length(comps))
  names(out) <- paste0("Comp", comps)

  for (j in comps) {
    v <- load[, j]
    ok <- is.finite(v) & (abs(v) > zero_tol)
    if (!any(ok)) next

    df <- data.frame(
      Feature = rownames(load)[ok],
      Component = j,
      Loading = unname(v[ok]),
      AbsLoading = abs(unname(v[ok])),
      stringsAsFactors = FALSE
    )
    df <- df[order(df$AbsLoading, decreasing = TRUE), , drop = FALSE]
    df$RankAbs <- seq_len(nrow(df))
    df$Type <- ifelse(df$Feature %in% met_features, "Metabolite",
      ifelse(df$Feature %in% tax_features, "Taxon", "Other")
    )

    if (!is.null(met_meta) && nrow(met_meta) > 0) {
      df$Super_Class <- NA_character_
      df$Class <- NA_character_
      df$Sub_Class <- NA_character_

      met_rows <- df$Type == "Metabolite" & df$Feature %in% rownames(met_meta)
      if (any(met_rows)) {
        df$Super_Class[met_rows] <- met_meta[df$Feature[met_rows], "Super_Class"]
        df$Class[met_rows] <- met_meta[df$Feature[met_rows], "Class"]
        df$Sub_Class[met_rows] <- met_meta[df$Feature[met_rows], "Sub_Class"]
      }
    }

    out[[j]] <- df
  }

  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) {
    out <- data.frame(
      Feature = character(0),
      Component = integer(0),
      Loading = numeric(0),
      AbsLoading = numeric(0),
      RankAbs = integer(0),
      Type = character(0),
      Super_Class = character(0),
      Class = character(0),
      Sub_Class = character(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out
}

compute_links <- function(met_mat, tax_mat, rho_cutoff = 0.6, q_cutoff = 0.05) {
  metabolites <- colnames(met_mat)
  taxa <- colnames(tax_mat)
  grid <- expand.grid(Metabolite = metabolites, Taxon = taxa, stringsAsFactors = FALSE)

  rho <- numeric(nrow(grid))
  pval <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    m <- grid$Metabolite[i]
    t <- grid$Taxon[i]
    x <- met_mat[, m]
    y <- tax_mat[, t]

    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) {
      rho[i] <- NA_real_
      pval[i] <- NA_real_
      next
    }
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    rho[i] <- unname(ct$estimate)
    pval[i] <- ct$p.value
  }

  grid$rho <- rho
  grid$p <- pval
  grid$q <- p.adjust(grid$p, method = "BH")

  keep <- is.finite(grid$rho) & is.finite(grid$q) & (abs(grid$rho) > rho_cutoff) & (grid$q < q_cutoff)
  grid[keep, , drop = FALSE]
}

make_link_colors <- function(rho_vec) {
  pos_pal <- grDevices::colorRampPalette(c("#FEE5D9", "#A50F15"))(100)
  neg_pal <- grDevices::colorRampPalette(c("#DEEBF7", "#08519C"))(100)

  cols <- character(length(rho_vec))
  for (i in seq_along(rho_vec)) {
    r <- rho_vec[i]
    s <- min(1, max(0, abs(r)))
    idx <- max(1, min(100, floor(s * 99) + 1))
    cols[i] <- if (is.na(r)) "grey80" else if (r >= 0) pos_pal[idx] else neg_pal[idx]
  }
  cols
}

plot_chord <- function(links, sector_order, grid_col, fc_vec, title_text, n_met, class_colors) {
  if (nrow(links) == 0) stop("No chord links to plot (after filtering).", call. = FALSE)

  require_pkg("circlize")
  suppressPackageStartupMessages(library(circlize))

  circos.clear()
  
  # Get actual sectors from links data
  actual_sectors <- unique(c(links$Metabolite, links$Taxon))
  
  if (length(actual_sectors) < 2) {
    stop("Too few sectors for chord diagram.", call. = FALSE)
  }
  
  circos.par(
    start.degree = 90,
    track.margin = c(0.005, 0.005),
    canvas.xlim = c(-1.6, 1.6),
    canvas.ylim = c(-1.25, 1.25)
  )

  link_cols <- make_link_colors(links$rho)
  link_lwd <- 0.6 + 2.2 * pmax(0, abs(links$rho) - rho_threshold) / (1 - rho_threshold)
  link_lwd[!is.finite(link_lwd)] <- 0.6
  link_lwd[abs(links$rho) >= highlight_abs_rho] <- link_lwd[abs(links$rho) >= highlight_abs_rho] * 1.6

  # Filter grid.col to only include actual sectors
  grid_col_filtered <- grid_col[names(grid_col) %in% actual_sectors]
  
  chordDiagram(
    x = links[, c("Metabolite", "Taxon", "rho")],
    grid.col = grid_col_filtered,
    col = link_cols,
    link.lwd = link_lwd,
    transparency = 0.12,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12),
    reduce = 0
  )

  # sector labels (bigger fonts)
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      nm <- get.cell.meta.data("sector.index")
      circos.text(
        x = CELL_META$xcenter,
        y = CELL_META$ylim[1] + mm_y(2.2),
        labels = nm,
        facing = "clockwise",
        niceFacing = TRUE,
        cex = label_cex
      )
    },
    bg.border = NA
  )

  # log2FC heat ring (Processing vs Control; metabolites only)
  fc_lim <- max(2, max(abs(fc_vec), na.rm = TRUE))
  fc_colors <- colorRamp2(c(-fc_lim, 0, fc_lim), c("#2C7BB6", "white", "#D7191C"))
  circos.track(
    ylim = c(0, 1),
    track.height = 0.06,
    bg.border = NA,
    panel.fun = function(x, y) {
      nm <- CELL_META$sector.index
      value <- fc_vec[nm]
      if (is.finite(value)) {
        circos.rect(
          CELL_META$xlim[1], 0,
          CELL_META$xlim[2], 1,
          col = fc_colors(value),
          border = NA
        )
      }
    }
  )

  title(title_text, line = -1)

  # Legends
  par(xpd = TRUE)
  legend(
    x = 1.05, y = 1.15,
    legend = c(
      sprintf("Positive (rho > %.1f)", rho_threshold),
      sprintf("Negative (rho < -%.1f)", rho_threshold),
      sprintf("Thicker: |rho| >= %.1f", highlight_abs_rho)
    ),
    col = c("#A50F15", "#08519C", "grey30"),
    lty = 1,
    lwd = c(2, 2, 3),
    bty = "n",
    cex = 0.85
  )
  legend(
    x = 1.05, y = 0.75,
    legend = c(sprintf("-%.1f", fc_lim), "0", sprintf("+%.1f", fc_lim)),
    fill = c(fc_colors(-fc_lim), fc_colors(0), fc_colors(fc_lim)),
    title = "log2FC\n(mean 25-100% vs 0%)",
    bty = "n",
    cex = 0.75
  )

  if (!is.null(class_colors) && length(class_colors) > 0) {
    legend(
      x = 1.05, y = 0.3,
      legend = names(class_colors),
      fill = unname(class_colors),
      border = NA,
      title = "Metabolite class",
      bty = "n",
      cex = 0.65
    )
    legend(
      x = 1.05, y = -0.65,
      legend = "Taxon",
      fill = "grey80",
      border = NA,
      bty = "n",
      cex = 0.65
    )
  }

  circos.clear()
}

## -------------------- packages --------------------
require_pkg("mixOmics")
require_pkg("circlize")

suppressPackageStartupMessages({
  library(mixOmics)
})

## -------------------- mixOmics sandbox compatibility --------------------
# In some sandboxed environments, `parallel::detectCores()` can return NA
# (e.g. sysctl blocked on macOS), which triggers an NA comparison inside
# mixOmics' internal `.check_cpus()` and aborts tuning/perf.
if (is.na(suppressWarnings(parallel::detectCores()))) {
  message("Note: parallel::detectCores() returned NA; patching mixOmics .check_cpus() for this session (cpus=1).")
  assignInNamespace(
    ".check_cpus",
    function(cpus) {
      if (is.null(cpus) || length(cpus) != 1 || !is.numeric(cpus) || is.na(cpus) || cpus <= 0) {
        stop("Number of CPUs should be a positive integer.", call. = FALSE)
      }
      cores <- suppressWarnings(parallel::detectCores())
      if (is.na(cores) || !is.finite(cores) || cores < 1) cores <- as.integer(cpus)
      if (cpus > cores) {
        message(sprintf("\nOnly %s CPUs available for parallel processing.\n", cores))
      }
      as.integer(cpus)
    },
    ns = "mixOmics"
  )
}

## -------------------- IO paths --------------------
script_dir <- get_script_dir()
project_root <- find_project_root(script_dir)
data_dir <- file.path(project_root, "data")
out_dir <- file.path(project_root, "figure9")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

exudate_path <- file.path(data_dir, "完整数据-分泌物.csv")
microbe_path <- file.path(data_dir, "完整数据-微生物.csv")

## -------------------- read + preprocess --------------------
met <- read_metabolome(exudate_path)
mic <- read_microbiome(microbe_path, rank = tax_rank)

common_samples <- intersect(rownames(met$x_log), rownames(mic$x_log))
if (length(common_samples) < 6) stop("Too few matched samples between metabolome and microbiome: ", length(common_samples), call. = FALSE)
common_samples <- sort(common_samples)

X_met <- met$x_log[common_samples, , drop = FALSE]
X_tax <- mic$x_log[common_samples, , drop = FALSE]
X_met_raw <- met$x_raw[common_samples, , drop = FALSE]

# Combine for sPLS-DA (metabolome + aggregated microbiome)
X <- cbind(X_met, X_tax)
Y <- make_group_ratio_5(rownames(X))
message("Group counts: ", paste(names(table(Y)), as.integer(table(Y)), sep = "=", collapse = ", "))

met_features <- colnames(X_met)
tax_features <- colnames(X_tax)

# sanity check: "citrate" naming consistency (Figure 7/9 note)
met_names_lower <- tolower(met$meta$Name)
has_citrate <- any(grepl("\\bcitrate\\b|\\bcitric acid\\b", met_names_lower, perl = TRUE))
if (!has_citrate) {
  message("Note: No metabolite annotated as 'citrate/citric acid' found in data/完整数据-分泌物.csv. ",
          "If manuscript text mentions citrate, please verify and correct accordingly.")
}

## -------------------- tune + fit + perf --------------------
keepX_grid <- unique(pmax(1, pmin(keepX_grid, ncol(X) - 1)))
keepX_grid <- keepX_grid[order(keepX_grid)]

tune <- tune.splsda(
  X,
  Y,
  ncomp = ncomp_max,
  validation = "Mfold",
  folds = cv_folds,
  nrepeat = cv_repeats,
  dist = cv_dist,
  measure = cv_measure,
  test.keepX = keepX_grid,
  progressBar = TRUE
)

if (!is.null(tune$choice.ncomp)) {
  if (is.list(tune$choice.ncomp) && "ncomp" %in% names(tune$choice.ncomp)) {
    best_ncomp <- tune$choice.ncomp$ncomp[1]
  } else {
    best_ncomp <- tune$choice.ncomp
  }
} else {
  best_ncomp <- ncomp_max
}

# Force ncomp=2 for visualization
if (best_ncomp < 2) {
  best_ncomp <- 2
}

best_keepX <- tune$choice.keepX
if (is.list(best_keepX)) best_keepX <- unlist(best_keepX)
if (length(best_keepX) < best_ncomp) {
  best_keepX <- c(best_keepX, rep(tail(best_keepX, 1), best_ncomp - length(best_keepX)))
}
best_keepX <- best_keepX[seq_len(best_ncomp)]

spls_model <- splsda(X, Y, ncomp = best_ncomp, keepX = best_keepX, scale = TRUE)

perf_res <- perf(
  spls_model,
  validation = "Mfold",
  folds = cv_folds,
  nrepeat = cv_repeats,
  dist = cv_dist,
  progressBar = TRUE
)

## -------------------- outputs: CV + tuning table --------------------
pdf(file.path(out_dir, "FigS9_splsda_CV.pdf"), width = 10, height = 8, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))
try({
  plot(tune)
  mtext(sprintf("M-fold CV: folds=%d, repeats=%d, measure=%s, dist=%s", cv_folds, cv_repeats, cv_measure, cv_dist), side = 3, line = 0.5, cex = 0.8)
}, silent = TRUE)
try({
  plot(perf_res)
  mtext(sprintf("M-fold CV: folds=%d, repeats=%d, dist=%s", cv_folds, cv_repeats, cv_dist), side = 3, line = 0.5, cex = 0.8)
}, silent = TRUE)
dev.off()

# Detailed tuning results table
tuning_tbl <- data.frame(
  choice_ncomp = best_ncomp,
  component = seq_len(best_ncomp),
  keepX = best_keepX,
  folds = cv_folds,
  repeats = cv_repeats,
  dist = cv_dist,
  measure = cv_measure,
  stringsAsFactors = FALSE
)
write.csv(tuning_tbl, file.path(out_dir, "splsda_tuning_results.csv"), row.names = FALSE)

# Export detailed CV tuning grid for Table S
cv_detail <- try({
  # Extract error rates from tune object
  err <- tune$error.rate
  if (!is.null(err)) {
    err_df <- as.data.frame(err)
    err_df$keepX <- rownames(err_df)
    err_df$cv_folds <- cv_folds
    err_df$cv_repeats <- cv_repeats
    err_df$cv_dist <- cv_dist
    err_df$cv_measure <- cv_measure
    err_df
  } else {
    NULL
  }
}, silent = TRUE)
if (!inherits(cv_detail, "try-error") && !is.null(cv_detail)) {
  write.csv(cv_detail, file.path(out_dir, "TableS_CV_tuning.csv"), row.names = FALSE)
}

## -------------------- outputs: loadings table (Table S) --------------------
loadings_tbl <- make_loadings_long(
  spls_model,
  met_features = met_features,
  tax_features = tax_features,
  met_meta = met$meta
)
write.csv(loadings_tbl, file.path(out_dir, "FigS9_splsda_loadings_nonzero.csv"), row.names = FALSE)

## -------------------- outputs: selected features --------------------
selected_tbl <- select_top_by_loading(
  spls_model,
  met_features = met_features,
  tax_features = tax_features,
  met_meta = met$meta,
  n_met = n_select_metabolites,
  n_tax = n_select_taxa
)
write.csv(selected_tbl, file.path(out_dir, "FigS9_selected_features.csv"), row.names = FALSE)

sel_met <- selected_tbl$Feature[selected_tbl$Type == "Metabolite"]
sel_tax <- selected_tbl$Feature[selected_tbl$Type == "Taxon"]

## -------------------- outputs: chord links table --------------------
met_sel <- X_met[common_samples, intersect(sel_met, colnames(X_met)), drop = FALSE]
tax_sel <- X_tax[common_samples, intersect(sel_tax, colnames(X_tax)), drop = FALSE]
if (ncol(met_sel) < 2) stop("Too few selected metabolites for chord after intersection.", call. = FALSE)
if (ncol(tax_sel) < 2) stop("Too few selected taxa for chord after intersection.", call. = FALSE)

links_sig <- compute_links(met_sel, tax_sel, rho_cutoff = rho_threshold, q_cutoff = q_threshold)
write.csv(links_sig, file.path(out_dir, "FigS9_chord_links.csv"), row.names = FALSE)

## -------------------- Chord label mapping (short codes) --------------------
# Create short codes for readability: M1...Mn for metabolites, T1...Tm for taxa
# Use actual column names from met_sel/tax_sel (which are used in links_sig)
met_names_chord <- colnames(met_sel)
tax_names_chord <- colnames(tax_sel)
met_codes <- setNames(paste0("M", seq_along(met_names_chord)), met_names_chord)
tax_codes <- setNames(paste0("T", seq_along(tax_names_chord)), tax_names_chord)
all_codes <- c(met_codes, tax_codes)

# Create label mapping table for Table S
label_map <- data.frame(
  ShortCode = c(met_codes, tax_codes),
  FullName = c(met_names_chord, tax_names_chord),
  Type = c(rep("Metabolite", length(met_names_chord)), rep("Taxon", length(tax_names_chord))),
  stringsAsFactors = FALSE
)
# Add metabolite class info
label_map$SuperClass <- NA_character_
met_in_meta <- met_names_chord[met_names_chord %in% rownames(met$meta)]
if (length(met_in_meta) > 0) {
  label_map$SuperClass[label_map$FullName %in% met_in_meta] <- met$meta[met_in_meta, "Super_Class"]
}
write.csv(label_map, file.path(out_dir, "TableS_chord_labels.csv"), row.names = FALSE)

# Create mapped links for chord (use short codes)
links_mapped <- links_sig
links_mapped$Metabolite <- met_codes[links_sig$Metabolite]
links_mapped$Taxon <- tax_codes[links_sig$Taxon]

# Debug: Check if any mapping failed (NA values)
if (any(is.na(links_mapped$Metabolite))) {
  message("WARNING: Some metabolite names not found in mapping:")
  message("  Not found: ", paste(links_sig$Metabolite[is.na(links_mapped$Metabolite)], collapse = ", "))
  message("  Available: ", paste(names(met_codes), collapse = ", "))
}
if (any(is.na(links_mapped$Taxon))) {
  message("WARNING: Some taxon names not found in mapping:")
  message("  Not found: ", paste(links_sig$Taxon[is.na(links_mapped$Taxon)], collapse = ", "))
  message("  Available: ", paste(names(tax_codes), collapse = ", "))
}
message("Chord sectors (first 5): ", paste(head(unique(c(links_mapped$Metabolite, links_mapped$Taxon)), 5), collapse = ", "))

## -------------------- main figure: A (scores) + B (loadings) + C (chord) --------------------
# sector colors: metabolites by Super_Class, taxa grey
# Note: Use met_names_chord/tax_names_chord for consistency with short codes
met_super <- rep("Unknown", length(met_names_chord))
names(met_super) <- met_names_chord
met_in_meta <- met_names_chord[met_names_chord %in% rownames(met$meta)]
if (length(met_in_meta) > 0) {
  met_super[met_in_meta] <- met$meta[met_in_meta, "Super_Class"]
  met_super[is.na(met_super) | !nzchar(met_super)] <- "Unknown"
}
classes <- unique(met_super)
class_cols <- setNames(grDevices::hcl.colors(length(classes), palette = "Dark 3"), classes)
# Use short codes for colors (met_codes values are M1, M2, etc.)
met_cols <- setNames(class_cols[met_super], met_codes[met_names_chord])
tax_cols <- setNames(rep("grey80", length(tax_names_chord)), tax_codes[tax_names_chord])
grid_col_all <- c(met_cols, tax_cols)

# log2FC for metabolites: mean(25-100%) vs 0% (Control)
ctrl_samples <- grep("^YS0_", common_samples, value = TRUE)
proc_samples <- setdiff(common_samples, ctrl_samples)
eps <- 1e-5
# Use met_sel which contains the actual data columns
mean_proc <- colMeans(X_met_raw[proc_samples, met_names_chord, drop = FALSE])
mean_ctrl <- colMeans(X_met_raw[ctrl_samples, met_names_chord, drop = FALSE])
log2fc <- log2((mean_proc + eps) / (mean_ctrl + eps))

# Use short codes for sector_order and fc_vec
sector_order_all <- c(met_codes[met_names_chord], tax_codes[tax_names_chord])
fc_vec_all <- rep(NA_real_, length(sector_order_all))
names(fc_vec_all) <- sector_order_all
fc_vec_all[met_codes[met_names_chord]] <- log2fc

# Filter to only sectors actually present in links_mapped (for chordDiagram)
chord_sectors <- unique(c(links_mapped$Metabolite, links_mapped$Taxon))
sector_order <- sector_order_all[sector_order_all %in% chord_sectors]
grid_col <- grid_col_all[sector_order]
fc_vec <- fc_vec_all[sector_order]


## -------------------- Plot A: sPLS-DA Scores --------------------
plot_A <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mar = c(4.5, 4.5, 3.5, 1))
  plotIndiv(
    spls_model,
    comp = if (best_ncomp >= 2) 1:2 else 1,
    group = Y,
    ind.names = FALSE,
    ellipse = if (best_ncomp >= 2) TRUE else FALSE,
    legend = TRUE,
    title = sprintf("sPLS-DA scores (ncomp=%d; keepX=%s)", best_ncomp, paste(best_keepX, collapse = ", "))
  )
}

pdf(file.path(out_dir, "Fig9A_scores.pdf"), width = 7, height = 6, useDingbats = FALSE)
plot_A()
dev.off()

tiff(file.path(out_dir, "Fig9A_scores.tiff"), width = 7, height = 6, units = "in", res = 600, compression = "lzw")
plot_A()
dev.off()

## -------------------- Plot B: Loadings (Comp1 + Comp2) --------------------
plot_B <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mfrow = c(1, 2))
  
  # Comp 1 loadings
  par(mar = c(5, 18, 4, 2))
  plotLoadings(
    spls_model,
    comp = 1,
    method = "mean",
    contrib = "max",
    ndisplay = min(loadings_top_n, best_keepX[1]),
    size.name = loadings_label_cex,
    name.var = substring(colnames(spls_model$X), 1, 30),
    legend = FALSE,  # Remove misleading Outcome legend
    title = sprintf("Loadings Comp 1 (keepX=%d; Top %d shown)", best_keepX[1], min(loadings_top_n, best_keepX[1]))
  )
  
  # Comp 2 loadings
  par(mar = c(5, 18, 4, 2))
  if (best_ncomp >= 2) {
    plotLoadings(
      spls_model,
      comp = 2,
      method = "mean",
      contrib = "max",
      ndisplay = min(loadings_top_n, best_keepX[2]),
      size.name = loadings_label_cex,
      name.var = substring(colnames(spls_model$X), 1, 30),
      legend = FALSE,  # Remove misleading Outcome legend
      title = sprintf("Loadings Comp 2 (keepX=%d; Top %d shown)", best_keepX[2], min(loadings_top_n, best_keepX[2]))
    )
  } else {
    plot.new()
    text(0.5, 0.5, "Comp 2 not used (choice.ncomp = 1)")
  }
}

pdf(file.path(out_dir, "Fig9B_loadings.pdf"), width = 14, height = 6, useDingbats = FALSE)
plot_B()
dev.off()

tiff(file.path(out_dir, "Fig9B_loadings.tiff"), width = 14, height = 6, units = "in", res = 600, compression = "lzw")
plot_B()
dev.off()

## -------------------- Plot C: Chord Diagram --------------------
plot_C <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mar = c(1, 1, 3, 1))
  plot_chord(
    links = links_mapped,  # Use short codes
    sector_order = sector_order,
    grid_col = grid_col,
    fc_vec = fc_vec,
    title_text = sprintf("Chord (Spearman; |rho|>%.1f, q<%.2f; see Table S for labels)", rho_threshold, q_threshold),
    n_met = length(sel_met),
    class_colors = class_cols
  )
}

pdf(file.path(out_dir, "Fig9C_chord.pdf"), width = 10, height = 8, useDingbats = FALSE)
plot_C()
dev.off()

tiff(file.path(out_dir, "Fig9C_chord.tiff"), width = 10, height = 8, units = "in", res = 600, compression = "lzw")
plot_C()
dev.off()

## -------------------- Main Fig.9 (A + B1/B2 + C) --------------------
plot_main <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  layout(
    matrix(
      c(
        1, 1,
        2, 3,
        4, 4
      ),
      nrow = 3,
      byrow = TRUE
    ),
    heights = c(1, 1, 1.6)
  )

  # A: scores
  par(mar = c(4.5, 4.5, 3.5, 1))
  plotIndiv(
    spls_model,
    comp = if (best_ncomp >= 2) 1:2 else 1,
    group = Y,
    ind.names = FALSE,
    ellipse = if (best_ncomp >= 2) TRUE else FALSE,
    legend = TRUE,
    title = sprintf("A  sPLS-DA scores (ncomp=%d; keepX=%s)", best_ncomp, paste(best_keepX, collapse = ", "))
  )

  # B1: loadings (Comp1)
  par(mar = c(5, 18, 4, 2))
  plotLoadings(
    spls_model,
    comp = 1,
    method = "mean",
    contrib = "max",
    ndisplay = min(loadings_top_n, best_keepX[1]),
    size.name = loadings_label_cex,
    name.var = substring(colnames(spls_model$X), 1, 30),
    title = "B1  Loadings (Comp 1)"
  )

  # B2: loadings (Comp2; if exists)
  par(mar = c(5, 18, 4, 2))
  if (best_ncomp >= 2) {
    plotLoadings(
      spls_model,
      comp = 2,
      method = "mean",
      contrib = "max",
      ndisplay = min(loadings_top_n, best_keepX[2]),
      size.name = loadings_label_cex,
      name.var = substring(colnames(spls_model$X), 1, 30),
      title = "B2  Loadings (Comp 2)"
    )
  } else {
    plot.new()
    text(0.5, 0.5, "Comp 2 not used (choice.ncomp = 1)")
  }

  # C: chord
  par(mar = c(1, 1, 3, 1))
  plot_chord(
    links = links_sig,
    sector_order = sector_order,
    grid_col = grid_col,
    fc_vec = fc_vec,
    title_text = sprintf("C  Chord (Spearman; |rho|>%.1f, q<%.2f)", rho_threshold, q_threshold),
    n_met = length(sel_met),
    class_colors = class_cols
  )
}

pdf(file.path(out_dir, "Fig9_main.pdf"), width = 12, height = 8.5, useDingbats = FALSE)
plot_main()
dev.off()

tiff(file.path(out_dir, "Fig9_main.tiff"), width = 12, height = 8.5, units = "in", res = 600, compression = "lzw")
plot_main()
dev.off()

message(
  "sPLS-DA: ncomp_max=", ncomp_max,
  "; choice.ncomp=", best_ncomp,
  "; keepX=", paste(best_keepX, collapse = ", "),
  "; CV folds=", cv_folds,
  "; repeats=", cv_repeats
)
message("Selected features: metabolites=", length(sel_met), "; taxa=", length(sel_tax))
message("Chord links: n=", nrow(links_sig), " (|rho|>", rho_threshold, ", q<", q_threshold, ")")
message("Done. Outputs written to: ", out_dir)
