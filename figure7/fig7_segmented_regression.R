# Figure 7 | Representative exudates: segmented regression (Muggeo) + breakpoint 95% CI
# R: 4.4.3 (as requested in ToDO_fig7.md)
#
# Input
# - data/完整数据-分泌物.csv
#
# Output (written into figure7/)
# - Fig7_main.pdf
# - Fig7_main.tiff
# - TableS_breakpoints.csv
# - FigS7_breakpoint_diagnostics.pdf (optional; set make_diagnostics <- TRUE)
#
# Notes
# - log2FC is computed relative to Control mean (RS = 0).
# - Segmented regression uses `segmented` (Muggeo class) with 1 breakpoint and 95% CI from `confint()`.

## ----------------------------- user config -----------------------------
input_csv <- file.path("data", "完整数据-分泌物.csv")
out_dir <- "figure7"

# Representative metabolites for Figure 7 (recommended: paste the exact "Name" from the CSV).
# If NULL, the script auto-selects `top_n` metabolites by max(abs(mean log2FC)) across RS levels.
representative_compounds <- NULL
top_n <- 6

# segmented regression settings
psi_guesses <- c(25, 50, 75) # initial breakpoint guesses (RS)
conf_level <- 0.95

# Optional outputs
make_diagnostics <- FALSE

# Figure layout
facet_ncol <- 3
fig_width_in <- 10
fig_height_in <- 7
fig_dpi <- 600
## ---------------------------------------------------------------------

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is required. Install it via install.packages('%s').", pkg, pkg),
      call. = FALSE
    )
  }
}

suggest_compounds <- function(query, available, n = 10) {
  hits <- available[grepl(query, available, ignore.case = TRUE)]
  if (length(hits) > 0) {
    return(utils::head(unique(hits), n))
  }
  hits2 <- agrep(query, available, ignore.case = TRUE, value = TRUE, max.distance = 0.2)
  utils::head(unique(hits2), n)
}

extract_est_ci <- function(seg_fit, conf_level = 0.95) {
  psi_mat <- seg_fit$psi
  if (is.null(psi_mat)) {
    psi_mat <- summary(seg_fit)$psi
  }
  if (is.null(psi_mat) || nrow(psi_mat) < 1) {
    return(list(est = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

est_col <- grep("Est", colnames(psi_mat), value = TRUE)
  if (length(est_col) == 0) {
    est <- as.numeric(psi_mat[1, 1])
  } else {
    est <- as.numeric(psi_mat[1, est_col[1]])
  }

  ci_low <- NA_real_
  ci_high <- NA_real_
  ci_mat <- tryCatch(stats::confint(seg_fit, level = conf_level), error = function(e) NULL)
  if (!is.null(ci_mat)) {
    if (is.vector(ci_mat) && length(ci_mat) >= 2) {
      ci_low <- as.numeric(ci_mat[1])
      ci_high <- as.numeric(ci_mat[2])
    } else if (is.matrix(ci_mat) && nrow(ci_mat) >= 1 && ncol(ci_mat) >= 2) {
      ci_low <- as.numeric(ci_mat[1, 1])
      ci_high <- as.numeric(ci_mat[1, 2])
    }
  }

  list(est = est, ci_low = ci_low, ci_high = ci_high)
}

fit_one_segmented <- function(df_one, psi_guesses = c(25, 50, 75)) {
  df_one <- df_one[is.finite(df_one$RS) & is.finite(df_one$log2FC), , drop = FALSE]
  if (nrow(df_one) < 6) {
    return(list(lm_fit = NULL, seg_fit = NULL))
  }

  lm_fit <- stats::lm(log2FC ~ RS, data = df_one)
  seg_fit <- NULL

  for (psi in psi_guesses) {
    seg_fit <- tryCatch(
      segmented::segmented(lm_fit, seg.Z = ~RS, psi = list(RS = psi)),
      error = function(e) NULL
    )
    if (!is.null(seg_fit)) break
  }

  list(lm_fit = lm_fit, seg_fit = seg_fit)
}

predict_with_ci <- function(model, rs_grid, conf_level = 0.95) {
  grid <- data.frame(RS = rs_grid)
  pred <- tryCatch(stats::predict(model, newdata = grid, se.fit = TRUE), error = function(e) NULL)

  if (!is.null(pred) && is.list(pred) && !is.null(pred$fit) && !is.null(pred$se.fit)) {
    fit <- as.numeric(pred$fit)
    se <- as.numeric(pred$se.fit)
    df_res <- tryCatch(stats::df.residual(model), error = function(e) NA_real_)
    alpha <- 1 - conf_level
    tcrit <- stats::qt(1 - alpha / 2, df = df_res)
    lwr <- fit - tcrit * se
    upr <- fit + tcrit * se
  } else {
    fit <- as.numeric(stats::predict(model, newdata = grid))
    lwr <- rep(NA_real_, length(fit))
    upr <- rep(NA_real_, length(fit))
  }

  data.frame(RS = rs_grid, fit = fit, lwr = lwr, upr = upr)
}

## ----------------------------- main script -----------------------------
require_pkg("ggplot2")
require_pkg("segmented")

if (!file.exists(input_csv)) {
  stop(sprintf("Input file not found: %s (run from repo root?)", input_csv), call. = FALSE)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_fig_pdf <- file.path(out_dir, "Fig7_main.pdf")
out_fig_tiff <- file.path(out_dir, "Fig7_main.tiff")
out_break_table <- file.path(out_dir, "TableS_breakpoints.csv")
out_diag_pdf <- file.path(out_dir, "FigS7_breakpoint_diagnostics.pdf")

exu <- utils::read.csv(
  input_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

if (!("Name" %in% colnames(exu))) {
  stop("Column 'Name' is required in 完整数据-分泌物.csv.", call. = FALSE)
}

sample_cols <- grep("^YS[0-9]+_F-[0-9]+$", colnames(exu), value = TRUE)
if (length(sample_cols) == 0) {
  stop("No sample columns matched pattern ^YS[0-9]+_F-[0-9]+$ in 完整数据-分泌物.csv.", call. = FALSE)
}

sample_meta <- data.frame(
  Sample = sample_cols,
  RS = as.numeric(sub("^YS([0-9]+)_F-[0-9]+$", "\\1", sample_cols)),
  Rep = as.integer(sub("^YS[0-9]+_F-([0-9]+)$", "\\1", sample_cols)),
  stringsAsFactors = FALSE
)
sample_meta <- sample_meta[order(sample_meta$RS, sample_meta$Rep), , drop = FALSE]
sample_cols <- sample_meta$Sample

if (!any(sample_meta$RS == 0)) {
  stop("No Control group detected (RS == 0). Expected columns like YS0_F-1.", call. = FALSE)
}

int_mat <- as.matrix(exu[, sample_cols, drop = FALSE])
mode(int_mat) <- "numeric"

control_cols <- sample_meta$Sample[sample_meta$RS == 0]
control_mean <- rowMeans(int_mat[, control_cols, drop = FALSE], na.rm = TRUE)

n_compounds <- nrow(exu)
n_samples <- length(sample_cols)

long_df <- data.frame(
  Compound = rep(exu$Name, each = n_samples),
  Sample = rep(sample_cols, times = n_compounds),
  RS = rep(sample_meta$RS, times = n_compounds),
  Rep = rep(sample_meta$Rep, times = n_compounds),
  Intensity = as.vector(t(int_mat)),
  ControlMean = rep(control_mean, each = n_samples),
  stringsAsFactors = FALSE
)
long_df$log2FC <- log2(long_df$Intensity / long_df$ControlMean)

available_compounds <- unique(long_df$Compound)

# Global citrate sanity check (ToDO_fig7.md: "where is citrate?")
if (!any(tolower(available_compounds) == "citrate")) {
  citrate_like <- available_compounds[grepl("citrate|citric", available_compounds, ignore.case = TRUE)]
  if (length(citrate_like) > 0) {
    message("NOTE: No exact metabolite named 'Citrate' was found in the CSV.")
    message("      Found citrate-like names: ", paste(utils::head(unique(citrate_like), 10), collapse = ", "))
  }
}

if (is.null(representative_compounds)) {
  mean_l2fc <- stats::aggregate(log2FC ~ Compound + RS, data = long_df, FUN = function(x) mean(x, na.rm = TRUE))
  mean_l2fc <- mean_l2fc[mean_l2fc$RS != 0, , drop = FALSE]
  max_abs <- tapply(abs(mean_l2fc$log2FC), mean_l2fc$Compound, max, na.rm = TRUE)
  max_abs <- sort(max_abs, decreasing = TRUE)
  representative_compounds <- names(utils::head(max_abs, top_n))
  message("Auto-selected representative compounds (top ", top_n, " by max|mean log2FC|):")
  message(paste0(" - ", representative_compounds, collapse = "\n"))
}

missing <- setdiff(representative_compounds, available_compounds)
if (length(missing) > 0) {
  if (any(tolower(missing) == "citrate")) {
    stop(
      "Requested compound 'Citrate' is NOT present in 完整数据-分泌物.csv.\n",
      "Please confirm the correct metabolite name/annotation (e.g., Triethylcitrate?) and fix the list.\n",
      call. = FALSE
    )
  }

  message("Requested compounds not found (exact match): ", paste(missing, collapse = ", "))
  for (q in missing) {
    sug <- suggest_compounds(q, available_compounds, n = 10)
    if (length(sug) > 0) {
      message("  Suggestions for '", q, "': ", paste(sug, collapse = ", "))
    }
  }
  stop("Please update `representative_compounds` with exact names from the CSV.", call. = FALSE)
}

plot_df <- long_df[long_df$Compound %in% representative_compounds, , drop = FALSE]
plot_df$Compound <- factor(plot_df$Compound, levels = representative_compounds)

bp_rows <- list()
pred_rows <- list()
fit_store <- list()

rs_grid <- seq(min(plot_df$RS, na.rm = TRUE), max(plot_df$RS, na.rm = TRUE), length.out = 200)

for (comp in representative_compounds) {
  df_one <- plot_df[plot_df$Compound == comp, c("RS", "log2FC"), drop = FALSE]
  fits <- fit_one_segmented(df_one, psi_guesses = psi_guesses)
  fit_store[[comp]] <- fits

  if (!is.null(fits$seg_fit)) {
    est_ci <- extract_est_ci(fits$seg_fit, conf_level = conf_level)
    pred <- predict_with_ci(fits$seg_fit, rs_grid = rs_grid, conf_level = conf_level)
    pred$Compound <- comp

    bp_rows[[comp]] <- data.frame(
      Compound = comp,
      breakpoint = est_ci$est,
      ci_lower = est_ci$ci_low,
      ci_upper = est_ci$ci_high,
      stringsAsFactors = FALSE
    )
    pred_rows[[comp]] <- pred
  } else if (!is.null(fits$lm_fit)) {
    pred <- predict_with_ci(fits$lm_fit, rs_grid = rs_grid, conf_level = conf_level)
    pred$Compound <- comp
    pred_rows[[comp]] <- pred

    bp_rows[[comp]] <- data.frame(
      Compound = comp,
      breakpoint = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    bp_rows[[comp]] <- data.frame(
      Compound = comp,
      breakpoint = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

bp_table <- do.call(rbind, bp_rows)
pred_df <- do.call(rbind, pred_rows)

utils::write.csv(bp_table, out_break_table, row.names = FALSE, fileEncoding = "UTF-8")

bp_vline <- bp_table[is.finite(bp_table$breakpoint), , drop = FALSE]
bp_band <- bp_table[is.finite(bp_table$ci_lower) & is.finite(bp_table$ci_upper), , drop = FALSE]

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = RS, y = log2FC)) +
  ggplot2::geom_point(
    ggplot2::aes(),
    position = ggplot2::position_jitter(width = 1.2, height = 0),
    alpha = 0.55,
    size = 1.6,
    color = "gray30"
  ) +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 2.1, color = "black") +
  ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 0.6, color = "black") +
  ggplot2::geom_ribbon(
    data = pred_df,
    ggplot2::aes(x = RS, ymin = lwr, ymax = upr, group = Compound),
    inherit.aes = FALSE,
    alpha = 0.15,
    fill = "#1f77b4"
  ) +
  ggplot2::geom_line(
    data = pred_df,
    ggplot2::aes(x = RS, y = fit, group = Compound),
    inherit.aes = FALSE,
    linewidth = 1.0,
    color = "#1f77b4"
  ) +
  ggplot2::geom_rect(
    data = bp_band,
    ggplot2::aes(xmin = ci_lower, xmax = ci_upper, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.06,
    fill = "red"
  ) +
  ggplot2::geom_vline(
    data = bp_vline,
    ggplot2::aes(xintercept = breakpoint),
    linetype = "dotted",
    color = "red",
    linewidth = 0.7
  ) +
  ggplot2::scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  ggplot2::facet_wrap(~Compound, scales = "free_y", ncol = facet_ncol) +
  ggplot2::labs(
    title = "Figure 7 | Segmented regression (1 breakpoint) of representative root exudates",
    x = "RS (%)",
    y = "log2FC vs Control mean (RS=0)"
  ) +
  ggplot2::theme_bw(base_size = 11, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(out_fig_pdf, p, width = fig_width_in, height = fig_height_in, units = "in")

save_tiff <- function(filename, plot_obj) {
  # Prefer ragg if available (better text rendering); fallback to grDevices::tiff.
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename,
      plot_obj,
      width = fig_width_in,
      height = fig_height_in,
      units = "in",
      device = ragg::agg_tiff,
      dpi = fig_dpi
    )
  } else {
    ggplot2::ggsave(
      filename,
      plot_obj,
      width = fig_width_in,
      height = fig_height_in,
      units = "in",
      device = "tiff",
      dpi = fig_dpi,
      compression = "lzw"
    )
  }
}
save_tiff(out_fig_tiff, p)

if (isTRUE(make_diagnostics)) {
  grDevices::pdf(out_diag_pdf, width = 8.5, height = 6.5)
  for (comp in representative_compounds) {
    fits <- fit_store[[comp]]
    df_one <- plot_df[plot_df$Compound == comp, c("RS", "log2FC"), drop = FALSE]
    df_one <- df_one[is.finite(df_one$RS) & is.finite(df_one$log2FC), , drop = FALSE]

    grDevices::par(mfrow = c(1, 2))
    graphics::plot(
      df_one$RS,
      df_one$log2FC,
      main = comp,
      xlab = "RS (%)",
      ylab = "log2FC vs Control mean",
      pch = 16,
      col = grDevices::adjustcolor("gray30", alpha.f = 0.7)
    )

    if (!is.null(fits$seg_fit)) {
      est_ci <- extract_est_ci(fits$seg_fit, conf_level = conf_level)
      pred <- predict_with_ci(fits$seg_fit, rs_grid = rs_grid, conf_level = conf_level)
      graphics::lines(pred$RS, pred$fit, col = "#1f77b4", lwd = 2)
      graphics::abline(v = est_ci$est, col = "red", lty = 2)
      graphics::mtext(sprintf("breakpoint = %.2f", est_ci$est), side = 3, line = 0.2, cex = 0.8)

      graphics::plot(
        stats::fitted(fits$seg_fit),
        stats::resid(fits$seg_fit),
        main = "Residuals vs fitted",
        xlab = "Fitted",
        ylab = "Residuals",
        pch = 16,
        col = grDevices::adjustcolor("gray30", alpha.f = 0.7)
      )
      graphics::abline(h = 0, lty = 2)
    } else if (!is.null(fits$lm_fit)) {
      pred <- predict_with_ci(fits$lm_fit, rs_grid = rs_grid, conf_level = conf_level)
      graphics::lines(pred$RS, pred$fit, col = "#1f77b4", lwd = 2)

      graphics::plot(
        stats::fitted(fits$lm_fit),
        stats::resid(fits$lm_fit),
        main = "Residuals vs fitted (LM fallback)",
        xlab = "Fitted",
        ylab = "Residuals",
        pch = 16,
        col = grDevices::adjustcolor("gray30", alpha.f = 0.7)
      )
      graphics::abline(h = 0, lty = 2)
    } else {
      graphics::mtext("Model fit failed.", side = 3, line = 0.2, cex = 0.9)
    }
  }
  grDevices::dev.off()
}

message("Done.")
message(" - ", out_fig_pdf)
message(" - ", out_fig_tiff)
message(" - ", out_break_table)
if (isTRUE(make_diagnostics)) message(" - ", out_diag_pdf)

