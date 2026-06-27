# ──────────────────────────────────────────────────────────────────────────────
# GAMM — Calibrated Height, Crown Area & Crown:Height Ratio
# Eucalyptus species × spacing trial | EucVision, IMPACT OAL, Stellenbosch
#
# Time series: 1 September 2025 onwards (t0 = 2025-09-01)
# Population:  Single-culture plots only | Living trees only
#              Toggle POPULATION_SUBSET for "All" or "Dominant" (top 20% height)
#
# Three response variables:
#   Response 1: Calibrated_Height_m  — Gaussian, identity link (raw scale)
#   Response 2: Crown_Area_m2        — Gaussian, identity link (log scale)
#   Response 3: CA:H ratio           — Gaussian, identity link (log scale)
#
# WHY LOG TRANSFORMATION FOR CROWN AND CA:H:
#   Gamma(log) failed to converge (4+ hrs) because log(Crown) and log(CA:H)
#   are LEFT-skewed (skew ~ -0.6 to -1.5 by spacing). Gaussian on the log
#   scale is a better fit than Gamma on the original scale.
#   Point predictions back-transformed via smearing correction:
#     E[exp(y)] ~ exp(fit + 0.5 * se^2)
#   Confidence intervals built on log scale then exponentiated:
#     [exp(fit - 1.96*se),  exp(fit + 1.96*se)]
#   This guarantees positive CI bounds and correct asymmetry.
#
# Model structure (per response, two models):
#   m_species: s(days, k=12, bs='cr')
#            + s(days, by=Species,   k=10, bs='tp')
#            + Spacing_f                              [fixed parametric]
#            + s(Plot_ID, bs='re') + s(Tree_ID, bs='re')
#
#   m_spacing: s(days, k=12, bs='cr')
#            + s(days, by=Spacing_f, k=10, bs='tp')
#            + Species                               [fixed parametric]
#            + s(Plot_ID, bs='re') + s(Tree_ID, bs='re')
#
# Global smooth bs='cr' (cubic regression spline) prevents boundary overshoot
# from the 59-day gap between Sep 1 and Oct 30 2025.
#
# Fitted with bam(), method='fREML', discrete=TRUE for computational efficiency.
# Baseline: Species = Grandis | Spacing = 1x1m
#
# Spacing codes:  1x1m =  1 m2/tree  |  2x2m =  4 m2/tree
#                 3x3m =  9 m2/tree  |  5x5m = 25 m2/tree
#
# Outputs saved to OUTPUT_DIR:
#   stats_table1_smooth_significance.csv
#   stats_table2_marginal_means.csv
#   stats_table3_pairwise_differences.csv
#   height / crown / cah  _curves_ species / spacing .png
#   height / crown / cah  _species / spacing _differences.png
#   combined_3x2_curves.png
#   diag_ height / crown / cah  _ species / spacing .png
#   Height/Crown/CAH_Species/Spacing_gamcheck.txt
# ──────────────────────────────────────────────────────────────────────────────


# ── 0. Packages ───────────────────────────────────────────────────────────────
# install.packages(c("mgcv", "tidyverse", "gratia", "patchwork"))
library(mgcv)
library(tidyverse)
library(gratia)
library(patchwork)
library(sf)
library(tictoc)
library(ggrepel)
library(multcompView)

# Add Calibri font
windowsFonts(Calibri = windowsFont("Calibri"))

# tictoc() function for runtime

# 1. Define the custom formatting function
toc_in_mins <- function(tic, toc, msg = "") {
  # Calculate elapsed minutes and round to 2 decimal places
  elapsed_mins <- round((toc - tic) / 60, 2)
  
  # Format the printed console message
  outmsg <- paste0(msg, ": ", elapsed_mins, " minutes elapsed")
  return(outmsg)
}

tic("Model Starts")


# ──────────────────────────────────────────────────────────────────────────────
# OUTPUT SETTINGS
# ──────────────────────────────────────────────────────────────────────────────

OUTPUT_DIR <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/10. GAMM"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# ANALYSIS SETTINGS
# ──────────────────────────────────────────────────────────────────────────────

BASELINE_SPECIES <- "Grandis"
BASELINE_SPACING <- "1x1m"

INCLUDE_SPECIES <- c(
  "Grandis",
  "Grandis clone",
  "Urophylla",
  "Cloeziana",
  "Cladocalyx"
)

# Toggle between "All" (whole population) and "Dominant" (top 20% by height)
POPULATION_SUBSET <- "All" 

cat("\n")
cat("=====================================================\n")
cat("GAMM ANALYSIS SETTINGS\n")
cat("=====================================================\n")
cat("Population subset: ", POPULATION_SUBSET, "\n")
cat("Species included:\n")
print(INCLUDE_SPECIES)
cat("\nBaseline species:", BASELINE_SPECIES, "\n")
cat("Baseline spacing:", BASELINE_SPACING, "\n")
cat("=====================================================\n\n")

# ── 1. Load & clean data ──────────────────────────────────────────────────────

# Define the dates where Crown Area polygons were borrowed
borrowed_dates <- as.Date(c("2025-11-14", "2026-03-16", "2026-04-08", 
                            "2026-04-13", "2026-04-29"))

df_raw <- read_csv("C:/Users/jakev/Downloads/UAV_Master_Dataset_25-05-2026.csv", 
                   show_col_types = FALSE) %>%
  mutate(Tree = round(as.numeric(Tree), 2))

df_base <- df_raw |>
  mutate(
    Date      = as.Date(Date),
    
    # If the flight is a borrowed date, replace Crown Area with NA
    Crown_Area_m2 = if_else(Date %in% borrowed_dates, NA_real_, Crown_Area_m2),
    
    t0        = as.Date("2025-09-01"),
    days      = as.numeric(Date - t0),
    Species   = factor(Species),
    Spacing_f = factor(Spacing,
                       levels = c(1, 2, 3, 5),
                       labels = c("1x1m", "2x2m", "3x3m", "5x5m")),
    Culture   = factor(Culture),
    Plot_ID   = factor(paste0(Compartment, "_", Plot)),
    Tree_ID   = factor(Tree_ID)
  ) |>
  filter(Death_Date == "Alive") |>
  filter(Date >= as.Date("2025-09-01")) |>
  filter(Culture == "Single")

# ── Dominant Tree Filter Switch ──
if (POPULATION_SUBSET == "Dominant") {
  df_base <- df_base |>
    group_by(Species, Spacing_f, Date) |>
    mutate(dom_threshold = quantile(Calibrated_Height_m, probs = 0.80, na.rm = TRUE)) |>
    filter(Calibrated_Height_m >= dom_threshold) |>
    ungroup() |>
    select(-dom_threshold)
}

# ── Height — Gaussian, raw scale ─────────────────────────────────────────────
df_h <- df_base |>
  filter(!is.na(Calibrated_Height_m)) |>
  rename(Height = Calibrated_Height_m)

# ── Crown area — Gaussian on LOG scale; back-transform predictions to m² ─────
# log(Crown) is left-skewed (-0.60 overall) → Gaussian on log scale fits well
# Gamma(log) was wrong: it assumes positive skewness on response scale
df_c <- df_base |>
  filter(!is.na(Crown_Area_m2), Crown_Area_m2 > 0) |>
  mutate(Crown = log(Crown_Area_m2))   # model on log scale

# ── CA:H ratio — Gaussian on LOG scale; back-transform predictions to m² m⁻¹ ─
# log(CA:H) is left-skewed (-0.69 overall) → same fix as crown
df_r <- df_base |>
  filter(!is.na(Crown_Area_m2), Crown_Area_m2 > 0,
         !is.na(Calibrated_Height_m), Calibrated_Height_m > 0) |>
  mutate(CAH = log(Crown_Area_m2 / Calibrated_Height_m))  # model on log scale

# Use Grandis as base reference
df_h$Species <- relevel(df_h$Species, ref = BASELINE_SPECIES)
df_c$Species <- relevel(df_c$Species, ref = BASELINE_SPECIES)
df_r$Species <- relevel(df_r$Species, ref = BASELINE_SPECIES)

df_h$Spacing_f <- relevel(df_h$Spacing_f, ref = BASELINE_SPACING)
df_c$Spacing_f <- relevel(df_c$Spacing_f, ref = BASELINE_SPACING)
df_r$Spacing_f <- relevel(df_r$Spacing_f, ref = BASELINE_SPACING)

cat("── Data summary ──────────────────────────────────────────────────────\n")
cat("Height dataset:     ", nrow(df_h), "obs |", n_distinct(df_h$Tree_ID), "trees\n")
cat("Crown area dataset: ", nrow(df_c), "obs |", n_distinct(df_c$Tree_ID), "trees\n")
cat("CA:H ratio dataset: ", nrow(df_r), "obs |", n_distinct(df_r$Tree_ID), "trees\n")
cat("Days range:         ", min(df_h$days), "to", max(df_h$days), "\n\n")


# ── 2. Shared plot settings ───────────────────────────────────────────────────
# Standardize colors across all potential plots (Matplotlib tab10 equivalent)
species_colors <- c(
  "Cladocalyx"    = "#1f77b4",  # Blue
  "Cloeziana"     = "#ff7f0e",  # Orange
  "Urophylla"     = "#9467bd",  # Purple
  "Grandis"       = "#2ca02c",  # Green
  "Grandis clone" = "#d62728",  # Red
  "Mixed"         = "black"     # Black
)

species_display <- c(
  "Cladocalyx"    = "Cladocalyx",
  "Cloeziana"     = "Cloeziana",
  "Grandis"       = "Grandis",
  "Grandis clone" = "Grandis clone",
  "Urophylla"     = "Urophylla",
  "Mixed"         = "Mixed"
)

spacing_colors <- c(
  "1x1m" = "#118AB2",
  "2x2m" = "#EF476F",
  "3x3m" = "#FFD166",
  "5x5m" = "#06D6A0"
)

spacing_display <- c(
  "1x1m" = "1m",
  "2x2m" = "2m",
  "3x3m" = "3m",
  "5x5m" = "5m"
)


# ──────────────────────────────────────────────────────────────────────────────
# ── SHARED FUNCTIONS ──────────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────

# ── Fit a pair of GAMMs (all three responses use gaussian() ) ─────────────────
fit_gamm_pair <- function(df, response_col) {
  
  cat("  Fitting species model...\n")
  m_sp <- bam(
    as.formula(paste0(response_col, " ~
      s(days, k = 12, bs = 'cr') +
      s(days, by = Species,   k = 10, bs = 'tp') +
      Spacing_f +
      s(Plot_ID, bs = 're') +
      s(Tree_ID, bs = 're')")),
    data     = df,
    family   = gaussian(),
    method   = "fREML",
    discrete = TRUE     
  )
  
  cat("  Fitting spacing model...\n")
  m_sc <- bam(
    as.formula(paste0(response_col, " ~
      s(days, k = 12, bs = 'cr') +
      s(days, by = Spacing_f, k = 10, bs = 'tp') +
      Species +
      s(Plot_ID, bs = 're') +
      s(Tree_ID, bs = 're')")),
    data     = df,
    family   = gaussian(),
    method   = "fREML",
    discrete = TRUE         
  )
  
  list(species = m_sp, spacing = m_sc)
}


# ── Predict, with optional delta-method back-transformation from log scale ────
# backtransform = TRUE  → use for Crown and CA:H (fitted on log scale)
# backtransform = FALSE → use for Height (fitted on raw scale)
predict_traj <- function(model, newdata, backtransform = FALSE) {
  
  preds <- predict(
    model,
    newdata = newdata,
    se.fit = TRUE,
    type = "response",
    exclude = c("s(Plot_ID)", "s(Tree_ID)")
  )
  
  fit_raw <- as.numeric(preds$fit)
  se_raw  <- as.numeric(preds$se.fit)
  
  if (backtransform) {
    fit_out <- exp(fit_raw + 0.5 * se_raw^2)  
    lwr_out <- exp(fit_raw - 1.96 * se_raw)   
    upr_out <- exp(fit_raw + 1.96 * se_raw)
  } else {
    fit_out <- fit_raw
    lwr_out <- fit_raw - 1.96 * se_raw
    upr_out <- fit_raw + 1.96 * se_raw
  }
  
  newdata |>
    mutate(fit = fit_out, lwr = lwr_out, upr = upr_out, se = se_raw)
}

# ── Marginal Means Table Function ─────────────────────────────────────────────
marginal_means <- function(model, group_var, group_levels,
                           fixed_var, fixed_level,
                           df_ref, key_days, key_labels, response_label,
                           backtransform = FALSE) {
  
  pred_base <- tibble(
    !!group_var := factor(group_levels, levels = levels(df_ref[[group_var]])),
    !!fixed_var := factor(fixed_level,  levels = levels(df_ref[[fixed_var]])),
    Culture = factor("Single", levels = levels(df_ref$Culture)),
    Plot_ID = levels(df_ref$Plot_ID)[1],
    Tree_ID = levels(df_ref$Tree_ID)[1]
  )
  
  map2_dfr(key_days, key_labels, function(d, lbl) {
    nd      <- pred_base |> mutate(days = d)
    
    preds   <- predict(model, newdata = nd, se.fit = TRUE,
                       type = "response",
                       exclude = c("s(Plot_ID)", "s(Tree_ID)"))
    
    fit_raw <- as.numeric(preds$fit)
    se_raw  <- as.numeric(preds$se.fit)
    
    if (backtransform) {
      fit_out <- exp(fit_raw + 0.5 * se_raw^2)
      se_out  <- exp(fit_raw) * se_raw
    } else {
      fit_out <- fit_raw
      se_out  <- se_raw
    }
    
    nd |>
      mutate(Response  = response_label,
             Timepoint = lbl,
             Mean      = round(fit_out, 3),
             SE        = round(se_out, 3),
             CI_lower  = round(Mean - 1.96 * SE, 3),
             CI_upper  = round(Mean + 1.96 * SE, 3)) |>
      select(Response, Timepoint, !!group_var, Mean, SE, CI_lower, CI_upper)
  })
}


# ── Rigorous lpmatrix Time-Series Differences ─────────────────────────────────
pairwise_diffs_rigorous <- function(model, df_ref, pred_grid, group_var, is_log_scale = FALSE) {
  
  grp_levels <- if (is.factor(pred_grid[[group_var]])) levels(pred_grid[[group_var]]) else unique(pred_grid[[group_var]])
  pairs <- combn(grp_levels, 2, simplify = FALSE)
  
  # Identify the fixed variable to hold constant
  fixed_var <- ifelse(group_var == "Species", "Spacing_f", "Species")
  fixed_level <- levels(df_ref[[fixed_var]])[1]
  
  map_dfr(pairs, function(pair) {
    lv1 <- pair[1]; lv2 <- pair[2]
    
    # 1. Extract the pre-calculated means from your prediction grid
    g1 <- pred_grid |> filter(.data[[group_var]] == lv1) |> arrange(days)
    g2 <- pred_grid |> filter(.data[[group_var]] == lv2) |> arrange(days)
    
    # 2. Build explicit newdata frames for lpmatrix extraction
    nd1 <- g1 |> mutate(!!fixed_var := factor(fixed_level, levels = levels(df_ref[[fixed_var]])),
                        Plot_ID = levels(df_ref$Plot_ID)[1], Tree_ID = levels(df_ref$Tree_ID)[1])
    
    nd2 <- g2 |> mutate(!!fixed_var := factor(fixed_level, levels = levels(df_ref[[fixed_var]])),
                        Plot_ID = levels(df_ref$Plot_ID)[1], Tree_ID = levels(df_ref$Tree_ID)[1])
    
    # 3. Extract lpmatrix (excluding plot/tree random effects)
    X1 <- predict(model, newdata = nd1, type = "lpmatrix", exclude = c("s(Plot_ID)", "s(Tree_ID)"))
    X2 <- predict(model, newdata = nd2, type = "lpmatrix", exclude = c("s(Plot_ID)", "s(Tree_ID)"))
    
    # 4. Covariance-aware standard error calculation
    Xdiff <- X1 - X2
    V     <- vcov(model, unconditional = TRUE) # unconditional accounts for smoothing parameter uncertainty
    
    # Calculate difference and SE strictly on the model's link scale
    fit_link <- as.numeric(Xdiff %*% coef(model))
    se_link  <- sqrt(rowSums((Xdiff %*% V) * Xdiff))
    
    # 95% CI on the link scale (using strict 1.96 multiplier)
    lwr_link <- fit_link - 1.96 * se_link
    upr_link <- fit_link + 1.96 * se_link
    
    if (is_log_scale) {
      
      # Absolute difference on the back-transformed scale (for the trendline)
      abs_diff <- g1$fit - g2$fit
      
      # --- Ribbon and rug: both derived from the lpmatrix CI on the log scale ---
      # fit_link  = log(group1) - log(group2)  =  log(group1 / group2)
      # lwr_link  = fit_link - 1.96 * se_link  (on log-ratio scale)
      # upr_link  = fit_link + 1.96 * se_link  (on log-ratio scale)
      #
      # To get the CI bounds on the ABSOLUTE difference scale we use:
      #   diff_upr = g2$fit * (exp(upr_link) - 1)   [upper bound on abs difference]
      #   diff_lwr = g2$fit * (exp(lwr_link) - 1)   [lower bound on abs difference]
      #
      # Derivation: if D = mu1 - mu2 and R = mu1/mu2 = exp(fit_link), then
      #   D = mu2 * (R - 1), so CI on D = mu2 * (exp(lwr/upr_link) - 1)
      
      lwr_abs <- g2$fit * (exp(lwr_link) - 1)
      upr_abs <- g2$fit * (exp(upr_link) - 1)
      
      # SE for the stats table (back-transformed via delta method)
      se_abs <- se_link * g2$fit   # delta method: se(D) ≈ se(log R) * mu2
      
      # Rug and ribbon are now fully consistent — both from lpmatrix CI
      sig <- (lwr_abs > 0) | (upr_abs < 0)
      
      tibble(comparison = paste0(lv1, " - ", lv2), days = g1$days,
             diff = abs_diff, se_diff = se_abs, lwr = lwr_abs, upr = upr_abs, sig = sig)
      
    } else {
      # For Height (raw scale): link scale is the absolute scale
      tibble(comparison = paste0(lv1, " - ", lv2), days = g1$days,
             diff = fit_link, se_diff = se_link, lwr = lwr_link, upr = upr_link, 
             sig = (lwr_link > 0) | (upr_link < 0))
    }
  })
}

# ── Custom Thesis Theme ───────────────────────────────────────────────────────
theme_thesis <- function() {
  theme_classic(base_size = 9, base_family = "Calibri") +
    theme(
      # Text and Titles
      plot.title       = element_text(size = 10, face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8),
      
      # Spines and Ticks
      axis.line        = element_line(colour = "black", linewidth = 1.2),
      axis.ticks       = element_line(colour = "black", linewidth = 1),
      axis.ticks.length = unit(4, "pt"),
      
      # Gridlines (y-axis only, dashed, grey, transparent)
      panel.grid.major.y = element_line(colour = alpha("#b0b0b0", 0.25), linewidth = 0.5, linetype = "solid"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      
      # Legend Settings (Moved to top left, border RESTORED)
      legend.position      = "top", 
      legend.justification = "left",
      legend.background    = element_rect(fill = "white", colour = "lightgray", linewidth = 0.5), # <-- Border is back!
      legend.title         = element_text(size = 10, face = "bold"),
      legend.text          = element_text(size = 8),
      legend.key.size      = unit(0.4, "cm"),
      legend.margin        = margin(t = 2, r = 5, b = 2, l = 5, unit = "pt"), # Adds a bit of breathing room inside the box
      
      # Trim outer margins to save vertical space
      plot.margin = margin(t = 2, r = 5, b = 2, l = 2, unit = "pt")
    )
}

# ── Shared curve plot builder ─────────────────────────────────────────────────
curve_plot <- function(curve_df, colour_var, colour_vals, colour_labels = NULL,
                       y_label, legend_title = NULL,
                       title = NULL, subtitle = NULL) {
  ggplot(curve_df, aes(x = days,
                       colour = .data[[colour_var]],
                       fill   = .data[[colour_var]])) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
    geom_line(aes(y = fit), linewidth = 0.9) +
    scale_colour_manual(values = colour_vals, labels = colour_labels, drop = FALSE) +
    scale_fill_manual(values = colour_vals, labels = colour_labels, drop = FALSE) +
    guides(colour = guide_legend(nrow = 2), fill = "none") +
    scale_x_continuous(breaks = seq(0, 270, by = 60), limits = c(0, 290)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "Days from 1 September 2025", y = y_label,
         colour = legend_title, fill = legend_title,
         title = title, subtitle = subtitle) +
    theme_thesis() +
    theme(
      plot.title    = element_text(size = 9, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.5, colour = "grey40", margin = margin(b = 3))
    )
}

# ── Plot pairwise differences ─────────────────────────────────────────────────
plot_diffs <- function(diff_df, y_label, fill_col = "steelblue", ncol = 3,
                       title = NULL, subtitle = NULL) {
  
  global_min_y <- min(diff_df$lwr, na.rm = TRUE)
  global_max_y <- max(diff_df$upr, na.rm = TRUE)
  rug_y_pos    <- global_min_y - ((global_max_y - global_min_y) * 0.05)
  
  ggplot(diff_df, aes(x = days, y = diff)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = fill_col) +
    geom_line(colour = fill_col, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_point(data = diff_df[diff_df$sig, ],
               aes(x = days, y = rug_y_pos),
               colour = "#ff3333", alpha = 0.8,
               shape = 124, size = 2) +
    facet_wrap(~ comparison, ncol = ncol) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    labs(x = "Days from 1 September 2025", y = y_label,
         title = title, subtitle = subtitle) +
    theme_thesis() +
    theme(
      plot.title    = element_text(size = 9, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.5, colour = "grey40", margin = margin(b = 3))
    )
}
  
  # ── Build main panel (9 comparisons, shared Y) ────────────────────────────
  p_main <- make_base(df_main, free_y = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = if (!is.null(df_outlier)) "" else "Days from 1 September 2025") +
    theme(
      axis.title.x = if (!is.null(df_outlier)) element_blank() else element_text(),
      axis.text.x  = if (!is.null(df_outlier)) element_blank() else element_text(),
      axis.ticks.x = if (!is.null(df_outlier)) element_blank() else element_line()
    )
  
  # ── If no outlier, return main plot directly ───────────────────────────────
  if (is.null(df_outlier)) return(p_main)
  
  # ── Build outlier panel (1 comparison, its own Y) ─────────────────────────
  p_outlier <- make_base(df_outlier, free_y = TRUE) +
    labs(title = NULL, subtitle = NULL,
         x = "Days from 1 September 2025", y = y_label)
  
  # ── Stack with patchwork: 9 panels tall, 1 panel below ────────────────────
  p_main / p_outlier +
    plot_layout(heights = c(3, 1))
}

# ── Statistics helpers ────────────────────────────────────────────────────────
smooth_sig_table <- function(model, response_label) {
  s      <- summary(model)
  sp_tbl <- as.data.frame(s$s.table)
  sp_tbl <- sp_tbl[
    grepl("^s\\(days\\):", rownames(sp_tbl)),
  ]
  sp_tbl$Term     <- rownames(sp_tbl)
  sp_tbl$Response <- response_label
  sp_tbl$Sig      <- ifelse(sp_tbl[["p-value"]] < 0.001, "***",
                            ifelse(sp_tbl[["p-value"]] < 0.01,  "**",
                                   ifelse(sp_tbl[["p-value"]] < 0.05,  "*",
                                          ifelse(sp_tbl[["p-value"]] < 0.1,   ".",  "ns"))))
  sp_tbl |>
    select(Response, Term, edf = edf, F = F, p = `p-value`, Sig) |>
    mutate(edf = round(edf, 2), F = round(F, 3), p = sprintf("%.2e", p))
}

pairwise_at_day <- function(diff_df, target_day, response_label, factor_label) {
  diff_df |>
    filter(abs(days - target_day) == min(abs(days - target_day))) |>
    slice(1, .by = comparison) |>
    mutate(
      z_stat     = diff / se_diff,
      p_val      = 2 * pnorm(-abs(z_stat)),          # two-tailed
      Response   = response_label,
      Factor     = factor_label,
      Timepoint  = paste0("Day ", target_day, " from 1 Sep 2025"),
      Difference = round(diff, 3),
      SE         = round(se_diff, 3),
      CI_lower   = round(lwr, 3),
      CI_upper   = round(upr, 3),
      Sig        = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        p_val < 0.10  ~ ".",
        TRUE          ~ "ns"
      )
    ) |>
    select(Response, Factor, Comparison = comparison,
           Difference, SE, CI_lower, CI_upper, Sig)
}


# ── Prediction grid helpers ───────────────────────────────────────────────────
make_full_grid <- function(days_seq, df_ref) {
  expand_grid(
    days = days_seq,
    Species   = levels(df_ref$Species),
    Spacing_f = levels(df_ref$Spacing_f)
  ) |>
    mutate(
      Species   = factor(Species, levels = levels(df_ref$Species)),
      Spacing_f = factor(Spacing_f, levels = levels(df_ref$Spacing_f)),
      Culture   = factor("Single", levels = levels(df_ref$Culture)),
      Plot_ID   = levels(df_ref$Plot_ID)[1],
      Tree_ID   = levels(df_ref$Tree_ID)[1]
    )
}

make_species_grid <- function(days_seq, df_ref) {
  expand_grid(
    days = days_seq,
    Species = levels(df_ref$Species)
  ) |>
    mutate(
      Species   = factor(Species, levels = levels(df_ref$Species)),
      Spacing_f = levels(df_ref$Spacing_f)[1],
      Culture   = factor("Single", levels = levels(df_ref$Culture)),
      Plot_ID   = levels(df_ref$Plot_ID)[1],
      Tree_ID   = levels(df_ref$Tree_ID)[1]
    )
}

make_spacing_grid <- function(days_seq, df_ref) {
  expand_grid(
    days = days_seq,
    Spacing_f = levels(df_ref$Spacing_f)
  ) |>
    mutate(
      Spacing_f = factor(Spacing_f, levels = levels(df_ref$Spacing_f)),
      Species   = levels(df_ref$Species)[1],
      Culture   = factor("Single", levels = levels(df_ref$Culture)),
      Plot_ID   = levels(df_ref$Plot_ID)[1],
      Tree_ID   = levels(df_ref$Tree_ID)[1]
    )
}

# ──────────────────────────────────────────────────────────────────────────────
# ── FIT MODELS ────────────────────────────────────────────────────────────────
# All three responses: Gaussian family (log-scale for Crown and CA:H)
# If models already in memory, skip to PREDICTION GRIDS
# Expected runtime: ~5-10 min per model pair (~30 min total)
# ──────────────────────────────────────────────────────────────────────────────

cat("══ RESPONSE 1: Calibrated Height (Gaussian, raw scale) ════════════════\n")
models_h <- fit_gamm_pair(df_h, "Height")
cat("\n── Height species model summary ─────────────────────────────────────\n")
print(summary(models_h$species))
cat("\n── Height spacing model summary ─────────────────────────────────────\n")
print(summary(models_h$spacing))

cat("\n══ RESPONSE 2: Crown Area (Gaussian, log scale) ════════════════════════\n")
cat("   Fitted on log(Crown_Area_m2); predictions back-transformed to m²\n\n")
models_c <- fit_gamm_pair(df_c, "Crown")
cat("\n── Crown species model summary ──────────────────────────────────────\n")
print(summary(models_c$species))
cat("\n── Crown spacing model summary ──────────────────────────────────────\n")
print(summary(models_c$spacing))

cat("\n══ RESPONSE 3: CA:H Ratio (Gaussian, log scale) ════════════════════════\n")
cat("   CA:H = Crown_Area_m2 / Height_m  (m2 m-1)\n")
cat("   Fitted on log(CA:H); predictions back-transformed to m2 m-1\n\n")
models_r <- fit_gamm_pair(df_r, "CAH")
cat("\n── CA:H species model summary ───────────────────────────────────────\n")
print(summary(models_r$species))
cat("\n── CA:H spacing model summary ───────────────────────────────────────\n")
print(summary(models_r$spacing))

# ──────────────────────────────────────────────────────────────────────────────
# MODEL DIAGNOSTICS
# ──────────────────────────────────────────────────────────────────────────────

cat("\n")
cat("=====================================================\n")
cat("MODEL DIAGNOSTICS\n")
cat("=====================================================\n")

# ──────────────────────────────────────────────────────────────────────────────
# Helper: run gam.check() and save output
# ──────────────────────────────────────────────────────────────────────────────

save_gam_check <- function(model, model_name, output_dir) {
  
  cat("\n-----------------------------------------------------\n")
  cat(model_name, "\n")
  cat("-----------------------------------------------------\n")
  
  txt_file <- file.path(
    output_dir,
    paste0(
      gsub("[^A-Za-z0-9]", "_", model_name),
      "_gamcheck.txt"
    )
  )
  
  capture.output(
    gam.check(model),
    file = txt_file
  )
  
  cat("Saved:", basename(txt_file), "\n")
  
  invisible(txt_file)
}

# ──────────────────────────────────────────────────────────────────────────────
# Run and save GAM checks
# ──────────────────────────────────────────────────────────────────────────────

save_gam_check(
  models_h$species,
  "Height Species",
  OUTPUT_DIR
)

save_gam_check(
  models_h$spacing,
  "Height Spacing",
  OUTPUT_DIR
)

save_gam_check(
  models_c$species,
  "Crown Area Species",
  OUTPUT_DIR
)

save_gam_check(
  models_c$spacing,
  "Crown Area Spacing",
  OUTPUT_DIR
)

save_gam_check(
  models_r$species,
  "CAH Ratio Species",
  OUTPUT_DIR
)

save_gam_check(
  models_r$spacing,
  "CAH Ratio Spacing",
  OUTPUT_DIR
)

cat("\n")
cat("All GAM diagnostic reports saved.\n")
cat("=====================================================\n")

# ──────────────────────────────────────────────────────────────────────────────
# ── PREDICTION GRIDS ──────────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────

days_h <- seq(min(df_h$days), max(df_h$days), length.out = 300)
days_c <- seq(min(df_c$days), max(df_c$days), length.out = 300)
days_r <- seq(min(df_r$days), max(df_r$days), length.out = 300)

cat("Predicting height trajectories...\n")
sp_diffs_h <- pairwise_diffs_rigorous(models_h$species, df_h, 
                                      predict_traj(models_h$species, make_species_grid(days_h, df_h)),
                                      "Species", is_log_scale = FALSE)
sc_diffs_h <- pairwise_diffs_rigorous(models_h$spacing, df_h,
                                      predict_traj(models_h$spacing, make_spacing_grid(days_h, df_h)),
                                      "Spacing_f", is_log_scale = FALSE)

cat("Predicting crown area trajectories (back-transforming to m2)...\n")
sp_diffs_c <- pairwise_diffs_rigorous(models_c$species, df_c,
                                      predict_traj(models_c$species, make_species_grid(days_c, df_c), backtransform = TRUE),
                                      "Species", is_log_scale = TRUE)
sc_diffs_c <- pairwise_diffs_rigorous(models_c$spacing, df_c,
                                      predict_traj(models_c$spacing, make_spacing_grid(days_c, df_c), backtransform = TRUE),
                                      "Spacing_f", is_log_scale = TRUE)

cat("Predicting CA:H ratio trajectories (back-transforming to m2 m-1)...\n")
sp_diffs_r <- pairwise_diffs_rigorous(models_r$species, df_r,
                                      predict_traj(models_r$species, make_species_grid(days_r, df_r), backtransform = TRUE),
                                      "Species", is_log_scale = TRUE)
sc_diffs_r <- pairwise_diffs_rigorous(models_r$spacing, df_r,
                                      predict_traj(models_r$spacing, make_spacing_grid(days_r, df_r), backtransform = TRUE),
                                      "Spacing_f", is_log_scale = TRUE)

cat("\nRange checks (all should be non-zero):\n")
cat("Height   species:", round(range(sp_diffs_h$diff), 3), "\n")
cat("Height   spacing:", round(range(sc_diffs_h$diff), 3), "\n")
cat("Crown    species:", round(range(sp_diffs_c$diff), 3), "\n")
cat("Crown    spacing:", round(range(sc_diffs_c$diff), 3), "\n")
cat("CA:H     species:", round(range(sp_diffs_r$diff), 3), "\n")
cat("CA:H     spacing:", round(range(sc_diffs_r$diff), 3), "\n")

# Growth curve grids
curve_h_sp <- predict_traj(models_h$species, make_species_grid(days_h, df_h))
curve_h_sc <- predict_traj(models_h$spacing, make_spacing_grid(days_h, df_h))

curve_c_sp <- predict_traj(models_c$species, make_species_grid(days_c, df_c), backtransform = TRUE)
curve_c_sc <- predict_traj(models_c$spacing, make_spacing_grid(days_c, df_c), backtransform = TRUE)

curve_r_sp <- predict_traj(models_r$species, make_species_grid(days_r, df_r), backtransform = TRUE)
curve_r_sc <- predict_traj(models_r$spacing, make_spacing_grid(days_r, df_r), backtransform = TRUE)


# ──────────────────────────────────────────────────────────────────────────────
# ── PLOTS ─────────────────────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# ── FIGURE CAPTION REFERENCE (Former Titles & Subtitles) ──────────────────────
# Use these to draft your Word document figure captions.
#
# ── DIFFERENCE PLOTS ──────────────────────────────────────────────────────────
# Species Pairwise Differences
#   Title:    Species pairwise [height / crown area / crown:height ratio] differences
#   Subtitle: Shaded = 95% CI  |  Red rug = significant period  |  Averaged over spacings
#
# Spacing Pairwise Differences
#   Title:    Spacing pairwise [height / crown area / crown:height ratio] differences
#   Subtitle: Shaded = 95% CI  |  Red rug = significant period  |  Averaged over species
#
# ── SINGLE CURVE PLOTS ────────────────────────────────────────────────────────
# Species Trajectories
#   Title:    GAMM-fitted [height / crown area / crown:height ratio] growth trajectories by species
#   Subtitle: Population mean, spacing controlled  |  Shaded = 95% CI
#
# Spacing Trajectories
#   Title:    GAMM-fitted [height / crown area / crown:height ratio] growth trajectories by spacing
#   Subtitle: Population mean, species controlled  |  Shaded = 95% CI
#
# ── COMBINED 3x2 GRID ─────────────────────────────────────────────────────────
#   Title:    GAMM-Fitted Growth Trajectories
#   Subtitle: Left column: Species-controlled  |  Right column: Spacing-controlled
# ──────────────────────────────────────────────────────────────────────────────

cat("\nGenerating plots...\n")

# ── HEIGHT plots ──────────────────────────────────────────────────────────────
p_h_sp_diff <- plot_diffs(sp_diffs_h,
                          y_label  = "Difference in height (m)",
                          fill_col = "steelblue", ncol = 3,
                          title    = "Species pairwise height differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over spacings")

p_h_sc_diff <- plot_diffs(sc_diffs_h,
                          y_label  = "Difference in height (m)",
                          fill_col = "darkorange", ncol = 3,
                          title    = "Spacing pairwise height differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over species")
ggsave(file.path(OUTPUT_DIR, "height_spacing_differences.png"), p_h_sc_diff,
       width = 6.30, height = 3.0, units = "in", dpi = 300)

p_h_sp_curves <- curve_plot(curve_h_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title  = "Species",
                            y_label       = "Calibrated height (m)",
                            title         = "GAMM-fitted height growth trajectories by species",
                            subtitle      = "Population mean, spacing controlled  |  Shaded = 95% CI")

p_h_sc_curves <- curve_plot(curve_h_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title  = "Spacing",
                            y_label       = "Calibrated height (m)",
                            title         = "GAMM-fitted height growth trajectories by spacing",
                            subtitle      = "Population mean, species controlled  |  Shaded = 95% CI")

# ── CROWN AREA plots ──────────────────────────────────────────────────────────
p_c_sp_diff <- plot_diffs(sp_diffs_c,
                          y_label  = "Difference in crown area (m\u00b2)",
                          fill_col = "#1b9e77", ncol = 3,
                          title    = "Species pairwise crown area differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over spacings")
ggsave(file.path(OUTPUT_DIR, "crown_species_differences.png"), p_c_sp_diff,
       width = 6.30, height = 5.5, units = "in", dpi = 300)

p_c_sc_diff <- plot_diffs(sc_diffs_c,
                          y_label  = "Difference in crown area (m\u00b2)",
                          fill_col = "#d95f02", ncol = 3,
                          title    = "Spacing pairwise crown area differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over species")
ggsave(file.path(OUTPUT_DIR, "crown_spacing_differences.png"), p_c_sc_diff,
       width = 6.30, height = 3.0, units = "in", dpi = 300)

p_c_sp_curves <- curve_plot(curve_c_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title  = "Species",
                            y_label       = "Crown area (m\u00b2)",
                            title         = "GAMM-fitted crown area growth trajectories by species",
                            subtitle      = "Population mean, spacing controlled  |  Shaded = 95% CI")

p_c_sc_curves <- curve_plot(curve_c_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title  = "Spacing",
                            y_label       = "Crown area (m\u00b2)",
                            title         = "GAMM-fitted crown area growth trajectories by spacing",
                            subtitle      = "Population mean, species controlled  |  Shaded = 95% CI")

# ── CA:H RATIO plots ──────────────────────────────────────────────────────────
p_r_sp_diff <- plot_diffs(sp_diffs_r,
                          y_label  = "Difference in CA:H ratio (m\u00b2 m\u207b\u00b9)",
                          fill_col = "#6a3d9a", ncol = 3,
                          title    = "Species pairwise CA:H ratio differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over spacings")
ggsave(file.path(OUTPUT_DIR, "cah_species_differences.png"), p_r_sp_diff,
       width = 6.30, height = 5.5, units = "in", dpi = 300)

p_r_sc_diff <- plot_diffs(sc_diffs_r,
                          y_label  = "Difference in CA:H ratio (m\u00b2 m\u207b\u00b9)",
                          fill_col = "#e31a1c", ncol = 3,
                          title    = "Spacing pairwise CA:H ratio differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over species")
ggsave(file.path(OUTPUT_DIR, "cah_spacing_differences.png"), p_r_sc_diff,
       width = 6.30, height = 3.0, units = "in", dpi = 300)

p_r_sp_curves <- curve_plot(curve_r_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title  = "Species",
                            y_label       = "CA:H ratio (m\u00b2 m\u207b\u00b9)",
                            title         = "GAMM-fitted CA:H ratio trajectories by species",
                            subtitle      = "Population mean, spacing controlled  |  Shaded = 95% CI")

p_r_sc_curves <- curve_plot(curve_r_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title  = "Spacing",
                            y_label       = "CA:H ratio (m\u00b2 m\u207b\u00b9)",
                            title         = "GAMM-fitted CA:H ratio trajectories by spacing",
                            subtitle      = "Population mean, species controlled  |  Shaded = 95% CI")

# ── AUTOMATED LABEL FUNCTION ──────────────────────────────────────────────────
attach_labels_auto <- function(curve_df, diff_df, group_var, target_day, y_positions) {
  
  # 1. Isolate the pairwise differences for the final day
  day_data <- diff_df |>
    filter(abs(days - target_day) == min(abs(days - target_day))) |>
    slice(1, .by = comparison)
  
  # 2. Create a logical vector of significance
  is_diff <- (day_data$lwr > 0) | (day_data$upr < 0)
  names(is_diff) <- gsub(" - ", "-", day_data$comparison)
  
  # 3. Generate the automated letters
  cld <- multcompView::multcompLetters(is_diff)$Letters
  letters_df <- tibble(
    !!group_var := names(cld),
    letter = as.character(cld)
  )
  
  # 4. Attach the letters to the endpoints of your curves
  curve_df |>
    filter(days == max(days)) |>
    arrange(desc(fit)) |> 
    left_join(letters_df, by = group_var) |>
    mutate(y_lab = y_positions) |>
    select(days, fit, y_lab, all_of(group_var), letter)
}

# ── SIGNIFICANCE LABELS (DAY 266) ─────────────────────────────────────────────

# 1. HEIGHT LABELS
lbl_h_sp <- attach_labels_auto(
  curve_df    = curve_h_sp,
  diff_df     = sp_diffs_h,
  group_var   = "Species",
  target_day  = 266,
  y_positions = c(4.1, 3.80, 3.5, 3.2, 2.85)
)

p_h_sp_curves <- p_h_sp_curves +
  geom_label(
    data = lbl_h_sp,
    aes(
      x = max(days) + 12,
      y = y_lab,
      label = letter,
      fill = Species
    ),
    color = "white",
    fontface = "bold",
    size = 3,
    label.r = unit(0, "lines"),
    show.legend = FALSE
  )

lbl_h_sc <- attach_labels_auto(
  curve_df    = curve_h_sc,
  diff_df     = sc_diffs_h,
  group_var   = "Spacing_f",
  target_day  = 266,
  y_positions = c(4.5, 4.2, 3.9, 3.6)
)

p_h_sc_curves <- p_h_sc_curves +
  geom_label(
    data = lbl_h_sc,
    aes(
      x = max(days) + 12,
      y = y_lab,
      label = letter,
      fill = Spacing_f
    ),
    color = "white",
    fontface = "bold",
    size = 3,
    label.r = unit(0, "lines"),
    show.legend = FALSE
  )

# 2. CROWN AREA LABELS
lbl_c_sp <- attach_labels_auto(
  curve_df    = curve_c_sp,
  diff_df     = sp_diffs_c,
  group_var   = "Species",
  target_day  = 266,
  y_positions = c(1.0, 0.92, 0.84, 0.76, 0.68)
)

p_c_sp_curves <- p_c_sp_curves +
  geom_label(
    data = lbl_c_sp,
    aes(
      x = max(days) + 12,
      y = y_lab,
      label = letter,
      fill = Species
    ),
    color = "white",
    fontface = "bold",
    size = 3,
    label.r = unit(0, "lines"),
    show.legend = FALSE
  )

lbl_c_sc <- attach_labels_auto(
  curve_df    = curve_c_sc,
  diff_df     = sc_diffs_c,
  group_var   = "Spacing_f",
  target_day  = 266,
  y_positions = c(2.35, 2.1, 1.8, 1.5)
)

p_c_sc_curves <- p_c_sc_curves +
  geom_label(
    data = lbl_c_sc,
    aes(
      x = max(days) + 12,
      y = y_lab,
      label = letter,
      fill = Spacing_f
    ),
    color = "white",
    fontface = "bold",
    size = 3,
    label.r = unit(0, "lines"),
    show.legend = FALSE
  )

# 3. CA:H RATIO LABELS
lbl_r_sp <- attach_labels_auto(
  curve_df    = curve_r_sp,
  diff_df     = sp_diffs_r,
  group_var   = "Species",
  target_day  = 266,
  y_positions = c(0.32, 0.295, 0.27, 0.245, 0.21)
)

p_r_sp_curves <- p_r_sp_curves +
  geom_label(
    data = lbl_r_sp,
    aes(
      x = max(days) + 12,
      y = y_lab,
      label = letter,
      fill = Species
    ),
    color = "white",
    fontface = "bold",
    size = 3,
    label.r = unit(0, "lines"),
    show.legend = FALSE
  )

lbl_r_sc <- attach_labels_auto(
  curve_df    = curve_r_sc,
  diff_df     = sc_diffs_r,
  group_var   = "Spacing_f",
  target_day  = 266,
  y_positions = c(0.58, 0.52, 0.46, 0.40)
)

p_r_sc_curves <- p_r_sc_curves +
  geom_label(
    data = lbl_r_sc,
    aes(
      x = max(days) + 12,
      y = y_lab,
      label = letter,
      fill = Spacing_f
    ),
    color = "white",
    fontface = "bold",
    size = 3,
    label.r = unit(0, "lines"),
    show.legend = FALSE
  )

# ── SAVE LABELLED CURVE PLOTS
ggsave(file.path(OUTPUT_DIR, "height_curves_species.png"), p_h_sp_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)
ggsave(file.path(OUTPUT_DIR, "height_curves_spacing.png"), p_h_sc_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)
ggsave(file.path(OUTPUT_DIR, "crown_curves_species.png"),  p_c_sp_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)
ggsave(file.path(OUTPUT_DIR, "crown_curves_spacing.png"),  p_c_sc_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)
ggsave(file.path(OUTPUT_DIR, "cah_curves_species.png"),    p_r_sp_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)
ggsave(file.path(OUTPUT_DIR, "cah_curves_spacing.png"),    p_r_sc_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)

# ── Combined 3x2 Figure (Max Width) ───────────────────────────────────────────

cat("\nAssembling 3x2 combined grid (independent Y-axes)...\n")

# 1. Helper function to strip axes and titles on inner plots
clean_panel <- function(p, keep_legend = FALSE, keep_x = FALSE, keep_y = TRUE,
                        keep_subtitle = FALSE) {
  p <- p + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
  
  if (!keep_legend)  p <- p + theme(legend.position  = "none")
  if (!keep_x)       p <- p + theme(axis.title.x     = element_blank(),
                                    axis.text.x      = element_blank(),
                                    axis.ticks.x     = element_blank())
  if (!keep_y)       p <- p + theme(axis.title.y     = element_blank())
  
  # Always strip the individual plot title; optionally strip the subtitle too
  p <- p + theme(plot.title = element_blank())
  if (!keep_subtitle) p <- p + theme(plot.subtitle = element_blank())
  
  return(p)
}

# 2. Apply cleaning to all 6 panels
#    Only c_sp_3x2 (top-left) keeps its subtitle — it reads "Left column: Species-controlled"
#    Re-write that subtitle to serve as a column-header hint for the reader
c_sp_3x2 <- clean_panel(p_c_sp_curves, keep_legend = TRUE,  keep_x = FALSE,
                        keep_y = TRUE,  keep_subtitle = FALSE)

c_sc_3x2 <- clean_panel(p_c_sc_curves, keep_legend = TRUE,  keep_x = FALSE,
                        keep_y = FALSE, keep_subtitle = FALSE)

h_sp_3x2 <- clean_panel(p_h_sp_curves, keep_legend = FALSE, keep_x = FALSE,
                        keep_y = TRUE,  keep_subtitle = FALSE)
h_sc_3x2 <- clean_panel(p_h_sc_curves, keep_legend = FALSE, keep_x = FALSE,
                        keep_y = FALSE, keep_subtitle = FALSE)

r_sp_3x2 <- clean_panel(p_r_sp_curves, keep_legend = FALSE, keep_x = TRUE,
                        keep_y = TRUE,  keep_subtitle = FALSE)
r_sc_3x2 <- clean_panel(p_r_sc_curves, keep_legend = FALSE, keep_x = TRUE,
                        keep_y = FALSE, keep_subtitle = FALSE)

# 3. Assemble and annotate with a single main title
p_combined_3x2 <- (c_sp_3x2 | c_sc_3x2) /
  (h_sp_3x2   | h_sc_3x2)   /
  (r_sp_3x2   | r_sc_3x2)   +
  plot_annotation(
    title    = "GAMM-fitted tree growth trajectories by species and spacing",
    subtitle = "Left column: spacing-controlled | Right column: species-controlled | Shaded = 95% CI",
    theme = theme(
      plot.title = element_text(
        size   = 10,
        face   = "bold",
        hjust  = 0.5,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        size   = 7.5,
        colour = "grey40",
        hjust  = 0.5,
        margin = margin(b = 4)
      )
    )
  )

# Height bumped by 0.2 in to give the main title breathing room
ggsave(file.path(OUTPUT_DIR, "combined_3x2_curves.png"), p_combined_3x2,
       width = 6.30, height = 6.7, units = "in", dpi = 300)

# ──────────────────────────────────────────────────────────────────────────────
# ── STATISTICS SUMMARY TABLES ─────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────

cat("\nGenerating statistics tables...\n")

key_days   <- c(59, 134, 203, 266)
key_labels <- c("2 months", "4.5 months", "6.5 months", "9 months")
final_day  <- 266

# TABLE 1: Smooth term significance — all six models
tbl1 <- bind_rows(
  smooth_sig_table(models_h$species, "Height (species model)"),
  smooth_sig_table(models_h$spacing, "Height (spacing model)"),
  smooth_sig_table(models_c$species, "Crown Area (species model)"),
  smooth_sig_table(models_c$spacing, "Crown Area (spacing model)"),
  smooth_sig_table(models_r$species, "CA:H Ratio (species model)"),
  smooth_sig_table(models_r$spacing, "CA:H Ratio (spacing model)")
)

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TABLE 1: GAMM Smooth Term Significance                             ║\n")
cat("║  EDF > 1 = non-linear  |  p < 0.05 = significant smooth            ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
print(as_tibble(tbl1), n = Inf)
write_csv(
  tbl1,
  file.path(
    OUTPUT_DIR,
    "stats_table1_smooth_significance.csv"
  )
)

# TABLE 2: Marginal means at key timepoints
# Height — raw scale (no back-transform)
tbl2_h_sp <- marginal_means(models_h$species, "Species",   levels(df_h$Species),
                            "Spacing_f", levels(df_h$Spacing_f)[1],
                            df_h, key_days, key_labels, "Height (m)")
tbl2_h_sc <- marginal_means(models_h$spacing, "Spacing_f", levels(df_h$Spacing_f),
                            "Species",   levels(df_h$Species)[1],
                            df_h, key_days, key_labels, "Height (m)")

# Crown — back-transform from log scale to m²
tbl2_c_sp <- marginal_means(models_c$species, "Species",   levels(df_c$Species),
                            "Spacing_f", levels(df_c$Spacing_f)[1],
                            df_c, key_days, key_labels, "Crown Area (m2)",
                            backtransform = TRUE)
tbl2_c_sc <- marginal_means(models_c$spacing, "Spacing_f", levels(df_c$Spacing_f),
                            "Species",   levels(df_c$Species)[1],
                            df_c, key_days, key_labels, "Crown Area (m2)",
                            backtransform = TRUE)

# CA:H — back-transform from log scale to m² m⁻¹
tbl2_r_sp <- marginal_means(models_r$species, "Species",   levels(df_r$Species),
                            "Spacing_f", levels(df_r$Spacing_f)[1],
                            df_r, key_days, key_labels, "CA:H Ratio (m2 m-1)",
                            backtransform = TRUE)
tbl2_r_sc <- marginal_means(models_r$spacing, "Spacing_f", levels(df_r$Spacing_f),
                            "Species",   levels(df_r$Species)[1],
                            df_r, key_days, key_labels, "CA:H Ratio (m2 m-1)",
                            backtransform = TRUE)

tbl2_list <- list(
  list(tbl2_h_sp, "2A: Heights by Species"),
  list(tbl2_h_sc, "2B: Heights by Spacing"),
  list(tbl2_c_sp, "2C: Crown Area by Species"),
  list(tbl2_c_sc, "2D: Crown Area by Spacing"),
  list(tbl2_r_sp, "2E: CA:H Ratio by Species"),
  list(tbl2_r_sc, "2F: CA:H Ratio by Spacing")
)

for (item in tbl2_list) {
  cat(paste0("\n╔══ TABLE ", item[[2]], " ══╗\n"))
  print(item[[1]] |>
          pivot_wider(names_from = Timepoint,
                      values_from = c(Mean, SE),
                      names_glue  = "{Timepoint} {.value}"), n = Inf)
}

write_csv(bind_rows(tbl2_h_sp, tbl2_h_sc, tbl2_c_sp,
                    tbl2_c_sc, tbl2_r_sp, tbl2_r_sc),
          file.path(OUTPUT_DIR,"stats_table2_marginal_means.csv"))

# TABLE 3: Pairwise differences at day 266
tbl3 <- bind_rows(
  pairwise_at_day(sp_diffs_h, final_day, "Height (m)",          "Species"),
  pairwise_at_day(sc_diffs_h, final_day, "Height (m)",          "Spacing"),
  pairwise_at_day(sp_diffs_c, final_day, "Crown Area (m2)",     "Species"),
  pairwise_at_day(sc_diffs_c, final_day, "Crown Area (m2)",     "Spacing"),
  pairwise_at_day(sp_diffs_r, final_day, "CA:H Ratio (m2 m-1)", "Species"),
  pairwise_at_day(sc_diffs_r, final_day, "CA:H Ratio (m2 m-1)", "Spacing")
)

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TABLE 3: All Pairwise Differences at Day 266 (25 May 2026)        ║\n")
cat("║  Difference = Group1 minus Group2 | Sig = CI excludes zero         ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
print(tbl3, n = Inf)
write_csv(tbl3, file.path(OUTPUT_DIR,"stats_table3_pairwise_differences.csv"))

cat("\n── Significant pairs at Day 266 ─────────────────────────────────────\n")
sig_only <- tbl3 |> filter(Sig %in% c("*", "**", "***"))
for (resp in c("Height (m)", "Crown Area (m2)", "CA:H Ratio (m2 m-1)")) {
  for (fac in c("Species", "Spacing")) {
    sub <- sig_only |> filter(Response == resp, Factor == fac)
    if (nrow(sub) > 0) {
      cat(paste0("\n", resp, " — ", fac, ":\n"))
      print(sub |> select(Comparison, Difference, CI_lower, CI_upper))
    }
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# ── DIAGNOSTICS ───────────────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────

cat("\n── k adequacy checks ────────────────────────────────────────────────\n")
model_list <- list(
  "Height species"  = models_h$species, "Height spacing"  = models_h$spacing,
  "Crown species"   = models_c$species, "Crown spacing"   = models_c$spacing,
  "CA:H species"    = models_r$species, "CA:H spacing"    = models_r$spacing
)
for (nm in names(model_list)) {
  cat(nm, "model:\n"); print(k.check(model_list[[nm]]))
}

diag_files <- list(
  list(models_h$species, file.path(OUTPUT_DIR, "diag_height_species.png"), "steelblue"),
  list(models_h$spacing, file.path(OUTPUT_DIR, "diag_height_spacing.png"), "darkorange"),
  list(models_c$species, file.path(OUTPUT_DIR, "diag_crown_species.png"), "#1b9e77"),
  list(models_c$spacing, file.path(OUTPUT_DIR, "diag_crown_spacing.png"), "#d95f02"),
  list(models_r$species, file.path(OUTPUT_DIR, "diag_cah_species.png"), "#6a3d9a"),
  list(models_r$spacing, file.path(OUTPUT_DIR, "diag_cah_spacing.png"), "#e31a1c")
)

for (item in diag_files) {
  p_diag <- appraise(item[[1]], 
                     point_col = item[[3]], 
                     point_alpha = 0.3, 
                     line_col = "black") & 
    theme_thesis()
  
  # Perfect 4.5 x 4.5 squares to be centered in Word
  ggsave(item[[2]], p_diag,
         width = 4.5, height = 4.5, units = "in", dpi = 300)
}


cat("\n── Saved outputs ─────────────────────────────────────────────────────\n")
cat("  STATISTICS TABLES\n")
cat("    stats_table1_smooth_significance.csv\n")
cat("    stats_table2_marginal_means.csv\n")
cat("    stats_table3_pairwise_differences.csv\n")
cat("  HEIGHT\n")
cat("    height_curves_species/spacing.png\n")
cat("    height_species/spacing_differences.png\n")
cat("  CROWN AREA  (back-transformed from log scale to m2)\n")
cat("    crown_curves_species/spacing.png\n")
cat("    crown_species/spacing_differences.png\n")
cat("  CA:H RATIO  (back-transformed from log scale to m2 m-1)\n")
cat("    cah_curves_species/spacing.png\n")
cat("    cah_species/spacing_differences.png\n")
cat("  COMBINED GRIDS\n")
cat("    combined_3x2_curves.png\n")
cat("  DIAGNOSTICS\n")
cat("    diag_height/crown/cah _species/_spacing .png\n")
cat("\n── Done ──────────────────────────────────────────────────────────────\n")

# End counter
toc(func.toc = toc_in_mins)
