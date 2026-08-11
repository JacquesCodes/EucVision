# ──────────────────────────────────────────────────────────────────────────────
# GAMM — Calibrated Height, Crown Area & Crown:Height Ratio
# ──────────────────────────────────────────────────────────────────────────────
# Author: Jacques Vermeulen
# Project: EucXylo (https://eucxylo.sun.ac.za/)
# ──────────────────────────────────────────────────────────────────────────────
# Eucalyptus species × spacing trial | IMPACT OAL, Stellenbosch
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
# Model structure (per response, ONE combined model):
#   model_x: s(days, k=12, bs='cr')
#          + s(days, by=Species,   k=10, bs='tp')
#          + s(days, by=Spacing_f, k=10, bs='tp')
#          + Species                              [fixed parametric — level shift]
#          + Spacing_f                            [fixed parametric — level shift]
#          + s(Plot_ID, bs='re') + s(Tree_ID, bs='re')
#
# WHY A SINGLE COMBINED MODEL (not separate species/spacing models):
#   A pilot comparison (AIC + curve back-prediction) showed the combined
#   model fits substantially better than either split model for all three
#   responses (AIC differences in the thousands), and that split-model curves
#   were measurably biased relative to the combined model — most severely for
#   the factor NOT being smoothed in that split model (e.g. spacing curves
#   from the species-only Height model were off by up to ~0.41 m). Height is
#   driven more by species identity; Crown Area and CA:H are driven more by
#   spacing — so omitting either factor's smooth term biases that model's
#   under-represented dimension. The combined model avoids this by letting
#   both factors have their own smooth AND parametric terms simultaneously.
#
# WHY THE BY-SMOOTH FACTOR IS ALSO PARAMETRIC:
#   s(days, by=factor) smooths carry a sum-to-zero constraint per level, so
#   they capture SHAPE only, relative to one shared intercept. Including the
#   same factor as a bare parametric term supplies the missing group-level
#   intercept shift, so species/spacing differences in overall magnitude are
#   correctly attributed rather than absorbed into a single global intercept.
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
#   diag_ height / crown / cah _combined.png
#   Height/Crown/CAH_Combined_gamcheck.txt
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
# OUTPUT SETTINGS ####
# ──────────────────────────────────────────────────────────────────────────────

OUTPUT_DIR <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/10. GAMM"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ──────────────────────────────────────────────────────────────────────────────
# ANALYSIS SETTINGS ####
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

# Human-readable label for plot subtitles and filenames
POPULATION_LABEL <- if (POPULATION_SUBSET == "Dominant") {
  "Dominant trees (top 20% by height)"
} else {
  "Full population"
}

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

df_raw <- read_csv("C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/01. Data Analysis/UAV_Master_Dataset_25-05-2026.csv", 
                   show_col_types = FALSE) %>%
  mutate(Tree = round(as.numeric(Tree), 2))

df_base <- df_raw |>
  mutate(
    Date      = as.Date(Date),
    Crown_Area_m2 = if_else(Date %in% borrowed_dates, NA_real_, Crown_Area_m2),
    t0        = as.Date("2025-09-01"),
    days      = as.numeric(Date - t0),
    
    # Define the exact order right here!
    Species   = factor(Species, levels = c(
      "Grandis",       # Baseline (Column 1)
      "Cladocalyx", 
      "Cloeziana", 
      "Urophylla", 
      "Grandis clone"  # Bottom row
    )),
    
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
# Use 1m by 1m as base reference
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
  "1x1m" = "1x1m",
  "2x2m" = "2x2m",
  "3x3m" = "3x3m",
  "5x5m" = "5x5m"
)

# ──────────────────────────────────────────────────────────────────────────────
# SHARED FUNCTIONS ####
# ──────────────────────────────────────────────────────────────────────────────

# ── Fit a single combined GAMM, optionally with group-specific residual variance ──
# var_groups = NULL           → original constant-variance fit (unchanged behaviour)
# var_groups = c("Species", "Spacing_f") → two-stage variance-weighted fit:
#   pass 1 estimates residual variance per group cell, pass 2 refits with
#   prior weights w = 1/sigma^2_group. bam() assumes Var(y) = sigma^2 / w,
#   so this gives each cell its own residual variance rather than pooling one
#   global value. This is the Welch-analogue for a GAMM: unequal-variance
#   groups no longer borrow each other's scatter when SEs are computed.
fit_gamm_combined <- function(df, response_col, var_groups = NULL, n_iter = 2) {
  
  form <- as.formula(paste0(response_col, " ~
      s(days, k = 12, bs = 'cr') +
      s(days, by = Species,   k = 10, bs = 'tp') +
      s(days, by = Spacing_f, k = 10, bs = 'tp') +
      Species +
      Spacing_f +
      s(Plot_ID, bs = 're') +
      s(Tree_ID, bs = 're')"))
  
  cat("  Pass 1: constant-variance fit...\n")
  df$.w <- 1
  m <- bam(form, data = df, family = gaussian(),
           method = "fREML", discrete = TRUE, weights = .w)
  
  if (is.null(var_groups)) {
    attr(m, "var_groups") <- NULL
    return(m)
  }
  
  grp <- interaction(df[var_groups], drop = TRUE)
  
  for (i in seq_len(n_iter)) {
    r <- residuals(m, type = "response")
    
    # Mean squared residual per group cell = that cell's residual variance.
    # Normalised by the geometric mean so weights centre near 1 and the
    # overall scale parameter stays interpretable.
    v <- tapply(r^2, grp, mean)
    v <- v / exp(mean(log(v)))
    
    df$.w <- 1 / as.numeric(v[as.character(grp)])
    
    cat("  Pass ", i + 1, ": variance-weighted refit ",
        "(weight range ", round(min(df$.w), 3), " to ",
        round(max(df$.w), 3), ")...\n", sep = "")
    
    m <- bam(form, data = df, family = gaussian(),
             method = "fREML", discrete = TRUE, weights = .w)
  }
  
  attr(m, "var_groups")  <- var_groups
  attr(m, "group_var")   <- v
  m
}

# ── Duan (1983) nonparametric smearing factors ────────────────────────────────
# For a model fitted on the log scale, E[exp(y)] = exp(mu) * mean(exp(resid)).
# Nonparametric: makes no lognormality assumption, which matters here because
# log(Crown) and log(CA:H) are left-skewed. Computed per variance group so that
# groups with wider residual spread get a correspondingly larger correction.
smearing_factors <- function(model, df) {
  r  <- residuals(model, type = "response")
  vg <- attr(model, "var_groups")
  
  if (is.null(vg)) {
    return(list(groups = NULL, factors = c(.global = mean(exp(r)))))
  }
  grp <- interaction(df[vg], drop = TRUE)
  list(groups = vg, factors = tapply(exp(r), grp, mean))
}

# Look up the right factor for each row of a prediction grid.
# Cells absent from the data (e.g. Grandis clone x 5x5m) fall back to the
# geometric mean of the observed factors.
smear_lookup <- function(smear, newdata) {
  if (is.null(smear))         return(rep(1, nrow(newdata)))
  if (is.null(smear$groups))  return(rep(as.numeric(smear$factors[1]), nrow(newdata)))
  
  g <- interaction(newdata[smear$groups], drop = FALSE)
  k <- as.numeric(smear$factors[as.character(g)])
  k[is.na(k)] <- exp(mean(log(smear$factors)))
  k
}

# ── Predict, with optional Duan smearing back-transformation from log scale ───
# backtransform = TRUE  → Crown and CA:H (fitted on log scale)
# backtransform = FALSE → Height (fitted on raw scale)
# The smearing factor scales the point estimate AND both CI bounds by the same
# constant, so the interval stays consistent with the point estimate and keeps
# its asymmetry on the response scale.
predict_traj <- function(model, newdata, backtransform = FALSE, smear = NULL) {
  
  preds <- predict(
    model,
    newdata = newdata,
    se.fit  = TRUE,
    type    = "link",
    exclude = c("s(Plot_ID)", "s(Tree_ID)")
  )
  
  fit_link <- as.numeric(preds$fit)
  se_link  <- as.numeric(preds$se.fit)
  
  if (backtransform) {
    k       <- smear_lookup(smear, newdata)
    fit_out <- k * exp(fit_link)
    lwr_out <- k * exp(fit_link - 1.96 * se_link)
    upr_out <- k * exp(fit_link + 1.96 * se_link)
  } else {
    fit_out <- fit_link
    lwr_out <- fit_link - 1.96 * se_link
    upr_out <- fit_link + 1.96 * se_link
  }
  
  newdata |>
    mutate(fit = fit_out, lwr = lwr_out, upr = upr_out, se = se_link)
}

# ── Marginal Means Table Function ─────────────────────────────────────────────
marginal_means <- function(model, group_var, group_levels,
                           fixed_var, fixed_level,
                           df_ref, key_days, key_labels, response_label,
                           backtransform = FALSE, smear = NULL) {
  
  pred_base <- tibble(
    !!group_var := factor(group_levels, levels = levels(df_ref[[group_var]])),
    !!fixed_var := factor(fixed_level,  levels = levels(df_ref[[fixed_var]])),
    Culture = factor("Single", levels = levels(df_ref$Culture)),
    Plot_ID = levels(df_ref$Plot_ID)[1],
    Tree_ID = levels(df_ref$Tree_ID)[1]
  )
  
  map2_dfr(key_days, key_labels, function(d, lbl) {
    nd <- pred_base |> mutate(days = d)
    
    preds <- predict(model, newdata = nd, se.fit = TRUE, type = "link",
                     exclude = c("s(Plot_ID)", "s(Tree_ID)"))
    
    fit_link <- as.numeric(preds$fit)
    se_link  <- as.numeric(preds$se.fit)
    
    if (backtransform) {
      k        <- smear_lookup(smear, nd)
      mean_out <- k * exp(fit_link)
      lwr_out  <- k * exp(fit_link - 1.96 * se_link)
      upr_out  <- k * exp(fit_link + 1.96 * se_link)
      # Approximate response-scale SE, for reporting only — the CI above is
      # the authoritative interval and is NOT derived from this value.
      se_out   <- mean_out * se_link
    } else {
      mean_out <- fit_link
      lwr_out  <- fit_link - 1.96 * se_link
      upr_out  <- fit_link + 1.96 * se_link
      se_out   <- se_link
    }
    
    nd |>
      mutate(Response  = response_label,
             Timepoint = lbl,
             Mean      = round(mean_out, 3),
             SE        = round(se_out,   3),
             CI_lower  = round(lwr_out,  3),
             CI_upper  = round(upr_out,  3)) |>
      select(Response, Timepoint, !!group_var, Mean, SE, CI_lower, CI_upper)
  })
}


# ── Rigorous lpmatrix Time-Series Differences ─────────────────────────────────
# Both branches now derive the point estimate and its SE from the SAME delta-
# method construction, so diff/se, the CI, and the significance rug are mutually
# consistent. For log-scale responses the difference is taken on the response
# scale (absolute units), with gradient d/dbeta [k1*exp(eta1) - k2*exp(eta2)]
#   = k1*exp(eta1)*x1 - k2*exp(eta2)*x2 = f1*x1 - f2*x2.
pairwise_diffs_rigorous <- function(model, df_ref, pred_grid, group_var,
                                    is_log_scale = FALSE, smear = NULL) {
  
  grp_levels <- if (is.factor(pred_grid[[group_var]])) {
    levels(pred_grid[[group_var]])
  } else {
    unique(pred_grid[[group_var]])
  }
  pairs <- combn(grp_levels, 2, simplify = FALSE)
  
  fixed_var   <- ifelse(group_var == "Species", "Spacing_f", "Species")
  fixed_level <- levels(df_ref[[fixed_var]])[1]
  
  V <- vcov(model, unconditional = TRUE)
  b <- coef(model)
  
  map_dfr(pairs, function(pair) {
    lv1 <- pair[1]; lv2 <- pair[2]
    
    g1 <- pred_grid |> filter(.data[[group_var]] == lv1) |> arrange(days)
    g2 <- pred_grid |> filter(.data[[group_var]] == lv2) |> arrange(days)
    
    nd1 <- g1 |> mutate(!!fixed_var := factor(fixed_level, levels = levels(df_ref[[fixed_var]])),
                        Plot_ID = levels(df_ref$Plot_ID)[1],
                        Tree_ID = levels(df_ref$Tree_ID)[1])
    nd2 <- g2 |> mutate(!!fixed_var := factor(fixed_level, levels = levels(df_ref[[fixed_var]])),
                        Plot_ID = levels(df_ref$Plot_ID)[1],
                        Tree_ID = levels(df_ref$Tree_ID)[1])
    
    X1 <- predict(model, newdata = nd1, type = "lpmatrix", exclude = c("s(Plot_ID)", "s(Tree_ID)"))
    X2 <- predict(model, newdata = nd2, type = "lpmatrix", exclude = c("s(Plot_ID)", "s(Tree_ID)"))
    
    eta1 <- as.numeric(X1 %*% b)
    eta2 <- as.numeric(X2 %*% b)
    
    if (is_log_scale) {
      k1 <- smear_lookup(smear, nd1)
      k2 <- smear_lookup(smear, nd2)
      
      f1 <- k1 * exp(eta1)
      f2 <- k2 * exp(eta2)
      
      diff_out <- f1 - f2
      G        <- X1 * f1 - X2 * f2          # row i scaled by f1[i] / f2[i]
      se_out   <- sqrt(rowSums((G %*% V) * G))
      
    } else {
      diff_out <- eta1 - eta2
      Xdiff    <- X1 - X2
      se_out   <- sqrt(rowSums((Xdiff %*% V) * Xdiff))
    }
    
    lwr_out <- diff_out - 1.96 * se_out
    upr_out <- diff_out + 1.96 * se_out
    
    tibble(group1 = lv1, group2 = lv2,
           comparison = paste0(lv1, " - ", lv2), days = g1$days,
           diff = diff_out, se_diff = se_out,
           lwr = lwr_out, upr = upr_out,
           sig = (lwr_out > 0) | (upr_out < 0))
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
plot_diffs <- function(diff_df, y_label, fill_col = "steelblue", ncol = NULL, 
                       title = NULL, subtitle = NULL) {
  
  # Lock in factor levels to ensure the baseline stays in the 1st column
  all_levels <- unique(c(diff_df$group1, diff_df$group2))
  diff_df <- diff_df |>
    mutate(
      group1 = factor(group1, levels = all_levels),
      group2 = factor(group2, levels = all_levels)
    )
  
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
    
    # facet_grid automatically generates the lower triangle matrix layout
    facet_grid(group2 ~ group1) +
    
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    labs(x = "Days from 1 September 2025", y = y_label,
         title = title, subtitle = subtitle) +
    theme_thesis() +
    theme(
      # Align titles to the absolute left edge of the plot, not just the panel
      plot.title.position = "plot", 
      
      # Ensure hjust = 0 for strict left-justification
      plot.title       = element_text(size = 9, face = "bold", margin = margin(b = 2), hjust = 0),
      plot.subtitle    = element_text(size = 7.5, colour = "grey40", margin = margin(b = 3), hjust = 0),
      
      strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.5),
      strip.text       = element_text(face = "bold", size = 8),
      panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.5)
    )
}

# ── Statistics helpers ────────────────────────────────────────────────────────


# Extract Model Fit Statistics
extract_fit_stats <- function(model, model_name) {
  s <- summary(model)
  tibble(
    Response = model_name,
    Adj_R_squared = round(s$r.sq, 4),
    Deviance_Explained = round(s$dev.expl, 4),
    REML_Score = round(s$sp.criterion, 1),
    N_Obs = s$n
  )
}


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
# FIT MODELS ####
# ──────────────────────────────────────────────────────────────────────────────
# All three responses: Gaussian family (log-scale for Crown and CA:H)
# Each response now gets ONE combined model (Species by-smooth + Spacing
# by-smooth together), replacing the earlier species/spacing model pair.
# If models already in memory, skip to PREDICTION GRIDS
# Expected runtime: ~10-15 min per model (~30-45 min total)
# ──────────────────────────────────────────────────────────────────────────────

cat("══ RESPONSE 1: Calibrated Height (Gaussian, raw scale) ════════════════\n")
model_h <- fit_gamm_combined(df_h, "Height",
                             var_groups = c("Species", "Spacing_f"))
cat("\n── Height combined model summary ────────────────────────────────────\n")
print(summary(model_h))

cat("\n══ RESPONSE 2: Crown Area (Gaussian, log scale) ════════════════════════\n")
cat("   Fitted on log(Crown_Area_m2); predictions back-transformed to m²\n\n")
model_c <- fit_gamm_combined(df_c, "Crown",
                             var_groups = c("Species", "Spacing_f"))
cat("\n── Crown combined model summary ─────────────────────────────────────\n")
print(summary(model_c))

cat("\n══ RESPONSE 3: CA:H Ratio (Gaussian, log scale) ════════════════════════\n")
cat("   CA:H = Crown_Area_m2 / Height_m  (m2 m-1)\n")
cat("   Fitted on log(CA:H); predictions back-transformed to m2 m-1\n\n")
model_r <- fit_gamm_combined(df_r, "CAH",
                             var_groups = c("Species", "Spacing_f"))
cat("\n── CA:H combined model summary ──────────────────────────────────────\n")
print(summary(model_r))

# ──────────────────────────────────────────────────────────────────────────────
# MODEL DIAGNOSTICS ####
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

save_gam_check(model_h, "Height Combined", OUTPUT_DIR)
save_gam_check(model_c, "Crown Area Combined", OUTPUT_DIR)
save_gam_check(model_r, "CAH Ratio Combined", OUTPUT_DIR)

cat("\n")
cat("All GAM diagnostic reports saved.\n")
cat("=====================================================\n")

# ──────────────────────────────────────────────────────────────────────────────
#  PREDICTION GRIDS ####
# ──────────────────────────────────────────────────────────────────────────────

days_h <- seq(min(df_h$days), max(df_h$days), length.out = 300)
days_c <- seq(min(df_c$days), max(df_c$days), length.out = 300)
days_r <- seq(min(df_r$days), max(df_r$days), length.out = 300)

# Smearing factors — computed once per log-scale model, reused everywhere
smear_c <- smearing_factors(model_c, df_c)
smear_r <- smearing_factors(model_r, df_r)

cat("\n── Duan smearing factors (Crown) ───────────────────────────────────\n")
print(round(smear_c$factors, 4))
cat("\n── Duan smearing factors (CA:H) ────────────────────────────────────\n")
print(round(smear_r$factors, 4))

cat("Predicting height trajectories...\n")
sp_diffs_h <- pairwise_diffs_rigorous(model_h, df_h,
                                      predict_traj(model_h, make_species_grid(days_h, df_h)),
                                      "Species", is_log_scale = FALSE)
sc_diffs_h <- pairwise_diffs_rigorous(model_h, df_h,
                                      predict_traj(model_h, make_spacing_grid(days_h, df_h)),
                                      "Spacing_f", is_log_scale = FALSE)

cat("Predicting crown area trajectories (back-transforming to m2)...\n")
sp_diffs_c <- pairwise_diffs_rigorous(model_c, df_c,
                                      predict_traj(model_c, make_species_grid(days_c, df_c),
                                                   backtransform = TRUE, smear = smear_c),
                                      "Species", is_log_scale = TRUE, smear = smear_c)
sc_diffs_c <- pairwise_diffs_rigorous(model_c, df_c,
                                      predict_traj(model_c, make_spacing_grid(days_c, df_c),
                                                   backtransform = TRUE, smear = smear_c),
                                      "Spacing_f", is_log_scale = TRUE, smear = smear_c)

cat("Predicting CA:H ratio trajectories (back-transforming to m2 m-1)...\n")
sp_diffs_r <- pairwise_diffs_rigorous(model_r, df_r,
                                      predict_traj(model_r, make_species_grid(days_r, df_r),
                                                   backtransform = TRUE, smear = smear_r),
                                      "Species", is_log_scale = TRUE, smear = smear_r)
sc_diffs_r <- pairwise_diffs_rigorous(model_r, df_r,
                                      predict_traj(model_r, make_spacing_grid(days_r, df_r),
                                                   backtransform = TRUE, smear = smear_r),
                                      "Spacing_f", is_log_scale = TRUE, smear = smear_r)

cat("\nRange checks (all should be non-zero):\n")
cat("Height   species:", round(range(sp_diffs_h$diff), 3), "\n")
cat("Height   spacing:", round(range(sc_diffs_h$diff), 3), "\n")
cat("Crown    species:", round(range(sp_diffs_c$diff), 3), "\n")
cat("Crown    spacing:", round(range(sc_diffs_c$diff), 3), "\n")
cat("CA:H     species:", round(range(sp_diffs_r$diff), 3), "\n")
cat("CA:H     spacing:", round(range(sc_diffs_r$diff), 3), "\n")

# Growth curve grids
curve_h_sp <- predict_traj(model_h, make_species_grid(days_h, df_h))
curve_h_sc <- predict_traj(model_h, make_spacing_grid(days_h, df_h))

curve_c_sp <- predict_traj(model_c, make_species_grid(days_c, df_c),
                           backtransform = TRUE, smear = smear_c)
curve_c_sc <- predict_traj(model_c, make_spacing_grid(days_c, df_c),
                           backtransform = TRUE, smear = smear_c)

curve_r_sp <- predict_traj(model_r, make_species_grid(days_r, df_r),
                           backtransform = TRUE, smear = smear_r)
curve_r_sc <- predict_traj(model_r, make_spacing_grid(days_r, df_r),
                           backtransform = TRUE, smear = smear_r)

# Pearson residuals rescale by the prior weights, so under a correct variance
# structure their SD should be ~1 in EVERY group. Deviation from 1 is what's
# left uncorrected. Response residuals are also reported so you can see the
# raw (uncorrected) spread that motivated the weighting.
check_heteroscedasticity <- function(model, df, response_label) {
  
  df$resid_resp <- residuals(model, type = "response")
  df$resid_pear <- residuals(model, type = "pearson")
  
  cat("\n====", response_label, "====\n")
  
  summarise_by <- function(d, gv) {
    d |> group_by(.data[[gv]]) |>
      summarise(n = n(),
                sd_response = round(sd(resid_resp), 4),
                sd_pearson  = round(sd(resid_pear), 4),
                .groups = "drop") |>
      arrange(desc(sd_response))
  }
  
  cat("\n-- Residuals by Species --\n");  print(summarise_by(df, "Species"))
  cat("\n-- Residuals by Spacing --\n");  print(summarise_by(df, "Spacing_f"))
  
  ratio <- function(x) round(max(x) / min(x), 2)
  sp <- summarise_by(df, "Species"); sc <- summarise_by(df, "Spacing_f")
  
  cat("\n-- SD ratio (max/min) --\n")
  cat("  Species  | response:", ratio(sp$sd_response),
      " pearson:", ratio(sp$sd_pearson), "\n")
  cat("  Spacing  | response:", ratio(sc$sd_response),
      " pearson:", ratio(sc$sd_pearson), "\n")
  cat("  Target: pearson ratio near 1.0 (say < 1.2). Response ratio is\n")
  cat("  expected to stay high — that is the real heterogeneity being modelled.\n")
}

check_heteroscedasticity(model_h, df_h, "Height")
check_heteroscedasticity(model_c, df_c, "Crown Area (log)")
check_heteroscedasticity(model_r, df_r, "CA:H Ratio (log)")

# ──────────────────────────────────────────────────────────────────────────────
# PLOTS ####
# ──────────────────────────────────────────────────────────────────────────────

cat("\nGenerating plots...\n")

# ── CROWN AREA plots ──────────────────────────────────────────────────────────
p_c_sp_diff <- plot_diffs(sp_diffs_c,
                          y_label  = "Difference in crown area per tree (m\u00b2)",
                          fill_col = "#1b9e77", ncol = 3,
                          title    = "Species pairwise crown area differences",
                          subtitle = paste0(POPULATION_LABEL, " | Simulated at ", BASELINE_SPACING, " | Shaded = 95% CI | Red rug = significant period"))
ggsave(file.path(OUTPUT_DIR, "crown_species_differences.png"), p_c_sp_diff,
       width = 6.30, height = 5, units = "in", dpi = 300)

p_c_sc_diff <- plot_diffs(sc_diffs_c,
                          y_label  = "Difference in crown area per tree (m\u00b2)",
                          fill_col = "#d95f02", ncol = 3,
                          title    = "Spacing pairwise crown area differences",
                          subtitle = paste0(POPULATION_LABEL, " | Simulated for ", BASELINE_SPECIES, " | Shaded = 95% CI | Red rug = significant period"))
ggsave(file.path(OUTPUT_DIR, "crown_spacing_differences.png"), p_c_sc_diff,
       width = 4, height = 3.0, units = "in", dpi = 300)

p_c_sp_curves <- curve_plot(curve_c_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title  = "Species",
                            y_label       = "Crown area per tree (m\u00b2)",
                            title         = "GAMM-fitted crown area growth trajectories by species",
                            subtitle = paste0(POPULATION_LABEL, " | Simulated at ", BASELINE_SPACING, " | Shaded = 95% CI"))

p_c_sc_curves <- curve_plot(curve_c_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title  = "Spacing",
                            y_label       = "Crown area per tree (m\u00b2)",
                            title         = "GAMM-fitted crown area growth trajectories by spacing",
                            subtitle = paste0(POPULATION_LABEL, " | Simulated for ", BASELINE_SPECIES, " | Shaded = 95% CI"))

# ── HEIGHT plots ──────────────────────────────────────────────────────────────
p_h_sp_diff <- plot_diffs(sp_diffs_h,
                          y_label  = "Difference in calibrated height per tree (m)",
                          fill_col = "steelblue", ncol = 3,
                          title    = "Species pairwise height differences",
                          subtitle = paste0(POPULATION_LABEL, " | Simulated at ", BASELINE_SPACING, " | Shaded = 95% CI | Red rug = significant period"))
ggsave(file.path(OUTPUT_DIR, "height_species_differences.png"), p_h_sp_diff,
       width = 6.30, height = 5, units = "in", dpi = 300)

p_h_sc_diff <- plot_diffs(sc_diffs_h,
                          y_label  = "Difference in calibrated height per tree (m)",
                          fill_col = "darkorange", ncol = 3,
                          title    = "Spacing pairwise height differences",
                          subtitle = paste0(POPULATION_LABEL, " | Simulated for ", BASELINE_SPECIES, " | Shaded = 95% CI | Red rug = significant period"))
ggsave(file.path(OUTPUT_DIR, "height_spacing_differences.png"), p_h_sc_diff,
       width = 4, height = 3.0, units = "in", dpi = 300)

p_h_sp_curves <- curve_plot(curve_h_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title  = "Species",
                            y_label       = "Calibrated height per tree (m)",
                            title         = "GAMM-fitted height growth trajectories by species",
                            subtitle = paste0(POPULATION_LABEL, " | Simulated at ", BASELINE_SPACING, " | Shaded = 95% CI"))

p_h_sc_curves <- curve_plot(curve_h_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title  = "Spacing",
                            y_label       = "Calibrated height per tree (m)",
                            title         = "GAMM-fitted height growth trajectories by spacing",
                            subtitle = paste0(POPULATION_LABEL, " | Simulated for ", BASELINE_SPECIES, " | Shaded = 95% CI"))


# ── CA:H RATIO plots ──────────────────────────────────────────────────────────
p_r_sp_diff <- plot_diffs(sp_diffs_r,
                          y_label  = "Difference in CA:H ratio per tree (m\u00b2 m\u207b\u00b9)",
                          fill_col = "#6a3d9a", ncol = 3,
                          title    = "Species pairwise CA:H ratio differences",
                          subtitle = paste0(POPULATION_LABEL, " | Simulated at ", BASELINE_SPACING, " | Shaded = 95% CI | Red rug = significant period"))
ggsave(file.path(OUTPUT_DIR, "cah_species_differences.png"), p_r_sp_diff,
       width = 6.30, height = 5, units = "in", dpi = 300)

p_r_sc_diff <- plot_diffs(sc_diffs_r,
                          y_label  = "Difference in CA:H ratio per tree (m\u00b2 m\u207b\u00b9)",
                          fill_col = "#e31a1c", ncol = 3,
                          title    = "Spacing pairwise CA:H ratio differences",
                          subtitle = paste0(POPULATION_LABEL, " | Simulated for ", BASELINE_SPECIES, " | Shaded = 95% CI | Red rug = significant period"))
ggsave(file.path(OUTPUT_DIR, "cah_spacing_differences.png"), p_r_sc_diff,
       width = 4, height = 3.0, units = "in", dpi = 300)

p_r_sp_curves <- curve_plot(curve_r_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title  = "Species",
                            y_label       = "CA:H ratio per tree (m\u00b2 m\u207b\u00b9)",
                            title         = "GAMM-fitted CA:H ratio trajectories by species",
                            subtitle = paste0(POPULATION_LABEL, " | Simulated at ", BASELINE_SPACING, " | Shaded = 95% CI"))

p_r_sc_curves <- curve_plot(curve_r_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title  = "Spacing",
                            y_label       = "CA:H ratio per tree (m\u00b2 m\u207b\u00b9)",
                            title         = "GAMM-fitted CA:H ratio trajectories by spacing",
                            subtitle = paste0(POPULATION_LABEL, " | Simulated for ", BASELINE_SPECIES, " | Shaded = 95% CI"))

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



# 1. CROWN AREA LABELS
lbl_c_sp <- attach_labels_auto(
  curve_df    = curve_c_sp,
  diff_df     = sp_diffs_c,
  group_var   = "Species",
  target_day  = 266,
  y_positions = c(1.15, 1.00, 0.85, 0.65, 0.50)
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
  y_positions = c(3.0, 2.5, 2.00, 0.90)
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

# 2. HEIGHT LABELS
lbl_h_sp <- attach_labels_auto(
  curve_df    = curve_h_sp,
  diff_df     = sp_diffs_h,
  group_var   = "Species",
  target_day  = 266,
  y_positions = c(4.4, 4, 3.6, 2.85, 2.05)
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

# 3. CA:H RATIO LABELS
lbl_r_sp <- attach_labels_auto(
  curve_df    = curve_r_sp,
  diff_df     = sp_diffs_r,
  group_var   = "Species",
  target_day  = 266,
  y_positions = c(0.36, 0.32, 0.28, 0.235, 0.195)
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
  y_positions = c(0.77, 0.67, 0.535, 0.23)
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
                        keep_subtitle = FALSE, tag = NULL) {
  p <- p + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
  
  if (!keep_legend)  p <- p + theme(legend.position  = "none")
  if (!keep_x)       p <- p + theme(axis.title.x     = element_blank(),
                                    axis.text.x      = element_blank(),
                                    axis.ticks.x     = element_blank())
  if (!keep_y)       p <- p + theme(axis.title.y     = element_blank())
  
  p <- p + theme(plot.title = element_blank())
  if (!keep_subtitle) p <- p + theme(plot.subtitle = element_blank())
  
  # Manual panel tag, anchored to the panel's own top-left corner via
  # -Inf/Inf — same annotate() + family/fontface combo already confirmed
  # to render bold correctly in font_test.png.
  if (!is.null(tag)) {
    p <- p + annotate("text", x = -Inf, y = Inf, label = tag,
                      family = "Calibri", fontface = "bold",
                      size = 4.5, hjust = -0.3, vjust = 1.5)
  }
  
  return(p)
}

# 2. Apply cleaning to all 6 panels
#    Only c_sp_3x2 (top-left) keeps its subtitle — it reads "Left column: Species-controlled"
#    Re-write that subtitle to serve as a column-header hint for the reader
c_sp_3x2 <- clean_panel(p_c_sp_curves, keep_legend = TRUE,  keep_x = FALSE,
                        keep_y = TRUE,  keep_subtitle = FALSE, tag = "(a)")
c_sc_3x2 <- clean_panel(p_c_sc_curves, keep_legend = TRUE,  keep_x = FALSE,
                        keep_y = FALSE, keep_subtitle = FALSE, tag = "(b)")
h_sp_3x2 <- clean_panel(p_h_sp_curves, keep_legend = FALSE, keep_x = FALSE,
                        keep_y = TRUE,  keep_subtitle = FALSE, tag = "(c)")
h_sc_3x2 <- clean_panel(p_h_sc_curves, keep_legend = FALSE, keep_x = FALSE,
                        keep_y = FALSE, keep_subtitle = FALSE, tag = "(d)")
r_sp_3x2 <- clean_panel(p_r_sp_curves, keep_legend = FALSE, keep_x = TRUE,
                        keep_y = TRUE,  keep_subtitle = FALSE, tag = "(e)")
r_sc_3x2 <- clean_panel(p_r_sc_curves, keep_legend = FALSE, keep_x = TRUE,
                        keep_y = FALSE, keep_subtitle = FALSE, tag = "(f)")

# 3. Assemble and annotate with a single main title
p_combined_3x2 <- (c_sp_3x2 | c_sc_3x2) /
  (h_sp_3x2   | h_sc_3x2)   /
  (r_sp_3x2   | r_sc_3x2)   +
  plot_annotation(
    title = paste0("GAMM-fitted growth trajectories \u2014 ", POPULATION_LABEL),
    subtitle = paste0("Left: Species comparisons (Simulated at ", BASELINE_SPACING, ")  |  Right: Spacing comparisons (Simulated for ", BASELINE_SPECIES, ")"),
    theme = theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.5, colour = "grey40", hjust = 0.5, margin = margin(b = 4))
    )
  )

# Height bumped by 0.2 in to give the main title breathing room
ggsave(
  file.path(OUTPUT_DIR, "combined_3x2_curves.png"), p_combined_3x2,
  width = 6.30, height = 6.7, units = "in", dpi = 300,
  device = ragg::agg_png
)

# ──────────────────────────────────────────────────────────────────────────────
# STATISTICS SUMMARY TABLES ####
# ──────────────────────────────────────────────────────────────────────────────

cat("\nGenerating statistics tables...\n")

key_days   <- c(59, 134, 203, 266)
key_labels <- c("2 months", "4.5 months", "6.5 months", "9 months")
final_day  <- 266

# TABLE 1: Smooth term significance — one combined model per response
# (each response's table now includes BOTH Species and Spacing by-smooth terms)
tbl1 <- bind_rows(
  smooth_sig_table(model_h, "Height"),
  smooth_sig_table(model_c, "Crown Area"),
  smooth_sig_table(model_r, "CA:H Ratio")
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
tbl2_h_sp <- marginal_means(model_h, "Species",   levels(df_h$Species),
                            "Spacing_f", levels(df_h$Spacing_f)[1],
                            df_h, key_days, key_labels, "Height (m)")
tbl2_h_sc <- marginal_means(model_h, "Spacing_f", levels(df_h$Spacing_f),
                            "Species",   levels(df_h$Species)[1],
                            df_h, key_days, key_labels, "Height (m)")

# Crown — back-transform from log scale to m²
tbl2_c_sp <- marginal_means(model_c, "Species",   levels(df_c$Species),
                            "Spacing_f", levels(df_c$Spacing_f)[1],
                            df_c, key_days, key_labels, "Crown Area (m2)",
                            backtransform = TRUE, smear = smear_c)
tbl2_c_sc <- marginal_means(model_c, "Spacing_f", levels(df_c$Spacing_f),
                            "Species",   levels(df_c$Species)[1],
                            df_c, key_days, key_labels, "Crown Area (m2)",
                            backtransform = TRUE, smear = smear_c)

# CA:H — back-transform from log scale to m² m⁻¹
tbl2_r_sp <- marginal_means(model_r, "Species",   levels(df_r$Species),
                            "Spacing_f", levels(df_r$Spacing_f)[1],
                            df_r, key_days, key_labels, "CA:H Ratio (m2 m-1)",
                            backtransform = TRUE, smear = smear_r)
tbl2_r_sc <- marginal_means(model_r, "Spacing_f", levels(df_r$Spacing_f),
                            "Species",   levels(df_r$Species)[1],
                            df_r, key_days, key_labels, "CA:H Ratio (m2 m-1)",
                            backtransform = TRUE, smear = smear_r)

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

# TABLE 4: Overall Model Fit Statistics
tbl_fit <- bind_rows(
  extract_fit_stats(model_h, "Height"),
  extract_fit_stats(model_c, "Crown Area"),
  extract_fit_stats(model_r, "CA:H Ratio")
)

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TABLE 4: Overall Model Fit Statistics                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
print(tbl_fit, n = Inf)

write_csv(
  tbl_fit, 
  file.path(OUTPUT_DIR, "stats_table4_model_fits.csv")
)


# ──────────────────────────────────────────────────────────────────────────────
# ── DIAGNOSTICS ───────────────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────

cat("\n── k adequacy checks ────────────────────────────────────────────────\n")
model_list <- list(
  "Height combined" = model_h,
  "Crown combined"  = model_c,
  "CA:H combined"   = model_r
)
for (nm in names(model_list)) {
  cat(nm, "model:\n"); print(k.check(model_list[[nm]]))
}

diag_files <- list(
  list(model_h, file.path(OUTPUT_DIR, "diag_height_combined.png"), "steelblue"),
  list(model_c, file.path(OUTPUT_DIR, "diag_crown_combined.png"), "#1b9e77"),
  list(model_r, file.path(OUTPUT_DIR, "diag_cah_combined.png"), "#6a3d9a")
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
cat("    diag_height/crown/cah _combined.png\n")
cat("\n── Done ──────────────────────────────────────────────────────────────\n")

# End counter
toc(func.toc = toc_in_mins)