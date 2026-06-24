# =============================================================================
# GAMM — Calibrated Height, Crown Area & Crown:Height Ratio
# Eucalyptus species × spacing trial
#
# Time series: 1 September 2025 onwards
#
# Three response variables:
#   Response 1: Calibrated_Height_m  — Gaussian  (identity link, raw scale)
#   Response 2: Crown_Area_m2        — Gaussian  (log-transformed; back-transformed for plots)
#   Response 3: CA:H ratio           — Gaussian  (log-transformed; back-transformed for plots)
#
# WHY LOG TRANSFORMATION FOR CROWN AND CA:H:
#   Gamma(log) failed to converge after 4+ hours because log(Crown) and
#   log(CA:H) are LEFT-skewed (skew ≈ -0.6 to -1.5 by spacing), meaning
#   the Gaussian on the log scale is a better fit than Gamma on original scale.
#   Predictions are back-transformed to m² and m² m⁻¹ using the delta method:
#   E[exp(y)] ≈ exp(fit + 0.5 * se²), SE ≈ exp(fit) * se
#
# For each response, two models:
#   m_species : s(days, by=Species)   + Spacing_f fixed + random effects
#   m_spacing : s(days, by=Spacing_f) + Species fixed   + random effects
#
# Global smooth uses bs='cr' (cubic regression spline) to prevent boundary
# overshoot caused by the 59-day gap between Sep 1 and Oct 30.
#
# Spacing codes:  1x1m = 1 m²/tree  |  2x2m = 4 m²/tree
#                 3x3m = 9 m²/tree  |  5x5m = 25 m²/tree
# =============================================================================


# ── 0. Packages ───────────────────────────────────────────────────────────────
# install.packages(c("mgcv", "tidyverse", "gratia", "patchwork"))
library(mgcv)
library(tidyverse)
library(gratia)
library(patchwork)
library(sf)
library(tictoc)

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


# =============================================================================
# OUTPUT SETTINGS
# =============================================================================

OUTPUT_DIR <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/10. GAMM"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# =============================================================================
# ANALYSIS SETTINGS
# =============================================================================

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
POPULATION_SUBSET <- "Dominant" 

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
# 1. Load one shapefile to serve as the master coordinate reference
base_crowns <- st_read("E:/Remote Sensing Media/02. 01 September 2025/09. Crown Metrics/Crown_Metrics_01_September_2025.shp", quiet = TRUE) %>%
  mutate(Tree = round(as.numeric(Tree), 2))

# 2. Extract X and Y centroids 
coords_df <- base_crowns %>%
  st_centroid() %>%
  mutate(
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(Compartment = Cmprtmn, Line, Plot, Tree, X, Y)

# 3. Merge X and Y into your main dataset
df_raw <- read_csv("C:/Users/jakev/Downloads/UAV_Master_Dataset_25-05-2026.csv", show_col_types = FALSE) %>%
  mutate(Tree = round(as.numeric(Tree), 2)) %>%
  left_join(coords_df, by = c("Compartment", "Line", "Plot", "Tree"))

df_base <- df_raw |>
  mutate(
    Date      = as.Date(Date),
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

# ── NEW: Dominant Tree Filter Switch ──
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


# =============================================================================
# ── SHARED FUNCTIONS ──────────────────────────────────────────────────────────
# =============================================================================

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
    exclude = c(
      "s(Plot_ID)",
      "s(Tree_ID)"
    )
  )
  
  fit_raw <- as.numeric(preds$fit)
  se_raw  <- as.numeric(preds$se.fit)
  
  if (backtransform) {
    
    fit_out <- exp(fit_raw + 0.5 * se_raw^2)
    se_out  <- exp(fit_raw) * se_raw
    
  } else {
    
    fit_out <- fit_raw
    se_out  <- se_raw
    
  }
  
  newdata |>
    mutate(
      fit = fit_out,
      se  = se_out
    )
}


# ── All pairwise differences ───────────────────────────────────────────────────
pairwise_diffs <- function(pred_grid, group_var, average_over) {
  grp_levels <- if (is.factor(pred_grid[[group_var]])) {
    levels(pred_grid[[group_var]])
  } else {
    unique(pred_grid[[group_var]])
  }
  pairs <- combn(grp_levels, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    lv1 <- pair[1]; lv2 <- pair[2]
    
    g1 <- pred_grid |> filter(.data[[group_var]] == lv1) |>
      group_by(days) |> summarise(fit = mean(fit), se = mean(se), .groups = "drop")
    g2 <- pred_grid |> filter(.data[[group_var]] == lv2) |>
      group_by(days) |> summarise(fit = mean(fit), se = mean(se), .groups = "drop")
    
    tibble(
      comparison = paste0(lv1, " - ", lv2),
      days       = g1$days,
      diff       = g1$fit - g2$fit,
      se_diff    = sqrt(g1$se^2 + g2$se^2),
      lwr        = diff - 1.96 * se_diff,
      upr        = diff + 1.96 * se_diff,
      sig        = (lwr > 0) | (upr < 0)
    )
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
      
      # Legend Settings (Upper left, framed, slightly transparent white background)
      legend.position      = c(0.02, 0.98), 
      legend.justification = c(0, 1),
      legend.background    = element_rect(fill = alpha("white", 0.95), colour = "lightgray", linewidth = 0.5),
      legend.title         = element_text(size = 10, face = "bold"),
      legend.text          = element_text(size = 8),
      legend.key.size      = unit(0.4, "cm")
    )
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
# Note: title and subtitle arguments have been removed to save vertical space
curve_plot <- function(curve_df, colour_var, colour_vals, colour_labels = NULL,
                       y_label, legend_title = NULL) {
  ggplot(curve_df, aes(x = days,
                       colour = .data[[colour_var]],
                       fill   = .data[[colour_var]])) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
    geom_line(aes(y = fit), linewidth = 0.9) +
    scale_colour_manual(
      values = colour_vals,
      labels = colour_labels,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = colour_vals,
      labels = colour_labels,
      drop = FALSE
    ) +
    
    # Keep 2-row legend so it doesn't stretch too wide across the top
    guides(
      colour = guide_legend(nrow = 2), 
      fill = "none"
    ) +
    scale_x_continuous(breaks = seq(0, 270, by = 60)) +
    
    # Padding returned to standard 5% since the legend is now outside
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    
    labs(x = "Days from 1 September 2025", y = y_label,
         colour = legend_title, fill = legend_title) +
    theme_thesis()
}

# ── Plot pairwise differences ─────────────────────────────────────────────────
plot_diffs <- function(diff_df, y_label, fill_col = "steelblue", ncol = 3) {
  
  # 1. Calculate a single global Y position so rugs align perfectly across all facets
  global_min_y <- min(diff_df$lwr, na.rm = TRUE)
  global_max_y <- max(diff_df$upr, na.rm = TRUE)
  
  # Set the rug position slightly below the lowest overall confidence interval
  rug_y_pos <- global_min_y - ((global_max_y - global_min_y) * 0.05)
  
  ggplot(diff_df, aes(x = days, y = diff)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = fill_col) +
    geom_line(colour = fill_col, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    
    # 2. The Floating Rug (using vertical lines, shape = 124)
    # Adjust `size = 6` up or down to make the rugs bigger or smaller!
    geom_point(data = diff_df[diff_df$sig, ],
               aes(x = days, y = rug_y_pos),
               colour = "#ff3333", alpha = 0.8,
               shape = 124, size = 2) +
    
    facet_wrap(~ comparison, ncol = ncol) +
    
    # 3. Crucial: Expand the bottom of the Y-axis by 15% so the floating rugs have room
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    
    labs(x = "Days from 1 September 2025", y = y_label) +
    theme_thesis()
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


marginal_means <- function(model, group_var, group_levels,
                           fixed_var, fixed_level,
                           df_ref, key_days, key_labels, response_label,
                           backtransform = FALSE) {
  
  # 1. Define reference X and Y coordinates
  x_ref <- mean(df_ref$X, na.rm = TRUE)
  y_ref <- mean(df_ref$Y, na.rm = TRUE)
  
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

# =============================================================================
# ── FIT MODELS ────────────────────────────────────────────────────────────────
# All three responses: Gaussian family (log-scale for Crown and CA:H)
# If models already in memory, skip to PREDICTION GRIDS
# Expected runtime: ~5-10 min per model pair (~30 min total)
# =============================================================================

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

# =============================================================================
# MODEL DIAGNOSTICS
# =============================================================================

cat("\n")
cat("=====================================================\n")
cat("MODEL DIAGNOSTICS\n")
cat("=====================================================\n")

# -----------------------------------------------------------------------------
# Helper: run gam.check() and save output
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# Run and save GAM checks
# -----------------------------------------------------------------------------

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

# =============================================================================
# ── PREDICTION GRIDS ──────────────────────────────────────────────────────────
# =============================================================================

days_h <- seq(min(df_h$days), max(df_h$days), length.out = 300)
days_c <- seq(min(df_c$days), max(df_c$days), length.out = 300)
days_r <- seq(min(df_r$days), max(df_r$days), length.out = 300)

cat("Predicting height trajectories...\n")
grid_h_sp_pred <- predict_traj(models_h$species, make_full_grid(days_h, df_h))
grid_h_sc_pred <- predict_traj(models_h$spacing, make_full_grid(days_h, df_h))

cat("Predicting crown area trajectories (back-transforming to m2)...\n")
grid_c_sp_pred <- predict_traj(models_c$species, make_full_grid(days_c, df_c),
                               backtransform = TRUE)
grid_c_sc_pred <- predict_traj(models_c$spacing, make_full_grid(days_c, df_c),
                               backtransform = TRUE)

cat("Predicting CA:H ratio trajectories (back-transforming to m2 m-1)...\n")
grid_r_sp_pred <- predict_traj(models_r$species, make_full_grid(days_r, df_r),
                               backtransform = TRUE)
grid_r_sc_pred <- predict_traj(models_r$spacing, make_full_grid(days_r, df_r),
                               backtransform = TRUE)

# Pairwise differences
cat("Computing pairwise differences...\n")
sp_diffs_h <- pairwise_diffs(grid_h_sp_pred, "Species",   "Spacing_f")
sc_diffs_h <- pairwise_diffs(grid_h_sc_pred, "Spacing_f", "Species")
sp_diffs_c <- pairwise_diffs(grid_c_sp_pred, "Species",   "Spacing_f")
sc_diffs_c <- pairwise_diffs(grid_c_sc_pred, "Spacing_f", "Species")
sp_diffs_r <- pairwise_diffs(grid_r_sp_pred, "Species",   "Spacing_f")
sc_diffs_r <- pairwise_diffs(grid_r_sc_pred, "Spacing_f", "Species")

cat("\nRange checks (all should be non-zero):\n")
cat("Height   species:", round(range(sp_diffs_h$diff), 3), "\n")
cat("Height   spacing:", round(range(sc_diffs_h$diff), 3), "\n")
cat("Crown    species:", round(range(sp_diffs_c$diff), 3), "\n")
cat("Crown    spacing:", round(range(sc_diffs_c$diff), 3), "\n")
cat("CA:H     species:", round(range(sp_diffs_r$diff), 3), "\n")
cat("CA:H     spacing:", round(range(sc_diffs_r$diff), 3), "\n")

# Growth curve grids
curve_h_sp <- predict_traj(models_h$species, make_species_grid(days_h, df_h)) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)
curve_h_sc <- predict_traj(models_h$spacing, make_spacing_grid(days_h, df_h)) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)

curve_c_sp <- predict_traj(models_c$species, make_species_grid(days_c, df_c),
                           backtransform = TRUE) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)
curve_c_sc <- predict_traj(models_c$spacing, make_spacing_grid(days_c, df_c),
                           backtransform = TRUE) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)

curve_r_sp <- predict_traj(models_r$species, make_species_grid(days_r, df_r),
                           backtransform = TRUE) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)
curve_r_sc <- predict_traj(models_r$spacing, make_spacing_grid(days_r, df_r),
                           backtransform = TRUE) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)


# =============================================================================
# ── PLOTS ─────────────────────────────────────────────────────────────────────
# =============================================================================

# =============================================================================
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
# =============================================================================

cat("\nGenerating plots...\n")

# ── HEIGHT plots ──────────────────────────────────────────────────────────────
p_h_sp_diff <- plot_diffs(sp_diffs_h,
                          y_label = "Difference in height (m)", fill_col = "steelblue", ncol = 3) 
ggsave(file.path(OUTPUT_DIR, "height_species_differences.png"), p_h_sp_diff,
       width = 6.30, height = 5.5, units = "in", dpi = 300)

p_h_sc_diff <- plot_diffs(sc_diffs_h,
                          y_label = "Difference in height (m)", fill_col = "darkorange", ncol = 3) 
ggsave(file.path(OUTPUT_DIR, "height_spacing_differences.png"), p_h_sc_diff,
       width = 6.30, height = 3.0, units = "in", dpi = 300)

p_h_sp_curves <- curve_plot(curve_h_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title = "Species",
                            y_label = "Calibrated height (m)")
ggsave(file.path(OUTPUT_DIR, "height_curves_species.png"), p_h_sp_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)

p_h_sc_curves <- curve_plot(curve_h_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title = "Spacing",
                            y_label = "Calibrated height (m)")
ggsave(file.path(OUTPUT_DIR, "height_curves_spacing.png"), p_h_sc_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)

# ── CROWN AREA plots ──────────────────────────────────────────────────────────
p_c_sp_diff <- plot_diffs(sp_diffs_c,
                          y_label = "Difference in crown area (m2)", fill_col = "#1b9e77", ncol = 3) 
ggsave(file.path(OUTPUT_DIR, "crown_species_differences.png"), p_c_sp_diff,
       width = 6.30, height = 5.5, units = "in", dpi = 300)

p_c_sc_diff <- plot_diffs(sc_diffs_c,
                          y_label = "Difference in crown area (m2)", fill_col = "#d95f02", ncol = 3) 
ggsave(file.path(OUTPUT_DIR, "crown_spacing_differences.png"), p_c_sc_diff,
       width = 6.30, height = 3.0, units = "in", dpi = 300)

p_c_sp_curves <- curve_plot(curve_c_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title = "Species",
                            y_label = "Crown area (m2)")
ggsave(file.path(OUTPUT_DIR, "crown_curves_species.png"), p_c_sp_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)

p_c_sc_curves <- curve_plot(curve_c_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title = "Spacing",
                            y_label = "Crown area (m2)")
ggsave(file.path(OUTPUT_DIR, "crown_curves_spacing.png"), p_c_sc_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)

# ── CA:H RATIO plots ──────────────────────────────────────────────────────────
p_r_sp_diff <- plot_diffs(sp_diffs_r,
                          y_label = "Difference in CA:H ratio (m2 m-1)", fill_col = "#6a3d9a", ncol = 3) 
ggsave(file.path(OUTPUT_DIR, "cah_species_differences.png"), p_r_sp_diff,
       width = 6.30, height = 5.5, units = "in", dpi = 300)

p_r_sc_diff <- plot_diffs(sc_diffs_r,
                          y_label = "Difference in CA:H ratio (m2 m-1)", fill_col = "#e31a1c", ncol = 3) 
ggsave(file.path(OUTPUT_DIR, "cah_spacing_differences.png"), p_r_sc_diff,
       width = 6.30, height = 3.0, units = "in", dpi = 300)

p_r_sp_curves <- curve_plot(curve_r_sp, "Species", species_colors,
                            colour_labels = species_display,
                            legend_title = "Species",
                            y_label = "CA:H ratio (m2 m-1)")
ggsave(file.path(OUTPUT_DIR, "cah_curves_species.png"), p_r_sp_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)

p_r_sc_curves <- curve_plot(curve_r_sc, "Spacing_f", spacing_colors,
                            colour_labels = spacing_display,
                            legend_title = "Spacing",
                            y_label = "CA:H ratio (m2 m-1)")
ggsave(file.path(OUTPUT_DIR, "cah_curves_spacing.png"), p_r_sc_curves,
       width = 6.30, height = 3.2, units = "in", dpi = 300)

# ── Combined 3x2 Figure (Max Width) ───────────────────────────────────────────

cat("\nAssembling 3x2 combined grid (independent Y-axes)...\n")

# 1. Helper function to strip X-axes on inner plots
clean_panel <- function(p, keep_legend = FALSE, keep_x = FALSE, keep_y = TRUE) {
  p <- p + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5)) 
  
  if (!keep_legend) p <- p + theme(legend.position = "none")
  
  if (!keep_x) {
    p <- p + theme(axis.title.x = element_blank(), 
                   axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank())
  }
  
  if (!keep_y) {
    p <- p + theme(axis.title.y = element_blank(), 
                   axis.text.y = element_blank(), 
                   axis.ticks.y = element_blank())
  }
  
  return(p)
}

# 2. Apply cleaning to all 6 panels
h_sp_3x2 <- clean_panel(p_h_sp_curves, keep_legend = TRUE, keep_x = FALSE, keep_y = TRUE)
h_sc_3x2 <- clean_panel(p_h_sc_curves, keep_legend = TRUE, keep_x = FALSE, keep_y = TRUE)

c_sp_3x2 <- clean_panel(p_c_sp_curves, keep_legend = FALSE, keep_x = FALSE, keep_y = TRUE)
c_sc_3x2 <- clean_panel(p_c_sc_curves, keep_legend = FALSE, keep_x = FALSE, keep_y = TRUE)

r_sp_3x2 <- clean_panel(p_r_sp_curves, keep_legend = FALSE, keep_x = TRUE,  keep_y = TRUE)
r_sc_3x2 <- clean_panel(p_r_sc_curves, keep_legend = FALSE, keep_x = TRUE,  keep_y = TRUE)

# 3. Assemble the 3x2 grid with patchwork (Removed plot_annotation titles)
p_combined_3x2 <- (h_sp_3x2 | h_sc_3x2) / 
  (c_sp_3x2 | c_sc_3x2) / 
  (r_sp_3x2 | r_sc_3x2)

# Height compressed to 6.5 inches
ggsave(file.path(OUTPUT_DIR, "combined_3x2_curves.png"), p_combined_3x2,
       width = 6.30, height = 6.5, units = "in", dpi = 300)

# =============================================================================
# ── STATISTICS SUMMARY TABLES ─────────────────────────────────────────────────
# =============================================================================

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


# =============================================================================
# ── DIAGNOSTICS ───────────────────────────────────────────────────────────────
# =============================================================================

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
