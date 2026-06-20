# =============================================================================
# GAMM — Calibrated Height, Crown Area & Crown:Height Ratio
# Eucalyptus species × spacing trial
#
# Time series: 1 September 2025 onwards
# Three response variables, each with two models:
#   Response 1: Calibrated_Height_m  — Gaussian, cubic regression spline
#   Response 2: Crown_Area_m2        — Gamma(log), crown expansion
#   Response 3: CA_H_ratio           — Gamma(log), crown architecture index
#
# For each response:
#   m_species : s(days, by=Species)   + Spacing_f fixed + random effects
#   m_spacing : s(days, by=Spacing_f) + Species fixed   + random effects
#
# Pairwise differences computed from predictions (bypasses gratia bug).
#
# Spacing codes:  1x1m = 1 m²/tree  |  2x2m = 4 m²/tree
#                 3x3m = 9 m²/tree  |  5x5m = 25 m²/tree
#
# CA:H ratio = Crown_Area_m2 / Calibrated_Height_m  (m² m⁻¹)
#   High ratio = wide spreading crown relative to height
#   Low ratio  = narrow columnar form
#   Expected gradient: wider spacing → higher CA:H (more crown per unit height)
# =============================================================================


# ── 0. Packages ───────────────────────────────────────────────────────────────
# install.packages(c("mgcv", "tidyverse", "gratia", "patchwork"))
library(mgcv)
library(tidyverse)
library(gratia)
library(patchwork)


# ── 1. Load & clean data ──────────────────────────────────────────────────────
df_raw <- read_csv("C:/Users/jakev/Downloads/UAV_Master_Dataset_25-05-2026.csv",
                   show_col_types = FALSE)

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
  filter(Date >= as.Date("2025-09-01"))

# Height dataset — Gaussian
df_h <- df_base |>
  filter(!is.na(Calibrated_Height_m)) |>
  rename(Height = Calibrated_Height_m)

# Crown area dataset — Gamma (strictly positive)
df_c <- df_base |>
  filter(!is.na(Crown_Area_m2)) |>
  filter(Crown_Area_m2 > 0) |>
  rename(Crown = Crown_Area_m2)

# CA:H ratio dataset — Gamma (strictly positive)
# Requires both crown area AND height to be non-missing
df_r <- df_base |>
  filter(!is.na(Crown_Area_m2), !is.na(Calibrated_Height_m)) |>
  filter(Crown_Area_m2 > 0, Calibrated_Height_m > 0) |>
  mutate(CAH = Crown_Area_m2 / Calibrated_Height_m)

cat("── Data summary ──────────────────────────────────────────────────────\n")
cat("Height dataset:     ", nrow(df_h), "obs |", n_distinct(df_h$Tree_ID), "trees\n")
cat("Crown area dataset: ", nrow(df_c), "obs |", n_distinct(df_c$Tree_ID), "trees\n")
cat("CA:H ratio dataset: ", nrow(df_r), "obs |", n_distinct(df_r$Tree_ID), "trees\n")
cat("Days range:         ", min(df_h$days), "to", max(df_h$days), "\n")
cat("CA:H range:         ", round(min(df_r$CAH), 3), "to", round(max(df_r$CAH), 3), "\n\n")


# ── 2. Shared plot settings ───────────────────────────────────────────────────
species_colors <- c(
  "Cladocalyx"    = "#336998",
  "Grandis"       = "#97dde3",
  "Cloeziana"     = "#ffffff",
  "Urophylla"     = "#e3acff",
  "Grandis clone" = "#ff7da0"
)

spacing_colors <- c(
  "1x1m" = "#118AB2",
  "2x2m" = "#EF476F",
  "3x3m" = "#FFD166",
  "5x5m" = "#06D6A0"
)

spacing_display <- c(
  "1x1m" = "1x1 m (1 m2/tree)",
  "2x2m" = "2x2 m (4 m2/tree)",
  "3x3m" = "3x3 m (9 m2/tree)",
  "5x5m" = "5x5 m (25 m2/tree)"
)


# =============================================================================
# ── SHARED FUNCTIONS ──────────────────────────────────────────────────────────
# =============================================================================

fit_gamm_pair <- function(df, response_col, family = gaussian()) {
  
  cat("  Fitting species model...\n")
  m_sp <- bam(
    as.formula(paste0(response_col, " ~
      s(days, k = 12, bs = 'cr') +
      s(days, by = Species,   k = 10, bs = 'tp') +
      Spacing_f + Culture +
      s(Plot_ID, bs = 're') +
      s(Tree_ID, bs = 're')")),
    data     = df,
    family   = family,
    method   = "fREML",
    discrete = FALSE
  )
  
  cat("  Fitting spacing model...\n")
  m_sc <- bam(
    as.formula(paste0(response_col, " ~
      s(days, k = 12, bs = 'cr') +
      s(days, by = Spacing_f, k = 10, bs = 'tp') +
      Species + Culture +
      s(Plot_ID, bs = 're') +
      s(Tree_ID, bs = 're')")),
    data     = df,
    family   = family,
    method   = "fREML",
    discrete = FALSE
  )
  
  list(species = m_sp, spacing = m_sc)
}


predict_traj <- function(model, newdata) {
  preds <- predict(model, newdata = newdata, se.fit = TRUE,
                   type    = "response",
                   exclude = c("s(Plot_ID)", "s(Tree_ID)"))
  newdata |>
    mutate(fit = as.numeric(preds$fit),
           se  = as.numeric(preds$se.fit))
}


pairwise_diffs <- function(pred_grid, group_var, average_over) {
  grp_levels <- if (is.factor(pred_grid[[group_var]])) {
    levels(pred_grid[[group_var]])
  } else {
    unique(pred_grid[[group_var]])
  }
  pairs <- combn(grp_levels, 2, simplify = FALSE)
  
  map_dfr(pairs, function(pair) {
    lv1 <- pair[1]; lv2 <- pair[2]
    
    g1 <- pred_grid |>
      filter(.data[[group_var]] == lv1) |>
      group_by(days) |>
      summarise(fit = mean(fit), se = mean(se), .groups = "drop")
    
    g2 <- pred_grid |>
      filter(.data[[group_var]] == lv2) |>
      group_by(days) |>
      summarise(fit = mean(fit), se = mean(se), .groups = "drop")
    
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


plot_diffs <- function(diff_df, title, subtitle, y_label,
                       fill_col = "steelblue", ncol = 3) {
  ggplot(diff_df, aes(x = days, y = diff)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                alpha = 0.2, fill = fill_col) +
    geom_line(colour = fill_col, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_rug(data = diff_df[diff_df$sig, ],
             sides = "b", colour = "red", alpha = 0.5) +
    facet_wrap(~ comparison, ncol = ncol) +
    labs(title = title, subtitle = subtitle,
         x = "Days from 1 September 2025",
         y = y_label) +
    theme_classic(base_size = 9) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(size = 7, face = "bold"),
      plot.title       = element_text(size = 9, face = "bold"),
      plot.subtitle    = element_text(size = 7, colour = "grey40")
    )
}


smooth_sig_table <- function(model, response_label) {
  s      <- summary(model)
  sp_tbl <- as.data.frame(s$s.table)
  sp_tbl <- sp_tbl[grepl("^s\\(days\\):", rownames(sp_tbl)), ]
  sp_tbl$Term     <- rownames(sp_tbl)
  sp_tbl$Response <- response_label
  sp_tbl$Sig      <- ifelse(sp_tbl[["p-value"]] < 0.001, "***",
                            ifelse(sp_tbl[["p-value"]] < 0.01,  "**",
                                   ifelse(sp_tbl[["p-value"]] < 0.05,  "*",
                                          ifelse(sp_tbl[["p-value"]] < 0.1,   ".",  "ns"))))
  sp_tbl |>
    select(Response, Term, edf = edf, F = F, p = `p-value`, Sig) |>
    mutate(edf = round(edf, 2), F = round(F, 3),
           p   = sprintf("%.2e", p))
}


marginal_means <- function(model, group_var, group_levels,
                           fixed_var, fixed_level,
                           df_ref, key_days, key_labels, response_label) {
  pred_base <- tibble(
    !!group_var := factor(group_levels, levels = levels(df_ref[[group_var]])),
    !!fixed_var := factor(fixed_level,  levels = levels(df_ref[[fixed_var]])),
    Culture = factor("Single", levels = levels(df_ref$Culture)),
    Plot_ID = levels(df_ref$Plot_ID)[1],
    Tree_ID = levels(df_ref$Tree_ID)[1]
  )
  map2_dfr(key_days, key_labels, function(d, lbl) {
    nd    <- pred_base |> mutate(days = d)
    preds <- predict(model, newdata = nd, se.fit = TRUE,
                     type    = "response",
                     exclude = c("s(Plot_ID)", "s(Tree_ID)"))
    nd |>
      mutate(Response  = response_label,
             Timepoint = lbl,
             Mean      = round(as.numeric(preds$fit), 3),
             SE        = round(as.numeric(preds$se.fit), 3),
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
      Response   = response_label,
      Factor     = factor_label,
      Timepoint  = paste0("Day ", target_day, " from 1 Sep 2025"),
      Difference = round(diff, 3),
      SE         = round(se_diff, 3),
      CI_lower   = round(lwr, 3),
      CI_upper   = round(upr, 3),
      Sig        = ifelse(sig, "YES ***", "no")
    ) |>
    select(Response, Factor, Comparison = comparison,
           Difference, SE, CI_lower, CI_upper, Sig)
}


# ── Shared prediction grid builder ───────────────────────────────────────────
make_full_grid <- function(days_seq, df_ref) {
  expand_grid(
    days = days_seq,
    Species   = levels(df_ref$Species),
    Spacing_f = levels(df_ref$Spacing_f)
  ) |>
    mutate(
      Species   = factor(Species,   levels = levels(df_ref$Species)),
      Spacing_f = factor(Spacing_f, levels = levels(df_ref$Spacing_f)),
      Culture   = factor("Single",  levels = levels(df_ref$Culture)),
      Plot_ID   = levels(df_ref$Plot_ID)[1],
      Tree_ID   = levels(df_ref$Tree_ID)[1]
    )
}

make_species_grid <- function(days_seq, df_ref) {
  expand_grid(days = days_seq, Species = levels(df_ref$Species)) |>
    mutate(
      Species   = factor(Species, levels = levels(df_ref$Species)),
      Spacing_f = levels(df_ref$Spacing_f)[1],
      Culture   = factor("Single", levels = levels(df_ref$Culture)),
      Plot_ID   = levels(df_ref$Plot_ID)[1],
      Tree_ID   = levels(df_ref$Tree_ID)[1]
    )
}

make_spacing_grid <- function(days_seq, df_ref) {
  expand_grid(days = days_seq, Spacing_f = levels(df_ref$Spacing_f)) |>
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
# If models already in memory, skip to PREDICTION GRIDS
# =============================================================================

cat("══ RESPONSE 1: Calibrated Height ══════════════════════════════════════\n")
models_h <- fit_gamm_pair(df_h, "Height", family = gaussian())
cat("\n── Height: Species model summary ────────────────────────────────────\n")
print(summary(models_h$species))
cat("\n── Height: Spacing model summary ────────────────────────────────────\n")
print(summary(models_h$spacing))

cat("\n══ RESPONSE 2: Crown Area ══════════════════════════════════════════════\n")
models_c <- fit_gamm_pair(df_c, "Crown", family = Gamma(link = "log"))
cat("\n── Crown: Species model summary ─────────────────────────────────────\n")
print(summary(models_c$species))
cat("\n── Crown: Spacing model summary ─────────────────────────────────────\n")
print(summary(models_c$spacing))

cat("\n══ RESPONSE 3: Crown:Height Ratio ══════════════════════════════════════\n")
cat("  CA:H = Crown_Area_m2 / Calibrated_Height_m  (m2 m-1)\n")
cat("  Gamma(log) family — strictly positive, right-skewed\n\n")
models_r <- fit_gamm_pair(df_r, "CAH", family = Gamma(link = "log"))
cat("\n── CA:H: Species model summary ──────────────────────────────────────\n")
print(summary(models_r$species))
cat("\n── CA:H: Spacing model summary ──────────────────────────────────────\n")
print(summary(models_r$spacing))


# =============================================================================
# ── PREDICTION GRIDS ──────────────────────────────────────────────────────────
# ▲▲▲ START HERE if all three model pairs already fitted ▲▲▲
# =============================================================================

days_h <- seq(min(df_h$days), max(df_h$days), length.out = 300)
days_c <- seq(min(df_c$days), max(df_c$days), length.out = 300)
days_r <- seq(min(df_r$days), max(df_r$days), length.out = 300)

cat("Predicting height trajectories...\n")
grid_h_sp_pred <- predict_traj(models_h$species, make_full_grid(days_h, df_h))
grid_h_sc_pred <- predict_traj(models_h$spacing, make_full_grid(days_h, df_h))

cat("Predicting crown area trajectories...\n")
grid_c_sp_pred <- predict_traj(models_c$species, make_full_grid(days_c, df_c))
grid_c_sc_pred <- predict_traj(models_c$spacing, make_full_grid(days_c, df_c))

cat("Predicting CA:H ratio trajectories...\n")
grid_r_sp_pred <- predict_traj(models_r$species, make_full_grid(days_r, df_r))
grid_r_sc_pred <- predict_traj(models_r$spacing, make_full_grid(days_r, df_r))

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

curve_c_sp <- predict_traj(models_c$species, make_species_grid(days_c, df_c)) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)
curve_c_sc <- predict_traj(models_c$spacing, make_spacing_grid(days_c, df_c)) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)

curve_r_sp <- predict_traj(models_r$species, make_species_grid(days_r, df_r)) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)
curve_r_sc <- predict_traj(models_r$spacing, make_spacing_grid(days_r, df_r)) |>
  mutate(lwr = fit - 1.96 * se, upr = fit + 1.96 * se)


# =============================================================================
# ── PLOTS ─────────────────────────────────────────────────────────────────────
# =============================================================================

cat("\nGenerating plots...\n")

# ── HEIGHT ────────────────────────────────────────────────────────────────────
p_h_sp_diff <- plot_diffs(sp_diffs_h,
                          title    = "Species pairwise height differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over spacings",
                          y_label  = "Difference in calibrated height (m)",
                          fill_col = "steelblue", ncol = 3)
ggsave("height_species_differences.pdf", p_h_sp_diff,
       width = 14, height = 9, units = "in", dpi = 300)

p_h_sc_diff <- plot_diffs(sc_diffs_h,
                          title    = "Spacing pairwise height differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over species",
                          y_label  = "Difference in calibrated height (m)",
                          fill_col = "darkorange", ncol = 3)
ggsave("height_spacing_differences.pdf", p_h_sc_diff,
       width = 10, height = 8, units = "in", dpi = 300)

p_h_sp_curves <- ggplot(curve_h_sp, aes(x = days, colour = Species, fill = Species)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 0.9) +
  scale_colour_manual(values = species_colors) +
  scale_fill_manual(values   = species_colors) +
  scale_x_continuous(breaks = seq(0, 270, by = 60)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title    = "GAMM-fitted height growth trajectories by species",
       subtitle = "Population mean, spacing controlled  |  Shaded = 95% CI",
       x = "Days from 1 September 2025",
       y = "Calibrated height (m)", colour = NULL, fill = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))
ggsave("height_curves_species.pdf", p_h_sp_curves,
       width = 7, height = 4.5, units = "in", dpi = 300)

p_h_sc_curves <- ggplot(curve_h_sc, aes(x = days, colour = Spacing_f, fill = Spacing_f)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 0.9) +
  scale_colour_manual(values = spacing_colors, labels = spacing_display) +
  scale_fill_manual(values   = spacing_colors, labels = spacing_display) +
  scale_x_continuous(breaks = seq(0, 270, by = 60)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title    = "GAMM-fitted height growth trajectories by spacing",
       subtitle = "Population mean, species controlled  |  Shaded = 95% CI",
       x = "Days from 1 September 2025",
       y = "Calibrated height (m)", colour = "Spacing", fill = "Spacing") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))
ggsave("height_curves_spacing.pdf", p_h_sc_curves,
       width = 7, height = 4.5, units = "in", dpi = 300)

# ── CROWN AREA ────────────────────────────────────────────────────────────────
p_c_sp_diff <- plot_diffs(sp_diffs_c,
                          title    = "Species pairwise crown area differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over spacings",
                          y_label  = "Difference in crown area (m2)",
                          fill_col = "steelblue", ncol = 3)
ggsave("crown_species_differences.pdf", p_c_sp_diff,
       width = 14, height = 9, units = "in", dpi = 300)

p_c_sc_diff <- plot_diffs(sc_diffs_c,
                          title    = "Spacing pairwise crown area differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over species",
                          y_label  = "Difference in crown area (m2)",
                          fill_col = "darkorange", ncol = 3)
ggsave("crown_spacing_differences.pdf", p_c_sc_diff,
       width = 10, height = 8, units = "in", dpi = 300)

p_c_sp_curves <- ggplot(curve_c_sp, aes(x = days, colour = Species, fill = Species)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 0.9) +
  scale_colour_manual(values = species_colors) +
  scale_fill_manual(values   = species_colors) +
  scale_x_continuous(breaks = seq(0, 270, by = 60)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title    = "GAMM-fitted crown area growth trajectories by species",
       subtitle = "Population mean, spacing controlled  |  Shaded = 95% CI",
       x = "Days from 1 September 2025",
       y = "Crown area (m2)", colour = NULL, fill = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))
ggsave("crown_curves_species.pdf", p_c_sp_curves,
       width = 7, height = 4.5, units = "in", dpi = 300)

p_c_sc_curves <- ggplot(curve_c_sc, aes(x = days, colour = Spacing_f, fill = Spacing_f)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 0.9) +
  scale_colour_manual(values = spacing_colors, labels = spacing_display) +
  scale_fill_manual(values   = spacing_colors, labels = spacing_display) +
  scale_x_continuous(breaks = seq(0, 270, by = 60)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title    = "GAMM-fitted crown area growth trajectories by spacing",
       subtitle = "Population mean, species controlled  |  Shaded = 95% CI",
       x = "Days from 1 September 2025",
       y = "Crown area (m2)", colour = "Spacing", fill = "Spacing") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))
ggsave("crown_curves_spacing.pdf", p_c_sc_curves,
       width = 7, height = 4.5, units = "in", dpi = 300)

# ── CA:H RATIO ────────────────────────────────────────────────────────────────
p_r_sp_diff <- plot_diffs(sp_diffs_r,
                          title    = "Species pairwise crown:height ratio differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over spacings",
                          y_label  = "Difference in CA:H ratio (m2 m-1)",
                          fill_col = "#6a3d9a", ncol = 3)
ggsave("cah_species_differences.pdf", p_r_sp_diff,
       width = 14, height = 9, units = "in", dpi = 300)

p_r_sc_diff <- plot_diffs(sc_diffs_r,
                          title    = "Spacing pairwise crown:height ratio differences",
                          subtitle = "Shaded = 95% CI  |  Red rug = significant period  |  Averaged over species",
                          y_label  = "Difference in CA:H ratio (m2 m-1)",
                          fill_col = "#e31a1c", ncol = 3)
ggsave("cah_spacing_differences.pdf", p_r_sc_diff,
       width = 10, height = 8, units = "in", dpi = 300)

p_r_sp_curves <- ggplot(curve_r_sp, aes(x = days, colour = Species, fill = Species)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 0.9) +
  scale_colour_manual(values = species_colors) +
  scale_fill_manual(values   = species_colors) +
  scale_x_continuous(breaks = seq(0, 270, by = 60)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title    = "GAMM-fitted crown:height ratio trajectories by species",
       subtitle = "Population mean, spacing controlled  |  Shaded = 95% CI",
       x = "Days from 1 September 2025",
       y = "CA:H ratio (m2 m-1)", colour = NULL, fill = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))
ggsave("cah_curves_species.pdf", p_r_sp_curves,
       width = 7, height = 4.5, units = "in", dpi = 300)

p_r_sc_curves <- ggplot(curve_r_sc, aes(x = days, colour = Spacing_f, fill = Spacing_f)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 0.9) +
  scale_colour_manual(values = spacing_colors, labels = spacing_display) +
  scale_fill_manual(values   = spacing_colors, labels = spacing_display) +
  scale_x_continuous(breaks = seq(0, 270, by = 60)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title    = "GAMM-fitted crown:height ratio trajectories by spacing",
       subtitle = "Population mean, species controlled  |  Shaded = 95% CI",
       x = "Days from 1 September 2025",
       y = "CA:H ratio (m2 m-1)", colour = "Spacing", fill = "Spacing") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))
ggsave("cah_curves_spacing.pdf", p_r_sc_curves,
       width = 7, height = 4.5, units = "in", dpi = 300)

# ── COMBINED PATCHWORK PANELS ─────────────────────────────────────────────────
# 3-panel: Height | Crown | CA:H ratio — by species
p_combined_species <- (p_h_sp_curves | p_c_sp_curves | p_r_sp_curves) +
  plot_annotation(
    title    = "GAMM-fitted growth trajectories by species",
    subtitle = "Left: Height  |  Centre: Crown area  |  Right: CA:H ratio  |  Spacing controlled",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, colour = "grey40")))
ggsave("combined_curves_species.pdf", p_combined_species,
       width = 21, height = 4.5, units = "in", dpi = 300)

# 3-panel: Height | Crown | CA:H ratio — by spacing
p_combined_spacing <- (p_h_sc_curves | p_c_sc_curves | p_r_sc_curves) +
  plot_annotation(
    title    = "GAMM-fitted growth trajectories by spacing",
    subtitle = "Left: Height  |  Centre: Crown area  |  Right: CA:H ratio  |  Species controlled",
    theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                     plot.subtitle = element_text(size = 9, colour = "grey40")))
ggsave("combined_curves_spacing.pdf", p_combined_spacing,
       width = 21, height = 4.5, units = "in", dpi = 300)


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

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TABLE 1: GAMM Smooth Term Significance                             ║\n")
cat("║  EDF > 1 = non-linear  |  p < 0.05 = significant smooth            ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
print(as_tibble(tbl1), n = Inf)
write_csv(tbl1, "stats_table1_smooth_significance.csv")

# TABLE 2: Marginal means at key timepoints
tbl2_h_sp <- marginal_means(models_h$species, "Species",   levels(df_h$Species),
                            "Spacing_f", levels(df_h$Spacing_f)[1],
                            df_h, key_days, key_labels, "Height (m)")
tbl2_h_sc <- marginal_means(models_h$spacing, "Spacing_f", levels(df_h$Spacing_f),
                            "Species",   levels(df_h$Species)[1],
                            df_h, key_days, key_labels, "Height (m)")
tbl2_c_sp <- marginal_means(models_c$species, "Species",   levels(df_c$Species),
                            "Spacing_f", levels(df_c$Spacing_f)[1],
                            df_c, key_days, key_labels, "Crown Area (m2)")
tbl2_c_sc <- marginal_means(models_c$spacing, "Spacing_f", levels(df_c$Spacing_f),
                            "Species",   levels(df_c$Species)[1],
                            df_c, key_days, key_labels, "Crown Area (m2)")
tbl2_r_sp <- marginal_means(models_r$species, "Species",   levels(df_r$Species),
                            "Spacing_f", levels(df_r$Spacing_f)[1],
                            df_r, key_days, key_labels, "CA:H Ratio (m2 m-1)")
tbl2_r_sc <- marginal_means(models_r$spacing, "Spacing_f", levels(df_r$Spacing_f),
                            "Species",   levels(df_r$Species)[1],
                            df_r, key_days, key_labels, "CA:H Ratio (m2 m-1)")

for (tbl in list(list(tbl2_h_sp, "2A: Heights by Species"),
                 list(tbl2_h_sc, "2B: Heights by Spacing"),
                 list(tbl2_c_sp, "2C: Crown Area by Species"),
                 list(tbl2_c_sc, "2D: Crown Area by Spacing"),
                 list(tbl2_r_sp, "2E: CA:H Ratio by Species"),
                 list(tbl2_r_sc, "2F: CA:H Ratio by Spacing"))) {
  cat(paste0("\n╔══ TABLE ", tbl[[2]], " ══╗\n"))
  print(tbl[[1]] |> pivot_wider(names_from = Timepoint,
                                values_from = c(Mean, SE), names_glue = "{Timepoint} {.value}"), n = Inf)
}

write_csv(bind_rows(tbl2_h_sp, tbl2_h_sc, tbl2_c_sp,
                    tbl2_c_sc, tbl2_r_sp, tbl2_r_sc),
          "stats_table2_marginal_means.csv")

# TABLE 3: Pairwise differences at day 266
tbl3 <- bind_rows(
  pairwise_at_day(sp_diffs_h, final_day, "Height (m)",         "Species"),
  pairwise_at_day(sc_diffs_h, final_day, "Height (m)",         "Spacing"),
  pairwise_at_day(sp_diffs_c, final_day, "Crown Area (m2)",    "Species"),
  pairwise_at_day(sc_diffs_c, final_day, "Crown Area (m2)",    "Spacing"),
  pairwise_at_day(sp_diffs_r, final_day, "CA:H Ratio (m2 m-1)", "Species"),
  pairwise_at_day(sc_diffs_r, final_day, "CA:H Ratio (m2 m-1)", "Spacing")
)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TABLE 3: All Pairwise Differences at Day 266 (25 May 2026)        ║\n")
cat("║  Difference = Group1 minus Group2 | Sig = CI excludes zero         ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
print(tbl3, n = Inf)
write_csv(tbl3, "stats_table3_pairwise_differences.csv")

cat("\n── Significant pairs at Day 266 ─────────────────────────────────────\n")
sig_only <- tbl3 |> filter(Sig == "YES ***")
for (resp in unique(tbl3$Response)) {
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
for (nm in c("Height species", "Height spacing",
             "Crown species",  "Crown spacing",
             "CA:H species",   "CA:H spacing")) {
  model <- switch(nm,
                  "Height species" = models_h$species, "Height spacing" = models_h$spacing,
                  "Crown species"  = models_c$species, "Crown spacing"  = models_c$spacing,
                  "CA:H species"   = models_r$species, "CA:H spacing"   = models_r$spacing)
  cat(nm, "model:\n"); print(k.check(model))
}

ggsave("diag_height_species.pdf", appraise(models_h$species),
       width = 10, height = 8, units = "in", dpi = 300)
ggsave("diag_height_spacing.pdf", appraise(models_h$spacing),
       width = 10, height = 8, units = "in", dpi = 300)
ggsave("diag_crown_species.pdf",  appraise(models_c$species),
       width = 10, height = 8, units = "in", dpi = 300)
ggsave("diag_crown_spacing.pdf",  appraise(models_c$spacing),
       width = 10, height = 8, units = "in", dpi = 300)
ggsave("diag_cah_species.pdf",    appraise(models_r$species),
       width = 10, height = 8, units = "in", dpi = 300)
ggsave("diag_cah_spacing.pdf",    appraise(models_r$spacing),
       width = 10, height = 8, units = "in", dpi = 300)


cat("\n── Saved outputs ─────────────────────────────────────────────────────\n")
cat("  STATISTICS TABLES\n")
cat("    stats_table1_smooth_significance.csv\n")
cat("    stats_table2_marginal_means.csv\n")
cat("    stats_table3_pairwise_differences.csv\n")
cat("  HEIGHT\n")
cat("    height_curves_species.pdf  |  height_curves_spacing.pdf\n")
cat("    height_species_differences.pdf  |  height_spacing_differences.pdf\n")
cat("  CROWN AREA\n")
cat("    crown_curves_species.pdf  |  crown_curves_spacing.pdf\n")
cat("    crown_species_differences.pdf  |  crown_spacing_differences.pdf\n")
cat("  CA:H RATIO\n")
cat("    cah_curves_species.pdf  |  cah_curves_spacing.pdf\n")
cat("    cah_species_differences.pdf  |  cah_spacing_differences.pdf\n")
cat("  COMBINED (3-panel)\n")
cat("    combined_curves_species.pdf  |  combined_curves_spacing.pdf\n")
cat("  DIAGNOSTICS\n")
cat("    diag_height/crown/cah _species/_spacing .pdf\n")
cat("\n── Done ──────────────────────────────────────────────────────────────\n")