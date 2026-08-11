# ──────────────────────────────────────────────────────────────────────────────
# LOCAL MORAN'S I (LISA) — EucVision GAMM Residuals
#
# Extends the global Moran's I analysis (10__EucVision_Moran_s_I.R) with a
# LOCAL indicator of spatial association (LISA) per tree, per date, per
# response. Global Moran's I tells you WHETHER spatial autocorrelation
# exists; local Moran's I tells you WHERE it concentrates — specifically,
# whether significant clustering aligns with plot boundaries (supporting a
# design-driven microsite/edge-effect interpretation) or crosses them freely
# (which would argue against a purely within-plot explanation).
#
# CLUSTER CLASSIFICATION (standard Anselin LISA quadrant convention):
#   High-High : residual above the sample mean, spatially lagged neighbours
#               also above the mean (positive cluster of over-predictions)
#   Low-Low   : residual below the mean, neighbours also below the mean
#               (positive cluster of under-predictions)
#   High-Low / Low-High: spatial outliers (value and neighbourhood disagree)
#   Not significant: local p >= 0.05, no reportable local pattern
#
# DEPENDENCIES (run after 10__EucVision_Moran_s_I.R, or ensure these exist
# in the current session):
#   - df_h$resid, df_c$resid, df_r$resid   (conditional residuals)
#   - flight_registry                       (26-date shapefile registry)
#   - T0                                    (as.Date("2025-09-01"))
#   - theme_thesis()                        (thesis ggplot theme)
#
# OUTPUT:
#   - Console summary (top clustering plots per response)
#   - Excel workbook: per-tree results, per-plot summary, notes
#   - Figure 1: lisa_temporal_trend.png — % significant trees over time
#   - Figure 2: lisa_spatial_snapshots.png — cluster maps at 3 representative
#     dates (first / middle / last flight), faceted by response
# ──────────────────────────────────────────────────────────────────────────────

library(sf)
library(spdep)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(openxlsx)


# ──────────────────────────────────────────────────────────────────────────────
# 1. LOCAL MORAN'S I FUNCTION
# ──────────────────────────────────────────────────────────────────────────────

run_local_moran <- function(data, target_date, shp_path, resid_col, response_name) {
  
  # Restrict crowns to ONLY the join-key columns (+ geometry, kept
  # automatically by sf::select) before joining. Crown shapefiles typically
  # carry their own "Species" (and sometimes "Spacing") attribute; if left
  # in, a name collision with df_sub's columns causes dplyr to silently
  # rename both to Species.x/Species.y, joined$Species then silently
  # returns NULL, and tibble() silently drops it — no error until much
  # later when something tries to use the missing column. Restricting to
  # join keys here removes the possibility of that collision entirely.
  crowns <- st_read(shp_path, quiet = TRUE) %>%
    mutate(Tree = round(as.numeric(Tree), 2)) %>%
    select(Cmprtmn, Line, Plot, Tree)
  
  df_sub <- data %>%
    filter(Date == as.Date(target_date)) %>%
    mutate(Tree = round(as.numeric(Tree), 2))
  
  df_sub$target_resid <- df_sub[[resid_col]]
  
  joined <- crowns %>%
    left_join(
      df_sub %>% select(Compartment, Line, Plot, Tree, target_resid,
                        Plot_ID, Species, Spacing_f),
      by = c("Cmprtmn" = "Compartment", "Line", "Plot", "Tree")
    ) %>%
    filter(!is.na(target_resid))
  
  # Need at least k+1 observations to build weights matrix
  if (nrow(joined) < 7) {
    warning(paste("Skipping", target_date, response_name, "- too few observations"))
    return(NULL)
  }
  
  coords <- st_coordinates(st_centroid(joined))
  nb     <- knn2nb(knearneigh(coords, k = 6))
  lw     <- nb2listw(nb, style = "W")
  
  lisa <- localmoran(joined$target_resid, lw)
  
  # Column names vary slightly across spdep versions — match by pattern
  # rather than hardcoding, so this doesn't silently break on an update.
  p_col <- grep("^Pr\\(", colnames(lisa))[1]
  z_col <- grep("^Z",     colnames(lisa))[1]
  
  local_p <- as.numeric(lisa[, p_col])
  local_z <- as.numeric(lisa[, z_col])
  local_I <- as.numeric(lisa[, "Ii"])
  local_sig <- local_p < 0.05
  
  # Standard Anselin LISA quadrant classification: mean-centre the value,
  # compute its spatial lag, and classify by the sign pair — this is more
  # reliable than using the sign of Ii alone.
  resid_c <- joined$target_resid - mean(joined$target_resid)
  lag_c   <- lag.listw(lw, resid_c)
  
  cluster_type <- case_when(
    !local_sig                ~ "Not significant",
    resid_c > 0 & lag_c > 0   ~ "High-High (positive cluster)",
    resid_c < 0 & lag_c < 0   ~ "Low-Low (positive cluster)",
    resid_c > 0 & lag_c < 0   ~ "High-Low (spatial outlier)",
    resid_c < 0 & lag_c > 0   ~ "Low-High (spatial outlier)",
    TRUE                      ~ "Not significant"
  )
  
  tibble(
    Date         = as.Date(target_date),
    days         = as.numeric(as.Date(target_date) - T0),
    Response     = response_name,
    Plot_ID      = joined$Plot_ID,
    Species      = joined$Species,
    Spacing_f    = joined$Spacing_f,
    Tree         = joined$Tree,
    X            = coords[, 1],
    Y            = coords[, 2],
    resid        = joined$target_resid,
    local_I      = local_I,
    local_z      = local_z,
    local_p      = local_p,
    local_sig    = local_sig,
    cluster_type = cluster_type
  )
}


# ──────────────────────────────────────────────────────────────────────────────
# 2. RUN ALL COMBINATIONS (26 dates x 3 responses)
# ──────────────────────────────────────────────────────────────────────────────

model_specs <- list(
  list(data = df_h, resid_col = "resid", response = "Height"),
  list(data = df_c, resid_col = "resid", response = "Crown Area"),
  list(data = df_r, resid_col = "resid", response = "CA:H Ratio")
)

total_tests <- length(flight_registry) * length(model_specs)
cat(paste0("Running ", total_tests, " local Moran's I (LISA) passes across ",
           length(flight_registry), " dates x ", length(model_specs), " responses...\n"))
cat("This will take several minutes (each pass classifies every tree). Progress below:\n\n")

lisa_list <- list()
counter    <- 0

for (flight in flight_registry) {
  for (spec in model_specs) {
    counter <- counter + 1
    cat(sprintf("  [%3d/%d] %s - %s\n",
                counter, total_tests,
                flight$date, spec$response))
    
    tab <- tryCatch(
      run_local_moran(spec$data, flight$date, flight$shp,
                      spec$resid_col, spec$response),
      error = function(e) {
        warning(paste("Error at", flight$date, spec$response, ":", e$message))
        NULL
      }
    )
    if (!is.null(tab)) lisa_list[[length(lisa_list) + 1]] <- tab
  }
}

lisa_results <- bind_rows(lisa_list)

lisa_results <- lisa_results %>%
  mutate(
    Response = factor(Response,
                      levels = c("Height", "Crown Area", "CA:H Ratio"),
                      labels = c("Calibrated height (m)",
                                 "Crown area (m\u00b2)",
                                 "CA:H ratio (m\u00b2 m\u207b\u00b9)")),
    cluster_type = factor(cluster_type,
                          levels = c("High-High (positive cluster)",
                                     "Low-Low (positive cluster)",
                                     "High-Low (spatial outlier)",
                                     "Low-High (spatial outlier)",
                                     "Not significant"))
  )

cat("\nTotal tree-date observations classified:", nrow(lisa_results), "\n")


# ──────────────────────────────────────────────────────────────────────────────
# 3. CONSOLE SUMMARY — which plots show persistent local clustering?
# ──────────────────────────────────────────────────────────────────────────────

plot_summary <- lisa_results %>%
  group_by(Response, Plot_ID) %>%
  summarise(
    Species    = first(Species),
    Spacing_f  = first(Spacing_f),
    N_obs      = n(),
    Pct_Sig    = round(100 * mean(local_sig), 1),
    Pct_HH     = round(100 * mean(cluster_type == "High-High (positive cluster)"), 1),
    Pct_LL     = round(100 * mean(cluster_type == "Low-Low (positive cluster)"), 1),
    Pct_Outlier = round(100 * mean(cluster_type %in% c("High-Low (spatial outlier)",
                                                       "Low-High (spatial outlier)")), 1),
    .groups = "drop"
  ) %>%
  arrange(Response, desc(Pct_Sig))

cat("\n================ TOP 5 PLOTS BY % SIGNIFICANT LOCAL CLUSTERING ================\n")
for (resp in levels(lisa_results$Response)) {
  cat(paste0("\n-- ", resp, " --\n"))
  print(
    as.data.frame(
      plot_summary %>% filter(Response == resp) %>% slice_head(n = 5) %>%
        select(Plot_ID, Species, Spacing_f, N_obs, Pct_Sig, Pct_HH, Pct_LL, Pct_Outlier)
    ),
    row.names = FALSE
  )
}
cat("\n=================================================================================\n")


# ──────────────────────────────────────────────────────────────────────────────
# 4. TEMPORAL TREND: % OF TREES SIGNIFICANT PER DATE, PER RESPONSE
# ──────────────────────────────────────────────────────────────────────────────

lisa_temporal <- lisa_results %>%
  group_by(Response, Date, days) %>%
  summarise(Pct_Sig = 100 * mean(local_sig), .groups = "drop")

make_lisa_panel <- function(data, response_label, keep_x = FALSE) {
  
  df_panel <- data %>% filter(Response == response_label)
  
  p <- ggplot(df_panel, aes(x = days, y = Pct_Sig)) +
    geom_line(colour = "#8E44AD", linewidth = 0.8) +
    geom_point(colour = "#8E44AD", size = 1.8, shape = 16) +
    scale_x_continuous(breaks = seq(0, 270, by = 60), limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.08))) +
    labs(y = paste0("% trees in significant\nlocal cluster  -  ", as.character(response_label))) +
    theme_thesis()
  
  if (!keep_x) {
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  } else {
    p <- p + labs(x = "Days from 1 September 2025")
  }
  p
}

resp_levels <- levels(lisa_results$Response)

pl_h <- make_lisa_panel(lisa_temporal, resp_levels[1], keep_x = FALSE)
pl_c <- make_lisa_panel(lisa_temporal, resp_levels[2], keep_x = FALSE)
pl_r <- make_lisa_panel(lisa_temporal, resp_levels[3], keep_x = TRUE)

p_lisa_temporal <- pl_h / pl_c / pl_r

ggsave(file.path(OUTPUT_DIR, "lisa_temporal_trend.png"), p_lisa_temporal,
       width = 8, height = 7, units = "in", dpi = 300)

cat("Saved: lisa_temporal_trend.png\n")


# ──────────────────────────────────────────────────────────────────────────────
# 5. SPATIAL SNAPSHOTS — cluster maps at 3 representative dates
# ──────────────────────────────────────────────────────────────────────────────

snapshot_dates <- as.Date(c(
  flight_registry[[1]]$date,
  flight_registry[[ceiling(length(flight_registry) / 2)]]$date,
  flight_registry[[length(flight_registry)]]$date
))

lisa_snap <- lisa_results %>%
  filter(Date %in% snapshot_dates) %>%
  mutate(Date_label = factor(format(Date, "%d %b %Y"),
                             levels = format(snapshot_dates, "%d %b %Y")))

lisa_colors <- c(
  "High-High (positive cluster)" = "#C0392B",
  "Low-Low (positive cluster)"   = "#2980B9",
  "High-Low (spatial outlier)"   = "#E67E22",
  "Low-High (spatial outlier)"   = "#8E44AD",
  "Not significant"              = "grey85"
)

p_lisa_snap <- ggplot(lisa_snap, aes(x = X, y = Y, colour = cluster_type)) +
  geom_point(size = 0.8, alpha = 0.85) +
  scale_colour_manual(values = lisa_colors, name = "LISA cluster") +
  coord_equal() +
  facet_grid(Response ~ Date_label) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.5))) +
  labs(x = NULL, y = NULL) +
  theme_thesis() +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 8)
  )

ggsave(file.path(OUTPUT_DIR, "lisa_spatial_snapshots.png"), p_lisa_snap,
       width = 9, height = 7, units = "in", dpi = 300)

cat("Saved: lisa_spatial_snapshots.png\n")


# ──────────────────────────────────────────────────────────────────────────────
# 6. EXCEL EXPORT
# ──────────────────────────────────────────────────────────────────────────────

output_xlsx <- file.path(OUTPUT_DIR, "EucVision_LISA_Results.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "LISA_Raw_Results")
addWorksheet(wb, "Plot_Summary")
addWorksheet(wb, "Notes")

style_title <- createStyle(
  fontName = "Arial", fontSize = 13, fontColour = "#2E4057",
  textDecoration = "bold", halign = "left", valign = "center"
)
style_header <- createStyle(
  fontName = "Arial", fontSize = 11, fontColour = "white",
  fgFill = "#2E4057", halign = "center", valign = "center",
  textDecoration = "bold", wrapText = TRUE,
  border = "Bottom", borderColour = "#FFFFFF", borderStyle = "medium"
)
style_center <- createStyle(
  fontName = "Arial", fontSize = 10,
  halign = "center", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC", borderStyle = "thin"
)

# Sheet 1: raw per-tree, per-date results
export_raw <- lisa_results %>%
  select(Date, days, Response, Plot_ID, Species, Spacing_f,
         X, Y, resid, local_I, local_z, local_p, local_sig, cluster_type) %>%
  arrange(Response, Date, Plot_ID)

writeData(wb, "LISA_Raw_Results",
          "EucVision GAMM - Local Moran's I (LISA) Cluster Classification, All Trees x All Dates",
          startRow = 1, startCol = 1)
addStyle(wb, "LISA_Raw_Results", style_title, rows = 1, cols = 1)
mergeCells(wb, "LISA_Raw_Results", cols = 1:14, rows = 1)
writeData(wb, "LISA_Raw_Results", export_raw,
          startRow = 2, startCol = 1, headerStyle = style_header, withFilter = TRUE)
setColWidths(wb, "LISA_Raw_Results", cols = 1:14,
             widths = c(13, 8, 24, 12, 16, 10, 12, 12, 10, 10, 10, 10, 10, 30))
freezePane(wb, "LISA_Raw_Results", firstActiveRow = 3)

# Sheet 2: per-plot summary
writeData(wb, "Plot_Summary",
          "Per-Plot Clustering Summary (aggregated across all flight dates)",
          startRow = 1, startCol = 1)
addStyle(wb, "Plot_Summary", style_title, rows = 1, cols = 1)
mergeCells(wb, "Plot_Summary", cols = 1:8, rows = 1)
writeData(wb, "Plot_Summary", plot_summary,
          startRow = 2, startCol = 1, headerStyle = style_header, withFilter = TRUE)
setColWidths(wb, "Plot_Summary", cols = 1:9,
             widths = c(24, 12, 16, 12, 10, 10, 10, 10, 10))
freezePane(wb, "Plot_Summary", firstActiveRow = 3)

# Sheet 3: notes
notes_df <- data.frame(
  Item = c(
    "Script version",
    "Purpose",
    "Cluster classification",
    "Spatial weights",
    "Significance threshold",
    "How to read Plot_Summary",
    "Key reference"
  ),
  Detail = c(
    paste("Generated:", Sys.time()),
    "Local Moran's I (LISA) per tree, per date, per response — identifies WHERE residual spatial autocorrelation concentrates, extending the global Moran's I test (10__EucVision_Moran_s_I.R) which only established THAT it exists.",
    "Standard Anselin LISA quadrant convention: High-High and Low-Low are positive local clusters (residual and its spatial lag agree in sign, mean-centred); High-Low and Low-High are spatial outliers (value and neighbourhood disagree). Not significant: local p >= 0.05.",
    "k-nearest neighbours, k = 6, row-standardised (style='W'), matching the global Moran's I test for consistency.",
    "p < 0.05 (two-sided), from spdep::localmoran().",
    "Pct_Sig = % of tree-date observations in that plot classified as any significant cluster type (HH/LL/HL/LH combined). High and persistent Pct_Sig across many dates for a plot supports a design-driven (microsite, edge, or grazing) explanation rather than a scattered/random pattern.",
    "Anselin, L. (1995) Local Indicators of Spatial Association-LISA. Geographical Analysis 27(2):93-115."
  )
)

writeData(wb, "Notes", notes_df, startRow = 1, startCol = 1,
          headerStyle = createStyle(fontName = "Arial", fontSize = 11,
                                    textDecoration = "bold",
                                    fgFill = "#2E4057", fontColour = "white"))
setColWidths(wb, "Notes", cols = 1:2, widths = c(28, 100))

saveWorkbook(wb, output_xlsx, overwrite = TRUE)
cat(paste0("Excel saved: ", output_xlsx, "\n"))


# ──────────────────────────────────────────────────────────────────────────────
# 7. SAVED OUTPUT SUMMARY
# ──────────────────────────────────────────────────────────────────────────────

cat("\n-- Saved outputs ----------------------------------------------------------\n")
cat("  FIGURES\n")
cat("    lisa_temporal_trend.png     (3-panel, 8 x 7 in, 300 dpi)\n")
cat("    lisa_spatial_snapshots.png  (3-response x 3-date grid, 9 x 7 in, 300 dpi)\n")
cat("  EXCEL\n")
cat("    EucVision_LISA_Results.xlsx\n")
cat("      Sheet 1: Raw per-tree, per-date results (", nrow(lisa_results), "rows)\n")
cat("      Sheet 2: Per-plot summary (", nrow(plot_summary), "rows)\n")
cat("      Sheet 3: Methodology notes\n")
cat("-- Done -------------------------------------------------------------------\n")

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(viridis)

# 1. Classify each tree's LISA behaviour across the season
#    (dominance-ratio version — a tree needs an 80/20+ split toward one
#    direction to be called "Dominant"; anything more balanced than that
#    is a genuine see-saw, not classification noise from a single blip)
persistence_df <- lisa_results %>%
  group_by(Plot_ID, Species, Spacing_f, Tree, Response) %>%
  summarise(
    total_flights = n(),
    n_HH      = sum(cluster_type == "High-High (positive cluster)", na.rm = TRUE),
    n_LL      = sum(cluster_type == "Low-Low (positive cluster)", na.rm = TRUE),
    n_outlier = sum(cluster_type %in% c("High-Low (spatial outlier)",
                                        "Low-High (spatial outlier)"), na.rm = TRUE),
    sig_count = n_HH + n_LL + n_outlier,
    persistence_pct = (sig_count / total_flights) * 100,
    # -1 = always Low-Low, +1 = always High-High, 0 = perfectly balanced
    net_bias  = if_else(n_HH + n_LL > 0, (n_HH - n_LL) / (n_HH + n_LL), NA_real_),
    # Static coords, but averaged defensively in case of any sub-pixel jitter
    X = mean(X),
    Y = mean(Y),
    .groups = "drop"
  ) %>%
  mutate(
    behaviour = case_when(
      sig_count == 0             ~ "Not significant",
      n_HH == 0 & n_LL == 0      ~ "Spatial outlier only",
      abs(net_bias) >= 0.6       ~ if_else(net_bias > 0, "Dominant High-High", "Dominant Low-Low"),
      TRUE                       ~ "True see-saw (balanced HH & LL)"
    ),
    behaviour = factor(behaviour, levels = c(
      "True see-saw (balanced HH & LL)",
      "Dominant High-High",
      "Dominant Low-Low",
      "Spatial outlier only",
      "Not significant"
    ))
  )

cat("\n-- Refined behaviour counts by response --\n")
print(persistence_df %>% count(Response, behaviour))

# 2. Generate the Master Behaviour Map
behaviour_colors <- c(
  "True see-saw (balanced HH & LL)" = "#7B2D8E",  # flags real model instability at that tree
  "Dominant High-High"              = "#C0392B",  # persistently under-predicted
  "Dominant Low-Low"                = "#2980B9",  # persistently over-predicted
  "Spatial outlier only"            = "#E67E22",
  "Not significant"                 = "grey88"
)

persistence_map <- ggplot(persistence_df,
                          aes(x = X, y = Y, colour = behaviour, alpha = persistence_pct)) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = behaviour_colors, name = "LISA behaviour\nacross season") +
  scale_alpha_continuous(range = c(0.15, 1), name = "Cluster\npersistence (%)", limits = c(0, 100)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  facet_wrap(~ Response, ncol = 3) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "EucVision: LISA Cluster Behaviour Map (2025-2026)",
    subtitle = "Colour = dominant cluster type across all flights | Opacity = how consistently (persistence %)",
    x = "Easting (Lo19)",
    y = "Northing (Lo19)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.background = element_rect(fill = "white", color = "grey80"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# 3. Save the output
ggsave(file.path(OUTPUT_DIR, "lisa_behaviour_map.png"), persistence_map,
       width = 16, height = 6, dpi = 300)

cat("Behaviour map saved as lisa_behaviour_map.png\n")

# ──────────────────────────────────────────────────────────────────────────────
# 8. TEMPORAL CHECK — do see-saw trees flip once, or oscillate randomly?
# ──────────────────────────────────────────────────────────────────────────────
# A per-tree spaghetti plot isn't viable here (700-850 see-saw trees per
# response would be unreadable). Instead:
#   (a) Population-level line chart: % of see-saw trees that are HH vs LL,
#       per date. A clean single flip (e.g. drought-driven) shows one line
#       dropping as the other rises around a specific date. Noisy
#       oscillation shows both jittering with no clear crossover.
#   (b) Tile/heatmap for the top plots by see-saw tree count: one row per
#       tree, one column per date, coloured by that date's cluster_type.
# ──────────────────────────────────────────────────────────────────────────────

# Identify which trees are "True see-saw" per response
seesaw_trees <- persistence_df %>%
  filter(behaviour == "True see-saw (balanced HH & LL)") %>%
  select(Response, Plot_ID, Tree)

# Pull their full per-date time series from lisa_results
seesaw_timeseries <- lisa_results %>%
  semi_join(seesaw_trees, by = c("Response", "Plot_ID", "Tree"))

# (a) Population-level HH vs LL balance over time, see-saw trees only
seesaw_temporal <- seesaw_timeseries %>%
  group_by(Response, Date, days) %>%
  summarise(
    pct_HH = 100 * mean(cluster_type == "High-High (positive cluster)"),
    pct_LL = 100 * mean(cluster_type == "Low-Low (positive cluster)"),
    .groups = "drop"
  ) %>%
  pivot_longer(c(pct_HH, pct_LL), names_to = "type", values_to = "pct") %>%
  mutate(type = recode(type, pct_HH = "High-High", pct_LL = "Low-Low"))

make_seesaw_panel <- function(data, response_label, keep_x = FALSE) {
  df_panel <- data %>% filter(Response == response_label)
  p <- ggplot(df_panel, aes(x = days, y = pct, colour = type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = c("High-High" = "#C0392B", "Low-Low" = "#2980B9"),
                        name = NULL) +
    scale_x_continuous(breaks = seq(0, 270, by = 60), limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.08))) +
    labs(y = paste0("% of see-saw trees\n-  ", as.character(response_label))) +
    theme_thesis()
  if (!keep_x) {
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  } else {
    p <- p + labs(x = "Days from 1 September 2025")
  }
  if (response_label != levels(data$Response)[1]) p <- p + theme(legend.position = "none")
  p
}

resp_levels <- levels(lisa_results$Response)
ss_h <- make_seesaw_panel(seesaw_temporal, resp_levels[1], keep_x = FALSE)
ss_c <- make_seesaw_panel(seesaw_temporal, resp_levels[2], keep_x = FALSE)
ss_r <- make_seesaw_panel(seesaw_temporal, resp_levels[3], keep_x = TRUE)

p_seesaw_temporal <- ss_h / ss_c / ss_r
ggsave(file.path(OUTPUT_DIR, "lisa_seesaw_temporal.png"), p_seesaw_temporal,
       width = 8, height = 7, units = "in", dpi = 300)
cat("Saved: lisa_seesaw_temporal.png\n")

# (b) Per-tree tile/heatmap for the top plots by see-saw tree count
top_seesaw_plots <- seesaw_trees %>%
  count(Response, Plot_ID, name = "n_seesaw_trees") %>%
  group_by(Response) %>%
  slice_max(n_seesaw_trees, n = 3) %>%
  ungroup()

heatmap_data <- seesaw_timeseries %>%
  semi_join(top_seesaw_plots, by = c("Response", "Plot_ID")) %>%
  mutate(Tree_label = paste0(Plot_ID, "_T", Tree))

# Order trees within each plot by Y position so spatial neighbours sit
# next to each other on the y-axis, making spatial sync easier to spot
tree_order <- heatmap_data %>%
  distinct(Response, Plot_ID, Tree_label, Y) %>%
  arrange(Response, Plot_ID, Y) %>%
  mutate(row_order = row_number())

heatmap_data <- heatmap_data %>%
  left_join(tree_order %>% select(Response, Tree_label, row_order),
            by = c("Response", "Tree_label"))

p_seesaw_heatmap <- ggplot(heatmap_data,
                           aes(x = days, y = row_order, fill = cluster_type)) +
  geom_tile() +
  scale_fill_manual(values = lisa_colors, name = "LISA cluster") +
  facet_grid(Plot_ID ~ Response, scales = "free_y", space = "free_y") +
  labs(x = "Days from 1 September 2025", y = "Tree (ordered by position within plot)") +
  theme_thesis() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8),
    strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.5)
  )

ggsave(file.path(OUTPUT_DIR, "lisa_seesaw_heatmap.png"), p_seesaw_heatmap,
       width = 10, height = 8, units = "in", dpi = 300)
cat("Saved: lisa_seesaw_heatmap.png\n")