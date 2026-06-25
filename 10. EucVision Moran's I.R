# =============================================================================
# MORAN'S I SPATIAL AUTOCORRELATION — EucVision GAMM Residuals
#
# Tests residual spatial autocorrelation across all 26 UAV flight dates
# for all 6 GAMM models (3 responses × 2 grouping factors).
#
# RESIDUAL TYPE: Conditional residuals (observed - full fitted values).
# predict(..., type = "response") includes all smooths AND random effects
# (s(Plot_ID), s(Tree_ID)), so spatial clustering driven by experimental
# design is properly accounted for before the test runs.
#
# Crown and CA:H fitted on log scale (Gaussian identity link), so
# predict() returns log-scale values — exp() applied before differencing
# against raw observed values.
#
# OUTPUT:
#   - Console table
#   - Excel workbook with results + notes sheet
#   - 3-panel line chart (Moran's I over time by response variable)
# =============================================================================

library(sf)
library(spdep)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(openxlsx)


# =============================================================================
# SETTINGS
# =============================================================================

OUTPUT_DIR <- "C:/Users/jakev/Stellenbosch University/JacquesV B.Sc. skripsie M.Sc. project - Documents/Processed Data/EucVision/10. GAMM"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Add Calibri font (matches thesis template)
windowsFonts(Calibri = windowsFont("Calibri"))

# Analysis start date (matches GAMM t0)
T0 <- as.Date("2025-09-01")


# =============================================================================
# 1. ATTACH CONDITIONAL RESIDUALS
# observed minus full model fitted values (all terms including random effects)
# =============================================================================

# Height — raw scale, column renamed to "Height" in df_h
df_h$resid_sp <- df_h$Height - predict(models_h$species, type = "response")
df_h$resid_sc <- df_h$Height - predict(models_h$spacing,  type = "response")

# Crown Area — modelled on log(Crown_Area_m2); exp() to back-transform
df_c$resid_sp <- df_c$Crown_Area_m2 - exp(predict(models_c$species, type = "response"))
df_c$resid_sc <- df_c$Crown_Area_m2 - exp(predict(models_c$spacing,  type = "response"))

# CA:H Ratio — modelled on log(Crown_Area_m2 / Calibrated_Height_m); exp() to back-transform
df_r$resid_sp <- (df_r$Crown_Area_m2 / df_r$Calibrated_Height_m) - exp(predict(models_r$species, type = "response"))
df_r$resid_sc <- (df_r$Crown_Area_m2 / df_r$Calibrated_Height_m) - exp(predict(models_r$spacing,  type = "response"))


# =============================================================================
# 2. FLIGHT DATE — SHAPEFILE REGISTRY (all 26 dates)
# =============================================================================

flight_registry <- list(
  list(date = "2025-09-01", shp = "E:/Remote Sensing Media/02. 01 September 2025/09. Crown Metrics/Crown_Metrics_01_September_2025.shp"),
  list(date = "2025-10-30", shp = "E:/Remote Sensing Media/03. 30 October 2025/09. Crown Metrics/Crown_Metrics_30_October_2025.shp"),
  list(date = "2025-11-07", shp = "E:/Remote Sensing Media/04. 07 November 2025/09. Crown Metrics/Crown_Metrics_07_November_2025.shp"),
  list(date = "2025-11-14", shp = "E:/Remote Sensing Media/05. 14 November 2025/09. Crown Metrics/Crown_Metrics_14_November_2025.shp"),
  list(date = "2025-11-17", shp = "E:/Remote Sensing Media/06. 17 November 2025/09. Crown Metrics/Crown_Metrics_17_November_2025.shp"),
  list(date = "2025-11-28", shp = "E:/Remote Sensing Media/07. 28 November 2025/09. Crown Metrics/Crown_Metrics_28_November_2025.shp"),
  list(date = "2025-12-08", shp = "E:/Remote Sensing Media/08. 08 December 2025/09. Crown Metrics/Crown_Metrics_08_December_2025.shp"),
  list(date = "2025-12-11", shp = "E:/Remote Sensing Media/09. 11 December 2025/09. Crown Metrics/Crown_Metrics_11_December_2025.shp"),
  list(date = "2025-12-22", shp = "E:/Remote Sensing Media/10. 22 December 2025/09. Crown Metrics/Crown_Metrics_22_December_2025.shp"),
  list(date = "2026-01-13", shp = "E:/Remote Sensing Media/11. 13 January 2026/09. Crown Metrics/Crown_Metrics_13_January_2026.shp"),
  list(date = "2026-01-22", shp = "E:/Remote Sensing Media/12. 22 January 2026/09. Crown Metrics/Crown_Metrics_22_January_2026.shp"),
  list(date = "2026-01-29", shp = "E:/Remote Sensing Media/13. 29 January 2026/09. Crown Metrics/Crown_Metrics_29_January_2026.shp"),
  list(date = "2026-02-06", shp = "E:/Remote Sensing Media/14. 06 February 2026/09. Crown Metrics/Crown_Metrics_06_February_2026.shp"),
  list(date = "2026-02-16", shp = "E:/Remote Sensing Media/15. 16 February 2026/09. Crown Metrics/Crown_Metrics_16_February_2026.shp"),
  list(date = "2026-02-23", shp = "E:/Remote Sensing Media/16. 23 February 2026/09. Crown Metrics/Crown_Metrics_23_February_2026.shp"),
  list(date = "2026-03-02", shp = "E:/Remote Sensing Media/17. 02 March 2026/09. Crown Metrics/Crown_Metrics_02_March_2026.shp"),
  list(date = "2026-03-09", shp = "E:/Remote Sensing Media/18. 09 March 2026/09. Crown Metrics/Crown_Metrics_09_March_2026.shp"),
  list(date = "2026-03-16", shp = "E:/Remote Sensing Media/19. 16 March 2026/09. Crown Metrics/Crown_Metrics_16_March_2026.shp"),
  list(date = "2026-03-23", shp = "E:/Remote Sensing Media/20. 23 March 2026/09. Crown Metrics/Crown_Metrics_23_March_2026.shp"),
  list(date = "2026-03-31", shp = "E:/Remote Sensing Media/21. 31 March 2026/09. Crown Metrics/Crown_Metrics_31_March_2026.shp"),
  list(date = "2026-04-08", shp = "E:/Remote Sensing Media/22. 08 April 2026/09. Crown Metrics/Crown_Metrics_08_April_2026.shp"),
  list(date = "2026-04-13", shp = "E:/Remote Sensing Media/23. 13 April 2026/09. Crown Metrics/Crown_Metrics_13_April_2026.shp"),
  list(date = "2026-04-23", shp = "E:/Remote Sensing Media/24. 23 April 2026/09. Crown Metrics/Crown_Metrics_23_April_2026.shp"),
  list(date = "2026-04-29", shp = "E:/Remote Sensing Media/25. 29 April 2026/09. Crown Metrics/Crown_Metrics_29_April_2026.shp"),
  list(date = "2026-05-08", shp = "E:/Remote Sensing Media/26. 08 May 2026/09. Crown Metrics/Crown_Metrics_08_May_2026.shp"),
  list(date = "2026-05-25", shp = "E:/Remote Sensing Media/27. 25 May 2026/09. Crown Metrics/Crown_Metrics_25_May_2026.shp")
)


# =============================================================================
# 3. MORAN'S I FUNCTION
# =============================================================================

run_moran <- function(data, target_date, shp_path, resid_col, response_name, model_type) {
  
  crowns <- st_read(shp_path, quiet = TRUE) %>%
    mutate(Tree = round(as.numeric(Tree), 2))
  
  df_sub <- data %>%
    filter(Date == as.Date(target_date)) %>%
    mutate(Tree = round(as.numeric(Tree), 2))
  
  df_sub$target_resid <- df_sub[[resid_col]]
  
  joined <- crowns %>%
    left_join(
      df_sub %>% select(Compartment, Line, Plot, Tree, target_resid),
      by = c("Cmprtmn" = "Compartment", "Line", "Plot", "Tree")
    ) %>%
    filter(!is.na(target_resid))
  
  # Need at least k+1 observations to build weights matrix
  if (nrow(joined) < 7) {
    warning(paste("Skipping", target_date, response_name, model_type, "- too few observations"))
    return(NULL)
  }
  
  coords  <- st_coordinates(st_centroid(joined))
  nb      <- knn2nb(knearneigh(coords, k = 6))
  lw      <- nb2listw(nb, style = "W")
  m_test  <- moran.test(joined$target_resid, lw)
  
  tibble(
    Date         = as.Date(target_date),
    days         = as.numeric(as.Date(target_date) - T0),
    Response     = response_name,
    Model        = model_type,
    Moran_I_Stat = round(m_test$estimate["Moran I statistic"], 4),
    Expected_I   = round(m_test$estimate["Expectation"],       6),
    Variance     = round(m_test$estimate["Variance"],          8),
    Z_Score      = round(m_test$statistic,                     4),
    p_value_raw  = m_test$p.value,
    p_value      = format.pval(m_test$p.value, eps = 0.001, digits = 3),
    Significant  = m_test$p.value < 0.05,
    Interpretation = case_when(
      m_test$p.value >= 0.05                                       ~ "Non-significant",
      abs(m_test$estimate["Moran I statistic"]) < 0.1              ~ "Significant - negligible effect",
      abs(m_test$estimate["Moran I statistic"]) < 0.3              ~ "Significant - weak effect",
      TRUE                                                          ~ "Significant - moderate/strong effect"
    )
  )
}


# =============================================================================
# 4. RUN ALL COMBINATIONS (26 dates x 6 models = 156 tests)
# =============================================================================

model_specs <- list(
  list(data = df_h, resid_col = "resid_sp", response = "Height",     model = "Species"),
  list(data = df_h, resid_col = "resid_sc", response = "Height",     model = "Spacing"),
  list(data = df_c, resid_col = "resid_sp", response = "Crown Area", model = "Species"),
  list(data = df_c, resid_col = "resid_sc", response = "Crown Area", model = "Spacing"),
  list(data = df_r, resid_col = "resid_sp", response = "CA:H Ratio", model = "Species"),
  list(data = df_r, resid_col = "resid_sc", response = "CA:H Ratio", model = "Spacing")
)

total_tests <- length(flight_registry) * length(model_specs)
cat(paste0("Running ", total_tests, " Moran's I tests across ",
           length(flight_registry), " dates x ", length(model_specs), " models...\n"))
cat("This will take several minutes. Progress below:\n\n")

results_list <- list()
counter      <- 0

for (flight in flight_registry) {
  for (spec in model_specs) {
    counter <- counter + 1
    cat(sprintf("  [%3d/%d] %s - %s (%s model)\n",
                counter, total_tests,
                flight$date, spec$response, spec$model))
    
    row <- tryCatch(
      run_moran(spec$data, flight$date, flight$shp,
                spec$resid_col, spec$response, spec$model),
      error = function(e) {
        warning(paste("Error at", flight$date, spec$response, spec$model, ":", e$message))
        NULL
      }
    )
    if (!is.null(row)) results_list[[length(results_list) + 1]] <- row
  }
}

moran_results <- bind_rows(results_list)

# Factor ordering for plots
moran_results <- moran_results %>%
  mutate(
    Response = factor(Response,
                      levels = c("Height", "Crown Area", "CA:H Ratio"),
                      labels = c("Calibrated height (m)",
                                 "Crown area (m2)",
                                 "CA:H ratio (m2 m-1)")),
    Model = factor(Model, levels = c("Species", "Spacing"))
  )


# =============================================================================
# 5. CONSOLE SUMMARY
# =============================================================================

cat("\n================ MORAN'S I RESULTS (CONDITIONAL RESIDUALS) ================\n")
print(
  as.data.frame(moran_results %>%
                  select(Date, Response, Model, Moran_I_Stat, Z_Score, p_value, Interpretation)),
  row.names = FALSE
)
cat("===========================================================================\n\n")

cat("-- Range summary ----------------------------------------------------------\n")
moran_results %>%
  group_by(Response, Model) %>%
  summarise(
    Min_I   = round(min(Moran_I_Stat),  4),
    Max_I   = round(max(Moran_I_Stat),  4),
    Mean_I  = round(mean(Moran_I_Stat), 4),
    N_sig   = sum(Significant),
    N_total = n(),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# =============================================================================
# 6. THESIS THEME (matches GAMM script exactly)
# =============================================================================

theme_thesis <- function() {
  theme_classic(base_size = 9, base_family = "Calibri") +
    theme(
      plot.title         = element_text(size = 10, face = "bold"),
      plot.subtitle      = element_text(size = 9,  colour = "grey40"),
      axis.title         = element_text(size = 9),
      axis.text          = element_text(size = 8),
      axis.line          = element_line(colour = "black", linewidth = 1.2),
      axis.ticks         = element_line(colour = "black", linewidth = 1),
      axis.ticks.length  = unit(4, "pt"),
      panel.grid.major.y = element_line(colour = alpha("#b0b0b0", 0.25),
                                        linewidth = 0.5, linetype = "solid"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position      = "top",
      legend.justification = "left",
      legend.background    = element_rect(fill = "white", colour = "lightgray",
                                          linewidth = 0.5),
      legend.title         = element_text(size = 10, face = "bold"),
      legend.text          = element_text(size = 8),
      legend.key.size      = unit(0.4, "cm"),
      legend.margin        = margin(t = 2, r = 5, b = 2, l = 5, unit = "pt"),
      plot.margin          = margin(t = 2, r = 5, b = 2, l = 2, unit = "pt"),
      strip.background     = element_blank(),
      strip.text           = element_text(size = 9, face = "bold")
    )
}

model_line_colors <- c(
  "Species" = "#2E4057",
  "Spacing" = "#E05A00"
)


# =============================================================================
# 7. MORAN'S I TEMPORAL LINE CHART
# =============================================================================

make_moran_panel <- function(data, response_label, keep_x = FALSE) {
  
  df_panel <- data %>% filter(Response == response_label)
  
  p <- ggplot(df_panel,
              aes(x = days, y = Moran_I_Stat,
                  colour = Model, group = Model)) +
    
    geom_hline(yintercept = 0,    linetype = "solid",  colour = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8, shape = 16) +
    
    geom_point(data = df_panel %>% filter(Significant),
               aes(x = days, y = Moran_I_Stat),
               shape = 1, size = 3.5, stroke = 0.8, colour = "#CC0000") +
    
    scale_colour_manual(values = model_line_colors, name = "Model") +
    
    scale_x_continuous(
      breaks = seq(0, 270, by = 60),
      limits = c(0, NA)
    ) +
    
    scale_y_continuous(
      limits = c(
        min(0, min(data$Moran_I_Stat, na.rm = TRUE) - 0.02),
        max(data$Moran_I_Stat, na.rm = TRUE) + 0.05
      ),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    
    labs(y = paste0("Moran's I  -  ", as.character(response_label))) +
    
    theme_thesis()
  
  if (!keep_x) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  } else {
    p <- p + labs(x = "Days from 1 September 2025")
  }
  
  if (response_label != levels(data$Response)[1]) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}

response_levels <- levels(moran_results$Response)

p_h <- make_moran_panel(moran_results, response_levels[1], keep_x = FALSE)
p_c <- make_moran_panel(moran_results, response_levels[2], keep_x = FALSE)
p_r <- make_moran_panel(moran_results, response_levels[3], keep_x = TRUE)

p_moran_temporal <- p_h / p_c / p_r +
  plot_layout(heights = c(1, 1, 1))

ggsave(
  file.path(OUTPUT_DIR, "moran_temporal_trend.png"),
  p_moran_temporal,
  width  = 8,
  height = 7,
  units  = "in",
  dpi    = 300
)

cat("Saved: moran_temporal_trend.png\n")


# =============================================================================
# 8. EXCEL EXPORT
# =============================================================================

output_xlsx <- file.path(OUTPUT_DIR, "EucVision_Morans_I_Results.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Morans_I_All_Dates")
addWorksheet(wb, "Range_Summary")
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
style_num4 <- createStyle(
  fontName = "Arial", fontSize = 10, numFmt = "0.0000",
  halign = "center", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC", borderStyle = "thin"
)
style_sig_strong <- createStyle(
  fontName = "Arial", fontSize = 10, fontColour = "#C0392B",
  halign = "left", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC", borderStyle = "thin"
)
style_sig_weak <- createStyle(
  fontName = "Arial", fontSize = 10, fontColour = "#E67E22",
  halign = "left", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC", borderStyle = "thin"
)
style_nonsig <- createStyle(
  fontName = "Arial", fontSize = 10, fontColour = "#27AE60",
  halign = "left", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC", borderStyle = "thin"
)

# Sheet 1
export_df <- moran_results %>%
  select(Date, days, Response, Model,
         Moran_I_Stat, Expected_I, Variance,
         Z_Score, p_value, Significant, Interpretation) %>%
  arrange(Response, Model, Date)

writeData(wb, "Morans_I_All_Dates",
          "EucVision GAMM - Moran's I Spatial Autocorrelation (All Dates, Conditional Residuals)",
          startRow = 1, startCol = 1)
addStyle(wb, "Morans_I_All_Dates", style_title, rows = 1, cols = 1)
mergeCells(wb, "Morans_I_All_Dates", cols = 1:11, rows = 1)

writeData(wb, "Morans_I_All_Dates", export_df,
          startRow = 2, startCol = 1,
          headerStyle = style_header, withFilter = TRUE)

n_data    <- nrow(export_df)
data_rows <- 3:(n_data + 2)

addStyle(wb, "Morans_I_All_Dates", style_center,
         rows = data_rows, cols = c(1:4, 9:10), gridExpand = TRUE, stack = FALSE)
addStyle(wb, "Morans_I_All_Dates", style_num4,
         rows = data_rows, cols = 5:8, gridExpand = TRUE, stack = FALSE)

for (i in seq_len(n_data)) {
  row_i    <- i + 2
  interp_i <- export_df$Interpretation[i]
  sty <- if (grepl("moderate|strong", interp_i)) style_sig_strong else
    if (grepl("weak",             interp_i)) style_sig_weak   else
      style_nonsig
  addStyle(wb, "Morans_I_All_Dates", sty, rows = row_i, cols = 11, stack = FALSE)
}

setColWidths(wb, "Morans_I_All_Dates", cols = 1:11,
             widths = c(13, 8, 24, 10, 13, 13, 13, 10, 10, 12, 36))
setRowHeights(wb, "Morans_I_All_Dates", rows = 1,         heights = 30)
setRowHeights(wb, "Morans_I_All_Dates", rows = 2,         heights = 22)
setRowHeights(wb, "Morans_I_All_Dates", rows = data_rows, heights = 18)
freezePane(wb, "Morans_I_All_Dates", firstActiveRow = 3)

# Sheet 2
summary_df <- moran_results %>%
  group_by(Response, Model) %>%
  summarise(
    Min_I   = round(min(Moran_I_Stat),  4),
    Max_I   = round(max(Moran_I_Stat),  4),
    Mean_I  = round(mean(Moran_I_Stat), 4),
    N_Sig   = sum(Significant),
    N_Total = n(),
    Pct_Sig = paste0(round(100 * sum(Significant) / n(), 1), "%"),
    .groups = "drop"
  )

writeData(wb, "Range_Summary",
          "Moran's I Range Summary by Response and Model",
          startRow = 1, startCol = 1)
addStyle(wb, "Range_Summary", style_title, rows = 1, cols = 1)
mergeCells(wb, "Range_Summary", cols = 1:8, rows = 1)
writeData(wb, "Range_Summary", summary_df,
          startRow = 2, startCol = 1, headerStyle = style_header)
setColWidths(wb, "Range_Summary", cols = 1:8,
             widths = c(28, 12, 10, 10, 10, 10, 10, 10))
setRowHeights(wb, "Range_Summary", rows = 1, heights = 30)
setRowHeights(wb, "Range_Summary", rows = 2, heights = 22)
setRowHeights(wb, "Range_Summary", rows = 3:(nrow(summary_df) + 2), heights = 18)

# Sheet 3
notes_df <- data.frame(
  Item = c(
    "Script version",
    "Residual type",
    "Why conditional residuals?",
    "Crown / CA:H back-transformation",
    "Spatial weights",
    "Dates tested",
    "Significance threshold",
    "Interpretation thresholds",
    "Key reference"
  ),
  Detail = c(
    paste("Generated:", Sys.time()),
    "Conditional: observed - predict(model, type='response'). All model terms included.",
    "Default residuals() on bam() are marginal and exclude random effect fitted values. Plot_ID encodes spatial clustering by design, inflating Moran's I if not accounted for.",
    "Crown and CA:H fitted on log scale (Gaussian identity link). predict() returns log-scale values; exp() applied before differencing against raw observed values.",
    "k-nearest neighbours, k = 6, row-standardised (style='W'). Based on crown polygon centroids extracted from UAV shapefiles.",
    paste("All", length(flight_registry), "UAV flight dates from 2025-09-01 to 2026-05-25"),
    "p < 0.05 (two-sided). Open red rings on temporal plot indicate significant dates.",
    "Non-sig: p >= 0.05 | Negligible: |I| < 0.1 | Weak: 0.1 <= |I| < 0.3 | Moderate/Strong: |I| >= 0.3",
    "Bivand et al. (2013) Applied Spatial Data Analysis with R. Legendre & Fortin (1989) Vegetatio 80:107-138."
  )
)

writeData(wb, "Notes", notes_df, startRow = 1, startCol = 1,
          headerStyle = createStyle(fontName = "Arial", fontSize = 11,
                                    textDecoration = "bold",
                                    fgFill = "#2E4057", fontColour = "white"))
setColWidths(wb, "Notes", cols = 1:2, widths = c(32, 100))

saveWorkbook(wb, output_xlsx, overwrite = TRUE)
cat(paste0("Excel saved: ", output_xlsx, "\n"))


# =============================================================================
# 9. SAVED OUTPUT SUMMARY
# =============================================================================

cat("\n-- Saved outputs ----------------------------------------------------------\n")
cat("  FIGURE\n")
cat("    moran_temporal_trend.png  (3-panel, 8 x 7 in, 300 dpi)\n")
cat("  EXCEL\n")
cat("    EucVision_Morans_I_Results.xlsx\n")
cat("      Sheet 1: All results (", nrow(moran_results), "rows)\n")
cat("      Sheet 2: Range summary by response + model\n")
cat("      Sheet 3: Methodology notes\n")
cat("-- Done -------------------------------------------------------------------\n")