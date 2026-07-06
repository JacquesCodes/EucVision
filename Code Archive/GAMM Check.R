# Strip out Species + Spacing_f from a single raw reference date — same
# logic your GAMM uses, just a plain lm() since we only have one date
build_raw_residual <- function(df_full, target_date, raw_col) {
  df_date <- df_full %>%
    filter(Date == as.Date(target_date)) %>%
    filter(!is.na(.data[[raw_col]]))
  
  m <- lm(as.formula(paste0(raw_col, " ~ Species + Spacing_f")), data = df_date)
  df_date$raw_resid <- residuals(m)
  df_date
}

# Run LISA on the reference-only residual (reuses your existing function
# unchanged — it mean-centres internally, so it works on any numeric column)
tls_check   <- build_raw_residual(df_h, "2025-11-28", "LiDAR_Height_m")
lisa_tls    <- run_local_moran(tls_check, "2025-11-28",
                               flight_registry[[which(sapply(flight_registry, function(f) f$date == "2025-11-28"))]]$shp,
                               "raw_resid", "TLS reference (28 Nov 2025)")

# field_check <- build_raw_residual(df_h, "2026-03-23", "Ground_Truth_Height")
# lisa_field  <- run_local_moran(field_check, "2026-03-23",
#                                flight_registry[[which(sapply(flight_registry, function(f) f$date == "2026-03-23"))]]$shp,
#                                "raw_resid", "Field reference (23 Mar 2026)")

# Compare directly against your existing GAMM-residual LISA results on
# these exact same calendar dates
gamm_tls   <- lisa_results %>% filter(Date == as.Date("2025-11-28"), Response == "Calibrated height (m)")
# gamm_field <- lisa_results %>% filter(Date == as.Date("2026-03-23"), Response == "Calibrated height (m)")

compare_tls <- gamm_tls %>%
  select(Plot_ID, Tree, gamm_cluster = cluster_type) %>%
  inner_join(lisa_tls %>% select(Plot_ID, Tree, ref_cluster = cluster_type),
             by = c("Plot_ID", "Tree"))

cat("Agreement rate, TLS reference vs GAMM residual (28 Nov):\n")
print(mean(compare_tls$gamm_cluster == compare_tls$ref_cluster))
print(table(compare_tls$gamm_cluster, compare_tls$ref_cluster))