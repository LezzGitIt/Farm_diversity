# Digitized land-cover metrics around each assemblage (EXPLORATORY) ----

### Exploratory test of whether the hand-digitized land-cover layers (Ch1-ssp-birds 00a_Extract_lcs.R output) add information beyond this repo's existing covariates. For each assemblage -- [data collector . farm . year-group . season] -- the 300 m buffers of its point counts are dissolved into one footprint (so overlapping buffers are not counted twice), the land-cover polygons inside it are rasterized at 1 m, and the class percentages (pland: forest / intpast / ssp / other) and total edge (te) are computed with landscapemetrics, matching Ch1's 00b_LSM.R. te scales with footprint area, so edge density (ed = te / covered ha) is the comparable metric across assemblages.

### Inputs read from Ch1 (not re-digitized, not re-run): `../Ch1-ssp-birds/Derived/Geospatial/shp/R_processed/Snapped_{middle,past,ubc,el_hatico}.gpkg` (2026-10-22 rerun). Each layer holds the polygons within 300 m of each point, copied per `Id_survey`; polygons are selected here spatially (those intersecting the footprint), so a point need not be keyed in the layer as long as a neighbour's copy covers it.

### One period layer per assemblage, used whole (no image_date filtering -- each layer is a mosaic of imagery years). Matches 00b: past for 2013 surveys (2016 for G-MB-G-CB), ubc for 2022, el_hatico for the 2025 UBC El Hatico survey, middle otherwise. Choices not covered by 00b's rule: the 2025-26 UBC assemblages and the 2024 UBC-Gaica re-surveys take ubc (the latest digitization; Aaron, 2026-10-22), and any assemblage with a point not in its period layer -- under any collector's id -- falls back to middle (the 2013 assemblages outside the past layer, UBC-Gaica 2024 on La Pradera 1 / La Cristalina / El Porvenir 1, and UBC 2025-26 on La Cristalina / El Porvenir 1, whose points have no digitization newer than middle).

### OUTPUTS
###   Derived/Excels/Landcover_assemblage_metrics.csv -- assemblage -> layer, footprint area, coverage, pland x 4, te, ed (cached; delete to rebuild)
###   Figures/Landcover_*.png -- coverage, correlation with existing covariates, descriptive plots

# Setup ----
library(tidyverse)
library(terra)
library(landscapemetrics)
library(cowplot)

source("Scripts/Shared_paths.R")   # WRANGLING_EXCELS, tax_div_path()
source("Scripts/Shared_fns.R")     # assemblage_habitat_mix()

ggplot2::theme_set(theme_cowplot())
terraOptions(progress = 0)

dir.create("Figures", showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

# Configuration ----

lc_dir <- "../Ch1-ssp-birds/Derived/Geospatial/shp/R_processed"
lc_periods <- c("middle", "past", "ubc", "el_hatico")

## Point-count buffer radius (m) -- the snapped layers only extend 300 m from each point, so this is also the maximum
buffer_radius_m <- 300

## Rasterization resolution (m), as in 00b
raster_res_m <- 1

## Land-cover classes and their integer codes in the rasterized landscape
lc_classes <- c("forest", "intpast", "other", "ssp")

## Projected CRS of the snapped layers (MAGNA-SIRGAS / Colombia Bogota)
lc_crs <- "EPSG:3116"

metrics_csv <- "Derived/Excels/Landcover_assemblage_metrics.csv"

# Point counts and assemblage membership ----

pc_coords <- read_csv(file.path(WRANGLING_EXCELS, "Site_covs.csv"), show_col_types = FALSE) %>%
  mutate(Id_scr = as.character(Id_scr)) %>%
  select(Id_survey_no_dc, Id_scr, Long, Lat, Habitat)

## The assemblages with a bird diversity estimate (the fuller all_farms export), as in 02a
target_assemblages <- read_csv(tax_div_path("all_farms"), show_col_types = FALSE) %>%
  distinct(Assemblage)

assemblage_pcs <- read_csv(file.path(WRANGLING_EXCELS, "Event_covs.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Year)) %>%
  left_join(pc_coords, by = "Id_survey_no_dc") %>%
  filter(!is.na(Id_scr)) %>%
  mutate(Assemblage = str_replace_all(paste(Uniq_db, Id_scr, Year_grp, Season, sep = "."), " |-", "_")) %>%
  semi_join(target_assemblages, by = "Assemblage") %>%
  distinct(Assemblage, Id_scr, Uniq_db, Year, Id_group, Id_survey, Id_survey_no_dc, Long, Lat, Habitat)

# Load the land-cover layers ----

## past: one polygon class was left as `setosfor`; recode to ssp, as 00b does
lc_layers <- map(set_names(lc_periods), \(period) {
  layer <- vect(file.path(lc_dir, paste0("Snapped_", period, ".gpkg")))
  if (period == "past") layer$lc_typ2 <- if_else(layer$lc_typ == "setosfor", "ssp", layer$lc_typ2)
  layer$lc_code <- match(layer$lc_typ2, lc_classes)
  stopifnot(!anyNA(layer$lc_code))
  layer[, c("Id_survey", "lc_typ2", "lc_code")]
})
## Point locations each layer covers: its Id_survey keys mapped to Id_survey_no_dc
survey_location <- read_csv(file.path(WRANGLING_EXCELS, "Event_covs.csv"), show_col_types = FALSE) %>%
  distinct(Id_survey, Id_survey_no_dc)
layer_locations <- map(lc_layers, \(layer) {
  survey_location$Id_survey_no_dc[survey_location$Id_survey %in% unique(layer$Id_survey)]
})

# Choose one period layer per assemblage ----

## Period by survey year (00b's rule, plus 2025-26 UBC and 2024 UBC-Gaica -> ubc), then fall back to middle unless every point of the assemblage is in that layer. A point counts as in a layer if its location (Id_survey_no_dc) is keyed there under any collector's Id_survey -- the UBC-Gaica ids are never keyed in ubc, but their locations are.
period_for_year <- function(uniq_db, year, id_group) {
  case_when(
    year == 2013 ~ "past",
    year == 2016 & str_detect(id_group, "G-MB-G-CB") ~ "past",
    uniq_db == "Ubc mbd" & year == 2022 ~ "ubc",
    uniq_db == "Ubc mbd" & str_detect(id_group, "-EH") ~ "el_hatico",
    uniq_db == "Ubc mbd" & year >= 2025 ~ "ubc",
    uniq_db == "Ubc gaica mbd" ~ "ubc",
    .default = "middle"
  )
}

assemblage_layer <- assemblage_pcs %>%
  mutate(period_rule = period_for_year(Uniq_db, Year, Id_group)) %>%
  mutate(keyed = map2_lgl(Id_survey_no_dc, period_rule, \(id, period) id %in% layer_locations[[period]])) %>%
  summarize(period_rule = first(period_rule), all_keyed = all(keyed),
            Num_pc = n_distinct(Id_survey_no_dc), .by = c(Assemblage, Id_scr, Uniq_db, Year)) %>%
  mutate(lc_period = if_else(all_keyed, period_rule, "middle"))

stopifnot(nrow(assemblage_layer) == n_distinct(assemblage_pcs$Assemblage))
print(count(assemblage_layer, period_rule, lc_period, all_keyed))

# Metric extraction ----

## Class percentages and total edge for one footprint. `points` = SpatVector of an assemblage's point counts (projected); `lc_layer` = land-cover polygons with an integer `lc_code` column. Overlapping buffers are dissolved; overlapping polygon copies rasterize to one cell value. Returns footprint area, the % of it covered by digitized polygons, pland per class (% of covered area) and te (m; edges with uncovered cells are not counted).
footprint_landcover_metrics <- function(points, lc_layer, radius_m, res_m, classes) {
  footprint <- aggregate(buffer(points, width = radius_m))
  polys <- lc_layer[footprint, ]
  template <- rast(ext(footprint), resolution = res_m, crs = crs(lc_layer))
  footprint_r <- rasterize(footprint, template)
  n_footprint <- global(footprint_r, "notNA")[[1]]
  out <- tibble(footprint_ha = n_footprint * res_m^2 / 1e4, covered_pct = 0, te = NA_real_)
  if (nrow(polys) == 0) return(out)
  lc_r <- rasterize(polys, template, field = "lc_code") %>%
    mask(footprint_r)
  n_covered <- global(lc_r, "notNA")[[1]]
  pland <- lsm_c_pland(lc_r) %>%
    mutate(class = classes[class]) %>%
    select(class, value) %>%
    pivot_wider(names_from = class, values_from = value)
  out %>%
    mutate(covered_pct = 100 * n_covered / n_footprint,
           te = lsm_l_te(lc_r)$value) %>%
    bind_cols(pland)
}

pts_proj <- vect(assemblage_pcs, geom = c("Long", "Lat"), crs = "EPSG:4326") %>%
  project(lc_crs)

extract_assemblage <- function(assemblage_id) {
  period <- assemblage_layer$lc_period[assemblage_layer$Assemblage == assemblage_id]
  points <- pts_proj[pts_proj$Assemblage == assemblage_id, ]
  footprint_landcover_metrics(points, lc_layers[[period]], buffer_radius_m, raster_res_m, lc_classes) %>%
    mutate(Assemblage = assemblage_id, lc_period = period, .before = 1)
}

## Cached: the extraction takes a few minutes. Delete the CSV to rebuild (e.g. after the digitization or the radius changes).
if (file.exists(metrics_csv)) {
  Lc_metrics <- read_csv(metrics_csv, show_col_types = FALSE)
} else {
  Lc_metrics <- map(assemblage_layer$Assemblage, extract_assemblage, .progress = TRUE) %>%
    list_rbind() %>%
    mutate(across(all_of(lc_classes), ~ replace_na(.x, 0)),
           ed = te / (footprint_ha * covered_pct / 100)) %>%
    left_join(assemblage_layer %>% select(Assemblage, Id_scr, Uniq_db, Year, Num_pc), by = "Assemblage") %>%
    select(Assemblage, Id_scr, Uniq_db, Year, Num_pc, lc_period, footprint_ha, covered_pct,
           all_of(lc_classes), te, ed) %>%
    mutate(across(c(footprint_ha, covered_pct, all_of(lc_classes), ed), ~ round(.x, 2)),
           te = round(te)) %>%
    arrange(Assemblage)
  write_csv(Lc_metrics, metrics_csv)
}

cat("\nLand-cover metrics:", nrow(Lc_metrics), "assemblages,", n_distinct(Lc_metrics$Id_scr), "farms.\n")
print(summary(Lc_metrics %>% select(footprint_ha, covered_pct, all_of(lc_classes), ed)))

# Cross-check against Ch1's point-level metrics ----

## Ch1's Event_covs_lsm.csv holds the same pland / te per point x year (300 m circle). The assemblage footprint value should track the mean of its points' values; it differs because shared buffer area is counted once here.
Ch1_points <- read_csv("../Ch1-ssp-birds/Derived/Excels/Event_covs_lsm.csv", show_col_types = FALSE) %>%
  distinct(Id_survey, Year, forest, intpast, other, ssp, te)

Point_mean <- assemblage_pcs %>%
  distinct(Assemblage, Id_survey, Year) %>%
  left_join(Ch1_points, by = c("Id_survey", "Year")) %>%
  summarize(across(c(forest, ssp, te), \(x) mean(x, na.rm = TRUE), .names = "{.col}_ptmean"),
            n_pts_ch1 = sum(!is.na(forest)), .by = Assemblage)

Cross_check <- Lc_metrics %>% left_join(Point_mean, by = "Assemblage")
cat("\nFootprint vs mean-of-points (Ch1) correlation: forest r =",
    round(cor(Cross_check$forest, Cross_check$forest_ptmean, use = "complete.obs"), 3),
    "; ssp r =", round(cor(Cross_check$ssp, Cross_check$ssp_ptmean, use = "complete.obs"), 3),
    "; ed vs te/28.3 ha r =", round(cor(Cross_check$ed, Cross_check$te_ptmean, use = "complete.obs"), 3), "\n")

# Compare with the existing covariates ----

## Assemblage covariates from 04a's persisted modelling frame (one row per assemblage; richness rows carry the q = 0 response), plus local-scale WVCC canopy (02a), Num.hab, and the share of surveyed points in each habitat
Model_frame <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
  mutate(Id_scr = as.character(Id_scr))

Asm_covs <- Model_frame %>%
  distinct(Assemblage, Ecoregion, CollectorXyear, dist_farm, Land_use_div, Water_mgmt_div,
           Pasture_mgmt_div, All_practices_div, Elev_mean, Tot_prec_mean, canopy_10k, pool_wes)

Canopy_local <- read_csv("Data/Geospatial/Canopy_by_scale_assemblage.csv", show_col_types = FALSE) %>%
  filter(radius_m %in% c(200, 1000)) %>%
  mutate(radius = paste0("canopy_", radius_m / 1000, "k")) %>%
  select(Assemblage, radius, canopy_cover) %>%
  pivot_wider(names_from = radius, values_from = canopy_cover)

Num_hab <- read_csv(tax_div_path("all_farms"), show_col_types = FALSE) %>%
  distinct(Assemblage, Num_hab = Num.hab)

Hab_share <- assemblage_habitat_mix(
  read_csv(file.path(WRANGLING_EXCELS, "Site_covs.csv"), show_col_types = FALSE),
  read_csv(file.path(WRANGLING_EXCELS, "Event_covs.csv"), show_col_types = FALSE)
) %>%
  select(Assemblage, pts_forest_pct, pts_ssp_pct)

## Log richness at Cmax (q = 0), the main-text response, for the descriptive plots
Richness <- Model_frame %>% filter(Hill == "richness") %>% select(Assemblage, log_richness = log_response)

Lc_compare <- Lc_metrics %>%
  left_join(Asm_covs, by = "Assemblage") %>%
  left_join(Canopy_local, by = "Assemblage") %>%
  left_join(Num_hab, by = "Assemblage") %>%
  left_join(Hab_share, by = "Assemblage") %>%
  left_join(Richness, by = "Assemblage")

## Coverage ----

## Which assemblages / farms have a usable value: digitized polygons over most of the footprint, and in the modelling frame (a farm with diversification indices)
coverage_min_pct <- 90

Lc_compare <- Lc_compare %>%
  mutate(in_model_frame = !is.na(Ecoregion),
         usable = covered_pct >= coverage_min_pct)

cat("\nCoverage of the digitized footprint (% of union of 300 m buffers):\n")
print(quantile(Lc_compare$covered_pct, c(0, .05, .25, .5, 1)))
cat("\nAssemblages by ecoregion (in the modelling frame) and coverage >=", coverage_min_pct, "%:\n")
Coverage_tbl <- Lc_compare %>%
  filter(in_model_frame) %>%
  summarize(assemblages = n(), assemblages_usable = sum(usable), farms = n_distinct(Id_scr),
            farms_usable = n_distinct(Id_scr[usable]), .by = Ecoregion) %>%
  arrange(Ecoregion)
print(Coverage_tbl)
cat("Not in the modelling frame (no diversification row):",
    paste(unique(Lc_compare$Id_scr[!Lc_compare$in_model_frame]), collapse = ", "), "\n")
cat("Low-coverage assemblages:\n")
print(Lc_compare %>% filter(!usable) %>% select(Assemblage, lc_period, covered_pct, footprint_ha))

p_coverage <- Lc_compare %>%
  ggplot(aes(x = covered_pct, y = fct_reorder(Uniq_db, covered_pct), colour = lc_period)) +
  geom_vline(xintercept = coverage_min_pct, linetype = 2, colour = "grey50") +
  geom_jitter(height = 0.15, width = 0, alpha = 0.7) +
  labs(x = "% of footprint digitized", y = NULL, colour = "Layer")
ggsave("Figures/Landcover_coverage.png", p_coverage, width = 7, height = 4, bg = "white")
print(p_coverage)

## Correlations with existing covariates ----

Lc_use <- Lc_compare %>% filter(in_model_frame, usable)

new_vars <- c(forest = "Forest %", ssp = "Silvopasture %", intpast = "Intensive pasture %", ed = "Edge density")
existing_vars <- c(
  canopy_0.2k = "WVCC 200 m", canopy_1k = "WVCC 1 km", canopy_10k = "WVCC 10 km",
  Elev_mean = "Elevation", Tot_prec_mean = "Precipitation", pool_wes = "Species pool (wes)",
  Land_use_div = "Land-use div.", Water_mgmt_div = "Water mgmt div.",
  Pasture_mgmt_div = "Pasture mgmt div.", All_practices_div = "All-practices div.",
  pts_forest_pct = "% points in forest", pts_ssp_pct = "% points in ssp",
  Num_pc = "Num point counts", dist_farm = "Dist. to farm", log_richness = "log richness (q0)"
)

Cor_long <- expand_grid(new = names(new_vars), existing = names(existing_vars)) %>%
  mutate(r = map2_dbl(new, existing, \(a, b) cor(Lc_use[[a]], Lc_use[[b]], method = "spearman", use = "pairwise.complete.obs")),
         new = factor(new_vars[new], levels = new_vars),
         existing = factor(existing_vars[existing], levels = rev(existing_vars)))

cat("\nSpearman correlations, new land-cover metrics vs existing covariates (n =", nrow(Lc_use), "assemblages):\n")
print(Cor_long %>% mutate(r = round(r, 2)) %>% pivot_wider(names_from = new, values_from = r))

cat("\nAmong the new metrics:\n")
print(round(cor(Lc_use[names(new_vars)], method = "spearman"), 2))

p_cor <- Cor_long %>%
  ggplot(aes(x = new, y = existing, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.2) +
  scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c0392b", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman r") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("Figures/Landcover_correlations.png", p_cor, width = 7, height = 6.5, bg = "white")
print(p_cor)

## Variance explained by the existing adjustment sets ----

## R2 of each new metric on (a) Ecoregion, (b) the component block in 04a (elev + elev2 + precip + precip2 + canopy_10k + pool_wes), (c) that block + local WVCC + surveyed-habitat share, (d) the collector-year (dataset) -- 1 - R2 is the variation the existing covariates do not already carry
r2_of <- function(response, rhs, data) {
  summary(lm(reformulate(rhs, response), data = data))$r.squared
}
rhs_sets <- list(
  Ecoregion = "Ecoregion",
  Component = c("poly(Elev_mean, 2)", "poly(Tot_prec_mean, 2)", "canopy_10k", "pool_wes"),
  Component_local = c("poly(Elev_mean, 2)", "poly(Tot_prec_mean, 2)", "canopy_10k", "pool_wes",
                      "canopy_0.2k", "pts_forest_pct", "pts_ssp_pct"),
  Collector_year = "CollectorXyear"
)
R2_tbl <- expand_grid(metric = names(new_vars), set = names(rhs_sets)) %>%
  mutate(R2 = map2_dbl(metric, set, \(m, s) r2_of(m, rhs_sets[[s]], Lc_use))) %>%
  pivot_wider(names_from = set, values_from = R2) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
cat("\nR2 of each new metric on the existing covariate sets:\n")
print(R2_tbl)

## Descriptive plots ----

Lc_long <- Lc_use %>%
  select(Assemblage, Ecoregion, log_richness, canopy_0.2k, pts_forest_pct, all_of(names(new_vars))) %>%
  pivot_longer(all_of(names(new_vars)), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(new_vars[metric], levels = new_vars))

p_by_eco <- Lc_long %>%
  ggplot(aes(x = Ecoregion, y = value, colour = Ecoregion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("Figures/Landcover_by_ecoregion.png", p_by_eco, width = 9, height = 6.5, bg = "white")
print(p_by_eco)

p_vs_richness <- Lc_long %>%
  ggplot(aes(x = value, y = log_richness, colour = Ecoregion)) +
  geom_point(alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6) +
  facet_wrap(~ metric, scales = "free_x") +
  labs(x = NULL, y = "log richness at Cmax (q = 0)")
ggsave("Figures/Landcover_vs_richness.png", p_vs_richness, width = 9, height = 6.5, bg = "white")
print(p_vs_richness)

p_vs_wvcc <- Lc_use %>%
  ggplot(aes(x = canopy_0.2k, y = forest, colour = Ecoregion)) +
  geom_point(alpha = 0.7) +
  labs(x = "WVCC canopy cover, 200 m around hull (%)", y = "Digitized forest, 300 m footprint (%)")
ggsave("Figures/Landcover_forest_vs_wvcc.png", p_vs_wvcc, width = 7, height = 5, bg = "white")
print(p_vs_wvcc)
