# Build the farm-level environmental covariate table from the point-count site covariates ----

### One-time prep step. Reads the unchanging point-count site covariates from the bird data-wrangling repo (`../Ssp-bird-data-wrangling/`), averages the climate and topography fields to the farm (`Id_gcs`), and writes the result to `Data/Farm_covariates.csv`. That file is then treated as a frozen raw input to this repo (like `Data/Complementary_Biodiversity_Paper_Birds_MJE_June_2026.xls`) -- re-run this script only when the wrangling repo's `Site_covs.csv` changes.

### Why farm-level averaging is safe here: climate and elevation barely vary among a farm's point counts (median within-farm SD ~0.02 C for temperature, ~16 mm for annual precipitation, ~7 m for elevation), so the farm mean loses almost nothing while giving one clean row per farm for the linking model.

# Setup ----
library(tidyverse)

wrangling_site_covs <- "../Ssp-bird-data-wrangling/Derived/Excels/Site_covs.csv"

# Load point-count site covariates ----

### `Site_covs.csv` is one row per point count (`Id_muestreo_no_dc`, 504 rows across 73 farms). `Avg_temp` (mean annual temperature, C) and `Tot_prec` (total annual precipitation, mm) are climate normals extracted at each point-count location; `Elev` is 90 m-resolution elevation (m). `Ecoregion` and `Departamento` are constant within a farm.
Site_covs <- read_csv(wrangling_site_covs, show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

# Aggregate to the farm ----

## Continuous fields: farm mean (centroid for coordinates, mean climate / elevation) plus the within-farm SD so downstream code can see how much spatial spread each mean hides
Farm_covariates <- Site_covs %>%
  summarize(
    Nombre_finca = paste(sort(unique(Nombre_finca)), collapse = "; "),
    Ecoregion = unique(Ecoregion),
    Departamento = unique(Departamento),
    N_point_counts = dplyr::n(),
    Long_mean = mean(Long),
    Lat_mean = mean(Lat),
    Elev_mean = mean(Elev),
    Elev_sd = sd(Elev),
    Avg_temp_mean = mean(Avg_temp),
    Avg_temp_sd = sd(Avg_temp),
    Tot_prec_mean = mean(Tot_prec),
    Tot_prec_sd = sd(Tot_prec),
    .by = Id_gcs
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  arrange(Id_gcs)

# Export as a frozen raw input ----

write_csv(Farm_covariates, "Data/Farm_covariates.csv")

cat(
  "Wrote Data/Farm_covariates.csv:", nrow(Farm_covariates), "farms,",
  sum(Farm_covariates$N_point_counts), "point counts.\n"
)
