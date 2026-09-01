# Farm-level data: environmental covariates + management diversification indices ----

### Two prep steps, was two scripts (01_farm_covariates + 02_match_farm_diversity):
###   1. Average the bird data-wrangling repo's point-count site covariates (climate normals, elevation, coordinates) to the farm -> Data/Farm_covariates.csv (a frozen raw input; re-run only when Site_covs.csv changes).
###   2. Read Maria Esquivel's four farm-management diversification indices, restrict to farms with a bird-diversity estimate (Scripts/00), and join the biodiversity summary + the covariates from step 1 -> Derived/Excels/Farm_diversity_matched.csv (the single farm-level analysis table).

# Setup ----
library(tidyverse)
library(readxl)
library(janitor)

source("Scripts/Model_fns.R")   # latest_file()

# ==================================================================== #
# 1. Farm-level environmental covariates                                #
# ==================================================================== #

### `Site_covs.csv` is one row per point count (`Id_muestreo_no_dc`, 504 rows across 73 farms). `Avg_temp` (mean annual temp, C) and `Tot_prec` (total annual precip, mm) are climate normals at each point-count location; `Elev` is 90 m-resolution elevation (m). Farm-level averaging is safe: climate / elevation barely vary among a farm's point counts (median within-farm SD ~0.02 C, ~16 mm, ~7 m).
Site_covs <- read_csv("../Ssp-bird-data-wrangling/Derived/Excels/Site_covs.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## Continuous fields: farm mean (centroid for coordinates) + the within-farm SD, so downstream code can see how much spatial spread each mean hides
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

## Frozen raw input to this repo (like the MJE xls); consumed by step 2 below, Scripts/03_exploratory.R, and Scripts/04b
write_csv(Farm_covariates, "Data/Farm_covariates.csv")
cat("Wrote Data/Farm_covariates.csv:", nrow(Farm_covariates), "farms,",
    sum(Farm_covariates$N_point_counts), "point counts.\n")

# ==================================================================== #
# 2. Match the management diversification indices to the diversity data #
# ==================================================================== #

# Load data ----

farm_div_xls <- "Data/Complementary_Biodiversity_Paper_Birds_MJE_June_2026.xls"
## The provider renamed the sheets Sheet1 -> "Data" and Sheet2 -> "Dictionary"; accept either
sheet_or <- function(preferred, fallback) if (preferred %in% excel_sheets(farm_div_xls)) preferred else fallback
Farm_div_raw <- read_excel(farm_div_xls, sheet = sheet_or("Data", "Sheet1"))
Farm_div_labels <- read_excel(farm_div_xls, sheet = sheet_or("Dictionary", "Sheet2"))

### Bird taxonomic-diversity Hill-number estimates per farm (Scripts/00; latest date-stamped export picked up automatically)
Tax_div_all_farms <- read_csv(latest_file("Derived/Excels", "^Tax_div_all_farms_.*\\.csv$"), show_col_types = FALSE)
Tax_div_coverage65 <- read_csv(latest_file("Derived/Excels", "^Tax_div_coverage65_.*\\.csv$"), show_col_types = FALSE)

# Clean the farm diversification data ----

## Rename indices to plain-English names; coerce the farm ID to character to match Id_gcs
Farm_div <- Farm_div_raw %>%
  rename(
    Id_gcs = ID,
    Land_use_div = `[0-1]_Indice_1`,
    Water_mgmt_div = `[0-1]_Indice_4`,
    Pasture_mgmt_div = `[0-1]_Indice_7`,
    All_practices_div = `[0-1]_Riqueza_Practicas_Al..`
  ) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  clean_names(case = "none")

# Match to biodiversity estimates ----

## A farm "has a biodiversity estimate" if it appears in Tax_div_all_farms (the more complete of the two exports)
Farm_div_matched <- Farm_div %>% semi_join(Tax_div_all_farms, by = "Id_gcs")
Farm_div_unmatched <- Farm_div %>% anti_join(Tax_div_all_farms, by = "Id_gcs")
cat(nrow(Farm_div_matched), "of", nrow(Farm_div), "diversification-index farms have a bird biodiversity estimate;",
    nrow(Farm_div_unmatched), "dropped.\n")

## Reverse direction: farms with a bird estimate but no diversification-index row (reference / control sites)
Bio_farms_unmatched <- Tax_div_all_farms %>%
  distinct(Id_gcs, Uniq_db, Ano_grp, Season) %>%
  anti_join(Farm_div, by = "Id_gcs") %>%
  arrange(Id_gcs)
cat(n_distinct(Bio_farms_unmatched$Id_gcs), "farms have a bird estimate but no diversification-index row.\n")

# Attach a farm-level biodiversity summary (context for later exploration) ----

Farm_richness <- Tax_div_coverage65 %>%
  filter(Order.q == 0) %>%
  summarize(Richness_mean = mean(No_Asy_TD), .by = Id_gcs)

Farm_shannon_simpson <- Tax_div_all_farms %>%
  filter(Order.q %in% c(1, 2)) %>%
  mutate(Hill_num = if_else(Order.q == 1, "Shannon_mean", "Simpson_mean")) %>%
  summarize(Estimate = mean(TD_asy), .by = c(Id_gcs, Hill_num)) %>%
  pivot_wider(names_from = Hill_num, values_from = Estimate)

Farm_n_assemblages <- Tax_div_all_farms %>%
  distinct(Id_gcs, Assemblage) %>%
  count(Id_gcs, name = "N_assemblages")

Farm_div_matched <- Farm_div_matched %>%
  left_join(Farm_richness, by = "Id_gcs") %>%
  left_join(Farm_shannon_simpson, by = "Id_gcs") %>%
  left_join(Farm_n_assemblages, by = "Id_gcs") %>%
  ## the farm-level environmental covariates from step 1 (drops the point-count coords + within-farm SDs the model does not need)
  left_join(
    Farm_covariates %>%
      select(Id_gcs, Nombre_finca, Ecoregion, Departamento,
             Elev_mean, Avg_temp_mean, Tot_prec_mean),
    by = "Id_gcs"
  )

stopifnot(sum(is.na(Farm_div_matched$Ecoregion)) == 0)

# Export ----

dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)
write_csv(Farm_div_matched, "Derived/Excels/Farm_diversity_matched.csv")
write_csv(Farm_div_unmatched, "Derived/Excels/Farm_diversity_unmatched.csv")
write_csv(Bio_farms_unmatched, "Derived/Excels/Bio_farms_unmatched.csv")
