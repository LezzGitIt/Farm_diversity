# Match farm management diversity metrics to bird biodiversity estimates ----

# Setup ----
library(tidyverse)
library(readxl)
library(janitor)

source("Scripts/Farm_diversity_fns.R")

# Load data ----

### Farm management diversity indices (Land use / water / pasture / all practices), one row per farm
farm_div_xls <- "Data/Complementary_Biodiversity_Paper_Birds_MJE_June_2026.xls"
## The provider renamed the sheets Sheet1 -> "Data" and Sheet2 -> "Dictionary"; accept either so the script works before and after that file update syncs
sheet_or <- function(preferred, fallback) if (preferred %in% excel_sheets(farm_div_xls)) preferred else fallback
Farm_div_raw <- read_excel(farm_div_xls, sheet = sheet_or("Data", "Sheet1"))
## Column-label lookup provided alongside the indices
Farm_div_labels <- read_excel(farm_div_xls, sheet = sheet_or("Dictionary", "Sheet2"))

### Bird taxonomic-diversity (Hill number) estimates per farm, produced by `Scripts/qmd/02_Analysis_iNEXT.qmd` into Derived/Excels/ (latest date-stamped export is picked up automatically)
Tax_div_all_farms <- read_csv(latest_file("Derived/Excels", "^Tax_div_all_farms_.*\\.csv$"), show_col_types = FALSE)
Tax_div_coverage65 <- read_csv(latest_file("Derived/Excels", "^Tax_div_coverage65_.*\\.csv$"), show_col_types = FALSE)

### Farm-level environmental covariates (ecoregion, climate normals, elevation), built from the point-count site covariates by `Scripts/00_farm_covariates.R`
Farm_covariates <- read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

# Clean farm diversity data ----

## Rename indices to plain-English names using the Sheet2 label lookup, and coerce the farm ID to
## character so it matches the type of Id_gcs in the biodiversity data below
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

## A farm "has an associated biodiversity estimate" if it appears in Tax_div_all_farms, the more
## complete of the two exports (Tax_div_coverage65 is a stricter subset, farms dropped there for low
## sample coverage; see Ch1-ssp-birds/Quarto_docs/02_Analysis_iNEXT.qmd, section "Remove farms low SC")
Farm_div_matched <- Farm_div %>%
  semi_join(Tax_div_all_farms, by = "Id_gcs")

Farm_div_unmatched <- Farm_div %>%
  anti_join(Tax_div_all_farms, by = "Id_gcs")

cat(
  nrow(Farm_div_matched), "of", nrow(Farm_div), "farms in the diversity-index spreadsheet",
  "have an associated bird biodiversity estimate and were kept;",
  nrow(Farm_div_unmatched), "were dropped (no matching Id_gcs in Tax_div_all_farms).\n"
)

## The reverse direction: farms with a bird biodiversity estimate but no matching row in the
## diversity-index spreadsheet (e.g. reference/control sites outside the SCR farm-diversification survey)
Bio_farms_unmatched <- Tax_div_all_farms %>%
  distinct(Id_gcs, Uniq_db, Ano_grp, Season) %>%
  anti_join(Farm_div, by = "Id_gcs") %>%
  arrange(Id_gcs)

cat(
  n_distinct(Bio_farms_unmatched$Id_gcs), "farms have a bird biodiversity estimate",
  "but no matching row in the diversity-index spreadsheet.\n"
)

# Attach a farm-level biodiversity summary for context in later exploration ----

## Non-asymptotic species richness (Order.q == 0) at Cmax, from the high-coverage subset
Farm_richness <- Tax_div_coverage65 %>%
  filter(Order.q == 0) %>%
  summarize(Richness_mean = mean(No_Asy_TD), .by = Id_gcs)

## Asymptotic Shannon/Simpson diversity (Order.q 1-2), from the full farm set
Farm_shannon_simpson <- Tax_div_all_farms %>%
  filter(Order.q %in% c(1, 2)) %>%
  mutate(Hill_num = if_else(Order.q == 1, "Shannon_mean", "Simpson_mean")) %>%
  summarize(Estimate = mean(TD_asy), .by = c(Id_gcs, Hill_num)) %>%
  pivot_wider(names_from = Hill_num, values_from = Estimate)

## Number of [farm x data collector x year x season] assemblages sampled per farm
Farm_n_assemblages <- Tax_div_all_farms %>%
  distinct(Id_gcs, Assemblage) %>%
  count(Id_gcs, name = "N_assemblages")

Farm_div_matched <- Farm_div_matched %>%
  left_join(Farm_richness, by = "Id_gcs") %>%
  left_join(Farm_shannon_simpson, by = "Id_gcs") %>%
  left_join(Farm_n_assemblages, by = "Id_gcs")

# Attach the farm-level environmental covariates ----

## Ecoregion, climate normals and elevation, so the matched file is a single farm-level analysis table for the linking model (drops the point-count coordinates and within-farm SDs, which the model does not need)
Farm_div_matched <- Farm_div_matched %>%
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
