# Match farm management diversity metrics to bird biodiversity estimates ----

# Setup ----
library(tidyverse)
library(readxl)
library(janitor)

# Load data ----

### Farm management diversity indices (Land use / water / pasture / all practices), one row per farm
Farm_div_raw <- read_excel("Data/Complementary_Biodiversity_Paper_Birds_MJE_June_2026.xls", sheet = "Sheet1")
## Column-label lookup provided alongside the indices
Farm_div_labels <- read_excel("Data/Complementary_Biodiversity_Paper_Birds_MJE_June_2026.xls", sheet = "Sheet2")

### Bird taxonomic-diversity (Hill number) estimates per farm, produced by the Ch1-ssp-birds repo's
### `Quarto_docs/02_Analysis_iNEXT.qmd` (see Data/Bird_biodiversity_estimates for provenance notes)
Tax_div_all_farms <- read_csv("Data/Bird_biodiversity_estimates/Tax_div_all_farms_06.04.26.csv", show_col_types = FALSE)
Tax_div_coverage65 <- read_csv("Data/Bird_biodiversity_estimates/Tax_div_coverage65_06.04.26.csv", show_col_types = FALSE)

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

# Export ----

dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)
write_csv(Farm_div_matched, "Derived/Excels/Farm_diversity_matched.csv")
write_csv(Farm_div_unmatched, "Derived/Excels/Farm_diversity_unmatched.csv")
