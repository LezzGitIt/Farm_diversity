# Bird taxonomic diversity vs farm management diversification -- the primary models ----

### Bayesian hierarchical measurement-error regressions of assemblage-level bird Hill-number diversity on the four farm management diversification indices -- one index at a time -- propagating each iNEXT estimate's bootstrap uncertainty into the response via `resp_se()`.

### PRIMARY DESIGN = two responses x two region-adjustment sets (see Scripts/dag.R):
###   Responses
###     * coverage / Cmax  -- abundance-based; keeps a survey-effort term, log(Num_pc). q = 1 / q = 2 use the asymptotic `TD_asy` (+ bootstrap SE), q = 0 the non-asymptotic `No_Asy_TD` at Cmax (+ coverage SE). Sections below: "Cmax arm".
###     * incidence at m* = 6  -- point-count-standardised (Scripts/00 ESTIMATE A); effort is already in the response, so NO Num_pc term. Sections below: "Incidence arm".
###   Adjustment sets
###     * "climate"   -- PRIMARY (adjustment set 1): elevation + precipitation (each x + x^2) + 10 km landscape forest cover (canopy_10k, Scripts/02a) + the range-rarity-weighted regional species pool (pool_wes, Scripts/01b). Adjusts Ecoregion's mechanisms without the factor soaking up management variation; pool_wes blocks the Ecoregion -> biogeographic history -> species pool backdoor.
###     * "ecoregion" -- ROBUSTNESS (adjustment set 2 = {Ecoregion} exactly): the 5-level factor alone. Closes the landscape-forest and species-pool backdoors itself; over-adjusts region-correlated management variation, so it brackets rather than replaces the primary.
### Both carry cyclic day-of-year + the farm / collector-year random effects. The habitat-count term was dropped 2026-08-31 (outcome-side only, coefficient ~ 0; the CollectorXyear RE carries the forest-only 2013 / 2016-17 surveys).

### For every (response x adjustment) a no-index BASELINE is fitted alongside the four index models, so the variance the indices add reads against what region + sampling already explain (`bayes_R2`).

### PRIMARY analysis = assemblages whose point counts average < 300 m from the farm (`dist_threshold`); the full set is refit as a distance-cutoff sensitivity for the index models, both responses.

### Shared helpers (priors, formula builder, compile-sharing grid fitter, tidiers) are in Scripts/Model_fns.R.

### OUTPUTS
###   Derived/Excels/Farm_mgmt_model_data.csv        -- the persisted Cmax modelling frame (raw + z-scored); reused by Scripts/04b and Scripts/05
###   Derived/Excels/Farm_mgmt_model_summaries.csv   -- Cmax arm: every fixed effect, 90% CrI, bayes_R2, convergence
###   Derived/Excels/Incidence_response_summaries.csv    -- incidence arm: focal coef + fit per model
###   Derived/Excels/Incidence_response_comparison.csv   -- focal coef, Cmax vs incidence, side by side
###   Derived/models/mod_*  (Cmax fits) and Derived/models/inc__*  (incidence fits)
### Figures are built from these tables by Scripts/05_farm_mgmt_plots.R.

# Setup ----
library(tidyverse)
library(brms)

source("Scripts/Model_fns.R")   # MGMT_PRIORS, mgmt_bf(), fit_model_grid(), tidy_model_fits(), focal_fit_summary(), latest_file()

### `conflicted` is deliberately not loaded: its symbol shims break rstan's Stan-model compilation. This script attaches no package that masks the dplyr verbs.

options(mc.cores = 4, brms.backend = "rstan")

dir.create("Derived/models", recursive = TRUE, showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

set.seed(1989)

# Modelling parameters ----

## Cmax arm: 4000 / 1500 / 0.995 -- three ~0.8-collinear climate terms (Elev, Elev^2, pool_wes) on ~80-90 rows; a few richness fits sat at R-hat 1.02-1.03 at 3000 / 0.99.
cmax_iter <- 4000; cmax_warmup <- 1500; cmax_adapt <- 0.995
## Incidence arm: 3000 / 1000 / 0.99.
inc_iter <- 3000; inc_warmup <- 1000; inc_adapt <- 0.99

dist_threshold <- 300

div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")

## Region-adjustment strings (Scripts/dag.R). `climate` = the DAG-sufficient primary; `ecoregion` = the proxy robustness check.
specs <- tribble(
  ~spec,        ~region_fixed,
  "climate",    "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z + pool_wes_z",
  "ecoregion",  "Ecoregion"
)

grid_for <- function(key_prefix) {
  expand_grid(
    hill = c("richness", "shannon", "simpson"),
    index = c("baseline", div_indices),
    specs,
    data_subset = c("primary", "full")
  ) %>%
    filter(!(data_subset == "full" & index == "baseline")) %>%   # full-data sensitivity: index models only
    arrange(data_subset, spec, hill, index) %>%
    mutate(key = paste0(key_prefix, paste(hill, index, spec, data_subset, sep = "__")),
           structure = paste(spec, if_else(index == "baseline", "baseline", "index"), sep = "__"))
}

# ==================================================================== #
# Cmax arm (abundance, coverage-standardised)                           #
# ==================================================================== #

# Load data ----

## Farm-level table: the four diversification indices, ecoregion, farm climate covariates; farms with a bird biodiversity estimate (Scripts/01a_farm_data.R)
Farm_level <- read_csv("Derived/Excels/Farm_diversity_matched.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## Assemblage diversity: q = 1 / q = 2 asymptotic from all_farms, q = 0 non-asymptotic from coverage65 (latest date-stamped exports from Scripts/00_bird_diversity_estimates.R)
Td_asy <- read_csv(latest_file("Derived/Excels", "^Tax_div_all_farms_.*\\.csv$"), show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  filter(Order.q %in% c(1, 2)) %>%
  transmute(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, Num.hab,
            Hill = if_else(Order.q == 1, "shannon", "simpson"),
            response = TD_asy, response_se = `s.e.`)

cov65 <- read_csv(latest_file("Derived/Excels", "^Tax_div_coverage65_.*\\.csv$"), show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## q = 0 uses the non-asymptotic estimate + its coverage-based SE (`No_Asy_TD_se`), which only the current Scripts/00 emits. If the export predates that, q = 0 is still fitted but with a negligible placeholder SE (1e-4) -- an ordinary model, no measurement error. TEMPORARY: re-run once 00 finishes.
q0_has_se <- "No_Asy_TD_se" %in% names(cov65)
if (!q0_has_se) {
  warning("coverage65 export has no No_Asy_TD_se -- q = 0 fitted WITHOUT measurement error ",
          "(placeholder SE). Re-run Scripts/00_bird_diversity_estimates.R and this script.")
}
Td_rich <- cov65 %>%
  filter(Order.q == 0) %>%
  transmute(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, Num.hab,
            Hill = "richness", response = No_Asy_TD,
            response_se = if (q0_has_se) No_Asy_TD_se else 1e-4 * No_Asy_TD)

Tax_div_long <- bind_rows(Td_rich, Td_asy)

## Assemblage sampling covariates from the point-count events: number of distinct point counts, mean day of year (migrant-season proxy), and mean distance of the point counts from the farm (Distancia_farm)
wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels/"

Assemblage_covs <- read_csv(paste0(wrangling_excels, "Site_covs.csv"), show_col_types = FALSE) %>%
  select(Id_muestreo_no_dc, Id_gcs, Distancia_farm) %>%
  left_join(read_csv(paste0(wrangling_excels, "Event_covs.csv"), show_col_types = FALSE),
            by = "Id_muestreo_no_dc") %>%
  filter(!is.na(Ano)) %>%
  mutate(Assemblage = str_replace_all(paste(Uniq_db, Id_gcs, Ano_grp, Season, sep = "."), " |-", "_")) %>%
  summarize(Num_pc = n_distinct(Id_muestreo),
            doy = mean(Julian_day, na.rm = TRUE),
            dist_farm = mean(Distancia_farm, na.rm = TRUE), .by = Assemblage)

## Canopy cover in the surrounding 10 km (Scripts/02a) -- the landscape scale that best explains diversity; a confounder per the DAG, in both specs
Canopy_10k <- read_csv("Data/Geospatial/Canopy_by_scale_assemblage.csv", show_col_types = FALSE) %>%
  filter(radius_m == 10000) %>%
  transmute(Assemblage, canopy_10k = canopy_cover)

## Range-rarity-weighted regional species pool (pool_wes = sum 1000 / sqrt(range); Scripts/01b) -- the DAG's SpeciesPool term, joined by farm; enters the CLIMATE spec only. Unlogged per Aaron (2026-08-31): a joint-environmental adjustment, not a separately interpretable coefficient.
Species_pool <- read_csv("Data/Farm_species_pool.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  select(Id_gcs, pool_wes)

# Build the Cmax modelling frame ----

## One row per [assemblage x response], with predictors, the log-scale response + SE, and the grouping factors
Model_data <- Tax_div_long %>%
  semi_join(Farm_level, by = "Id_gcs") %>%
  left_join(Assemblage_covs, by = "Assemblage") %>%
  left_join(Canopy_10k, by = "Assemblage") %>%
  left_join(Species_pool, by = "Id_gcs") %>%
  left_join(
    Farm_level %>% select(Id_gcs, Ecoregion, all_of(div_indices),
                          Elev_mean, Avg_temp_mean, Tot_prec_mean),
    by = "Id_gcs"
  ) %>%
  mutate(
    CollectorXyear = paste(Uniq_db, Ano_grp, sep = "_"),
    Num_hab_num = as.numeric(as.character(Num.hab)),
    log_response = log(response),
    ## delta-method SE of log(response): SE(log X) ~= SE(X) / X
    se_log = response_se / response,
    Num_pc_log = log(Num_pc),
    ## cyclic day of year -- captures the migrant-season wave without assuming linearity
    doy_sin = sin(2 * pi * doy / 365),
    doy_cos = cos(2 * pi * doy / 365)
  )

## z-score continuous predictors across the modelling rows. Num_hab_num is scaled and carried in the persisted frame for Scripts/04b even though it is no longer a model term.
predictors_to_scale <- c(div_indices, "Elev_mean", "Avg_temp_mean", "Tot_prec_mean",
                         "Num_pc_log", "Num_hab_num", "canopy_10k", "pool_wes")
Model_data <- Model_data %>%
  mutate(across(all_of(predictors_to_scale), ~ as.numeric(scale(.x)), .names = "{.col}_z")) %>%
  rename(Num_hab_z = Num_hab_num_z)

## Persist the frame: Scripts/05 maps standardised axes back to raw units; Scripts/04b (Piedemonte cut, spec checks) reuses it
Model_data %>%
  select(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, CollectorXyear, Hill, Ecoregion,
         response, response_se, log_response, se_log, Num_pc, doy, doy_sin, doy_cos, dist_farm,
         all_of(predictors_to_scale), ends_with("_z")) %>%
  write_csv("Derived/Excels/Farm_mgmt_model_data.csv")

# Fit the Cmax grid ----

## The focal index is carried in a generic column `focal_z`; index == "baseline" drops it. Row filters and the fixed term order below are held stable so unchanged fits reuse their cached `.rds`.
cmax_frame <- function(hill, index, region_fixed, data_subset) {
  df <- Model_data %>%
    filter(Hill == hill, !is.na(se_log), !is.na(Num_pc_log_z), !is.na(doy_sin))
  if (data_subset == "primary") df <- df %>% filter(!is.na(dist_farm), dist_farm < dist_threshold)
  if (index != "baseline") {
    df <- df %>% mutate(focal_z = .data[[paste0(index, "_z")]]) %>% filter(!is.na(focal_z))
  }
  if (str_detect(region_fixed, "prec|Elev")) {
    df <- df %>% filter(!is.na(Tot_prec_mean_z), !is.na(Elev_mean_z))
  }
  if (str_detect(region_fixed, "canopy_10k")) df <- df %>% filter(!is.na(canopy_10k_z))
  if (str_detect(region_fixed, "pool_wes")) df <- df %>% filter(!is.na(pool_wes_z))
  df
}

cmax_grid <- grid_for("mod_")

cmax_fits <- fit_model_grid(
  cmax_grid,
  build_bf = function(row) mgmt_bf(
    "log_response | resp_se(se_log, sigma = TRUE)",
    if (row$index != "baseline") "focal_z" else NULL,
    c(row$region_fixed, "Num_pc_log_z", "doy_sin", "doy_cos", "(1 | Id_gcs)", "(1 | CollectorXyear)")
  ),
  build_data = function(row) cmax_frame(row$hill, row$index, row$region_fixed, row$data_subset),
  iter = cmax_iter, warmup = cmax_warmup, adapt_delta = cmax_adapt
)

Model_summaries <- tidy_model_fits(cmax_fits, cmax_grid) %>%
  select(hill, index, spec, data_subset, term, estimate, conf_low, conf_high,
         p_direction_pos, n_obs, bayes_R2, max_rhat, n_divergent) %>%
  mutate(across(c(estimate, conf_low, conf_high, p_direction_pos, bayes_R2), ~ round(.x, 4)),
         note = if_else(hill == "richness" & !q0_has_se, "q0 placeholder SE -- rerun after 00", ""))

write_csv(Model_summaries, "Derived/Excels/Farm_mgmt_model_summaries.csv")
if (!q0_has_se) message("NOTE: q = 0 fitted without measurement error (placeholder SE). Re-run after Scripts/00.")

# ==================================================================== #
# Incidence arm (point-count-standardised, m* = 6, no Num_pc)           #
# ==================================================================== #

# Load the incidence response + reuse the Cmax predictors ----

## Predictors + grouping factors from the persisted Cmax frame (Hill-invariant per assemblage), so the incidence coefficients are on the same z-scale as the Cmax ones -- only the response and the dropped Num_pc term differ
Predictors <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  distinct(Assemblage, Id_gcs, Ecoregion, CollectorXyear, dist_farm, doy_sin, doy_cos,
           Elev_mean_z, Tot_prec_mean_z, canopy_10k_z, pool_wes_z,
           Land_use_div_z, Water_mgmt_div_z, Pasture_mgmt_div_z, All_practices_div_z)

Incidence <- read_csv(latest_file("Derived/Excels", "^Tax_div_incidence_.*\\.csv$"), show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  transmute(Assemblage, Hill,
            response = qD, response_se = qD_se,
            log_response = log(qD), se_log = qD_se / qD)

Inc_data <- Incidence %>%
  inner_join(Predictors, by = "Assemblage") %>%
  filter(is.finite(se_log), se_log > 0)

# Fit the incidence grid ----

## No Num_pc term -- effort is already standardised in the response
inc_frame <- function(hill, index, region_fixed, data_subset) {
  df <- Inc_data %>% filter(Hill == hill, !is.na(doy_sin))
  if (data_subset == "primary") df <- df %>% filter(!is.na(dist_farm), dist_farm < dist_threshold)
  if (str_detect(region_fixed, "Elev|prec")) {
    df <- df %>% filter(!is.na(Elev_mean_z), !is.na(Tot_prec_mean_z))
  }
  if (str_detect(region_fixed, "canopy_10k")) df <- df %>% filter(!is.na(canopy_10k_z))
  if (str_detect(region_fixed, "pool_wes")) df <- df %>% filter(!is.na(pool_wes_z))
  if (index != "baseline") {
    df <- df %>% mutate(focal_z = .data[[paste0(index, "_z")]]) %>% filter(!is.na(focal_z))
  }
  df
}

inc_grid <- grid_for("inc__")

inc_fits <- fit_model_grid(
  inc_grid,
  build_bf = function(row) mgmt_bf(
    "log_response | resp_se(se_log, sigma = TRUE)",
    if (row$index != "baseline") "focal_z" else NULL,
    c(row$region_fixed, "doy_sin", "doy_cos", "(1 | Id_gcs)", "(1 | CollectorXyear)")
  ),
  build_data = function(row) inc_frame(row$hill, row$index, row$region_fixed, row$data_subset),
  iter = inc_iter, warmup = inc_warmup, adapt_delta = inc_adapt
)

Inc_summaries <- focal_fit_summary(inc_fits, inc_grid) %>%
  select(hill, index, spec, data_subset, n_obs, bayes_R2,
         focal_est, focal_lo, focal_hi, p_direction_pos, max_rhat, n_divergent) %>%
  mutate(across(c(bayes_R2, focal_est, focal_lo, focal_hi), ~ round(.x, 4)),
         p_direction_pos = round(p_direction_pos, 3))

write_csv(Inc_summaries, "Derived/Excels/Incidence_response_summaries.csv")

# ==================================================================== #
# Compare the two responses                                            #
# ==================================================================== #

Abund_focal <- Model_summaries %>%
  filter(term == "focal_z", data_subset == "primary") %>%
  transmute(hill, index, spec, resp = "abundance-Cmax (with Num_pc)",
            focal_est = estimate, focal_lo = conf_low, focal_hi = conf_high,
            p_direction_pos, n_obs)

Inc_focal <- Inc_summaries %>%
  filter(index != "baseline", data_subset == "primary") %>%
  transmute(hill, index, spec, resp = "incidence-at-6-PC (no Num_pc)",
            focal_est, focal_lo, focal_hi, p_direction_pos, n_obs)

Comparison <- bind_rows(Abund_focal, Inc_focal) %>%
  arrange(spec, hill, index, resp)

write_csv(Comparison, "Derived/Excels/Incidence_response_comparison.csv")

### The figures built from these tables -- the Cmax-vs-incidence forest plot and the distance-cutoff sensitivity -- are in Scripts/05_farm_mgmt_plots.R (every management figure lives there).

# ==================================================================== #
# Report                                                               #
# ==================================================================== #

cat("\n== Cmax arm: bayes_R2, baseline vs + each index (primary, dist < ", dist_threshold, " m) ==\n", sep = "")
Model_summaries %>%
  filter(data_subset == "primary") %>%
  distinct(hill, index, spec, n_obs, bayes_R2) %>%
  pivot_wider(names_from = index, values_from = bayes_R2) %>%
  print(n = Inf)

cat("\n== Focal coefficient (focal_z), primary: Cmax vs incidence ==\n")
Comparison %>%
  select(spec, hill, index, resp, focal_est, focal_lo, focal_hi, p_direction_pos) %>%
  print(n = Inf)

cat("\n== Distance-cutoff sensitivity: focal_z, primary (< ", dist_threshold, " m) vs full set ==\n", sep = "")
bind_rows(
  Model_summaries %>% filter(term == "focal_z") %>%
    transmute(response = "cmax", hill, index, spec, data_subset,
              est = estimate, lo = conf_low, hi = conf_high),
  Inc_summaries %>% filter(index != "baseline") %>%
    transmute(response = "incidence", hill, index, spec, data_subset,
              est = focal_est, lo = focal_lo, hi = focal_hi)
) %>%
  pivot_wider(names_from = data_subset, values_from = c(est, lo, hi)) %>%
  arrange(response, spec, hill, index) %>%
  print(n = Inf)

cat("\n== Convergence: fits with max R-hat > 1.01 or divergences ==\n")
bind_rows(
  Model_summaries %>% mutate(arm = "cmax") %>%
    distinct(arm, hill, index, spec, data_subset, max_rhat, n_divergent),
  Inc_summaries %>% mutate(arm = "incidence") %>%
    distinct(arm, hill, index, spec, data_subset, max_rhat, n_divergent)
) %>%
  filter(max_rhat > 1.01 | n_divergent > 0) %>%
  print(n = Inf)
