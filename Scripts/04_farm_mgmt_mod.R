# Bird taxonomic diversity vs farm management diversification (Bayesian measurement-error models) ----

### Fits Bayesian hierarchical regressions of assemblage-level bird Hill-number diversity on the four farm management diversification indices -- one index at a time -- propagating each iNEXT diversity estimate's bootstrap uncertainty into the response via `resp_se()`.

### For every (response x region-adjustment) a **no-index baseline** is fitted alongside the four index models, so the variance the indices add can be read against the variance region + sampling already explain (`bayes_R2`).

### Region adjustment, two model versions per fit (see Scripts/dag.R for the causal rationale):
###   * "climate" -- PRIMARY. The DAG-sufficient set: farm elevation + precipitation (each as x + x^2) + surrounding 10 km canopy cover. Adjusts the mechanisms Ecoregion stands for without the 5-level factor soaking up management variation.
###   * "ecoregion" -- ROBUSTNESS. Ecoregion as a 5-level fixed effect + canopy. The factor is a coarse proxy for climate (R^2 ~ 0.8) and the only handle on the regional species pool the climate version misses; it also over-adjusts on any region-correlated management variation, so it brackets rather than replaces the primary.
### Both versions carry the sampling covariates below, including Num.hab (number of habitat types surveyed -- a sampling-scope control: 2013 & 2016/17 surveys were forest-only, 2019+ spanned the land-use gradient).

### Responses (all log-transformed, modelled with `resp_se(se, sigma = TRUE)`):
###   * richness  -- q = 0 non-asymptotic `No_Asy_TD` at Cmax + its coverage-based SE (from the coverage65 export)
###   * shannon   -- q = 1 asymptotic `TD_asy` + bootstrap SE (from the all_farms export)
###   * simpson   -- q = 2 asymptotic `TD_asy` + bootstrap SE
### q = 1 / q = 2 use the asymptotic estimate because those profiles asymptote; q = 0 does not, so it uses the coverage-standardised estimate.

### Model per fit:
###   log(response) | resp_se(se_log, sigma = TRUE) ~ [index_z +] <region adjustment> + Num_pc_log_z + Num_hab_z + doy_sin + doy_cos + (1 | Id_gcs) + (1 | CollectorXyear)
### `Id_gcs` is nested in Ecoregion automatically. `CollectorXyear` is the [dataset x year-group] batch (field protocol + which farms + which habitats surveyed). `doy_sin` / `doy_cos` are a cyclic term on the assemblage's mean day of year -- a control for the Nearctic migrant influx (Oct-Mar). Fixed `Year` is omitted (folded into CollectorXyear).

### PRIMARY analysis = assemblages whose point counts average < 300 m from the farm (`dist_threshold`); the full set is refit as a sensitivity for the index models. See Scripts/dag.R and Project_notes.md for the distance-cutoff decision.

# Setup ----
library(tidyverse)
library(brms)

source("Scripts/Farm_diversity_fns.R")

### `conflicted` is deliberately not loaded: its symbol shims break rstan's Stan-model compilation. This script attaches no package that masks the dplyr verbs.

options(mc.cores = 4, brms.backend = "rstan")

dir.create("Derived/models", recursive = TRUE, showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

set.seed(1989)

# Modelling parameters ----

chains <- 4
iter <- 3000
warmup <- 1000
## 0.99 (not 0.97): the small-n ecoregion fits (5-level factor + canopy + sampling + 2 REs on ~75 rows) threw tens of divergences at 0.97
adapt_delta <- 0.99

div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")

## Region-adjustment model versions (Scripts/dag.R). Both carry `canopy_10k_z` -- the
## 10 km scale that best explains diversity (Scripts/06b_scale_of_effect.R) and, per the
## DAG, a confounder to adjust (it blocks the unobserved farmer-values -> landscape-forest
## backdoor), not just a precision term. `climate` is the DAG-sufficient primary; `ecoregion`
## is the proxy robustness check.
specs <- tribble(
  ~spec,        ~region_fixed,                                                                              ~role,
  "climate",    "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",    "primary",
  "ecoregion",  "Ecoregion + canopy_10k_z",                                                                  "robustness"
)

# Load data ----

## Farm-level table: the four diversification indices, ecoregion, farm climate covariates; farms with a bird biodiversity estimate (Scripts/02_match_farm_diversity.R)
Farm_level <- read_csv("Derived/Excels/Farm_diversity_matched.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## Assemblage diversity: q = 1 / q = 2 asymptotic from all_farms, q = 0 non-asymptotic from coverage65 (latest date-stamped exports from Scripts/00_bird_diversity_estimates.R)
Td_asy <- read_csv(latest_file("Derived/Excels", "^Tax_div_all_farms_.*\\.csv$"), show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  filter(Order.q %in% c(1, 2)) %>%
  transmute(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, Num.hab,
            Hill = if_else(Order.q == 1, "shannon", "simpson"),
            response = TD_asy, response_se = `s.e.`)

## q = 0 needs the coverage-based SE (`No_Asy_TD_se`), which only the current
## `Scripts/00_bird_diversity_estimates.R` emits. If the export predates that,
## richness is skipped with a message rather than erroring.
cov65 <- read_csv(latest_file("Derived/Excels", "^Tax_div_coverage65_.*\\.csv$"), show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## q = 0 uses the non-asymptotic estimate. Its measurement-error SE (`No_Asy_TD_se`,
## coverage-based) only comes from the current `Scripts/00_bird_diversity_estimates.R`.
## If the export predates that, q = 0 is still fitted but with a negligible placeholder
## SE (1e-4) -- i.e. an ordinary model, no measurement error. TEMPORARY: re-run once
## `00` finishes so the q = 0 CrIs reflect the estimation uncertainty.
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
responses <- sort(unique(Tax_div_long$Hill))

## Assemblage sampling covariates from the point-count events: number of distinct point counts, mean day of year (migrant-season proxy), and mean distance of the point counts from the farm (Distancia_farm)
wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels/"

## PRIMARY analysis keeps assemblages whose point counts average < 300 m from the farm; the full set is a sensitivity run. (Aaron's choice: drops the clearly-displaced survey groups, keeps GPS jitter; the coefficients barely move -- see Project_notes.md.)
dist_threshold <- 300

Assemblage_covs <- read_csv(paste0(wrangling_excels, "Site_covs.csv"), show_col_types = FALSE) %>%
  select(Id_muestreo_no_dc, Id_gcs, Distancia_farm) %>%
  left_join(read_csv(paste0(wrangling_excels, "Event_covs.csv"), show_col_types = FALSE),
            by = "Id_muestreo_no_dc") %>%
  filter(!is.na(Ano)) %>%
  mutate(Assemblage = str_replace_all(paste(Uniq_db, Id_gcs, Ano_grp, Season, sep = "."), " |-", "_")) %>%
  summarize(Num_pc = n_distinct(Id_muestreo),
            doy = mean(Julian_day, na.rm = TRUE),
            dist_farm = mean(Distancia_farm, na.rm = TRUE), .by = Assemblage)

## Canopy cover in the surrounding 10 km (Scripts/06a_Extract_cc_buff.R) -- the landscape scale that best explains diversity; in both model versions (a confounder per the DAG, not just precision)
Canopy_10k <- read_csv("Data/Geospatial/Canopy_by_scale_assemblage.csv", show_col_types = FALSE) %>%
  filter(radius_m == 10000) %>%
  transmute(Assemblage, canopy_10k = canopy_cover)

# Build the modelling frame ----

## One row per [assemblage x response], with predictors, the log-scale response + SE, and the grouping factors
Model_data <- Tax_div_long %>%
  semi_join(Farm_level, by = "Id_gcs") %>%
  left_join(Assemblage_covs, by = "Assemblage") %>%
  left_join(Canopy_10k, by = "Assemblage") %>%
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

## z-score continuous predictors across the modelling rows (per-response scaling would fragment the interpretation; the row set is near-identical across responses)
predictors_to_scale <- c(div_indices, "Elev_mean", "Avg_temp_mean", "Tot_prec_mean",
                         "Num_pc_log", "Num_hab_num", "canopy_10k")
Model_data <- Model_data %>%
  mutate(across(all_of(predictors_to_scale), ~ as.numeric(scale(.x)), .names = "{.col}_z")) %>%
  rename(Num_hab_z = Num_hab_num_z)

## Persist the frame so Scripts/05_farm_mgmt_plots.R can map standardized axes to raw units, and Scripts/04c_farm_mgmt_piedemonte.R can reuse it (dist_farm is carried for the < 300 m cutoff)
Model_data %>%
  select(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, CollectorXyear, Hill, Ecoregion,
         response, response_se, log_response, se_log, Num_pc, doy, doy_sin, doy_cos, dist_farm,
         all_of(predictors_to_scale), ends_with("_z")) %>%
  write_csv("Derived/Excels/Farm_mgmt_model_data.csv")

# Priors (weakly informative, log-diversity scale) ----

mod_priors <- c(
  prior(student_t(3, 3, 2.5), class = "Intercept"),
  prior(normal(0, 0.75), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# Model formula + data frame for one fit ----

## The focal index is carried in a generic column `focal_z`; index == "baseline" drops it. Num.hab is a standing sampling control in every model. Formulas with identical structure share a compiled Stan model.
formula_for <- function(index, region_fixed) {
  terms <- c(if (index != "baseline") "focal_z",
             region_fixed,
             "Num_pc_log_z", "Num_hab_z", "doy_sin", "doy_cos",
             "(1 | Id_gcs)", "(1 | CollectorXyear)")
  bf(as.formula(paste0("log_response | resp_se(se_log, sigma = TRUE) ~ ",
                       paste(terms, collapse = " + "))))
}

frame_for <- function(hill, index, region_fixed, data_subset) {
  df <- Model_data %>%
    filter(Hill == hill, !is.na(se_log), !is.na(Num_pc_log_z), !is.na(Num_hab_z), !is.na(doy_sin))
  if (data_subset == "primary") df <- df %>% filter(!is.na(dist_farm), dist_farm < dist_threshold)
  if (index != "baseline") {
    df <- df %>% mutate(focal_z = .data[[paste0(index, "_z")]]) %>% filter(!is.na(focal_z))
  }
  if (str_detect(region_fixed, "prec|Elev")) {
    df <- df %>% filter(!is.na(Tot_prec_mean_z), !is.na(Elev_mean_z))
  }
  if (str_detect(region_fixed, "canopy_10k")) df <- df %>% filter(!is.na(canopy_10k_z))
  df
}

# Fit the grid ----

## data_subset: "primary" (point counts average < dist_threshold m from the farm) is the analysis; "full" (all assemblages) is a sensitivity run, index models only.
fit_grid <- expand_grid(
  hill = responses,
  index = c("baseline", div_indices),
  specs,
  data_subset = c("primary", "full")
) %>%
  filter(!(data_subset == "full" & index == "baseline")) %>%   # full-data sensitivity: index models only
  arrange(data_subset, spec, hill, index) %>%
  mutate(key = paste(hill, index, spec, data_subset, sep = "__"),
         structure = paste(spec, if_else(index == "baseline", "baseline", "index"), sep = "__"))

fit_one <- function(hill, index, spec, region_fixed, role, data_subset, key, structure, base_fit) {
  file <- sprintf("Derived/models/mod_%s", key)
  frame <- frame_for(hill, index, region_fixed, data_subset)
  common <- list(chains = chains, iter = iter, warmup = warmup, seed = 1989,
                 control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
                 file = file, file_refit = "on_change")
  if (is.null(base_fit)) {
    do.call(brm, c(list(formula = formula_for(index, region_fixed),
                        data = frame, prior = mod_priors), common))
  } else {
    do.call(update, c(list(object = base_fit, newdata = frame, recompile = FALSE), common))
  }
}

mod_fits <- vector("list", nrow(fit_grid))
names(mod_fits) <- fit_grid$key
base_by_structure <- list()

for (i in seq_len(nrow(fit_grid))) {
  row <- as.list(fit_grid[i, ])
  message(sprintf("[%d/%d] %s", i, nrow(fit_grid), row$key))
  base <- base_by_structure[[row$structure]]
  fit <- fit_one(row$hill, row$index, row$spec, row$region_fixed,
                 row$role, row$data_subset, row$key, row$structure, base_fit = base)
  if (is.null(base)) base_by_structure[[row$structure]] <- fit
  mod_fits[[row$key]] <- fit
}

# Collect coefficients + Bayesian R-squared ----

tidy_fixef <- function(fit) {
  draws <- as_draws_matrix(fit, variable = "^b_", regex = TRUE)
  imap(asplit(draws, 2), ~ tibble(
    term = str_remove(.y, "^b_"),
    estimate = median(.x), conf_low = quantile(.x, 0.05),
    conf_high = quantile(.x, 0.95), p_direction_pos = mean(.x > 0)
  )) %>% list_rbind()
}

Model_summaries <- pmap(fit_grid, function(hill, index, spec, data_subset, key, ...) {
  fit <- mod_fits[[key]]
  r2 <- bayes_R2(fit)[, "Estimate"]
  tidy_fixef(fit) %>%
    mutate(hill = hill, index = index, spec = spec, data_subset = data_subset,
           n_obs = nobs(fit), bayes_R2 = r2,
           max_rhat = round(max(rhat(fit), na.rm = TRUE), 3),
           n_divergent = sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)) %>%
    relocate(hill, index, spec, data_subset)
}) %>%
  list_rbind() %>%
  mutate(across(c(estimate, conf_low, conf_high, p_direction_pos, bayes_R2), ~ round(.x, 4)),
         ## flag the q = 0 rows fitted without a real measurement-error SE
         note = if_else(hill == "richness" & !q0_has_se, "q0 placeholder SE -- rerun after 00", ""))

write_csv(Model_summaries, "Derived/Excels/Farm_mgmt_model_summaries.csv")

if (!q0_has_se) message("NOTE: q = 0 fitted without measurement error (placeholder SE). Re-run after Scripts/00_bird_diversity_estimates.R.")

# Report ----

cat("\nBayesian R-squared: baseline (region + sampling only) vs + each index (primary analysis, dist < ", dist_threshold, " m)\n", sep = "")
Model_summaries %>%
  filter(data_subset == "primary") %>%
  distinct(hill, index, spec, n_obs, bayes_R2) %>%
  pivot_wider(names_from = index, values_from = bayes_R2) %>%
  print(n = Inf)

cat("\nManagement-diversification coefficient (focal_z), primary analysis\n")
Model_summaries %>%
  filter(term == "focal_z", data_subset == "primary") %>%
  select(hill, index, spec, estimate, conf_low, conf_high, p_direction_pos, n_obs) %>%
  arrange(spec, hill, index) %>%
  print(n = Inf)

cat("\nFull-data sensitivity: focal_z, primary (dist < ", dist_threshold, " m) vs full set\n", sep = "")
Model_summaries %>%
  filter(term == "focal_z") %>%
  select(hill, index, spec, data_subset, estimate, conf_low, conf_high, n_obs) %>%
  pivot_wider(names_from = data_subset, values_from = c(estimate, conf_low, conf_high, n_obs)) %>%
  print(n = Inf)

cat("\nConvergence: max R-hat and divergences by fit\n")
Model_summaries %>%
  distinct(hill, index, spec, data_subset, max_rhat, n_divergent) %>%
  filter(max_rhat > 1.01 | n_divergent > 0) %>%
  print(n = Inf)
