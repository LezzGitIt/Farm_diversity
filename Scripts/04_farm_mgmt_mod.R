# Bird taxonomic diversity vs farm management diversification (Bayesian measurement-error models) ----

### Fits Bayesian hierarchical regressions of assemblage-level bird Hill-number diversity on the four farm management diversification indices -- one index at a time -- propagating each iNEXT diversity estimate's bootstrap uncertainty into the response via `resp_se()`.

### For every (response x region-adjustment) a **no-index baseline** is fitted alongside the four index models, so the variance the indices add can be read against the variance region + sampling already explain (`bayes_R2`).

### Region adjustment, two ways per fit:
###   * "ecoregion" -- Ecoregion as a 5-level fixed effect (conservative; absorbs all between-region variation).
###   * "climate"   -- farm mean precipitation + elevation instead of the label (mechanistic; leaves more between-region variation for management, at the risk of residual confounding).
### A third "ecoregion_numhab" spec adds the number of habitat types sampled -- a SENSITIVITY check only: Num.hab is partly a mediator of the land-use index (a more land-use-diversified farm has its points across more habitats), so it is not in the primary models.

### Responses (all log-transformed, modelled with `resp_se(se, sigma = TRUE)`):
###   * richness  -- q = 0 non-asymptotic `No_Asy_TD` at Cmax + its coverage-based SE (from the coverage65 export)
###   * shannon   -- q = 1 asymptotic `TD_asy` + bootstrap SE (from the all_farms export)
###   * simpson   -- q = 2 asymptotic `TD_asy` + bootstrap SE
### q = 1 / q = 2 use the asymptotic estimate because those profiles asymptote; q = 0 does not, so it uses the coverage-standardised estimate.

### Model per fit:
###   log(response) | resp_se(se_log, sigma = TRUE) ~ [index_z +] <region adjustment> + log(Num_pc)_z + doy_sin + doy_cos + (1 | Id_gcs) + (1 | CollectorXyear)
### `Id_gcs` is nested in Ecoregion automatically. `CollectorXyear` is a batch effect for the [dataset x year-group] cohorts (different field protocols). `doy_sin` / `doy_cos` are a cyclic term on the assemblage's mean day of year -- a control for the Nearctic migrant influx (Oct-Mar). Fixed `Year` is omitted (near-collinear with data collector).

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
adapt_delta <- 0.97

div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")

## Region-adjustment specs: the non-focal fixed effects and whether this is a sensitivity spec
specs <- tribble(
  ~spec,               ~region_fixed,                    ~extra_fixed,   ~sensitivity,
  "ecoregion",         "Ecoregion",                      "",             FALSE,
  "climate",           "Tot_prec_mean_z + Elev_mean_z",  "",             FALSE,
  "ecoregion_numhab",  "Ecoregion",                      "Num_hab_z",    TRUE
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

Td_rich <- if ("No_Asy_TD_se" %in% names(cov65)) {
  cov65 %>%
    filter(Order.q == 0) %>%
    transmute(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, Num.hab,
              Hill = "richness", response = No_Asy_TD, response_se = No_Asy_TD_se)
} else {
  message("coverage65 export has no No_Asy_TD_se -- skipping q = 0 richness ",
          "(re-run Scripts/00_bird_diversity_estimates.R to add it).")
  NULL
}

Tax_div_long <- bind_rows(Td_rich, Td_asy)
responses <- sort(unique(Tax_div_long$Hill))

## Assemblage sampling covariates from the point-count events: number of distinct point counts, and mean day of year (migrant-season proxy)
wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels/"
Assemblage_covs <- read_csv(paste0(wrangling_excels, "Site_covs.csv"), show_col_types = FALSE) %>%
  select(Id_muestreo_no_dc, Id_gcs) %>%
  left_join(read_csv(paste0(wrangling_excels, "Event_covs.csv"), show_col_types = FALSE),
            by = "Id_muestreo_no_dc") %>%
  filter(!is.na(Ano)) %>%
  mutate(Assemblage = str_replace_all(paste(Uniq_db, Id_gcs, Ano_grp, Season, sep = "."), " |-", "_")) %>%
  summarize(Num_pc = n_distinct(Id_muestreo),
            doy = mean(Julian_day, na.rm = TRUE), .by = Assemblage)

# Build the modelling frame ----

## One row per [assemblage x response], with predictors, the log-scale response + SE, and the grouping factors
Model_data <- Tax_div_long %>%
  semi_join(Farm_level, by = "Id_gcs") %>%
  left_join(Assemblage_covs, by = "Assemblage") %>%
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
                         "Num_pc_log", "Num_hab_num")
Model_data <- Model_data %>%
  mutate(across(all_of(predictors_to_scale), ~ as.numeric(scale(.x)), .names = "{.col}_z")) %>%
  rename(Num_hab_z = Num_hab_num_z)

## Persist the frame so Scripts/05_farm_mgmt_plots.R can map standardized axes to raw units
Model_data %>%
  select(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, CollectorXyear, Hill, Ecoregion,
         response, response_se, log_response, se_log, Num_pc, doy,
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

## The focal index is carried in a generic column `focal_z`; index == "baseline" drops it. Formulas with identical structure share a compiled Stan model.
formula_for <- function(index, region_fixed, extra_fixed) {
  terms <- c(if (index != "baseline") "focal_z",
             region_fixed,
             if (nzchar(extra_fixed)) extra_fixed,
             "Num_pc_log_z", "doy_sin", "doy_cos",
             "(1 | Id_gcs)", "(1 | CollectorXyear)")
  bf(as.formula(paste0("log_response | resp_se(se_log, sigma = TRUE) ~ ",
                       paste(terms, collapse = " + "))))
}

frame_for <- function(hill, index, region_fixed, extra_fixed) {
  df <- Model_data %>%
    filter(Hill == hill, !is.na(se_log), !is.na(Num_pc_log_z), !is.na(doy_sin))
  if (index != "baseline") {
    df <- df %>% mutate(focal_z = .data[[paste0(index, "_z")]]) %>% filter(!is.na(focal_z))
  }
  if (str_detect(region_fixed, "prec|Elev")) {
    df <- df %>% filter(!is.na(Tot_prec_mean_z), !is.na(Elev_mean_z))
  }
  if (nzchar(extra_fixed)) df <- df %>% filter(!is.na(Num_hab_z))
  df
}

# Fit the grid ----

fit_grid <- expand_grid(
  hill = responses,
  index = c("baseline", div_indices),
  specs
) %>%
  ## sensitivity spec: skip the baseline (only need it to sit the index coefficients against)
  filter(!(sensitivity & index == "baseline")) %>%
  arrange(spec, hill, index) %>%
  mutate(key = paste(hill, index, spec, sep = "__"),
         structure = paste(spec, if_else(index == "baseline", "baseline", "index"), sep = "__"))

fit_one <- function(hill, index, spec, region_fixed, extra_fixed, sensitivity, key, structure, base_fit) {
  file <- sprintf("Derived/models/mod_%s", key)
  frame <- frame_for(hill, index, region_fixed, extra_fixed)
  common <- list(chains = chains, iter = iter, warmup = warmup, seed = 1989,
                 control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
                 file = file, file_refit = "on_change")
  if (is.null(base_fit)) {
    do.call(brm, c(list(formula = formula_for(index, region_fixed, extra_fixed),
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
  fit <- fit_one(row$hill, row$index, row$spec, row$region_fixed, row$extra_fixed,
                 row$sensitivity, row$key, row$structure, base_fit = base)
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

Model_summaries <- pmap(fit_grid, function(hill, index, spec, key, ...) {
  fit <- mod_fits[[key]]
  r2 <- bayes_R2(fit)[, "Estimate"]
  tidy_fixef(fit) %>%
    mutate(hill = hill, index = index, spec = spec,
           n_obs = nobs(fit), bayes_R2 = r2,
           max_rhat = round(max(rhat(fit), na.rm = TRUE), 3),
           n_divergent = sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)) %>%
    relocate(hill, index, spec)
}) %>%
  list_rbind() %>%
  mutate(across(c(estimate, conf_low, conf_high, p_direction_pos, bayes_R2), ~ round(.x, 4)))

write_csv(Model_summaries, "Derived/Excels/Farm_mgmt_model_summaries.csv")

# Report ----

cat("\nBayesian R-squared: baseline (region + sampling only) vs + each index\n")
Model_summaries %>%
  distinct(hill, index, spec, n_obs, bayes_R2) %>%
  filter(!str_detect(spec, "numhab")) %>%
  pivot_wider(names_from = index, values_from = bayes_R2) %>%
  print(n = Inf)

cat("\nManagement-diversification coefficient (focal_z), by response x index x spec\n")
Model_summaries %>%
  filter(term == "focal_z") %>%
  select(hill, index, spec, estimate, conf_low, conf_high, p_direction_pos, n_obs) %>%
  arrange(spec, hill, index) %>%
  print(n = Inf)

cat("\nConvergence: max R-hat and divergences by fit\n")
Model_summaries %>%
  distinct(hill, index, spec, max_rhat, n_divergent) %>%
  filter(max_rhat > 1.01 | n_divergent > 0) %>%
  print(n = Inf)
