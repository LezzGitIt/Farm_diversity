# Piedemonte-only cut of the bird-diversity ~ farm-management analysis ----

### The full-data analysis (Scripts/04_farm_mgmt_mod.R) is dominated by between-ecoregion differences, and the one management signal (pasture / water) is confounded with region. Restricting to a single ecoregion holds the regional species pool (approximately) constant -- the strongest available test for a real management effect. See Scripts/qmd/Piedemonte_proposal.qmd for the full rationale.

### Piedemonte-specific choices:
###   * WITHIN Piedemonte elevation, precipitation and 10 km canopy are r ~ 0.91-0.94 collinear (one foothill gradient), so they cannot go in the model separately. They are replaced by env_pc1: the first principal component of the three (captures ~95% of their joint variance).
###   * NO "ecoregion" model version -- region is fixed, nothing to bracket against.
###   * env_pc1 is linear only (elevation range ~270-760 m; 25 farms will not support a quadratic).
###   * All continuous predictors are re-standardised WITHIN the Piedemonte set, so a "1 SD" effect is 1 SD among Piedemonte farms.

### Model per fit (same measurement-error structure as 04):
###   log(diversity) | resp_se(se_log, sigma = TRUE) ~ [focal_z +] env_pc1_z + Num_pc_log_z + Num_hab_z + doy_sin + doy_cos + <random effects>
### Primary random effects: (1 | Id_gcs) + (1 | CollectorXyear). A sensitivity drops the CollectorXyear term -- all 8 collector-year batches are present within Piedemonte and the forest-only (2013 / 2016-17) vs land-use-gradient (2019+) sampling change is the main thing that grouping carries, so this checks whether it re-absorbs a management signal.

### Reads Scripts/04_farm_mgmt_mod.R's persisted modelling frame (Farm_mgmt_model_data.csv). Run 04 first.

# Setup ----
library(tidyverse)
library(brms)

source("Scripts/Farm_diversity_fns.R")

## conflicted is deliberately not loaded (its shims break rstan's Stan compilation).

options(mc.cores = 4, brms.backend = "rstan")
dir.create("Derived/models", recursive = TRUE, showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)
set.seed(1989)

# Modelling parameters (matched to 04) ----

chains <- 4
## tighter than 04 (0.99 / 3000): the Piedemonte fits are small-n (17-48 assemblages) so the sampler needs more care, especially the q = 0 richness fits
iter <- 4000
warmup <- 1500
adapt_delta <- 0.999

div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")
dist_threshold <- 300

mod_priors <- c(
  prior(student_t(3, 3, 2.5), class = "Intercept"),
  prior(normal(0, 0.75), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# Build the Piedemonte modelling frame ----

## 04's frame: one row per [assemblage x Hill]. Drop 04's full-data z-scores -- everything is re-standardised within Piedemonte below.
Pied <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  filter(Ecoregion == "Piedemonte") %>%
  select(-ends_with("_z"))

## env_pc1: PC1 of {elevation, precipitation, canopy_10k}, fitted on the distinct Piedemonte farms then applied to every assemblage row.
env_vars <- c("Elev_mean", "Tot_prec_mean", "canopy_10k")
pied_farms <- Pied %>% distinct(Id_gcs, .keep_all = TRUE) %>% select(Id_gcs, all_of(env_vars))
env_pca <- prcomp(pied_farms %>% select(all_of(env_vars)), scale. = TRUE)
pc1_var <- summary(env_pca)$importance["Proportion of Variance", "PC1"]
message(sprintf("env_pc1 captures %.1f%% of the joint variance of elevation / precip / canopy within Piedemonte", 100 * pc1_var))
message("env_pc1 loadings: ", paste(sprintf("%s %.2f", env_vars, env_pca$rotation[, "PC1"]), collapse = "  "))

pied_farms <- pied_farms %>%
  mutate(env_pc1 = as.numeric(predict(env_pca, newdata = select(pied_farms, all_of(env_vars)))[, 1]))

Pied <- Pied %>%
  left_join(pied_farms %>% select(Id_gcs, env_pc1), by = "Id_gcs") %>%
  ## 04's CSV carries `doy` (mean Julian day) but not the cyclic terms -- rebuild them (same as 04)
  mutate(doy_sin = sin(2 * pi * doy / 365),
         doy_cos = cos(2 * pi * doy / 365))

## re-standardise every continuous predictor within the Piedemonte set (1 SD = 1 SD among Piedemonte assemblages)
predictors_to_scale <- c(div_indices, "env_pc1", "Num_pc_log", "Num_hab_num")
Pied <- Pied %>%
  mutate(across(all_of(predictors_to_scale), ~ as.numeric(scale(.x)), .names = "{.col}_z"))

write_csv(Pied, "Derived/Excels/Farm_mgmt_piedemonte_data.csv")

responses <- sort(unique(Pied$Hill))

# Model formula + data frame for one fit ----

## re_spec: "full_re" = both grouping factors (primary); "no_collector" = drop the CollectorXyear batch RE (sensitivity)
re_terms <- c(full_re = "(1 | Id_gcs) + (1 | CollectorXyear)",
              no_collector = "(1 | Id_gcs)")

formula_for <- function(index, re_spec) {
  terms <- c(if (index != "baseline") "focal_z",
             "env_pc1_z", "Num_pc_log_z", "Num_hab_num_z", "doy_sin", "doy_cos",
             re_terms[[re_spec]])
  bf(as.formula(paste0("log_response | resp_se(se_log, sigma = TRUE) ~ ",
                       paste(terms, collapse = " + "))))
}

frame_for <- function(hill, index, data_subset) {
  df <- Pied %>%
    filter(Hill == hill, !is.na(se_log), !is.na(env_pc1_z), !is.na(Num_pc_log_z),
           !is.na(Num_hab_num_z), !is.na(doy_sin))
  if (data_subset == "primary") df <- df %>% filter(!is.na(dist_farm), dist_farm < dist_threshold)
  if (index != "baseline") {
    df <- df %>% mutate(focal_z = .data[[paste0(index, "_z")]]) %>% filter(!is.na(focal_z))
  }
  df
}

# Fit the grid ----

## primary  = < 300 m, both REs, baseline + 4 indices
## sensitivity 1 = full assemblage set, both REs, index models
## sensitivity 2 = < 300 m, no CollectorXyear RE, index models
fit_grid <- bind_rows(
  expand_grid(hill = responses, index = c("baseline", div_indices), data_subset = "primary", re_spec = "full_re"),
  expand_grid(hill = responses, index = div_indices,                data_subset = "full",    re_spec = "full_re"),
  expand_grid(hill = responses, index = div_indices,                data_subset = "primary", re_spec = "no_collector")
) %>%
  mutate(key = paste(hill, index, data_subset, re_spec, sep = "__"),
         structure = paste(re_spec, if_else(index == "baseline", "baseline", "index"), sep = "__"))

fit_one <- function(hill, index, data_subset, re_spec, key, structure, base_fit) {
  ## "pied_" prefix (not "mod_") so Scripts/05_farm_mgmt_plots.R's mod_* glob ignores these
  file <- sprintf("Derived/models/pied_%s", key)
  frame <- frame_for(hill, index, data_subset)
  common <- list(chains = chains, iter = iter, warmup = warmup, seed = 1989,
                 control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
                 file = file, file_refit = "on_change")
  if (is.null(base_fit)) {
    do.call(brm, c(list(formula = formula_for(index, re_spec), data = frame, prior = mod_priors), common))
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
  fit <- fit_one(row$hill, row$index, row$data_subset, row$re_spec, row$key, row$structure, base_fit = base)
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

Model_summaries <- pmap(fit_grid, function(hill, index, data_subset, re_spec, key, ...) {
  fit <- mod_fits[[key]]
  tidy_fixef(fit) %>%
    mutate(hill = hill, index = index, data_subset = data_subset, re_spec = re_spec,
           n_obs = nobs(fit), n_farm = length(unique(fit$data$Id_gcs)),
           bayes_R2 = bayes_R2(fit)[, "Estimate"],
           max_rhat = round(max(rhat(fit), na.rm = TRUE), 3),
           n_divergent = sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)) %>%
    relocate(hill, index, data_subset, re_spec)
}) %>%
  list_rbind() %>%
  mutate(across(c(estimate, conf_low, conf_high, p_direction_pos, bayes_R2), ~ round(.x, 4)))

write_csv(Model_summaries, "Derived/Excels/Farm_mgmt_piedemonte_summaries.csv")

# Report ----

cat("\n== Piedemonte-only farm-management models ==\n")
cat(sprintf("env_pc1: %.0f%% of elev/precip/canopy joint variance; loadings %s\n",
            100 * pc1_var, paste(sprintf("%s %.2f", env_vars, env_pca$rotation[, "PC1"]), collapse = " ")))

cat("\nFocal index coefficient (focal_z), primary analysis (< ", dist_threshold, " m, both REs)\n", sep = "")
Model_summaries %>%
  filter(term == "focal_z", data_subset == "primary", re_spec == "full_re") %>%
  select(hill, index, estimate, conf_low, conf_high, p_direction_pos, n_obs, n_farm) %>%
  arrange(hill, index) %>%
  print(n = Inf)

cat("\nBayesian R-squared: baseline vs + each index (primary)\n")
Model_summaries %>%
  filter(data_subset == "primary", re_spec == "full_re") %>%
  distinct(hill, index, bayes_R2) %>%
  pivot_wider(names_from = index, values_from = bayes_R2) %>%
  print(n = Inf)

cat("\nSensitivity: focal_z across the three fits (primary / full set / no-collector RE)\n")
Model_summaries %>%
  filter(term == "focal_z") %>%
  mutate(fit = case_when(data_subset == "full" ~ "full_set",
                         re_spec == "no_collector" ~ "no_collector_re",
                         TRUE ~ "primary")) %>%
  select(hill, index, fit, estimate, conf_low, conf_high) %>%
  pivot_wider(names_from = fit, values_from = c(estimate, conf_low, conf_high)) %>%
  arrange(hill, index) %>%
  print(n = Inf)

cat("\nConvergence (fits with max R-hat > 1.01 or any divergence)\n")
Model_summaries %>%
  distinct(hill, index, data_subset, re_spec, max_rhat, n_divergent) %>%
  filter(max_rhat > 1.01 | n_divergent > 0) %>%
  print(n = Inf)
