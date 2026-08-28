# Link bird taxonomic diversity to farm management diversification (Bayesian measurement-error models) ----

### Fits Bayesian hierarchical regressions of assemblage-level bird Hill-number diversity on the four farm management diversification indices, one index at a time, propagating the iNEXT bootstrap uncertainty of each diversity estimate into the response via `resp_se()`.

### Two adjustment sets are fitted for every (Hill number x index) pair, so the management coefficient can be read against how region is controlled for:
###   * "ecoregion" -- Ecoregion as a 5-level fixed effect (categorical adjustment; conservative, absorbs all between-region variation).
###   * "climate"   -- farm mean annual precipitation + elevation instead of the label (mechanistic adjustment; leaves more between-region variation for management to attach to, at the risk of residual confounding).

### Response: log of the asymptotic Hill number (`TD_asy`) for q = 1 (Shannon) and q = 2 (Simpson), which the iNEXT4steps profiles in `00_bird_diversity_estimates.R` show do asymptote. Non-asymptotic richness (q = 0) is deferred until `02` re-exports a coverage-based SE for it (bundle with the nboot 100 -> 500 bump).

### Model per fit:
###   log(TD_asy) | resp_se(se_log, sigma = TRUE) ~ index_z + <adjustment> + log(Num_pc)_z + (1 | Id_gcs) + (1 | CollectorXyear)
### `Id_gcs` is nested in Ecoregion automatically (each farm sits in one ecoregion, farm IDs are globally unique). `CollectorXyear` is a batch effect for the 8 [dataset x year-group] cohorts, which used materially different field protocols (see the data-paper summary in Project_notes.md). Fixed `Year` is omitted -- near-collinear with data collector.

# Setup ----
library(tidyverse)
library(brms)
library(tidybayes)

source("Scripts/Farm_diversity_fns.R")

### `conflicted` is deliberately not loaded here: its symbol shims break rstan's Stan-model compilation (rstan scans every package with `apropos()`/`exists()`, tripping on ambiguous names like `ar` / `lag`). This script attaches no package that masks the dplyr verbs, so plain tidyverse ordering is enough.

options(mc.cores = 4, brms.backend = "rstan")

dir.create("Derived/models", recursive = TRUE, showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

set.seed(1989)

# Modelling parameters ----

chains <- 4
iter <- 3000
warmup <- 1000
adapt_delta <- 0.97

hill_numbers <- c(Shannon = 1, Simpson = 2)
div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")
adjustments <- c("ecoregion", "climate")

# Load data ----

## Farm-level table: the four diversification indices, ecoregion, and farm climate covariates, restricted to farms with a bird biodiversity estimate (Scripts/02_match_farm_diversity.R)
Farm_level <- read_csv("Derived/Excels/Farm_diversity_matched.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## Assemblage-level bird diversity estimates with their iNEXT bootstrap SE (latest date-stamped export from 00_bird_diversity_estimates.R)
Tax_div <- read_csv(latest_file("Derived/Excels", "^Tax_div_all_farms_.*\\.csv$"), show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## Number of distinct point-count locations contributing to each assemblage -- a sampling-effort / spatial-heterogeneity covariate that was a strong predictor in Scripts/qmd/03_Farm_diversity.qmd. Built the same way there: distinct Id_muestreo per [Uniq_db . Id_gcs . Ano_grp . Season].
wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels/"
Assemblage_effort <- read_csv(paste0(wrangling_excels, "Site_covs.csv"), show_col_types = FALSE) %>%
  select(Id_muestreo_no_dc, Id_gcs) %>%
  left_join(read_csv(paste0(wrangling_excels, "Event_covs.csv"), show_col_types = FALSE),
            by = "Id_muestreo_no_dc") %>%
  mutate(Assemblage = str_replace_all(
    paste(Uniq_db, Id_gcs, Ano_grp, Season, sep = "."), " |-", "_"
  )) %>%
  summarize(Num_pc = n_distinct(Id_muestreo), .by = Assemblage)

# Build the assemblage-level modelling frame ----

## One row per [assemblage x Hill number], carrying the response, its log-scale SE, the farm predictors, and the grouping factors. Restricted to the matched farms and to q = 1 / q = 2.
Model_data <- Tax_div %>%
  filter(Order.q %in% hill_numbers) %>%
  semi_join(Farm_level, by = "Id_gcs") %>%
  left_join(Assemblage_effort, by = "Assemblage") %>%
  left_join(
    Farm_level %>% select(Id_gcs, Ecoregion, all_of(div_indices),
                          Elev_mean, Avg_temp_mean, Tot_prec_mean),
    by = "Id_gcs"
  ) %>%
  mutate(
    Hill = fct_recode(factor(Order.q), Shannon = "1", Simpson = "2"),
    CollectorXyear = paste(Uniq_db, Ano_grp, sep = "_"),
    log_TD = log(TD_asy),
    ## Delta-method SE of log(TD_asy): SE(log X) ~= SE(X) / X. Adequate where the CV is small (median ~0.07-0.10 here); the few high-CV assemblages get down-weighted, which is the point.
    se_log = `s.e.` / TD_asy,
    Num_pc_log = log(Num_pc)
  )

## z-score the predictors (over the modelling rows) so the priors below are on a common scale and the intercept is interpretable at the mean farm
predictors_to_scale <- c(div_indices, "Elev_mean", "Avg_temp_mean", "Tot_prec_mean", "Num_pc_log")
Model_data <- Model_data %>%
  mutate(across(all_of(predictors_to_scale), ~ as.numeric(scale(.x)), .names = "{.col}_z"))

## Persist the modelling frame (raw + z-scored predictors, response, SE) so Scripts/05_farm_mgmt_plots.R can map standardized axes back to raw index units and draw partial residuals without re-deriving the scaling
Model_data %>%
  select(Assemblage, Id_gcs, Uniq_db, Ano_grp, Season, CollectorXyear, Hill,
         Ecoregion, TD_asy, `s.e.`, log_TD, se_log, Num_pc,
         all_of(predictors_to_scale), ends_with("_z")) %>%
  write_csv("Derived/Excels/Linking_model_data.csv")

# Priors (weakly informative, on the log-diversity scale) ----

link_priors <- c(
  prior(student_t(3, 3, 2.5), class = "Intercept"),
  prior(normal(0, 0.75), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# Model formulas ----

## The focal index is always carried in a generic column `focal_z`, so all "ecoregion" fits share one compiled Stan model and all "climate" fits share another
formula_for <- function(adjustment) {
  rhs <- switch(
    adjustment,
    ecoregion = "focal_z + Ecoregion + Num_pc_log_z + (1 | Id_gcs) + (1 | CollectorXyear)",
    climate   = "focal_z + Tot_prec_mean_z + Elev_mean_z + Num_pc_log_z + (1 | Id_gcs) + (1 | CollectorXyear)"
  )
  bf(as.formula(paste0("log_TD | resp_se(se_log, sigma = TRUE) ~ ", rhs)))
}

## Rows for one (Hill number, index, adjustment) fit: drop assemblages missing the focal index (or the climate covariates for the climate adjustment)
frame_for <- function(hill, index, adjustment) {
  df <- Model_data %>%
    filter(Hill == hill) %>%
    mutate(focal_z = .data[[paste0(index, "_z")]]) %>%
    filter(!is.na(focal_z), !is.na(se_log), !is.na(Num_pc_log_z))
  if (adjustment == "climate") {
    df <- df %>% filter(!is.na(Tot_prec_mean_z), !is.na(Elev_mean_z))
  }
  df
}

# Fit the grid ----

## First fit of each adjustment structure compiles; the rest reuse that compiled model via update(recompile = FALSE). file = caches every fit to Derived/models/ so re-runs load from disk.
fit_grid <- expand_grid(
  hill = names(hill_numbers),
  index = div_indices,
  adjustment = adjustments
) %>%
  arrange(adjustment, hill, index)

fit_link_model <- function(hill, index, adjustment, base_fit) {
  file <- sprintf("Derived/models/link_%s_%s_%s", tolower(hill), index, adjustment)
  frame <- frame_for(hill, index, adjustment)
  if (is.null(base_fit)) {
    brm(
      formula = formula_for(adjustment), data = frame, prior = link_priors,
      chains = chains, iter = iter, warmup = warmup, seed = 1989,
      control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
      file = file, file_refit = "on_change"
    )
  } else {
    ## reuse the already-compiled Stan model of the first fit with this adjustment structure
    update(
      base_fit, newdata = frame, recompile = FALSE,
      chains = chains, iter = iter, warmup = warmup, seed = 1989,
      control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
      file = file, file_refit = "on_change"
    )
  }
}

link_fits <- vector("list", nrow(fit_grid))
names(link_fits) <- with(fit_grid, paste(tolower(hill), index, adjustment, sep = "_"))
base_fits <- list()

for (i in seq_len(nrow(fit_grid))) {
  row <- fit_grid[i, ]
  key <- row$adjustment
  message(sprintf("[%d/%d] %s / %s / %s", i, nrow(fit_grid), row$hill, row$index, row$adjustment))
  base <- base_fits[[key]]
  fit <- fit_link_model(row$hill, row$index, row$adjustment, base_fit = base)
  if (is.null(base)) base_fits[[key]] <- fit
  link_fits[[i]] <- fit
}

# Collect coefficients ----

## Posterior summary of every fixed effect from every fit, plus the posterior probability the effect is positive; the `focal_z` rows are the management-diversification result
tidy_fixef <- function(fit, hill, index, adjustment) {
  draws <- as_draws_matrix(fit, variable = "^b_", regex = TRUE)
  imap(
    asplit(draws, 2),
    ~ tibble(
      term = str_remove(.y, "^b_"),
      estimate = median(.x),
      conf_low = quantile(.x, 0.05),
      conf_high = quantile(.x, 0.95),
      p_direction_pos = mean(.x > 0)
    )
  ) %>%
    list_rbind() %>%
    mutate(Hill = hill, Index = index, Adjustment = adjustment) %>%
    select(Hill, Index, Adjustment, term, estimate, conf_low, conf_high, p_direction_pos)
}

Link_coefficients <- pmap(fit_grid, function(hill, index, adjustment) {
  key <- paste(tolower(hill), index, adjustment, sep = "_")
  tidy_fixef(link_fits[[key]], hill, index, adjustment)
}) %>%
  list_rbind() %>%
  mutate(across(c(estimate, conf_low, conf_high, p_direction_pos), ~ round(.x, 4)))

write_csv(Link_coefficients, "Derived/Excels/Linking_model_coefficients.csv")

# Convergence diagnostics ----

Link_diagnostics <- pmap(fit_grid, function(hill, index, adjustment) {
  key <- paste(tolower(hill), index, adjustment, sep = "_")
  fit <- link_fits[[key]]
  s <- summary(fit)
  tibble(
    Hill = hill, Index = index, Adjustment = adjustment,
    n_obs = nobs(fit),
    max_rhat = max(rhat(fit), na.rm = TRUE),
    min_bulk_ess = min(c(s$fixed[, "Bulk_ESS"], s$random$Id_gcs[, "Bulk_ESS"]), na.rm = TRUE),
    n_divergent = sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
  )
}) %>%
  list_rbind() %>%
  mutate(across(c(max_rhat, min_bulk_ess), ~ round(.x, 3)))

print(Link_diagnostics, n = Inf)
write_csv(Link_diagnostics, "Derived/Excels/Linking_model_diagnostics.csv")

cat("\nManagement-diversification coefficients (focal_z), by spec:\n")
Link_coefficients %>%
  filter(term == "focal_z") %>%
  arrange(Hill, Index, Adjustment) %>%
  print(n = Inf)
