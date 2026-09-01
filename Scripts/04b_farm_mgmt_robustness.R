# Farm-management models: robustness / specification checks ----

### Every section below stress-tests the Scripts/04a primary result (bird diversity ~ farm-management diversification) by changing ONE thing and checking the focal (management) coefficient. Full rationale + findings for each are in Project_notes.md and the section headers. None of these changes the substantive conclusion.

### SECTIONS (was one script each; set `run_sections` to fit a subset):
###   piedemonte    (was 04c) -- single-ecoregion cut: region fixed, species pool ~ constant. The strongest test for a real effect. Both responses.
###   pool_blocks   (was 04d) -- does the range-map species pool count (pool_point) add anything? Env-adjustment blocks that never combine all three axes; + a precipitation functional-form sub-check.
###   spec_checks   (was 04e) -- drop (1|CollectorXyear); drop resp_se(); + the response scale / likelihood family (log / raw / sqrt / Student-t) as a % effect.
###   response_dist (was 04f) -- which response distribution fits the iNEXT estimates best (LOO on a common scale, posterior-predictive shape checks, how the SE scales with the estimate).
###   ideal_adj     (was 04h) -- the DAG-ideal adjustment (climate + climate^2 + species pool + landscape forest, no Ecoregion factor) + a full collinearity assessment.
###   endemism_pool (was 04i) -- swap the raw pool count for the range-rarity-weighted pool_cwe; does an endemism-weighted pool change the estimate?
###   spline_env    (was 04j) -- replace the elevation / precipitation quadratics with rigid k = 4 splines; is the quadratic's symmetry distorting the focal coefficient?

### `pool_blocks` / `spec_checks` / `response_dist` / `ideal_adj` / `endemism_pool` still carry the PRE-2026-08-31 climate block (Num_hab_z, no pool_wes_z; ecoregion has canopy) -- deliberately frozen (Aaron): their conclusions are about the response / RE / pool, not the exact climate parameterisation. `piedemonte` and `spline_env` are on the current spec.

### Shared helpers: Scripts/Model_fns.R. Reads Scripts/04a's persisted frame Derived/Excels/Farm_mgmt_model_data.csv (+ the incidence export, + Data/Farm_species_pool.csv). Run 04a first.

# Setup ----
library(tidyverse)
library(brms)
library(cowplot)
library(loo)
library(mgcv)

source("Scripts/Model_fns.R")   # MGMT_PRIORS, mgmt_bf(), fit_model_grid(), tidy_model_fits(), focal_fit_summary(), latest_file()

### `conflicted` deliberately not loaded -- its shims break rstan's Stan-model compilation (see Scripts/04a).

options(mc.cores = 4, brms.backend = "rstan")
ggplot2::theme_set(theme_cowplot())

dir.create("Derived/models", recursive = TRUE, showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

set.seed(1989)

## Which sections to run this pass (each is independently cached via file_refit = "on_change")
run_sections <- c("piedemonte", "pool_blocks", "spec_checks", "response_dist",
                  "ideal_adj", "endemism_pool", "spline_env")

# Shared constants + inputs ----

div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")
dist_threshold <- 300
model_data_csv <- "Derived/Excels/Farm_mgmt_model_data.csv"

## Scripts/04a's persisted [assemblage x Hill] frame + the cyclic day-of-year terms (04a stores `doy`, not doy_sin/doy_cos)
Base_frame <- read_csv(model_data_csv, show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs),
         doy_sin = sin(2 * pi * doy / 365),
         doy_cos = cos(2 * pi * doy / 365))

## Per-farm potential species pool metrics (Scripts/01b): raw range-map count + range-rarity-weighted variants
Pool_raw <- read_csv("Data/Farm_species_pool.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

index_lab <- c(Land_use_div = "Land use", Water_mgmt_div = "Water mgmt",
               Pasture_mgmt_div = "Pasture mgmt", All_practices_div = "All practices")
hill_lab  <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)", simpson = "Simpson (q = 2)")


# ================================================================== #
# Section: piedemonte  (was Scripts/04c)                              #
# ================================================================== #
### Restricting to Piedemonte holds the regional species pool ~ constant. WITHIN Piedemonte elevation / precipitation / canopy are r ~ 0.91-0.94 collinear -> replaced by env_pc1 (PC1 of the three, ~95% of their joint variance). No "ecoregion" version (region is fixed). Predictors re-standardised WITHIN Piedemonte. Two responses (Cmax with Num_pc, incidence without); a Cmax-only sensitivity drops the CollectorXyear RE.

if ("piedemonte" %in% run_sections) {
  message("\n=== Section: piedemonte ===")

  pied_iter <- 4000; pied_warmup <- 1500; pied_adapt <- 0.999

  ## Piedemonte assemblages, 04a full-data z-scores dropped (everything re-standardised within Piedemonte below)
  Pied <- Base_frame %>%
    filter(Ecoregion == "Piedemonte") %>%
    select(-ends_with("_z"))

  ## env_pc1: PC1 of {elevation, precipitation, canopy_10k}, fitted on the distinct Piedemonte farms then applied to every assemblage row
  env_vars <- c("Elev_mean", "Tot_prec_mean", "canopy_10k")
  pied_farms <- Pied %>% distinct(Id_gcs, .keep_all = TRUE) %>% select(Id_gcs, all_of(env_vars))
  env_pca <- prcomp(pied_farms %>% select(all_of(env_vars)), scale. = TRUE)
  pc1_var <- summary(env_pca)$importance["Proportion of Variance", "PC1"]
  message(sprintf("env_pc1 captures %.1f%% of the joint variance of elev / precip / canopy within Piedemonte", 100 * pc1_var))
  message("env_pc1 loadings: ", paste(sprintf("%s %.2f", env_vars, env_pca$rotation[, "PC1"]), collapse = "  "))

  pied_farms <- pied_farms %>%
    mutate(env_pc1 = as.numeric(predict(env_pca, newdata = select(pied_farms, all_of(env_vars)))[, 1]))

  Pied <- Pied %>%
    left_join(pied_farms %>% select(Id_gcs, env_pc1), by = "Id_gcs") %>%
    mutate(doy_sin = sin(2 * pi * doy / 365), doy_cos = cos(2 * pi * doy / 365))

  ## Second response: the point-count-standardised incidence estimate (Scripts/00 ESTIMATE A); its models drop Num_pc_log_z
  Incidence <- read_csv(latest_file("Derived/Excels", "^Tax_div_incidence_.*\\.csv$"), show_col_types = FALSE) %>%
    transmute(Assemblage, Hill, inc_log_response = log(qD), inc_se_log = qD_se / qD)
  Pied <- Pied %>% left_join(Incidence, by = c("Assemblage", "Hill"))

  ## re-standardise every continuous predictor within the Piedemonte set
  pied_scale <- c(div_indices, "env_pc1", "Num_pc_log")
  Pied <- Pied %>% mutate(across(all_of(pied_scale), ~ as.numeric(scale(.x)), .names = "{.col}_z"))

  write_csv(Pied, "Derived/Excels/Farm_mgmt_piedemonte_data.csv")

  pied_responses <- sort(unique(Pied$Hill))
  re_terms <- c(full_re = "(1 | Id_gcs) + (1 | CollectorXyear)", no_collector = "(1 | Id_gcs)")

  pied_bf <- function(row) {
    lhs <- if (row$response_type == "incidence") "inc_log_response | resp_se(inc_se_log, sigma = TRUE)"
           else "log_response | resp_se(se_log, sigma = TRUE)"
    mgmt_bf(lhs,
            if (row$index != "baseline") "focal_z" else NULL,
            c("env_pc1_z",
              if (row$response_type == "cmax") "Num_pc_log_z",
              "doy_sin", "doy_cos", re_terms[[row$re_spec]]))
  }

  pied_frame <- function(row) {
    df <- Pied %>% filter(Hill == row$hill, !is.na(env_pc1_z), !is.na(doy_sin))
    df <- if (row$response_type == "incidence") df %>% filter(!is.na(inc_se_log))
          else df %>% filter(!is.na(se_log), !is.na(Num_pc_log_z))
    if (row$data_subset == "primary") df <- df %>% filter(!is.na(dist_farm), dist_farm < dist_threshold)
    if (row$index != "baseline") df <- df %>% mutate(focal_z = .data[[paste0(row$index, "_z")]]) %>% filter(!is.na(focal_z))
    df
  }

  ## grid order (NOT sorted) fixes which fit compiles Stan per structure -- keep the bind_rows order
  pied_grid <- bind_rows(
    expand_grid(hill = pied_responses, index = c("baseline", div_indices), data_subset = "primary", re_spec = "full_re",     response_type = "cmax"),
    expand_grid(hill = pied_responses, index = div_indices,                data_subset = "full",    re_spec = "full_re",     response_type = "cmax"),
    expand_grid(hill = pied_responses, index = div_indices,                data_subset = "primary", re_spec = "no_collector", response_type = "cmax"),
    expand_grid(hill = pied_responses, index = c("baseline", div_indices), data_subset = "primary", re_spec = "full_re",     response_type = "incidence"),
    expand_grid(hill = pied_responses, index = div_indices,                data_subset = "full",    re_spec = "full_re",     response_type = "incidence")
  ) %>%
    mutate(key = paste("pied", hill, index, data_subset, re_spec, response_type, sep = "__"),
           structure = paste(re_spec, response_type, if_else(index == "baseline", "baseline", "index"), sep = "__"))

  pied_fits <- fit_model_grid(pied_grid, build_bf = pied_bf, build_data = pied_frame,
                              iter = pied_iter, warmup = pied_warmup, adapt_delta = pied_adapt)

  Pied_summaries <- tidy_model_fits(pied_fits, pied_grid) %>%
    mutate(n_farm = map_int(key, ~ length(unique(pied_fits[[.x]]$data$Id_gcs)))) %>%
    select(hill, index, data_subset, re_spec, response_type, term, estimate, conf_low, conf_high,
           p_direction_pos, n_obs, n_farm, bayes_R2, max_rhat, n_divergent) %>%
    mutate(across(c(estimate, conf_low, conf_high, p_direction_pos, bayes_R2), ~ round(.x, 4)))

  write_csv(Pied_summaries, "Derived/Excels/Farm_mgmt_piedemonte_summaries.csv")

  cat("\n[piedemonte] env_pc1:", sprintf("%.0f%%", 100 * pc1_var),
      "of elev/precip/canopy joint variance; loadings",
      paste(sprintf("%s %.2f", env_vars, env_pca$rotation[, "PC1"]), collapse = " "), "\n")
  cat("[piedemonte] focal_z, primary (< 300 m, both REs), by response\n")
  Pied_summaries %>%
    filter(term == "focal_z", data_subset == "primary", re_spec == "full_re") %>%
    select(response_type, hill, index, estimate, conf_low, conf_high, p_direction_pos, n_obs, n_farm) %>%
    arrange(response_type, hill, index) %>% print(n = Inf)
  cat("[piedemonte] convergence (max R-hat > 1.01 or divergences)\n")
  Pied_summaries %>% distinct(response_type, hill, index, data_subset, re_spec, max_rhat, n_divergent) %>%
    filter(max_rhat > 1.01 | n_divergent > 0) %>% print(n = Inf)
}


# ================================================================== #
# Section: pool_blocks  (was Scripts/04d)                             #
# ================================================================== #
### Does the range-map species-pool count (pool_point, Scripts/01b) earn a place in 04a? pool_point ~ poly(Elev,2)+poly(precip,2) R^2 = 0.86, so fit FOUR env-adjustment blocks that never put all three axes together and check the focal coefficient + fit. Then a precipitation functional-form sub-check (quadratic vs spline, with / without elevation). Outcome: pool does not move the management coefficient; the pasture/water lean is residual precipitation confounding.

if ("pool_blocks" %in% run_sections) {
  message("\n=== Section: pool_blocks ===")

  pb_iter <- 3000; pb_warmup <- 1000; pb_adapt <- 0.99

  Md <- Base_frame %>%
    left_join(Pool_raw %>% select(Id_gcs, pool_point), by = "Id_gcs") %>%
    mutate(pool_point_z = as.numeric(scale(pool_point)))

  env_blocks <- c(
    clim        = "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
    pool        = "pool_point_z + canopy_10k_z",
    pool_elev   = "pool_point_z + Elev_mean_z + I(Elev_mean_z^2) + canopy_10k_z",
    pool_precip = "pool_point_z + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z"
  )
  pb_sampling <- c("Num_pc_log_z", "Num_hab_z", "doy_sin", "doy_cos", "(1 | Id_gcs)", "(1 | CollectorXyear)")

  pb_frame <- function(row) {
    Md %>%
      filter(Hill == row$hill, !is.na(se_log), !is.na(canopy_10k_z), !is.na(pool_point_z),
             !is.na(Num_pc_log_z), !is.na(Num_hab_z), !is.na(Tot_prec_mean_z), !is.na(Elev_mean_z),
             !is.na(dist_farm), dist_farm < dist_threshold) %>%
      mutate(focal_z = .data[[paste0(row$index, "_z")]]) %>%
      filter(!is.na(focal_z))
  }

  pb_grid <- expand_grid(hill = c("richness", "shannon"), index = div_indices, env = names(env_blocks)) %>%
    mutate(key = paste("pooltest", hill, index, env, sep = "__"), structure = env)

  pb_fits <- fit_model_grid(
    pb_grid,
    build_bf = function(row) mgmt_bf("log_response | resp_se(se_log, sigma = TRUE)", "focal_z",
                                     c(env_blocks[[row$env]], pb_sampling)),
    build_data = pb_frame,
    iter = pb_iter, warmup = pb_warmup, adapt_delta = pb_adapt
  )

  summ <- pmap_dfr(pb_grid, function(hill, index, env, key, ...) {
    f <- pb_fits[[key]]
    b <- fixef(f)["focal_z", ]
    pool_row <- if ("pool_point_z" %in% rownames(fixef(f))) fixef(f)["pool_point_z", ] else rep(NA, 4)
    ess <- tryCatch(min(summary(f)$fixed[, "Bulk_ESS"], na.rm = TRUE), error = function(e) NA)
    tibble(hill, index, env,
           focal = sprintf("%+.3f [%+.3f, %+.3f]", b["Estimate"], b["Q2.5"], b["Q97.5"]),
           focal_est = round(b["Estimate"], 3),
           pool_coef = if (is.na(pool_row[1])) "" else sprintf("%+.3f [%+.3f, %+.3f]", pool_row[1], pool_row[3], pool_row[4]),
           bayes_R2 = round(bayes_R2(f)[, "Estimate"], 3),
           min_bulk_ess = round(ess),
           max_rhat = round(max(rhat(f), na.rm = TRUE), 3),
           n_div = sum(subset(nuts_params(f), Parameter == "divergent__")$Value))
  })
  write_csv(summ, "Derived/Excels/Species_pool_model_test.csv")

  ## precipitation functional-form sub-check (quadratic vs spline; with / without elevation), pasture / water only
  precip_specs <- c(
    "poly2"    = "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
    "splP"     = "Elev_mean_z + I(Elev_mean_z^2) + s(Tot_prec_mean_z, k = 5) + canopy_10k_z",
    "precOnly" = "Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
    "splPonly" = "s(Tot_prec_mean_z, k = 5) + canopy_10k_z"
  )
  precip_lab <- c(poly2 = "Elev^2 + precip^2", splP = "Elev^2 + spline(precip)",
                  precOnly = "precip^2 only (no elevation)", splPonly = "spline(precip) only")

  pflex_grid <- expand_grid(hill = c("richness", "shannon"),
                            index = c("Pasture_mgmt_div", "Water_mgmt_div"),
                            pspec = names(precip_specs)) %>%
    mutate(key = paste("precipflex", hill, index, pspec, sep = "__"), structure = key)

  pflex_fits <- fit_model_grid(
    pflex_grid,
    build_bf = function(row) mgmt_bf("log_response | resp_se(se_log, sigma = TRUE)", "focal_z",
                                     c(precip_specs[[row$pspec]], "Num_pc_log_z", "Num_hab_z",
                                       "doy_sin", "doy_cos", "(1 | Id_gcs)", "(1 | CollectorXyear)")),
    build_data = pb_frame,
    iter = pb_iter, warmup = pb_warmup, adapt_delta = pb_adapt
  )

  pflex <- pmap_dfr(pflex_grid, function(hill, index, pspec, key, ...) {
    b <- fixef(pflex_fits[[key]])["focal_z", ]
    tibble(index = recode(index, Pasture_mgmt_div = "Pasture mgmt", Water_mgmt_div = "Water mgmt"),
           hill = recode(hill, richness = "Richness (q=0)", shannon = "Shannon (q=1)"),
           block = precip_lab[[pspec]],
           focal = sprintf("%+.3f [%+.3f, %+.3f]", b["Estimate"], b["Q2.5"], b["Q97.5"]),
           focal_est = round(b["Estimate"], 3))
  }) %>%
    arrange(index, hill, factor(block, levels = precip_lab))
  write_csv(pflex, "Derived/Excels/Precip_vs_elev_sensitivity.csv")

  # figures (was Scripts/05d) ----

  ## Fig 1: the pasture / water coefficient across every environmental spec, ordered by how completely precipitation is controlled
  pb_spec_lab <- c(pool_elev = "species pool + elev²", clim = "elev² + precip²",
                   splP = "elev² + s(precip)", pool = "species pool",
                   pool_precip = "species pool + precip²", precOnly = "precip² only",
                   splPonly = "s(precip) only")
  pb_hill_lab <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)")

  pull_focal <- function(fit, spec) {
    b <- fixef(fit)["focal_z", ]
    tibble(spec = spec, est = b["Estimate"], lo = b["Q2.5"], hi = b["Q97.5"])
  }
  env_coefs <- bind_rows(
    pmap_dfr(pb_grid, function(hill, index, env, key, ...)
      bind_cols(tibble(hill, index), pull_focal(pb_fits[[key]], env))),
    pmap_dfr(pflex_grid, function(hill, index, pspec, key, ...)
      bind_cols(tibble(hill, index), pull_focal(pflex_fits[[key]], pspec)))
  ) %>%
    filter(index %in% c("Pasture_mgmt_div", "Water_mgmt_div"), spec %in% names(pb_spec_lab)) %>%
    mutate(Index = factor(recode(index, !!!index_lab), levels = c("Pasture mgmt", "Water mgmt")),
           Hill  = factor(recode(hill, !!!pb_hill_lab), levels = unname(pb_hill_lab)),
           Spec  = factor(recode(spec, !!!pb_spec_lab), levels = pb_spec_lab))

  p_env_specs <- env_coefs %>%
    ggplot(aes(est, fct_rev(Spec))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#1b7837", fatten = 2) +
    facet_grid(Index ~ Hill) +
    labs(x = "Diversification-index effect on log diversity (posterior median, 95% CrI)", y = NULL,
         title = "The pasture / water coefficient tracks how completely precipitation is adjusted",
         subtitle = "Environmental term varied; 10 km landscape forest cover + sampling controls held fixed.") +
    theme(axis.text.y = element_text(size = 8.5), plot.subtitle = element_text(size = 9.5))
  ggsave("Figures/Species_pool_forest.png", p_env_specs, width = 10, height = 6, bg = "white")
  print(p_env_specs)

  ## Fig 2: variance partition of the range-map species-pool count onto the elevation / precipitation gradients
  canopy_by_farm <- Base_frame %>% summarize(canopy_10k = mean(canopy_10k, na.rm = TRUE), .by = Id_gcs)
  pe_vp <- Pool_raw %>%
    left_join(read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>% mutate(Id_gcs = as.character(Id_gcs)), by = "Id_gcs") %>%
    left_join(canopy_by_farm, by = "Id_gcs")
  r2v <- function(f) summary(lm(f, pe_vp))$r.squared
  full_vp <- r2v(pool_point ~ poly(Elev_mean, 2) + poly(Tot_prec_mean, 2))
  onlyE <- r2v(pool_point ~ poly(Elev_mean, 2)); onlyP <- r2v(pool_point ~ poly(Tot_prec_mean, 2))
  varpart <- tibble(
    component = c("Unique to precipitation", "Unique to elevation", "Shared elev + precip", "Unexplained"),
    value = c(full_vp - onlyE, full_vp - onlyP, onlyE + onlyP - full_vp, 1 - full_vp)
  ) %>% mutate(component = factor(component, levels = rev(component)))

  p_pool_varpart <- ggplot(varpart, aes(value, "spp_pool", fill = component)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * value)), position = position_stack(vjust = 0.5), size = 3.2, colour = "white") +
    scale_fill_manual(values = c("Unique to precipitation" = "#2166ac", "Unique to elevation" = "#b2182b",
                                 "Shared elev + precip" = "#7f7f7f", "Unexplained" = "grey85"), name = NULL) +
    scale_x_continuous(labels = scales::percent, expand = expansion(0)) +
    labs(x = sprintf("Variance of the species-pool count across the %d farms", nrow(pe_vp)), y = NULL,
         title = "The range-map species-pool count is a composite of the two environmental gradients",
         subtitle = sprintf("~ poly(Elev,2): R² = %.2f   |   ~ poly(precip,2): R² = %.2f   |   both: R² = %.2f", onlyE, onlyP, full_vp)) +
    theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    guides(fill = guide_legend(nrow = 2))
  ggsave("Figures/Species_pool_varpart.png", p_pool_varpart, width = 9, height = 2.8, bg = "white")
  print(p_pool_varpart)

  cat("\n[pool_blocks] variance partition: pool_point ~ climate (farm level, lm)\n")
  pe <- Pool_raw %>%
    left_join(read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>% mutate(Id_gcs = as.character(Id_gcs)),
              by = "Id_gcs")
  r2 <- function(f) round(summary(lm(f, pe))$r.squared, 3)
  full <- r2(pool_point ~ poly(Elev_mean, 2) + poly(Tot_prec_mean, 2))
  cat("  ~ poly(Elev,2):", r2(pool_point ~ poly(Elev_mean, 2)),
      " ~ poly(precip,2):", r2(pool_point ~ poly(Tot_prec_mean, 2)),
      " both:", full, "\n")
  cat("[pool_blocks] focal coefficient by environmental block (primary)\n")
  summ %>% select(hill, index, env, focal) %>% pivot_wider(names_from = env, values_from = focal) %>%
    arrange(hill, index) %>% print(width = Inf)
  cat("[pool_blocks] precipitation functional-form sub-check (pasture / water)\n")
  pflex %>% select(index, hill, block, focal) %>% print(n = Inf)
}


# ================================================================== #
# Section: spec_checks  (was Scripts/04e)                             #
# ================================================================== #
### (a) drop (1 | CollectorXyear) -- is it acting like (1|Ecoregion)? (b) drop resp_se() -- how much does the measurement-error weighting matter? (c) response scale / likelihood family (Gaussian log/raw/sqrt, Student-t log) as a % change in diversity per +1 SD. Outcome: focal coefficient robust to all three.

if ("spec_checks" %in% run_sections) {
  message("\n=== Section: spec_checks ===")

  sc_iter <- 3000; sc_warmup <- 1000; sc_adapt <- 0.99

  Md <- Base_frame %>%
    mutate(sqrt_response = sqrt(response), se_sqrt = response_se / (2 * sqrt(response)))

  sc_env <- c(
    climate   = "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
    ecoregion = "Ecoregion + canopy_10k_z"
  )
  variant_lhs <- c(base = "log_response | resp_se(se_log, sigma = TRUE)",
                   no_collector = "log_response | resp_se(se_log, sigma = TRUE)",
                   no_resp_se = "log_response")
  variant_re  <- c(base = "(1 | Id_gcs) + (1 | CollectorXyear)",
                   no_collector = "(1 | Id_gcs)",
                   no_resp_se = "(1 | Id_gcs) + (1 | CollectorXyear)")

  sc_frame <- function(row) {
    df <- Md %>%
      filter(Hill == row$hill, !is.na(se_log), !is.na(canopy_10k_z), !is.na(Num_pc_log_z),
             !is.na(Num_hab_z), !is.na(dist_farm), dist_farm < dist_threshold)
    if (str_detect(sc_env[[row$env]], "Elev|prec")) df <- df %>% filter(!is.na(Elev_mean_z), !is.na(Tot_prec_mean_z))
    df %>% mutate(focal_z = .data[[paste0(row$index, "_z")]]) %>% filter(!is.na(focal_z))
  }

  sc_grid <- expand_grid(hill = c("richness", "shannon"), index = div_indices,
                         env = names(sc_env), variant = names(variant_lhs)) %>%
    mutate(key = paste("speccheck", hill, index, env, variant, sep = "__"),
           structure = paste(env, variant, sep = "__"))

  sc_fits <- fit_model_grid(
    sc_grid,
    build_bf = function(row) brms::bf(stats::as.formula(paste0(
      variant_lhs[[row$variant]], " ~ focal_z + ", sc_env[[row$env]],
      " + Num_pc_log_z + Num_hab_z + doy_sin + doy_cos + ", variant_re[[row$variant]]))),
    build_data = sc_frame,
    iter = sc_iter, warmup = sc_warmup, adapt_delta = sc_adapt
  )

  summ <- pmap_dfr(sc_grid, function(hill, index, env, variant, key, ...) {
    f <- sc_fits[[key]]
    b <- fixef(f)["focal_z", ]
    re_sd <- tryCatch({
      s <- summary(f)$random
      if (!is.null(s$CollectorXyear)) round(s$CollectorXyear["sd(Intercept)", "Estimate"], 3) else NA_real_
    }, error = function(e) NA_real_)
    tibble(hill, index, env, variant,
           focal = sprintf("%+.3f [%+.3f, %+.3f]", b["Estimate"], b["Q2.5"], b["Q97.5"]),
           focal_est = round(b["Estimate"], 3),
           sigma = round(mean(as_draws_matrix(f, variable = "sigma")), 3),
           collectorXyear_sd = re_sd,
           max_rhat = round(max(rhat(f), na.rm = TRUE), 3))
  })
  write_csv(summ, "Derived/Excels/Spec_checks.csv")

  ## (c) response scale / likelihood family -- pasture / water only, climate spec
  lik_priors_raw <- c(
    prior(student_t(3, 20, 15), class = "Intercept"), prior(normal(0, 5), class = "b"),
    prior(exponential(0.15), class = "sd"), prior(exponential(0.15), class = "sigma"))
  lik_priors_sqrt <- c(
    prior(student_t(3, 5, 4), class = "Intercept"), prior(normal(0, 2), class = "b"),
    prior(exponential(0.5), class = "sd"), prior(exponential(0.5), class = "sigma"))
  lik_scale <- function(lik) dplyr::case_when(lik == "gaussian_raw" ~ "raw", lik == "gaussian_sqrt" ~ "sqrt", TRUE ~ "log")
  lik_rhs <- paste("focal_z +", sc_env[["climate"]],
                   "+ Num_pc_log_z + Num_hab_z + doy_sin + doy_cos + (1 | Id_gcs) + (1 | CollectorXyear)")

  lik_grid <- expand_grid(hill = c("richness", "shannon"),
                          index = c("Pasture_mgmt_div", "Water_mgmt_div"),
                          lik = c("gaussian_log", "student_log", "gaussian_raw", "gaussian_sqrt")) %>%
    mutate(key = paste("likcheck", hill, index, lik, sep = "__"), structure = key)

  lik_lhs <- function(scale) switch(scale,
    raw  = "response | resp_se(response_se, sigma = TRUE)",
    sqrt = "sqrt_response | resp_se(se_sqrt, sigma = TRUE)",
    log  = "log_response | resp_se(se_log, sigma = TRUE)")

  lik_fits <- fit_model_grid(
    lik_grid,
    build_bf = function(row) brms::bf(stats::as.formula(paste(lik_lhs(lik_scale(row$lik)), "~", lik_rhs))),
    build_data = function(row) sc_frame(modifyList(row, list(env = "climate"))),
    prior = function(row, data) switch(lik_scale(row$lik), raw = lik_priors_raw, sqrt = lik_priors_sqrt, MGMT_PRIORS),
    family = function(row) if (row$lik == "student_log") student() else gaussian(),
    iter = sc_iter, warmup = sc_warmup, adapt_delta = sc_adapt
  )

  lik_summ <- pmap_dfr(lik_grid, function(hill, index, lik, key, ...) {
    f <- lik_fits[[key]]
    scale <- lik_scale(lik)
    b <- fixef(f)["focal_z", ]
    ## per-scale % change in diversity per +1 SD of the index (raw/sqrt need the frame's own means)
    df_full <- sc_frame(list(hill = hill, index = index, env = "climate"))
    to_pct <- switch(scale,
                     raw  = function(x) 100 * x / mean(df_full$response),
                     sqrt = function(x) 100 * 2 * mean(df_full$sqrt_response) * x / mean(df_full$response),
                     log  = function(x) 100 * (exp(x) - 1))
    tibble(index = recode(index, Pasture_mgmt_div = "Pasture mgmt", Water_mgmt_div = "Water mgmt"),
           hill = recode(hill, richness = "Richness (q=0)", shannon = "Shannon (q=1)"),
           likelihood = recode(lik, gaussian_log = "Gaussian on log", student_log = "Student-t on log",
                               gaussian_raw = "Gaussian on raw", gaussian_sqrt = "Gaussian on sqrt"),
           pct_effect = sprintf("%+.1f%% [%+.1f, %+.1f]", to_pct(b["Estimate"]), to_pct(b["Q2.5"]), to_pct(b["Q97.5"])),
           pct_est = round(to_pct(b["Estimate"]), 1),
           max_rhat = round(max(rhat(f), na.rm = TRUE), 3))
  })
  write_csv(lik_summ, "Derived/Excels/Likelihood_sensitivity.csv")

  cat("\n[spec_checks] focal coefficient: base vs drop (1|CollectorXyear)\n")
  summ %>% filter(variant %in% c("base", "no_collector")) %>%
    select(hill, index, env, variant, focal_est) %>%
    pivot_wider(names_from = variant, values_from = focal_est) %>%
    mutate(shift = round(no_collector - base, 3)) %>% arrange(env, hill, index) %>% print(n = Inf)
  cat("[spec_checks] focal coefficient: base (resp_se) vs no_resp_se\n")
  summ %>% filter(variant %in% c("base", "no_resp_se")) %>%
    select(hill, index, env, variant, focal_est) %>%
    pivot_wider(names_from = variant, values_from = focal_est) %>%
    mutate(shift = round(no_resp_se - base, 3)) %>% arrange(env, hill, index) %>% print(n = Inf)
  cat("[spec_checks] % change per +1 SD of the index, by likelihood (climate spec)\n")
  lik_summ %>% select(index, hill, likelihood, pct_effect) %>%
    pivot_wider(names_from = likelihood, values_from = pct_effect) %>% arrange(index, hill) %>% print(width = Inf)
}


# ================================================================== #
# Section: response_dist  (was Scripts/04f)                           #
# ================================================================== #
### Which response distribution fits the iNEXT estimates best? Gaussian on {log, raw, sqrt} + Student-t on {log, raw}, all with resp_se(sigma = TRUE), on the baseline + All-practices models, q = 0/1/2, both adjustments. Evaluated by LOO on a common raw-y scale, posterior-predictive shape statistics, and how the bootstrap SE scales with the estimate. Outcome: SE ~ sqrt(estimate), so sqrt fits best; the focal effect is invariant; log kept for interpretability + positivity.

if ("response_dist" %in% run_sections) {
  message("\n=== Section: response_dist ===")

  rd_iter <- 3000; rd_warmup <- 1000; rd_adapt <- 0.99
  rd_chains <- 4; rd_draws_kept <- rd_chains * (rd_iter - rd_warmup)

  rd_specs <- c(
    climate   = "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
    ecoregion = "Ecoregion + canopy_10k_z"
  )
  rd_families <- tribble(
    ~family,          ~fam_obj,   ~scale,
    "gaussian_log",   "gaussian", "log",
    "student_log",    "student",  "log",
    "gaussian_raw",   "gaussian", "raw",
    "student_raw",    "student",  "raw",
    "gaussian_sqrt",  "gaussian", "sqrt"
  )
  fam_lab <- c(gaussian_log = "Gaussian / log", student_log = "Student-t / log",
               gaussian_raw = "Gaussian / raw", student_raw = "Student-t / raw",
               gaussian_sqrt = "Gaussian / sqrt")
  rd_scale_of <- function(family) rd_families$scale[match(family, rd_families$family)]

  to_raw  <- function(x, scale) switch(scale, log = exp(x), sqrt = pmax(x, 0)^2, x)
  obs_raw <- function(fit, scale) switch(scale,
                                         log  = exp(fit$data$log_response),
                                         sqrt = fit$data$sqrt_response^2,
                                         fit$data$response)

  Md <- Base_frame %>%
    mutate(sqrt_response = sqrt(response), se_sqrt = response_se / (2 * sqrt(response)))

  rd_frame <- function(row) {
    df <- Md %>%
      filter(Hill == row$hill, response > 0, !is.na(se_log),
             !is.na(Num_pc_log_z), !is.na(Num_hab_z), !is.na(doy_sin),
             !is.na(canopy_10k_z), !is.na(dist_farm), dist_farm < dist_threshold)
    if (str_detect(rd_specs[[row$spec]], "Elev|prec")) df <- df %>% filter(!is.na(Elev_mean_z), !is.na(Tot_prec_mean_z))
    if (row$index != "baseline") df <- df %>% mutate(focal_z = .data[[paste0(row$index, "_z")]]) %>% filter(!is.na(focal_z))
    df
  }
  rd_rhs <- function(spec, index) paste(c(if (index != "baseline") "focal_z", rd_specs[[spec]],
                                          "Num_pc_log_z", "Num_hab_z", "doy_sin", "doy_cos",
                                          "(1 | Id_gcs)", "(1 | CollectorXyear)"), collapse = " + ")
  rd_lhs <- function(scale) switch(scale,
    log  = "log_response | resp_se(se_log, sigma = TRUE)",
    raw  = "response | resp_se(response_se, sigma = TRUE)",
    sqrt = "sqrt_response | resp_se(se_sqrt, sigma = TRUE)")
  rd_prior <- function(frame, scale) {
    if (scale == "log") return(MGMT_PRIORS)
    y <- if (scale == "sqrt") frame$sqrt_response else frame$response
    m <- stats::median(y); s <- stats::sd(y)
    c(set_prior(sprintf("student_t(3, %.3f, %.3f)", m, 2.5 * s), class = "Intercept"),
      set_prior(sprintf("normal(0, %.3f)", s), class = "b"),
      set_prior(sprintf("exponential(%.5f)", 1 / s), class = "sd"),
      set_prior(sprintf("exponential(%.5f)", 1 / s), class = "sigma"))
  }

  rd_grid <- expand_grid(hill = names(hill_lab), spec = names(rd_specs),
                         index = c("baseline", "All_practices_div"), rd_families) %>%
    mutate(key = paste("respdist", hill, spec, index, family, sep = "__"),
           structure = paste(family, spec, index, sep = "__")) %>%
    arrange(structure, hill)

  rd_fits <- fit_model_grid(
    rd_grid,
    build_bf = function(row) brms::bf(stats::as.formula(paste(rd_lhs(rd_scale_of(row$family)), "~",
                                                              rd_rhs(row$spec, row$index)))),
    build_data = rd_frame,
    prior = function(row, data) rd_prior(data, rd_scale_of(row$family)),
    family = function(row) if (row$fam_obj == "student") student() else gaussian(),
    iter = rd_iter, warmup = rd_warmup, adapt_delta = rd_adapt
  )

  ## (1) LOO on a common raw-y scale (change-of-variables Jacobian)
  pointwise_raw_ll <- function(fit, scale) {
    ll <- log_lik(fit)
    if (scale == "log")  ll <- sweep(ll, 2, log(obs_raw(fit, "log")), "-")
    if (scale == "sqrt") ll <- sweep(ll, 2, log(2 * sqrt(obs_raw(fit, "sqrt"))), "-")
    ll
  }
  loo_one <- function(fit, scale) {
    ll <- pointwise_raw_ll(fit, scale)
    r_eff <- loo::relative_eff(exp(ll), chain_id = rep(seq_len(rd_chains), each = rd_draws_kept / rd_chains))
    suppressWarnings(loo::loo(ll, r_eff = r_eff))
  }
  rd_loo <- imap(rd_fits, function(fit, key) loo_one(fit, rd_scale_of(rd_grid$family[match(key, rd_grid$key)])))

  Loo_table <- rd_grid %>%
    group_by(hill, spec, index) %>%
    group_modify(function(g, ...) {
      cmp <- loo::loo_compare(rd_loo[g$key])
      pos <- match(rownames(cmp), g$key)
      tibble(family = g$family[pos], elpd_loo = cmp[, "elpd_loo"], se_elpd = cmp[, "se_elpd_loo"],
             elpd_diff = cmp[, "elpd_diff"], se_diff = cmp[, "se_diff"],
             p_loo = map_dbl(g$key[pos], ~ rd_loo[[.x]]$estimates["p_loo", "Estimate"]),
             n_pareto_gt_07 = map_int(g$key[pos], ~ sum(loo::pareto_k_values(rd_loo[[.x]]) > 0.7)))
    }) %>%
    ungroup() %>%
    mutate(across(c(elpd_loo, se_elpd, elpd_diff, se_diff, p_loo), ~ round(.x, 2)))
  write_csv(Loo_table, "Derived/Excels/Response_distribution_loo.csv")

  ## (2) posterior-predictive shape statistics (raw scale)
  skewness <- function(x) { x <- x[is.finite(x)]; mean((x - mean(x))^3) / (mean((x - mean(x))^2))^1.5 }
  stat_fns <- list(skew = skewness, sd = sd, min = min, max = max,
                   q05 = function(x) quantile(x, 0.05, names = FALSE),
                   q95 = function(x) quantile(x, 0.95, names = FALSE))
  Ppc_table <- rd_grid %>%
    mutate(scale = rd_scale_of(family)) %>%
    pmap_dfr(function(hill, spec, index, family, key, scale, ...) {
      fit <- rd_fits[[key]]
      yrep <- to_raw(posterior_predict(fit, ndraws = 2000), scale)
      y_obs <- obs_raw(fit, scale)
      imap_dfr(stat_fns, function(fn, nm) {
        t_obs <- fn(y_obs); t_rep <- apply(yrep, 1, fn)
        tibble(stat = nm, obs = t_obs, bayes_p = mean(t_rep >= t_obs))
      }) %>% mutate(hill = hill, spec = spec, index = index, family = family, .before = 1)
    }) %>%
    mutate(obs = round(obs, 3), bayes_p = round(bayes_p, 3))
  write_csv(Ppc_table, "Derived/Excels/Response_distribution_ppc_stats.csv")

  ## (3) how the iNEXT SE scales with the estimate
  Scale_diag <- Md %>%
    filter(response > 0, !is.na(response_se), !is.na(dist_farm), dist_farm < dist_threshold) %>%
    group_by(Hill) %>%
    group_modify(function(g, ...) {
      m <- stats::lm(log(response_se) ~ log(response), data = g)
      ci <- stats::confint(m)["log(response)", ]
      tibble(n = nrow(g), se_scaling_exponent = coef(m)[["log(response)"]],
             exp_ci_low = ci[[1]], exp_ci_high = ci[[2]],
             cor_se_vs_estimate = cor(g$response_se, g$response),
             cor_selog_vs_estimate = cor(g$response_se / g$response, g$response),
             cv_se = sd(g$response_se) / mean(g$response_se),
             cv_selog = sd(g$response_se / g$response) / mean(g$response_se / g$response))
    }) %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  write_csv(Scale_diag, "Derived/Excels/Response_distribution_scale_diag.csv")

  ## figure: posterior-predictive density by candidate (baseline, climate spec)
  ppc_dens_df <- rd_grid %>%
    filter(index == "baseline", spec == "climate") %>%
    mutate(scale = rd_scale_of(family)) %>%
    pmap_dfr(function(hill, family, key, scale, ...) {
      fit <- rd_fits[[key]]
      yrep <- to_raw(posterior_predict(fit, ndraws = 60), scale)
      y_obs <- obs_raw(fit, scale)
      bind_rows(tibble(hill = hill, family = family, kind = "rep", draw = as.vector(row(yrep)), value = as.vector(yrep)),
                tibble(hill = hill, family = family, kind = "obs", draw = 0L, value = y_obs))
    }) %>%
    mutate(Hill = factor(recode(hill, !!!hill_lab), levels = unname(hill_lab)),
           Family = factor(recode(family, !!!fam_lab), levels = unname(fam_lab)))
  p_ppc <- ggplot(ppc_dens_df %>% filter(kind == "rep"), aes(value, group = draw)) +
    geom_density(colour = "#9ecae1", alpha = 0.5, linewidth = 0.25) +
    geom_density(data = ppc_dens_df %>% filter(kind == "obs"), aes(value), inherit.aes = FALSE,
                 colour = "#08306b", linewidth = 0.9) +
    facet_grid(Hill ~ Family, scales = "free") + coord_cartesian(xlim = c(0, NA)) +
    labs(x = "Diversity estimate (effective no. of species)", y = "Density",
         title = "Posterior-predictive density by candidate response distribution",
         subtitle = "Baseline model, climate + canopy adjustment; dark = observed, light = 60 predictive draws (raw scale)") +
    theme_minimal_grid(11) + theme(strip.text = element_text(size = 9))
  ggsave("Figures/Response_distribution_ppc.png", p_ppc, bg = "white", width = 14, height = 8, dpi = 150)
  print(p_ppc)

  ## figure: the scale diagnostic
  diag_df <- Md %>%
    filter(response > 0, !is.na(response_se), !is.na(dist_farm), dist_farm < dist_threshold) %>%
    mutate(Hill = factor(recode(Hill, !!!hill_lab), levels = unname(hill_lab)), se_log = response_se / response)
  exp_lab <- Scale_diag %>%
    mutate(Hill = factor(recode(Hill, !!!hill_lab), levels = unname(hill_lab)),
           lab = sprintf("log-log slope = %.2f [%.2f, %.2f]", se_scaling_exponent, exp_ci_low, exp_ci_high))
  p_se_raw <- ggplot(diag_df, aes(response, response_se)) +
    geom_point(alpha = 0.5, size = 1.4) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#08519c") +
    facet_wrap(~ Hill, scales = "free") +
    labs(x = "Diversity estimate", y = "iNEXT bootstrap SE", title = "(A) The bootstrap SE increases with the estimate (raw scale)")
  p_se_log <- ggplot(diag_df, aes(response, se_log)) +
    geom_point(alpha = 0.5, size = 1.4) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#08519c") +
    geom_text(data = exp_lab, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 3) +
    facet_wrap(~ Hill, scales = "free") +
    labs(x = "Diversity estimate", y = "SE of log(estimate) = SE / estimate",
         title = "(B) SE / estimate falls as the estimate rises: dividing by the estimate over-corrects")
  resid_keys <- c("Gaussian / log"  = "respdist__richness__climate__baseline__gaussian_log",
                  "Gaussian / sqrt" = "respdist__richness__climate__baseline__gaussian_sqrt",
                  "Gaussian / raw"  = "respdist__richness__climate__baseline__gaussian_raw")
  resid_df <- imap_dfr(resid_keys, function(key, lab) {
    fit <- rd_fits[[key]]
    tibble(candidate = lab, fitted = fitted(fit)[, "Estimate"],
           resid_pearson = residuals(fit, type = "pearson")[, "Estimate"])
  }) %>% mutate(candidate = factor(candidate, levels = names(resid_keys)))
  p_resid <- ggplot(resid_df, aes(fitted, resid_pearson)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(alpha = 0.6, size = 1.4) + geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#08519c") +
    facet_wrap(~ candidate, scales = "free_x") +
    labs(x = "Fitted value (model scale)", y = "Pearson residual", title = "(C) Pearson residuals vs fitted -- richness baseline, climate spec")
  p_scale_diag <- plot_grid(p_se_raw, p_se_log, p_resid, ncol = 1)
  ggsave("Figures/Response_distribution_scale_diag.png", p_scale_diag, bg = "white", width = 10, height = 11, dpi = 150)
  print(p_scale_diag)

  cat("\n[response_dist] LOO best candidate per cell:\n")
  Loo_table %>% group_by(hill, spec, index) %>% slice_max(elpd_loo, n = 1) %>% ungroup() %>%
    count(family, name = "n_cells_won") %>% arrange(desc(n_cells_won)) %>% print()
  cat("[response_dist] SE-scaling exponent (0 = raw ideal, 1 = log ideal, ~0.5 = sqrt):\n")
  print(Scale_diag %>% select(Hill, n, se_scaling_exponent, exp_ci_low, exp_ci_high))
}


# ================================================================== #
# Section: ideal_adj  (was Scripts/04h)                               #
# ================================================================== #
### The DAG-ideal adjustment: focal + elev + elev^2 + precip + precip^2 + pool_point + canopy + sampling + REs -- all three collinear environmental axes at once (04d deliberately avoided this). Run it anyway + assess the collinearity properly (predictor VIF, joint posterior correlations, focal stability vs the 04a specs, prior sensitivity). Outcome: environmental block badly collinear, but the management coefficient (VIF ~ 1.5) is clean and prior-stable.

if ("ideal_adj" %in% run_sections) {
  message("\n=== Section: ideal_adj ===")

  ia_iter <- 3000; ia_warmup <- 1000; ia_adapt <- 0.99

  Md <- Base_frame %>%
    left_join(Pool_raw %>% select(Id_gcs, pool_point), by = "Id_gcs") %>%
    mutate(pool_point_z = as.numeric(scale(pool_point)))

  ia_env  <- c("Elev_mean_z", "I(Elev_mean_z^2)", "Tot_prec_mean_z", "I(Tot_prec_mean_z^2)", "pool_point_z", "canopy_10k_z")
  ia_samp <- c("Num_pc_log_z", "Num_hab_z", "doy_sin", "doy_cos")
  ia_re   <- c("(1 | Id_gcs)", "(1 | CollectorXyear)")

  ia_frame <- function(row) {
    df <- Md %>%
      filter(Hill == row$hill, !is.na(se_log), !is.na(canopy_10k_z), !is.na(pool_point_z),
             !is.na(Num_pc_log_z), !is.na(Num_hab_z), !is.na(Tot_prec_mean_z), !is.na(Elev_mean_z),
             !is.na(dist_farm), dist_farm < dist_threshold)
    if (row$index != "baseline") df <- df %>% mutate(focal_z = .data[[paste0(row$index, "_z")]]) %>% filter(!is.na(focal_z))
    df
  }
  ia_bf <- function(row) mgmt_bf("log_response | resp_se(se_log, sigma = TRUE)",
                                 if (row$index != "baseline") "focal_z" else NULL,
                                 c(ia_env, ia_samp, ia_re))

  ## (1) design-matrix collinearity
  design_df <- ia_frame(list(hill = "richness", index = "Water_mgmt_div")) %>%
    transmute(focal_z, Elev_z = Elev_mean_z, Elev_z2 = Elev_mean_z^2,
              prec_z = Tot_prec_mean_z, prec_z2 = Tot_prec_mean_z^2,
              pool_z = pool_point_z, canopy_z = canopy_10k_z,
              Num_pc_z = Num_pc_log_z, Num_hab_z, doy_sin, doy_cos)
  Cor_pred <- cor(design_df, use = "pairwise")
  vif_manual <- function(dat) map_dbl(names(dat), function(v) {
    r2 <- summary(lm(reformulate(setdiff(names(dat), v), response = v), data = dat))$r.squared
    1 / (1 - r2)
  }) %>% set_names(names(dat))
  VIF <- vif_manual(design_df)
  Collinearity <- tibble(predictor = names(VIF), VIF = round(VIF, 2)) %>% arrange(desc(VIF))
  high_cor_pairs <- which(upper.tri(Cor_pred) & abs(Cor_pred) > 0.6, arr.ind = TRUE) %>%
    as_tibble() %>%
    mutate(pair = sprintf("%s ~ %s", rownames(Cor_pred)[row], colnames(Cor_pred)[col]),
           r = round(Cor_pred[cbind(row, col)], 2)) %>%
    select(pair, r) %>% arrange(desc(abs(r)))
  write_csv(Collinearity, "Derived/Excels/Ideal_adjustment_collinearity.csv")

  ## (2) fit the grid
  ia_grid <- expand_grid(hill = c("richness", "shannon", "simpson"), index = c("baseline", div_indices)) %>%
    arrange(hill, index) %>%
    mutate(key = paste("ideal", hill, index, sep = "__"),
           structure = if_else(index == "baseline", "baseline", "index"))

  ia_fits <- fit_model_grid(ia_grid, build_bf = ia_bf, build_data = ia_frame,
                            iter = ia_iter, warmup = ia_warmup, adapt_delta = ia_adapt)

  Ideal_summaries <- pmap_dfr(ia_grid, function(hill, index, key, ...) {
    fit <- ia_fits[[key]]
    draws <- as_draws_matrix(fit, variable = "^b_", regex = TRUE)
    focal <- if ("b_focal_z" %in% colnames(draws)) draws[, "b_focal_z"] else NULL
    cc <- cor(draws[, setdiff(colnames(draws), "b_Intercept"), drop = FALSE])
    cc[!upper.tri(cc)] <- NA
    worst_ix <- which.max(abs(cc))
    focal_cor <- if (is.null(focal)) NA_real_ else {
      fc <- abs(cc["b_focal_z", ]); fc <- fc[!is.na(fc) & names(fc) != "b_focal_z"]
      if (length(fc)) round(max(fc), 2) else NA_real_
    }
    tibble(hill, index, n_obs = nobs(fit),
           focal_est = if (is.null(focal)) NA_real_ else round(median(focal), 4),
           focal_lo  = if (is.null(focal)) NA_real_ else round(quantile(focal, 0.05), 4),
           focal_hi  = if (is.null(focal)) NA_real_ else round(quantile(focal, 0.95), 4),
           p_direction_pos = if (is.null(focal)) NA_real_ else round(mean(focal > 0), 3),
           bayes_R2 = round(bayes_R2(fit)[, "Estimate"], 4),
           min_bulk_ess = round(min(summary(fit)$fixed[, "Bulk_ESS"], na.rm = TRUE)),
           max_post_cor = round(cc[worst_ix], 2),
           max_post_cor_pair = paste(rownames(cc)[row(cc)[worst_ix]], colnames(cc)[col(cc)[worst_ix]], sep = " ~ "),
           focal_max_post_cor = focal_cor,
           max_rhat = round(max(rhat(fit), na.rm = TRUE), 3),
           n_divergent = sum(subset(nuts_params(fit), Parameter == "divergent__")$Value))
  })
  write_csv(Ideal_summaries, "Derived/Excels/Ideal_adjustment_summaries.csv")

  ## (3) focal stability: ideal vs 04a climate / ecoregion
  Spec_compare <- read_csv("Derived/Excels/Farm_mgmt_model_summaries.csv", show_col_types = FALSE) %>%
    filter(data_subset == "primary", term == "focal_z", spec %in% c("climate", "ecoregion")) %>%
    transmute(hill, index, spec, focal_est = estimate, focal_lo = conf_low, focal_hi = conf_high) %>%
    bind_rows(Ideal_summaries %>% filter(index != "baseline") %>%
                transmute(hill, index, spec = "ideal (climate^2 + pool + canopy)", focal_est, focal_lo, focal_hi)) %>%
    arrange(hill, index, spec)

  ## (4) prior sensitivity on b (pasture / water)
  ps_grid <- expand_grid(hill = c("richness", "shannon"),
                         index = c("Pasture_mgmt_div", "Water_mgmt_div"),
                         prior_sd = c(0.5, 0.75, 1.5)) %>%
    mutate(key = paste("idealPS", hill, index, prior_sd, sep = "__"), structure = key)
  ps_fits <- fit_model_grid(
    ps_grid, build_bf = ia_bf, build_data = ia_frame,
    prior = function(row, data) c(prior(student_t(3, 3, 2.5), class = "Intercept"),
                                  set_prior(sprintf("normal(0, %.2f)", row$prior_sd), class = "b"),
                                  prior(exponential(1), class = "sd"), prior(exponential(1), class = "sigma")),
    iter = ia_iter, warmup = ia_warmup, adapt_delta = ia_adapt
  )
  Prior_sens <- pmap_dfr(ps_grid, function(hill, index, prior_sd, key, ...) {
    b <- fixef(ps_fits[[key]])["focal_z", ]
    tibble(hill, index, prior_sd,
           focal = sprintf("%+.3f [%+.3f, %+.3f]", b["Estimate"], b["Q2.5"], b["Q97.5"]),
           focal_est = round(b["Estimate"], 3))
  })
  write_csv(Prior_sens, "Derived/Excels/Ideal_adjustment_prior_sens.csv")

  ## figure
  p_effects <- Spec_compare %>%
    mutate(Index = factor(recode(index, !!!index_lab), levels = rev(unname(index_lab))),
           Hill = factor(recode(hill, !!!hill_lab), levels = unname(hill_lab)),
           Spec = factor(recode(spec,
                                `ideal (climate^2 + pool + canopy)` = "ideal: climate² + pool + canopy",
                                climate = "04a climate", ecoregion = "04a Ecoregion"),
                         levels = c("ideal: climate² + pool + canopy", "04a climate", "04a Ecoregion"))) %>%
    ggplot(aes(focal_est, Index, colour = Spec)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = focal_lo, xmax = focal_hi), position = position_dodge(width = 0.6), size = 0.4) +
    scale_colour_manual(values = c("ideal: climate² + pool + canopy" = "#08519c",
                                   "04a climate" = "#d95f02", "04a Ecoregion" = "#238b45"), name = NULL) +
    facet_wrap(~ Hill) +
    labs(x = "Management-diversification coefficient (focal_z, 90% CrI)", y = NULL,
         title = "Management effect under the DAG-ideal adjustment vs the Scripts/04a specs",
         subtitle = "Response: abundance / Cmax, primary subset (< 300 m)") +
    theme(legend.position = "bottom")
  ggsave("Figures/Ideal_adjustment_effects.png", p_effects, bg = "white", width = 11, height = 6, dpi = 150)
  print(p_effects)

  cat("\n[ideal_adj] VIF (design matrix incl. quadratics):\n"); print(Collinearity, n = Inf)
  cat("[ideal_adj] pairs |r| > 0.6:\n"); print(high_cor_pairs, n = Inf)
  cat("[ideal_adj] focal coefficient + worst joint posterior correlation:\n")
  Ideal_summaries %>% select(hill, index, focal_est, focal_lo, focal_hi, bayes_R2, focal_max_post_cor, max_rhat, n_divergent) %>%
    print(n = Inf, width = Inf)
}


# ================================================================== #
# Section: endemism_pool  (was Scripts/04i)                           #
# ================================================================== #
### Swap the raw pool count for the range-rarity-weighted pool_cwe (log then z-scored) in three specs that vary whether it competes with the climate axes. Reports the management coefficient, the pool coefficient + its sign, fit, collinearity. Outcome: no change to the management conclusion; on the log scale pool_cwe is r ~ 0.85 with elevation -- a collinear substitute, not a new axis.

if ("endemism_pool" %in% run_sections) {
  message("\n=== Section: endemism_pool ===")

  ep_iter <- 3000; ep_warmup <- 1000; ep_adapt <- 0.99
  pool_var <- "pool_cwe"

  Md <- Base_frame %>%
    left_join(Pool_raw %>% select(Id_gcs, pool_point, pool_we, pool_cwe, pool_wes), by = "Id_gcs") %>%
    mutate(pool_z = as.numeric(scale(log(.data[[pool_var]]))))

  ep_sampling <- c("Num_pc_log_z", "Num_hab_z", "doy_sin", "doy_cos", "(1 | Id_gcs)", "(1 | CollectorXyear)")
  ep_specs <- c(
    swap      = "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + pool_z + canopy_10k_z",
    we_precip = "pool_z + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
    we_canopy = "pool_z + canopy_10k_z"
  )

  ep_frame <- function(row) {
    df <- Md %>%
      filter(Hill == row$hill, !is.na(se_log), !is.na(canopy_10k_z), !is.na(pool_z),
             !is.na(Num_pc_log_z), !is.na(Num_hab_z), !is.na(Tot_prec_mean_z), !is.na(Elev_mean_z),
             !is.na(dist_farm), dist_farm < dist_threshold)
    if (row$index != "baseline") df <- df %>% mutate(focal_z = .data[[paste0(row$index, "_z")]]) %>% filter(!is.na(focal_z))
    df
  }

  ep_grid <- expand_grid(hill = c("richness", "shannon", "simpson"),
                         index = c("baseline", div_indices), spec = names(ep_specs)) %>%
    arrange(spec, hill, index) %>%
    mutate(key = paste("poolwe", hill, index, spec, sep = "__"),
           structure = paste(spec, if_else(index == "baseline", "baseline", "index"), sep = "__"))

  ep_fits <- fit_model_grid(
    ep_grid,
    build_bf = function(row) mgmt_bf("log_response | resp_se(se_log, sigma = TRUE)",
                                     if (row$index != "baseline") "focal_z" else NULL,
                                     c(ep_specs[[row$spec]], ep_sampling)),
    build_data = ep_frame,
    iter = ep_iter, warmup = ep_warmup, adapt_delta = ep_adapt
  )

  coef_ci <- function(draws, v) {
    if (!v %in% colnames(draws)) return(tibble(est = NA, lo = NA, hi = NA, pd_pos = NA))
    x <- draws[, v]
    tibble(est = median(x), lo = quantile(x, 0.05), hi = quantile(x, 0.95), pd_pos = mean(x > 0))
  }
  Pool_summaries <- pmap_dfr(ep_grid, function(hill, index, spec, key, ...) {
    fit <- ep_fits[[key]]
    draws <- as_draws_matrix(fit, variable = "^b_", regex = TRUE)
    slopes <- setdiff(colnames(draws), "b_Intercept")
    cc <- cor(draws[, slopes, drop = FALSE]); diag(cc) <- 0
    we_cor <- if ("b_pool_z" %in% slopes) round(max(abs(cc["b_pool_z", ])), 2) else NA_real_
    focal_cor <- if ("b_focal_z" %in% slopes) round(max(abs(cc["b_focal_z", ])), 2) else NA_real_
    fo <- coef_ci(draws, "b_focal_z"); we <- coef_ci(draws, "b_pool_z")
    tibble(hill, index, spec, n_obs = nobs(fit),
           focal_est = round(fo$est, 4), focal_lo = round(fo$lo, 4), focal_hi = round(fo$hi, 4),
           focal_pd_pos = round(fo$pd_pos, 3),
           pool_est = round(we$est, 4), pool_lo = round(we$lo, 4), pool_hi = round(we$hi, 4),
           pool_pd_pos = round(we$pd_pos, 3),
           bayes_R2 = round(bayes_R2(fit)[, "Estimate"], 4),
           pool_max_post_cor = we_cor, focal_max_post_cor = focal_cor,
           min_bulk_ess = round(min(summary(fit)$fixed[, "Bulk_ESS"], na.rm = TRUE)),
           max_rhat = round(max(rhat(fit), na.rm = TRUE), 3),
           n_divergent = sum(subset(nuts_params(fit), Parameter == "divergent__")$Value))
  })
  write_csv(Pool_summaries, "Derived/Excels/Pool_endemism_summaries.csv")

  ref <- read_csv("Derived/Excels/Farm_mgmt_model_summaries.csv", show_col_types = FALSE) %>%
    filter(data_subset == "primary", term == "focal_z", spec == "climate") %>%
    transmute(hill, index, spec = "04a climate (no pool)", focal_est = estimate, focal_lo = conf_low, focal_hi = conf_high)
  ref_ideal <- tryCatch(
    read_csv("Derived/Excels/Ideal_adjustment_summaries.csv", show_col_types = FALSE) %>%
      filter(index != "baseline") %>%
      transmute(hill, index, spec = "ideal_adj (pool_point)", focal_est, focal_lo, focal_hi),
    error = function(e) tibble())
  Comparison <- Pool_summaries %>%
    filter(index != "baseline") %>%
    transmute(hill, index, spec = recode(spec, swap = "pool_cwe swap (ideal spec)",
                                         we_precip = "pool_cwe + precip", we_canopy = "pool_cwe + canopy"),
              focal_est, focal_lo, focal_hi) %>%
    bind_rows(ref, ref_ideal) %>%
    arrange(hill, index, spec)
  write_csv(Comparison, "Derived/Excels/Pool_endemism_comparison.csv")

  p_effects <- Comparison %>%
    mutate(Index = factor(recode(index, !!!index_lab), levels = rev(unname(index_lab))),
           Hill = factor(recode(hill, !!!hill_lab), levels = unname(hill_lab)),
           Spec = factor(spec, levels = c("pool_cwe swap (ideal spec)", "pool_cwe + precip", "pool_cwe + canopy",
                                          "ideal_adj (pool_point)", "04a climate (no pool)"))) %>%
    ggplot(aes(focal_est, Index, colour = Spec)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = focal_lo, xmax = focal_hi), position = position_dodge(width = 0.7), size = 0.35) +
    scale_colour_brewer(palette = "Set1", name = NULL) +
    facet_wrap(~ Hill) +
    labs(x = "Management coefficient (focal_z, 90% CrI)", y = NULL,
         title = "Does the range-rarity-weighted pool change the management estimate?",
         subtitle = "Response: abundance / Cmax, primary subset (< 300 m)") +
    theme(legend.position = "bottom")
  ggsave("Figures/Pool_endemism_effects.png", p_effects, bg = "white", width = 12, height = 6, dpi = 150)
  print(p_effects)

  cat("\n[endemism_pool] management coefficient: pool_cwe specs vs ideal_adj (pool_point) vs 04a climate\n")
  Comparison %>%
    mutate(cell = sprintf("%+.3f [%+.3f, %+.3f]", focal_est, focal_lo, focal_hi)) %>%
    select(hill, index, spec, cell) %>%
    pivot_wider(names_from = spec, values_from = cell) %>% print(n = Inf, width = Inf)
  cat("[endemism_pool] baseline bayes_R2 by spec\n")
  Pool_summaries %>% filter(index == "baseline") %>% select(hill, spec, bayes_R2) %>%
    pivot_wider(names_from = spec, values_from = bayes_R2) %>% print()
}


# ================================================================== #
# Section: spline_env  (was Scripts/04j)                              #
# ================================================================== #
### Replace the elevation / precipitation quadratics in the primary `climate` spec with rigid k = 4 penalized splines (both responses, baseline + 4 indices, primary subset). Compare quadratic vs spline by focal coefficient, LOO, bayes_R2, PPC, fitted env shapes, and mgcv edf. Outcome: no substantive change; LOO / bayes_R2 mildly favour the quadratic -- no evidence it is too restrictive.

if ("spline_env" %in% run_sections) {
  message("\n=== Section: spline_env ===")

  se_iter <- 4000; se_warmup <- 1500; se_adapt <- 0.999

  env_forms <- tribble(
    ~form,       ~env_terms,
    "quadratic", "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z + pool_wes_z",
    "spline",    "s(Elev_mean_z, k = 4) + s(Tot_prec_mean_z, k = 4) + canopy_10k_z + pool_wes_z"
  )

  Predictors <- Base_frame %>%
    distinct(Assemblage, Id_gcs, CollectorXyear, Ecoregion, dist_farm, doy_sin, doy_cos,
             Elev_mean_z, Tot_prec_mean_z, canopy_10k_z, pool_wes_z, Num_pc_log_z,
             Land_use_div_z, Water_mgmt_div_z, Pasture_mgmt_div_z, All_practices_div_z)
  Cmax_resp <- Base_frame %>% transmute(Assemblage, Hill, response_arm = "cmax", log_response, se_log)
  Inc_resp <- read_csv(latest_file("Derived/Excels", "^Tax_div_incidence_.*\\.csv$"), show_col_types = FALSE) %>%
    transmute(Assemblage, Hill, response_arm = "incidence", log_response = log(qD), se_log = qD_se / qD)
  SE_data <- bind_rows(Cmax_resp, Inc_resp) %>%
    inner_join(Predictors, by = "Assemblage") %>%
    filter(is.finite(se_log), se_log > 0)

  se_frame <- function(row) {
    df <- SE_data %>%
      filter(response_arm == row$response_arm, Hill == row$hill,
             !is.na(doy_sin), !is.na(Elev_mean_z), !is.na(Tot_prec_mean_z),
             !is.na(canopy_10k_z), !is.na(pool_wes_z), !is.na(dist_farm), dist_farm < dist_threshold)
    if (row$response_arm == "cmax") df <- df %>% filter(!is.na(Num_pc_log_z))
    if (row$index != "baseline") df <- df %>% mutate(focal_z = .data[[paste0(row$index, "_z")]]) %>% filter(!is.na(focal_z))
    df
  }
  se_bf <- function(row) mgmt_bf(
    "log_response | resp_se(se_log, sigma = TRUE)",
    if (row$index != "baseline") "focal_z" else NULL,
    c(env_forms$env_terms[match(row$form, env_forms$form)],
      if (row$response_arm == "cmax") "Num_pc_log_z",
      "doy_sin", "doy_cos", "(1 | Id_gcs)", "(1 | CollectorXyear)")
  )

  se_grid <- expand_grid(response_arm = c("cmax", "incidence"), env_forms,
                         hill = c("richness", "shannon", "simpson"),
                         index = c("baseline", div_indices)) %>%
    mutate(structure = paste(response_arm, form, if_else(index == "baseline", "baseline", "index"), sep = "__"),
           key = paste("spl", response_arm, form, hill, index, sep = "__")) %>%
    arrange(response_arm, form, structure, hill, index)

  se_fits <- fit_model_grid(se_grid, build_bf = se_bf, build_data = se_frame,
                            iter = se_iter, warmup = se_warmup, adapt_delta = se_adapt)
  se_loo <- imap(se_fits, ~ suppressWarnings(loo(.x)))

  smooth_sds <- function(fit) {
    v <- brms::variables(fit); sds <- v[startsWith(v, "sds_")]
    if (!length(sds)) return(tibble(sds_elev = NA_real_, sds_prec = NA_real_))
    med <- apply(as_draws_matrix(fit, variable = sds), 2, median)
    tibble(sds_elev = unname(med[str_detect(names(med), "Elev")][1]),
           sds_prec = unname(med[str_detect(names(med), "prec|Prec")][1]))
  }

  SE_summaries <- pmap(se_grid, function(response_arm, form, hill, index, key, ...) {
    fit <- se_fits[[key]]
    has_focal <- "b_focal_z" %in% brms::variables(fit)
    focal <- if (has_focal) as.numeric(as_draws_matrix(fit, variable = "b_focal_z")) else NA_real_
    lo <- se_loo[[key]]
    bind_cols(
      tibble(response_arm, form, hill, index, n_obs = nobs(fit),
             focal_est = if (has_focal) median(focal) else NA_real_,
             focal_lo  = if (has_focal) quantile(focal, 0.05) else NA_real_,
             focal_hi  = if (has_focal) quantile(focal, 0.95) else NA_real_,
             pd_pos    = if (has_focal) mean(focal > 0) else NA_real_,
             bayes_R2  = bayes_R2(fit)[, "Estimate"],
             elpd_loo  = lo$estimates["elpd_loo", "Estimate"],
             p_loo     = lo$estimates["p_loo", "Estimate"],
             n_pareto_bad = sum(lo$diagnostics$pareto_k > 0.7),
             max_rhat  = max(rhat(fit), na.rm = TRUE),
             n_divergent = sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)),
      smooth_sds(fit))
  }) %>%
    list_rbind() %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  write_csv(SE_summaries, "Derived/Excels/Spline_env_check_summaries.csv")

  Focal_compare <- SE_summaries %>%
    filter(index != "baseline") %>%
    select(response_arm, hill, index, form, focal_est, focal_lo, focal_hi, pd_pos) %>%
    pivot_wider(names_from = form, values_from = c(focal_est, focal_lo, focal_hi, pd_pos)) %>%
    mutate(delta_est = round(focal_est_spline - focal_est_quadratic, 4),
           quad_excl0 = sign(focal_lo_quadratic) == sign(focal_hi_quadratic),
           spline_excl0 = sign(focal_lo_spline) == sign(focal_hi_spline),
           crI_conclusion_changes = quad_excl0 != spline_excl0) %>%
    arrange(response_arm, hill, index)
  write_csv(Focal_compare, "Derived/Excels/Spline_env_check_focal_compare.csv")

  Loo_compare <- se_grid %>%
    distinct(response_arm, hill, index) %>%
    pmap(function(response_arm, hill, index) {
      kq <- paste("spl", response_arm, "quadratic", hill, index, sep = "__")
      ks <- paste("spl", response_arm, "spline", hill, index, sep = "__")
      cmp <- suppressWarnings(loo_compare(list(quadratic = se_loo[[kq]], spline = se_loo[[ks]])))
      d_spline <- se_loo[[ks]]$estimates["elpd_loo", "Estimate"] - se_loo[[kq]]$estimates["elpd_loo", "Estimate"]
      se_diff <- cmp[2, "se_diff"]
      tibble(response_arm, hill, index,
             elpd_diff_spline_minus_quad = round(d_spline, 3), se_diff = round(se_diff, 3),
             favours = case_when(d_spline >  2 * se_diff ~ "spline",
                                 d_spline < -2 * se_diff ~ "quadratic", TRUE ~ "negligible"))
    }) %>%
    list_rbind() %>%
    arrange(response_arm, hill, index)
  write_csv(Loo_compare, "Derived/Excels/Spline_env_check_loo.csv")

  ## approximate frequentist edf of each env smooth (mgcv, unweighted)
  edf_one <- function(response_arm, hill) {
    df <- se_frame(list(response_arm = response_arm, hill = hill, index = "baseline")) %>%
      mutate(Id_gcs = factor(Id_gcs), CollectorXyear = factor(CollectorXyear))
    rhs <- c("s(Elev_mean_z, k = 4)", "s(Tot_prec_mean_z, k = 4)", "canopy_10k_z", "pool_wes_z",
             if (response_arm == "cmax") "Num_pc_log_z", "doy_sin", "doy_cos",
             "s(Id_gcs, bs = 're')", "s(CollectorXyear, bs = 're')")
    m <- mgcv::gam(as.formula(paste("log_response ~", paste(rhs, collapse = " + "))), data = df, method = "REML")
    st <- summary(m)$s.table
    tibble(response_arm, hill,
           edf_elev = round(st["s(Elev_mean_z)", "edf"], 2), p_elev = signif(st["s(Elev_mean_z)", "p-value"], 2),
           edf_prec = round(st["s(Tot_prec_mean_z)", "edf"], 2), p_prec = signif(st["s(Tot_prec_mean_z)", "p-value"], 2))
  }
  Edf_tbl <- expand_grid(response_arm = c("cmax", "incidence"), hill = c("richness", "shannon", "simpson")) %>%
    pmap(edf_one) %>% list_rbind()
  write_csv(Edf_tbl, "Derived/Excels/Spline_env_check_edf.csv")

  arm_lab <- c(cmax = "Coverage / Cmax", incidence = "Incidence (m* = 6)")

  p_focal <- SE_summaries %>%
    filter(index != "baseline") %>%
    mutate(Index = factor(recode(index, !!!index_lab), levels = rev(unname(index_lab))),
           Hill  = factor(recode(hill, !!!hill_lab), levels = unname(hill_lab)),
           Arm   = factor(recode(response_arm, !!!arm_lab), levels = unname(arm_lab)),
           Form  = str_to_title(form)) %>%
    ggplot(aes(focal_est, Index, colour = Form)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = focal_lo, xmax = focal_hi), position = position_dodge(width = 0.5), size = 0.4) +
    scale_colour_manual(values = c(Quadratic = "#d95f02", Spline = "#1b9e77"), name = "Elev / precip form") +
    facet_grid(Arm ~ Hill) +
    labs(x = "Management-diversification coefficient (focal_z, 90% CrI)", y = NULL,
         title = "Focal management effect: quadratic vs rigid-spline elevation / precipitation",
         subtitle = "Primary climate spec, primary distance subset.") +
    theme(legend.position = "bottom")
  ggsave("Figures/Spline_env_check_focal.png", p_focal, bg = "white", width = 11, height = 7, dpi = 150)
  print(p_focal)

  ce_curve <- function(key, var) {
    fit <- se_fits[[key]]
    ce <- conditional_effects(fit, effects = var, resolution = 120)[[var]]
    meta <- se_grid %>% filter(key == !!key)
    tibble(x = ce[[var]], est = ce$estimate__,
           variable = if (str_detect(var, "Elev")) "Elevation (z)" else "Precipitation (z)",
           response_arm = meta$response_arm, form = meta$form, hill = meta$hill) %>%
      mutate(est_c = est - mean(est), lo_c = ce$lower__ - mean(est), hi_c = ce$upper__ - mean(est))
  }
  Shape_df <- se_grid %>%
    filter(index == "baseline") %>%
    pull(key) %>%
    map(~ bind_rows(ce_curve(.x, "Elev_mean_z"), ce_curve(.x, "Tot_prec_mean_z"))) %>%
    list_rbind() %>%
    mutate(Hill = factor(recode(hill, !!!hill_lab), levels = unname(hill_lab)),
           Arm  = recode(response_arm, !!!arm_lab), Form = str_to_title(form))
  p_shapes <- ggplot(Shape_df, aes(x, est_c, colour = Form, fill = Form)) +
    geom_ribbon(aes(ymin = lo_c, ymax = hi_c), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = Arm), linewidth = 0.7) +
    scale_colour_manual(values = c(Quadratic = "#d95f02", Spline = "#1b9e77"), name = "Form") +
    scale_fill_manual(values = c(Quadratic = "#d95f02", Spline = "#1b9e77"), name = "Form") +
    scale_linetype_manual(values = c(`Coverage / Cmax` = "solid", `Incidence (m* = 6)` = "22"), name = "Response") +
    facet_grid(Hill ~ variable, scales = "free_x") +
    labs(x = "Standardised predictor", y = "Centred partial effect on log(diversity)",
         title = "Fitted elevation / precipitation response: quadratic vs rigid spline (k = 4)") +
    theme(legend.position = "bottom")
  ggsave("Figures/Spline_env_check_shapes.png", p_shapes, bg = "white", width = 10, height = 9, dpi = 150)
  print(p_shapes)

  ppc_panels <- se_grid %>%
    filter(index == "baseline") %>%
    pmap(function(response_arm, form, hill, key, ...) {
      pp_check(se_fits[[key]], ndraws = 60) +
        ggtitle(sprintf("%s | %s | %s", recode(response_arm, !!!arm_lab), str_to_title(form), hill)) +
        theme(legend.position = "none", plot.title = element_text(size = 9))
    })
  p_ppc <- plot_grid(plotlist = ppc_panels, ncol = 3)
  ggsave("Figures/Spline_env_check_ppc.png", p_ppc, bg = "white", width = 13, height = 12, dpi = 150)
  print(p_ppc)

  cat("\n[spline_env] approximate env-smooth edf (mgcv; 1 = linear, ~2 = quadratic-equivalent, 3 = max):\n")
  print(Edf_tbl, n = Inf)
  cat("[spline_env] LOO spline - quadratic (favours 'spline' only if diff > 2 SE):\n")
  print(Loo_compare, n = Inf)
  cat("[spline_env] largest |change| in the focal coefficient:",
      round(max(abs(Focal_compare$delta_est), na.rm = TRUE), 4), "SD;",
      sum(Focal_compare$crI_conclusion_changes, na.rm = TRUE), "of", nrow(Focal_compare), "CrI zero-exclusion flips\n")
}

message("\n=== 04b done: ", paste(run_sections, collapse = ", "), " ===")
