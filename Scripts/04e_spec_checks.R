# Sensitivity of the farm-management model to specification choices ----

### (a) the (1 | CollectorXyear) random effect -- CollectorXyear (dataset x year-group, 8 levels) is partly confounded with Ecoregion: 4 of 8 levels are 100% Piedemonte, though the 3 largest datasets (Cipav, Gaica 2013/14, Gaica 2016/17; 72 of 108 assemblages) span 4-5 ecoregions each (Cramer's V = 0.38). So the RE absorbs some region structure but is not equivalent to (1 | Ecoregion). Does dropping it change the management coefficients on the full 5-ecoregion data?
### (b) resp_se() -- the measurement-error response term. What happens if the iNEXT point estimates are used as-is, with no per-assemblage SE?
### (c) response scale / likelihood -- the main models are Gaussian on log(diversity). The raw estimates are close to symmetric and logging mildly over-corrects (see qmd/Farm_mgmt_summary.qmd fig-loglik), so is the near-null result specific to the log scale? Refit pasture/water x richness/shannon (climate, < 300 m) as Gaussian-on-log, Student-t-on-log, and Gaussian-on-raw, all with resp_se, and compare the focal effect expressed as a % change in diversity per SD of the index.

### Refits the climate and ecoregion primary models (< 300 m) for each index x {richness, shannon} under three variants: base / no CollectorXyear RE / no resp_se. Reads Scripts/04_farm_mgmt_mod.R's persisted frame. Outputs Derived/Excels/Spec_checks.csv and Derived/Excels/Likelihood_sensitivity.csv.

### OUTCOME (2026-08-29):
### (a) CollectorXyear RE sd is 0.18-0.31 (log scale) -- real batch/protocol structure, not negligible. But dropping the RE shifts the focal coefficients by <= 0.022 (mostly <= 0.01): pasture climate +0.080 -> +0.077, water climate +0.074 -> +0.080. So it is NOT acting like (1|Ecoregion) -- if it were absorbing region confounding, dropping it from the climate spec would inflate pasture/water, and it doesn't. It captures field-protocol differences among the 6 datasets, as intended. Keep it. (Also tested in Piedemonte, 04c -- same, <= 0.01.)
### (b) Dropping resp_se gives LARGER |focal| (mean 0.050 vs 0.040, ~25% bigger; pasture richness +0.080 -> +0.106) and WIDER CrIs (mean 0.174 vs 0.162). resp_se turns the fit into a precision-weighted regression -- the well-sampled (low-SE) assemblages get more influence, and they show a weaker management association than the noisy ones. So resp_se pulls the coefficients down and tightens them by extracting a cleaner signal from the precise points. It is both the honest choice (the iNEXT SEs are real) and the conservative one -- without it, pasture/water richness CrIs would just barely exclude zero.

# Setup ----
library(tidyverse)
library(brms)
source("Scripts/Farm_diversity_fns.R")
options(mc.cores = 4, brms.backend = "rstan")
dir.create("Derived/models", recursive = TRUE, showWarnings = FALSE)
set.seed(1989)

chains <- 4; iter <- 3000; warmup <- 1000; adapt_delta <- 0.99
div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")
dist_threshold <- 300
mod_priors <- c(
  prior(student_t(3, 3, 2.5), class = "Intercept"),
  prior(normal(0, 0.75), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

Model_data <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs),
         doy_sin = sin(2 * pi * doy / 365), doy_cos = cos(2 * pi * doy / 365))

env_blocks <- c(
  climate   = "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
  ecoregion = "Ecoregion + canopy_10k_z"
)

## variant -> (response LHS, random-effect block)
variant_lhs <- c(base = "log_response | resp_se(se_log, sigma = TRUE)",
                 no_collector = "log_response | resp_se(se_log, sigma = TRUE)",
                 no_resp_se = "log_response")
variant_re  <- c(base = "(1 | Id_gcs) + (1 | CollectorXyear)",
                 no_collector = "(1 | Id_gcs)",
                 no_resp_se = "(1 | Id_gcs) + (1 | CollectorXyear)")

formula_for <- function(env, variant) {
  bf(as.formula(paste0(
    variant_lhs[[variant]], " ~ focal_z + ", env_blocks[[env]],
    " + Num_pc_log_z + Num_hab_z + doy_sin + doy_cos + ", variant_re[[variant]]
  )))
}

frame_for <- function(hill, index, env) {
  df <- Model_data %>%
    filter(Hill == hill, !is.na(se_log), !is.na(canopy_10k_z), !is.na(Num_pc_log_z),
           !is.na(Num_hab_z), !is.na(dist_farm), dist_farm < dist_threshold)
  if (str_detect(env_blocks[[env]], "Elev|prec")) df <- df %>% filter(!is.na(Elev_mean_z), !is.na(Tot_prec_mean_z))
  df %>% mutate(focal_z = .data[[paste0(index, "_z")]]) %>% filter(!is.na(focal_z))
}

grid <- expand_grid(hill = c("richness", "shannon"), index = div_indices,
                    env = names(env_blocks), variant = names(variant_lhs)) %>%
  mutate(key = paste("speccheck", hill, index, env, variant, sep = "__"),
         structure = paste(env, variant, sep = "__"))

base_by_structure <- list()
fits <- vector("list", nrow(grid)); names(fits) <- grid$key
for (i in seq_len(nrow(grid))) {
  row <- as.list(grid[i, ])
  message(sprintf("[%d/%d] %s", i, nrow(grid), row$key))
  frame <- frame_for(row$hill, row$index, row$env)
  common <- list(chains = chains, iter = iter, warmup = warmup, seed = 1989,
                 control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
                 file = sprintf("Derived/models/%s", row$key), file_refit = "on_change")
  base <- base_by_structure[[row$structure]]
  fit <- if (is.null(base)) {
    do.call(brm, c(list(formula = formula_for(row$env, row$variant), data = frame, prior = mod_priors), common))
  } else {
    do.call(update, c(list(object = base, newdata = frame, recompile = FALSE), common))
  }
  if (is.null(base)) base_by_structure[[row$structure]] <- fit
  fits[[row$key]] <- fit
}

# Summarise ----

summ <- pmap_dfr(grid, function(hill, index, env, variant, key, ...) {
  f <- fits[[key]]
  b <- fixef(f)["focal_z", ]
  ## sd of the CollectorXyear RE, where present
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

cat("\n== (a) CollectorXyear random effect ==\n")
cat("CollectorXyear RE sd in the base models (0 => it explains nothing):\n")
summ %>% filter(variant == "base") %>% distinct(env, hill, index, collectorXyear_sd) %>%
  pivot_wider(names_from = env, values_from = collectorXyear_sd) %>% print(n = Inf)

cat("\nFocal coefficient: base vs dropping (1 | CollectorXyear)\n")
summ %>% filter(variant %in% c("base", "no_collector")) %>%
  select(hill, index, env, variant, focal_est) %>%
  pivot_wider(names_from = variant, values_from = focal_est) %>%
  mutate(shift = round(no_collector - base, 3)) %>%
  arrange(env, hill, index) %>% print(n = Inf)

cat("\n== (b) measurement-error response (resp_se) ==\n")
cat("Focal coefficient: base (resp_se) vs no_resp_se; and sigma\n")
summ %>% filter(variant %in% c("base", "no_resp_se")) %>%
  select(hill, index, env, variant, focal_est, sigma) %>%
  pivot_wider(names_from = variant, values_from = c(focal_est, sigma)) %>%
  mutate(coef_shift = round(focal_est_no_resp_se - focal_est_base, 3)) %>%
  arrange(env, hill, index) %>% print(n = Inf)

# (c) Response scale / likelihood family ----

### Gaussian-on-log (the main model), Student-t-on-log (robust to the mild left skew of log diversity), Gaussian-on-raw (the raw estimates are closer to symmetric). All keep resp_se(., sigma = TRUE); resp_se only supports gaussian + student, which is why Gamma / skew-normal are not options here. Priors follow the response scale: log models keep the 04 priors; the raw model uses response-scale priors. The focal effect is reported as a % change in diversity per +1 SD of the index -- exp(b) - 1 for the log models, b / mean(response) for the raw model -- so the three are directly comparable.

lik_indices <- c("Pasture_mgmt_div", "Water_mgmt_div")

lik_priors_log <- mod_priors   # student_t(3,3,2.5) intercept is on the log scale
lik_priors_raw <- c(
  prior(student_t(3, 20, 15), class = "Intercept"),
  prior(normal(0, 5), class = "b"),
  prior(exponential(0.15), class = "sd"),
  prior(exponential(0.15), class = "sigma")
)

lik_grid <- expand_grid(hill = c("richness", "shannon"), index = lik_indices,
                        lik = c("gaussian_log", "student_log", "gaussian_raw")) %>%
  mutate(key = paste("likcheck", hill, index, lik, sep = "__"))

lik_rhs <- paste("focal_z +", env_blocks[["climate"]],
                 "+ Num_pc_log_z + Num_hab_z + doy_sin + doy_cos + (1 | Id_gcs) + (1 | CollectorXyear)")

lik_summ <- pmap_dfr(lik_grid, function(hill, index, lik, key) {
  df <- frame_for(hill, index, "climate")
  raw <- lik == "gaussian_raw"
  lhs <- if (raw) "response | resp_se(response_se, sigma = TRUE)" else "log_response | resp_se(se_log, sigma = TRUE)"
  fam <- if (lik == "student_log") student() else gaussian()
  f <- brm(bf(as.formula(paste(lhs, "~", lik_rhs))), data = df, family = fam,
           prior = if (raw) lik_priors_raw else lik_priors_log,
           chains = chains, iter = iter, warmup = warmup, seed = 1989,
           control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
           file = sprintf("Derived/models/%s", key), file_refit = "on_change")
  b <- fixef(f)["focal_z", ]
  to_pct <- function(x) if (raw) 100 * x / mean(df$response) else 100 * (exp(x) - 1)
  tibble(index = recode(index, Pasture_mgmt_div = "Pasture mgmt", Water_mgmt_div = "Water mgmt"),
         hill = recode(hill, richness = "Richness (q=0)", shannon = "Shannon (q=1)"),
         likelihood = recode(lik, gaussian_log = "Gaussian on log", student_log = "Student-t on log",
                             gaussian_raw = "Gaussian on raw"),
         pct_effect = sprintf("%+.1f%% [%+.1f, %+.1f]", to_pct(b["Estimate"]), to_pct(b["Q2.5"]), to_pct(b["Q97.5"])),
         pct_est = round(to_pct(b["Estimate"]), 1),
         max_rhat = round(max(rhat(f), na.rm = TRUE), 3))
})
write_csv(lik_summ, "Derived/Excels/Likelihood_sensitivity.csv")

cat("\n== (c) response scale / likelihood family ==\n")
cat("Focal effect as % change in diversity per +1 SD of the index (climate spec, < 300 m):\n")
lik_summ %>% select(index, hill, likelihood, pct_effect) %>%
  pivot_wider(names_from = likelihood, values_from = pct_effect) %>%
  arrange(index, hill) %>% print(width = Inf)
