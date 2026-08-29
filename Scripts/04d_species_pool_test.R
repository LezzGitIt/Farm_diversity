# Does the range-map species pool add anything to the farm-management model? ----

### The DAG (Scripts/dag.R) flags the regional species pool as the one confounder the "climate" model version might not fully capture with elevation + precipitation + canopy. Scripts/01b_species_pool.R builds a per-farm potential species pool (`pool_point`, from Ayerbe ranges + Suarez-Castro 2024 elevational limits). This script tests whether it earns a place in Scripts/04_farm_mgmt_mod.R.

### pool_point is collinear with the climate terms (pool ~ poly(Elev,2) + poly(precip,2) R^2 = 0.86), so putting all three in one model is uninformative (that fit throws ESS / collinearity warnings). Instead, following Aaron's suggestion, fit the climate model with FOUR environmental-adjustment blocks -- never all three axes at once -- and compare the focal-index coefficient and fit:
###   clim        : poly(Elev,2) + poly(precip,2) + canopy            [the current 04 "climate" spec]
###   pool        : pool_point + canopy                                [pool replaces both climate axes]
###   pool_elev   : pool_point + poly(Elev,2) + canopy                 [pool + elevation]
###   pool_precip : pool_point + poly(precip,2) + canopy               [pool + precipitation]
### All four also carry the standard sampling controls + REs. If the focal coefficient is stable across blocks and the pool-based blocks fit about as well, pool is a valid (more parsimonious) substitute; if the focal coefficient moves, the blocks disagree about the confounder and that is worth knowing.

### Reads Scripts/04_farm_mgmt_mod.R's persisted frame + Data/Farm_species_pool.csv. Primary analysis (< 300 m). Outputs Derived/Excels/Species_pool_model_test.csv.

# Setup ----
library(tidyverse)
library(brms)
source("Scripts/Farm_diversity_fns.R")
options(mc.cores = 4, brms.backend = "rstan")
dir.create("Derived/models", recursive = TRUE, showWarnings = FALSE)
set.seed(1989)

chains <- 4; iter <- 3000; warmup <- 1000; adapt_delta <- 0.99
mod_priors <- c(
  prior(student_t(3, 3, 2.5), class = "Intercept"),
  prior(normal(0, 0.75), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

div_indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")
dist_threshold <- 300

# Data ----

pool <- read_csv("Data/Farm_species_pool.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>% select(Id_gcs, pool_point)

Model_data <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  left_join(pool, by = "Id_gcs") %>%
  mutate(
    pool_point_z = as.numeric(scale(pool_point)),
    doy_sin = sin(2 * pi * doy / 365),   # 04's CSV carries doy, not the cyclic terms
    doy_cos = cos(2 * pi * doy / 365)
  )

# Model blocks ----

env_blocks <- c(
  clim        = "Elev_mean_z + I(Elev_mean_z^2) + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z",
  pool        = "pool_point_z + canopy_10k_z",
  pool_elev   = "pool_point_z + Elev_mean_z + I(Elev_mean_z^2) + canopy_10k_z",
  pool_precip = "pool_point_z + Tot_prec_mean_z + I(Tot_prec_mean_z^2) + canopy_10k_z"
)
sampling <- "Num_pc_log_z + Num_hab_z + doy_sin + doy_cos + (1 | Id_gcs) + (1 | CollectorXyear)"

formula_for <- function(env) {
  bf(as.formula(paste0("log_response | resp_se(se_log, sigma = TRUE) ~ focal_z + ",
                       env_blocks[[env]], " + ", sampling)))
}

frame_for <- function(hill, index) {
  Model_data %>%
    filter(Hill == hill, !is.na(se_log), !is.na(canopy_10k_z), !is.na(pool_point_z),
           !is.na(Num_pc_log_z), !is.na(Num_hab_z), !is.na(Tot_prec_mean_z), !is.na(Elev_mean_z),
           !is.na(dist_farm), dist_farm < dist_threshold) %>%
    mutate(focal_z = .data[[paste0(index, "_z")]]) %>%
    filter(!is.na(focal_z))
}

# Fit the grid ----

grid <- expand_grid(hill = c("richness", "shannon"), index = div_indices, env = names(env_blocks)) %>%
  mutate(key = paste("pooltest", hill, index, env, sep = "__"))

base_by_env <- list()
fits <- vector("list", nrow(grid)); names(fits) <- grid$key

for (i in seq_len(nrow(grid))) {
  row <- as.list(grid[i, ])
  message(sprintf("[%d/%d] %s", i, nrow(grid), row$key))
  frame <- frame_for(row$hill, row$index)
  common <- list(chains = chains, iter = iter, warmup = warmup, seed = 1989,
                 control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
                 file = sprintf("Derived/models/%s", row$key), file_refit = "on_change")
  base <- base_by_env[[row$env]]
  fit <- if (is.null(base)) {
    do.call(brm, c(list(formula = formula_for(row$env), data = frame, prior = mod_priors), common))
  } else {
    do.call(update, c(list(object = base, newdata = frame, recompile = FALSE), common))
  }
  if (is.null(base)) base_by_env[[row$env]] <- fit
  fits[[row$key]] <- fit
}

# Summarise ----

summ <- pmap_dfr(grid, function(hill, index, env, key) {
  f <- fits[[key]]
  b <- fixef(f)["focal_z", ]
  pool_row <- if ("pool_point_z" %in% rownames(fixef(f))) fixef(f)["pool_point_z", ] else rep(NA, 4)
  ess <- tryCatch({
    s <- summary(f)$fixed
    min(s[, "Bulk_ESS"], na.rm = TRUE)
  }, error = function(e) NA)
  tibble(
    hill, index, env,
    focal = sprintf("%+.3f [%+.3f, %+.3f]", b["Estimate"], b["Q2.5"], b["Q97.5"]),
    focal_est = round(b["Estimate"], 3),
    pool_coef = if (is.na(pool_row[1])) "" else sprintf("%+.3f [%+.3f, %+.3f]", pool_row[1], pool_row[3], pool_row[4]),
    bayes_R2 = round(bayes_R2(f)[, "Estimate"], 3),
    min_bulk_ess = round(ess),
    max_rhat = round(max(rhat(f), na.rm = TRUE), 3),
    n_div = sum(subset(nuts_params(f), Parameter == "divergent__")$Value)
  )
})

write_csv(summ, "Derived/Excels/Species_pool_model_test.csv")

# Report ----

cat("\n== Variance partitioning: pool_point ~ climate (farm level, lm) ==\n")
pe <- read_csv("Data/Farm_species_pool.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  left_join(read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>% mutate(Id_gcs = as.character(Id_gcs)),
            by = "Id_gcs")
r2 <- function(f) round(summary(lm(f, pe))$r.squared, 3)
full <- r2(pool_point ~ poly(Elev_mean, 2) + poly(Tot_prec_mean, 2))
cat("  ~ poly(Elev,2)            :", r2(pool_point ~ poly(Elev_mean, 2)), "\n")
cat("  ~ poly(precip,2)          :", r2(pool_point ~ poly(Tot_prec_mean, 2)), "\n")
cat("  ~ both                    :", full, "\n")
cat("  unique to elev            :", full - r2(pool_point ~ poly(Tot_prec_mean, 2)), "\n")
cat("  unique to precip          :", full - r2(pool_point ~ poly(Elev_mean, 2)), "\n")

cat("\n== Focal coefficient by environmental block (primary analysis) ==\n")
summ %>%
  select(hill, index, env, focal) %>%
  pivot_wider(names_from = env, values_from = focal) %>%
  arrange(hill, index) %>%
  print(width = Inf)

cat("\n== Fit (bayes_R2) and diagnostics by block ==\n")
summ %>%
  select(hill, index, env, bayes_R2, min_bulk_ess, max_rhat, n_div) %>%
  arrange(hill, index, env) %>%
  print(n = Inf)
