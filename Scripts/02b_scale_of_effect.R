# Scale of effect: which canopy-cover radius best explains bird diversity ----

### For each concentric-disc radius from Scripts/02a_extract_canopy_buffers.R, fit a mixed model of assemblage-level diversity on standardized canopy cover at that radius, and read off which radius adds the most explained variance (largest gain in marginal R^2 over a matching null without the canopy term; equivalently the largest AIC improvement). That radius is the "scale of effect" for canopy cover, and it is the radius Scripts/04a uses for the landscape-forest covariate.

### ONE adjustment set -- the same covariates the Scripts/04a 'component' spec uses (minus the landscape-forest term, which is what is being varied): elevation + elevation^2 + precipitation + precipitation^2 + the range-rarity-weighted species pool (pool_wes, Scripts/01b) + sampling effort (log number of point counts) + cyclic day-of-year + (1 | CollectorXyear) + (1 | Id_gcs). So the chosen radius is optimal for the model that is actually fitted downstream, not for some other adjustment set. (An earlier version of this script compared five exploratory adjustment sets -- raw / topography / Mundlak -- showing the canopy effect flips sign and scale with the adjustment; that was a diagnostic, not needed for the radius decision, and was removed 2026-09-02.)

### Primary response: q = 0 non-asymptotic richness (`No_Asy_TD`, coverage65 export). q = 1 / q = 2 asymptotic diversity are run alongside. The 3 Ref_* reference sites (not SCR farms) are dropped.

# Setup ----
library(tidyverse)
library(lme4)
library(MuMIn)
library(cowplot)

source("Scripts/Model_fns.R")   # latest_file()
ggplot2::theme_set(theme_cowplot())

dir.create("Figures", showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels/"

radii_m <- c(seq(200, 2000, by = 200), seq(3000, 10000, by = 1000))

# Load data ----

## Canopy cover per [assemblage x radius], frozen by Scripts/02a_extract_canopy_buffers.R (mean canopy within r m of the convex hull of the assemblage's point counts, from its own survey-year raster)
Canopy_by_scale <- read_csv("Data/Geospatial/Canopy_by_scale_assemblage.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

## Farm-level covariates: ecoregion, elevation, precipitation (Scripts/01a_farm_data.R) + range-rarity-weighted species pool (Scripts/01b_species_pool.R)
Farm_covariates <- read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  select(Id_gcs, Ecoregion, Elev_mean, Tot_prec_mean) %>%
  left_join(read_csv("Data/Farm_species_pool.csv", show_col_types = FALSE) %>%
              mutate(Id_gcs = as.character(Id_gcs)) %>% select(Id_gcs, pool_wes),
            by = "Id_gcs")

Site_covs <- read_csv(paste0(wrangling_excels, "Site_covs.csv"), show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  select(Id_muestreo_no_dc, Id_gcs)

## Assemblage = [data collector . farm . year-group . season]; one Event_covs row per point-count survey
Assemblage_surveys <- read_csv(paste0(wrangling_excels, "Event_covs.csv"), show_col_types = FALSE) %>%
  left_join(Site_covs, by = "Id_muestreo_no_dc") %>%
  filter(!is.na(Ano)) %>%
  mutate(
    Id_gcs = as.character(Id_gcs),
    Assemblage = str_replace_all(paste(Uniq_db, Id_gcs, Ano_grp, Season, sep = "."), " |-", "_"),
    CollectorXyear = paste(Uniq_db, Ano_grp, sep = "_")
  )

## Diversity estimates: q = 0 from the coverage65 export, q = 1 / q = 2 (asymptotic) from all_farms
all_farms_export <- latest_file("Derived/Excels", "^Tax_div_all_farms_.*\\.csv$")

Td_q0 <- read_csv(latest_file("Derived/Excels", "^Tax_div_coverage65_.*\\.csv$"), show_col_types = FALSE) %>%
  filter(Order.q == 0) %>%
  transmute(Assemblage, richness_q0 = No_Asy_TD)

Td_q12 <- read_csv(all_farms_export, show_col_types = FALSE) %>%
  filter(Order.q %in% c(1, 2)) %>%
  mutate(metric = if_else(Order.q == 1, "shannon_q1", "simpson_q2")) %>%
  select(Assemblage, metric, TD_asy) %>%
  pivot_wider(names_from = metric, values_from = TD_asy)

# Assemble the modelling frame ----

## Canopy is already per assemblage (02a); attach the collector-x-year batch key, sampling effort, mean day of year, farm covariates, and the diversity responses
Assemblage_keys <- Assemblage_surveys %>%
  summarize(CollectorXyear = first(CollectorXyear),
            Num_pc = n_distinct(Id_muestreo),
            doy = mean(Julian_day, na.rm = TRUE), .by = c(Assemblage, Id_gcs))

Scale_data <- Canopy_by_scale %>%
  filter(!str_detect(Id_gcs, "^Ref")) %>%
  left_join(Assemblage_keys, by = c("Assemblage", "Id_gcs")) %>%
  left_join(Farm_covariates, by = "Id_gcs") %>%
  left_join(Td_q0, by = "Assemblage") %>%
  left_join(Td_q12, by = "Assemblage") %>%
  mutate(
    num_pc_log = log(Num_pc),
    doy_sin = sin(2 * pi * doy / 365),
    doy_cos = cos(2 * pi * doy / 365)
  )

stopifnot(sum(is.na(Scale_data$Ecoregion)) == 0)

cat("Assemblages with canopy at all radii:", n_distinct(Scale_data$Assemblage),
    "| with q0 richness:", n_distinct(Scale_data$Assemblage[!is.na(Scale_data$richness_q0)]),
    "| with q1/q2:", n_distinct(Scale_data$Assemblage[!is.na(Scale_data$shannon_q1)]), "\n")

# Model specification ----

## The Scripts/04a 'component' adjustment set (minus the landscape-forest term, which is the focal canopy term here)
other_fixed <- paste("elev_z + I(elev_z^2) + precip_z + I(precip_z^2) + pool_wes_z",
                     "num_pc_log_z + doy_sin + doy_cos", sep = " + ")
re    <- "(1 | CollectorXyear) + (1 | Id_gcs)"
focal <- "canopy_z"

# Fit one radius x response ----

## MuMIn::r.squaredGLMM returns a named-row matrix ("delta" etc.) for GLMMs and a single unnamed row for Gaussian lmer -- take the right one either way
grab_r2 <- function(model) {
  x <- suppressWarnings(MuMIn::r.squaredGLMM(model))
  if ("delta" %in% rownames(x)) x["delta", ] else x[1, ]
}

fit_lmer <- function(formula_text, data) {
  tryCatch(
    suppressMessages(lmer(as.formula(formula_text), data, REML = FALSE)),
    error = function(e) NULL
  )
}

fit_scale_model <- function(radius, response_col) {
  d <- Scale_data %>%
    filter(radius_m == radius, !is.na(.data[[response_col]]), !is.na(canopy_cover),
           !is.na(doy_sin), !is.na(pool_wes)) %>%
    mutate(
      y = log(.data[[response_col]]),
      canopy_z = as.numeric(scale(canopy_cover)),
      num_pc_log_z = as.numeric(scale(num_pc_log)),
      elev_z = as.numeric(scale(Elev_mean)),
      precip_z = as.numeric(scale(Tot_prec_mean)),
      pool_wes_z = as.numeric(scale(pool_wes))
    )

  full <- fit_lmer(paste("y ~", focal, "+", other_fixed, "+", re), d)
  null <- fit_lmer(paste("y ~", other_fixed, "+", re), d)
  if (is.null(full) || is.null(null)) {
    return(tibble(response = response_col, radius_m = radius, n = nrow(d),
                  canopy_beta = NA_real_, canopy_lo = NA_real_, canopy_hi = NA_real_,
                  R2m_full = NA_real_, R2c_full = NA_real_, dR2m_canopy = NA_real_,
                  AIC_full = NA_real_, dAIC_vs_null = NA_real_, singular = NA))
  }

  r2_full <- grab_r2(full)
  r2_null <- grab_r2(null)
  ci <- tryCatch(confint(full, parm = focal, method = "Wald"),
                 error = function(e) c(NA_real_, NA_real_))
  tibble(
    response = response_col, radius_m = radius, n = nobs(full),
    canopy_beta = fixef(full)[[focal]], canopy_lo = ci[1], canopy_hi = ci[2],
    R2m_full = r2_full[["R2m"]], R2c_full = r2_full[["R2c"]],
    dR2m_canopy = r2_full[["R2m"]] - r2_null[["R2m"]],
    AIC_full = AIC(full), dAIC_vs_null = AIC(full) - AIC(null),
    singular = isSingular(full)
  )
}

# Run the grid ----

model_grid <- expand_grid(
  response = c("richness_q0", "shannon_q1", "simpson_q2"),
  radius = radii_m
)

Scale_results <- pmap(model_grid, function(response, radius) fit_scale_model(radius, response)) %>%
  list_rbind() %>%
  mutate(across(c(canopy_beta, canopy_lo, canopy_hi, R2m_full, R2c_full,
                  dR2m_canopy, AIC_full, dAIC_vs_null), ~ round(.x, 4)))

write_csv(Scale_results, "Derived/Excels/Scale_effect_results.csv")

## The scale of effect: radius maximising the canopy term's marginal-R^2 gain, per response
Best_scales <- Scale_results %>%
  filter(!is.na(dR2m_canopy)) %>%
  slice_max(dR2m_canopy, n = 1, by = response)

cat("\nBest-supported canopy radius (max gain in marginal R^2 over the matching null):\n")
print(Best_scales %>% select(response, radius_m, dR2m_canopy, canopy_beta,
                             canopy_lo, canopy_hi, dAIC_vs_null, singular), n = Inf)

# Plot: fit and effect size across radii ----

response_labels <- c(richness_q0 = "Richness (q = 0)",
                     shannon_q1 = "Shannon (q = 1)", simpson_q2 = "Simpson (q = 2)")

plot_data <- Scale_results %>%
  mutate(Response = factor(recode(response, !!!response_labels), levels = unname(response_labels)))

best_primary <- Best_scales %>% filter(response == "richness_q0")

p_scale_r2 <- ggplot(plot_data, aes(radius_m, dR2m_canopy)) +
  geom_vline(data = best_primary, aes(xintercept = radius_m), linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(linewidth = 0.9, colour = "#1b7837") +
  geom_point(size = 1.8, colour = "#1b7837") +
  facet_wrap(~Response, scales = "free_y") +
  scale_x_continuous(breaks = c(200, 600, 1200, 2000, 4000, 6000, 8000, 10000), trans = "sqrt") +
  labs(
    x = "Disc radius (m, sqrt scale)", y = expression(Delta ~ "marginal " * R^2 ~ "from canopy cover"),
    title = "Scale of effect: variance in diversity explained by canopy cover vs radius",
    subtitle = "Gain in marginal R-squared from adding canopy cover to the 'component' adjustment set. Dashed line = best radius for q = 0 richness."
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Scale_effect_canopy_r2.png", p_scale_r2, bg = "white", width = 12, height = 5.5)
print(p_scale_r2)

p_scale_beta <- ggplot(plot_data, aes(radius_m, canopy_beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = canopy_lo, ymax = canopy_hi), alpha = 0.15, fill = "#1b7837") +
  geom_line(linewidth = 0.9, colour = "#1b7837") +
  geom_point(size = 1.8, colour = "#1b7837") +
  facet_wrap(~Response) +
  scale_x_continuous(breaks = c(200, 600, 1200, 2000, 4000, 6000, 8000, 10000), trans = "sqrt") +
  labs(
    x = "Disc radius (m, sqrt scale)", y = "Standardized canopy-cover effect on log diversity",
    title = "Canopy-cover effect size across spatial scales",
    subtitle = "lmer fixed effect with 95% Wald interval, from the 'component' adjustment set."
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Scale_effect_canopy_beta.png", p_scale_beta, bg = "white", width = 12, height = 5.5)
print(p_scale_beta)
