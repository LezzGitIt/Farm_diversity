# Scale of effect: which canopy-cover radius best explains bird diversity ----

### For each concentric-disc radius from Scripts/02a_extract_canopy_buffers.R, fit a mixed model of assemblage-level diversity on standardized canopy cover at that radius, and read off which radius adds the most explained variance (largest gain in marginal R^2 over a matching null without the canopy term; equivalently the largest AIC improvement). That radius is the "scale of effect" for canopy cover.

### Every model also controls for sampling effort (log number of point counts), the number of habitat types sampled, and a cyclic term on the assemblage's mean day of year (a Nearctic-migrant-season proxy). Five region-adjustment specifications:
###   sampling        canopy                                            + (1|CollectorXyear)
###   sampling_farm   "                                                  + (1|CollectorXyear) + (1|Id_gcs)
###   topography      canopy + elevation + precipitation                 + (1|CollectorXyear)
###   topography_farm "                                                  + (1|CollectorXyear) + (1|Id_gcs)
###   within_eco      canopy_within_ecoregion + Ecoregion                + (1|CollectorXyear)
### `within_eco` is a Mundlak-style within-region check: canopy is centred on its ecoregion mean, and Ecoregion (fixed) absorbs all between-region variation, so the canopy coefficient is the effect of a farm being more wooded than its regional peers.

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

## Farm-level covariates: ecoregion, elevation, precipitation (Scripts/01a_farm_data.R)
Farm_covariates <- read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  select(Id_gcs, Ecoregion, Elev_mean, Tot_prec_mean)

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

## Number of habitat types sampled per assemblage (constant across Order.q)
Assemblage_numhab <- read_csv(all_farms_export, show_col_types = FALSE) %>%
  distinct(Assemblage, Num.hab)

# Assemble the modelling frame ----

## Canopy is already per assemblage (06a); attach the collector-x-year batch key, sampling effort, mean day of year, farm covariates, habitat count, and the diversity responses
Assemblage_keys <- Assemblage_surveys %>%
  summarize(CollectorXyear = first(CollectorXyear),
            Num_pc = n_distinct(Id_muestreo),
            doy = mean(Julian_day, na.rm = TRUE), .by = c(Assemblage, Id_gcs))

Scale_data <- Canopy_by_scale %>%
  filter(!str_detect(Id_gcs, "^Ref")) %>%
  left_join(Assemblage_keys, by = c("Assemblage", "Id_gcs")) %>%
  left_join(Farm_covariates, by = "Id_gcs") %>%
  left_join(Assemblage_numhab, by = "Assemblage") %>%
  left_join(Td_q0, by = "Assemblage") %>%
  left_join(Td_q12, by = "Assemblage") %>%
  mutate(
    num_pc_log = log(Num_pc),
    num_hab_num = as.numeric(as.character(Num.hab)),
    doy_sin = sin(2 * pi * doy / 365),
    doy_cos = cos(2 * pi * doy / 365)
  )

stopifnot(sum(is.na(Scale_data$Ecoregion)) == 0)

cat("Assemblages with canopy at all radii:", n_distinct(Scale_data$Assemblage),
    "| with q0 richness:", n_distinct(Scale_data$Assemblage[!is.na(Scale_data$richness_q0)]),
    "| with q1/q2:", n_distinct(Scale_data$Assemblage[!is.na(Scale_data$shannon_q1)]), "\n")

# Model specifications ----

## Each spec: the fixed effects besides the focal canopy term, the random-effect terms, and which column is the focal canopy term
## Every spec carries the same sampling / seasonal controls; they differ only in how region is adjusted for
sampling_controls <- "num_pc_log_z + num_hab_z + doy_sin + doy_cos"
specs <- tribble(
  ~spec,             ~region_fixed,          ~re,                                    ~focal,
  "sampling",        NA,                     "(1 | CollectorXyear)",                  "canopy_z",
  "sampling_farm",   NA,                     "(1 | CollectorXyear) + (1 | Id_gcs)",   "canopy_z",
  "topography",      "elev_z + precip_z",    "(1 | CollectorXyear)",                  "canopy_z",
  "topography_farm", "elev_z + precip_z",    "(1 | CollectorXyear) + (1 | Id_gcs)",   "canopy_z",
  "within_eco",      "Ecoregion",            "(1 | CollectorXyear)",                  "canopy_within_z"
) %>%
  mutate(other_fixed = if_else(is.na(region_fixed), sampling_controls,
                               paste(sampling_controls, region_fixed, sep = " + ")))

# Fit one radius x response x spec ----

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

fit_scale_model <- function(radius, response_col, spec, other_fixed, re, focal) {
  d <- Scale_data %>%
    filter(radius_m == radius, !is.na(.data[[response_col]]), !is.na(canopy_cover),
           !is.na(doy_sin), !is.na(num_hab_num)) %>%
    mutate(
      y = log(.data[[response_col]]),
      canopy_z = as.numeric(scale(canopy_cover)),
      canopy_within_z = as.numeric(scale(canopy_cover - ave(canopy_cover, Ecoregion))),
      num_pc_log_z = as.numeric(scale(num_pc_log)),
      num_hab_z = as.numeric(scale(num_hab_num)),
      elev_z = as.numeric(scale(Elev_mean)),
      precip_z = as.numeric(scale(Tot_prec_mean)),
      Ecoregion = factor(Ecoregion)
    )

  full <- fit_lmer(paste("y ~", focal, "+", other_fixed, "+", re), d)
  null <- fit_lmer(paste("y ~", other_fixed, "+", re), d)
  if (is.null(full) || is.null(null)) {
    return(tibble(response = response_col, spec = spec, radius_m = radius, n = nrow(d),
                  canopy_beta = NA_real_, canopy_lo = NA_real_, canopy_hi = NA_real_,
                  R2m_full = NA_real_, R2c_full = NA_real_, dR2m_canopy = NA_real_,
                  AIC_full = NA_real_, dAIC_vs_null = NA_real_, singular = NA))
  }

  r2_full <- grab_r2(full)
  r2_null <- grab_r2(null)
  ci <- tryCatch(confint(full, parm = focal, method = "Wald"),
                 error = function(e) c(NA_real_, NA_real_))
  tibble(
    response = response_col, spec = spec, radius_m = radius, n = nobs(full),
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
  specs,
  radius = radii_m
)

Scale_results <- pmap(model_grid, function(response, spec, other_fixed, re, focal, radius, ...) {
  fit_scale_model(radius, response, spec, other_fixed, re, focal)
}) %>%
  list_rbind() %>%
  mutate(across(c(canopy_beta, canopy_lo, canopy_hi, R2m_full, R2c_full,
                  dR2m_canopy, AIC_full, dAIC_vs_null), ~ round(.x, 4)))

write_csv(Scale_results, "Derived/Excels/Scale_effect_results.csv")

## The scale of effect: radius maximising the canopy term's marginal-R^2 gain, per response x spec
Best_scales <- Scale_results %>%
  filter(!is.na(dR2m_canopy)) %>%
  slice_max(dR2m_canopy, n = 1, by = c(response, spec))

cat("\nBest-supported canopy radius (max gain in marginal R^2 over the matching null):\n")
print(Best_scales %>% select(response, spec, radius_m, dR2m_canopy, canopy_beta,
                             canopy_lo, canopy_hi, dAIC_vs_null, singular), n = Inf)

# Plot: fit and effect size across radii ----

response_labels <- c(richness_q0 = "Richness (q = 0)",
                     shannon_q1 = "Shannon (q = 1)", simpson_q2 = "Simpson (q = 2)")
spec_order <- c("sampling", "sampling_farm", "topography", "topography_farm", "within_eco")

plot_data <- Scale_results %>%
  mutate(
    Response = factor(recode(response, !!!response_labels), levels = unname(response_labels)),
    spec = factor(spec, levels = spec_order)
  )

best_primary <- Best_scales %>% filter(response == "richness_q0", spec == "sampling")

p_scale_r2 <- ggplot(plot_data, aes(radius_m, dR2m_canopy, colour = spec)) +
  geom_vline(data = best_primary, aes(xintercept = radius_m), linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~Response, scales = "free_y") +
  scale_x_continuous(breaks = radii_m, trans = "sqrt") +
  scale_colour_brewer(palette = "Dark2", name = "Model") +
  labs(
    x = "Disc radius (m, sqrt scale)", y = expression(Delta ~ "marginal " * R^2 ~ "from canopy cover"),
    title = "Scale of effect: variance in diversity explained by canopy cover vs radius",
    subtitle = "Gain in marginal R-squared over the matching null. Dashed line = best radius for q = 0 richness (sampling model)."
  ) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Scale_effect_canopy_r2.png", p_scale_r2, bg = "white", width = 12, height = 5.5)
print(p_scale_r2)

p_scale_beta <- ggplot(plot_data, aes(radius_m, canopy_beta, colour = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = canopy_lo, ymax = canopy_hi, fill = spec), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~Response) +
  scale_x_continuous(breaks = radii_m, trans = "sqrt") +
  scale_colour_brewer(palette = "Dark2", name = "Model", aesthetics = c("colour", "fill")) +
  labs(
    x = "Disc radius (m, sqrt scale)", y = "Standardized canopy-cover effect on log diversity",
    title = "Canopy-cover effect size across spatial scales",
    subtitle = "lmer fixed effect with 95% Wald interval. within_eco: coefficient is the within-ecoregion (Mundlak) canopy effect."
  ) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Scale_effect_canopy_beta.png", p_scale_beta, bg = "white", width = 12, height = 5.5)
print(p_scale_beta)
