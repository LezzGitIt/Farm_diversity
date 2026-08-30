# Figures for the species-pool / precipitation-confounding analysis ----

### Reads the fits from Scripts/04d_species_pool_test.R and draws:
###   p_env_specs     -- the management-index (pasture / water) coefficient across every environmental-adjustment spec, labelled by the variables it contains, ordered by how completely precipitation is controlled
###   p_pool_varpart  -- how the range-map species pool decomposes onto the elevation and precipitation gradients
### Safe to source() -- reads cached models + a farm-level table, never refits.

### p_env_specs merges what were two report objects (the pooltest forest + the precipitation functional-form table): both show the same thing -- the pasture / water coefficient tracks precipitation adjustment. The spp_pool coefficient itself is not plotted; the variance partition below already shows spp_pool is ~89% elevation + precipitation, and 04d shows adding it does not move the focal coefficient.

# Setup ----
library(tidyverse)
library(brms)
library(cowplot)
ggplot2::theme_set(theme_cowplot(11))
dir.create("Figures", showWarnings = FALSE)

index_lab <- c(Land_use_div = "Land use", Water_mgmt_div = "Water mgmt",
               Pasture_mgmt_div = "Pasture mgmt", All_practices_div = "All practices")
hill_lab  <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)")

### every spec also carries 10 km canopy + the sampling controls + REs; the label is just the elevation / precipitation / species-pool part
### `spec` is the key that appears in the model filename (pooltest__* env, or precipflex__* pspec); poly2 is the same formula as clim, so it is dropped in favour of clim
spec_lab <- c(
  pool_elev = "species pool + elev²",
  clim      = "elev² + precip²",
  splP      = "elev² + s(precip)",
  pool      = "species pool",
  pool_precip = "species pool + precip²",
  precOnly  = "precip² only",
  splPonly  = "s(precip) only"
)
spec_order <- names(spec_lab)   # top -> bottom: precipitation least controlled -> most controlled

# Focal coefficient from every 04d environmental spec ----

pull_focal <- function(f, spec) {
  b <- fixef(readRDS(f))["focal_z", ]
  tibble(spec = spec, est = b["Estimate"], lo = b["Q2.5"], hi = b["Q97.5"])
}

pooltest <- list.files("Derived/models", pattern = "^pooltest__(richness|shannon)__.+__(clim|pool|pool_elev|pool_precip)\\.rds$", full.names = TRUE) %>%
  map(function(f) {
    p <- str_match(basename(f), "^pooltest__(richness|shannon)__(.+)__(clim|pool|pool_elev|pool_precip)\\.rds$")
    bind_cols(tibble(hill = p[2], index = p[3]), pull_focal(f, p[4]))
  }) %>% list_rbind()

precipflex <- list.files("Derived/models", pattern = "^precipflex__(richness|shannon)__.+__(splP|precOnly|splPonly)\\.rds$", full.names = TRUE) %>%
  map(function(f) {
    p <- str_match(basename(f), "^precipflex__(richness|shannon)__(.+)__(splP|precOnly|splPonly)\\.rds$")
    bind_cols(tibble(hill = p[2], index = p[3]), pull_focal(f, p[4]))
  }) %>% list_rbind()

coefs <- bind_rows(pooltest, precipflex) %>%
  filter(index %in% c("Pasture_mgmt_div", "Water_mgmt_div"), spec %in% spec_order) %>%
  mutate(
    Index = factor(recode(index, !!!index_lab), levels = unname(index_lab)),
    Hill  = factor(recode(hill, !!!hill_lab), levels = unname(hill_lab)),
    Spec  = factor(recode(spec, !!!spec_lab), levels = spec_lab[spec_order])
  )

# Figure: focal coefficient across environmental specs ----

p_env_specs <- coefs %>%
  ggplot(aes(est, fct_rev(Spec))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#1b7837", fatten = 2) +
  facet_grid(Index ~ Hill) +
  labs(x = "Standardized effect of the management index on log diversity (posterior median, 95% CrI)", y = NULL,
       title = "The pasture / water coefficient tracks how completely precipitation is adjusted",
       subtitle = "Environmental term varied; 10 km canopy + sampling controls held fixed in every spec.\nTop rows: elevation competes with precipitation, or precipitation enters only weakly. Bottom rows: precipitation is the sole / dominant term.") +
  theme(axis.text.y = element_text(size = 8.5), plot.subtitle = element_text(size = 9.5))
ggsave("Figures/Species_pool_forest.png", p_env_specs, width = 10, height = 6, bg = "white")
print(p_env_specs)

# Figure: variance partitioning of spp_pool onto the two gradients ----

### one row per farm: pool + climate + a single canopy value (Farm_mgmt_model_data carries canopy per assemblage-year, so collapse to the farm mean before joining or the R2s are computed on duplicated rows)
canopy_by_farm <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  summarize(canopy_10k = mean(canopy_10k, na.rm = TRUE), .by = Id_gcs)

pe <- read_csv("Data/Farm_species_pool.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  left_join(read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>% mutate(Id_gcs = as.character(Id_gcs)),
            by = "Id_gcs") %>%
  left_join(canopy_by_farm, by = "Id_gcs")

r2 <- function(f) summary(lm(f, pe))$r.squared
full  <- r2(pool_point ~ poly(Elev_mean, 2) + poly(Tot_prec_mean, 2))
onlyE <- r2(pool_point ~ poly(Elev_mean, 2))
onlyP <- r2(pool_point ~ poly(Tot_prec_mean, 2))

varpart <- tibble(
  component = c("Unique to precipitation", "Unique to elevation", "Shared elev + precip", "Unexplained"),
  value = c(full - onlyE, full - onlyP, onlyE + onlyP - full, 1 - full)
) %>%
  mutate(component = factor(component, levels = rev(component)))

p_pool_varpart <- ggplot(varpart, aes(value, "spp_pool", fill = component)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * value)), position = position_stack(vjust = 0.5), size = 3.2, colour = "white") +
  scale_fill_manual(values = c("Unique to precipitation" = "#2166ac", "Unique to elevation" = "#b2182b",
                               "Shared elev + precip" = "#7f7f7f", "Unexplained" = "grey85"), name = NULL) +
  scale_x_continuous(labels = scales::percent, expand = expansion(0)) +
  labs(x = sprintf("Variance of the range-map species pool (spp_pool) across the %d farms", nrow(pe)), y = NULL,
       title = "spp_pool is a composite of the two environmental gradients",
       subtitle = sprintf("spp_pool ~ poly(Elev,2): R² = %.2f   |   ~ poly(precip,2): R² = %.2f   |   both: R² = %.2f", onlyE, onlyP, full)) +
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill = guide_legend(nrow = 2))
ggsave("Figures/Species_pool_varpart.png", p_pool_varpart, width = 9, height = 2.8, bg = "white")
print(p_pool_varpart)

# spp_pool correlations, for the text ----

pool_cor <- tibble(
  gradient = c("Elevation", "Precipitation", "Canopy cover (10 km)"),
  r = c(cor(pe$pool_point, pe$Elev_mean), cor(pe$pool_point, pe$Tot_prec_mean),
        cor(pe$pool_point, pe$canopy_10k, use = "pairwise"))
) %>% mutate(r = round(r, 2))
print(pool_cor)
