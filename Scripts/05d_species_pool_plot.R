# Figures for the species-pool / precipitation-confounding analysis ----

### Reads the fits from Scripts/04d_species_pool_test.R and draws:
###   p_pool_forest   -- the management-index coefficient and the spp_pool coefficient across the four environmental-adjustment blocks
###   p_pool_varpart  -- how the range-map species pool decomposes onto the elevation and precipitation gradients
### Safe to source() -- reads cached models + a farm-level table, never refits.

# Setup ----
library(tidyverse)
library(brms)
library(cowplot)
ggplot2::theme_set(theme_cowplot(11))
dir.create("Figures", showWarnings = FALSE)

index_lab <- c(Land_use_div = "Land use", Water_mgmt_div = "Water mgmt",
               Pasture_mgmt_div = "Pasture mgmt", All_practices_div = "All practices")
hill_lab  <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)")
env_lab   <- c(clim = "clim\nElev² + precip²",
               pool = "pool\nspp_pool",
               pool_elev = "pool_elev\nspp_pool + Elev²",
               pool_precip = "pool_precip\nspp_pool + precip²")
env_order <- names(env_lab)

# Coefficients from the 04d fits ----

fit_files <- list.files("Derived/models", pattern = "^pooltest__.*\\.rds$", full.names = TRUE)

coefs <- map(fit_files, function(f) {
  fit <- readRDS(f)
  parts <- str_match(basename(f), "^pooltest__(richness|shannon)__(.+)__(clim|pool|pool_elev|pool_precip)\\.rds$")
  fx <- fixef(fit)
  terms <- intersect(c("focal_z", "pool_point_z"), rownames(fx))
  map(terms, ~ tibble(
    hill = parts[2], index = parts[3], env = parts[4],
    term = if (.x == "focal_z") "Management index" else "spp_pool",
    est = fx[.x, "Estimate"], lo = fx[.x, "Q2.5"], hi = fx[.x, "Q97.5"]
  )) %>% list_rbind()
}) %>% list_rbind() %>%
  mutate(
    Index = factor(recode(index, !!!index_lab), levels = unname(index_lab)),
    Hill  = factor(recode(hill, !!!hill_lab), levels = unname(hill_lab)),
    Env   = factor(recode(env, !!!env_lab), levels = env_lab[env_order])
  )

# Figure: coefficient forest plot ----

p_pool_forest <- coefs %>%
  filter(index %in% c("Pasture_mgmt_div", "Water_mgmt_div")) %>%
  ggplot(aes(est, fct_rev(Env), colour = term, shape = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), fatten = 2) +
  facet_grid(Index ~ Hill) +
  scale_colour_manual(values = c("Management index" = "#1b7837", "spp_pool" = "#762a83"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(x = "Standardized effect on log diversity (posterior median, 95% CrI)", y = NULL,
       title = "Management effect and species-pool effect across environmental-adjustment blocks",
       subtitle = "Each block controls a different subset of {spp_pool, Elev², precip²} -- never all three (they are collinear)") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
ggsave("Figures/Species_pool_forest.png", p_pool_forest, width = 10, height = 6.5, bg = "white")
print(p_pool_forest)

# Figure: variance partitioning of spp_pool onto the two gradients ----

pe <- read_csv("Data/Farm_species_pool.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs)) %>%
  left_join(read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>% mutate(Id_gcs = as.character(Id_gcs)),
            by = "Id_gcs") %>%
  left_join(read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
              mutate(Id_gcs = as.character(Id_gcs)) %>% distinct(Id_gcs, canopy_10k),
            by = "Id_gcs")

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
  labs(x = "Variance of the range-map species pool (spp_pool) across the 73 farms", y = NULL,
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
