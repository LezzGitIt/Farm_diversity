# Figures for the bird-diversity ~ management-diversification linking models ----

### Reads the fitted models and summary tables written by Scripts/03_linking_model.R and builds the figures. Each plot object is both saved to Figures/ and printed, and this script is safe to source() -- it never refits. The plot objects are named so they can be pulled into .qmd reports and the manuscript later.

# Setup ----
library(tidyverse)
library(brms)
library(cowplot)

ggplot2::theme_set(theme_cowplot())

dir.create("Figures", showWarnings = FALSE)

# Labels and shared scales ----

div_index_labels <- c(
  Land_use_div = "Land use",
  Water_mgmt_div = "Water management",
  Pasture_mgmt_div = "Pasture management",
  All_practices_div = "All practices"
)
div_level_order <- unname(div_index_labels)

adjustment_labels <- c(ecoregion = "Ecoregion (categorical)", climate = "Precip + elevation")
adjustment_colours <- c("Ecoregion (categorical)" = "#1b7837", "Precip + elevation" = "#762a83")

identifiable_indices <- c("Land_use_div", "All_practices_div")

# Load 03 outputs ----

Link_coefficients <- read_csv("Derived/Excels/Linking_model_coefficients.csv", show_col_types = FALSE) %>%
  mutate(
    Index = factor(recode(Index, !!!div_index_labels), levels = div_level_order),
    Adjustment = recode(Adjustment, !!!adjustment_labels)
  )

Model_data <- read_csv("Derived/Excels/Linking_model_data.csv", show_col_types = FALSE)

## brms fits, keyed by the 03 filename convention link_<hill>_<index>_<adjustment>.rds
link_fit_files <- list.files("Derived/models", pattern = "^link_.*\\.rds$", full.names = TRUE)
link_fits <- set_names(
  map(link_fit_files, readRDS),
  str_remove(basename(link_fit_files), "\\.rds$")
)

# Figure 1: management-diversification coefficient across specs ----

## The focal index effect (standardized, log-diversity scale) with its 90% credible interval, for every Hill number x index x region-adjustment combination. This is the headline: does farm management diversification track bird diversity once region is controlled for, and how much does that answer depend on how region is controlled for.
p_index_effects <- Link_coefficients %>%
  filter(term == "focal_z") %>%
  ggplot(aes(x = estimate, y = fct_rev(Index), colour = Adjustment)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(
    aes(xmin = conf_low, xmax = conf_high),
    position = position_dodge(width = 0.5), fatten = 2.5
  ) +
  facet_wrap(~Hill) +
  scale_colour_manual(values = adjustment_colours, name = "Region adjustment") +
  labs(
    x = "Standardized effect on log diversity (posterior median, 90% CrI)",
    y = NULL,
    title = "Farm management diversification vs bird diversity",
    subtitle = "Effect of a 1-SD increase in each index, after adjusting for region and sampling effort"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Linking_index_effects.png", p_index_effects, bg = "white", width = 10, height = 5)
print(p_index_effects)

# Figure 2: full fixed-effect structure ----

## Every fixed effect from every model, so the index effect can be read against the effort and region terms that dominate. Ecoregion / climate terms only appear for the adjustment set that uses them.
coef_term_labels <- c(
  focal_z = "Focal index",
  Num_pc_log_z = "log(point counts)",
  Tot_prec_mean_z = "Precipitation",
  Elev_mean_z = "Elevation",
  EcoregionBoyacasantander = "Ecoregion: Boyaca-Santander",
  EcoregionCafetera = "Ecoregion: Cafetera",
  EcoregionPiedemonte = "Ecoregion: Piedemonte",
  EcoregionRiocesar = "Ecoregion: Rio Cesar"
)

p_all_coefficients <- Link_coefficients %>%
  filter(term != "Intercept") %>%
  mutate(term = factor(recode(term, !!!coef_term_labels), levels = rev(unname(coef_term_labels)))) %>%
  ggplot(aes(x = estimate, y = term, colour = Adjustment)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(
    aes(xmin = conf_low, xmax = conf_high),
    position = position_dodge(width = 0.5), fatten = 1.8, size = 0.4
  ) +
  facet_grid(Hill ~ Index) +
  scale_colour_manual(values = adjustment_colours, name = "Region adjustment") +
  labs(
    x = "Standardized effect on log diversity (posterior median, 90% CrI)", y = NULL,
    title = "All fixed effects, by model"
  ) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9))
ggsave("Figures/Linking_all_coefficients.png", p_all_coefficients, bg = "white", width = 12, height = 7)
print(p_all_coefficients)

# Figure 3: posterior predictive checks ----

## Density overlay of the observed log-diversity against draws from the posterior predictive, for the "All practices" models (representative of all specs)
ppcheck_keys <- names(link_fits) %>% keep(~ str_detect(.x, "All_practices_div"))
ppcheck_panels <- map(ppcheck_keys, function(key) {
  parts <- str_match(key, "^link_(shannon|simpson)_(.+)_(ecoregion|climate)$")
  pp_check(link_fits[[key]], ndraws = 100) +
    labs(subtitle = paste0(str_to_title(parts[2]), " -- ", adjustment_labels[[parts[4]]])) +
    theme(legend.position = "none")
})
p_ppcheck <- plot_grid(plotlist = ppcheck_panels, ncol = 2)
p_ppcheck <- plot_grid(
  ggdraw() + draw_label("Posterior predictive checks (All practices models)", fontface = "bold", x = 0.02, hjust = 0),
  p_ppcheck, ncol = 1, rel_heights = c(0.06, 1)
)
ggsave("Figures/Linking_ppcheck.png", p_ppcheck, bg = "white", width = 10, height = 8)
print(p_ppcheck)

# Figure 4: conditional effect of the identifiable indices, on the diversity scale ----

## Predicted diversity across the observed range of each index (back-transformed from log), holding the other predictors at their means, for land use and all practices -- the two indices least confounded with ecoregion. Points are the raw assemblage estimates.
index_scaling <- Model_data %>%
  summarize(across(all_of(names(div_index_labels)),
                   list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(),
               names_to = c("Index", "stat"), names_pattern = "^(.*_div)_(mean|sd)$") %>%
  pivot_wider(names_from = stat, values_from = value)

conditional_index_effect <- function(key) {
  parts <- str_match(key, "^link_(shannon|simpson)_(.+)_(ecoregion|climate)$")
  index <- parts[3]
  sc <- index_scaling %>% filter(Index == index)
  as_tibble(conditional_effects(link_fits[[key]], effects = "focal_z")[["focal_z"]]) %>%
    transmute(
      Hill = str_to_title(parts[2]),
      Index = factor(div_index_labels[[index]], levels = div_level_order),
      Adjustment = adjustment_labels[[parts[4]]],
      index_value = focal_z * sc$sd + sc$mean,
      diversity = exp(estimate__),
      lower = exp(lower__), upper = exp(upper__)
    )
}

identifiable_keys <- names(link_fits) %>%
  keep(~ str_detect(.x, paste(identifiable_indices, collapse = "|")))
Conditional_lines <- map(identifiable_keys, conditional_index_effect) %>% list_rbind()

Raw_points <- Model_data %>%
  select(Hill, all_of(identifiable_indices), TD_asy) %>%
  pivot_longer(all_of(identifiable_indices), names_to = "Index", values_to = "index_value") %>%
  filter(!is.na(index_value)) %>%
  mutate(Index = factor(recode(Index, !!!div_index_labels), levels = div_level_order))

p_index_conditional <- ggplot(Conditional_lines, aes(index_value, diversity)) +
  geom_point(
    data = Raw_points, aes(x = index_value, y = TD_asy), inherit.aes = FALSE,
    alpha = 0.25, size = 1
  ) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Adjustment), alpha = 0.2) +
  geom_line(aes(colour = Adjustment), linewidth = 1) +
  facet_grid(Hill ~ Index, scales = "free") +
  scale_colour_manual(values = adjustment_colours, name = "Region adjustment", aesthetics = c("colour", "fill")) +
  labs(
    x = "Diversification index [0-1]", y = "Asymptotic Hill-number diversity",
    title = "Predicted bird diversity across the identifiable diversification indices",
    subtitle = "Other predictors held at their means; points are the raw assemblage estimates"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Linking_index_conditional.png", p_index_conditional, bg = "white", width = 10, height = 7)
print(p_index_conditional)
