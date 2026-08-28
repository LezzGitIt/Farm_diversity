# Figures for the bird-diversity ~ farm-management-diversification models ----

### Reads the fitted models and the summary table written by Scripts/04_farm_mgmt_mod.R and builds the figures. Every plot object is saved to Figures/ and printed, and the script is safe to source() -- it never refits. Plot objects are named for reuse in the .qmd report and the manuscript.

# Setup ----
library(tidyverse)
library(brms)
library(cowplot)

ggplot2::theme_set(theme_cowplot())
dir.create("Figures", showWarnings = FALSE)

# Labels ----

index_labels <- c(
  baseline = "(baseline: no index)",
  Land_use_div = "Land use", Water_mgmt_div = "Water management",
  Pasture_mgmt_div = "Pasture management", All_practices_div = "All practices"
)
index_order <- unname(index_labels)

spec_labels <- c(
  ecoregion = "Ecoregion (fixed)", climate = "Precip + elevation",
  ecoregion_numhab = "Ecoregion + Num.hab (sensitivity)"
)
spec_colours <- c("Ecoregion (fixed)" = "#1b7837", "Precip + elevation" = "#762a83",
                  "Ecoregion + Num.hab (sensitivity)" = "#d95f02")

hill_labels <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)", simpson = "Simpson (q = 2)")
hill_order <- unname(hill_labels)

identifiable_indices <- c("Land_use_div", "All_practices_div")

# Load 04 outputs ----

Model_summaries <- read_csv("Derived/Excels/Farm_mgmt_model_summaries.csv", show_col_types = FALSE) %>%
  mutate(
    Hill = factor(recode(hill, !!!hill_labels), levels = hill_order),
    Index = factor(recode(index, !!!index_labels), levels = index_order),
    Spec = recode(spec, !!!spec_labels)
  )

Model_data <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE)

## brms fits, keyed by the 04 filename convention mod_<hill>__<index>__<spec>.rds
mod_fit_files <- list.files("Derived/models", pattern = "^mod_.*\\.rds$", full.names = TRUE)
mod_fits <- set_names(map(mod_fit_files, readRDS), str_remove(basename(mod_fit_files), "\\.rds$"))

# Figure 1: management-diversification coefficient across specs ----

p_index_effects <- Model_summaries %>%
  filter(term == "focal_z", spec %in% c("ecoregion", "climate")) %>%
  ggplot(aes(estimate, fct_rev(Index), colour = Spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high),
                  position = position_dodge(width = 0.5), fatten = 2.5) +
  facet_wrap(~Hill) +
  scale_colour_manual(values = spec_colours, name = "Region adjustment") +
  labs(
    x = "Standardized effect on log diversity (posterior median, 90% CrI)", y = NULL,
    title = "Farm management diversification vs bird diversity",
    subtitle = "Effect of a 1-SD increase in each index, adjusting for region + sampling + season"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Farm_mgmt_index_effects.png", p_index_effects, bg = "white", width = 10, height = 4)
print(p_index_effects)

# Figure 2: variance explained -- baseline vs + each index ----

## bayes_R2 of the no-index baseline (region + sampling + season) and of each index model, so the variance an index adds is visible.
p_bayes_r2 <- Model_summaries %>%
  distinct(Hill, Index, Spec, spec, bayes_R2) %>%
  filter(spec %in% c("ecoregion", "climate")) %>%
  ggplot(aes(bayes_R2, fct_rev(Index), colour = Spec)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  facet_wrap(~Hill) +
  scale_colour_manual(values = spec_colours, name = "Region adjustment") +
  labs(
    x = "Bayesian R²", y = NULL,
    title = "Variance explained: baseline vs adding each diversification index",
    subtitle = "Baseline = region + sampling effort + habitat count + day of year, no index"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Farm_mgmt_bayes_r2.png", p_bayes_r2, bg = "white", width = 10, height = 4)
print(p_bayes_r2)

# Figure 3: posterior predictive checks ----

ppcheck_keys <- names(mod_fits) %>% keep(~ str_detect(.x, "__All_practices_div__(ecoregion|climate)$"))
ppcheck_panels <- map(ppcheck_keys, function(key) {
  parts <- str_match(key, "^mod_(richness|shannon|simpson)__.+__(ecoregion|climate)$")
  suppressWarnings(pp_check(mod_fits[[key]], ndraws = 100)) +
    labs(subtitle = paste0(hill_labels[[parts[2]]], " -- ", spec_labels[[parts[3]]])) +
    theme(legend.position = "none")
})
p_ppcheck <- plot_grid(
  ggdraw() + draw_label("Posterior predictive checks (All practices models)", fontface = "bold", x = 0.02, hjust = 0),
  plot_grid(plotlist = ppcheck_panels, ncol = 2),
  ncol = 1, rel_heights = c(0.05, 1)
)
ggsave("Figures/Farm_mgmt_ppcheck.png", p_ppcheck, bg = "white", width = 10, height = 9)
print(p_ppcheck)

# Figure 4: conditional effect of the identifiable indices, on the diversity scale ----

index_scaling <- Model_data %>%
  summarize(across(all_of(setdiff(names(index_labels), "baseline")),
                   list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("index", "stat"),
               names_pattern = "^(.*_div)_(mean|sd)$") %>%
  pivot_wider(names_from = stat, values_from = value)

conditional_index_effect <- function(key) {
  parts <- str_match(key, "^mod_(richness|shannon|simpson)__(.+)__(ecoregion|climate)$")
  index <- parts[3]
  sc <- index_scaling %>% filter(index == !!index)
  as_tibble(conditional_effects(mod_fits[[key]], effects = "focal_z")[["focal_z"]]) %>%
    transmute(
      Hill = factor(hill_labels[[parts[2]]], levels = hill_order),
      Index = factor(index_labels[[index]], levels = index_order),
      Spec = spec_labels[[parts[4]]],
      index_value = focal_z * sc$sd + sc$mean,
      diversity = exp(estimate__), lower = exp(lower__), upper = exp(upper__)
    )
}

identifiable_keys <- names(mod_fits) %>%
  keep(~ str_detect(.x, paste0("__(", paste(identifiable_indices, collapse = "|"), ")__(ecoregion|climate)$")))
Conditional_lines <- map(identifiable_keys, conditional_index_effect) %>% list_rbind()

Raw_points <- Model_data %>%
  select(Hill, all_of(identifiable_indices), response) %>%
  mutate(Hill = factor(recode(Hill, !!!hill_labels), levels = hill_order)) %>%
  pivot_longer(all_of(identifiable_indices), names_to = "index", values_to = "index_value") %>%
  filter(!is.na(index_value)) %>%
  mutate(Index = factor(recode(index, !!!index_labels), levels = index_order))

p_index_conditional <- ggplot(Conditional_lines, aes(index_value, diversity)) +
  geom_point(data = Raw_points, aes(index_value, response), inherit.aes = FALSE, alpha = 0.25, size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Spec), alpha = 0.2) +
  geom_line(aes(colour = Spec), linewidth = 1) +
  facet_grid(Hill ~ Index, scales = "free_y") +
  scale_colour_manual(values = spec_colours, name = "Region adjustment", aesthetics = c("colour", "fill")) +
  labs(
    x = "Diversification index [0-1]", y = "Hill-number diversity",
    title = "Predicted bird diversity across the identifiable diversification indices",
    subtitle = "Other predictors at their means; points are the raw assemblage estimates"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Farm_mgmt_index_conditional.png", p_index_conditional, bg = "white", width = 10, height = 8)
print(p_index_conditional)
