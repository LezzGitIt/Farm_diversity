# Figures for the bird-diversity ~ farm-management-diversification analysis ----

### The manuscript / summary-report figures for the question of interest. Reads the summary tables + fitted models from Scripts/04a (primary: both responses x both adjustment sets) and Scripts/04b (the Piedemonte cut; the forest_only survey-scope interaction). Safe to source() -- never refits. Every plot object is saved, printed, and named for reuse in the .qmd reports.

### (Was three scripts: 05 = the primary forest / R2 / pp-check / conditional-effect figures; 05c = the Piedemonte forest plot; the incidence-response + distance-sensitivity figures used to be built in 04g. The species-pool diagnostic figures moved to Scripts/04b's pool_blocks section.)

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

## climate (= the 'component' spec in the write-up) = the DAG-sufficient primary (temp² + precip² + landscape forest + endemism index); ecoregion = the proxy robustness check, Ecoregion alone (Scripts/dag.R)
spec_labels <- c(climate = "Component", ecoregion = "Ecoregion")
spec_colours <- c(Component = "#762a83", Ecoregion = "#1b7837")

hill_labels <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)", simpson = "Simpson (q = 2)")
hill_order <- unname(hill_labels)

identifiable_indices <- c("Land_use_div", "All_practices_div")

# Load 04 outputs ----

Model_summaries <- read_csv("Derived/Excels/Farm_mgmt_model_summaries.csv", show_col_types = FALSE) %>%
  ## tolerate a summaries file written before the primary/full split
  { if ("data_subset" %in% names(.)) . else mutate(., data_subset = "primary") } %>%
  mutate(
    Hill = factor(recode(hill, !!!hill_labels), levels = hill_order),
    Index = factor(recode(index, !!!index_labels), levels = index_order),
    Spec = recode(spec, !!!spec_labels)
  )

## the primary figures use the < 300 m analysis set; the full set feeds the sensitivity figure only
Primary <- Model_summaries %>% filter(data_subset == "primary")

Model_data <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE)

## brms fits, keyed by the 04 filename convention mod_<hill>__<index>__<spec>__<data_subset>.rds
mod_fit_files <- list.files("Derived/models", pattern = "^mod_.*\\.rds$", full.names = TRUE)
mod_fits <- set_names(map(mod_fit_files, readRDS), str_remove(basename(mod_fit_files), "\\.rds$"))

## the fit-dependent figures (pp-checks, conditional effects) are skipped when no .rds are present, e.g. an interim render before 04 has run
have_fits <- length(mod_fits) > 0
if (!have_fits) message("05: no fitted models in Derived/models/ -- skipping pp-check and conditional-effect figures.")

# Figure 1: management-diversification coefficient, primary vs robustness ----

p_index_effects <- Primary %>%
  filter(term == "focal_z") %>%
  ggplot(aes(estimate, fct_rev(Index), colour = Spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high),
                  position = position_dodge(width = 0.5), fatten = 2.5) +
  facet_wrap(~Hill) +
  scale_colour_manual(values = spec_colours, name = NULL) +
  labs(
    x = "Standardized effect on log diversity (posterior median, 90% CrI)", y = NULL,
    title = "Farm management diversification vs bird diversity",
    subtitle = "1-SD increase in each index; assemblages < 300 m from the farm.\nComponent (temperature + precipitation + landscape forest + endemism index) = primary, Ecoregion = robustness"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Farm_mgmt_index_effects.png", p_index_effects, bg = "white", width = 10, height = 4)
print(p_index_effects)

# Figure 2: variance explained -- baseline vs + each index ----

p_bayes_r2 <- Primary %>%
  distinct(Hill, Index, Spec, spec, bayes_R2) %>%
  ggplot(aes(bayes_R2, fct_rev(Index), colour = Spec)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  facet_wrap(~Hill) +
  scale_colour_manual(values = spec_colours, name = NULL) +
  labs(
    x = "Bayesian R²", y = NULL,
    title = "Variance explained: baseline vs adding each diversification index",
    subtitle = "Baseline = region + sampling effort + day of year + random effects, no index"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Farm_mgmt_bayes_r2.png", p_bayes_r2, bg = "white", width = 10, height = 4)
print(p_bayes_r2)

## The distance-to-farm cutoff sensitivity (both responses) is Figure 6 below.

# Figure 3: posterior predictive checks ----

if (have_fits) {

ppcheck_keys <- names(mod_fits) %>% keep(~ str_detect(.x, "__All_practices_div__(ecoregion|climate)__primary$"))
ppcheck_panels <- map(ppcheck_keys, function(key) {
  parts <- str_match(key, "^mod_(richness|shannon|simpson)__.+__(ecoregion|climate)__primary$")
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
  parts <- str_match(key, "^mod_(richness|shannon|simpson)__(.+)__(ecoregion|climate)__primary$")
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
  keep(~ str_detect(.x, paste0("__(", paste(identifiable_indices, collapse = "|"), ")__(ecoregion|climate)__primary$")))
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
  scale_colour_manual(values = spec_colours, name = NULL, aesthetics = c("colour", "fill")) +
  labs(
    x = "Diversification index [0-1]", y = "Hill-number diversity",
    title = "Predicted bird diversity across the identifiable diversification indices",
    subtitle = "Other predictors at their means; points are the raw assemblage estimates"
  ) +
  theme(legend.position = "bottom")
ggsave("Figures/Farm_mgmt_index_conditional.png", p_index_conditional, bg = "white", width = 10, height = 8)
print(p_index_conditional)

} else {
  p_ppcheck <- NULL
  p_index_conditional <- NULL
}

# Figure 5: coverage/Cmax vs point-count-standardised (incidence) response ----

### was Scripts/04g's p_effects. The focal management coefficient under each response, both adjustment sets.
inc_index_lab <- c(Land_use_div = "Land use", Water_mgmt_div = "Water mgmt",
                   Pasture_mgmt_div = "Pasture mgmt", All_practices_div = "All practices")

focal_both_responses <- bind_rows(
  Model_summaries %>%
    filter(term == "focal_z", data_subset == "primary") %>%
    transmute(hill, index, spec, resp = "Coverage / Cmax",
              est = estimate, lo = conf_low, hi = conf_high),
  read_csv("Derived/Excels/Incidence_response_summaries.csv", show_col_types = FALSE) %>%
    filter(index != "baseline", data_subset == "primary") %>%
    transmute(hill, index, spec, resp = "Incidence (m* = 6)",
              est = focal_est, lo = focal_lo, hi = focal_hi)
) %>%
  mutate(Index = factor(recode(index, !!!inc_index_lab), levels = rev(unname(inc_index_lab))),
         Hill  = factor(recode(hill, !!!hill_labels), levels = hill_order),
         Spec  = recode(spec, !!!spec_labels))

p_incidence_effects <- focal_both_responses %>%
  ggplot(aes(est, Index, colour = resp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.4) +
  scale_colour_manual(values = c("Coverage / Cmax" = "#d95f02", "Incidence (m* = 6)" = "#08519c"), name = NULL) +
  facet_grid(Spec ~ Hill) +
  labs(x = "Standardized diversification-index effect on log diversity (90% CrI)", y = NULL,
       title = "Management effect: coverage/Cmax response vs point-count-standardised response") +
  ## extra right margin so the rotated facet_grid strip labels ("Component" / "Ecoregion") are not clipped when scaled into the PDF text column
  theme(legend.position = "bottom", plot.margin = margin(6, 16, 6, 6))
ggsave("Figures/Incidence_response_effects.png", p_incidence_effects, bg = "white", width = 11, height = 7, dpi = 150)
print(p_incidence_effects)

# Figure 6: distance-to-farm cutoff sensitivity, both responses ----

### was Scripts/04g's p_dist. Focal coefficient on the primary (< 300 m) set vs every assemblage, both specs, both responses.
subset_lab <- c(primary = "Primary (< 300 m)", full = "All assemblages")
resp_pal   <- c("Coverage / Cmax" = "#d95f02", "Incidence (m* = 6)" = "#08519c")

dist_both <- bind_rows(
  Model_summaries %>%
    filter(term == "focal_z", index != "baseline") %>%
    transmute(hill, index, spec, data_subset, response = "Coverage / Cmax",
              est = estimate, lo = conf_low, hi = conf_high),
  read_csv("Derived/Excels/Incidence_response_summaries.csv", show_col_types = FALSE) %>%
    filter(index != "baseline") %>%
    transmute(hill, index, spec, data_subset, response = "Incidence (m* = 6)",
              est = focal_est, lo = focal_lo, hi = focal_hi)
) %>%
  mutate(Index = factor(recode(index, !!!inc_index_lab), levels = rev(unname(inc_index_lab))),
         Hill  = factor(recode(hill, !!!hill_labels), levels = hill_order),
         Spec  = recode(spec, climate = "Component spec", ecoregion = "Ecoregion spec"),
         Subset   = factor(recode(data_subset, !!!subset_lab), levels = unname(subset_lab)),
         Response = factor(response, levels = names(resp_pal)))

p_dist_sensitivity <- ggplot(dist_both, aes(est, Index, colour = Response, shape = Subset,
                                            group = interaction(Response, Subset))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), size = 0.35) +
  scale_colour_manual(values = resp_pal, name = NULL) +
  scale_shape_manual(values = c("Primary (< 300 m)" = 16, "All assemblages" = 1), name = NULL) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1)) +
  facet_grid(Spec ~ Hill) +
  labs(x = "Standardized diversification-index effect on log diversity (90% CrI)", y = NULL,
       title = "Distance-to-farm cutoff sensitivity, both responses",
       subtitle = "Point counts averaging < 300 m from the farm (filled) vs all assemblages (open).") +
  ## extra right margin so the rotated facet_grid strip labels ("Component spec" / "Ecoregion spec") are not clipped when scaled into the PDF text column
  theme(legend.position = "bottom", panel.spacing.x = unit(1, "lines"), plot.margin = margin(6, 16, 6, 6))
ggsave("Figures/Dist_sensitivity_combined.png", p_dist_sensitivity, bg = "white", width = 11, height = 7.5, dpi = 150)
print(p_dist_sensitivity)

# Figure 7: Piedemonte-only cut ----

### was Scripts/05c. Focal-index forest plot for the single-ecoregion analysis (Scripts/04b piedemonte section): both responses, primary (< 300 m) + all-assemblages, full random effects.
pied_resp_lab <- c(cmax = "Coverage / Cmax", incidence = "Incidence (m* = 6)")

Pied_summaries <- read_csv("Derived/Excels/Farm_mgmt_piedemonte_summaries.csv", show_col_types = FALSE) %>%
  filter(term == "focal_z", re_spec == "full_re") %>%
  { if ("response_type" %in% names(.)) . else mutate(., response_type = "cmax") } %>%
  mutate(fit = if_else(data_subset == "full", "All assemblages", "Primary (< 300 m)"),
         Hill = factor(recode(hill, !!!hill_labels), levels = hill_order),
         Index = factor(recode(index, !!!inc_index_lab), levels = unname(inc_index_lab)),
         Response = factor(recode(response_type, !!!pied_resp_lab), levels = unname(pied_resp_lab)),
         fit = factor(fit, levels = c("Primary (< 300 m)", "All assemblages")))

p_pied_effects <- ggplot(Pied_summaries, aes(estimate, fct_rev(Index), colour = fit, shape = fit)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high), position = position_dodge(width = 0.6), fatten = 2) +
  facet_grid(Response ~ Hill) +
  scale_colour_manual(values = c("Primary (< 300 m)" = "#d95f02", "All assemblages" = "grey30"), name = NULL) +
  scale_shape_manual(values = c(17, 16), name = NULL) +
  labs(x = "Standardized effect on log diversity (posterior median, 90% CrI)", y = NULL,
       title = "Piedemonte-only: management diversification vs bird diversity",
       subtitle = "Region fixed, species pool ~ constant; env_pc1 + sampling adjusted.\nBoth responses; primary (< 300 m) vs all assemblages.") +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 10))
ggsave("Figures/Farm_mgmt_piedemonte.png", p_pied_effects, bg = "white", width = 10, height = 6.5)
print(p_pied_effects)

# Figure 8: forest-only interaction -- does the slope differ by survey scope? ----

### Scripts/04b forest_only section, `forest_int` variant (Cmax response): focal * All_forest, where All_forest = 1 for assemblages whose point counts were all in forest (the Gaica 2013 / 2016-17 surveys and their UBC / UBC-Gaica re-surveys) and 0 for those that sampled the land-use gradient. Main focal term = gradient slope; focal:All_forest = forest-only minus gradient; forest-only slope = their sum (Derived/Excels/Forest_int_slopes.csv). The land-cover fits (silvopasture % + edge density together) also include the 3 reference farms.
scope_labels <- c(`0` = "Gradient-surveyed", `1` = "Forest-only")
scope_colours <- c(`Gradient-surveyed` = "#d95f02", `Forest-only` = "#1b9e77")
fo_term_labels <- c(index_labels[-1], ssp_z = "Silvopasture %", ed_z = "Edge density")
fo_term_order <- unname(fo_term_labels)

## Slopes: gradient, forest-only, and their difference, every focal predictor x spec x Hill
Fo_int_slopes <- read_csv("Derived/Excels/Forest_int_slopes.csv", show_col_types = FALSE) %>%
  mutate(Predictor = factor(if_else(focal == "landcover", fo_term_labels[term], fo_term_labels[focal]), levels = fo_term_order),
         Hill = factor(recode(hill, !!!hill_labels), levels = hill_order),
         Spec = recode(spec, !!!spec_labels),
         Slope = factor(recode(slope, gradient = "Gradient-surveyed", forest_only = "Forest-only", difference = "Difference (forest-only - gradient)"),
                        levels = c("Gradient-surveyed", "Forest-only", "Difference (forest-only - gradient)")))

p_fo_int_slopes <- ggplot(Fo_int_slopes, aes(estimate, fct_rev(Predictor), colour = Slope, shape = Slope)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high), position = position_dodge(width = 0.7), fatten = 2) +
  facet_grid(Spec ~ Hill) +
  scale_colour_manual(values = c(scope_colours, `Difference (forest-only - gradient)` = "grey30"), name = NULL) +
  scale_shape_manual(values = c(16, 16, 4), name = NULL) +
  labs(x = "Standardized effect on log diversity (posterior median, 90% CrI)", y = NULL,
       title = "Management / land-cover slope by survey scope (focal x All_forest interaction)",
       subtitle = "Cmax response; each index fitted alone, silvopasture % + edge density fitted together (+ reference farms)") +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 10))
ggsave("Figures/Forest_int_slopes.png", p_fo_int_slopes, bg = "white", width = 11, height = 6.5)
print(p_fo_int_slopes)

## Conditional effects: predicted diversity along the focal predictor, one line per survey scope (other predictors -- including the other land-cover metric -- at their means, Ecoregion at its reference level, random effects excluded), with the assemblages each fit used; the land-cover fits include the 3 reference farms (all gradient-surveyed), drawn as labelled triangles
fo_int_files <- list.files("Derived/models", pattern = "^forestonly__forest_int__.*\\.rds$", full.names = TRUE)

if (length(fo_int_files) > 0) {
  fo_int_fits <- set_names(map(fo_int_files, readRDS), str_remove(basename(fo_int_files), "\\.rds$"))

  ## land-cover back-transform: 04b z-scores silvopasture % and edge density across the 04a frame's assemblages
  Lc_frame <- read_csv("Derived/Excels/Landcover_assemblage_metrics.csv", show_col_types = FALSE) %>%
    semi_join(Model_data, by = "Assemblage")
  lc_scaling <- tibble(index = c("ssp_z", "ed_z"),
                       mean = c(mean(Lc_frame$ssp, na.rm = TRUE), mean(Lc_frame$ed, na.rm = TRUE)),
                       sd = c(sd(Lc_frame$ssp, na.rm = TRUE), sd(Lc_frame$ed, na.rm = TRUE)))
  fo_scaling <- bind_rows(index_scaling, lc_scaling)

  ## one fit x one focal term (an index fit has focal_z; a land-cover fit has ssp_z and ed_z)
  fo_int_conditional <- function(key, term) {
    parts <- str_match(key, "^forestonly__forest_int__(richness|shannon|simpson)__(.+)__(ecoregion|climate)$")
    focal <- parts[3]
    sc <- fo_scaling %>% filter(index == if (focal == "landcover") term else focal)
    fit <- fo_int_fits[[key]]
    labels <- tibble(Hill = factor(hill_labels[[parts[2]]], levels = hill_order),
                     Predictor = factor(fo_term_labels[[if (focal == "landcover") term else focal]], levels = fo_term_order),
                     Spec = spec_labels[[parts[4]]])
    lines <- as_tibble(conditional_effects(fit, effects = paste0(term, ":All_forest"),
                                           int_conditions = list(All_forest = c(0, 1)))[[1]]) %>%
      mutate(x = .data[[term]] * sc$sd + sc$mean, Scope = scope_labels[as.character(All_forest)],
             diversity = exp(estimate__), lower = exp(lower__), upper = exp(upper__)) %>%
      select(x, Scope, diversity, lower, upper) %>%
      bind_cols(labels)
    points <- as_tibble(fit$data) %>%
      mutate(x = .data[[term]] * sc$sd + sc$mean, Scope = scope_labels[as.character(All_forest)],
             diversity = exp(log_response),
             Farm = if_else(str_detect(Id_scr, "^Ref_"), "Reference farm", "SCR farm"),
             ## "Ref_El hatico" + "Cipav_mbd_16_17" -> "Hatico 16-17" (El Hatico was surveyed twice)
             ref_label = if_else(Farm == "Reference farm",
                                 paste(str_to_title(str_remove(str_remove(Id_scr, "^Ref_"), "^El ")),
                                       str_replace(str_extract(CollectorXyear, "\\d+_\\d+$"), "_", "-")),
                                 NA_character_)) %>%
      select(x, Scope, diversity, Farm, ref_label) %>%
      bind_cols(labels)
    list(lines = lines, points = points)
  }

  fo_int_terms <- tibble(key = names(fo_int_fits)) %>%
    mutate(term = map(key, \(k) if (str_detect(k, "__landcover__")) c("ssp_z", "ed_z") else "focal_z")) %>%
    unnest(term)
  Fo_int_ce <- map2(fo_int_terms$key, fo_int_terms$term, fo_int_conditional)
  Fo_int_lines <- map(Fo_int_ce, "lines") %>% list_rbind()
  Fo_int_points <- map(Fo_int_ce, "points") %>% list_rbind()

  plot_fo_int_conditional <- function(lines, points, facets, x_label, title) {
    ggplot(lines, aes(x, diversity, colour = Scope, fill = Scope)) +
      geom_point(data = points, aes(shape = Farm), alpha = 0.45, size = 1.3) +
      ggrepel::geom_text_repel(data = filter(points, !is.na(ref_label)), aes(label = ref_label),
                               size = 2.6, show.legend = FALSE, min.segment.length = 0, seed = 1989) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 1) +
      facets +
      scale_colour_manual(values = scope_colours, name = NULL, aesthetics = c("colour", "fill")) +
      scale_shape_manual(values = c(`SCR farm` = 16, `Reference farm` = 17), name = NULL) +
      labs(x = x_label, y = "Hill-number diversity", title = title,
           subtitle = "Lines: posterior median + 90% CrI, other predictors at their means (Ecoregion at its reference level), random effects excluded.\nPoints: assemblage estimates (unadjusted).") +
      theme(legend.position = "bottom", plot.subtitle = element_text(size = 9), panel.spacing.x = unit(1.5, "lines"))
  }

  p_fo_int_ssp <- plot_fo_int_conditional(
    filter(Fo_int_lines, Predictor == "Silvopasture %"), filter(Fo_int_points, Predictor == "Silvopasture %"),
    facet_grid(Hill ~ Spec, scales = "free_y"), "Silvopasture (% of the 300 m point-count footprint)",
    "Silvopasture % x survey scope (edge density held at its mean)")
  ggsave("Figures/Forest_int_ssp_conditional.png", p_fo_int_ssp, bg = "white", width = 8, height = 8)
  print(p_fo_int_ssp)

  p_fo_int_ed <- plot_fo_int_conditional(
    filter(Fo_int_lines, Predictor == "Edge density"), filter(Fo_int_points, Predictor == "Edge density"),
    facet_grid(Hill ~ Spec, scales = "free_y"), "Edge density (m / ha of the 300 m point-count footprint)",
    "Edge density x survey scope (silvopasture % held at its mean)")
  ggsave("Figures/Forest_int_ed_conditional.png", p_fo_int_ed, bg = "white", width = 8, height = 8)
  print(p_fo_int_ed)

  p_fo_int_index_component <- plot_fo_int_conditional(
    filter(Fo_int_lines, Spec == "Component", Predictor %in% index_labels), filter(Fo_int_points, Spec == "Component", Predictor %in% index_labels),
    facet_grid(Hill ~ Predictor, scales = "free"), "Diversification index [0-1]",
    "Diversification index x survey scope -- Component adjustment set")
  ggsave("Figures/Forest_int_index_conditional_component.png", p_fo_int_index_component, bg = "white", width = 11, height = 8)
  print(p_fo_int_index_component)

  p_fo_int_index_ecoregion <- plot_fo_int_conditional(
    filter(Fo_int_lines, Spec == "Ecoregion", Predictor %in% index_labels), filter(Fo_int_points, Spec == "Ecoregion", Predictor %in% index_labels),
    facet_grid(Hill ~ Predictor, scales = "free"), "Diversification index [0-1]",
    "Diversification index x survey scope -- Ecoregion adjustment set")
  ggsave("Figures/Forest_int_index_conditional_ecoregion.png", p_fo_int_index_ecoregion, bg = "white", width = 11, height = 8)
  print(p_fo_int_index_ecoregion)
} else {
  message("05: no forest_int fits in Derived/models/ -- skipping the forest-only conditional-effect figures.")
}
