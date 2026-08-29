# Figure for the Piedemonte-only farm-management models ----

### Reads Scripts/04c_farm_mgmt_piedemonte.R's summary table and draws the focal-index forest plot (primary fit + the two sensitivities). Safe to source() -- never refits.

# Setup ----
library(tidyverse)
library(cowplot)
ggplot2::theme_set(theme_cowplot())
dir.create("Figures", showWarnings = FALSE)

index_labels <- c(Land_use_div = "Land use", Water_mgmt_div = "Water management",
                  Pasture_mgmt_div = "Pasture management", All_practices_div = "All practices")
index_order <- unname(index_labels)
hill_labels <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)", simpson = "Simpson (q = 2)")

Summ <- read_csv("Derived/Excels/Farm_mgmt_piedemonte_summaries.csv", show_col_types = FALSE) %>%
  filter(term == "focal_z") %>%
  mutate(
    fit = case_when(data_subset == "full" ~ "All assemblages",
                    re_spec == "no_collector" ~ "No collector RE",
                    TRUE ~ "Primary (< 300 m)"),
    Hill = factor(recode(hill, !!!hill_labels), levels = unname(hill_labels)),
    Index = factor(recode(index, !!!index_labels), levels = index_order),
    fit = factor(fit, levels = c("Primary (< 300 m)", "All assemblages", "No collector RE"))
  )

# Forest plot ----

p_pied_effects <- ggplot(Summ, aes(estimate, fct_rev(Index), colour = fit, shape = fit)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high),
                  position = position_dodge(width = 0.6), fatten = 2) +
  facet_wrap(~Hill) +
  scale_colour_manual(values = c("Primary (< 300 m)" = "#d95f02",
                                 "All assemblages" = "grey30",
                                 "No collector RE" = "#1b7837"), name = NULL) +
  scale_shape_manual(values = c(17, 16, 15), name = NULL) +
  labs(x = "Standardized effect on log diversity (posterior median, 90% CrI)", y = NULL,
       title = "Piedemonte-only: management diversification vs bird diversity",
       subtitle = "Region fixed, species pool ~ constant; env_pc1 + sampling adjusted. n = 17-24 farms per index") +
  theme(legend.position = "bottom")

ggsave("Figures/Farm_mgmt_piedemonte.png", p_pied_effects, bg = "white", width = 10, height = 4.2)
print(p_pied_effects)
