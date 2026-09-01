# Figure for the Piedemonte-only farm-management models ----

### Reads Scripts/04c_farm_mgmt_piedemonte.R's summary table and draws the focal-index forest plot (primary fit + the all-assemblages sensitivity). Safe to source() -- never refits.

### The drop-CollectorXyear-RE sensitivity is fit by 04c but not shown here -- it moves the focal coefficient by <= 0.01 (see the 04c header) and the extra series clutters the plot.

# Setup ----
library(tidyverse)
library(cowplot)
ggplot2::theme_set(theme_cowplot())
dir.create("Figures", showWarnings = FALSE)

index_labels <- c(Land_use_div = "Land use", Water_mgmt_div = "Water management",
                  Pasture_mgmt_div = "Pasture management", All_practices_div = "All practices")
index_order <- unname(index_labels)
hill_labels <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)", simpson = "Simpson (q = 2)")

resp_labels <- c(cmax = "Coverage / Cmax", incidence = "Incidence (m* = 6)")

Summ <- read_csv("Derived/Excels/Farm_mgmt_piedemonte_summaries.csv", show_col_types = FALSE) %>%
  filter(term == "focal_z", re_spec == "full_re") %>%
  ## a summaries file written before the two-response split has no response_type column
  { if ("response_type" %in% names(.)) . else mutate(., response_type = "cmax") } %>%
  mutate(
    fit = if_else(data_subset == "full", "All assemblages", "Primary (< 300 m)"),
    Hill = factor(recode(hill, !!!hill_labels), levels = unname(hill_labels)),
    Index = factor(recode(index, !!!index_labels), levels = index_order),
    Response = factor(recode(response_type, !!!resp_labels), levels = unname(resp_labels)),
    fit = factor(fit, levels = c("Primary (< 300 m)", "All assemblages"))
  )

# Forest plot ----

p_pied_effects <- ggplot(Summ, aes(estimate, fct_rev(Index), colour = fit, shape = fit)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf_low, xmax = conf_high),
                  position = position_dodge(width = 0.6), fatten = 2) +
  facet_grid(Response ~ Hill) +
  scale_colour_manual(values = c("Primary (< 300 m)" = "#d95f02",
                                 "All assemblages" = "grey30"), name = NULL) +
  scale_shape_manual(values = c(17, 16), name = NULL) +
  labs(x = "Standardized effect on log diversity (posterior median, 90% CrI)", y = NULL,
       title = "Piedemonte-only: management diversification vs bird diversity",
       subtitle = "Region fixed, species pool ~ constant; env_pc1 + sampling adjusted.\nBoth responses; primary (< 300 m) vs all assemblages. n = 17-24 farms per index.") +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 10))

ggsave("Figures/Farm_mgmt_piedemonte.png", p_pied_effects, bg = "white", width = 10, height = 6.5)
print(p_pied_effects)
