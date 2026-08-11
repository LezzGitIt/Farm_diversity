# Summary plots exploring the farm management diversity metrics ----

# Setup ----
library(tidyverse)
library(GGally)
library(cowplot)

ggplot2::theme_set(theme_cowplot())

dir.create("Figures", showWarnings = FALSE)

# Load data ----

## Matched farm diversity metrics (only farms with an associated bird biodiversity estimate; see
## Scripts/01_Match_farm_diversity.R) and the un-matched farms, kept for the representativeness check below
Farm_div_matched <- read_csv("Derived/Excels/Farm_diversity_matched.csv", show_col_types = FALSE)
Farm_div_unmatched <- read_csv("Derived/Excels/Farm_diversity_unmatched.csv", show_col_types = FALSE)

div_index_names <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")
div_index_labels <- c(
  Land_use_div = "Land use",
  Water_mgmt_div = "Water management",
  Pasture_mgmt_div = "Pasture management",
  All_practices_div = "All practices"
)

# Distribution of the four diversity indices ----

p_distributions <- Farm_div_matched %>%
  select(Id_gcs, all_of(div_index_names)) %>%
  pivot_longer(cols = all_of(div_index_names), names_to = "Index", values_to = "Value") %>%
  mutate(Index = recode(Index, !!!div_index_labels)) %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 15, fill = "#2a78d6", color = "white") +
  facet_wrap(~Index) +
  labs(
    x = "Diversification index [0-1]", y = "Number of farms",
    title = "Distribution of farm management diversification indices",
    subtitle = paste0("n = ", nrow(Farm_div_matched), " farms with an associated bird biodiversity estimate")
  )
ggsave("Figures/Farm_diversity_index_distributions.png", p_distributions, bg = "white", width = 9, height = 7)

# Correlation among diversity indices and farm-level environmental covariates ----

cor_vars <- c(div_index_names, "dist_predio_cercano", "WVCC_mean", "DEM_mean", "Biomasa_mean")
cor_labels <- c(
  div_index_labels,
  dist_predio_cercano = "Dist. nearest\nSCR farm",
  WVCC_mean = "Canopy\ncover",
  DEM_mean = "Elevation",
  Biomasa_mean = "Biomass"
)

p_corr <- Farm_div_matched %>%
  select(all_of(cor_vars)) %>%
  rename(!!!setNames(names(cor_labels), unname(cor_labels))) %>%
  ggcorr(label = TRUE, label_size = 3, label_round = 2, hjust = 0.75, size = 3.5, layout.exp = 2) +
  labs(title = "Correlation among diversification indices and farm covariates")
ggsave("Figures/Farm_diversity_correlations.png", p_corr, bg = "white", width = 8, height = 7)

# Pairwise relationships among the four diversity indices ----

p_pairs <- Farm_div_matched %>%
  select(all_of(div_index_names)) %>%
  rename(!!!setNames(names(div_index_labels), unname(div_index_labels))) %>%
  ggpairs(
    lower = list(continuous = wrap("points", alpha = 0.6, color = "#2a78d6")),
    diag = list(continuous = wrap("densityDiag", fill = "#2a78d6", alpha = 0.4))
  ) +
  labs(title = "Pairwise relationships among diversification indices")
ggsave("Figures/Farm_diversity_index_pairs.png", p_pairs, bg = "white", width = 9, height = 9)

# Matched vs. unmatched farms: are the indices similar for farms dropped due to no biodiversity match? ----

p_matched_compare <- bind_rows(
  Farm_div_matched %>% mutate(Match_status = "Matched"),
  Farm_div_unmatched %>% mutate(Match_status = "No biodiversity estimate")
) %>%
  select(Id_gcs, Match_status, all_of(div_index_names)) %>%
  pivot_longer(cols = all_of(div_index_names), names_to = "Index", values_to = "Value") %>%
  mutate(Index = recode(Index, !!!div_index_labels)) %>%
  ggplot(aes(x = Match_status, y = Value, fill = Match_status)) +
  geom_boxplot(outliers = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  facet_wrap(~Index) +
  labs(
    x = NULL, y = "Diversification index [0-1]",
    title = "Farm diversification indices: matched vs. unmatched farms",
    subtitle = "Checking whether farms dropped for lacking a biodiversity estimate look systematically different"
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("Figures/Farm_diversity_matched_vs_unmatched.png", p_matched_compare, bg = "white", width = 8, height = 7)
