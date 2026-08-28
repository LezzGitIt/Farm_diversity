# Exploratory plots of the farm management diversification indices ----

### Two groups of plots: (1) the four indices on their own -- distributions, mutual correlation, correlation with farm covariates, and a matched-vs-unmatched representativeness check; (2) how the indices vary among the 5 ecoregions, contrasted with how bird diversity varies, to gauge how far ecoregion confounds any management-bird-diversity relationship.

# Setup ----
library(tidyverse)
library(GGally)
library(cowplot)
library(conflicted)

## Pin the dplyr verbs so this script still runs when sourced into a session that already has MASS / car attached (e.g. after 03_Farm_diversity.qmd)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::recode)

ggplot2::theme_set(theme_cowplot())

dir.create("Figures", showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

# Load data ----

## Matched farm diversification indices (only farms with an associated bird biodiversity estimate; see Scripts/02_match_farm_diversity.R) and the un-matched farms, kept for the representativeness check below. The matched file now carries Ecoregion and the farm-level climate covariates (added in 01 from Data/Farm_covariates.csv).
Farm_div_matched <- read_csv("Derived/Excels/Farm_diversity_matched.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))
Farm_div_unmatched <- read_csv("Derived/Excels/Farm_diversity_unmatched.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))

Farm_div_eco <- Farm_div_matched

stopifnot(sum(is.na(Farm_div_eco$Ecoregion)) == 0)

# Labels ----

div_index_names <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")
div_index_labels <- c(
  Land_use_div = "Land use",
  Water_mgmt_div = "Water management",
  Pasture_mgmt_div = "Pasture management",
  All_practices_div = "All practices"
)

bird_metric_names <- c("Richness_mean", "Shannon_mean", "Simpson_mean")
bird_metric_labels <- c(
  Richness_mean = "Richness (q = 0)",
  Shannon_mean = "Shannon (q = 1)",
  Simpson_mean = "Simpson (q = 2)"
)

## Farm-level environmental covariates for the by-ecoregion comparison. Canopy / elevation / biomass are the MJE per-farm medians; temperature and precipitation are the point-count means from Data/Farm_covariates.csv (joined in 01). Temperature is included to show it is ~ -elevation (r about -0.99); precipitation is a more independent axis.
covar_names <- c("WVCC_median", "DEM_median", "Biomasa_median", "Avg_temp_mean", "Tot_prec_mean")
covar_labels <- c(
  WVCC_median = "Canopy cover (WVCC)",
  DEM_median = "Elevation (DEM)",
  Biomasa_median = "Biomass",
  Avg_temp_mean = "Mean temperature",
  Tot_prec_mean = "Total precipitation"
)

## Panel order: specific practice groups first, overall index last
div_level_order <- unname(div_index_labels)

# Distribution of the four diversity indices ----

p_distributions <- Farm_div_matched %>%
  select(Id_gcs, all_of(div_index_names)) %>%
  pivot_longer(cols = all_of(div_index_names), names_to = "Index", values_to = "Value") %>%
  mutate(Index = factor(recode(Index, !!!div_index_labels), levels = div_level_order)) %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 15, fill = "#2a78d6", color = "white") +
  facet_wrap(~Index) +
  labs(
    x = "Diversification index [0-1]", y = "Number of farms",
    title = "Distribution of farm management diversification indices",
    subtitle = paste0("n = ", nrow(Farm_div_matched), " farms with an associated bird biodiversity estimate")
  )
ggsave("Figures/Farm_diversity_index_distributions.png", p_distributions, bg = "white", width = 9, height = 7)
print(p_distributions)

# Correlation among diversity indices and farm-level environmental covariates ----

cor_vars <- c(div_index_names, "dist_predio_cercano", "WVCC_mean", "DEM_mean", "Biomasa_mean",
              "Avg_temp_mean", "Tot_prec_mean")
cor_labels <- c(
  div_index_labels,
  dist_predio_cercano = "Dist. nearest\nSCR farm",
  WVCC_mean = "Canopy\ncover",
  DEM_mean = "Elevation",
  Biomasa_mean = "Biomass",
  Avg_temp_mean = "Temperature",
  Tot_prec_mean = "Precipitation"
)

p_corr <- Farm_div_matched %>%
  select(all_of(cor_vars)) %>%
  rename(!!!setNames(names(cor_labels), unname(cor_labels))) %>%
  ggcorr(label = TRUE, label_size = 3, label_round = 2, hjust = 0.75, size = 3.5, layout.exp = 2) +
  labs(title = "Correlation among diversification indices and farm covariates")
ggsave("Figures/Farm_diversity_correlations.png", p_corr, bg = "white", width = 8, height = 7)
print(p_corr)

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
print(p_pairs)

# Representativeness: are the indices similar for farms dropped due to no biodiversity match? ----

p_matched_compare <- bind_rows(
  Farm_div_matched %>% mutate(Match_status = "Matched"),
  Farm_div_unmatched %>% mutate(Match_status = "No biodiversity estimate")
) %>%
  select(Id_gcs, Match_status, all_of(div_index_names)) %>%
  pivot_longer(cols = all_of(div_index_names), names_to = "Index", values_to = "Value") %>%
  mutate(Index = factor(recode(Index, !!!div_index_labels), levels = div_level_order)) %>%
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
print(p_matched_compare)

# By ecoregion: ordering and long format ----

### Motivation: at this spatial scale bird diversity is mostly explained by ecoregion (see Scripts/qmd/03_Farm_diversity.qmd). If the management indices also track ecoregion, then Ecoregion in a bird-diversity model soaks up variation that farm management could otherwise explain.

## Order ecoregions by ascending mean farm-level bird richness, so every plot reads left (low-diversity ecoregion) to right (high) and a management gradient, if present, shows as a trend
ecoregion_order <- Farm_div_eco %>%
  summarize(Richness_mean = mean(Richness_mean, na.rm = TRUE), .by = Ecoregion) %>%
  arrange(Richness_mean) %>%
  pull(Ecoregion)

Farm_div_eco <- Farm_div_eco %>%
  mutate(Ecoregion = factor(Ecoregion, levels = ecoregion_order))

Div_long <- Farm_div_eco %>%
  select(Id_gcs, Ecoregion, all_of(div_index_names)) %>%
  pivot_longer(all_of(div_index_names), names_to = "Index", values_to = "Value") %>%
  mutate(Index = factor(recode(Index, !!!div_index_labels), levels = div_level_order))

Bird_long <- Farm_div_eco %>%
  select(Id_gcs, Ecoregion, all_of(bird_metric_names)) %>%
  pivot_longer(all_of(bird_metric_names), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(recode(Metric, !!!bird_metric_labels), levels = unname(bird_metric_labels)))

Covar_long <- Farm_div_eco %>%
  select(Id_gcs, Ecoregion, all_of(covar_names)) %>%
  pivot_longer(all_of(covar_names), names_to = "Covariate", values_to = "Value") %>%
  mutate(Covariate = factor(recode(Covariate, !!!covar_labels), levels = unname(covar_labels)))

# By ecoregion: group summaries and variance explained ----

## One reusable summary so the index and bird-metric plots share identical error-bar logic
summarise_by_group <- function(df, group_col) {
  df %>%
    filter(!is.na(Value)) %>%
    summarize(
      n = dplyr::n(),
      Mean = mean(Value),
      SD = sd(Value),
      .by = c(Ecoregion, {{ group_col }})
    ) %>%
    mutate(
      SE = SD / sqrt(n),
      CI_lower = Mean - qt(0.975, pmax(n - 1, 1)) * SE,
      CI_upper = Mean + qt(0.975, pmax(n - 1, 1)) * SE
    )
}

Div_summary <- summarise_by_group(Div_long, Index)
Bird_summary <- summarise_by_group(Bird_long, Metric)
Covar_summary <- summarise_by_group(Covar_long, Covariate)

## For each index / bird metric: R-squared of `response ~ Ecoregion` (share of variation that is between-ecoregion) plus a Kruskal-Wallis test (rank-based, robust to the bounded [0-1] indices and small groups)
variance_explained_by_ecoregion <- function(df, group_col) {
  df %>%
    filter(!is.na(Value)) %>%
    nest(data = -{{ group_col }}) %>%
    mutate(
      n = map_int(data, nrow),
      R2_ecoregion = map_dbl(data, ~ summary(lm(Value ~ Ecoregion, data = .x))$r.squared),
      KW_p_value = map_dbl(data, ~ kruskal.test(Value ~ Ecoregion, data = .x)$p.value)
    ) %>%
    select(-data) %>%
    rename(Variable = {{ group_col }})
}

Ecoregion_variance <- bind_rows(
  variance_explained_by_ecoregion(Div_long, Index) %>% mutate(Type = "Management diversification"),
  variance_explained_by_ecoregion(Bird_long, Metric) %>% mutate(Type = "Bird diversity"),
  variance_explained_by_ecoregion(Covar_long, Covariate) %>% mutate(Type = "Farm covariate")
) %>%
  relocate(Type) %>%
  mutate(across(c(R2_ecoregion, KW_p_value), ~ round(.x, 3)))

print(Ecoregion_variance)
write_csv(Ecoregion_variance, "Derived/Excels/Ecoregion_variance_explained.csv")

# By ecoregion: shared plot elements ----

ecoregion_fill <- scale_fill_brewer(palette = "Dark2", name = "Ecoregion")
tilted_x <- theme(axis.text.x = element_text(angle = 30, hjust = 1))
roomy_margin <- theme(plot.margin = margin(5.5, 5.5, 5.5, 28))
index_y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.25))

# By ecoregion: diversification indices (bar + CI + jittered farms) ----

## R-squared label per facet, positioned near the top of the [0-1] panel
div_r2_labels <- Ecoregion_variance %>%
  filter(Type == "Management diversification") %>%
  mutate(Index = factor(Variable, levels = div_level_order),
         label = paste0("R2(ecoregion) = ", sprintf("%.2f", R2_ecoregion)))

p_div_bars <- ggplot(Div_summary, aes(x = Ecoregion, y = Mean, fill = Ecoregion)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(CI_lower, 0), ymax = pmin(CI_upper, 1)), width = 0.25) +
  geom_jitter(
    data = Div_long %>% filter(!is.na(Value)),
    aes(x = Ecoregion, y = Value), inherit.aes = FALSE,
    width = 0.12, height = 0, alpha = 0.35, size = 1
  ) +
  geom_text(aes(y = -0.06, label = paste0("n=", n)), size = 2.7, color = "grey30") +
  geom_text(
    data = div_r2_labels, aes(x = 0.6, y = 1.04, label = label),
    inherit.aes = FALSE, hjust = 0, size = 3, fontface = "italic"
  ) +
  facet_wrap(~Index) +
  ecoregion_fill +
  index_y_scale +
  coord_cartesian(ylim = c(-0.08, 1.06), clip = "off") +
  labs(
    x = NULL, y = "Diversification index [0-1]",
    title = "Farm management diversification by ecoregion",
    subtitle = "Ecoregions ordered by ascending bird richness. Bars = mean, whiskers = 95% CI, points = farms."
  ) +
  theme(legend.position = "none") +
  tilted_x +
  roomy_margin
ggsave("Figures/Diversification_by_ecoregion_bars.png", p_div_bars, bg = "white", width = 10, height = 8)
print(p_div_bars)

# By ecoregion: diversification indices (boxplot view of the spread) ----

p_div_box <- Div_long %>%
  filter(!is.na(Value)) %>%
  ggplot(aes(x = Ecoregion, y = Value, fill = Ecoregion)) +
  geom_boxplot(outliers = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.4, size = 1) +
  facet_wrap(~Index) +
  ecoregion_fill +
  index_y_scale +
  labs(
    x = NULL, y = "Diversification index [0-1]",
    title = "Farm management diversification by ecoregion (distribution)"
  ) +
  theme(legend.position = "none") +
  tilted_x +
  roomy_margin
ggsave("Figures/Diversification_by_ecoregion_box.png", p_div_box, bg = "white", width = 10, height = 8)
print(p_div_box)

# By ecoregion: bird diversity on the same farms (for the confounding comparison) ----

p_bird_bars <- ggplot(Bird_summary, aes(x = Ecoregion, y = Mean, fill = Ecoregion)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.25) +
  geom_jitter(
    data = Bird_long %>% filter(!is.na(Value)),
    aes(x = Ecoregion, y = Value), inherit.aes = FALSE,
    width = 0.12, height = 0, alpha = 0.35, size = 1
  ) +
  facet_wrap(~Metric, scales = "free_y") +
  ecoregion_fill +
  labs(
    x = NULL, y = "Bird diversity (farm mean)",
    title = "Bird diversity by ecoregion, on the farms with a management index",
    subtitle = "Same ecoregion order as above; compare the shape of this pattern with the management pattern"
  ) +
  theme(legend.position = "none") +
  tilted_x +
  roomy_margin
ggsave("Figures/Bird_diversity_by_ecoregion.png", p_bird_bars, bg = "white", width = 10, height = 4.5)
print(p_bird_bars)

# By ecoregion: farm environmental covariates (canopy cover, elevation, biomass, temperature, precipitation) ----

## These are the mechanisms Ecoregion partly stands in for; showing they are also ecoregion-structured motivates swapping Ecoregion for climate/topography in the linking model
p_covar_bars <- ggplot(Covar_summary, aes(x = Ecoregion, y = Mean, fill = Ecoregion)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.25) +
  geom_jitter(
    data = Covar_long %>% filter(!is.na(Value)),
    aes(x = Ecoregion, y = Value), inherit.aes = FALSE,
    width = 0.12, height = 0, alpha = 0.35, size = 1
  ) +
  facet_wrap(~Covariate, scales = "free_y") +
  ecoregion_fill +
  labs(
    x = NULL, y = "Farm value",
    title = "Farm environmental covariates by ecoregion",
    subtitle = "Per-farm median canopy cover / elevation / biomass (MJE); point-count-mean temperature and precipitation"
  ) +
  theme(legend.position = "none") +
  tilted_x +
  roomy_margin
ggsave("Figures/Covariates_by_ecoregion.png", p_covar_bars, bg = "white", width = 11, height = 8)
print(p_covar_bars)

# By ecoregion: standardized signal, management vs birds vs covariates side by side ----

## z-score each variable across farms, average by ecoregion, and show management indices, bird metrics and farm covariates as panels of lines on a common axis; parallel up/down movement across ecoregions = the two are confounded
Standardized_means <- bind_rows(
  Div_long %>% rename(Variable = Index) %>% mutate(Type = "Management diversification"),
  Bird_long %>% rename(Variable = Metric) %>% mutate(Type = "Bird diversity"),
  Covar_long %>% rename(Variable = Covariate) %>% mutate(Type = "Farm covariate")
) %>%
  filter(!is.na(Value)) %>%
  mutate(Value_z = as.numeric(scale(Value)), .by = Variable) %>%
  summarize(Mean_z = mean(Value_z), .by = c(Type, Variable, Ecoregion)) %>%
  mutate(Type = factor(Type, levels = c("Management diversification", "Bird diversity", "Farm covariate")))

## One line per variable, coloured by family; each family gets its own panel so lines do not overplot and the panel-to-panel shape comparison is the whole point
p_std <- ggplot(Standardized_means, aes(x = Ecoregion, y = Mean_z, group = Variable, color = Type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~Type) +
  scale_color_manual(values = c(
    "Management diversification" = "#1b7837",
    "Bird diversity" = "#762a83",
    "Farm covariate" = "#b35806"
  )) +
  labs(
    x = NULL, y = "Ecoregion mean (z-scored across farms)",
    title = "Ecoregion signal: management diversification vs bird diversity vs farm covariates",
    subtitle = "Ecoregions ordered by ascending bird richness. One line per variable. Parallel shapes across panels = confounding."
  ) +
  theme(legend.position = "none") +
  tilted_x
ggsave("Figures/Ecoregion_signal_management_vs_birds.png", p_std, bg = "white", width = 12, height = 5)
print(p_std)

# Console summary ----

cat("\nFarms per ecoregion (matched set), ordered by ascending bird richness:\n")
Farm_div_eco %>% count(Ecoregion) %>% print()

cat("\nShare of each variable's variation that is between-ecoregion (R^2 of ~ Ecoregion):\n")
Ecoregion_variance %>% arrange(desc(R2_ecoregion)) %>% print()
