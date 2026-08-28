# Causal DAG for the bird-diversity ~ farm-management-diversification analysis ----

### Reference artifact, not a pipeline stage (like Scripts/Farm_diversity_fns.R). It encodes the assumed causal structure behind Scripts/04_farm_mgmt_mod.R (Aaron's hand-drawn DAG, 2026-08-28), derives the adjustment set the model should target, checks the implied independencies against the data, and maps the 04 specs onto it. Consult before trusting a 04 coefficient. Produces Figures/DAG.png and Derived/Excels/DAG_adjustment_sets.csv.

# Setup ----

library(tidyverse)
library(dagitty)
library(ggdag)
library(cowplot)

ggplot2::theme_set(theme_dag())
dir.create("Figures", showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

# Variable glossary ----

### Every node, what it actually is, and which column carries it. Read this before editing the DAG -- several names are easy to misread.
# FarmDiv        EXPOSURE. One of the four [0-1] farm-management diversification indices from Maria Esquivel (Land_use_div / Water_mgmt_div / Pasture_mgmt_div / All_practices_div). One model per index.
# BirdDiv        OUTCOME. The iNEXT Hill-number diversity estimate for the assemblage (q = 0/1/2), i.e. the *observed / estimated* diversity -- so every sampling variable points into it. Point-count duration differences across datasets are NOT a node: iNEXT's individual-based accumulation + coverage standardisation already absorb them.
# Ecoregion      5-level biogeographic region. A coarse categorical summary of Climate + LandForest + regional species pool; the SCR programme also rolled out differently by region, so it also sits upstream of FarmDiv.
# Climate        Farm-mean precipitation + elevation (confirmed: elevation IS part of this node; temperature is r = -0.99 with elevation so dropped). Two columns (Elev_mean, Tot_prec_mean), entered as poly(.,2) in 04.
# LandForest     Amount of woody-vegetation cover in the surrounding landscape (WVCC canopy at the 06b scale-of-effect radius, ~10 km). "Landscape forest" in the drawing. R^2 ~ 0.45 with Ecoregion -- NOT covered by it, ~58% of its SD is within-region.
# ForestConfig   Configuration / fragmentation of that landscape forest (connectivity, edge, patch size) -- confirmed, distinct from the amount. UNOBSERVED. Not confounding (no arrow to FarmDiv), just unmodelled outcome noise.
# FarmerValues   UNOBSERVED farmer priorities / attitudes. Drives how much a farmer diversifies management AND the land-use context they maintain. Confirmed: NO direct arrow to BirdDiv -> conditioning on LandForest fully blocks it *as drawn*. The reason LandForest must be adjusted, not just added for precision.
# NumHab         Number of DISTINCT HABITAT TYPES SURVEYED on the farm. A SAMPLING-SCOPE covariate (distinct from NumPC = effort). CRITICAL: the sampling scope changed over time -- 2013 and 2016/17 surveys covered FOREST ONLY; 2019 onward covered the land-use gradient (forest + silvopasture + pasture). So early-year assemblages have low NumHab and their diversity estimate reflects forest birds, not the whole farm. NumHab + the Observer(CollectorXyear) RE + Year together are what control for this. Not a mediator of FarmDiv.
# NumPC          Number of distinct point counts in the assemblage (per [farm x team x year], confirmed -- not per farm). Sampling effort.
# Season         Calendar timing of the surveys within the season (Nearctic-migrant influx Oct-Mar). All modelled data is "Early" season, so this is purely the within-Early calendar position. Column: mean Julian day -> doy_sin/doy_cos.
# Year           Survey year-group (Ano_grp: 13_14, 16_17, 19, 22, 24, 25_26). Fixed Year was dropped from 04 (near-collinear with Observer); folded into the Observer/CollectorXyear RE. Kept as a node because the sampling-scope change above is a Year effect.
# Observer       The data collector x year-group batch = one of 6 field teams / 8 CollectorXyear levels (Gaica-mbd, Cipav, Gaica-distancia, Unillanos, UBC-mbd, UBC-gaica). Protocol + which farms + which habitats surveyed. NOT the individual observer (not in the data). Region-bound (4 of 6 teams are Piedemonte-only). = the CollectorXyear random effect.

# The DAG ----

dag <- dagitty('dag {
  FarmDiv [exposure]
  BirdDiv [outcome]
  FarmerValues [latent]
  ForestConfig [latent]

  Ecoregion -> Climate
  Ecoregion -> FarmDiv
  Ecoregion -> LandForest
  Ecoregion -> Observer
  Ecoregion -> Season

  Climate -> FarmDiv
  Climate -> LandForest
  Climate -> BirdDiv

  FarmerValues -> FarmDiv
  FarmerValues -> LandForest

  LandForest -> BirdDiv
  LandForest -> ForestConfig
  ForestConfig -> BirdDiv

  FarmDiv -> BirdDiv

  Year -> Observer
  Year -> BirdDiv
  Year -> NumHab
  Observer -> Season
  Season -> BirdDiv
  Observer -> NumPC
  Observer -> NumHab
  Observer -> BirdDiv
  NumPC -> BirdDiv
  NumHab -> BirdDiv
}')

### Year -> NumHab: the sampling-SCOPE norm shifted around 2019 -- 2013 & 2016/17 surveys were forest-only, 2019+ spanned the land-use gradient. So NumHab (habitat types surveyed) is low in the early years by protocol, not by farm. This is a Year effect on top of the team effect (Observer -> NumHab). Controlling NumHab + the CollectorXyear RE handles it, but the early forest-only assemblages are, in effect, estimating forest-bird diversity rather than whole-farm diversity -- a caveat for interpreting those rows. See Project_notes.md "Background: the source bird dataset".

### NOT drawn, and why:
# - Ecoregion -> BirdDiv direct: omitted per Aaron\'s DAG -- Ecoregion acts on birds only THROUGH Climate + LandForest (+ the unmeasured species pool). If a residual species-pool effect exists, the climate spec under-adjusts; the ecoregion robustness spec is the guard against that.
# - FarmDiv -> LandForest: the 10 km landscape is dominated by non-farm land, so on-farm tree planting barely moves it. (At a small radius this arrow would exist.)
# - FarmDiv -> NumHab: per Aaron, NumHab is what was surveyed, not a farm property caused by diversification.

coordinates(dag) <- list(
  x = c(FarmDiv = 1.7, BirdDiv = 2.3, FarmerValues = 3.8, ForestConfig = 4.1,
        Ecoregion = 0, Climate = 1.1, LandForest = 3.3,
        Year = 0, Season = 4.2, NumPC = 0.7, NumHab = 1.5, Observer = 1.1),
  y = c(FarmDiv = 4.2, BirdDiv = 2, FarmerValues = 4.6, ForestConfig = 2.7,
        Ecoregion = 3.4, Climate = 3.1, LandForest = 3.6,
        Year = 2.2, Season = 2.2, NumPC = 0.2, NumHab = 0.2, Observer = 0.9)
)

# Adjustment set ----

### Total effect of FarmDiv on BirdDiv. NumHab is a sampling covariate here (not a mediator), so it is a normal control, not something to hold out.
total_min <- adjustmentSets(dag, exposure = "FarmDiv", outcome = "BirdDiv",
                            effect = "total", type = "minimal")
total_can <- adjustmentSets(dag, exposure = "FarmDiv", outcome = "BirdDiv",
                            effect = "total", type = "canonical")

cat("== Total-effect adjustment (minimal) ==\n"); print(total_min)
cat("\n== Total-effect adjustment (canonical) ==\n"); print(total_can)

### READING (2026-08-28, Aaron's DAG)
# dagitty gives two minimal sufficient sets:
#   (1) { Climate, LandForest, Observer, Season, Year }   -- block every Ecoregion channel to BirdDiv without conditioning on Ecoregion itself
#   (2) { Climate, Ecoregion, LandForest }                -- condition on Ecoregion instead of its Observer/Season channels
# BOTH require Climate and LandForest directly. Points:
#   - Climate is a fork with its own arrows to FarmDiv and BirdDiv, so conditioning on Ecoregion (its parent) does NOT substitute for it -- Climate has to be in the model explicitly.
#   - There is no minimal set that is "Ecoregion + sampling, no Climate". So the pure "ecoregion" spec (5-level factor, no elev/precip) is NOT dag-exact; it works only insofar as the factor PROXIES Climate (R^2 ~ 0.8 on elev/precip). It also mops up any residual species-pool effect the climate spec misses. => keep it as a ROBUSTNESS check, not the primary.
#   - A combined "Ecoregion + Climate" model is set (2) plus redundancy -- unnecessary AND unidentifiable at ~0.8 collinearity. Skip it (confirms Aaron's read).
#   - LandForest (canopy) must be ADJUSTED, not merely added for precision: it blocks the FarmerValues backdoor  FarmDiv <- FarmerValues -> LandForest -> BirdDiv. => canopy belongs in BOTH specs.
#   - NumHab / NumPC are sampling covariates (per Aaron), on the outcome only -- kept for precision, not needed to block a backdoor. Watch: NumHab ~ Land_use / All_practices sits at p ~ 0.05-0.08 (below), so for those two indices NumHab may absorb a sliver of real signal; a "drop NumHab" sensitivity for them is on the backlog.
#   - FarmerValues is FULLY blocked by conditioning on LandForest *as this DAG is drawn* (its only path to BirdDiv is via LandForest). The limitation is the assumption itself: if bird-minded farmers also do in-field things not captured by canopy or the four indices (retain snags, less pesticide, wider fencerows), FarmerValues -> BirdDiv would need its own arrow and residual confounding returns. ForestConfig (landscape configuration) is unmeasured -> its contribution to BirdDiv is unmodelled noise, not confounding (no arrow to FarmDiv).
# => 04 PRIMARY = climate spec = poly(Elev,2) + poly(Precip,2) + canopy + CollectorXyear RE + doy + NumPC + NumHab  (== minimal set 1).
#    04 ROBUSTNESS = ecoregion spec = Ecoregion + canopy + CollectorXyear RE + doy + NumPC + NumHab.
#    No combined spec. LOO for the doy / RE nuisance terms: deferred (backlog).

# Testable implications ----

ici <- impliedConditionalIndependencies(dag)
cat("\n== Implied conditional independencies (first 20) ==\n")
print(head(ici, 20))

### The checks that bear on IDENTIFICATION of the management effect: is any covariate we treat as a confounder / precision term actually associated with the exposure beyond the region+climate we condition on? (Interdependencies purely among the nuisance covariates -- Season~Observer etc. -- are not tested: they are all conditioned on, so they do not move the FarmDiv coefficient.) Small p = tension.
data_path <- "Derived/Excels/Farm_mgmt_model_data.csv"
if (file.exists(data_path)) {
  d <- read_csv(data_path, show_col_types = FALSE) %>%
    distinct(Id_gcs, .keep_all = TRUE)   # farm-level covariates: one row per farm

  p_of <- function(fit, term) tryCatch(summary(fit)$coefficients[term, "Pr(>|t|)"], error = function(e) NA_real_)
  indices <- c("Land_use_div", "Water_mgmt_div", "Pasture_mgmt_div", "All_practices_div")

  test_indep <- function(y, index, given) {
    f <- reformulate(c(index, given), response = y)
    p_of(lm(f, d), index)
  }

  checks <- expand_grid(
    target = c("canopy_10k", "Num_pc_log", "Num_hab_num"),
    index = indices
  ) %>%
    mutate(
      given = if_else(target == "canopy_10k", "Ecoregion + Elev_mean + Tot_prec_mean", "Ecoregion"),
      implication = paste0(target, " _||_ ", index, " | ", given),
      p_value = pmap_dbl(list(target, index, given), ~ test_indep(..1, ..2, ..3)),
      flag = if_else(p_value < 0.05, "** exposure-covariate tension", "")
    ) %>%
    select(implication, p_value, flag)
  cat("\n== Data spot-checks: covariate _||_ exposure | conditioning set (12 tests) ==\n")
  print(checks, width = Inf, n = 12)
  cat("canopy: a real FarmDiv->canopy arrow at 10 km is implausible -> a flag here is most likely the FarmerValues path (already blocked by conditioning on canopy) or chance. NumPC/NumHab flags would mean sampling effort tracks management beyond region -- worth a closer look if strong.\n")
}

# How the 04 specs map onto the DAG ----

spec_map <- tribble(
  ~spec,        ~adjusts,                                                              ~role,
  "climate",    "poly(Elev,2) + poly(Precip,2) + canopy + CollectorXyear RE + doy + NumPC + NumHab", "PRIMARY -- equals minimal set {Climate, LandForest, Observer, Season, Year}",
  "ecoregion",  "Ecoregion + canopy + CollectorXyear RE + doy + NumPC + NumHab",               "ROBUSTNESS -- 5-level factor proxies Climate (R2 ~ 0.8); also catches any residual species-pool effect; canopy still needed for the FarmerValues backdoor",
  "(combined)", "Ecoregion + Climate together",                                               "NOT USED -- minimal set 2 plus redundancy; unidentifiable at ~0.8 collinearity"
)
cat("\n== 04 spec <-> DAG ==\n"); print(spec_map, width = Inf)

adj_out <- tibble(
  effect = "total",
  minimal_set = map_chr(as.list(total_min), ~ paste(sort(.x), collapse = " + "))
)
write_csv(adj_out, "Derived/Excels/DAG_adjustment_sets.csv")

# Figure ----

dag_tidy <- dag %>%
  tidy_dagitty() %>%
  mutate(role = case_when(
    name == "FarmDiv" ~ "exposure",
    name == "BirdDiv" ~ "outcome",
    name %in% c("FarmerValues", "ForestConfig") ~ "unobserved",
    name %in% c("NumPC", "NumHab", "Season", "Year", "Observer") ~ "sampling",
    TRUE ~ "confounder / covariate"
  ))

role_cols <- c(exposure = "#1b7837", outcome = "#762a83", unobserved = "grey72",
               sampling = "#4393c3", "confounder / covariate" = "grey30")

p_dag <- ggplot(dag_tidy, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_colour = "grey55") +
  geom_dag_point(aes(colour = role), size = 21) +
  geom_dag_text(colour = "white", size = 2.3) +
  scale_colour_manual(values = role_cols, name = NULL) +
  expand_limits(x = c(-0.5, 4.7), y = c(-0.2, 5)) +
  labs(title = "Assumed causal structure: farm management diversification -> bird diversity",
       subtitle = "Minimal adjustment set: Climate + landscape forest + Observer/Season/Year (the climate spec). NumHab/NumPC = sampling covariates. FarmerValues + forest configuration unobserved.") +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 9))

ggsave("Figures/DAG.png", p_dag, width = 12, height = 8.5, bg = "white")
print(p_dag)
