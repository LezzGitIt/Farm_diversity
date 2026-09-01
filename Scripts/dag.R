# Causal DAG for the bird-diversity ~ farm-management-diversification analysis ----

### Reference artifact, not a pipeline stage (like Scripts/Model_fns.R). It encodes the assumed causal structure behind Scripts/04a_farm_mgmt_models.R (Aaron's hand-drawn DAG, 2026-08-28; simplified 2026-08-31), derives the adjustment set the model should target, checks the implied independencies against the data, and maps the 04 specs onto it. Consult before trusting a 04 coefficient. Produces Figures/DAG.png and Derived/Excels/DAG_adjustment_sets.csv.

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
# Ecoregion      5-level biogeographic region. A coarse categorical summary of Climate + LandForest + BiogeoHistory; the SCR programme also rolled out differently by region, so it also sits upstream of FarmDiv.
# Climate        Farm-mean precipitation + elevation (confirmed: elevation IS part of this node; temperature is r = -0.99 with elevation so dropped). Two columns (Elev_mean, Tot_prec_mean), entered as poly(.,2) in 04.
# BiogeoHistory  UNOBSERVED biogeographic / evolutionary history of the region (Andean uplift, Pleistocene refugia, dispersal barriers). A cause of the regional SpeciesPool that Climate does not capture. No arrow to FarmDiv except through Ecoregion.
# SpeciesPool    Regional potential species pool -- how many / which species could occur at the farm. Caused by BiogeoHistory AND Climate. Partially MEASURED: Scripts/01b_species_pool.R builds a range-map richness count (pool_point) plus range-rarity-weighted variants (pool_wes / pool_cwe); Scripts/04b (pool_blocks / endemism_pool sections) tests them. The compositional / endemism part is only partly captured. On the path Ecoregion -> BiogeoHistory -> SpeciesPool -> BirdDiv, so the climate spec (no Ecoregion factor) blocks that backdoor only via the endemism-index proxy -- this is why the Ecoregion spec stays as a DAG-valid adjustment, not just a robustness check.
# EndemismIndex  OBSERVED. Range-rarity-weighted endemism index (pool_wes = Sum 1000/sqrt(global range km^2), Scripts/01b_species_pool.R -- Ayerbe range polygons + AVONET ranges), z-scored. A measured PROXY for BiogeoHistory: the range-restricted part of the regional pool tracks biogeographic history, not contemporary climate (the raw range-map count pool_point tracks climate and is ~redundant with elev + precip). Entered in the 04a `climate` spec to help block Ecoregion -> BiogeoHistory -> SpeciesPool -> BirdDiv. PARTIAL only -- pool_wes is r ~ 0.8 with elevation, R^2 ~ 0.5 with Ecoregion, so it is not a stand-in for conditioning on SpeciesPool itself; the `ecoregion` spec stays as the check. A leaf node (no arrow to BirdDiv or FarmDiv), so adding it does not change the adjustment sets.
# LandForest     Amount of woody-vegetation cover in the surrounding *landscape* -- WVCC canopy in the ~10 km buffer (the 02b scale-of-effect radius). At that scale it is set by topography, protected areas and hundreds of landholders, NOT by any one farmer -- so it takes no arrow from FarmDiv and is a PRECISION covariate (predicts BirdDiv, reduces residual variance), not a confounder to adjust. R^2 ~ 0.45 with Ecoregion; ~58% of its SD is within-region.
# NumPC          Number of distinct point counts in the assemblage (per [farm x team x year], confirmed -- not per farm). Sampling effort.
# Season         Calendar timing of the surveys within the season (Nearctic-migrant influx Oct-Mar). All modelled data is "Early" season, so this is purely the within-Early calendar position. Column: mean Julian day -> doy_sin/doy_cos.
# Year           Survey year-group (Ano_grp: 13_14, 16_17, 19, 22, 24, 25_26). Fixed Year was dropped from 04 (near-collinear with Observer); folded into the Observer/CollectorXyear RE. Kept as a node because the sampling-scope change (see NumHab under "Variables to consider") is a Year effect.
# Observer       The data collector x year-group batch = one of 6 field teams / 8 CollectorXyear levels (Gaica-mbd, Cipav, Gaica-distancia, Unillanos, UBC-mbd, UBC-gaica). Protocol + which farms + which habitats surveyed. NOT the individual observer (not in the data). Region-bound (4 of 6 teams are Piedemonte-only). = the CollectorXyear random effect.

# The DAG ----

dag <- dagitty('dag {
  FarmDiv [exposure]
  BirdDiv [outcome]
  BiogeoHistory [latent]

  Ecoregion -> Climate
  Ecoregion -> FarmDiv
  Ecoregion -> LandForest
  Ecoregion -> Observer
  Ecoregion -> Season
  Ecoregion -> BiogeoHistory

  Climate -> LandForest
  Climate -> BirdDiv
  Climate -> SpeciesPool

  BiogeoHistory -> SpeciesPool
  BiogeoHistory -> EndemismIndex
  SpeciesPool -> BirdDiv

  LandForest -> BirdDiv

  FarmDiv -> BirdDiv

  Year -> Observer
  Year -> BirdDiv
  Observer -> Season
  Season -> BirdDiv
  Observer -> NumPC
  Observer -> BirdDiv
  NumPC -> BirdDiv
}')

### NOT drawn, and why:
# - Ecoregion -> BirdDiv direct: omitted per Aaron's DAG -- Ecoregion acts on birds only THROUGH Climate, LandForest and BiogeoHistory -> SpeciesPool. The SpeciesPool branch is explicit: it is what makes the climate spec insufficient on its own and the Ecoregion spec a DAG-valid adjustment.
# - FarmDiv -> LandForest: LandForest is the ~10 km buffer, dominated by non-farm land -- a single farmer's planting does not move it. (At a farm-scale radius this arrow would exist; that is the on-farm forest effect under "Variables to consider".)
# - Climate -> FarmDiv: removed by Aaron 2026-08-30 -- the management indices are structured by region (culture, extension networks) rather than directly by rainfall / elevation, and Ecoregion -> FarmDiv already carries that.

# Variables to consider (removed or unmeasured) ----

### The DAG above is deliberately minimal. These are variables that were dropped for simplicity, or are unmeasured, and why each could matter. The two unmeasured ones (FarmerValues, FarmSize) are why the honest reading is "consistent with no effect", not "no effect": both would bias a true effect toward the positive and neither is in the identified adjustment set.
#
# NumHab      Number of distinct habitat types surveyed on the farm. A sampling-SCOPE covariate (distinct from NumPC = effort). The sampling scope shifted around 2019: 2013 & 2016/17 surveys were FOREST-ONLY, 2019+ spanned the land-use gradient (forest + silvopasture + pasture). So early-year assemblages estimate forest-bird diversity, not whole-farm diversity. DROPPED from the DAG because it is outcome-side only (child of Year + Observer, parent of BirdDiv) -- never on a path between FarmDiv and BirdDiv, so it is in neither adjustment set -- and in the fitted 04 models its coefficient is -0.03 to -0.09 with every CrI spanning zero: NumPC + the CollectorXyear RE + Year already absorb the scope shift. The early forest-only rows remain a caveat for interpretation, not an adjustment need. (04 still carries Num_hab_z as a legacy precision term; harmless, could be dropped for consistency.)
#
# FarmSize    UNMEASURED farm area (hectares). Bigger farms have more room / units to diversify (-> FarmDiv) and get more point counts (-> NumPC); both legs are blocked by conditioning on NumPC. The unblocked concern is a DIRECT FarmSize -> BirdDiv (species-area: more area -> larger populations, more interior habitat, less edge) that NumPC does not stand in for. Measuring farm hectares would close this.
#
# FarmerValues  UNMEASURED farmer priorities / attitudes. Drives how much a farmer diversifies management (-> FarmDiv). As a cause of the exposure only, it needs no adjustment. BUT if bird-minded farmers also do in-field things the four indices miss (retain snags, less pesticide, wider fencerows), FarmerValues gains a second path to BirdDiv and becomes unblockable confounding.
#
# On-farm / boundary woody cover  The farm-scale version of LandForest (fencerow trees, live fences, scattered pasture trees, farm canopy). At a farm-scale radius, FarmDiv -> woody cover and FarmerValues -> woody cover both exist -- this is the concrete mechanism behind the FarmerValues caveat. A farm-level canopy / tree-density measure would let this path be modelled explicitly rather than left as confounding.
#
# ForestConfig  Landscape forest CONFIGURATION (edge density, patch size, connectivity) independent of total cover. Could predict BirdDiv beyond LandForest amount. Takes no arrow from FarmDiv at the 10 km scale, so it is outcome noise / precision, not confounding. Removed by Aaron 2026-08-30 for simplicity.
#
# Realized / compositional species pool  pool_point is range-map RICHNESS; the compositional & endemism part of SpeciesPool (the BiogeoHistory branch that climate proxies miss) is only partly captured by pool_wes / pool_cwe, which are themselves ~0.8 correlated with elevation. An eBird-derived realized-composition or phylogenetic-turnover measure would be a cleaner instrument for that branch.

### Two visual blocks: Ecoregion + the region / environment drivers it summarises on the left, the sampling covariates on the right, exposure far left and outcome in the middle. The Ecoregion -> Observer / -> Season edges cross the figure by necessity (field teams are region-bound).
coordinates(dag) <- list(
  x = c(FarmDiv = 0.5, BirdDiv = 11.2,
        Ecoregion = 3.3, Climate = 6.8, BiogeoHistory = 3.0, EndemismIndex = 3.2, SpeciesPool = 6.4, LandForest = 10.0,
        Year = 14.2, Observer = 14.2, Season = 14.2, NumPC = 14.2),
  y = c(FarmDiv = 4.2, BirdDiv = 4.2,
        Ecoregion = 8.4, Climate = 8.4, BiogeoHistory = 6.2, EndemismIndex = 3.2, SpeciesPool = 5.6, LandForest = 8.2,
        Year = 8.6, Observer = 6.5, Season = 4.4, NumPC = 2.3)
)

# Adjustment set ----

### Total effect of FarmDiv on BirdDiv.
total_min <- adjustmentSets(dag, exposure = "FarmDiv", outcome = "BirdDiv",
                            effect = "total", type = "minimal")
total_can <- adjustmentSets(dag, exposure = "FarmDiv", outcome = "BirdDiv",
                            effect = "total", type = "canonical")

cat("== Total-effect adjustment (minimal) ==\n"); print(total_min)
cat("\n== Total-effect adjustment (canonical) ==\n"); print(total_can)

### READING (2026-08-31, Aaron's simplified DAG -- Climate -> FarmDiv removed, so the sets are smaller than before; 2026-09-01 -- observed EndemismIndex proxy added off BiogeoHistory, adjustment sets unchanged)
# dagitty gives two minimal sufficient sets:
#   (1) { Climate, LandForest, SpeciesPool, Observer, Season, Year }   -- block every Ecoregion channel to BirdDiv without conditioning on Ecoregion itself
#   (2) { Ecoregion }                                                  -- one node closes all of them at the source
# Points:
#   - FarmDiv's ONLY parent is Ecoregion, so every backdoor runs FarmDiv <- Ecoregion -> ... -> BirdDiv. Conditioning on Ecoregion (set 2) closes all of them -- nothing else is strictly required. The Ecoregion spec is exactly the minimal adjustment.
#   - Since Climate -> FarmDiv was removed, Climate is no longer a fork on the exposure. In set (1) it is still needed -- to block FarmDiv <- Ecoregion -> Climate -> BirdDiv (and the Climate -> LandForest / -> SpeciesPool legs) without conditioning on Ecoregion -- but it is NOT needed on top of Ecoregion in set (2).
#   - Set (1) REQUIRES a SpeciesPool term: FarmDiv <- Ecoregion -> BiogeoHistory -> SpeciesPool -> BirdDiv is a backdoor Climate + LandForest do NOT block. The 04a `climate` spec carries the endemism index (pool_wes = the observed EndemismIndex proxy for BiogeoHistory) for exactly this. PARTIAL block only -- pool_wes is r ~ 0.8 with elevation, R^2 ~ 0.5 with Ecoregion -- so the Ecoregion spec is kept as the check. Scripts/04b (pool_blocks: the range-map count pool_point; endemism_pool: the range-rarity-weighted pool_cwe) confirms no pool variant moves the management coefficient.
#   - Observer + Year together handle the sampling side of set (1): conditioning on Observer blocks FarmDiv <- Ecoregion -> Observer -> {NumPC, Season, BirdDiv}, but it also opens the collider path FarmDiv <- Ecoregion -> [Observer] <- Year -> BirdDiv, so Year is then required too. In 04 both are folded into the CollectorXyear RE (Year was near-collinear with Observer).
#   - NumPC is NOT a required confounder in this DAG: its only backdoor (Ecoregion -> Observer -> NumPC -> BirdDiv) is already blocked by conditioning on Observer. It is kept in 04 for PRECISION and to guard the unmeasured FarmDiv <- FarmSize -> NumPC -> BirdDiv path (see "Variables to consider"). Canopy was ALSO pure precision in the Ecoregion spec and was dropped 2026-08-31 -- that spec is now Ecoregion alone (exactly set 2).
#   - A combined "Ecoregion + Climate + pool" model is over-specified and unidentifiable at this collinearity. Skip it.
#   - Two unmeasured paths are NOT closed by any set here (see "Variables to consider"): FarmerValues -> on-farm practices the 4 indices miss -> BirdDiv, and FarmSize -> BirdDiv directly (species-area). Both push a null toward positive, so "consistent with no effect" is the honest reading.
# => 04a `climate` spec   = poly(Elev,2) + poly(Precip,2) + canopy_10k + pool_wes + CollectorXyear RE + doy + NumPC  (= minimal set 1: Climate + LandForest + SpeciesPool [via the pool_wes proxy] + Observer/Year [via the RE] + Season [via doy]; PLUS NumPC for precision).
#    04a `ecoregion` spec = Ecoregion + CollectorXyear RE + doy + NumPC  (= minimal set 2 = {Ecoregion} exactly; PLUS doy / NumPC / RE as nuisance / precision).
#    Report both; they agree the effect is ~0 once precipitation is well controlled. LOO for the doy / RE nuisance terms: deferred (backlog).
#    (Writeup names for the two specs are still being chosen -- e.g. "Component" / "Ecoregion" -- the code ids stay `climate` / `ecoregion`.)

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
    target = c("canopy_10k", "Num_pc_log"),
    index = indices
  ) %>%
    mutate(
      given = if_else(target == "canopy_10k", "Ecoregion + Elev_mean + Tot_prec_mean", "Ecoregion"),
      implication = paste0(target, " _||_ ", index, " | ", given),
      p_value = pmap_dbl(list(target, index, given), ~ test_indep(..1, ..2, ..3)),
      flag = if_else(p_value < 0.05, "** exposure-covariate tension", "")
    ) %>%
    select(implication, p_value, flag)
  cat("\n== Data spot-checks: covariate _||_ exposure | conditioning set (8 tests) ==\n")
  print(checks, width = Inf, n = 8)
  cat("canopy: a real FarmDiv->canopy arrow at 10 km is implausible -> a flag here just means the 10 km canopy still tracks management within region (shared topography / programme rollout), not a causal path; conditioning on it is still fine. A NumPC flag is expected -- that is the FarmSize -> NumPC confounder we now adjust for.\n")
}

# How the 04 specs map onto the DAG ----

spec_map <- tribble(
  ~spec,        ~adjusts,                                                              ~role,
  "climate",    "poly(Elev,2) + poly(Precip,2) + canopy_10k + pool_wes + CollectorXyear RE + doy + NumPC", "minimal set 1, fully observed: Climate + LandForest + SpeciesPool (via the pool_wes endemism-index proxy) + Observer/Year (via the RE) + Season (via doy); PLUS NumPC for precision. The primary spec",
  "ecoregion",  "Ecoregion + CollectorXyear RE + doy + NumPC",                            "minimal set 2 = {Ecoregion} exactly -- the factor closes the Climate, LandForest and BiogeoHistory->SpeciesPool backdoors at the source. PLUS doy / NumPC / RE as nuisance / precision. Canopy dropped 2026-08-31 (pure precision here). The conservative bracket",
  "(combined)", "Ecoregion + Climate + pool together",                                   "NOT USED -- over-specified, unidentifiable at ~0.8 collinearity. The canonical adjustment set adds nothing over {Ecoregion} for the focal coefficient"
)
cat("\n== 04 spec <-> DAG ==\n"); print(spec_map, width = Inf)

adj_out <- tibble(
  effect = "total",
  minimal_set = map_chr(as.list(total_min), ~ paste(sort(.x), collapse = " + "))
)
write_csv(adj_out, "Derived/Excels/DAG_adjustment_sets.csv")

# Figure ----

### Short node names drive dagitty / the adjustment sets; these are the readable labels for the plot only.
node_labels <- c(
  FarmDiv       = "Farm \ndiversification",
  BirdDiv       = "Bird\ndiversity",
  Ecoregion     = "Ecoregion",
  Climate       = "Climate\n(elev + precip)",
  BiogeoHistory = "Biogeographic\nhistory",
  EndemismIndex = "Endemism\nindex",
  SpeciesPool   = "Regional\nspecies pool",
  LandForest    = "Landscape\nforest",
  Year          = "Survey year",
  Observer      = "Field team",
  Season        = "Survey timing",
  NumPC         = "Number\nof surveys"
)

dag_tidy <- dag %>%
  tidy_dagitty() %>%
  mutate(
    role = case_when(
      name == "FarmDiv" ~ "exposure",
      name == "BirdDiv" ~ "outcome",
      name == "BiogeoHistory" ~ "unobserved",
      name %in% c("NumPC", "Season", "Year", "Observer") ~ "sampling",
      TRUE ~ "region / environment"
    ),
    label = node_labels[name]
  )

role_cols <- c(exposure = "#1b7837", outcome = "#762a83", unobserved = "grey72",
               sampling = "#4393c3", "region / environment" = "grey30")

### Draw edges as manually shortened, gently curved segments (not geom_dag_edges) so the arrowheads clear the node discs instead of hiding under them and near-parallel edges fan apart. r_start / r_end are node radii in data units, tuned to the point size below under coord_equal().
gd    <- dag_tidy$data
nodes <- gd %>% distinct(name, x, y, role, label)
r_start <- 0.45
r_end   <- 0.60
edges <- gd %>%
  filter(!is.na(to)) %>%
  mutate(
    dx = xend - x, dy = yend - y, len = sqrt(dx^2 + dy^2),
    x1 = x    + dx * r_start / len, y1 = y    + dy * r_start / len,
    x2 = xend - dx * r_end   / len, y2 = yend - dy * r_end   / len
  )
### The two Ecoregion -> sampling edges cross the whole figure; give them a deeper bow so they sweep clear of the environment cluster instead of grazing its nodes.
is_long   <- edges$name == "Ecoregion" & edges$to %in% c("Observer", "Season")
edge_arrow <- grid::arrow(length = grid::unit(11, "pt"), type = "closed")

p_dag <- ggplot() +
  geom_curve(data = edges[!is_long, ], aes(x = x1, y = y1, xend = x2, yend = y2),
             curvature = -0.10, colour = "grey40", linewidth = 0.7, lineend = "round",
             arrow = edge_arrow) +
  geom_curve(data = edges[is_long, ], aes(x = x1, y = y1, xend = x2, yend = y2),
             curvature = 0.32, colour = "grey55", linewidth = 0.6, lineend = "round",
             arrow = edge_arrow) +
  geom_point(data = nodes, aes(x, y, colour = role), size = 30) +
  geom_text(data = nodes, aes(x, y, label = label), colour = "white",
            fontface = "bold", size = 3.35, lineheight = 0.9) +
  scale_colour_manual(values = role_cols, name = NULL) +
  coord_equal(clip = "off") +
  expand_limits(x = c(-0.4, 15.0), y = c(1.7, 9.0)) +
  labs(title = "Assumed causal structure: farm management diversification -> bird diversity",
       subtitle = "Left: Ecoregion and the region / environment drivers it summarises. Right: the sampling covariates.\nTwo DAG-valid adjustments: {Climate, LandForest, SpeciesPool, Observer, Season, Year} or {Ecoregion} alone.\nGrey = unobserved. The endemism index is the measured proxy for biogeographic history; see the script for variables held out.") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 12),
        plot.margin = margin(8, 8, 8, 8))

ggsave("Figures/DAG.png", p_dag, width = 19, height = 8, bg = "white")
print(p_dag)
