# Causal DAG for the bird-diversity ~ farm-management-diversification analysis ----

### Reference artifact, not a pipeline stage (like Scripts/Farm_diversity_fns.R). It encodes the assumed causal structure behind Scripts/04_farm_mgmt_mod.R, derives the adjustment sets the model specs should target, and lists the testable implied independencies. Consult it before trusting a coefficient from 04. Produces Figures/DAG.png and Derived/Excels/DAG_adjustment_sets.csv.

# Setup ----

library(tidyverse)
library(dagitty)
library(ggdag)
library(cowplot)

ggplot2::theme_set(theme_dag())
dir.create("Figures", showWarnings = FALSE)
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

# The DAG ----

### Node meanings
# Mgmt        one of the four [0-1] diversification indices (exposure, one model per index).
# BirdDiv     iNEXT Hill-number diversity of the assemblage (outcome, as measured -- so sampling nuisances point into it).
# Ecoregion   5-level biogeographic region; the SCR programme also rolled out differently by region, so it causes Mgmt too.
# Elev        farm-mean elevation (temperature is r = -0.99 with it, so Temp is the same node -- use one).
# Precip      farm-mean annual precipitation.
# Canopy10k   woody-vegetation cover in the surrounding 10 km (the scale-of-effect peak from 06b).
# NumHab      count of habitat types recorded across the farm's point counts.
# NumPC       number of point counts in the assemblage (sampling effort).
# DOY         mean day of year of the surveys (Nearctic-migrant season proxy).
# Collector   [data collector x year] batch = survey protocol (radius, duration, single vs repeat visits).
# FarmSize    unobserved farm area / capacity; larger operations get more point counts and may diversify management more. Deliberately NOT drawn into NumHab: habitat *count* on a working cattle farm is a land-use decision and a regional-template feature, not a function of raw hectares (a large monoculture pasture still scores NumHab = 1).

dag <- dagitty('dag {
  Mgmt [exposure]
  BirdDiv [outcome]
  FarmSize [latent]

  Ecoregion -> Elev
  Ecoregion -> Precip
  Ecoregion -> Canopy10k
  Ecoregion -> Mgmt
  Ecoregion -> BirdDiv
  Ecoregion -> Collector

  Elev -> Canopy10k
  Elev -> Mgmt
  Elev -> BirdDiv

  Precip -> Canopy10k
  Precip -> Mgmt
  Precip -> BirdDiv

  Canopy10k -> BirdDiv

  Mgmt -> NumHab
  Mgmt -> BirdDiv
  NumHab -> BirdDiv

  FarmSize -> Mgmt
  FarmSize -> NumPC

  Ecoregion -> DOY
  Collector -> DOY
  Collector -> Canopy10k

  NumPC -> BirdDiv
  DOY -> BirdDiv
  Collector -> BirdDiv
}')

### Edges added 2026-08-28 after data spot-checks. Collector = [dataset x year] here also stands in for WHICH farm subset that dataset covered -- so it drives the survey window (Collector -> DOY) and the surrounding landscape of those farms (Collector -> Canopy10k). This clears the strong flags (DOY ~ Ecoregion/Elev p < 1e-4; Canopy10k ~ DOY | Ecoregion p = 2e-7). None of these touch identification of the Mgmt -> BirdDiv effect (Collector, DOY, Canopy10k have no arrow into Mgmt), only whether they are clean precision controls -- and they are all in the model anyway. Two residual weak flags left as-is (p ~ 0.02-0.03, 6 tests): Canopy10k ~ Pasture_mgmt_div and Collector ~ Land_use_div, both given region -- possibly which-farms-each-collector-visited confounding; revisit if they strengthen.

coordinates(dag) <- list(
  x = c(Mgmt = 0, BirdDiv = 5, Ecoregion = 2.4, Elev = 1.5, Precip = 2.6,
        Canopy10k = 4.0, NumHab = 2.5, FarmSize = 0, NumPC = 1, DOY = 2.6, Collector = 3.9),
  y = c(Mgmt = 1.6, BirdDiv = 1.6, Ecoregion = 4.2, Elev = 3.2, Precip = 2.7,
        Canopy10k = 3.2, NumHab = 0.7, FarmSize = 3.2, NumPC = -0.6, DOY = -0.6, Collector = -0.6)
)

# Adjustment sets ----

### The total effect of Mgmt on BirdDiv. NumHab is a mediator (Mgmt -> NumHab -> BirdDiv), so it must stay OUT of the total-effect model -- that is the 04 "ecoregion_numhab" spec, which by conditioning on the mediator targets the DIRECT effect instead.
total_sets <- adjustmentSets(dag, exposure = "Mgmt", outcome = "BirdDiv",
                             effect = "total", type = "canonical")

### The direct effect (blocking the NumHab path deliberately).
direct_sets <- adjustmentSets(dag, exposure = "Mgmt", outcome = "BirdDiv",
                              effect = "direct", type = "minimal")

### Minimal sufficient sets for the total effect (usually the practical choice).
total_min <- adjustmentSets(dag, exposure = "Mgmt", outcome = "BirdDiv",
                            effect = "total", type = "minimal")

cat("== Total-effect adjustment (canonical) ==\n"); print(total_sets)
cat("\n== Total-effect adjustment (minimal) ==\n"); print(total_min)
cat("\n== Direct-effect adjustment (minimal) ==\n"); print(direct_sets)

### READING (2026-08-28)
# The minimal total-effect set is { Ecoregion, Elev, Precip, NumPC }. Elev and Precip carry arrows into BOTH Mgmt and BirdDiv that Ecoregion does not mediate (Elev's and Precip's only parent is Ecoregion, but as forks on the Mgmt <- Elev -> BirdDiv path they stay open unless conditioned directly). Symmetrically, Ecoregion -> BirdDiv (biogeographic species pool) is not mediated by climate. So:
#  - the 04 "ecoregion" spec (Ecoregion only) leaves the Elev/Precip forks open  -> mild UNDER-adjustment on climate, but in practice Elev/Precip are ~0.8 R^2 with Ecoregion so little is left;
#  - the 04 "climate" spec (Elev + Precip, no Ecoregion) leaves Ecoregion -> BirdDiv open -> residual biogeographic confounding.
# The fully sufficient model would include the 5-level factor AND both continuous axes, but they are ~80% collinear -> not identifiable here. So the two specs are BOUNDING CASES that happen to agree (both ~0 under ecoregion; weak positive water/pasture under climate we attribute to the residual Ecoregion path). Adding poly(Elev,2)+poly(Precip,2) to the climate spec makes it a closer Ecoregion substitute and is the agreed next step; also worth trying the combined Ecoregion + poly(Elev,2) + poly(Precip,2) spec and checking posterior-SD inflation before treating it as primary.
# Canopy10k is in the canonical set but NOT the minimal one -> it is a precision covariate (predicts BirdDiv, no arrow to Mgmt), fine to keep in the climate spec, optional in the ecoregion spec.
# FarmSize (latent) confounds via FarmSize -> Mgmt and FarmSize -> NumPC -> BirdDiv; conditioning on NumPC closes that path. Residual worry only if FarmSize -> BirdDiv directly (large farms = more interior habitat) -- a stated limitation; measuring farm area would settle it.

### FarmSize is latent, so any set that requires it is not feasible -- flag which minimal sets are actually estimable.
feasible <- Filter(function(s) !("FarmSize" %in% s), as.list(total_min))
cat("\n== Feasible minimal sets for the total effect (no latent FarmSize) ==\n")
if (length(feasible) == 0) {
  cat("  none listed -- see notes; the working set is {Ecoregion, NumPC} (or {Elev, Precip, NumPC})\n")
} else {
  print(feasible)
}

# Testable implications ----

### Implied conditional independencies -- check the strong ones against the data as a DAG sanity test.
ici <- impliedConditionalIndependencies(dag)
cat("\n== Implied conditional independencies (first 15) ==\n")
print(head(ici, 15))

### Spot-check a few against the 04 modelling frame (if it exists). A small p-value = the data contradict that implied independence = the DAG is missing an edge.
data_path <- "Derived/Excels/Farm_mgmt_model_data.csv"
if (file.exists(data_path)) {
  d <- read_csv(data_path, show_col_types = FALSE) %>%
    distinct(Id_gcs, .keep_all = TRUE)   # farm-level covariates: one row per farm

  p_of <- function(fit, term) tryCatch(summary(fit)$coefficients[term, "Pr(>|t|)"], error = function(e) NA_real_)

  checks <- tibble(
    implication = c(
      "Canopy10k _||_ NumPC",
      "DOY _||_ Elev | Ecoregion",
      "Canopy10k _||_ DOY | Ecoregion, Collector",
      "Canopy10k _||_ Land_use_div | Ecoregion, Elev, Precip",
      "Canopy10k _||_ Pasture_mgmt_div | Ecoregion, Elev, Precip",
      "Collector _||_ Land_use_div | Ecoregion"
    ),
    p_value = c(
      cor.test(d$canopy_10k, d$Num_pc_log)$p.value,
      p_of(lm(Elev_mean ~ doy + Ecoregion, d), "doy"),
      p_of(lm(canopy_10k ~ doy + Ecoregion + CollectorXyear, d), "doy"),
      p_of(lm(canopy_10k ~ Land_use_div + Ecoregion + Elev_mean + Tot_prec_mean, d), "Land_use_div"),
      p_of(lm(canopy_10k ~ Pasture_mgmt_div + Ecoregion + Elev_mean + Tot_prec_mean, d), "Pasture_mgmt_div"),
      anova(lm(Land_use_div ~ Ecoregion + CollectorXyear, d), lm(Land_use_div ~ Ecoregion, d))$`Pr(>F)`[2]
    )
  ) %>% mutate(flag = if_else(p_value < 0.05, "** tension with the DAG (6 tests -- weigh multiple-testing)", ""))
  cat("\n== Data spot-checks of implied independencies ==\n")
  print(checks, width = Inf)
  cat("Any remaining flags are among precision covariates (Canopy10k / DOY / Collector); none has an arrow into Mgmt, so they do not affect identification of the management effect, only how clean those controls are.\n")
}

# How the 04 specs map onto the DAG ----

spec_map <- tribble(
  ~spec,              ~adjusts,                                                      ~targets,        ~dag_verdict,
  "ecoregion",        "Ecoregion + NumPC + DOY + Collector",                          "total (bound)", "closes the Ecoregion backdoor; Elev/Precip forks left open but they are ~0.8 R2 with Ecoregion so little residual -- upper bound on adjustment (risks absorbing region-tracking management signal)",
  "climate",          "poly(Elev,2) + poly(Precip,2) + Canopy10k + NumPC + DOY + Collector", "total (bound)", "mechanistic Ecoregion substitute; quadratics + canopy tighten the stand-in; Ecoregion->BirdDiv (species pool) left open -> lower bound (residual biogeographic confounding)",
  "ecoregion + climate", "Ecoregion + poly(Elev,2) + poly(Precip,2) + sampling",      "total (target)","the DAG-sufficient set; try it, but check posterior-SD inflation from Ecoregion/Elev/Precip collinearity before making it primary",
  "ecoregion_numhab", "ecoregion spec + NumHab",                                      "direct",        "conditions on the Mgmt->NumHab->BirdDiv mediator on purpose -- interpret as direct, not total"
)
cat("\n== 04 spec <-> DAG ==\n"); print(spec_map, width = Inf)

adj_out <- bind_rows(
  tibble(effect = "total",  set = map_chr(as.list(total_min),  ~ paste(sort(.x), collapse = " + "))),
  tibble(effect = "direct", set = map_chr(as.list(direct_sets), ~ paste(sort(.x), collapse = " + ")))
)
write_csv(adj_out, "Derived/Excels/DAG_adjustment_sets.csv")

# Figure ----

dag_tidy <- dag %>%
  tidy_dagitty() %>%
  mutate(role = case_when(
    name == "Mgmt" ~ "exposure",
    name == "BirdDiv" ~ "outcome",
    name == "NumHab" ~ "mediator",
    name == "FarmSize" ~ "unobserved",
    name %in% c("NumPC", "DOY", "Collector") ~ "sampling nuisance",
    TRUE ~ "confounder / covariate"
  ))

role_cols <- c(exposure = "#1b7837", outcome = "#762a83", mediator = "#d95f02",
               unobserved = "grey70", "sampling nuisance" = "#4393c3",
               "confounder / covariate" = "grey30")

p_dag <- ggplot(dag_tidy, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_colour = "grey55") +
  geom_dag_point(aes(colour = role), size = 20) +
  geom_dag_text(colour = "white", size = 2.5) +
  scale_colour_manual(values = role_cols, name = NULL) +
  expand_limits(x = c(-0.4, 5.4)) +
  labs(title = "Assumed causal structure: farm management diversification -> bird diversity",
       subtitle = "NumHab is a mediator (out of the total-effect model); Canopy10k is a precision covariate; sampling nuisances adjusted for precision") +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 9))

ggsave("Figures/DAG.png", p_dag, width = 12, height = 8.5, bg = "white")
print(p_dag)
