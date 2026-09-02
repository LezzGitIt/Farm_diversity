# Does farm-management diversification benefit bird diversity?

## Overview

This repository contains the analysis code linking farm-level **management
diversification** on Colombian Sustainable Cattle Ranching (SCR) farms to
**bird taxonomic diversity**.

Bird diversity is summarised as iNEXT Hill numbers (richness q = 0, Shannon
q = 1, Simpson q = 2) estimated per `[collector . farm . year-group . season]`
assemblage from standardised point counts. Farm management is summarised by four
`[0-1]` diversification indices — land use, water management, pasture management,
and an overall "all practices" composite — compiled by Maria Jimena Esquivel. The
two datasets are joined on farm ID (69 farms have both).

The full writeup is `Scripts/qmd/Farm_mgmt_summary.qmd` (renders to PDF):
question, confounding structure, approach, preliminary results, interpretation.

## Repository structure

```
Scripts/
  Model_fns.R                    # Shared helpers (file picker + the brms model machinery)
  dag.R                          # Assumed causal DAG behind the models -> Figures/DAG.png

  00_bird_diversity_estimates.R  # Point counts -> per-assemblage Hill-number diversity (two estimates:
                                 #   point-count-standardised incidence, and abundance/coverage-standardised)
  01a_farm_data.R                # Farm-level environmental covariates + the four diversification indices,
                                 #   matched to farms that have a bird-diversity estimate
  01b_species_pool.R             # Per-farm potential species pool from range maps (raw count + endemism-weighted)
  02a_extract_canopy_buffers.R   # Woody-vegetation canopy cover around each assemblage at 18 radii (200 m - 10 km)
  02b_scale_of_effect.R          # Which radius best explains diversity -> picks the 10 km landscape-forest covariate
  03_exploratory.R               # Exploratory plots: the indices, and how everything varies among ecoregions

  04a_farm_mgmt_models.R         # PRIMARY analysis: bird diversity ~ each diversification index,
                                 #   two responses x two adjustment sets, Bayesian measurement-error models (brms)
  04b_farm_mgmt_robustness.R     # Specification checks (toggleable sections): single-ecoregion cut,
                                 #   species-pool terms, measurement-error / likelihood checks,
                                 #   DAG-ideal adjustment + collinearity, spline vs quadratic environment
  05_farm_mgmt_plots.R           # Every figure for the management question (reads 04a/04b output; never refits)

  qmd/
    Farm_mgmt_summary.qmd         # The whole-story synthesis: question, approach, results
```

`00`–`02b` are data prep plus the landscape-forest scale-of-effect analysis;
`03` is exploratory work to understand the data;
`04`–`05` are the management-diversification analysis. Generated data, models and
figures are written to `Data/`, `Derived/` and `Figures/` and are not tracked.

## Running the analysis

Scripts use paths relative to the project root and are numbered in run order:

``` r
source("Scripts/00_bird_diversity_estimates.R")   # slow (~1 h); re-run only when the point-count data changes
source("Scripts/01a_farm_data.R")
source("Scripts/01b_species_pool.R")
source("Scripts/02a_extract_canopy_buffers.R")    # heavy (reads national rasters); outputs are cached
source("Scripts/02b_scale_of_effect.R")
source("Scripts/03_exploratory.R")
source("Scripts/04a_farm_mgmt_models.R")          # ~100 brms fits; cached, so re-runs are fast
source("Scripts/04b_farm_mgmt_robustness.R")      # the sensitivity checks
source("Scripts/05_farm_mgmt_plots.R")
```

Then render the report from the project root:

``` bash
quarto render Scripts/qmd/Farm_mgmt_summary.qmd
```

## Data availability

The bird point-count data come from the data paper *"Bird diversity in working
landscapes of Colombia"* (Skinner et al.), currently in prep; they are not yet
public. The farm-management diversification indices were provided by Maria Jimena
Esquivel. Neither dataset is distributed in this repository at this stage.

## Dependencies

R packages: `tidyverse`, `brms` (with `rstan`), `iNEXT`, `iNEXT.4steps`, `sf`,
`terra`, `tidyterra`, `lme4`, `MuMIn`, `mgcv`, `loo`, `posterior`, `cowplot`,
`GGally`, `dagitty`, `ggdag`, `readxl`, `janitor`, `knitr`.

## Citation

[To be added upon publication]
