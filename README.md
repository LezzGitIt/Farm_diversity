# Does farm-management diversification benefit bird diversity?

## Overview

This repository contains the analysis code linking farm-level **management diversification** on Colombian Sustainable Cattle Ranching (SCR) farms to **bird taxonomic diversity**.

Bird diversity is summarised as iNEXT Hill numbers (richness q = 0, Shannon q = 1, Simpson q = 2) estimated per `[collector . farm . year-group . season]` assemblage from standardised point counts. Farm management is summarised by four `[0-1]` diversification indices — land use, water management, pasture management, and an overall "all practices" composite — compiled by Maria Jimena Esquivel. The two datasets are joined on farm ID (69 farms have both).

The full writeup is `Scripts/qmd/Farm_mgmt_summary.qmd` (renders to PDF): question, confounding structure, approach, preliminary results, interpretation.

## Repository structure

All code lives in `Scripts/`.

**Shared**

- `Model_fns.R` — Shared helpers: input paths (the data-paper folder, and which dated diversity-estimate file the downstream scripts read) + the brms model machinery
- `dag.R` — Assumed causal DAG behind the models → `Figures/DAG.png`

**Data prep and exploration**

- `00_bird_diversity_estimates.R` — Point counts → per-assemblage Hill-number diversity, two estimates: point-count-standardised incidence, and abundance/coverage-standardised
- `01a_farm_data.R` — Farm-level environmental covariates + the four diversification indices, matched to farms with at least one assemblage-level bird-diversity estimate
- `01b_species_pool.R` — Per-farm potential species pool from range maps (raw count + endemism-weighted)
- `02a_extract_canopy_buffers.R` — Woody-vegetation canopy cover (WVCC) around each assemblage at 18 radii (200 m – 10 km)
- `02b_scale_of_effect.R` — Which WVCC radius best explains diversity 
- `03_exploratory.R` — Exploratory plots: the indices, and how everything varies among ecoregions

**Management-diversification analysis**

- `04a_farm_mgmt_models.R` — PRIMARY analysis: bird diversity ~ each diversification index, two responses × two adjustment sets, Bayesian measurement-error models (brms)
- `04b_farm_mgmt_robustness.R` — Model specification checks (set `run_sections` to run a subset): single-ecoregion cut, species-pool terms, measurement-error / likelihood checks, DAG-ideal adjustment + collinearity, spline vs quadratic environment
- `05_farm_mgmt_plots.R` — Every figure for the management question (reads 04a/04b output)

**Reports (`qmd/`)**

- `Farm_mgmt_summary.qmd` — Synthesis report for collaborators: question, approach, results
- `_archive/` — Superseded write-ups, kept for methods history only (not part of the run order): `02_Analysis_iNEXT.qmd` = the Cmax / iNEXT-4steps rationale behind `00`; `03_Farm_diversity.qmd` = the old negative-binomial GLMM that motivated the analysis

`00`–`02b` are data prep plus the landscape-forest scale-of-effect analysis; `03` is exploratory work to understand the data; `04`–`05` are the management-diversification analysis. Generated tables, models and figures are written to `Derived/` and `Figures/` (plus the cached canopy extraction in `Data/Geospatial/`); none are tracked.

## Running the analysis

Scripts use paths relative to the project root and are numbered in run order. The bird point-count tables are read from the sibling data-paper repo (`../Ssp-bird-data-wrangling/Derived/Excels/`, English column names; set `WRANGLING_EXCELS` in `Scripts/Model_fns.R`). Every script downstream of `00` reads the diversity export whose date stamp is pinned by `TAX_DIV_VERSION` in `Model_fns.R` -- after re-running `00`, update that value deliberately (re-running `00` alone does not change what the models see).

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

The bird point-count data come from the data paper *"Bird diversity in working landscapes of Colombia"* (Skinner et al.), currently in prep; they are not yet public. The farm-management diversification indices were provided by Maria Jimena Esquivel. Neither dataset is distributed in this repository at this stage.

## Dependencies

R packages: `tidyverse`, `brms` (with `rstan`), `iNEXT`, `iNEXT.4steps`, `sf`, `terra`, `tidyterra`, `lme4`, `MuMIn`, `mgcv`, `loo`, `posterior`, `cowplot`, `GGally`, `dagitty`, `ggdag`, `readxl`, `janitor`, `knitr`.

## Citation

\[To be added upon publication\]
