# Re-analysis checklist — farm management diversification vs bird diversity

When a better dataset lands (more complete diversification indices, additional
farm IDs, etc.), work through this. It records what has already been decided or
tested so nothing is re-litigated from scratch. Full detail is in the synthesis
report `Scripts/qmd/Farm_mgmt_summary.qmd` and `Project_notes.md` (local only,
chronological).

---

## 1. Rebuild the inputs

- [ ] `00_bird_diversity_estimates.R` — re-run at `nboot = 500` **only if the
      point-count data changed**. Both passes; pass 2's low-coverage removal set
      comes from pass 1's own `SC_obs < 0.65`.
- [ ] `01_farm_covariates.R` — re-run **only if the wrangling repo's
      `Site_covs.csv` changed** (climate normals, elevation, coordinates,
      `Distancia_farm`).
- [ ] `01b_species_pool.R` — re-run **only if the farm set changed** (the
      1 788-species Ayerbe range layer is cached in
      `Data/Geospatial/Ayerbe_ranges.gpkg`; elevational limits are layered
      Suárez-Castro 2024 → Ayerbe-guide → Hilty).
- [ ] `02_match_farm_diversity.R` — check the new match rate (was 69 / 94).
      Note which farms gained water/pasture-management values. Drop the
      `sheet_or("Data","Sheet1")` fallback once the renamed xlsx sheets sync.
- [ ] Confirm `Data/Complementary_Biodiversity_Paper_Birds_MJE_June_2026.xls` is
      the current file; check `dist_predio_cercano` and the WVCC / DEM / Biomasa
      columns are still there.

## 2. Model specification — settled, keep as-is

- Bayesian measurement-error regression (`brms`, `resp_se(se, sigma = TRUE)`),
  one index at a time; responses = log Hill numbers q = 0 / 1 / 2.
- **Two environmental-adjustment versions** (`Scripts/dag.R`):
  - `climate` (primary) — `Elev_z + I(Elev_z^2) + prec_z + I(prec_z^2) + canopy_10k_z`
  - `ecoregion` (robustness) — 5-level factor + `canopy_10k_z`
- Sampling controls throughout — `log(Num_pc)_z`, `Num_hab_z`, cyclic
  day-of-year, `(1 | Id_gcs)`, `(1 | CollectorXyear)`.
- **300 m** distance-to-farm cutoff = primary; full set = sensitivity.
- No-index baseline per response for the `bayes_R2` comparison.
- `adapt_delta = 0.99` (0.999 + `iter = 4000` for the small-n Piedemonte fits).
- Priors: `student_t(3, 3, 2.5)` intercept, `normal(0, 0.75)` b,
  `exponential(1)` sd and sigma.

## 3. Specification checks — re-run each with the new data

- [ ] **DAG** (`dag.R`) — re-render; re-run the implied-conditional-independence
      spot-checks; confirm the minimal adjustment set is still
      `{Climate, LandForest, NumPC, Observer, Season, Year}` (NumPC now enters as
      a genuine confounder via the latent `FarmSize -> NumPC` leg, not just
      precision; `LandForest` is an Ecoregion proxy in the climate spec, pure
      precision under the ecoregion factor).
- [ ] **Distance cutoff** — regenerate the `Distancia_farm` histogram; confirm
      300 m still drops only the clearly-displaced survey groups; compare
      300 / 500 / full.
- [ ] **Precipitation functional form** (`04d` + `Precip_vs_elev_sensitivity.csv`)
      — the key confounding finding: the pasture/water coefficient is ~+0.08
      when precipitation is weakly controlled and ~+0.04 (CrI includes 0) when
      it is the dominant environmental control. Re-check; consider making
      `s(prec)` the primary climate term.
- [x] **Drop `(1 | CollectorXyear)`** — `04e_spec_checks.R` (full 5-ecoregion) +
      `04c` (Piedemonte). RE sd 0.18-0.31 (real batch structure) but dropping it
      shifts the focal coefficients by ≤ 0.02 — it is *not* acting like
      `(1 | Ecoregion)` (Cramér's V with Ecoregion only 0.38; the 3 biggest
      datasets span 4-5 regions). Keep it. Re-check with new data.
- [x] **`resp_se` vs point estimates** — `04e_spec_checks.R`. Dropping it gives
      ~25 % *larger* |focal| and ~7 % wider CrIs — `resp_se` is a precision
      weight that trusts the well-sampled assemblages (which show a weaker
      management signal). It is the honest *and* conservative choice; without it
      pasture/water richness CrIs would barely exclude zero. Keep it.
- [ ] **Species pool** (`01b` + `04d_species_pool_test.R`) — `spp_pool` is
      ~89 % explained by elevation + precipitation and adds nothing to the
      management coefficients. Re-check only if the farm set shifts substantially.
- [ ] **Piedemonte-only cut** (`04c_farm_mgmt_piedemonte.R` +
      `05c_...plot.R` + `qmd/Piedemonte_report.qmd`) — region fixed; the
      collinear elev/precip/canopy replaced by `env_pc1` (PC1, ~95 % of joint
      variance); no ecoregion version; sensitivities = full set + drop-collector.
      More farms here would raise the power (n = 17 for water/pasture).
- [ ] **Canopy scale of effect** (`06a_Extract_cc_buff.R` +
      `06b_scale_of_effect.R`) — re-run if the point counts changed; the CSV
      cache means `06a` won't re-extract otherwise. Peak was ~8-9 km,
      topo-adjusted.
- [ ] **Baseline vs index `bayes_R2`** — confirm the indices still add ~0.

## 4. Flagged, not yet done — revisit with better data

- [ ] **`s(precip)` as the primary climate term** so `climate` and `ecoregion`
      converge (the DAG-consistent precipitation adjustment).
- [ ] **`dist_predio_cercano`** (distance to nearest neighbouring property, from
      the MJE xls; range 0-36 km, median ~1.4 km) as a covariate — a farm
      *isolation* measure, **uncorrelated with `Distancia_farm`** (r ≈ 0;
      `Distancia_farm` is instead a data-quality flag for off-farm surveys).
      Test as a predictor and/or a spatial-autocorrelation control.
- [ ] **Index-level Mundlak / within-between decomposition** — ecoregion mean +
      farm deviation per index; the within coefficient is the confound-free
      estimate. `06b` already has the canopy analogue (`within_eco`).
- [ ] **Functional / phylogenetic diversity** responses, or
      species-of-conservation-concern subsets — may respond where taxonomic Hill
      numbers do not.
- [ ] **eBird composition** (`auk`) or an **endemism-weighted pool** — for the
      *compositional* biogeographic residual the range-map richness count misses.
      Needs the Colombia EBD download; `auk` not yet installed.
- [ ] **Farm area (hectares)** as a covariate — the DAG's latent `FarmSize` has
      an unblockable `FarmSize -> BirdDiv` path (species-area / interior habitat /
      edge) on top of the `-> FarmDiv` and `-> NumPC` legs that conditioning on
      `NumPC` already closes. Measuring hectares would close the last leg and
      remove a positive bias on the management coefficient. Check whether the MJE
      xls or the wrangling repo carries a farm-area column.
- [ ] **On-farm woody cover / practice detail the four indices miss** — the
      `FarmerValues -> {snag retention, pesticide use, fencerow width, ...} ->
      BirdDiv` path is a *farm-scale* forest/practice effect that `canopy_10k`
      (10 km buffer, set by topography and hundreds of landholders — a scale
      mismatch with farmer decisions) cannot stand in for. A farm-boundary canopy
      or hedgerow-length layer, or a finer practice inventory, would let this be
      adjusted rather than left as a caveat.
- [ ] **Prior sensitivity on `sd` / `sigma`** given many single-assemblage farms.
- [ ] **Investigate the missing** water/pasture-management values (9 farms,
      mostly Piedemonte) — may be resolved by the new data.
- [ ] **Retire `Data/Excels/*_06.04.26.csv`** (unread since the in-repo regen).
- [ ] **Update `Ch1-ssp-birds/Project_notes.md`** (stale re: the diversity
      sub-pipeline move).

## 5. Interpretation anchors — test the new data against these

1. **No detectable effect of management diversification on bird taxonomic
   diversity.** Land use and All practices (the region-separable indices) are
   flat null under every spec.
2. **The pasture / water positive under the parametric `climate` spec is
   precipitation confounding** — `Pasture_mgmt_div ~ precip` r = 0.30, and it
   vanishes whenever precipitation is adjusted well (including under the
   `Ecoregion` factor, R² 0.83 with precip).
3. **The `Ecoregion` adjustment is legitimate confounding control, not
   over-adjustment** — the attenuation is reproducible with no region factor at
   all, purely by controlling precipitation better. (Caveat: with everything this
   collinear this is an inference from the pattern, not a proof.)
4. **Surrounding landscape forest cover at ~8-9 km is positively associated with
   bird richness**, net of climate — a landscape-context effect, not on-farm
   management.
5. Everything is **limited by the collinearity** of elevation, precipitation,
   canopy, ecoregion and the region-bound indices — identification of a
   management effect from observational data is fragile by construction.

## 6. Branch

All of the above lives on branch **`refactor-pipeline-scripts`** — not yet
merged or pushed. Review and merge before the re-analysis, or rebase the new
work onto it.
