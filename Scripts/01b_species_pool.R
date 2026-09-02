# Regional species-pool covariate from Ayerbe range maps + Suarez-Castro elevational limits ----

### Per-farm "how many bird species could occur here" -- a potential-species-pool covariate. Range maps and published elevational limits are external to the point-count data, so this is NOT circular: it predicts the regional avifauna available to a farm without any knowledge of what was surveyed there. It is the one thing `Ecoregion` carries that elevation + climate + canopy in Scripts/04a_farm_mgmt_models.R do not (biogeographic species pool); see Scripts/dag.R (the SpeciesPool branch).

### Two-step filter per farm:
###   1. GEOSPATIAL -- the species' Ayerbe Colombian range polygon contains the farm (point-in-polygon).
###   2. ELEVATIONAL -- the farm elevation falls within the species elevational range (no buffer). Elevational limits are layered by source priority: Suarez-Castro et al. (2024) AOH table S3 (1,652 spp), then the manual book digitisations from ../Ssp-bird-data-wrangling/Scripts/03_FT_elev.R -- Hazen's read-out of the Ayerbe-Quinones (2018) field guide, then Hazen's read-out of Hilty. Species still without a limit (~13%, mostly Nearctic migrants / seabirds) pass this step unfiltered.

### Also builds RANGE-RARITY-WEIGHTED pool metrics (Aaron's idea) from AVONET global range sizes: `pool_we` (weighted endemism, Sum 1/range), `pool_wes` (softer, Sum 1/sqrt range), `pool_cwe` (corrected -- proportion range-restricted), `pool_rr` (count of restricted-range species, < 50,000 km2 BirdLife EBA criterion). The hypothesis: total pool richness tracks contemporary climate (why `pool_point` is redundant), but the range-restricted component tracks biogeographic history, which `Ecoregion` carries and climate does not -- so a weighted metric may be the non-redundant compositional axis. `Data/Farm_species_pool.csv` carries all of them.

### Outputs -> Data/Farm_species_pool.csv, treated as raw input like Data/Farm_covariates.csv. The combined Ayerbe range layer is cached in Data/Geospatial/Ayerbe_ranges.gpkg (the 1,890 shapefile reads are slow on a cold OneDrive).

### OUTCOME (2026-08-29): `pool_point ~ poly(Elev,2) + poly(precip,2)` R^2 = 0.87 (unique: elev 0.29, precip 0.45, shared 0.12). Scripts/04b_farm_mgmt_robustness.R (pool_blocks section) tested it properly (blocks that never combine all three axes) and it does not change the management coefficients -> NOT added to 04a, kept for the record. Along the way the pool_blocks section + Derived/Excels/Precip_vs_elev_sensitivity.csv showed the pasture/water "signal" in the climate spec is residual PRECIPITATION confounding (Pasture_mgmt_div ~ precip r = 0.30), not a missing species-pool axis; see Project_notes.md.
### OUTCOME (weighted-endemism metrics, 2026-08-31): the hypothesis holds. `pool_we` / `pool_cwe` are much LESS climate-redundant than the raw count: R^2(~ poly(elev,2)+poly(precip,2)) drops from 0.87 (pool_point) to 0.43-0.45, and the within-Ecoregion fraction rises from 0.22 to 0.72. They are negatively correlated with `pool_point` (r ~ -0.5) -- a genuinely different axis: the Andean coffee region (Cafetera) and Boyaca-Santander are the endemism hotspots, the species-rich lowlands (Piedemonte, Bajo Magdalena) are endemism-poor. `pool_rr` (raw count < 50k km2) stays climate-bound (R^2 0.85, r_elev 0.90 -- restricted-range birds are overwhelmingly Andean, so counting them tracks elevation); `pool_wes` (softer weight) is intermediate (0.66). So the full 1/range weighting (`pool_we`, or `pool_cwe`) is the one worth testing in 04d / 04h as the compositional biogeographic axis the raw count misses. Caveats: still r ~ 0.6 with elevation and ~0.48 R^2 with Ecoregion (not clean); the metric has a huge dynamic range driven by a handful of ultra-narrow-range species, so check robustness (drop-one-species).

# Setup ----
library(tidyverse)
library(sf)
library(readxl)

sf::sf_use_s2(TRUE)
dir.create("Data/Geospatial", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

ayerbe_dir   <- "../Geospatial_data/Ayerbe_shapefiles_1890spp"
## LOCAL external-dataset path -- not in this repo. Repoint to your copy of the
## Ayerbe (2018) per-species elevational-range tables.
elev_dir     <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Elev_ranges"
ranges_cache <- "Data/Geospatial/Ayerbe_ranges.gpkg"

# Combined Ayerbe range layer (cached) ----

if (file.exists(ranges_cache)) {
  message("loading cached range layer: ", ranges_cache)
  ranges <- st_read(ranges_cache, quiet = TRUE)
} else {
  shp <- list.files(ayerbe_dir, pattern = "[.]shp$", full.names = TRUE)
  message(length(shp), " Ayerbe shapefiles -- reading (slow on a cold OneDrive) ...")
  ranges <- map(seq_along(shp), function(i) {
    if (i %% 200 == 0) message("  ", i, " / ", length(shp))
    g <- tryCatch(st_read(shp[i], quiet = TRUE), error = function(e) NULL)
    if (is.null(g) || nrow(g) == 0) return(NULL)
    g <- st_zm(g)
    tibble(species = gsub("[.]shp$", "", basename(shp[i])), geometry = st_geometry(g)) %>%
      st_as_sf(crs = st_crs(g))
  }) %>%
    compact() %>%
    map(~ st_transform(.x, 4326)) %>%
    bind_rows() %>%
    st_make_valid()
  ranges <- ranges[st_is(ranges, c("POLYGON", "MULTIPOLYGON")), ] %>%
    st_cast("MULTIPOLYGON", warn = FALSE)
  st_write(ranges, ranges_cache, delete_dsn = TRUE, quiet = TRUE)
  message("wrote ", ranges_cache, ": ", n_distinct(ranges$species), " species, ", nrow(ranges), " features")
}

# Elevational limits: layered by source priority ----

suarez <- read_csv(file.path(elev_dir, "Suarez_castro_AOH_birds_table_S3_V3.csv"), show_col_types = FALSE)

## Hazen's manual read-outs from the two field guides (03_FT_elev.R), keyed by Ayerbe binomial -- the range-polygon species names
hazen_ayerbe <- read_xlsx(file.path(elev_dir, "Hazen_Elev_ranges_Ayerbe.xlsx")) %>%
  transmute(name = Species_ayerbe, elev_lo = Min_ayerbe, elev_hi = Max_ayerbe)
hazen_hilty <- suppressWarnings(read_xlsx(file.path(elev_dir, "Hazen_Elev_ranges_Hilty.xlsx"))) %>%
  transmute(name = Species_ayerbe, elev_lo = Min_Hilty, elev_hi = Max_Hilty)

## bind in priority order (Suarez-Castro first, on either its BirdLife or Clements name), keep the first non-NA per name
elev_by_name <- bind_rows(
  suarez %>% transmute(name = Scientific.Name,       elev_lo = Minimum.elevation, elev_hi = Maximum.elevation, src = "Suarez-Castro 2024"),
  suarez %>% transmute(name = Name.Clements.eBird.,  elev_lo = Minimum.elevation, elev_hi = Maximum.elevation, src = "Suarez-Castro 2024"),
  hazen_ayerbe %>% mutate(src = "Ayerbe 2018 guide"),
  hazen_hilty  %>% mutate(src = "Hilty guide")
) %>%
  filter(!is.na(name), !is.na(elev_lo), !is.na(elev_hi)) %>%
  distinct(name, .keep_all = TRUE)

elev_lookup <- tibble(species = sort(unique(ranges$species))) %>%
  left_join(elev_by_name, by = c("species" = "name"))
message(sum(!is.na(elev_lookup$elev_lo)), " of ", nrow(elev_lookup),
        " Ayerbe species have an elevational limit (",
        round(100 * mean(!is.na(elev_lookup$elev_lo))), "%); the rest pass the filter unfiltered. Sources: ",
        paste(sprintf("%s %d", names(table(elev_lookup$src)), table(elev_lookup$src)), collapse = ", "))

# Global range size per species (AVONET) -- for the range-rarity-weighted pool ----

### `pool_point` counts every species equally and is ~87% a composite of elevation + precipitation (it tracks contemporary climate / productivity, as total richness does). The RANGE-RESTRICTED component of the pool is instead set by biogeographic history (isolation, orogeny, refugia), which `Ecoregion` carries and climate does not. So a range-rarity-weighted pool should be the compositional axis that is NOT redundant with the climate terms.
### `Range.Size` (global extent-of-occurrence area, km^2) from AVONET1 (Tobias et al. 2022, BirdLife taxonomy). Direct name match to the Ayerbe range-polygon species covers ~95%; a genus-median fills most of the rest; the last ~2% are excluded from the weighted metrics (they still count in pool_point).

## LOCAL external-dataset path -- not in this repo. Repoint to your copy of
## AVONET1_BirdLife.csv (Tobias et al. 2022).
avonet_path <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Avonet_Data/TraitData/AVONET1_BirdLife.csv"
avonet <- read_csv(avonet_path, show_col_types = FALSE) %>%
  transmute(name = Species1, genus = word(Species1, 1), range_km2 = as.numeric(Range.Size))

avo_by_name  <- set_names(avonet$range_km2, avonet$name)
avo_by_genus <- avonet %>% filter(!is.na(range_km2)) %>%
  summarize(gmed = median(range_km2), .by = genus) %>% { set_names(.$gmed, .$genus) }

## keyed by the RAW range-polygon species names (as elo/ehi are); match AVONET on a space-normalised copy, genus-median fallback
range_species <- sort(unique(ranges$species))
range_norm    <- str_squish(str_replace_all(range_species, "_", " "))
range_km2 <- coalesce(avo_by_name[range_norm], avo_by_genus[word(range_norm, 1)]) %>%
  set_names(range_species)

## BirdLife "restricted-range species" threshold (Stattersfield et al. 1998, Endemic Bird Areas): breeding range < 50,000 km^2
rr_threshold <- 50000

message(sum(!is.na(range_km2)), " of ", length(range_km2), " range-polygon species have a global range size (",
        round(100 * mean(!is.na(range_km2))), "%); ",
        sum(range_km2 < rr_threshold, na.rm = TRUE), " are restricted-range (< ", rr_threshold, " km2).")

# Farms ----

fc <- read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))
farms <- st_as_sf(fc, coords = c("Long_mean", "Lat_mean"), crs = 4326)

elo <- set_names(elev_lookup$elev_lo, elev_lookup$species)
ehi <- set_names(elev_lookup$elev_hi, elev_lookup$species)

# Pool per farm (point-in-polygon) ----

geo_hits <- st_intersects(farms, ranges)   # step 1: which range polygons contain each farm point

Farm_species_pool <- map_dfr(seq_len(nrow(farms)), function(i) {
  sp_geo <- unique(ranges$species[geo_hits[[i]]])
  fe <- fc$Elev_mean[i]
  lo <- elo[sp_geo]; hi <- ehi[sp_geo]
  in_band <- is.na(lo) | (hi >= fe & lo <= fe)   # step 2: farm elevation within [species_min, species_max]; NA elev -> keep
  sp_pool <- sp_geo[in_band]                      # the species pool for this farm
  rk <- range_km2[sp_pool]; rk <- rk[!is.na(rk)]  # their global range sizes (drop the ~2% without one)
  tibble(
    Id_gcs = fc$Id_gcs[i],
    pool_geo_point = length(sp_geo),
    pool_point     = length(sp_pool),
    ## weighted endemism (Williams et al. 1996; Crisp et al. 2001): sum of inverse range sizes, in units of 1 / million km^2
    pool_we  = sum(1e6 / rk),
    ## a softer version -- sum of 1 / sqrt(range), less dominated by the very narrowest ranges
    pool_wes = sum(1000 / sqrt(rk)),
    ## corrected weighted endemism -- the PROPORTION of the pool that is range-restricted (decoupled from raw richness)
    pool_cwe = sum(1e6 / rk) / length(rk),
    ## count of restricted-range species (BirdLife EBA criterion, < 50,000 km^2)
    pool_rr  = sum(rk < rr_threshold),
    n_range_known = length(rk)
  )
})

write_csv(Farm_species_pool, "Data/Farm_species_pool.csv")

# Diagnostics ----

pool_eco <- Farm_species_pool %>%
  left_join(fc %>% select(Id_gcs, Ecoregion, Elev_mean, Tot_prec_mean), by = "Id_gcs")

cat("\n== pool_point per farm (mean [range]) ==\n")
pool_eco %>%
  summarize(pool_point = sprintf("%.0f [%.0f-%.0f]", mean(pool_point), min(pool_point), max(pool_point)),
            pre_elev_filter = sprintf("%.0f", mean(pool_geo_point)), .by = Ecoregion) %>%
  arrange(Ecoregion) %>% print(width = Inf)

cat("\nelevation filter keeps", round(100 * mean(pool_eco$pool_point / pool_eco$pool_geo_point)), "% of the geospatial hits on average\n")

# How climate-redundant is each pool metric? ----

### The go / no-go test: `pool_point` is ~87% a climate composite and adds nothing to Scripts/04a (tested in Scripts/04d, to be folded into 04b). A range-rarity-weighted metric earns a place only if it carries variation that elevation + precipitation do NOT -- i.e. a LOWER `r2_climate_poly` and a HIGHER within-Ecoregion fraction than `pool_point`.

pool_metrics <- c("pool_point", "pool_we", "pool_wes", "pool_cwe", "pool_rr")

pool_diag <- map_dfr(pool_metrics, function(v) {
  y <- pool_eco[[v]]
  tibble(
    metric          = v,
    r_precip        = cor(y, pool_eco$Tot_prec_mean),
    r_elevation     = cor(y, pool_eco$Elev_mean),
    r2_climate_poly = summary(lm(y ~ poly(pool_eco$Elev_mean, 2) + poly(pool_eco$Tot_prec_mean, 2)))$r.squared,
    r2_ecoregion    = summary(lm(y ~ pool_eco$Ecoregion))$r.squared,
    frac_within_eco = sd(y - ave(y, pool_eco$Ecoregion)) / sd(y),
    r_with_pool_point = cor(y, pool_eco$pool_point)
  )
}) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))

cat("\n== Pool metrics: climate-redundancy check ==\n")
cat("(want a range-rarity metric to have LOWER r2_climate_poly and HIGHER frac_within_eco than pool_point)\n")
print(pool_diag, width = Inf)

cat("\n== Each metric per Ecoregion (mean) ==\n")
pool_eco %>%
  summarize(across(all_of(pool_metrics), ~ round(mean(.x), 2)), .by = Ecoregion) %>%
  arrange(Ecoregion) %>% print(width = Inf)

# Figure ----

pool_lab <- c(pool_point = "Species count", pool_we = "Weighted endemism (Sum 1/range)",
              pool_wes = "Weighted endemism (Sum 1/sqrt range)", pool_cwe = "Corrected weighted endemism",
              pool_rr = "Restricted-range species (< 50,000 km2)")

p_pool <- pool_eco %>%
  pivot_longer(all_of(pool_metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(recode(metric, !!!pool_lab), levels = unname(pool_lab))) %>%
  ggplot(aes(fct_reorder(Ecoregion, value), value)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Potential species-pool metrics by ecoregion (farm point)",
       subtitle = "Species whose Ayerbe range + elevational range include the farm. Range sizes from AVONET1 (BirdLife).") +
  theme_minimal(11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave("Figures/Species_pool_by_ecoregion.png", p_pool, width = 11, height = 6.5, bg = "white")
print(p_pool)
