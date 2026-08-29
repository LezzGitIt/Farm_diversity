# Regional species-pool covariate from Ayerbe range maps + Suarez-Castro elevational limits ----

### Per-farm "how many bird species could occur here" -- a potential-species-pool covariate. Range maps and published elevational limits are external to the point-count data, so this is NOT circular: it predicts the regional avifauna available to a farm without any knowledge of what was surveyed there. It is the one thing `Ecoregion` carries that elevation + climate + canopy in Scripts/04_farm_mgmt_mod.R do not (biogeographic species pool); see Scripts/dag.R and Scripts/qmd/Species_pool_proposal.qmd.

### Two-step filter per farm:
###   1. GEOSPATIAL -- the species' Ayerbe Colombian range polygon contains the farm (point-in-polygon).
###   2. ELEVATIONAL -- the farm's elevation falls within the species' elevational range +/- `elev_margin_m`, using Minimum/Maximum elevation from Suarez-Castro et al. (2024) AOH table S3 (1,652 Colombian species). Species with no elevational limit (~17%, mostly Nearctic migrants and waterbirds) pass this step.

### `pool_point` -> Data/Farm_species_pool.csv, treated as raw input like Data/Farm_covariates.csv. The combined Ayerbe range layer is cached in Data/Geospatial/Ayerbe_ranges.gpkg (the 1,890 shapefile reads are slow on a cold OneDrive).

### OUTCOME (2026-08-28): `pool_point ~ poly(Elev,2) + poly(precip,2)` R^2 = 0.89 -- the pool is ~redundant with the climate terms already in 04. The decisive test (add pool_point_z to the climate models; Derived/Excels/Species_pool_decisive_test.csv) confirmed it: the pasture/water coefficients move by <= 0.006 and pool_point_z is null with a huge CrI. So this covariate is NOT added to 04 -- it's kept for the record. The ecoregion<->climate difference is compositional biogeography or Ecoregion over-adjustment, not a missing richness axis; an endemism-weighted or eBird measure would be the next thing to try.

# Setup ----
library(tidyverse)
library(sf)

sf::sf_use_s2(TRUE)
dir.create("Data/Geospatial", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

ayerbe_dir   <- "../Geospatial_data/Ayerbe_shapefiles_1890spp"
suarez_path  <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Elev_ranges/Suarez_castro_AOH_birds_table_S3_V3.csv"
ranges_cache <- "Data/Geospatial/Ayerbe_ranges.gpkg"

elev_margin_m <- 250

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

# Elevational limits: Suarez-Castro et al. (2024) AOH table S3 ----

suarez <- read_csv(suarez_path, show_col_types = FALSE) %>%
  select(Scientific.Name, Name.Clements.eBird., Minimum.elevation, Maximum.elevation)

## one row per name (BirdLife or Clements), so an Ayerbe species can match on either
suarez_by_name <- bind_rows(
  suarez %>% transmute(name = Scientific.Name, elev_lo = Minimum.elevation, elev_hi = Maximum.elevation),
  suarez %>% transmute(name = Name.Clements.eBird., elev_lo = Minimum.elevation, elev_hi = Maximum.elevation)
) %>%
  filter(!is.na(name), !is.na(elev_lo)) %>%
  distinct(name, .keep_all = TRUE)

elev_lookup <- tibble(species = sort(unique(ranges$species))) %>%
  left_join(suarez_by_name, by = c("species" = "name"))
n_with_elev <- sum(!is.na(elev_lookup$elev_lo))
message(n_with_elev, " of ", nrow(elev_lookup), " Ayerbe species matched to a Suarez-Castro elevational range (",
        round(100 * n_with_elev / nrow(elev_lookup)), "%); the rest pass the elevation filter unfiltered")

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
  in_band <- is.na(lo) | (hi >= fe - elev_margin_m & lo <= fe + elev_margin_m)   # step 2 (NA elev -> keep)
  tibble(Id_gcs = fc$Id_gcs[i], pool_geo_point = length(sp_geo), pool_point = sum(in_band))
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

cat("\npool_point diagnostics:\n")
tibble(
  r2_ecoregion    = summary(lm(pool_point ~ Ecoregion, pool_eco))$r.squared,
  r_precip        = cor(pool_eco$pool_point, pool_eco$Tot_prec_mean),
  r_elevation     = cor(pool_eco$pool_point, pool_eco$Elev_mean),
  r2_climate_poly = summary(lm(pool_point ~ poly(Elev_mean, 2) + poly(Tot_prec_mean, 2), pool_eco))$r.squared,
  frac_within_eco = sd(pool_eco$pool_point - ave(pool_eco$pool_point, pool_eco$Ecoregion)) / sd(pool_eco$pool_point)
) %>% mutate(across(everything(), ~ round(.x, 2))) %>% print(width = Inf)

# Figure ----

p_pool <- pool_eco %>%
  ggplot(aes(fct_reorder(Ecoregion, pool_point), pool_point)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.7) +
  labs(x = NULL, y = "Potential species pool (# species)",
       title = "Range-map species pool by ecoregion (farm point)",
       subtitle = "Species whose Ayerbe range includes the farm and whose elevational range (Suarez-Castro 2024) overlaps farm elevation +/- 250 m") +
  theme_minimal(11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("Figures/Species_pool_by_ecoregion.png", p_pool, width = 7, height = 4.5, bg = "white")
print(p_pool)
