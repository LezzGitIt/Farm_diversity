# Regional species-pool covariate from Ayerbe range maps + Suarez-Castro elevational limits ----

### Per-farm "how many bird species could occur here" -- a potential-species-pool covariate. Range maps and published elevational limits are external to the point-count data, so this is NOT circular: it predicts the regional avifauna available to a farm without any knowledge of what was surveyed there. It is the one thing `Ecoregion` carries that elevation + climate + canopy in Scripts/04_farm_mgmt_mod.R do not (biogeographic species pool); see Scripts/dag.R and Scripts/qmd/Species_pool_proposal.qmd.

### Two-step filter per farm:
###   1. GEOSPATIAL -- the species' Ayerbe Colombian range polygon contains the farm (point-in-polygon; also computed at a 10 km and 25 km buffer for robustness to farm-centroid / polygon-edge imprecision).
###   2. ELEVATIONAL -- the farm's elevation falls within the species' elevational range +/- `elev_margin_m`, using Minimum/Maximum elevation from Suarez-Castro et al. (2024) AOH table S3 (1,652 Colombian species). Species with no elevational limit (~17%, mostly Nearctic migrants and waterbirds) pass this step.

### `pool_point` (no buffer) is the primary metric. Frozen covariate -> Data/Farm_species_pool.csv, treated as raw input like Data/Farm_covariates.csv. The combined Ayerbe range layer is cached in Data/Geospatial/Ayerbe_ranges.gpkg (the 1,890 shapefile reads are slow on a cold OneDrive).

# Setup ----
library(tidyverse)
library(sf)

sf::sf_use_s2(TRUE)
dir.create("Data/Geospatial", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

ayerbe_dir   <- "../Geospatial_data/Ayerbe_shapefiles_1890spp"
suarez_path  <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Datasets_external/Elev_ranges/Suarez_castro_AOH_birds_table_S3_V3.csv"
ranges_cache <- "Data/Geospatial/Ayerbe_ranges.gpkg"

buffer_km    <- c(0, 10, 25)   # 0 = point-in-polygon (primary)
elev_margin_m <- 250
metric_crs   <- 3116           # MAGNA-SIRGAS / Colombia Bogota -- national, metres

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
farms   <- st_as_sf(fc, coords = c("Long_mean", "Lat_mean"), crs = 4326)
ranges_m <- st_transform(ranges, metric_crs)
farms_m  <- st_transform(farms, metric_crs)

elo <- set_names(elev_lookup$elev_lo, elev_lookup$species)
ehi <- set_names(elev_lookup$elev_hi, elev_lookup$species)

# Pool per farm x buffer ----

pool_long <- map(buffer_km, function(bk) {
  geom <- if (bk == 0) st_geometry(farms_m) else st_buffer(st_geometry(farms_m), bk * 1000)
  map_dfr(seq_len(nrow(farms)), function(i) {
    sp_geo <- unique(ranges_m$species[st_intersects(geom[i], ranges_m)[[1]]])   # step 1
    fe <- fc$Elev_mean[i]
    lo <- elo[sp_geo]; hi <- ehi[sp_geo]
    in_band <- is.na(lo) | (hi >= fe - elev_margin_m & lo <= fe + elev_margin_m) # step 2 (NA elev -> keep)
    tibble(Id_gcs = fc$Id_gcs[i], buffer_km = bk,
           pool_geo = length(sp_geo), pool = sum(in_band))
  })
}) %>% list_rbind()

Farm_species_pool <- pool_long %>%
  select(Id_gcs, buffer_km, pool_geo, pool) %>%
  pivot_wider(names_from = buffer_km, values_from = c(pool_geo, pool),
              names_glue = "{.value}_{buffer_km}km") %>%
  rename_with(~ str_replace(.x, "_0km$", "_point"))

write_csv(Farm_species_pool, "Data/Farm_species_pool.csv")

# Diagnostics ----

pool_eco <- Farm_species_pool %>%
  left_join(fc %>% select(Id_gcs, Ecoregion, Elev_mean, Tot_prec_mean), by = "Id_gcs")

cat("\n== Species pool per farm (mean [range]) ==\n")
pool_eco %>%
  summarize(across(c(pool_point, pool_10km, pool_25km, pool_geo_point),
                   ~ sprintf("%.0f [%.0f-%.0f]", mean(.x), min(.x), max(.x))), .by = Ecoregion) %>%
  arrange(Ecoregion) %>% print(width = Inf)

cat("\nelevation filter effect (pool / pool_geo):\n")
pool_eco %>%
  summarize(point = round(mean(pool_point / pool_geo_point), 2),
            b25 = round(mean(pool_25km / pool_geo_25km), 2)) %>% print()

cat("\npool ~ Ecoregion R^2, and correlation with precipitation / elevation:\n")
tibble(metric = c("pool_point", "pool_10km", "pool_25km")) %>%
  mutate(
    r2_ecoregion = map_dbl(metric, ~ summary(lm(pool_eco[[.x]] ~ pool_eco$Ecoregion))$r.squared),
    r_precip     = map_dbl(metric, ~ cor(pool_eco[[.x]], pool_eco$Tot_prec_mean)),
    r_elevation  = map_dbl(metric, ~ cor(pool_eco[[.x]], pool_eco$Elev_mean)),
    frac_within_eco = map_dbl(metric, ~ sd(pool_eco[[.x]] - ave(pool_eco[[.x]], pool_eco$Ecoregion)) / sd(pool_eco[[.x]]))
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% print()

# Figure ----

p_pool <- pool_eco %>%
  select(Id_gcs, Ecoregion, pool_point, pool_25km) %>%
  pivot_longer(c(pool_point, pool_25km), names_to = "metric", values_to = "pool") %>%
  mutate(metric = recode(metric, pool_point = "Farm point (primary)", pool_25km = "25 km buffer")) %>%
  ggplot(aes(fct_reorder(Ecoregion, pool), pool)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = NULL, y = "Potential species pool (# species)",
       title = "Range-map species pool by ecoregion",
       subtitle = "Species whose Ayerbe range includes the farm and whose elevational range (Suarez-Castro 2024) overlaps farm elevation +/- 250 m") +
  theme_minimal(11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("Figures/Species_pool_by_ecoregion.png", p_pool, width = 10, height = 4.5, bg = "white")
print(p_pool)
