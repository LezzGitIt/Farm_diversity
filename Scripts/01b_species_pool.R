# Regional species-pool covariate from Ayerbe range maps ----

### Builds a per-farm "regional species pool" covariate: the number of bird species whose Colombian range (Ayerbe et al. range maps) overlaps a buffer around the farm. This is the one thing `Ecoregion` carries that elevation + climate + canopy do not (the biogeographic species pool) -- see Scripts/dag.R and Scripts/qmd/Species_pool_proposal.qmd. Intended as a continuous covariate for the "climate" model version in Scripts/04_farm_mgmt_mod.R, so that version can stand without the `Ecoregion` factor.

### Two counts per farm, at three buffer radii (25 / 50 / 100 km):
###   * pool_all       -- every species whose range polygon intersects the buffer
###   * pool_elevband  -- only those species whose Colombian elevational range also overlaps the farm's elevation +/- 300 m, so the covariate is a *biogeographic residual* rather than elevational turnover (which the model's Elev term already handles)

### Frozen covariate, written to Data/Farm_species_pool.csv and treated as raw input (like Data/Farm_covariates.csv). Re-run only when the range maps or the farm set change. The combined range layer and per-species elevation envelope are cached in Data/Geospatial/.

# Setup ----
library(tidyverse)
library(sf)
library(terra)

sf::sf_use_s2(TRUE)
dir.create("Data/Geospatial", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

ayerbe_dir <- "../Geospatial_data/Ayerbe_shapefiles_1890spp"
dem_path   <- "../Geospatial_data/Environmental/elevation/COL_elv_msk.tif"

ranges_cache <- "Data/Geospatial/Ayerbe_ranges.gpkg"

radii_km <- c(25, 50, 100)
elev_band_m <- 300
metric_crs <- 3116   # MAGNA-SIRGAS / Colombia Bogota -- national, metres

# Combined range layer (cached) ----

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
  ## keep only polygonal geometry, cast to a single type for the gpkg
  ranges <- ranges[st_is(ranges, c("POLYGON", "MULTIPOLYGON")), ] %>%
    st_cast("MULTIPOLYGON", warn = FALSE)
  st_write(ranges, ranges_cache, delete_dsn = TRUE, quiet = TRUE)
  message("wrote ", ranges_cache, ": ", n_distinct(ranges$species), " species, ", nrow(ranges), " features")
}

# Farms ----

fc <- read_csv("Data/Farm_covariates.csv", show_col_types = FALSE) %>%
  mutate(Id_gcs = as.character(Id_gcs))
farms <- st_as_sf(fc, coords = c("Long_mean", "Lat_mean"), crs = 4326)

# Project farms + ranges to a metric CRS ----

ranges_m <- st_transform(ranges, metric_crs)
farms_m  <- st_transform(farms, metric_crs)

# Per-farm elevation-band region (DEM cells within farm elevation +/- 300 m, inside the 100 km buffer) ----

## An earlier version used each species' global elevational min/max from the DEM, but Ayerbe range polygons are single blobs, so a montane species' [min, max] spans the whole gradient and the +/- 300 m filter removes nothing. Instead: build the actual elevation-band terrain around each farm and count only species whose range overlaps it.

dem <- rast(dem_path)
band_poly <- map(seq_len(nrow(farms)), function(i) {
  fe <- fc$Elev_mean[i]
  buf100 <- st_buffer(farms_m[i, ], max(radii_km) * 1000) %>% st_transform(st_crs(dem))
  dcrop <- crop(dem, vect(buf100))
  band <- (dcrop >= fe - elev_band_m) & (dcrop <= fe + elev_band_m)
  band[band == 0] <- NA
  bp <- tryCatch(st_as_sf(as.polygons(band, dissolve = TRUE)), error = function(e) NULL)
  if (is.null(bp) || nrow(bp) == 0) return(NULL)
  st_transform(st_union(st_geometry(bp)), metric_crs)
})

# Pool counts per farm x radius ----

pool_long <- map(radii_km, function(rk) {
  map_dfr(seq_len(nrow(farms)), function(i) {
    buf_i <- st_buffer(farms_m[i, ], rk * 1000)
    hit_i <- ranges_m[st_intersects(buf_i, ranges_m)[[1]], ]
    n_all <- n_distinct(hit_i$species)
    bp <- band_poly[[i]]
    if (is.null(bp) || length(bp) == 0) {
      n_band <- NA_integer_
    } else {
      region_i <- suppressWarnings(st_intersection(st_geometry(buf_i), bp))
      n_band <- if (length(region_i) == 0) 0L
                else n_distinct(hit_i$species[lengths(st_intersects(hit_i, region_i)) > 0])
    }
    tibble(Id_gcs = farms$Id_gcs[i], radius_km = rk, pool_all = n_all, pool_elevband = n_band)
  })
}) %>% list_rbind()

message(pool_long %>% filter(radius_km == max(radii_km)) %>% pull(pool_all) %>% max(),
        " species in the largest single-farm pool; ",
        round(100 * mean(pool_long$pool_elevband / pool_long$pool_all, na.rm = TRUE)),
        "% survive the elevation-band filter on average")

Farm_species_pool <- pool_long %>%
  select(Id_gcs, radius_km, pool_all, pool_elevband) %>%
  pivot_wider(names_from = radius_km, values_from = c(pool_all, pool_elevband),
              names_glue = "{.value}_{radius_km}k") %>%
  left_join(fc %>% select(Id_gcs, Nombre_finca, Ecoregion, Elev_mean), by = "Id_gcs")

write_csv(Farm_species_pool %>% select(-Nombre_finca, -Ecoregion, -Elev_mean),
          "Data/Farm_species_pool.csv")

# Diagnostics ----

cat("\n== Regional species pool per farm ==\n")
Farm_species_pool %>%
  select(Ecoregion, starts_with("pool_")) %>%
  summarize(across(starts_with("pool_"), ~ sprintf("%.0f (%.0f-%.0f)", mean(.x), min(.x), max(.x))), .by = Ecoregion) %>%
  arrange(Ecoregion) %>%
  print(width = Inf)

r2_eco <- Farm_species_pool %>%
  select(Ecoregion, starts_with("pool_")) %>%
  pivot_longer(-Ecoregion, names_to = "metric") %>%
  summarize(r2_ecoregion = round(summary(lm(value ~ Ecoregion))$r.squared, 2), .by = metric)
cat("\nR^2 of  pool ~ Ecoregion  (should be high -- the pool is meant to carry the between-region signal):\n")
print(r2_eco)

## how much of each metric survives within-ecoregion (the part that could add signal beyond the Ecoregion factor)
cat("\nwithin-ecoregion SD as a fraction of total SD, by metric:\n")
Farm_species_pool %>%
  select(Ecoregion, starts_with("pool_")) %>%
  pivot_longer(-Ecoregion, names_to = "metric") %>%
  summarize(frac_within = round(sd(value - ave(value, Ecoregion)) / sd(value), 2), .by = metric) %>%
  print()

# Figure ----

p_pool <- Farm_species_pool %>%
  select(Id_gcs, Ecoregion, pool_all_50k, pool_elevband_50k) %>%
  pivot_longer(c(pool_all_50k, pool_elevband_50k), names_to = "metric", values_to = "pool") %>%
  mutate(metric = recode(metric, pool_all_50k = "All species (50 km)",
                         pool_elevband_50k = "Elevation-band matched (50 km)")) %>%
  ggplot(aes(fct_reorder(Ecoregion, pool), pool)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = NULL, y = "Regional species pool (# species)",
       title = "Range-map species pool by ecoregion",
       subtitle = "50 km buffer around each farm; elevation-band = species whose elevational range overlaps farm elevation +/- 300 m") +
  theme_minimal(11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("Figures/Species_pool_by_ecoregion.png", p_pool, width = 10, height = 4.5, bg = "white")
print(p_pool)
