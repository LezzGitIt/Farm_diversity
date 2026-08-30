# Regional species-pool covariate from Ayerbe range maps + Suarez-Castro elevational limits ----

### Per-farm "how many bird species could occur here" -- a potential-species-pool covariate. Range maps and published elevational limits are external to the point-count data, so this is NOT circular: it predicts the regional avifauna available to a farm without any knowledge of what was surveyed there. It is the one thing `Ecoregion` carries that elevation + climate + canopy in Scripts/04_farm_mgmt_mod.R do not (biogeographic species pool); see Scripts/dag.R and Scripts/qmd/_archive/Species_pool_proposal.qmd.

### Two-step filter per farm:
###   1. GEOSPATIAL -- the species' Ayerbe Colombian range polygon contains the farm (point-in-polygon).
###   2. ELEVATIONAL -- the farm elevation falls within the species elevational range (no buffer). Elevational limits are layered by source priority: Suarez-Castro et al. (2024) AOH table S3 (1,652 spp), then the manual book digitisations from ../Ssp-bird-data-wrangling/Scripts/03_FT_elev.R -- Hazen's read-out of the Ayerbe-Quinones (2018) field guide, then Hazen's read-out of Hilty. Species still without a limit (~13%, mostly Nearctic migrants / seabirds) pass this step unfiltered.

### `pool_point` -> Data/Farm_species_pool.csv, treated as raw input like Data/Farm_covariates.csv. The combined Ayerbe range layer is cached in Data/Geospatial/Ayerbe_ranges.gpkg (the 1,890 shapefile reads are slow on a cold OneDrive).

### OUTCOME (2026-08-29): `pool_point ~ poly(Elev,2) + poly(precip,2)` R^2 = 0.87 (unique: elev 0.29, precip 0.45, shared 0.12). Scripts/04d_species_pool_test.R tested it properly (blocks that never combine all three axes) and it does not change the management coefficients -> NOT added to 04, kept for the record. Along the way 04d + Derived/Excels/Precip_vs_elev_sensitivity.csv showed the pasture/water "signal" in the climate spec is residual PRECIPITATION confounding (Pasture_mgmt_div ~ precip r = 0.30), not a missing species-pool axis; see Project_notes.md.

# Setup ----
library(tidyverse)
library(sf)
library(readxl)

sf::sf_use_s2(TRUE)
dir.create("Data/Geospatial", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

ayerbe_dir   <- "../Geospatial_data/Ayerbe_shapefiles_1890spp"
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
       subtitle = "Species whose Ayerbe range includes the farm and whose elevational range (Suarez-Castro 2024 + guides) includes the farm elevation") +
  theme_minimal(11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("Figures/Species_pool_by_ecoregion.png", p_pool, width = 7, height = 4.5, bg = "white")
print(p_pool)
