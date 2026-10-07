# Woody vegetation canopy cover at multiple spatial scales around each assemblage ----

### Scale-of-effect preparation. An assemblage -- [data collector . farm . year-group . season] -- is the unit the iNEXT diversity estimates are computed at. For each assemblage this takes the convex hull of the point counts that assemblage sampled, buffers that hull outward at 18 radii (200 m to 2000 m in 200 m steps, then 3000 m to 10 km in 1 km steps), and extracts the mean woody-vegetation canopy cover inside each buffered hull from the WVCC raster year matching the assemblage's survey year. Every assemblage has exactly one Id_scr and one survey year, so the output feeds Scripts/02b_scale_of_effect.R directly with no farm / group roll-up.

### Self-contained: needs only the point-count coordinates (`Site_covs.csv`), the assemblage membership (`Event_covs.csv`), and the Colombia Woody Vegetation Structure and Change rasters (~25 m, Zenodo 18154841; canopy cover of trees / shrubs > 2 m) in `../Geospatial_data/Environmental/`. The heavy step is clipping the ~500 MB national rasters; outputs are frozen into `Data/Geospatial/` so `08` and any re-analysis do not need them:
###   * `Assemblage_hulls.gpkg`           -- convex hull of each assemblage's point counts
###   * `canopy_rasters/farm_<Id_scr>_<year>.tif` -- national WVCC clipped to each farm's 5.5 km extent (assemblages of the same farm x year share the clip)
###   * `Canopy_by_scale_assemblage.csv`  -- assemblage x radius -> mean canopy cover

# Setup ----
library(tidyverse)
library(terra)

source("Scripts/Model_fns.R")   # WRANGLING_EXCELS, tax_div_path()

Sys.setenv(GDAL_NUM_THREADS = "ALL_CPUS")
terraOptions(memfrac = 0.6, progress = 0)

wvcc_dir <- "../Geospatial_data/Environmental"

out_dir <- "Data/Geospatial"
raster_out_dir <- file.path(out_dir, "canopy_rasters")
dir.create(raster_out_dir, recursive = TRUE, showWarnings = FALSE)

## Buffered-hull radii (metres): the model in 08 asks which of these best explains diversity. Extended to 10 km because the topography-adjusted scale-of-effect curve had not peaked by 5 km.
radii_m <- c(seq(200, 2000, by = 200), seq(3000, 10000, by = 1000))

## Crop the national rasters to this much around a farm's points -- 500 m past the widest radius so edge pixels are never clipped
clip_buffer_m <- 10500

## WVCC raster years available on disk
raster_years <- c(2013, 2014, 2016, 2017, 2019, 2022, 2024)
nearest_raster_year <- function(year) raster_years[which.min(abs(raster_years - year))]

# Point counts and assemblage membership ----

pc_coords <- read_csv(file.path(WRANGLING_EXCELS, "Site_covs.csv"), show_col_types = FALSE) %>%
  mutate(Id_scr = as.character(Id_scr)) %>%
  select(Id_survey_no_dc, Id_scr, Long, Lat)

## The assemblages we need canopy for: those with a bird diversity estimate (the fuller all_farms export). This also drops the 9 Otun Quimbaya mature-forest reference point counts (`MB-R-OQ_*`), which are not on an SCR farm and are absent from Site_covs.
target_assemblages <- read_csv(
  tax_div_path("all_farms"), show_col_types = FALSE
) %>% distinct(Assemblage)

## One Event_covs row per point-count survey; build the assemblage id and snap the survey year to the nearest WVCC raster year
assemblage_surveys <- read_csv(file.path(WRANGLING_EXCELS, "Event_covs.csv"), show_col_types = FALSE) %>%
  filter(!is.na(Year)) %>%
  left_join(pc_coords %>% select(Id_survey_no_dc, Id_scr), by = "Id_survey_no_dc") %>%
  filter(!is.na(Id_scr)) %>%
  mutate(
    Id_scr = as.character(Id_scr),
    Assemblage = str_replace_all(paste(Uniq_db, Id_scr, Year_grp, Season, sep = "."), " |-", "_"),
    raster_year = map_dbl(Year, nearest_raster_year)
  ) %>%
  semi_join(target_assemblages, by = "Assemblage")

## Assemblage -> its distinct point counts
assemblage_pcs <- assemblage_surveys %>%
  distinct(Assemblage, Id_scr, Id_survey_no_dc, raster_year)

## Assemblage -> its farm and raster year (verified one of each per assemblage; Id_scr is part of the key, and no assemblage's surveys span two WVCC years)
assemblage_meta <- assemblage_pcs %>% distinct(Assemblage, Id_scr, raster_year)
stopifnot(nrow(assemblage_meta) == n_distinct(assemblage_meta$Assemblage))

cat(n_distinct(assemblage_pcs$Assemblage), "assemblages across",
    n_distinct(assemblage_pcs$Id_scr), "farms and",
    n_distinct(assemblage_pcs$Id_survey_no_dc), "point counts.\n")

# Per-assemblage convex hulls ----

pts <- vect(pc_coords, geom = c("Long", "Lat"), crs = "EPSG:4326")

assemblage_ids <- sort(unique(assemblage_pcs$Assemblage))

## Convex hull of an assemblage's point counts; 1-2 point assemblages give a point / line, so pad by 0.5 m (negligible) to keep every hull a polygon
assemblage_hull_of <- function(assemblage_id) {
  pc_ids <- assemblage_pcs %>% filter(Assemblage == assemblage_id) %>% pull(Id_survey_no_dc)
  hull <- convHull(pts[pts$Id_survey_no_dc %in% pc_ids, ])
  if (geomtype(hull) != "polygons") hull <- buffer(hull, width = 0.5)
  hull$Assemblage <- assemblage_id
  hull
}
assemblage_hulls <- vect(map(assemblage_ids, assemblage_hull_of))
writeVector(assemblage_hulls, file.path(out_dir, "Assemblage_hulls.gpkg"), overwrite = TRUE)

# Clip a needed year's national raster to a farm's 5.5 km extent (cached) ----

## Cache unit is [farm x raster year]: assemblages of the same farm and year share a clip, and a farm's clip (all its points + 5.5 km) contains every buffered hull of its assemblages. Returns instantly if the clip already exists on disk.
clip_farm_year <- function(id_scr, year) {
  out_tif <- file.path(raster_out_dir, sprintf("farm_%s_%d.tif", id_scr, year))
  if (file.exists(out_tif)) return(out_tif)
  message("clip: farm ", id_scr, " / ", year)
  farm_pts <- pts[pts$Id_scr == id_scr, ]
  clip_extent <- ext(buffer(farm_pts, width = clip_buffer_m))
  national <- rast(file.path(wvcc_dir, sprintf("Colombia_WVCC_%d.tif", year)))
  clipped <- crop(national, clip_extent, snap = "out")
  writeRaster(clipped, out_tif, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  out_tif
}

# Extract mean canopy cover within each buffered hull (cached to CSV) ----

canopy_csv <- file.path(out_dir, "Canopy_by_scale_assemblage.csv")

extract_assemblage <- function(assemblage_id) {
  message("extract: ", assemblage_id)
  meta <- assemblage_meta %>% filter(Assemblage == assemblage_id)
  clip_farm_year(meta$Id_scr, meta$raster_year)
  hull <- assemblage_hulls[assemblage_hulls$Assemblage == assemblage_id, ]
  canopy_r <- rast(file.path(raster_out_dir, sprintf("farm_%s_%d.tif", meta$Id_scr, meta$raster_year)))
  map(radii_m, function(radius) {
    buffered_hull <- buffer(hull, width = radius)
    mean_cover <- terra::extract(canopy_r, buffered_hull, fun = "mean", exact = TRUE,
                                 ID = FALSE, na.rm = TRUE)[[1]]
    tibble(Assemblage = assemblage_id, Id_scr = meta$Id_scr, raster_year = meta$raster_year,
           radius_m = radius, canopy_cover = mean_cover)
  }) %>% list_rbind()
}

## The extraction loop is the slow step, so the CSV is an incremental cache: only assemblages not already in it are extracted, then appended. The cache is keyed by assemblage id alone, so if an existing assemblage's point counts change (new or dropped points upstream), delete that assemblage's rows from the CSV to re-extract it; delete the whole CSV (or change the radii) to force a full rebuild.
Canopy_cached <- if (file.exists(canopy_csv)) {
  read_csv(canopy_csv, show_col_types = FALSE) %>% mutate(Id_scr = as.character(Id_scr))
} else {
  tibble(Assemblage = character(), Id_scr = character(), raster_year = double(), radius_m = double(), canopy_cover = double())
}
to_extract <- setdiff(assemblage_ids, Canopy_cached$Assemblage)
message(length(assemblage_ids) - length(to_extract), " assemblages loaded from the canopy cache; extracting ", length(to_extract), ".")
Canopy_by_scale <- bind_rows(Canopy_cached, map(to_extract, extract_assemblage) %>% list_rbind()) %>%
  filter(Assemblage %in% assemblage_ids) %>%
  arrange(Assemblage, radius_m)
if (length(to_extract) > 0) write_csv(Canopy_by_scale, canopy_csv)

cat("\nCanopy_by_scale:", nrow(Canopy_by_scale), "rows;",
    sum(is.na(Canopy_by_scale$canopy_cover)), "missing canopy values.\n")
print(
  Canopy_by_scale %>%
    summarize(mean = mean(canopy_cover, na.rm = TRUE),
              min = min(canopy_cover, na.rm = TRUE),
              max = max(canopy_cover, na.rm = TRUE), .by = radius_m)
)

# Diagnostic plot: one assemblage's point counts, hull, and buffered rings on the canopy raster ----

## Change this to inspect a different assemblage; default is a 12-point Bajo Magdalena assemblage with a strong canopy gradient across radii
example_assemblage <- "Gaica_mbd.1053.16_17.Early"

library(tidyterra)

ex_meta <- assemblage_meta %>% filter(Assemblage == example_assemblage)
ex_pc_ids <- assemblage_pcs %>% filter(Assemblage == example_assemblage) %>% pull(Id_survey_no_dc)
ex_pts <- pts[pts$Id_survey_no_dc %in% ex_pc_ids, ]
ex_hull <- assemblage_hulls[assemblage_hulls$Assemblage == example_assemblage, ]
ex_rings <- vect(map(radii_m, function(radius) {
  ring <- buffer(ex_hull, width = radius)
  ring$radius_m <- radius
  ring
}))

## Crop the farm-year canopy raster to just the widest ring for a tight view (ensures the clip exists even when extraction was skipped)
ex_canopy <- rast(clip_farm_year(ex_meta$Id_scr, ex_meta$raster_year)) %>%
  crop(ext(ex_rings[ex_rings$radius_m == max(radii_m), ]), snap = "out")

p_example <- ggplot() +
  geom_spatraster(data = ex_canopy) +
  scale_fill_viridis_c(name = "Canopy\ncover (%)", limits = c(0, 100), na.value = "transparent") +
  geom_spatvector(data = ex_rings, fill = NA, colour = "white", linewidth = 0.3) +
  geom_spatvector(data = ex_hull, fill = NA, colour = "red", linewidth = 0.7) +
  geom_spatvector(data = ex_pts, colour = "red", size = 1.3) +
  coord_sf(expand = FALSE) +
  labs(
    title = paste0("Buffered-hull scales: ", example_assemblage, " (WVCC ", ex_meta$raster_year, ")"),
    subtitle = paste0(nrow(ex_pts), " point counts + convex hull (red); ",
                      length(radii_m), " rings, ", min(radii_m), " m to ", max(radii_m) / 1000, " km (white)")
  ) +
  theme_minimal()
ggsave("Figures/Scale_effect_example_buffers.png", p_example, bg = "white", width = 8.5, height = 8)
print(p_example)
