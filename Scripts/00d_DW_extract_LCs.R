## PhD birds in silvopastoral landscapes##
# Mathilde Manual Digitization & GIS work 
# The purpose of this script is to generate files that Mathilde will need to digitize landcover from satellite imagery in QGIS, as well as to examine the files that she digitized & compare them to other publicly available spatial layers

## NOTE:: Actively developing this file, at present it is in a state of disarray

# All files can be found on Box: SCR Box -> Geodatabase -> Mathilde spatial

# Contents
# 1)
# 2)
# 3)
# 4)
# 5)
# 6)

# TO DO:
# 1) Could it make sense to just combine all these files right away? Search LCs_comb
# 2) Once we've buffered live fences, could pull data_year(s) from the surrounding polygon(s)
# 3) If there are overlaps, could union (aggregate I think?) polygons of the same landcover & data_year.

# Make sure to tidy up script before you stop working on this.
# NOTE:: I copied & pasted most of the code from section 3 to section 2, so will likely want to go back in and clean up. Probably best to make a new script for the comparison of Mathilde's data & Brandt's data.

# Load libraries & data ---------------------------------------------------
#library(geobuffer)
library(terra)
library(tidyterra)
library(sf)
library(tidyverse)
library(janitor)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(conflicted)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")
load("Rdata/the_basics_11.21.24.Rdata")

# 1. Create files to aid in digitization -----------------------------------------------------
# Create spatial file of point counts -- note that there are multiple coordinates for a single point count in a few cases

# load("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/R_Files/PhD/Manual_dig.Rdata")
# load("Rdata/Manual_dig_2.20.24.Rdata")

# Add habitat to point counts and pivot_wider to have single row per PC (there are still multiple rows for some PCs b/c there are multiple coordinates in some cases)
Bird_pcs_hab_sf <- Bird_pcs %>%
  distinct(
    Id_muestreo, Id_group, Departamento, Uniq_db,
    Id_gcs, Nombre_finca_mixed, Nombre_finca, Habitat_og,
    Habitat_ut, Longitud_decimal, Latitud_decimal
  ) %>%
  group_by(Id_muestreo) %>% # Habitat_og, Habitat_ut
  mutate(Habitat = row_number()) %>% # Create numbers to add to columns
  pivot_wider(
    names_from = Habitat,
    values_from = c(Habitat_ut, Habitat_og)
  ) %>%
  ungroup() %>%
  st_as_sf(coords = c("Longitud_decimal", "Latitud_decimal"), crs = 4326)

# Combine Habitat_ut & Habitat_og into a single column with unite()
endings <- 1:7
endings_ <- c(paste0("_", 1:7))
for (i in 1:length(endings)) {
  print(i)
  Pc_locs_sf <- Bird_pcs_hab_sf %>%
    unite(
      col = "Hab_UT_OG",
      ends_with(endings_[i]), sep = "_"
    ) %>%
    rename_with(~ paste0("Hab_UT_OG", endings[i]), Hab_UT_OG)
}

Bird_pcs_hab_sf %>%
  distinct(Id_muestreo, Uniq_db) %>%
  nrow() - 54 # %>% count(Uniq_db)

# >Buffers ---------------------------------------------------------------
# Create buffer sizes vector and prepare loop
# NOTE:: Even though in terra::buffer the argument is called 'width' it still refers to the radius
buffer_dists <- c(25, 50, 100, 300, 500, 1000, 5000)
buffers <- union <- v_union <- list()
area_km2 <- as.numeric()

# Create buffers of varying sizes to determine how much area would need to be digitized.. Time is limited so need to pick a reasonable buffer size for Mathilde
sf_use_s2(TRUE) # Otherwise buffer_dists gives weird error about arc degrees

# Function geobuffer_pts to make a geodesic buffer doesn't seem to work anymore in 2024. See github for previous versions of code as I've deleted several small things (including what is needed to make the st_union().

cents_sv <- Pc_locs_sf %>% vect() %>% # sv = spatvector
  project("epsg:4686")
for (i in 1:length(buffer_dists)) {
  print(i)
  buffers[[i]] <- cents_sv %>%
    buffer(width = buffer_dists[i]) %>%
    st_as_sf()
  # 2 options here.. Can combine overlapping polygons st_union -> (direct to, I think) st_cast, or have separate polygons for each Id_group buffers -> group_by -> summarize -> st_cast
  union[[i]] <- buffers[[i]] %>% # st_union() %>% #st_union merges all overlapping polygons
    group_by(Id_group) %>%
    summarize(geometry = st_union(geometry)) %>%
    st_cast("MULTIPOLYGON")
  v_union[[i]] <- st_make_valid(union[[i]])
  area_km2[[i]] <- sum(st_area(v_union[[i]])) / 1000000
}
warnings() # If 'polygon from first part only' this is important means it's only taking a single polygon when not overlapping

# Nested list of buffers #DELETE UNNECESSARY
buffersLL <- map(.x = list(buffers, union), .f = \(buffL)
                 setNames(buffL, paste0("b", buffer_dists)))

# Name buffers & unioned buffers
names(buffers) <- paste0("b", buffer_dists)
names(union) <- paste0("b", buffer_dists)

# Convert buffers & unioned buffers to vector and set crs to MAGNA-SIRGAS
buffers_ms <- map(buffers, .f = \(x) vect(x) %>% project("epsg:4686"))
union_ms <- map(union, .f = \(x) vect(x) %>% project("epsg:4686"))

# I think numbers have differed in terms of area at least in part b/c initially you did dissolve all polygons (this provides the accurate area)
buffers_df <- data.frame(buffer_dists, area_km2) %>%
  mutate(area3times = area_km2 * 3)
sapply(buffers_df, round, 2)

# >Export files ----------------------------------------------------------
# Create kmz & shapefiles of buffers + points. WILL WANT TO CHANGE PROJECTION
Buffers_cipav <- buffers_ms$b25 %>% filter(Uniq_db == "CIPAV MBD")

# Put items that need to be exported in a list
files.exp <- list(PC_Cents = cents_sv, Buffers_cipav = Buffers_cipav, Buffers_50m = buffers_ms$b50, Buffers_300m = union_ms$b300)
files.exp.sf <- map(files.exp, \(x) st_as_sf(x))

# Rename Centroids file & arrange so they're in order -- if get error may need to tibble and then rename
files.exp.sf$PC_Cents <- files.exp.sf$PC_Cents %>%
  # tibble() %>%
  rename(name = Id_muestreo) %>%
  arrange(name) # %>%
# st_as_sf()

dir <- "Intermediate_products_geospatial/"
filepath <- c("kml", "shp")
filetype <- c("kml", "ESRI Shapefile")

# Create folders if necessary
map(names(files.exp.sf), \(x) dir.create(path = paste0(dir, "shp/", x)))
# Sf
map(names(files.exp.sf), \(name)
    map2(
      .x = filepath, .y = filetype,
      .f = \(path, type) {
        path_shp <- paste0(dir, path, "/", name, "/", name, ".", path)
        path_kml <- paste0(dir, path, "/", name, ".", path)
        st_write(
          obj = files.exp.sf[[name]],
          dsn = if (path == "shp") path_shp else path_kml,
          driver = type,
          layer = name
        )
      }
    ))

# >Example farm figure ----------------------------------------------------
Bird_pcs %>%
  filter(Nombre_finca == "La herradura" | Id_gcs == 4121) %>%
  distinct(Id_gcs, Id_group, Nombre_finca, Nombre_finca_mixed, Latitud_decimal, Longitud_decimal)

Ex <- buffers[[4]] %>% filter(Uniq_db == "GAICA MBD" & Id_group == "G-MB-M-LH")
Ex_cents <- Pc_locs_sf %>% filter(Id_group == "G-MB-M-LH")
Ex500 <- Ex %>% st_union()
Ex500 %>% ggplot() +
  geom_sf() +
  geom_sf(data = Ex_cents) +
  ggspatial::annotation_scale(location = "br")


# >Metadata files --------------------------------------------------------
## Create file w/ appropriate dates for extraction of Planet imagery for Seth
Pc_date4 %>%
  distinct(
    Id_muestreo, Nombre_institucion,
    Uniq_db, Departamento, Nombre_finca, Ano, Mes
  ) %>%
  distinct(Departamento, Mes, Ano) %>%
  group_by(Departamento, Ano) %>%
  summarize(Month = round(mean(Mes), 0)) %>%
  mutate(
    Month = ifelse(Departamento == "Meta" & Ano == 2017, 10, Month),
    Month2 = ifelse(Month >= 3 & Month <= 6, 12, 5)
  ) %>%
  as.data.frame() # %>%
# write.xlsx("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Mentorship/Digitization_Mathilde/Digitization_Dep_month_year2.xlsx", row.names = F)

## Metadata file for Mathilde to fill out
Pc_date4 %>%
  mutate(Id_muestreo = str_split_i(Id_muestreo, "_", i = 1)) %>%
  group_by(Id_muestreo, Ano) %>%
  summarize(across(), Month.min = min(Mes), Month.max = max(Mes)) %>%
  distinct(
    Id_muestreo, Nombre_institucion,
    Uniq_db, Departamento, Nombre_finca, Ano, Month.min, Month.max
  ) %>%
  group_by(Id_muestreo) %>%
  group_modify(~ add_row(.x, Ano = 2023, .before = 0)) %>%
  View()
mutate(Closest_date = NA, Source_closest = NA, N.photos = NA, Notes = NA, Questions = NA) %>%
  as.data.frame() # %>%
# write.xlsx("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Mentorship/Digitization_Mathilde/Digitization_metadata.xlsx", row.names = F)

# save.image("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/R_Files/PhD/ColombiaPhD.Rdata")

# 2. Extract LC covs ---------------------------------------------------------
## TO TRY: If Mathilde doesn't have old shapefiles (from July), can skip over problematic polygons (e.g. topology Error) & get the indices of all problematic coordinates / polygons for Mathilde to redraw
# Aggregate polygons by lc & id group (or just lc?) first, & then try erasing? I TRIED AGGREGATING JUST BY LC_TYP AND HAD NO LUCK
# Focus on Id_group first, try to aggregate by Id_group as well
# Figure out a workflow with a single Id_group (can more easily plot), then expand
# Once Mathilde uploads original polygons, bring this in & remove all invalid polygons & then try workflow in terra or (especially sf)
# Even without the old polygons, could try workflow in sf again, using st_union()
# Go back to step 3
# Go back to remove polygon step -- are the rows in the dataframe same as # of polygons ?
# Plot LC_holes2, are there holes in it? If LC_middle is good , you could continue trying to iron out the rest of work flow, something might occur to you?
# See erase() help file, is there an example?
# Look back at old 'Misc' section

# Terra approach -------------------------------------------------------------------
## Goal:: Each point count needs landcover values in 300m radius that are temporally matched with when the point count was sampled
path <- "../../Mentorship/Digitization_Mathilde/Final_docs"
files.shp <- list.files(path = path, , pattern = "shapefile")
files.shp <- files.shp[-2] # Remove old lc_middle file

shp.lc.L <- map(.x = files.shp, \(shp)
                vect(paste0(path, "/", shp)) %>%
                  clean_names())
names(shp.lc.L) <- files.shp

# Add rownumber for easy identification & communication with Mathilde
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% mutate(poly_num = row_number())
})

# Could it make sense to just combine all these files right away? Polys can overlap, & we have data_year. Within a data_year (e.g., 2013) we shouldn't have overlap..
## NOTE:: TESTING FOR OVERLAP IS A GOOD THING TO DO, ALTHOUGH MATHILDE PROBABLY NEEDS TO FIX THE FEW PROBLEMATIC POLYGONS BEFORE WE CAN EVEN TEST THIS
LCs_comb <- pmap(
  .l = list(shp.lc.L[2], shp.lc.L[3], shp.lc.L[4]),
  \(polys, trees, lf){
    rbind(polys, trees, lf)
  }
)
sum(sapply(shp.lc.L[2:4], nrow))

# >Explore dataframes  ----------------------------------------------------------
## Polygons
map(shp.lc.L[c(2:4)], \(polys){
  names(polys)
})

map(shp.lc.L[c(2:4)], \(polys){
  table(polys$trees_perc)
})

# Trees percent
# LCs with trees_perc != NA
map(shp.lc.L[c(2:4)], \(polys){
  polys %>%
    filter(!is.na(trees_perc)) %>%
    pull(lc_typ) %>%
    table()
})

# Number of polygons
map(shp.lc.L[c(2:4)], \(polys){
  polys %>% filter(lc_typ == "distrees" & is.na(trees_perc))
})

# Forest type
map(shp.lc.L[c(2:4)], \(polys){
  polys %>%
    filter(!is.na(forest_typ)) %>%
    pull(lc_typ) %>%
    table()
})

# Crop type
map(shp.lc.L[c(2:4)], \(polys){
  polys %>%
    filter(!is.na(crop_typ)) %>%
    pull(lc_typ) %>%
    table()
})

# LC type
map(shp.lc.L[c(2:4)], \(polys){
  table(polys$lc_typ)
})

# Data_years#
# Data_year
map(shp.lc.L[c(2:4)], \(polys){
  table(polys$data_year)
})

# Data_year2
map(shp.lc.L[c(2:4)], \(polys){
  table(polys$data_year2)
})

# Data_year3
map(shp.lc.L[c(2:4)], \(polys){
  table(polys$data_year3)
})

# NOTE:: With data_year the first value is not always the earliest.. For example in some data_year1 is 2017 & data_year2 is 2016
shp.lc.L$landcovers_middle_shapefile %>% filter(data_year == 2016 & !is.na(data_year2))
shp.lc.L$landcovers_middle_shapefile %>%
  filter(data_year2 == 2016) %>%
  pull(data_year)

# Image date
map(shp.lc.L[c(2:4)], \(polys){
  table(polys$data_year)
})

# DOES IT MAKE SENSE THAT THERE ARE POLYGONS WITH DATA YEAR = 2022 IN THIS SHAPEFILE?
shp.lc.L$landcovers_middle_shapefile %>% filter(data_year == 2022)
shp.lc.L$landcovers_middle_shapefile %>% filter(data_year2 == 2022)
shp.lc.L$landcovers_middle_shapefile %>% filter(!is.na(data_year3))

# Or do it by database instead of by column across all databases
lapply(shp.lc.L$landcovers_middle_shapefile, table)
lapply(shp.lc.L$landcovers_middle_UBC_shapefile, table)
lapply(shp.lc.L$landcovers_past_shapefile, table)

# Check NAs and NaN (not a number)
lapply(shp.lc.L$landcovers_middle_shapefile, \(col) table(is.na(col)))
lapply(shp.lc.L$landcovers_middle_shapefile, \(col) table(is.nan(col)))

## Linear features
map(shp.lc.L[c(5:7)], \(lf){
  names(lf)
})

# CHECK:: Ensure tree_shrub category has only F (forest remnant), S (shrub), and T (tree)
map(shp.lc.L[c(5:7)], \(lf){
  table(lf$tree_shrub)
})


# Is this the year the linear feature was digitized?
map(shp.lc.L[c(5:7)], \(lf){
  table(lf$layer)
})

# Is this the year the linear feature was digitized? Can it be deleted?
map(shp.lc.L[c(5:7)], \(lf){
  table(lf$path)
})

# >Explore invalid polygons -----------------------------------------
# NOTE:: several polygons are not valid
map(shp.lc.L, \(shp){
  table(is.valid(shp))
})

## We lose 4 polygons in the middle shapefile
# Before makeValid
map(shp.lc.L[c(2:4)], \(polys){
  polys %>% nrow()
})

# Now only 3439 polys using makeValid
map(shp.lc.L[c(2:4)], \(polys){
  polys %>%
    makeValid() %>%
    nrow()
})

# These are the 4 polygon numbers that are lost, & causing mismatch with number of rows & number of polygons
valid_nums <- shp.lc.L$landcovers_middle_shapefile %>%
  makeValid() %>%
  pull(poly_num)
prob_polys <- setdiff(1:3443, valid_nums) # problematic polygons
# Note extent is NA, geometry is none, causing major issues.
shp.lc.L$landcovers_middle_shapefile[prob_polys, ] %>% head(n = 4)

# Can't plot these 4 polygons
ggplot() +
  geom_spatvector(data = shp.lc.L$landcovers_middle_shapefile[prob_polys, ])

# Generate list of invalid polygons
invalid <- map(shp.lc.L[2:4], \(polys){
  polys %>%
    mutate(is_valid = is.valid(.)) %>%
    filter(!is_valid)
})

# Make polygons valid & then plot
# It seems like the polygons still have weird stuff going on, may be best to have Mathilde go in & fix these issues?
valid <- map(invalid, \(polys_df){
  makeValid(polys_df)
})

# Plot invalid polygons, the # corresponds to the row numbers in the attribute table of QGIS
Invalid_poly_plots <- map2(invalid, valid, \(invalid_polys, valid_polys)
                           map(seq_len(nrow(valid_polys)), function(i) {
                             invalid_plot <- ggplot() +
                               geom_spatvector(data = invalid_polys[i, ], fill = "red", color = "black") +
                               labs( # title = list_name,
                                 subtitle = paste("Invalid Polygon", pull(invalid_polys[i, ], poly_num))
                               ) +
                               theme_min
                             valid_plot <- ggplot() +
                               geom_spatvector(data = valid_polys[i, ], fill = "red", color = "black") +
                               labs( # title = list_name,
                                 subtitle = paste("Valid Polygon", pull(valid_polys[i, ], poly_num))
                               ) +
                               theme_min
                             ggarrange(invalid_plot, valid_plot)
                           }))

# Print 3 PDFs (for each shapefile) with the invalid polygons in the left column, and the valid column (after using makeValid) on the right side
imap(Invalid_poly_plots, \(poly_plots, name){
  pdf(file = paste0("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Mentorship/Digitization_Mathilde/Feedback/invalid_polygons_", name, format(Sys.Date(), "%m.%d.%y"), ".pdf"), width = 8.5, height = 11, bg = "white")
  print(marrangeGrob(grobs = poly_plots, ncol = 3, nrow = 3, layout_matrix = matrix(1:3, 3, 1, TRUE)))
  dev.off()
})

# >Make polygons valid -----------------------------------------------------
# SHORT TERM SOLUTION -- Delete all invalid polygons
# LONG TERM SOLUTION -- Need to redraw

# Identify the problematic polygons
TF <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% is.valid()
})

# Remove all non-valid polygons to iron out the rest of the workflow
shp.lc.L[c(2:4)] <- map2(shp.lc.L[c(2:4)], TF, \(polys, TF){
  polys[TF, ]
})

# Make polygon files valid
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% makeValid()
})

# NOTE:: All polygons are valid at this point
map(shp.lc.L, \(shp){
  table(is.valid(shp))
})

# >Overlapping polygons ----------------------------------------------------
## Mathilde says, 'In one shapefile, polygons cannot overlap, it is not something that Q allows you to do'. However, some polygons DO overlap.
# Create a safe version of the `relate()` function with `possibly`
safe_relate <- possibly(
  ~relate(.x, relation = "overlaps", pairs = TRUE),
  otherwise = NULL,
  quiet = FALSE
)

# Map over the input with safe_relate function
overlap_matrix <- map(
  shp.lc.L[c(2:4)], safe_relate, .progress = TRUE
)

# Problematic coordinates given in the error message
problematic_coords <- matrix(c(
  -73.793410454502194, 3.6727849683307507 #,
  #-75.07352327711844, 10.800882487507739 # fixed
), ncol = 2, byrow = TRUE)
problematic_points <- terra::vect(problematic_coords, type = "points", crs = crs(shp.lc.L[[2]]))

# Find the polygon containing this point
problematic_polygons <- terra::relate(shp.lc.L[[2]], problematic_points,
                                      pairs = TRUE, relation = "intersects") %>%
  as_tibble()
problematic_polygons2 <- c(problematic_polygons$id.x, problematic_polygons$id.y) %>% unique()

# Manipulate problematic_polygons2 list if necessary? 
#problematic_polygons2 <- c(1087, 1100, 1102, 3089)
problematic_polygons2 <- setNames(problematic_polygons2, problematic_polygons2)

# Examine these polygons 
shp.lc.L[[2]] %>% filter(poly_num %in% problematic_polygons2) %>% 
  select(where(~ !all(is.na(.))))

# NOTE:: Several of these do have very small areas & would be removed if small areas were removed
shp.lc.L$landcovers_middle_shapefile_no_small_plg %>%
  filter(poly_num %in% problematic_polygons2) %>%
  pull(area)

# Remove the problematic polygons & run relate() again to see if it works
shp.lc.L[[2]] %>%
  filter(!poly_num %in% problematic_polygons2) %>%
  relate(relation = "overlaps", pairs = TRUE)

prob_poly_plots <- imap(problematic_polygons2, \(prob_poly, names){
  shp.lc.L$landcovers_middle_shapefile %>%
    filter(poly_num %in% prob_poly) %>%
    ggplot() +
    geom_spatvector() +
    labs(title = names) +
    theme_min
})
ggarrange(prob_poly_plots[[1]], prob_poly_plots[[2]])
#ggarrange(prob_poly_plots[[1]], prob_poly_plots[[2]], prob_poly_plots[[3]], prob_poly_plots[[4]])

# Plot problematic polygons
shp.lc.L$landcovers_middle_shapefile %>%
  filter(poly_num %in% problematic_polygons2) %>%
  ggplot() +
  geom_spatvector(aes(fill = factor(poly_num)), alpha = .5) +
  theme_min

# NOTE:: With these few polygons removed ALL relationships are valid!
relationships <- c("intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint")
map(relationships, \(relat){
  shp.lc.L[[2]] %>%
    filter(!poly_num %in% problematic_polygons2) %>%
    relate(relation = relat, pairs = TRUE)
}, .progress = TRUE)

# Sort each row and remove duplicates
overlap_indices <- map(overlap_matrix, \(om){
  unique(t(apply(om, 1, sort))) %>% as.data.frame()
})

# Generate overlap 'area' vector & plots for Mathilde
# For middle_UBC
area_mid_UBC <- as.numeric()
plots_mid_UBC <- list()
for (i in 1:nrow(overlap_indices$landcovers_middle_UBC_shapefile)) { #
  print(i)
  index <- overlap_indices$landcovers_middle_UBC_shapefile %>% slice(i)
  overlapping <- shp.lc.L$landcovers_middle_UBC_shapefile[c(index[, 1], index[, 2])]
  intersection <- terra::intersect(overlapping[1, ], overlapping[2, ])
  if (!is.empty(intersection)) {
    area_mid_UBC[i] <- round(expanse(intersection, unit = "m"), 2)
  } else {
    area_mid_UBC[i] <- 0
  }
  if (area_mid_UBC[i] > 50) {
    p <- ggplot() +
      geom_spatvector(
        data = overlapping,
        fill = c("red", "green"), alpha = .3
      ) +
      geom_spatvector(data = intersection, fill = "black") +
      labs(
        title = paste0("Polygons ", index[, 1], ", ", index[, 2]),
        subtitle = paste0("Overlap area = ", area_mid_UBC[i], " meters sq")
      ) +
      theme_min
    # Save plots to list
    plots_mid_UBC[[length(plots_mid_UBC) + 1]] <- p
  }
}
sort(area_mid_UBC)

# Confirm # of plots
length(plots_mid_UBC)

# Plots middle_UBC
ggarrange(plots_mid_UBC[[1]], plots_mid_UBC[[2]], plots_mid_UBC[[3]], plots_mid_UBC[[4]], ncol = 2, nrow = 2)

# For 'past' shapefile
area_past <- as.numeric()
plots_past <- list()
for (i in 1:nrow(overlap_indices$landcovers_past_shapefile)) { #
  print(i)
  index <- overlap_indices$landcovers_past_shapefile %>% slice(i)
  overlapping <- shp.lc.L$landcovers_past_shapefile[c(index[, 1], index[, 2])]
  intersection <- terra::intersect(overlapping[1, ], overlapping[2, ])
  if (!is.empty(intersection)) {
    area_past[i] <- round(expanse(intersection, unit = "m"), 2)
  } else {
    area_past[i] <- 0
  }
  if (area_past[i] > 50) {
    p <- ggplot() +
      geom_spatvector(
        data = overlapping,
        fill = c("red", "green"), alpha = .3
      ) +
      geom_spatvector(data = intersection, fill = "black") +
      labs(
        title = paste0("Polygons ", index[, 1], ", ", index[, 2]),
        subtitle = paste0("Overlap area = ", area_past[i], " meters sq")
      ) +
      theme_min
    # Save plots to list
    plots_past[[length(plots_past) + 1]] <- p
  }
}
sort(area_past)

# Confirm # of plots
length(plots_past)
# Plots past
ggarrange(plots_past[[1]], plots_past[[2]], plots_past[[3]], plots_past[[4]], ncol = 2, nrow = 2)

## Delete overlapping polygons
no_overlap <- map(shp.lc.L[c(3:4)], \(polys){
  erase(polys, sequential = TRUE)
})

# Compare areas
map_dbl(shp.lc.L[c(3:4)], \(polys){
  sum(expanse(polys))
})
map_dbl(no_overlap, \(polys){
  sum(expanse(polys))
})

# Determine if any polygons are completely isolated (not touching any other polygons)
touching <- map(shp.lc.L[c(3:4)], \(polys){
  relate(polys, relation = "intersects", pairs = TRUE) # intersects == if they touch, overlap, or are within
})

# Detect the polygons that are completely isolated
setdiff(1:475, unique(touching$landcovers_middle_UBC_shapefile[, 1]))
setdiff(1:634, unique(touching$landcovers_past_shapefile[, 1]))

# DELETE
# Try snapping to see if it fixes validity -- snap may be making invalid geometries that cannot be fixed with makeValid()
# snapped <- map(shp.lc.L[2:4], \(sf){
# terra::snap(sf, tolerance = 1)
# })

map(snapped, \(sf){
  table(is.valid(sf))
})
# Through here

# >Format -----------------------------------------------------------------
# Formatting of dataframe
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% mutate(
    image_date = as.Date(image_date),
    img_date_yr = chron::years(image_date),
    trees_perc = as.numeric(trees_perc)
  )
})

# Format linear features & individual trees
# Add LC
shp.lc.L[c(1, 5:7)] <- map2(
  shp.lc.L[c(1, 5:7)], c("distrees", rep("livefence", 3)),
  \(lf, lc_typ){
    lf %>% mutate(lc_typ = lc_typ)
  }
)

# Remove empty id column & set the trees percentage = 100
shp.lc.L$individual_trees_shapefile <- shp.lc.L$individual_trees_shapefile %>%
  select(-id) %>%
  mutate(trees_perc = 100)

# >Buffer -----------------------------------------------------------------
# Buffer trees & live fences (lf) by appropriate distances.
# Discuss buffers with Mathilde - "S" for shrub should be 2m radius (NOT diameter), "T" should be 5m wide, and "F" should be ~Xm radius, indiv trees should be radius 5m at low elevation & 3m at high elevation? See Mathilde's metadata file
# STILL NEED TO DIFFERENTIATE BUFFER SIZE BY TREE_SHRUB COLUMN
trees_lf_buff <- pmap(
  .l = list(shp.lc.L[c(1, 5:7)], buf.rad = c(4, 5, 5, 5)), # Buffer radii
  .f = \(shpL, buf.rad) buffer(shpL, buf.rad)
)

map(trees_lf_buff, \(shp){
  sum(expanse(shp))
})

# After buffering there may still be overlap between individual trees & live fences, but there shouldn't be. Cut out live fences using the individual trees polygons
trees_lf_buff[2:4] <- map2(trees_lf_buff[2:4], trees_lf_buff[1], \(lfs, trees){
  terra::erase(lfs, trees)
})

# >Combine polygons -------------------------------------------------------
## Cut out the buffered individual trees & live fences from the LC polygons
# Live fences
LCs_holes1 <- map2(shp.lc.L[2:4], trees_lf_buff[2:4], \(polys, lf){
  polys %>% terra::erase(lf)
})

# Individual trees
LCs_holes2 <- map2(LCs_holes1, trees_lf_buff[1], \(polys, trees){
  polys %>% terra::erase(trees)
})

# Combine layers & visualize an example
LCs_comb_l <- pmap(
  .l = list(LCs_holes2, trees_lf_buff[1], trees_lf_buff[2:4]),
  \(polys, trees, lf){
    rbind(polys, trees, lf)
  }
)

LCs_comb_l$landcovers_middle_shapefile %>%
  filter(lc_typ == "livefence" & !is.na(data_year)) %>%
  data.frame()


# CONFIRMED that LCs_comb_l is the sum of all of these SpatVectors
460 + 1255 + 1066
622 + 1255 + 1023

# CHECK:: Should be no NAs in lc_typ
map(LCs_comb_l, \(polys){
  polys %>%
    filter(is.na(lc_typ)) %>%
    nrow()
})

# Visualize an example buffer
ggplot() +
  geom_spatvector(data = LCs_comb_l[[1]], aes(fill = lc_typ), alpha = .4) +
  geom_spatvector(data = union_ms$b100, alpha = .3) + # Add buffers
  geom_sf(data = Pc_locs_sf, size = .5) + # Add centroids
  coord_sf(xlim = c(-73.683, -73.666), ylim = c(3.323, 3.343))

# ISSUE with middle_UBC file
# The terra::erase() function modifies the geometries, and if there are any invalid geometries or empty geometries after the erase operation, this could lead to inconsistencies in the number of rows when you later try to visualize the data.
ggplot() +
  geom_spatvector(data = LCs_comb_l[[2]], aes(fill = lc_typ), alpha = .4) +
  geom_spatvector(data = union_ms$b100, alpha = .3) + # Add buffers
  geom_sf(data = Pc_locs_sf, size = .5) + # Add centroids
  coord_sf(xlim = c(-73.683, -73.666), ylim = c(3.323, 3.343))

# Difference in number of polygons & number of lc types
nrow(shp.lc.L$landcovers_middle_UBC_shapefile) # 471 rows before erase function
nrow(LCs_holes1$landcovers_middle_UBC_shapefile) # 458 rows after
length(LCs_holes1$landcovers_middle_UBC_shapefile$lc_typ) # 437


# NOTE:: It seems like Mathilde left the linear features in the 'past' shapefile. Maybe this doesn't matter?
ggplot() +
  geom_spatvector(data = LCs_comb_l[[3]], aes(fill = lc_typ), alpha = .4) +
  geom_spatvector(data = union_ms$b100, alpha = .3) + # Add 50m radius buffers
  geom_sf(data = Pc_locs_sf, size = .5) + # Add centroids
  coord_sf(xlim = c(-73.683, -73.666), ylim = c(3.323, 3.343))


# >LC within 300m buffer ---------------------------------------
# See notes in notebook & at the top of script for ideas on how to do this
buffers_ms$b300

UBC_pts <- Pc_locs_sf %>% filter(Nombre_institucion == "UBC")

# Finca El Porvenir (both UniLlanos & I surveyed) have the data_year == 2019 but are not in the ‘middle_ubc’ file as well, or they don’t have data_year2 == 2022 (if there is no more recent imagery).
map(LCs_comb_l[c(1, 2)], \(polys){
  polys %>%
    filter(data_year == 2022) %>%
    ggplot() +
    geom_spatvector(aes(fill = lc_typ), alpha = .4) +
    geom_sf(data = UBC_pts, size = .5)
})

# >Id_group ----------------------------------------------------------------
## THIS IS AN OPTION TO CONSIDER, BUT IT MAY NOT BE NECESSARY. NEED TO TALK WITH MATHILDE & BETTER UNDERSTAND WHAT'S GOING ON. THINK BIG PICTURE, DON'T SPEND A BUNCH OF TIME ON THIS UNTIL WE'RE CONFIDENT IT IS THE BEST WAY FORWARD. Think about definition of 'intersects', is that what you want?

## Append the Id_group of each polygon to the spatvector data frame##
# relate() produces a matrix where column 1 is the Id # of the LC polygons & column 2 is the Id number of the buffer in union b300.
ints.df <- LCs_comb_l[[1]] %>%
  relate(union_ms$b300, "intersects", pairs = TRUE) %>%
  data.frame() %>%
  rename_with(~ c("id.lc", "id.buff")) %>%
  tibble()

# NOTE:: There are some digitized polygons that are contained in multiple buffers
overlaps4 <- ints.df %>%
  count(id.lc, sort = T) %>%
  filter(n == 4) %>%
  pull(id.lc)

## Let's visualize some example polygons + buffers that have 4 overlaps
test_lcs <- filter(LCs_comb_l[[1]], poly_num %in% as.character(overlaps4[1:20]))
extent <- ext(test_lcs) #+ c(-.0003, .0003, -.0000, .0000)

# Subset appropriate buffer IDs
ints.df %>% filter(id.lc == test_lcs$poly_num[1])
buff <- union_ms$b300[c(71, 75, 73, 78), ]


# Plot to better understand -- Andorra? Most of these polygons should only have 2 overlaps not 4 according to plot
ggplot() +
  geom_spatvector(
    data = test,
    aes(fill = lc_typ), alpha = .4
  ) +
  geom_spatvector(data = buff, aes(color = Id_group), alpha = .3) #+ #Add buffers
# geom_sf(data = Pc_locs_sf, size = .5) #+ #Add centroids
# coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = TRUE) #+
# ggrepel::geom_text_repel(data = union_ms$b300, aes(label = Id_group))


### NEXT STEP:: See if all polygons are contained within union_ms$b300
LCs_comb_l[[1]]

# Easy solution for now, just take the first Id_group, otherwise will need multiple columns for each Id_group.
# Kind of confusing b/c the 2020 file is formatted differently than the other files will be... Could also figure out which years were sampled when & select Id_group based on year somehow (think about this in more detail with Mathilde).
id.buffs <- ints.df %>%
  slice_head(by = id.lc) %>%
  pull(id.buff)
length(id.buffs)
nrow(LCs_comb_l[[1]])
LCs_comb_l[[1]]$Id_group <- union_ms$b300[id.buffs, ]$Id_group
unique(LCs_comb_l[[1]]$Id_group)

# Visualize the different Id_groups
# NOTE: The 3 of these overlap , examine with Mathilde
LCs_comb_l[[1]] %>%
  filter(Id_group %in% c("G-AD-M-CO", "G-MB-M-EA", "G-AD-M-LP")) %>%
  ggplot() +
  geom_spatvector(aes(fill = Id_group), alpha = .4)

# >Aggregate --------------------------------------------------------------
# NOTE:: Better to aggregate second because aggregate does not weigh the 'trees percent' column by the size of the polygon if you set fun = mean. Instead calculate the 'effective number of treed hectares' for each polygon & then sum them in the aggregate call.
# ALTERNATIVE:: Could probably just avoid aggregating at all for these calculations...
# To add 'crops_type' to aggregate call I think would need to paste lc & cropstype (e.g., crops_fruit). Doesn't really matter too much as there are very few crop polygons that are trees
LCs_comb %>%
  distinct(Id_group, lc, crops_type) %>%
  data.frame()
# The agg_tree_cover_ha is equal to the 'effective number of treed hectares'
LCs_agg_ID_lc <- aggregate(LCs_comb_a,
                           by = c("Id_group", "lc"), # "crops_type"
                           fun = function(x) {
                             sum(x, na.rm = TRUE)
                           }
)

# Visualize just the predominant landcovers to make legend easier to read
LCs_agg_ID_lc %>%
  filter(lc %in% c("intpast", "forest", "crops", "distrees", "livefence", "shrub")) %>%
  ggplot() +
  geom_spatvector(aes(fill = lc), alpha = .4) +
  coord_sf(xlim = c(-73.683, -73.666), ylim = c(3.323, 3.343))

# Inspect results
LCs_comb_a %>%
  data.frame() %>%
  dplyr::select(Id_group, lc, area, trees_prop, tree_cover_ha) %>%
  mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>%
  # filter(!is.na(crops_type)) %>%
  head(n = 30)

# LC areas by category
LCs_comb_a %>%
  as.data.frame() %>%
  group_by(lc) %>% # Id_group
  summarize(sum.area.ha = sum(area)) %>%
  arrange(desc(sum.area.ha))


# >Extract LC -------------------------------------------------------------
test <- terra::crop(LCs_comb_l[[3]], union_ms$b300)


# >Exploring options -------------------------------------------------------
# Calculate and add area attribute for each polygon
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], function(polys) {
  areas <- expanse(polys, unit = "m") # Calculate area in square meters
  polys$area <- round(areas, 5) # Add area as a new attribute
  return(polys)
})

# Filter out small polygons
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], function(polys) {
  polys %>% filter(area > 0.001)
})

# Apply a small buffer
small_buffer <- map(shp.lc.L[c(2:4)], function(polys) {
  polys %>% buffer(width = .00001)
})

# Try disaggregating
disagg_ubc <- shp.lc.L$landcovers_middle_UBC_shapefile %>% disagg()

# None of these potential solutions fix the issue. What else is there other than erase?
nrow(no_overlap$landcovers_middle_UBC_shapefile)
length(no_overlap$landcovers_middle_UBC_shapefile$lc_typ)
nrow(shp.lc.L$landcovers_middle_UBC_shapefile)

# Tried aggregating before using erase, but this didn't work, not sure why
# shp.lc.L[[1]] <- buffer(shp.lc.L[[1]], width = 4)

# shp.agg <- map(shp.lc.L, \(shp){
#  shp %>% aggregate(by = "lc_typ", dissolve = TRUE)
# })

# Sf approach -------------------------------------------------------------
# Original motivation is that st_distinct() might handle things differently than terra's erase().
# This is mostly the same code as above just using sf syntax

# ENSURE you are paying attention to sf_use_s2... Standard geometries are flat, while S2 geometries account for the Earth's curvature. Generally..
# Flat Geometries: Suitable for local studies, small regions, or when working with data that does not span large areas.
# S2 Geometries: Ideal for global datasets, spatial analyses over long distances, and complex geometry manipulations.
sf_use_s2(TRUE)

shp.lc.L <- map(.x = files.shp, \(shp){
  st_read(paste0(path, "/", shp)) %>%
    clean_names()
})
names(shp.lc.L) <- files.shp

# Add rownumber for easy identification & communication with Mathilde
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% mutate(poly_num = row_number())
})

# >Explore invalid polygons -----------------------------------------
# NOTE:: several polygons are not valid
map(shp.lc.L, \(polys){
  polys %>%
    st_is_valid() %>%
    table()
})

## We don't lose any polygons if working in sf
# Before makeValid
map(shp.lc.L[c(2:4)], \(polys){
  polys %>% nrow()
})

# But there are still a few polygons that are not actually valid
map(shp.lc.L[c(2:4)], \(polys){
  polys %>%
    st_is_valid() %>%
    table()
})

# Identify the problematic polygons
TF <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% st_is_valid()
})
# Visualize
map2(shp.lc.L[c(2:4)], TF, \(polys, TF){
  polys[!TF, ] %>% ggplot() +
    geom_sf()
})
# Subset just valid polygons for now to continue hammering out workflow
shp.lc.L[c(2:4)] <- map2(shp.lc.L[c(2:4)], TF, \(polys, TF){
  polys[TF, ]
})

# Generate list of invalid polygons
invalid <- map(shp.lc.L[2:4], \(polys){
  polys %>%
    mutate(is_valid = st_is_valid(.)) %>%
    filter(!is_valid)
})

# Make polygons valid & then plot
# It seems like the polygons still have weird stuff going on, may be best to have Mathilde go in & fix these issues?
valid <- map(invalid, \(polys_df){
  st_make_valid(polys_df)
})

# ERROR
# Plot invalid polygons, the # corresponds to the row numbers in the attribute table of QGIS
Invalid_poly_plots <- map2(invalid, valid, \(invalid_polys, valid_polys)
                           map(seq_len(nrow(invalid_polys)), function(i) {
                             invalid_plot <- ggplot() +
                               geom_sf(data = invalid_polys[i, ], fill = "red", color = "black") +
                               labs( # title = list_name,
                                 subtitle = paste("Invalid Polygon", pull(invalid_polys[i, ], poly_num))
                               ) +
                               theme_min
                             valid_plot <- ggplot() +
                               geom_sf(data = valid_polys[i, ], fill = "red", color = "black") +
                               labs( # title = list_name,
                                 subtitle = paste("Valid Polygon", pull(valid_polys[i, ], poly_num))
                               ) +
                               theme_min
                             ggarrange(invalid_plot, valid_plot)
                           }))

# Troubleshoot
for (i in c(1:25)) {
  p <- ggplot() +
    geom_sf(data = valid$landcovers_middle_UBC_shapefile[i, ], fill = "red", color = "black") +
    labs( # title = list_name,
      subtitle = paste("Valid Polygon", pull(valid$landcovers_middle_UBC_shapefile[i, ], poly_num))
    )
  print(p)
}

# Print 3 PDFs (for each shapefile) with the invalid polygons in the left column, and the valid column (after using makeValid) on the right side
imap(Invalid_poly_plots, \(poly_plots, name){
  pdf(file = paste0("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Mentorship/Digitization_Mathilde/Feedback/invalid_polygons_", name, format(Sys.Date(), "%m.%d.%y"), ".pdf"), width = 8.5, height = 11, bg = "white")
  print(marrangeGrob(grobs = poly_plots, ncol = 3, nrow = 3, layout_matrix = matrix(1:3, 3, 1, TRUE)))
  dev.off()
})

# >Make polygons valid -----------------------------------------------------
# SHORT TERM SOLUTION -- Delete these 4 corrupt polygons
# LONG TERM SOLUTION -- Need to redraw
shp.lc.L <- map(shp.lc.L, \(polys){
  polys %>% st_make_valid()
})


# NOTE:: Not all polygons are valid at this point
map(shp.lc.L, \(shp){
  table(st_is_valid(shp))
})

# Remove the single polygon that wasn't fixed
shp.lc.L[c(2:4)] <- map2(shp.lc.L[c(2:4)], TF, \(polys, TF){
  polys %>% filter(TF)
})

# >Buffer -----------------------------------------------------------
# Define the buffer radii for each type
buffer_sizes <- c(4, 5, 5, 5)

# Buffer trees and live fences by appropriate distances
trees_lf_buff <- pmap(
  .l = list(shp.lc.L[c(1, 5:7)], buffer_sizes),
  .f = \(shpL, buf_rad) st_buffer(shpL, dist = buf_rad)
)

# Define the helper function
st_erase <- function(x, y) {
  st_difference(x, st_union(st_combine(y)))
}

# Notice this produces error related to the trees! I don't know why the trees would be problematic, but after playing with different subsets of both sets of polygons, it is these trees that are causing problems.
trees_lf_buff[2:4] <- map2(trees_lf_buff[2:4], trees_lf_buff[1], \(lfs, trees){
  st_erase(lfs, trees)
})

# Join all polygons
union_trees <- st_union(trees_lf_buff[[1]])

# Notice there is no error when trees are unioned
map2(trees_lf_buff[2], union_trees, \(lfs, trees){
  lfs %>% st_erase(union_trees)
})

# Understanding the error message... I think this means
# Visualize the problematic trees
trees_lf_buff[[1]][c(261:262), ] %>%
  ggplot() +
  geom_sf()

trees_lf_buff[[1]][c(922:923), ] %>%
  ggplot() +
  geom_sf()


# >Combine polygons -------------------------------------------------------
## Cut out the buffered individual trees & live fences from the LC polygons
# Live fences
LCs_holes1 <- map2(shp.lc.L[2:4], trees_lf_buff[2:4], \(polys, lf){
  polys %>% st_erase(lf)
})

# ChatGPT
# Erase live fences from LC polygons
LCs_holes1 <- Map(function(lc, lf) {
  st_erase(lc, lf) # Erase all of lf from lc
}, sf_lc_layers, sf_lf_buffers)

# Erase individual trees from the resulting LC polygons
sf_trees_buffer <- st_as_sf(trees_lf_buff[[1]])
LCs_holes2 <- lapply(LCs_holes1, function(lc) {
  st_erase(lc, sf_trees_buffer) # Erase all of trees from each lc layer
})
