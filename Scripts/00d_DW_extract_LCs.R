## PhD birds in silvopastoral landscapes##
## Data wrangling 00d -- Extract landcover associated with each point 
## This script extracts relevant landcover and landscape metrics for each point count location

## NOTE:: Anything relying on Lsm_df_rast (which still contains the 'raster_sample_plots' column) will not export. If changes need to be made to exported files, will need to rerun the sample_lsm() function

# Contents
# 1) Load libraries and data
# 2) Identify problematic polygons 
#   a) 1 polygon that is invalid even after st_make_valid
#   b) Another polygon causing a topology exception at specific coords
#   c) (Later) 5 polygons that are lost during erase 
# 3) Format: including buffering of trees & live fences
# 4) Erase: Cut out the buffered trees & live fences from the landcover polygons
# 5) Combine: Combine the holey landcover with the buffered trees + live fences
# 6) Extract LC: Use buffers to extract the landcover polygons associated with each 
# 7) landscapemetrics: Rasterize to then use landscapemetrics package

## TO DO: See little journal
# Something weird going on with image_date -- trace backwards to understand where date errors are coming from

# Load libraries & data --------------------------------------------------------
pkgs <- c(
  "terra", "tidyterra", "sf", "tidyverse", "cowplot", "maptiles",
  "xlsx", "readxl", "gridExtra", "ggpubr", "conflicted", "janitor"
)

# Install any necessary packages 
#lapply(pkgs, install.packages)

# Load packages
lapply(pkgs, library, character.only = TRUE)

# Set basic themes, bring in functions, & data
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")
load("Rdata/the_basics_02.27.25.Rdata")

# Bring in data -----------------------------------------------------------
#Mathilde shapefiles
path <- "../../Mentorship/Digitization_Mathilde/Final_docs/final_shp_for_natalia"
files.shp <- list.files(path = path, pattern = "final") #, recursive = TRUE)
#files.shp <- files.shp[-2] # Remove old lc_middle file

shp.lc.L <- map(.x = files.shp, \(shp)
                vect(paste0(path, "/", shp)) %>%
                  clean_names())
names(shp.lc.L) <- files.shp

# Formatting LC dataframes
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% mutate(
    #poly_num = row_number(),
    image_date = as.Date(image_date),
    img_date_yr = lubridate::years(image_date),
    trees_perc = as.numeric(trees_perc)
  )
})

# Export with poly number for Natalia
#shp.lc.L[[2]] %>% writeVector("../../Mentorship/Digitization_Mathilde/Final_docs/landcovers_middle_final_poly_num/landcovers_middle_final_poly_num.shp", overwrite = TRUE)

# Buffer -----------------------------------------------------------------
## Format linear features & individual trees
# Add LC
shp.lc.L[c(1, 5:7)] <- map2(
  shp.lc.L[c(1, 5:7)], c("distrees", rep("livefence", 3)),
  \(lc, lc_typ){
    lc %>% mutate(lc_typ = lc_typ)
  }
)

## Buffer trees & live fences (lf) by appropriate distances.
# Discuss buffers with Mathilde - "S" for shrub should be 2m radius (NOT diameter), "T" should be 5m wide, and "F" should be ~Xm radius, indiv trees should be radius 5m at low elevation & 3m at high elevation? See Mathilde's metadata file
# STILL NEED TO DIFFERENTIATE BUFFER SIZE BY TREE_SHRUB COLUMN
trees_lf_buff <- pmap(
  .l = list(shp.lc.L[c(1, 5:7)], buf.rad = c(4, 5, 5, 5)), # Buffer radii
  .f = \(shpL, buf.rad) buffer(shpL, buf.rad) %>% 
    select(-contains("tree_shrub"))
)

# Merge polygons so there are no weird overlaps within a shapefile
t_lf_agg <- map(trees_lf_buff, \(t_lf){
  aggregate(t_lf, by = "lc_typ", dissolve = TRUE) %>% #  by = names(t_lf)
    select(-agg_n)
})

# Remove empty id column & set the trees percentage = 100
t_lf_agg[[1]] <- t_lf_agg[[1]] %>%
  select(-mean_id) %>%
  mutate(trees_perc = 100)

# Erase -------------------------------------------------------------------
## First erase (live fences) to create LC_holes1
# After buffering there may still be overlap between individual trees & live fences, but there shouldn't be. Cut out live fences using the individual trees polygons
trees_lf_buff[2:4] <- map2(t_lf_agg[2:4], t_lf_agg[1], \(lfs, trees){
  terra::erase(lfs, trees)
})

## Cut out the buffered individual trees & live fences from the LC polygons
# Live fences
Lc_holes1 <- erase(shp.lc.L[[2]], trees_lf_buff[[2]])

# CHECK lengths:: Has any metadata been deleted with erase()? If so, important to remove the missing IDs
nrow(Lc_holes1)
length(Lc_holes1$poly_num)

# >Rm miss_ids ------------------------------------------------------------
# Remove the polygons that are lost during erase so metadata aligns
miss_ids <- setdiff(shp.lc.L[[2]]$poly_num, Lc_holes1$poly_num) 
miss_ids <- setNames(miss_ids, miss_ids)
Lc_holes2 <- Lc_holes1 %>% filter(!poly_num %in% miss_ids)

# NOTE:: Don't really understand why removing 1 polygon solves the issue with the difference in lengths for Lc_holes1 (difference of 4)
nrow(Lc_holes2)
length(Lc_holes2$poly_num)

# Visualize polygons that are lost in the series of erases
#if(FALSE){
miss_id_plots <- imap(miss_ids, \(id, name){
  miss_poly <- shp.lc.L[[2]] %>% filter(poly_num == id)
  ggplot() + 
    geom_spatvector(data = miss_poly) + 
    labs(title = name) + 
    theme_min
})
if(length(miss_id_plots) < 15){
  miss_id_plots
}
#}
#--- End >Rm miss_ids

## Second erase (individual trees) to create LC_holes2
Lc_holes3 <- erase(Lc_holes2, trees_lf_buff[[1]])

# Combine landcovers ------------------------------------------------------
# Combine the holey polygon landcover, the live fences, and the individual trees
Lcs_comb <- rbind(Lc_holes3, trees_lf_buff[[2]], trees_lf_buff[[1]])
unique(Lcs_comb$lc_typ)

Lcs_comb2 <- Lcs_comb %>% mutate(lc_typ2 = case_when(
  lc_typ %in% c("bareground", "crop", "crops", "water", "shrub", "devland") ~ "other",
  lc_typ %in% c("distrees", "forbanks", "silvopast", "livefence") ~ "ssp",
  .default = lc_typ
))

## Visualize an example buffer
if(FALSE){
  # Buffers from Landcover_digitization.R script
  path <- "Derived_geospatial/shp/"
  Buff_50m <- vect(paste0(path, "Buffers_50m/Buffers_50m.shp"))
  #Buff_300m <- vect(paste0(path, "Buffers_300m/Buffers_300m.shp"))
  ggplot() +
    geom_spatvector(data = Lcs_comb2, aes(fill = lc_typ2), alpha = .4) +
    #geom_spatvector(data = cropped_lcs, aes(fill = lc_typ2), alpha = .8) +
    geom_spatvector(data = Buff_50m, alpha = .3) +
    geom_sf(data = Pc_locs_sf, size = .5) + 
    coord_sf(xlim = c(-73.683, -73.666), ylim = c(3.323, 3.343))
}

# Extract LC  -------------------------------------------------------------
Pc_vect <- Pc_locs_sf %>% vect() %>% 
  project("EPSG:4686")
Pc_vect_proj <- Pc_vect %>% project("EPSG:3116")

# ALTERNATIVE APPROACH NEW DOESN'T WORK
if(FALSE){
  Pc_vect <- Pc_locs_sf %>% 
    distinct(Ecoregion, Departamento, Id_group, Id_muestreo, geometry) %>%
    vect() %>% 
    project("EPSG:4686")
}

Buffer_rad_nmr <- seq(from = 300, to = 50, by = -50)
Buffer_rad_nmr <- setNames(Buffer_rad_nmr, Buffer_rad_nmr)
Buffers <- map(Buffer_rad_nmr, \(rad){
  Pc_vect %>% buffer(rad) %>% 
    select(Uniq_db, Ecoregion, Departamento, Id_group, Id_muestreo)
})

##  Intersect each buffer with the landcover object
# NEW APPROACH (FAIL)
#if(FALSE){
Buff_300_l <- Buffers$`300` %>% terra::split("Id_muestreo")
Intersect_lcs <- map(Buff_300_l, \(buff_id){
  Lcs_comb2 %>% terra::intersect(buff_id) %>% 
    project("EPSG:3116")
})
#}

# OLD APPROACH 
if(FALSE){
  Intersect_lcs <- Lcs_comb2 %>% terra::intersect(Buffers$`300`) %>% 
    project("EPSG:3116") %>%
    mutate(poly_num = row_number())
}

# Perf circle - DELETE
ggplot() + geom_spatvector(data = Intersect_lcs$`G-AD-M-LC_01`)

# Snap vertices -----------------------------------------------------------
# There are small gaps in the data, which will elevate some landscape metrics (e.g. total edge). Fill in these gaps. 
Id_groups <- unique(Intersect_lcs$Id_group)
start <- Sys.time()
# NOTE:: tidyterra::filter() is not working in map function? Weird
Snapped_lcs_l <- map(Intersect_lcs, \(Lc_id) {
  Lc_id %>% snap(tolerance = 5)
})
Sys.time() - start 

# Return back to a single SpatVector
Snapped_lcs <- Snapped_lcs_l %>% vect()

# Troubleshoot -----------------------------------------------------------------
# Identify problematic polys & Id groups
dissappear_polys <- c(1559, 1569)
Prob_id_groups <- c("G-AD-M-EPO3", "U-MB-M-EPO3")

# Select problematic buffers
Id_epo3 <- Buffers$`300` %>% filter(Id_group %in% Prob_id_groups )

# Visualize problem
Id_epo3_cropped <- Lcs_comb2 %>% terra::intersect(Id_epo3)
ggplot() + geom_spatvector(data = Id_epo3_cropped) 

# I confirmed the problem begins with Lc_holes1 object (after first erase), but don't remember how I determined this.
Lc_holes1 

## Intersect epo3 buffer with vp2
# NOTE:: It is not the act of intersect that is removing the polygons in question
Cropped_lcs_vp2 <- shp.lc.L[[2]] %>% terra::intersect(Id_epo3) %>%
  project("EPSG:3116")

EPO3_07 <- Cropped_lcs_vp2 %>% filter(Id_muestreo == "G-AD-M-EPO3_07")

ggplot() + geom_spatvector(data = Cropped_lcs_vp2, alpha = .7) + 
  geom_spatvector(data = EPO3_07, aes(fill = lc_typ), alpha = .6)

## Highlight polygons that are disappearing w/ erase in blue. There are livefences likely causing at least some of the problems. Please redraw these livefences as polygons (use imagery to estimate the appropriate width & shape)
miss_polys_vp2 <- shp.lc.L[[2]] %>% filter(poly_num %in% dissappear_polys)
ggplot() + 
  geom_spatvector(data = miss_polys_vp2, fill = "pink", alpha = .7) +
  geom_spatvector(data = miss_polys_vp2, 
                  color = "blue", alpha = 1, linewidth = 3) +
  geom_spatvector(data = Id_epo3_cropped, aes(fill = lc_typ2), alpha = .3) 

# Export ------------------------------------------------------------------
## Export cropped landcover shapefile - The polygons that mathilde digitized that have been processed in this R script
Snapped_lcs %>% 
  terra::writeVector("Derived_geospatial/shp/Snapped_lcs.gpkg") #overwrite = TRUE

# Extras -----------------------------------------------------------------
# >Overlapping polys ------------------------------------------------------
# Create a safe version of the `relate()` function with `possibly`
safe_relate <- possibly(
  ~relate(.x, relation = "overlaps", pairs = TRUE),
  otherwise = NULL,
  quiet = FALSE
)

# Map over the input with safe_relate function
overlap_matrix <- map(
  shp.lc.L[[2]], safe_relate, .progress = TRUE
)

# Sort each row and remove duplicates
overlap_indices <- map(overlap_matrix, \(om){
  unique(t(apply(om, 1, sort))) %>% as.data.frame()
})

## Generate overlap 'area' vector & plots for Natalia
# For 'middle_final' shapefile
area_middle <- as.numeric()
plots_middle <- list()
for (i in 1:nrow(overlap_indices$landcovers_middle_final)) { #
  print(i)
  index <- overlap_indices$landcovers_middle_final %>% slice(i)
  overlapping <- shp.lc.L$landcovers_middle_final[c(index[, 1], index[, 2])]
  intersection <- terra::intersect(overlapping[1, ], overlapping[2, ])
  if (!is.empty(intersection)) {
    area_middle[i] <- round(expanse(intersection, unit = "m"), 2)
  } else {
    area_middle[i] <- 0
  }
  if (area_middle[i] > 10) {
    p <- ggplot() +
      geom_spatvector(
        data = overlapping,
        fill = c("red", "green"), alpha = .3
      ) +
      geom_spatvector(data = intersection, fill = "black") +
      labs(
        title = paste0("Polygons ", index[, 1], ", ", index[, 2]),
        subtitle = paste0("Overlap area = ", area_middle[i], " meters sq")
      ) +
      theme_min
    # Save plots to list
    plots_middle[[length(plots_middle) + 1]] <- p
  }
}
sort(area_middle)

# Confirm # of plots
length(plots_middle) # Only plots > 150 meters squared
# Plots middle
ggarrange(plots_middle[[1]], plots_middle[[2]], plots_middle[[3]], plots_middle[[4]], plots_middle[[5]], ncol = 2, nrow = 3)

# >Problematic polygons ---------------------------------------------------
# >>Invalid ---------------------------------------------------------------
# NOTE:: There are 45 polygons that are invalid by sf standards. These don't seem to influence the workflow at present but important to keep in mind if things start to go wrong. 
# See below for plotting of these polygons
if(FALSE){
  Lcs_comb_sf <- Lcs_comb2 %>% st_as_sf()
  TF <- Lcs_comb_sf %>% st_is_valid()
  corrected <- Lcs_comb_sf[!TF,] %>% st_as_sf() %>%
    st_make_valid()
  corrected
}

# NOTE:: these are different areas, indicating that holes are being filled in
Lcs_comb2 %>% filter(poly_num == 2904) %>% expanse()
corrected  %>% filter(poly_num == 2904) %>% st_area()

# Create a list of individual polygons
polygon_list <- map(1:nrow(corrected), ~ corrected [.x, ])

# Generate individual plots
map(polygon_list, ~ ggplot() +
      geom_spatvector(data = .x, fill = "blue", alpha = 0.5) +
      labs(title = paste("Polygon", as.character(.x$poly_num))) +
      theme_minimal())
##

# >>Topology exception-------------------------------------------------
## There were several polygons that were causing problems executing the script. Identify and remove for now until Mathilde / Natalia can fix 

## At present, the relate function works! (no Topology exception), thus the rest of this code is unnecessary 
terra::relate(shp.lc.L[[2]], pairs = FALSE, relation = "overlaps")

## Remove invalid polygons identified in sf
# For whatever reason sf is more strict than terra -- some polygons that are valid in terra, are not valid in sf
map(shp.lc.L, ~ table(is.valid(.x))) 
Tf_l <- map(shp.lc.L, \(shp){
  sf_shp <- sf::st_as_sf(shp) 
  sf_shp %>% st_is_valid() %>% 
    table()
  sf_shp2 <- sf_shp %>% st_make_valid()
  TF <- sf_shp2 %>% st_is_valid()
  return(TF)
}) 
valid_polys <- shp.lc.L[[2]][Tf_l[[2]], ] #%>% mutate(poly_num = row_number())

# Identify the polygons that overlap with the given points 
prob_coords <- matrix(c(-73.793538544036281, 3.6730780784730706), ncol = 2, byrow = TRUE)
prob_point <- terra::vect(prob_coords, type = "points", crs = crs(shp.lc.L[[2]]))
TF <- terra::relate(valid_polys, prob_point, 
                    pairs = FALSE, relation = "intersects") 
table(TF[,1])
Prob_poly_id <- which(TF[,1])
Prob_poly <- valid_polys %>% filter(poly_num %in% Prob_poly_id)

# Use this polygon to bring in additional surrounding polygons:
TF <- terra::relate(valid_polys, Prob_poly, 
                    pairs = FALSE, relation = "intersects") 
table(TF[,1])
Prob_poly_ids <- which(TF[,1])
Prob_polys <- valid_polys %>% filter(poly_num %in% Prob_poly_ids)

# Visualize the problematic point and the surrounding polygons 
ggplot() +
  geom_spatvector(data = Prob_polys) +
  geom_spatvector(data = Prob_poly, fill = "red") +
  geom_spatvector(data = prob_point, color = "blue", size = 4, alpha = 0.5) + 
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +  # Round to 1 decimal place
  labs(x = "Longitude",
       y = "Latitude") +
  theme_minimal() 

## REMOVE: Prob_polygon
valid_polys2 <- valid_polys #%>% filter(!poly_num %in% Prob_poly_id)

# >Subsetting -------------------------------------------------------------
# Large rasters, can be helpful to subset to more manageable size
join_group <- Pc_locs_sf %>% distinct(Id_muestreo, Id_group)
Id_buff <- expand_grid(Id_muestreo = Pc_locs_sf$Id_muestreo, 
                       Buffer_rad = names(Buffers)) %>% 
  left_join(join_group)

# Subset
cropped_sub <- cropped_lcs %>% filter(Id_group == "UBC-MB-M-A") 
subset <- unique(Id_buff$Id_group)[1:3]
Lcs_sub <- cropped_lcs %>% filter(Id_group %in% subset)
Id_buff_sub <- Id_buff %>% filter(Id_group %in% subset)
Pc_sub <- Pc_vect_proj %>% filter(Id_group %in% subset)


# >Option 2 ---------------------------------------------------------------
## DELETE 
# Instead of using built in lsm function.. I don't think this is necessary , but keep around for a while in case
Id_muestreo <- unique(Id_buff_sub$Id_muestreo)
Id_muestreo <- setNames(Id_muestreo, Id_muestreo)

Buffers_proj <- map(Buffers, ~.x %>% project("EPSG:3116"))
Buffer_rad_chr <- setNames(as.character(Buffer_rad_nmr), Buffer_rad_nmr)
mask_rast <- map(Id_muestreo, \(id){
  map(Buffer_rad_chr, \(rad){
    Buff_id <- Buffers_proj[[rad]] %>% filter(Id_muestreo == id)
    crop_rast <- crop(rast, Buff_id, mask = TRUE)
    return(crop_rast)
  })
})

pland <- map_depth(mask_rast, 2, ~lsm_c_pland(.x))
lsm_df3 <- list_flatten(pland, name_spec = "{outer}_buff{inner}") %>% 
  list_rbind(names_to = "Id_buff") %>% 
  mutate(Id_muestreo = str_split_i(Id_buff, "_buff", 1), 
         Buff_rad = str_split_i(Id_buff, "_buff", 2)) %>% 
  select(-Id_buff)
lsm_df3
## 

# OLD, no longer need to split into a list. DELETE?
if(FALSE){
  # Split into a list so all polygon cells are split up by Id_muestreo. 
  cropped_sub_l <- cropped_lcs %>% filter(Id_group == "UBC-MB-M-A") %>% 
    terra::split("Id_muestreo")
  names(cropped_sub_l) <- map_chr(cropped_sub_l, \(crop_l){
    crop_l %>% pull(Id_muestreo) %>% 
      unique()
  })
  
  # Plot 
  imap(cropped_sub_l, \(cropped, names){
    cropped %>%
      ggplot() +
      geom_spatvector(aes(fill = lc_typ2)) +
      labs(title = names)
  })
}