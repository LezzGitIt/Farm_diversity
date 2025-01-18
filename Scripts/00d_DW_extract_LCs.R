## PhD birds in silvopastoral landscapes##
## Data wrangling 00d -- Extract landcover associated with each point 
## This script generates the outputs related to functional traits (Avo_traits_final), & elevational ranges (Elev_ranges) that will be used in future scripts 

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

# Load libraries & data --------------------------------------------------------
library(terra)
library(tidyterra)
library(sf)
library(tidyverse)
library(janitor)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(conflicted)
library(landscapemetrics)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")
load("Rdata/the_basics_01.09.25.Rdata")

# Bring in data -----------------------------------------------------------
#Mathilde shapefiles
path <- "../../Mentorship/Digitization_Mathilde/Final_docs/final_shp_for_natalia"
files.shp <- list.files(path = path, pattern = "final") #, recursive = TRUE)
#files.shp <- files.shp[-2] # Remove old lc_middle file

shp.lc.L <- map(.x = files.shp, \(shp)
                vect(paste0(path, "/", shp)) %>%
                  clean_names())
names(shp.lc.L) <- files.shp

# Buffers from Landcover_digitization.R script
path <- "Derived_geospatial/shp/"
Buff_50m <- vect(paste0(path, "Buffers_50m/Buffers_50m.shp"))
#Buff_300m <- vect(paste0(path, "Buffers_300m/Buffers_300m.shp"))

# Formatting LC dataframes
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], \(polys){
  polys %>% mutate(
    image_date = as.Date(image_date),
    img_date_yr = lubridate::years(image_date),
    trees_perc = as.numeric(trees_perc)
  )
})

# Prob polys: identify and remove ---------------------------------------------------------
## There are several polygons that are causing problems executing the script. Identify and remove for now until Mathilde / Natalia can fix 
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
valid_polys <- shp.lc.L[[2]][Tf_l[[2]], ] %>% mutate(poly_num = row_number())

# Visualize the polygon that is still invalid after st_make_valid
invalid <- shp.lc.L[[2]] %>% filter(poly_num %in% which(!Tf_l[[2]])) 
ggplot() + 
  geom_spatvector(data = invalid) + 
  labs(title = "poly_num 3177")

## Relate still causes Topology exception
#terra::relate(valid_polys, pairs = FALSE, relation = "overlaps")

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
valid_polys2 <- valid_polys %>% filter(!poly_num %in% Prob_poly_id)

# Formatting --------------------------------------------------------------
## Format linear features & individual trees
# Add LC
shp.lc.L[c(1, 5:7)] <- map2(
  shp.lc.L[c(1, 5:7)], c("distrees", rep("livefence", 3)),
  \(lf, lc_typ){
    lf %>% mutate(lc_typ = lc_typ)
  }
)

# Remove empty id column & set the trees percentage = 100
shp.lc.L[[1]] <- shp.lc.L[[1]] %>%
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


# Erase -------------------------------------------------------------------
# After buffering there may still be overlap between individual trees & live fences, but there shouldn't be. Cut out live fences using the individual trees polygons
trees_lf_buff[2:4] <- map2(trees_lf_buff[2:4], trees_lf_buff[1], \(lfs, trees){
  terra::erase(lfs, trees)
})

## Cut out the buffered individual trees & live fences from the LC polygons
# Live fences
Lc_holes1 <- erase(valid_polys2, trees_lf_buff[[2]])

# NOTE:: Different lengths, certain metadata is being deleted with erase()
nrow(Lc_holes1)
length(Lc_holes1$poly_num)

# Remove the 5 polygons that are lost during erase
miss_ids <- setdiff(valid_polys2$poly_num, Lc_holes1$poly_num) 
Lc_holes2 <- Lc_holes1 %>% filter(!poly_num %in% miss_ids)

# Erase the individual trees
Lc_holes3 <- erase(Lc_holes2, trees_lf_buff[[1]])

# Combine landcovers ------------------------------------------------------
# Combine the holey polygon landcover, the live fences, and the individual trees
Lcs_comb <- rbind(Lc_holes3, trees_lf_buff[[2]], trees_lf_buff[[1]])
unique(Lcs_comb$lc_typ)

Lcs_comb2 <- Lcs_comb %>% mutate(lc_typ2 = case_when(
  lc_typ %in% c("bareground", "crops", "water", "shrub", "devland") ~ "other",
  lc_typ %in% c("distrees", "forbanks", "silvopast", "livefence") ~ "ssp",
  .default = lc_typ
))

# Visualize an example buffer
ggplot() +
  geom_spatvector(data = Lcs_comb2, aes(fill = lc_typ2), alpha = .4) +
  #geom_spatvector(data = cropped_buffers, aes(fill = lc_typ2), alpha = .8) +
  geom_spatvector(data = Buff_50m, alpha = .3) +
  geom_sf(data = Pc_locs_sf, size = .5) + 
  coord_sf(xlim = c(-73.683, -73.666), ylim = c(3.323, 3.343))

# Visualize polygons that are lost in the 
if(FALSE){
  miss_id_plots <- imap(miss_ids, \(id, name){
    miss_poly <- valid_polys %>% filter(poly_num == id)
    ggplot() + 
      geom_spatvector(data = miss_poly) + 
      labs(title = name) + 
      theme_min
  })
  ggarrange(miss_id_plots[[1]], miss_id_plots[[2]], miss_id_plots[[3]], miss_id_plots[[4]], miss_id_plots[[5]])
}

# Extract LC  -------------------------------------------------------------
Buff_100m <- Pc_locs_sf %>% vect() %>% 
  buffer(100) %>% 
  project("epsg:4686") %>% 
  select(Uniq_db, Ecoregion, Departamento, Id_group, Id_muestreo)

# NOTE:: There are 45 polygons that are invalid by sf standards. These don't seem to influence the workflow at present but important to keep in mind if things start to go wrong. 
# See below for plotting of these polygons
if(FALSE){
  Lcs_comb_sf <- Lcs_comb2 %>% st_as_sf()
  TF <- Lcs_comb_sf %>% st_is_valid()
  corrected <- Lcs_comb_sf[!TF,] %>% st_as_sf() %>%
    st_make_valid()
  corrected
}

# Intersect each buffer with the landcover object
cropped_buffers <- terra::intersect(Lcs_comb2, Buff_100m) %>%  
  #aggregate(by = c("lc_typ", "Id_muestreo"))
  project("EPSG:3116") %>% 
  mutate(poly_num = row_number())

# Split into a list so all polygon cells are split up by Id_muestreo. 
cropped_sub_l <- cropped_buffers %>% filter(Id_group == "UBC-MB-M-A") %>% 
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

# landscapemetrics --------------------------------------------------------
# >Rasterize --------------------------------------------------------------
# Rasterize ensuring landcover is the value of each cell
rast_l <- map(cropped_sub_l, \(cropped){
  rast <- cropped %>% rast(resolution = 1)
  terra::rasterize(cropped, rast, field = "lc_typ2")
})
plot(rast_l[[1]])
# Ensure that landscapes pass the check 
map(rast_l, check_landscape)

# >Calc metrics ------------------------------------------------------------
# Change class to say lc_typ, and include the id (name of the list)  # Look up from EWPW code
landscapemetrics::lsm_c_pland(rast_l)

# Working -----------------------------------------------------------------
## Invalid polygons from above 
# Note these are different areas, indicating that holes are being filled in
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

##NOTE:: May be able to keep as a single raster and later intersect with buffers 
cropped_sub <- cropped_buffers %>% filter(Id_group == "UBC-MB-M-A")
rast <- cropped_sub %>% rast(resolution = 1)
rast <- terra::rasterize(cropped_sub, rast, field = "lc_typ2")

terra::intersect(rast, cropped_sub)
##
