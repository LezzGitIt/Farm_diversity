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

# Load libraries & data --------------------------------------------------------
library(terra)
library(tidyterra)
library(sf)
library(tidyverse)
library(janitor)
library(cowplot)
library(xlsx)
library(readxl)
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
  #geom_spatvector(data = cropped_lcs, aes(fill = lc_typ2), alpha = .8) +
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
Pc_vect <- Pc_locs_sf %>% vect() %>% 
  project("EPSG:4686")
Pc_vect_proj <- Pc_vect %>% project("EPSG:3116")

Buffer_rad_nmr <- seq(from = 300, to = 50, by = -50)
Buffer_rad_nmr <- setNames(Buffer_rad_nmr, Buffer_rad_nmr)
Buffers <- map(Buffer_rad_nmr, \(rad){
  Pc_vect %>% buffer(rad) %>% 
    select(Uniq_db, Ecoregion, Departamento, Id_group, Id_muestreo)
})

# Intersect each buffer with the landcover object
cropped_lcs <- Lcs_comb2 %>% terra::intersect(Buffers$`300`) %>%
  project("EPSG:3116") %>%
  mutate(poly_num = row_number())

# landscapemetrics --------------------------------------------------------
# >Rasterize --------------------------------------------------------------
# Rasterize ensuring landcover is the value of each cell
# I think using the 300m buffers is better, and maybe one at a time if needed for speed
rast <- cropped_lcs %>%  #Lcs_sub
  rast(resolution  = 2)

# Rasterize
terraOptions(threads = 4) #About 42 mins
if(FALSE){
  start <- Sys.time()
  Lc_rast <- terra::rasterize(cropped_lcs, rast, field = "lc_typ2") # Lcs_sub
  Sys.time() - start
}
# Ultimately don't need to save the 1.5 GB large landcover raster once lsm are finalized
Lc_rast <- rast("Derived_geospatial/tif/Landcovers_mathilde_2m.tif")

# >Calc metrics ------------------------------------------------------------
# OPTION 1 - Use built in lsm function 
# Provides convenient % inside column
if(FALSE){
  # About 10 minutes
  Lsm_df <- map(Buffer_rad_nmr, \(rad){
    sample_lsm(Lc_rast, y = Pc_vect_proj, 
               shape = "circle", size = rad, 
               plot_id = Pc_vect_proj$Id_muestreo,
               what = c("lsm_c_pland", "lsm_c_te"),
               #level = "landscape", type = "diversity metric",
               return_raster = TRUE)
  }) %>% list_rbind(names_to = "Buffer_rad") %>% 
    mutate(Buffer_rad = as.numeric(Buffer_rad))

  ## Create tbl to map lc_typ onto classes 
  # Take a single plot with all classes present
  Sing_plot <- Lsm_df %>% 
    filter(plot_id == "G-MB-Q-PORT_01" & Buffer_rad == 300) %>%
    pull(raster_sample_plots) %>% 
    .[[1]]
  map_class <- tibble(
    lc_typ = terra::unique(Sing_plot$lc_typ2)[[1]],
    class = 0:3 # Adjust if needed
  )
  
  # Merge with map_class
  Lsm_df_rast <- Lsm_df %>% left_join(map_class) %>% 
    rename(Id_muestreo = plot_id) %>% 
    mutate(across(where(is.numeric), \(col) round(col, 2))) %>%
    select(-c(id)) # What is id?
}

# Instead of recalculating (slow), load in lsm file
Lsm_df2 <- read_xlsx("Derived/Excels/Lsm_df_01.20.25.xlsx")

# # NOTE:: There are 'percentage_inside' over 100 and under 90?! 
To_digitize <- Lsm_df2 %>% filter(percentage_inside < 90) %>% #> 100
  arrange(desc(percentage_inside)) %>% 
  distinct(Id_muestreo, Buffer_rad, percentage_inside) %>% 
  #filter(str_detect(Id_muestreo, "\\(1\\)|OQ")) %>% 
  arrange(Id_muestreo, Buffer_rad)
To_digitize

# Save and export ---------------------------------------------------------
# Export files to digitize for Natalia
To_digitize %>% left_join(Pc_locs_sf) %>%
  slice_max(Buffer_rad) %>%
  filter(Id_muestreo != "OQ_Practica") #%>% 
#st_write("Derived_geospatial/shp/To_digitize/To_digitize.shp")

# Export Excel of lsm
Lsm_df_rast %>% select(-raster_sample_plots) %>%
  as.data.frame() %>%
if(FALSE){
  write.xlsx(
    file = paste0("Derived/Excels/Lsm_df_", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), 
                  showNA = FALSE, row.names = FALSE
    )
}

## Export landcover shapefiles for Diana Laura 
cropped_lcs %>% 
  select(Ecoregion, Departamento, lc_typ, trees_perc, image_date, alt_source) #%>%
  #terra::writeVector("Derived_geospatial/Diana_laura_digitization/Training_polygons/Training_polygons.shp")

## Export rasters of buffer sizes
# Create dataframe to contain plots
Lsm_rast_buffs <- Lsm_df_rast %>%
  group_by(Id_muestreo, Buffer_rad) %>%
  slice(1) %>% # Take the first row for each group
  ungroup() %>%
  rename(Plots = raster_sample_plots) %>% 
  select(Id_muestreo, Buffer_rad, Plots, percentage_inside)

# Export rasters: 1 folder per Id_muestreo, 6 buffer sizes per folder 
path <- "Derived_geospatial/tif/Id_buffers_tifs/"
Lsm_rast_export <- Lsm_rast_buffs %>% 
  mutate(
    dir_path = paste0(path, Id_muestreo), # Create directory path
    filename = paste0(dir_path, "/Landcover_", Id_muestreo, "_buff", Buffer_rad, ".tif") 
  ) %>% arrange(Id_muestreo, Buffer_rad)
Lsm_rast_export %>% rowwise() %>% 
  group_split() %>% 
  purrr::walk(function(row) {
    if (!dir.exists(row$dir_path)) {
      dir.create(row$dir_path, recursive = TRUE, showWarnings = FALSE)
      writeRaster(row$Plots[[1]], filename = row$filename, overwrite = FALSE)
  }})

# Import rasters ------------------------------------------------------------
## If you want to do any visualization or additional calculation of metrics, can bring rasters back in
# Create a nested list with ID at top level and buffer size as next level
path <- "Derived_geospatial/tif/Id_buffers_tifs/"
Ids_files <- list.files(path = path, full.names = FALSE)
Lc_rast_l <- map(Ids_files, function(id) {
  # Get the files for this ID
  buffer_files <- list.files(
    path = file.path(path, id), 
    pattern = "\\.tif$", 
    full.names = TRUE
  )
  
  # Extract buffer sizes from filenames and sort them
  buffer_sizes <- str_extract(buffer_files, "(?<=_buff)\\d+(?=\\.tif$)") %>% as.numeric()
  sorted_files <- buffer_files[order(buffer_sizes)]
  
  rasters <- map(sorted_files, ~ terra::rast(.x))
  
  # Name each raster by its buffer size extracted from the filename
  names(rasters) <- str_extract(buffer_files, "(?<=_buff)\\d+(?=\\.tif$)")
  return(rasters)
})
names(Lc_rast_l) <- Ids_files

# Ensure that landscapes pass the check 
#check_landscape(Lc_rast_l)

# Visualize some examples
map(Lc_rast_l$`C-MB-A-ED_03`, \(sr){ #sr = spatraster
    plot(sr, main = dim(sr))
  }) 

# Extras -----------------------------------------------------------------
# >Working ----------------------------------------------------------------

# >Interactive tmap -------------------------------------------------------
library(tmap)
library(tmaptools)

# Ensure that shp_test and trees_lf_buff[[2]] are SpatVectors
# Convert to sf for better compatibility with tmap (optional)
shp_test_sf <- sf::st_as_sf(shp_test) %>% st_make_valid()
trees_sf <- sf::st_as_sf(trees_lf_buff[[2]])

TF <- shp_test_sf %>% 
  st_is_valid()
shp_test_sf <- shp_test_sf[TF,]

# Create an interactive map
tmap_mode("view")  # Switch to interactive mode

shp_test_sf %>% filter()
tm_shape() +
  tm_polygons(col = "blue", alpha = 0.5, border.alpha = 0.7) +
  #tm_shape(trees_sf) +
  #tm_borders(col = "red", lwd = 2, alpha = 0.8) +
  tm_layout(title = "Interactive Map of Problematic Areas",
            legend.outside = TRUE)

# >Invalid polygons  -------------------------------------------------------
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