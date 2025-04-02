## PhD birds in silvopastoral landscapes##
## Data wrangling 00d -- Extract landcover associated with each point 
## This script extracts relevant landcover and landscape metrics for each point count location

## Inputs & Outputs
# Inputs: Mathilde / Natalia's digitized polygons & linear features (3 layers representing different temporal periods), a single layer of dispersed trees (no temporal component)
# Output: Snapped landcovers file

# Instructions
# Go to 'Which shapefile?' section and set tbl_row <- # 1 = middle, 2 = past, 3 = ubc. This determines which file is processed & ultimately exported 

# Script contents ---------------------------------------------------------
# 1) Load libraries and data
# 2) Subset individual trees: Relevant for past & ubc files
# 3) Buffer: Buffering of trees & live fences
# 4) Erase: Cut out the buffered trees & live fences from the landcover polygons
# 5) Combine: Combine the holey landcover with the buffered trees + live fences
# 6) Extract LC: Use buffers to extract the landcover polygons associated with each point count
# 7) Snap: Fill in any small gaps in the LC polygons
# Extras:: Several things associated with problem solving 

## TO DO: See little journal
# Something weird going on with image_date -- trace backwards to understand where date errors are coming from

# Load libraries & Rdata --------------------------------------------------
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
conflicted::conflicts_prefer(terra::intersect)
source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")
load("Rdata/the_basics_02.27.25.Rdata")

# Which shapefile?  -------------------------------------------------------
# This determines which file is processed & ultimately exported 
tbl_row <- 2 # 1 = middle, 2 = past, 3 = ubc

# Create Index_tbl
Uniq_db <- unique(Bird_pcs$Uniq_db)
Index_tbl <- tibble(name = c("middle", "past", "ubc"), index = c(2,3,4), data_collector = list(Uniq_db[c(1:3,6)], c("Gaica mbd", "Cipav mbd"), "Ubc mbd"))

# Bring in data -----------------------------------------------------------
#Mathilde shapefiles
path <- "../../Mentorship/Digitization_Mathilde/Final_docs/final_shp_for_natalia"
files.shp <- list.files(path = path, pattern = "final") #, recursive = TRUE)
#files.shp <- files.shp[-2] # Remove old lc_middle file

shp.lc.L <- map(.x = files.shp, \(shp)
                vect(paste0(path, "/", shp)) %>%
                  clean_names())
files.shp2 <- str_remove(files.shp, ".gpkg")
names(shp.lc.L) <- files.shp2

# Formatting LC dataframes
shp.lc.L[c(2:4)] <- map(shp.lc.L[c(2:4)], \(polys){
  polys_df <- polys %>% mutate(
    #poly_num = row_number(),
    image_date = as.Date(image_date),
    img_date_yr = lubridate::years(image_date),
    trees_perc = as.numeric(trees_perc)
  )
})

# Subset individual trees -------------------------------------------------
## The individual trees file was not differentiated by the data collection year, given that individual trees left in pasture did not change much over the project. 
# Thus, our objective is to subset the individual trees to only the relevant IDs if we are working with UBC or past shapefiles
if((tbl_row + 1) %in% c(3,4)){
  shp.lc.L[[1]] <- intersect(shp.lc.L[[tbl_row + 1]], shp.lc.L[[1]])
}

# Visualize 
ggplot() + geom_spatvector(data = shp.lc.L[[tbl_row + 1]]) + 
  geom_spatvector(data = shp.lc.L[[1]]) + # Individual trees
  geom_spatvector(data = shp.lc.L[[tbl_row + 1 + 3]]) # Live fences

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
Lc_holes1 <- erase(shp.lc.L[[tbl_row + 1]], trees_lf_buff[[tbl_row + 1]])

## Second erase (individual trees) to create LC_holes2
Lc_holes2 <- erase(Lc_holes1, trees_lf_buff[[1]])

# Combine landcovers ------------------------------------------------------
# Combine the holey polygon landcover, the live fences, and the individual trees
Lcs_comb <- rbind(Lc_holes2, trees_lf_buff[[tbl_row + 1]], trees_lf_buff[[1]])
nrow(Lcs_comb)
length(Lcs_comb$lc_typ)

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
Pc_vect <- Pc_locs_sf %>% 
  filter(Uniq_db %in% Index_tbl[[tbl_row, "data_collector"]][[1]]) %>%
  vect() %>% 
  project("EPSG:4686") 
Pc_vect_proj <- Pc_vect %>% project("EPSG:3116")

# ALTERNATIVE APPROACH NEW DOESN'T WORK
if(FALSE){
  Pc_vect <- Pc_locs_sf %>% 
    distinct(Ecoregion, Departamento, Id_group, Id_muestreo, geometry) %>%
    vect() %>% 
    project("EPSG:4686")
}

## Cycle through years? 
if(FALSE){
  Pc_date8 %>% distinct(Ecoregion, Id_muestreo, Ano) %>% 
    pivot_wider(id.cols = Id_muestreo, 
                names_prefix = "Ano", names_glue = 1:3,
                values_from = Ano) 
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

# Export ------------------------------------------------------------------
## Export cropped landcover shapefile - The polygons that mathilde digitized that have been processed in this R script

Snapped_lcs %>% 
  terra::writeVector(paste0("Derived_geospatial/shp/R_processed/Snapped_", Index_tbl[[tbl_row, "name"]], ".gpkg"), overwrite = TRUE) #, overwrite = TRUE
stop()

# Extras -----------------------------------------------------------------
# >Zoom in og LC file -------------------------------------------
# Use the extent of a problematic set of polygons to examine the original shapefile
bbox <- terra::ext(Intersect_lcs$`G-MB-M-A_01-B`)  # Extract extent

shp.3116 <- shp.lc.L[[3]] %>% project("EPSG:3116")
ggplot() +
  geom_spatvector(data = shp.3116) +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4])) +
  ggtitle("G-MB-M-A_01-B")

# >Rm miss_ids ------------------------------------------------------------
# CHECK lengths:: Has any metadata been deleted with erase()? If so, important to remove the missing IDs
nrow(Lc_holes1)
length(Lc_holes1$poly_num)

# Remove the polygons that are lost during erase so metadata aligns
miss_ids <- setdiff(shp.lc.L[[tbl_row + 1]]$poly_num, Lc_holes1$poly_num) 
miss_ids <- setNames(miss_ids, miss_ids)
Lc_holes2 <- Lc_holes1 %>% filter(!poly_num %in% miss_ids)

# NOTE:: Don't really understand why removing 1 polygon solves the issue with the difference in lengths for Lc_holes1 (difference of 4)
nrow(Lc_holes2)
length(Lc_holes2$poly_num)

## There are new missing polygons created in this erase! 
# Examine and remove
miss_ids3 <- setdiff(shp.lc.L[[tbl_row + 1]]$poly_num, Lc_holes3$poly_num)
miss_ids3 <- setNames(miss_ids3, miss_ids3)
Lc_holes4 <- Lc_holes3 %>% filter(!poly_num %in% miss_ids3)

## Visualize polygons that are lost in the series of erases
#if(FALSE){
miss_id_plots <- imap(miss_ids3, \(id, name){
  miss_poly <- shp.lc.L[[tbl_row + 1]] %>% filter(poly_num == id)
  ggplot() + 
    geom_spatvector(data = miss_poly) + 
    labs(title = name) + 
    theme_min
})
if(length(miss_id_plots) < 15){
  miss_id_plots
}
#}

# >Troubleshoot -----------------------------------------------------------------
# Identify problematic polys & Id groups -- change depending on which groups are
dissappear_polys_id <- miss_ids3
dissappear_polys <- Snapped_lcs %>% filter(poly_num %in% dissappear_polys_id)
Prob_id_groups <- unique(dissappear_polys$Id_group)
#Prob_id_groups <- c("G-AD-M-EPO3", "U-MB-M-EPO3")
#Prob_id_muestreo <- unique(dissappear_polys$Id_muestreo)

# Select problematic buffers
Prob_groups <- Buffers$`300` %>% filter(Id_group %in% Prob_id_groups )

# Visualize problem
Prob_cropped <- Lcs_comb2 %>% terra::intersect(Prob_groups)
ggplot() + geom_spatvector(data = Prob_cropped) 

# I confirmed the problem begins with Lc_holes1 object (after first erase), but don't remember how I determined this.
Lc_holes1 

## Intersect prob buffer with vp2
# NOTE:: It is not the act of intersect that is removing the polygons in question
Cropped_lcs_vp2 <- shp.lc.L[[tbl_row + 1]] %>% terra::intersect(Prob_groups) %>%
  project("EPSG:3116")

#Prob_muestreo <- Cropped_lcs_vp2 %>% filter(Id_muestreo %in% Prob_id_muestreo)

# Visualize
ggplot() + geom_spatvector(data = Cropped_lcs_vp2, alpha = .7) + 
  geom_spatvector(data = Prob_muestreo, aes(fill = lc_typ), alpha = .6)

## Highlight polygons that are disappearing w/ erase() in blue. There are livefences likely causing at least some of the problems.
miss_polys_vp2 <- shp.lc.L[[tbl_row + 1]] %>% filter(poly_num %in% dissappear_polys_id)
miss_polys_vp3 <- miss_polys_vp2 %>% filter(poly_num %in% c(471, 456))
# Visualize
ggplot() + 
  geom_spatvector(data = miss_polys_vp3, fill = "pink", alpha = .7) +
  geom_spatvector(data = miss_polys_vp3, 
                  color = "blue", alpha = 1, linewidth = 3) +
  geom_spatvector(data = Prob_cropped, aes(fill = lc_typ2), alpha = .3) 


# >Overlapping polys ------------------------------------------------------
# NOTE:: Did the overlapping polygons piece with test object (created next line) & found that there is no overlap
#test <- list(landcovers_middle_final = Intersect_lcs) #

# Create a safe version of the `relate()` function with `possibly`
safe_relate <- possibly(
  ~terra::relate(.x, relation = "overlaps", pairs = TRUE),
  otherwise = NULL,
  quiet = FALSE
)

# Map over the input with safe_relate function
overlap_matrix <- purrr::map(
  shp.lc.L[tbl_row + 1], safe_relate, .progress = TRUE #shp.lc.L[2]
)

# Sort each row and remove duplicates
overlap_indices <- map(overlap_matrix, \(om){
  unique(t(apply(om, 1, sort))) %>% as.data.frame()
})

## Generate overlap 'area' vector & plots for Natalia
# For 'middle_final' shapefile
area_middle <- as.numeric()
plots_middle <- list()
for (i in 1:nrow(overlap_indices[[1]])) { #
  print(i)
  index_loop <- overlap_indices[[1]] %>% slice(i)
  overlapping <- shp.lc.L[[tbl_row + 1]][c(index_loop[, 1], index_loop[, 2])]
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
        title = paste0("Polygons ", index_loop[, 1], ", ", index_loop[, 2]),
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
ggarrange(plots_middle[[1]], plots_middle[[2]], plots_middle[[3]], plots_middle[[4]], plots_middle[[5]], plots_middle[[6]], ncol = 2, nrow = 3)

# >Problematic polygons ---------------------------------------------------
# Identify problematic polygons, including polygons: That were invalid even after st_make_valid, polygon causing a topology exception at specific coords, polygons that are lost during erase 

# >>Invalid ---------------------------------------------------------------
# NOTE:: There are 45 polygons that are invalid by sf standards. These don't seem to influence the workflow at present but important to keep in mind if things start to go wrong. 
# See below for plotting of these polygons

# Convert the polygons to sf
Lcs_comb_sf <- shp.lc.L[[tbl_row + 1]] %>% st_as_sf()

# Identify problematic polygons
TF <- Lcs_comb_sf %>% st_is_valid()
problematic <- Lcs_comb_sf[!TF,]
corrected <- problematic %>% st_make_valid()

#NOTE:: Difference in terra & sf
table(TF)
shp.lc.L[[tbl_row + 1]] %>% is.valid() %>% table() # makeValid() %>%

# NOTE:: these are different areas, indicating that holes are being filled in
#Lcs_comb2 %>% filter(poly_num == 2904) %>% expanse()
#corrected  %>% filter(poly_num == 2904) %>% st_area()

# Create individual plots
plot_list <- map2(1:nrow(problematic[1:76,]), 1:nrow(corrected[1:76,]), ~ {
  prob_plot <- ggplot() +
    geom_sf(data = problematic[.x, ], fill = "red", alpha = 0.5) +
    labs(title = paste("Problematic Polygon", as.character(problematic$poly_num[.x]))) +
    theme_minimal()
  
  corr_plot <- ggplot() +
    geom_sf(data = corrected[.y, ], fill = "blue", alpha = 0.5) +
    labs(title = paste("Corrected Polygon", as.character(corrected$poly_num[.y]))) +
    theme_minimal()
  
  ggarrange(prob_plot, corr_plot, ncol = 2)  # Arrange side by side
})

# Save all plots in a single PDF file
#pdf("Figures/Prob_polys/Invalid_Corrected_vs_Problematic_unaltered.pdf", width = 12, height = 8)
walk(plot_list, print)
#dev.off()

# >>Topology exception-------------------------------------------------
## There were several polygons that were causing problems executing the script. Identify and remove for now until Mathilde / Natalia can fix 

## At present, the relate function works! (no Topology exception), thus the rest of this code is unnecessary 
terra::relate(shp.lc.L[[tbl_row + 1]], pairs = FALSE, relation = "overlaps")

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
valid_polys <- shp.lc.L[[tbl_row + 1]][Tf_l[[tbl_row + 1]], ] #%>% mutate(poly_num = row_number())

# Identify the polygons that overlap with the given points 
prob_coords <- matrix(c(-73.793538544036281, 3.6730780784730706), ncol = 2, byrow = TRUE)
prob_point <- terra::vect(prob_coords, type = "points", crs = crs(shp.lc.L[[index]]))
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