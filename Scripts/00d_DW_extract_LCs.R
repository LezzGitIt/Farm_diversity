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
  "terra", "tidyterra", "sf", "tidyverse", "janitor", "cowplot", "maptiles",
  "xlsx", "readxl", "gridExtra", "ggpubr", "conflicted", "landscapemetrics"
)

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Set basic themes, bring in functions, & data
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
    poly_num = row_number(),
    image_date = as.Date(image_date),
    img_date_yr = lubridate::years(image_date),
    trees_perc = as.numeric(trees_perc)
  )
})

# Export with poly number for Natalia
#shp.lc.L[[2]] %>% writeVector("../../Mentorship/Digitization_Mathilde/Final_docs/landcovers_middle_final_poly_num/landcovers_middle_final_poly_num.shp", overwrite = TRUE)

# Prob polys: identify and remove--------------------------------------------------
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
valid_polys <- shp.lc.L[[2]][Tf_l[[2]], ] #%>% mutate(poly_num = row_number())

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
valid_polys2 <- valid_polys #%>% filter(!poly_num %in% Prob_poly_id)

# >Buffer -----------------------------------------------------------------
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
# After buffering there may still be overlap between individual trees & live fences, but there shouldn't be. Cut out live fences using the individual trees polygons
trees_lf_buff[2:4] <- map2(t_lf_agg[2:4], t_lf_agg[1], \(lfs, trees){
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
miss_ids <- setNames(miss_ids, miss_ids)
Lc_holes2 <- Lc_holes1 %>% filter(!poly_num %in% miss_ids)

# NOTE:: Don't really understand why removing 1 polygon solves the issue with the difference in lengths for Lc_holes1 (difference of 4)
nrow(Lc_holes2)
length(Lc_holes2$poly_num)

# Erase the individual trees
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
    miss_poly <- valid_polys2 %>% filter(poly_num == id)
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


# LANDSCAPEMETRICS --------------------------------------------------------
# Rasterize --------------------------------------------------------------

# Rasterize ensuring landcover is the value of each cell
rast <- cropped_lcs %>%  #Lcs_sub
  rast(resolution = 1)

# Generate a list with each Id_muestreo 
Lcs_id_muestreo <- cropped_lcs %>% terra::split("Id_muestreo")

# Create extent with desired resolution
rast_l <- map(Lcs_id_muestreo, \(polys_group){
  polys_group %>% rast(resolution = 2)
})

# Rasterize polygons -- very quick 
Lc_rast_l <- map2(Lcs_id_muestreo, rast_l, \(polys_id, rast) {
  Lc_rast <- terra::rasterize(polys_id, rast, field = "lc_typ2") # Lcs_sub
})

# join_lc_class ---------------------------------------------------------------

# Extract 
uniq_classes <- get_unique_values(Lc_rast_l)

over4 <- keep(uniq_classes, ~ .x %>% length() > 4)
under4 <- keep(uniq_classes, ~ .x %>% length() < 4)
not4 <- keep(uniq_classes, ~ .x %>% length() != 4)
is4 <- keep(uniq_classes, ~ .x %>% length() == 4)

# Turn numeric for subsetting
is4_i <- str_split_i(names(is4), "_", 2) %>% as.numeric()
not4_i <- str_split_i(names(not4), "_", 2) %>% as.numeric()

## Create tbl to map lc_typ onto classes 
# Select any point count with all landcover types
head(is4_i) # Many to choose from

join_lc_class_all <- tibble(
  lc_typ2 = terra::unique(Lc_rast_l[[1]]$lc_typ2)[[1]],
  class = 0:3 # Adjust if needed
)

# There are several point counts that have < 4 landcover classes present.
join_lc_class <- map2(Lc_rast_l, uniq_classes, \(rast, class) {
  if(length(class) != 4){
    lc_typ2 <- rast$lc_typ2 %>% unique()
    tibble(lc_typ2) %>% 
      mutate(class = row_number()-1) %>% 
      right_join(tibble(class)) # This adds NA where length(class) == 5
  } else{
    join_lc_class_all
  }
}) %>% list_rbind(names_to = "Id_muestreo") %>% 
  # Replace NA with empty
  mutate(lc_typ2 = ifelse(is.na(lc_typ2), "empty", lc_typ2)) 

# Visualize rasters -------------------------------------------------------

## Visualize buffers if helpful (e.g. polygons with raster cells == NA)
# NOTE:: All Id_groups will show NAs due to the raster cells OUTSIDE of landcovers
ggplot() +
  geom_spatraster(data = Lc_rast_l$`C-MB-B-LM_04`) # 30, 37

# NOTE:: G-MB-M-EPO1 is error with coordinates that still needs to be fixed 
plots_u4 <- imap(Lc_rast_l[not4_i[1:3]], \(rast, group){
  ggplot() +
    geom_spatraster(data = rast) +
    labs(title = group)
})
#plots_u4

# Visualize with satellite imagery 
sat_bg <- get_tiles(Lc_rast_l[[1]], crop = TRUE, provider = "Esri.WorldImagery")
plot_tiles(sat_bg) 
plot(Lc_rast_l[[1]], add = TRUE, alpha = .5)

# Functions that could be helpful for exploration
#zoom(new = FALSE)
#draw(x = "points")
#click(id = TRUE)

# Calc metrics ------------------------------------------------------------

# Calculate lsm using Landscapemetrics package
Pc_cents <- Pc_vect_proj %>% filter(!Id_muestreo %in% paste0("OQ_0", 7:9)) %>%
  terra::split("Id_muestreo")

# Generate lsm for all items in these lists
start <- Sys.time()
Lsm_l <- map2(Lc_rast_l, Pc_cents, \(raster, points) {
  scale_sample(
    landscape = raster, 
    y = points, 
    what = c("lsm_c_pland", "lsm_c_te"), 
    shape = "circle", 
    size = Buffer_rad_nmr
  )
}, .progress = TRUE)
Sys.time() - start

# Rbind our lsm list, format, and join with join_lc_class for lc_typ
Lsm_df <- Lsm_l %>% list_rbind(names_to = "Id_muestreo") %>% 
  rename(buffer = size) %>% 
  select(-c(plot_id, id, layer)) %>% 
  left_join(join_lc_class)

# Should be no NAs if merge worked appropriately
Lsm_df %>% filter(is.na(lc_typ2))

# Inspect Lsm_df ----------------------------------------------------------

# Instead of recalculating (slow), load in lsm file
#Lsm_df <- read_xlsx("Derived/Excels/Lsm/Lsm_df_02.08.25.xlsx")

## NOTE:: There are 'percentage_inside' over 100 and under 90. Generate 1 row per Id_muestreo to handle more easily 
To_digitize <- Lsm_df %>% filter(percentage_inside < 100) %>% #> 100
  summarize(percent_inside = min(percentage_inside), .by = Id_muestreo) %>%
  mutate(sum_fun = "min") %>%
  arrange(percent_inside) #%>%
#filter(str_detect(Id_muestreo, "\\(1\\)|OQ")) 
To_digitize

over100 <- Lsm_df %>% filter(percentage_inside >= 104.7) %>% 
  summarize(percent_inside = max(percentage_inside), .by = Id_muestreo) %>%
  mutate(sum_fun = "max") %>% 
  arrange(desc(percent_inside))
over100 

# Create dataframe of ids that are problematic due to percentage inside (pi)
prob_pi_ids <- rbind(To_digitize, over104) %>% 
  mutate(across(where(is_double), ~ round(.x, 2)))

# Plot problematic pc buffers ------------------------------------------------

# Plot, using pmap to iterate over dataframe in rowwise fashion
prob_pi_plots <- pmap(prob_pi_ids[,1:2], \(Id_muestreo, percent_inside) {
  ggplot() +
    geom_spatraster(data = Lc_rast_l[[Id_muestreo]]) + 
    labs(subtitle = Id_muestreo, 
         caption = paste("Percent Inside:", percent_inside)) +
    theme_min + 
    theme(legend.position = "none")
})
# Visualize example
prob_pi_plots[[2]]

## Print PDF file with 9 plots per page
pdf("Figures/Prob_polys/Polygons_percent_inside.pdf")
print(marrangeGrob(grobs = prob_pi_plots, ncol = 3, nrow = 3, 
                   layout_matrix = matrix(1:9, 3, 3, TRUE)))
dev.off()

# Correlations ------------------------------------------------------------
# Not sure exactly what this is doing, may be easier to just use cor()
corr_l <- map(Lsm_l, \(lsm_tbl){
  lsm_tbl %>% filter(size == 300) %>%
    calculate_correlation(simplify = TRUE)
}) 
corr_df <- corr_l %>% list_rbind(names_to = "Id_muestreo") %>% 
  filter(metric_1 != metric_2) %>% 
  arrange(value)

# show_lsm() #Show landscape metrics on patch level printed in their corresponding patch.

# Save and export ---------------------------------------------------------
# Export centroids of point counts still to be digitized for Natalia
To_digitize %>% left_join(Pc_locs_sf) %>%
  slice_max(Buffer_rad) %>%
  filter(Id_muestreo != "OQ_Practica") #%>% 
#st_write("Derived_geospatial/shp/To_digitize_natalia/To_digitize.shp")

# Export Excel of lsm
Lsm_df %>% as.data.frame() %>%
  #if(FALSE){
  write.xlsx(
    file = paste0("Derived/Excels/Lsm/Lsm_df_", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), 
    showNA = FALSE, row.names = FALSE
  )
#}

# Extras -----------------------------------------------------------------
# >Working ----------------------------------------------------------------

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