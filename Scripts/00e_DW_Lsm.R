## PhD birds in silvopastoral landscapes##
## Data wrangling 00e -- Generate landscape metrics

## Inputs & Outputs
# Inputs: Snapped landcovers file from OOd script
# Output: Lsm_df files with landscapemetrics

# Instructions
# Go to 'Which shapefile?' section and set filename to middle, ubc, or past. This determines which shapefile is processed & the landscapemetrics that are ultimately exported 

# Script contents ---------------------------------------------------------
# 1) Load libraries and data
# 2) Set shapefile ('Which shapefile?')
# 3) Rasterize: landscapemetrics package requires rasters as input
# 4) Join landcover classes: We process the landcovers in a list based on the Id_muestreo (point count ID), & there may be differing numbers of landcovers in each point count. This step ensures that the class definitions (e.g., 0 = "forest) are equivalent for each point count
# 4) Visualize: any unusual rasters to help diagnose problems
# 5) Calc landscapemetrics: Use scale_sample() function to extract the landscapemetrics for each buffer
# 6) Inspect 'percentage inside': Examine point counts that have 'percentage inside' significantly < 100 or > 100
# 7) Plot: Plot point counts with low or high 'percentage inside' values
# 8) Pivot long to wide: Landscapemetrics defaults to a long dataframe, but we ultimately want a single row per Id_muestreo X buffer size. 
# 9) Date_year manual: Update Date_year column for "past" & "ubc" files
# 10) Export landscape metric dataframes: Export dataframes that are then joined in future script

# Load libraries --------------------------------------------------------
pkgs <- c(
  "terra", "tidyterra", "sf", "tidyverse", "janitor", "cowplot", "maptiles",
  "xlsx", "readxl", "gridExtra", "ggpubr", "conflicted", "landscapemetrics"
)

# Load packages
lapply(pkgs, library, character.only = TRUE)

# Set basic themes, conflicts preference, bring in handy functions
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")

# Load data
load("Rdata/the_basics_02.27.25.Rdata")

# Which shapefile?  --------------------------------------------------------
file_name <- "past" # middle, ubc, past

Snapped_lcs <- vect(paste0("Derived_geospatial/shp/R_processed/Snapped_", file_name, ".gpkg"))

# Still to do -------------------------------------------------------------
Metadata <- Snapped_lcs %>% data.frame()
# Many polygons have image date of 1970-01-01
Metadata %>% tabyl(image_date)
Metadata %>% filter(image_date == lubridate::as_date("1970-01-01")) %>% 
  tabyl(lc_typ) 

## CLEAN UP FILES 

# Correct a mistake ---------------------------------------------------------
if(file_name == "past"){
  Snapped_lcs <- Snapped_lcs %>% 
           mutate(lc_typ2 = if_else(lc_typ == "setosfor", "ssp", lc_typ2))
}

# Rasterize --------------------------------------------------------------
# Rasterize ensuring landcover is the value of each cell
rast <- Snapped_lcs %>%  #Lcs_sub
  rast(resolution = 1)

# Generate a list with each Id_muestreo 
Lcs_id_muestreo <- Snapped_lcs %>% terra::split("Id_muestreo")

# Create extent with desired resolution
rast_l <- map(Lcs_id_muestreo, \(polys_group){
  polys_group %>% rast(resolution = 1)
})

# Rasterize polygons -- very quick 
Lc_rast_l <- map2(Lcs_id_muestreo, rast_l, \(polys_id, rast) {
  Lc_rast <- terra::rasterize(polys_id, rast, field = "lc_typ2") # Lcs_sub
})

# Mask raster
#Lc_rast_crop <- map2(Lcs_id_muestreo, Lc_rast_l, \(polys_id, rast) {
 # terra::crop(rast, polys_id, mask = TRUE, touches = TRUE, snap = "out")
#})

# join_lc_class ---------------------------------------------------------------

# Extract 
uniq_classes <- get_unique_values(Lc_rast_l)

map(Lc_rast_l[c(15:17)], \(rast){
  ggplot() + geom_spatraster(data = rast)
})

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
  lc_typ2 = terra::unique(Lc_rast_l[[is4_i[1]]]$lc_typ2)[[1]],
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
# NOTE:: All Id_groups will show NAs due to the raster cells OUTSIDE of landcovers UBC−MB−M−LRO1_01

# NOTE:: In past file the 'setosfor' lc is causing lc > 4. Maybe fixing this will be enough
plots_u4 <- imap(Lc_rast_l[not4_i[1:4]], \(rast, group){
  ggplot() +
    geom_spatraster(data = rast) +
    labs(title = group)
})
#plots_u4

# Visualize with satellite imagery 
if(FALSE){
  sat_bg <- get_tiles(Lc_rast_l[[1]], crop = TRUE, provider = "Esri.WorldImagery")
  plot_tiles(sat_bg) 
  plot(Lc_rast_l[[1]], add = TRUE, alpha = .5)
}

# Functions that could be helpful for exploration
#zoom(new = FALSE)
#draw(x = "points")
#click(id = TRUE)

# Calc metrics ------------------------------------------------------------
## Calculate lsm using Landscapemetrics package

# Generate objects going into scale_sample() function
Pc_vect <- Pc_locs_sf %>% select(-c(Uniq_db, Nombre_institucion)) %>%
  distinct() %>%
  filter(Id_muestreo %in% names(Lc_rast_l)) %>%
  vect() %>% 
  project("EPSG:4686")

Pc_vect_proj <- Pc_vect %>% project("EPSG:3116")
Pc_cents <- Pc_vect_proj %>%
  terra::split("Id_muestreo")

Buffer_rad_nmr <- seq(from = 300, to = 50, by = -50)
Buffer_rad_nmr <- setNames(Buffer_rad_nmr, Buffer_rad_nmr)
# Generate lsm for all points / rasterized polys in these lists
start <- Sys.time()
safe_scale_sample <- safely(scale_sample)
Lsm_l_safe <- map2(Lc_rast_l, Pc_cents, \(raster, points) {
  safe_scale_sample(
    landscape = raster, 
    y = points, 
    what = c("lsm_c_pland", "lsm_l_te"), 
    shape = "circle", 
    size = Buffer_rad_nmr
  )
}, .progress = TRUE)
Sys.time() - start

# Extract results & errors, identify problematic points 
Lsm_l <- map(Lsm_l_safe, ~ .x$result) %>% 
  discard(is.null)
errors <- map(Lsm_l_safe, ~ .x$error)
prob_points <- names(discard(errors, is.null))
prob_points <- setNames(prob_points, prob_points)
prob_points 

# Rbind our lsm list, format, and join with join_lc_class for lc_typ
Lsm_df_long <- Lsm_l %>% list_rbind(names_to = "Id_muestreo") %>% 
  rename(buffer = size) %>% 
  select(-c(plot_id, id, layer)) %>% 
  left_join(join_lc_class)

# Should be no PLAND NAs if merge worked appropriately.
Lsm_df_long %>% filter(is.na(lc_typ2) & metric != "te")

# Percent_inside inspection------------------------------------------------------
## NOTE:: There are 'percentage_inside' (pi) over 100 and under 98. Generate 1 row per Id_muestreo to handle more easily 
Low_pi <- Lsm_df_long %>% filter(percentage_inside < 100) %>% #> 100
  summarize(percent_inside = min(percentage_inside), .by = Id_muestreo) %>%
  mutate(sum_fun = "min") %>%
  arrange(percent_inside)# %>%
  #filter(str_detect(Id_muestreo, "\\(1\\)|OQ")) 
Low_pi

High_pi <- Lsm_df_long %>% filter(percentage_inside >= 100) %>% 
  summarize(percent_inside = max(percentage_inside), .by = Id_muestreo) %>%
  mutate(sum_fun = "max") %>% 
  arrange(desc(percent_inside))
High_pi

# Create dataframe of ids that are problematic due to percentage inside
prob_pi_ids <- rbind(Low_pi, High_pi) %>% 
  mutate(across(where(is_double), ~ round(.x, 2)))

# Subset the 20 ids with the lowest & highest percent inside
prob_ids_min_max <- prob_pi_ids %>% slice(1:10, .by = sum_fun)

# >Plot problematic pi buffers ------------------------------------------------
## Plot, using pmap to iterate over dataframe in rowwise fashion
# Can use prob_pi_ids if that is better, but the highest percentage_inside values aren't very informative visually
prob_pi_plots <- pmap(prob_ids_min_max[,1:2], \(Id_muestreo, percent_inside) {
  ggplot() +
    geom_spatraster(data = Lc_rast_l[[Id_muestreo]]) + 
    labs(subtitle = Id_muestreo, 
         caption = paste("Percent Inside:", percent_inside)) +
    theme_min + 
    theme(legend.position = "none")
})
# Visualize example
# prob_pi_plots

## Print PDF file with 9 plots per page
if(FALSE){
  pdf(paste0("Figures/Prob_polys/Polygons_percent_inside_", file_name, ".pdf"))
  print(marrangeGrob(grobs = prob_pi_plots, ncol = 3, nrow = 3, 
                     layout_matrix = matrix(1:9, 3, 3, TRUE)))
  dev.off()
}

# Pivot long to wide ------------------------------------------------------
# Landscapemetrics package default returns long tibble, but we need a wide tibble for analyses 

# PLAND df
Lsm_pland <- Lsm_df_long %>% filter(metric == "pland") %>% 
  pivot_wider(id_cols = c(Id_muestreo, metric, buffer),
              names_from = lc_typ2,
              values_from = value, values_fill = 0) 

## CHECK:: All should equal 100
Check_add_tbl <- Lsm_pland %>% 
  mutate(Total = rowSums(across(4:7)))
table(near(Check_add_tbl$Total, 100))

# Total edge df
Lsm_te <- Lsm_df_long %>% filter(metric == "te") %>% 
  pivot_wider(id_cols = c(Id_muestreo, buffer),
              names_from = metric,
              values_from = value)

# Join tibbles - 484 rows
Lsm_df <- Lsm_pland %>% full_join(Lsm_te) %>% 
  select(-metric) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Date_year manual --------------------------------------------------------
## Make manual changes to data_year in ubc & past files. 
# The idea was to have the data_year column specify which point counts a set of polygons correspond to (& it could be multiple polygons, hence 'data_year2' & columns). Unfortunately, there were errors in Mathilde's metadata, so we forego correcting the 'middle'shapefile, & instead focus on the past & ubc files, which we use to overwrite the middle file in the correct locations.
if(file_name == "past"){
  Lsm_df <- Lsm_df %>% mutate(
    data_year = 2013, # General
    # "G-MB-G-CB" is an exception
    data_year = ifelse(str_detect(Id_muestreo, "G-MB-G-CB"), 2016, data_year)
  )
}

if(file_name == "ubc"){
  Lsm_df <- Lsm_df %>% mutate(data_year = 2022)
}

# Save and export ---------------------------------------------------------
stop()
# Export centroids of point counts still to be digitized for Natalia
Low_pi %>% left_join(Pc_locs_sf) %>%
  filter(Id_muestreo != "OQ_Practica") #%>% 
#st_write("Derived_geospatial/shp/To_digitize_natalia/To_digitize.shp")

# Export 
prob_pi_ids %>% as.data.frame() #%>%
  #write.xlsx(file = paste0("Derived/Excels/Lsm/Polys_percent_inside.xlsx"), 
  #showNA = FALSE, row.names = FALSE
#)

# Export Excel of lsm
Lsm_df_exp <- Lsm_df %>% mutate(lc_file = file_name) %>%
  as.data.frame() 
Lsm_df_exp %>%
  #if(FALSE){
    write.xlsx(
    file = paste0("Derived/Excels/Lsm/Lsm_df_", file_name, "_", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), showNA = FALSE, row.names = FALSE
  )
  #}

warnings() # 50+ warnings: Double values will be converted to integer.

# Extras ------------------------------------------------------------------
# >Correlations ------------------------------------------------------------
# Not sure exactly what this is doing, may be easier to just use cor()
corr_l <- map(Lsm_l, \(lsm_tbl){
  lsm_tbl %>% filter(size == 300) %>%
    calculate_correlation(simplify = TRUE)
}) 
corr_df <- corr_l %>% list_rbind(names_to = "Id_muestreo") %>% 
  filter(metric_1 != metric_2) %>% 
  arrange(value)

# show_lsm() #Show landscape metrics on patch level printed in their corresponding patch.

# >Troubleshooting ---------------------------------------------------------
# NOTE:: data_year is nan for livefence in past & middle_ubc file
tabyl(Snapped_lcs$data_year)
Snapped_lcs %>% filter(data_year == 2016)
Snapped_lcs %>% filter(is.nan(data_year)) %>% 
  pull(lc_typ)

Id_prob <- c("G-MB-A-LRE", "G-MB-M-LP1", "G-MB-G-CB") #"G-MB-A-LRE", G-MB-M-LP1, "G-MB-G-CB"
Id_prob <- setNames(Id_prob, Id_prob)
#"G-MB-G-CB" change to 2013
imap(Id_prob, \(id, name){
  Snapped_lcs %>% 
    filter(data_year %in% c(2016, 2017, 2019) & Id_group == id) %>% 
    ggplot() + 
    geom_spatvector() +
    labs(title = name)
})

prob <- "G-MB-G-CB"

Snapped_lcs %>% filter(Id_muestreo == "G-MB-M-A_01-B") %>%
  ggplot() + 
  geom_spatvector() +
  labs(title = prob)

# In middle: "G-MB-A-LRE" = 2016, G-MB-G-CB = 2017, "G-MB-M-LP1" = 2017
# G-MB-G-CB SHOULD ACTUALLY BE 2016 IN PAST FILE!? 
Bird_pcs %>% filter(Id_group %in% Id_prob) %>% 
  distinct(Id_group, Ano) %>% 
  arrange(Id_group)


## PART 2 -- TROUBLESHOOT
Prob_ids <- Lsm_pland %>% filter(empty > 0) %>% 
  pull(Id_muestreo) %>% 
  unique() 
Prob_ids <- setNames(Prob_ids, Prob_ids) # Does G-MB-A-CL_01 still show up? 

df <- Snapped_lcs %>% filter(Id_muestreo == Prob_ids[1])
df2 <- Snapped_lcs %>% 
  filter(Id_muestreo == Prob_ids[1] & lc_typ2 == "other")
ggplot() +
  geom_spatvector(data = df, aes(fill = lc_typ2), alpha = .1) +
  geom_spatvector(data = df2, aes(fill = lc_typ2), color = "red", alpha = 1) +
  labs(title = Prob_ids[1])

Join_data_year <- Snapped_lcs %>% as_tibble() %>%
  select(Id_muestreo, Uniq_db, starts_with("data_year")) %>% 
  distinct() %>% 
  filter(!if_all(starts_with("data_year"), ~ is.nan(.x)))

# Some manual cleaning
Join_data_year %>% filter(data_year = case_when(
  Uniq_db == "Ubc mbd" ~ 2022
))

# Identify points with multiple rows
Ids_mult_rows <- Join_data_year %>% mutate(n = n(), .by = Id_muestreo) %>% 
  filter(n > 1) %>% 
  arrange(desc(n)) %>% view()
pull(Id_muestreo)
Join_data_year %>% filter(Id_muestreo == "G-MB-M-EA_01")

## Investigating ##
Snapped_lcs %>% as_tibble() %>% 
  filter(Id_group %in% c("U-MB-M-EPO3", "UBC-MB-M-EPO3", "G-AD-M-EPO3")) %>% 
  select(Id_group, starts_with("data_year"), image_date) %>% 
  distinct()

Join_data_year %>% 
  filter(Uniq_db %in% c("Unillanos mbd", "Ubc mbd") & Id_muestreo == "U-MB-M-EPO3_04") %>%
  filter(Id_muestreo %in% Ids_mult_rows) %>% 
  tabyl(data_year)

## 

EPO3 <- c("U-MB-M-EPO3", "UBC-MB-M-EPO3", "G-AD-M-EPO3")
CO_EA <- c("G-MB-M-EA", "G-AD-M-CO")

Snapped_lcs %>% data.frame() %>%
  filter(Id_group %in% CO_EA) %>% 
  distinct(Id_group, data_year, data_year2, image_date)

Prob_poly <- Snapped_lcs %>% filter(Id_muestreo == "G-MB-M-EA_02") %>% 
  filter(data_year == 2019)
Prob_poly <- Snapped_lcs %>% filter(Id_muestreo == "G-AD-M-CO_02") %>% 
  filter(data_year == 2017 & is.na(data_year2))
Snapped_lcs %>% filter(Id_group %in% CO_EA) %>%
  #filter(data_year == 2019) %>% #c(2016, 2017) 
  ggplot() + 
  geom_spatvector(aes(fill = Id_group), alpha = .5) + 
  geom_spatvector(data = Prob_poly, fill = "black", alpha = 1, size = 5)

Bird_pcs %>% filter(Id_group %in% c("G-AD-M-CO", "G-MB-M-EA")) %>% 
  distinct(Id_muestreo, Ano) %>% view()

