## PhD birds in silvopastoral landscapes##
## Data wrangling 00e -- Generate landscape metrics

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
Snapped_lcs <- vect("Derived_geospatial/shp/Snapped_lcs.gpkg")

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

# Or any specific problems 
# G-AD-M-LC_01 has several bumps sticking out, several internal holes as well
ggplot() + geom_spatraster(data = Lc_rast_l[["G-AD-M-LC_01"]])

# Visualize with satellite imagery 
sat_bg <- get_tiles(Lc_rast_l[[1]], crop = TRUE, provider = "Esri.WorldImagery")
plot_tiles(sat_bg) 
plot(Lc_rast_l[[1]], add = TRUE, alpha = .5)

# Functions that could be helpful for exploration
#zoom(new = FALSE)
#draw(x = "points")
#click(id = TRUE)

# Calc metrics ------------------------------------------------------------
## Calculate lsm using Landscapemetrics package

# Generate objects going into scale_sample() function
Pc_vect <- Pc_locs_sf %>% 
  vect() %>% 
  project("EPSG:4686")
Pc_vect_proj <- Pc_vect %>% project("EPSG:3116")
Pc_cents <- Pc_vect_proj %>%
  terra::split("Id_muestreo")

Buffer_rad_nmr <- seq(from = 300, to = 50, by = -50)
Buffer_rad_nmr <- setNames(Buffer_rad_nmr, Buffer_rad_nmr)

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

# Should be no NAs if merge worked appropriately.
Lsm_df %>% filter(is.na(lc_typ2))

# On 2.28.25 there started being a problematic class 4 here in 3 point counts, not sure what that's about. It seems like all others have class starting at 0, and this starts at 1
join_lc_class %>% filter(Id_muestreo == "C-MB-VC-LCH_04")
Lsm_df %>% filter(Id_muestreo == "C-MB-VC-LCH_04") %>% 
  pull(class)

# Inspect Lsm_df ----------------------------------------------------------

# Instead of recalculating (slow), load in lsm file
#Lsm_df <- read_xlsx("Derived/Excels/Lsm/Lsm_df_02.08.25.xlsx")

## NOTE:: There are 'percentage_inside' (pi) over 100 and under 98. Generate 1 row per Id_muestreo to handle more easily 
Low_pi <- Lsm_df %>% filter(percentage_inside < 100) %>% #> 100
  summarize(percent_inside = min(percentage_inside), .by = Id_muestreo) %>%
  mutate(sum_fun = "min") %>%
  arrange(percent_inside)# %>%
  #filter(str_detect(Id_muestreo, "\\(1\\)|OQ")) 
Low_pi

High_pi <- Lsm_df %>% filter(percentage_inside >= 100) %>% 
  summarize(percent_inside = max(percentage_inside), .by = Id_muestreo) %>%
  mutate(sum_fun = "max") %>% 
  arrange(desc(percent_inside))
High_pi

# Create dataframe of ids that are problematic due to percentage inside
prob_pi_ids <- rbind(Low_pi, High_pi) %>% 
  mutate(across(where(is_double), ~ round(.x, 2)))

# Subset the 20 ids with the lowest & highest percent inside
prob_ids_min_max <- prob_pi_ids %>% slice(1:10, .by = sum_fun)

# Plot problematic pc buffers ------------------------------------------------

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
  pdf("Figures/Prob_polys/Polygons_percent_inside.pdf")
  print(marrangeGrob(grobs = prob_pi_plots, ncol = 3, nrow = 3, 
                     layout_matrix = matrix(1:9, 3, 3, TRUE)))
  dev.off()
}

# There are several missing polygons in these Id_groups
Prob_id_groups <- Snapped_lcs %>% 
  filter(Id_group %in% c("G-AD-M-EPO3", "U-MB-M-EPO3"))  
ggplot() + geom_spatvector(data = Prob_id_groups, aes(fill = lc_typ2)) 

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
Low_pi %>% left_join(Pc_locs_sf) %>%
  filter(Id_muestreo != "OQ_Practica") #%>% 
#st_write("Derived_geospatial/shp/To_digitize_natalia/To_digitize.shp")

# Export 
prob_pi_ids %>% as.data.frame() #%>%
  #write.xlsx(file = paste0("Derived/Excels/Lsm/Polys_percent_inside.xlsx"), 
  #showNA = FALSE, row.names = FALSE
#)

# Export Excel of lsm
Lsm_df_exp <- Lsm_df %>% as.data.frame() 
Lsm_df_exp %>%
  if(FALSE){
    write.xlsx(
    file = paste0("Derived/Excels/Lsm/Lsm_df_", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), showNA = FALSE, row.names = FALSE
  )
  }

stop()
warnings() # 50+ warnings: Double values will be converted to integer.
