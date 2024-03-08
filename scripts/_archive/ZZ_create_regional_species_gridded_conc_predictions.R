# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: November 14, 2023
# Description: This script uses our estimated coefficients to predict the concentration of species across a grid


# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(c(regionalPMcoeffs_normalized, us_region_map, full_samplePM_df), cache = drake_cache)

# to do - get grid cells in each region and then predict
create_regional_species_conc_predictions <- function(grid_fp, regionalPMcoeffs_normalized, us_region_map, full_samplePM_df) {

  # COMBINE GRID CELLS WITH REGIONS
  #read in gridded smoke data
  grid_10km <- st_read(grid_fp) %>%
    st_transform(4326) 
  
  # group geometry into regions + remake polygons
  regions_sf <- us_region_map %>% 
    group_by(region) %>% 
    dplyr::summarise(reg_geom = st_union(geometry)) %>% 
    ungroup() %>% 
    st_make_valid()
  
  # determine which grid cells are in which regions by doing a spatial join
  grid_10km_region_sf <- grid_10km %>% 
    st_join(regions_sf) %>% 
    filter(!is.na(region)) %>% 
    rename(grid_id_10km = 'ID') %>% 
    dplyr::select(-COORDX, -COORDY) 
  

  # create species list
  species_list <- c('AS', 'NI', 'OC', 'PB') 
  # species_list <- c('PB') 
  # current_species <- species_list[1]
  
  # map over each species
  all_species_regions_preds_df <- map_df(species_list, function(current_species) {
   
     print(current_species)
   
  # get a list of all the regions  
  region_list <- unique(grid_10km_region_sf$region)
  
  # current_region <- region_list[1]

  # map over each region, match grid cells w/ smoke data + predict over our sample but with regional betas
  all_regions_preds_df <- map_df(region_list, function(current_region) {

    print(current_region)
  
  # filter geometry to current region
  current_regions_sf <- grid_10km_region_sf %>%
    filter(region == current_region) 
  
  # filter dataframe of betas to current species + grab estimate
  current_species_betaS <- regionalPMcoeffs_normalized %>%
    filter(species == current_species) %>% 
    filter(region == current_region) %>% 
    filter(pm_type == 'smokePM') %>% 
    distinct(Estimate) %>% 
    as.numeric()
  
  # grab long name
  current_species_name <- regionalPMcoeffs_normalized %>%
    filter(species == current_species) %>% 
    filter(region == current_region) %>% 
    distinct(species_long) %>% 
    as.character()
  
  # grab short name
  current_species_code <- regionalPMcoeffs_normalized %>%
    filter(species == current_species) %>% 
    filter(region == current_region) %>% 
    distinct(species) %>% 
    as.character()

  # predict species concentration due to smoke in each grid cell
  predicted_gridded_species_conc_df <- full_samplePM_df %>%
    filter(region == current_region) %>% 
    mutate(pred_grid_conc = current_species_betaS*smokePM,
           species_long = paste0(current_species_name),
           species = paste0(current_species_code))
  
  write_fst(predicted_gridded_species_conc_df, 
            file.path(data_fp, 
                      paste0("intermediate/regional_gridded", current_species_code, "_", current_region, ".fst")))

 # summarize across full sample to get avg exposure over 15 yrs
 sample_avg_predictions <- predicted_gridded_species_conc_df %>%
   group_by(grid_id_10km, region, species_long, species) %>%
   dplyr::summarise(sample_avg_pred_conc = mean(pred_grid_conc, na.rm = TRUE)) %>%
   ungroup() 
 
}) %>%
  bind_rows() # bind for all species within a given region

}) %>%
  bind_rows() # bind all predictions together for all regions
  
  all_species_regions_preds_sf <- all_species_regions_preds_df %>% 
    left_join(grid_10km_region_sf, by = c('region','grid_id_10km'))
  
  return(all_species_regions_preds_sf)

}





