# Emma Krasovich Southworth, emmars@stanford.edu
# Created: Jan 16, 2023 | Last Updated: Jan 16, 2023
# Description: merge speciation data with plume data, basic cleaning

# loadd(c(cleaned_spec_df, improve_smoke_dens_fire_dist), cache = drake_cache)
# pm_fp = file.path(wip_gdrive_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')
# grid_fp = file.path(wip_gdrive_fp, 'intermediate/10km_grid_wgs84.shp')

merge_speciation_plumes_and_pm <- function(cleaned_spec_df, improve_smoke_dens_fire_dist, pm_fp, grid_fp) {

  # MERGE PLUMES W SPECIATION
  plumes_spec_df <- improve_smoke_dens_fire_dist %>% 
    left_join(cleaned_spec_df, 
              by = c('Dataset', 'SiteCode', 'Date', 'Elevation', 'Longitude', 'Latitude', 'epsg')) 
  
  # NOW PREPARE MARISSA's PREDICTED SMOKE PM TO BE MERGED WITH SPECIATION DATA
  #### read in PM data ####
  # wildfire pm data
  pm <- readRDS(pm_fp) %>% 
    rename(Date = date)
  
  # Load 10 km grid
  grid_10km <- st_read(grid_fp) %>% 
    st_transform(unique(plumes_spec_df$epsg))
  
  #### Match monitor locations in speciation data with grid 10km ID ####
  # make points spatial data for monitors in df
  points <- plumes_spec_df %>% 
    # get a unique list of sites
    distinct(Dataset, SiteCode, Longitude, Latitude) %>% 
    # convert so sf object
    st_as_sf(coords = c('Longitude',"Latitude"), crs = unique(plumes_spec_df$epsg)) 
  

  # pull id cells that interesect with the monitoring sites
  points$grid_id_10km <- apply(st_intersects(grid_10km, points, sparse = FALSE), 2, 
                               function(col) {grid_10km[which(col), ]$ID})
  

  # Add grid_id_10km to the speciation df
  speciation_w_pm_df <- left_join(plumes_spec_df, points, 
                          by = c("Dataset", "SiteCode")) %>% 
    mutate(grid_id_10km = as.integer(grid_id_10km)) %>%
    # match PM2.5 from the pm dataset by gridcell
    left_join(pm, 
              by = c("grid_id_10km", "Date")) %>% 
    st_drop_geometry() %>% 
    dplyr::select(-c(LandUseCode, geometry, grid_id_10km, EPACode)) %>% 
    # Replace NA values (non-smoke days) to 0 (no smoke detected)
    mutate(smokePM_pred = ifelse(is.na(smokePM_pred), 0, smokePM_pred)) %>% 
    dplyr::select(-c(low_count, med_count, high_count, density_missing))
  
  
  return(speciation_w_pm_df)
}
