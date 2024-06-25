# Emma Krasovich Southworth, emmars@stanford.edu
# Created: Jan 16, 2023 | Last Updated: Jan 16, 2023
# Description: merge speciation data with plume data, basic cleaning

# loadd(CONUS_spec_df, cache = drake_cache)
# pm_fp = file.path(data_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')
# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')

join_speciation_w_gridded_smokePM <- function(CONUS_spec_df, pm_fp, grid_fp) {

  # NOW PREPARE MARISSA's PREDICTED SMOKE PM TO BE MERGED WITH SPECIATION DATA
  #### read in PM data ####
  # wildfire pm data
  smoke_pm <- readRDS(pm_fp) %>% 
    rename(Date = 'date') %>% 
    distinct()
  
  # Load 10 km grid
  grid_10km <- st_read(grid_fp) %>% 
    st_transform(unique(CONUS_spec_df$epsg))
  
  #### Match monitor locations in speciation data with grid 10km ID ####
  # make points spatial data for monitors in df
  points <- CONUS_spec_df %>% 
    # get a unique list of sites
    distinct(Dataset, site_id, long, lat) %>% 
    # convert so sf object
    st_as_sf(coords = c(x ='long', y="lat"), crs = st_crs(grid_10km)) 
  

  # pull id cells that intersect with the monitoring sites
  pts_in_grid_df <- points %>% 
    st_join(grid_10km) %>% 
    st_drop_geometry() %>% 
    rename(grid_id_10km = 'ID') 
  
  # Add grid_id_10km to the speciation df, join the smoke PM data by grid cell
  spec_w_smoke_pm_df <- CONUS_spec_df %>% 
    mutate(Date = as.Date(Date, format = '%Y-%m-%d')) %>% 
    left_join(pts_in_grid_df, 
              by = c("Dataset", "site_id")) %>% 
    # match PM2.5 from the smoke_pm dataset by gridcell
    left_join(smoke_pm, 
              by = c("grid_id_10km", "Date")) %>% 
    # Replace NA values (non-smoke days) to 0 (no smoke detected)
    mutate(smokePM2.5 = ifelse(is.na(smokePM_pred), 0, smokePM_pred)) %>% 
    # drop extra vars
    dplyr::select(-c(AQSCode, grid_id_10km)) %>% 
    distinct()


  return(spec_w_smoke_pm_df)
}
