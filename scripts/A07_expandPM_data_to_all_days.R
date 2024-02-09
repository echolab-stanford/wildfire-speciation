# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 25, 2025
# Description: This script creates the full universe of days in our sample from 2006-2020, adds in smoke PM
# and then matches each grid cell with region

# us_states_fp = file.path(boundaries_fp, 'all_national_states.rds')
# pm_fp = file.path(data_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')
# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(us_region_map, cache = drake_cache)

# to do - get grid cells in each region and then predict
expandPM_data_to_all_days <- function(pm_fp, grid_fp, us_states_fp) {

# COMBINE GRID CELLS + SMOKE PREDICTIONS:
# smoke pm 
smoke_pm <- readRDS(pm_fp) %>% 
  rename(Date = 'date')

#read in gridded smoke data
grid_10km <- st_read(grid_fp) %>%
  st_transform(4326) 

us_regions <- readRDS(us_states_fp) %>% 
  st_as_sf() %>% 
  st_transform(crs = 4326) %>% 
  filter(!STUSPS %in% c('AK', 'HI', "PR", 'VI', 'AS', 'GU', 'MP')) %>% 
  mutate(region = case_when(
    STUSPS %in% c('OR', 'CA', 'WA') ~ 'pacific',
    STUSPS %in% c('NV', 'UT', 'CO', 'WY', 'MT', 'ID') ~ 'rocky_mountain',
    STUSPS %in% c('AZ', 'NM', 'TX', 'OK') ~ 'southwest',
    STUSPS %in% c('ND', 'SD', 'NE', 'KS', 'MN', 'IA', 'MO', 'WI', 'IL', 'IN', 'MI', 'OH') ~ 'midwest',
    STUSPS %in% c('KY', 'WV', 'VA', 'TN', 'NC', 'MS', 'AL', 'GA', 'SC', 'FL', 'LA', 'AR') ~ 'southeast',
    STUSPS %in% c('ME', 'VT', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'PA', 'DE', 'MD', 'DC') ~ 'northeast',
    STUSPS %in% c('HI', 'AK', 'VI') ~ 'noncontiguous',
    is.na(STUSPS)  ~ 'outside US',
    STUSPS %in% c('AB', 'ON') ~ 'outside US'
  )) %>%  
  # Exclude regions that are  'noncontiguous' and 'outside US'
  filter(!region %in% c('noncontiguous', 'outside US')) %>% 
  dplyr::select(st_ab = 'STUSPS', st_name = 'NAME', region)

# group geometry into regions + remake polygons
regions_sf <- us_regions %>% 
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

# create a dataframe that has all the dates between 2006-2020, not just smoke days, 
# and then create 0 values on non-smoke days
full_days_df <- tibble(
  Date = seq(as.Date("2006-01-01"), as.Date("2020-12-31"), by = "days")
) %>% 
  # create a date for each grid cell for all days in sample
  crossing(unique(grid_10km_region_sf$grid_id_10km)) %>% 
  rename(grid_id_10km = `unique(grid_10km_region_sf$grid_id_10km)`) %>% 
  # merge smoke date
  left_join(smoke_pm, by = c('Date','grid_id_10km')) %>% 
  # for any days that are NA, turn to 0
  mutate(smokePM_pred = ifelse(is.na(smokePM_pred), 0, smokePM_pred)) %>% 
  rename(smokePM = 'smokePM_pred') %>% 
  # add in which region the grid cell is in
  left_join(grid_10km_region_sf %>% 
              st_drop_geometry(), by = 'grid_id_10km')

return(full_days_df)
}