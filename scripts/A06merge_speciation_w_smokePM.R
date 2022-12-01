# description: merge speciation data with smoke PM data 

# pm_fp = file.path(wip_gdrive_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')
# grid_fp = file.path(wip_gdrive_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(cleaned_spec_df, cache = drake_cache)

merge_speciation_w_smokePM <- function(cleaned_spec_df, pm_fp, grid_fp) {

#### read in PM data ####
# wildfire pm data
pm <- readRDS(pm_fp) %>% 
  rename(Date = date)
# Load 10 km grid
grid_10km <- st_read(grid_fp)


#### Match monitor locations in speciation data with grid 10km ID ####
# make points spatial data for monitors in df
points <- cleaned_spec_df %>% 
  # get a unique list of sites
  distinct(Dataset, SiteCode, Longitude, Latitude) %>% 
  # convert so sf object
  st_as_sf(coords = c('Longitude',"Latitude"), crs = unique(cleaned_spec_df$epsg)) 


points$grid_id_10km <- apply(st_intersects(grid_10km, points, sparse = FALSE), 2, 
                             function(col) {grid_10km[which(col), ]$ID})

# Add grid_id_10km to the speciation df
speciation <- left_join(cleaned_spec_df, 
                        points, by = c("Dataset", "SiteCode")) %>% 
  mutate(grid_id_10km = as.integer(grid_id_10km))

# Merge files to get pm prediction and add to speciation df
pm_speciation_df <- left_join(speciation, pm, 
                        by = c("grid_id_10km", "Date")) %>% 
  # get duration that a site is online
  mutate(year = year(Date),
         month = month(Date),
         doy = yday(Date), 
         spec_units = 'ug_m3') %>% 
  filter(year > 2005) %>%  # we dont have smokePM data before this
  mutate(season = case_when(
    month %in% c(12,1,2) ~ 'winter',
    month %in% c(3,4,5) ~ 'spring',
    month %in% c(6,7,8) ~ 'summer',
    month %in% c(9,10,11) ~ 'autumn',
  )) %>% 
  mutate(region = case_when(
    State %in% c('OR', 'CA', 'WA') ~ 'pacific',
    State %in% c('NV', 'UT', 'CO', 'WY', 'MT', 'ID') ~ 'rocky_mountain',
    State %in% c('AZ', 'NM', 'TX', 'OK') ~ 'southwest',
    State %in% c('ND', 'SD', 'NE', 'KS', 'MN', 'IA', 'MO', 'WI', 'IL', 'IN', 'MI', 'OH') ~ 'midwest',
    State %in% c('KY', 'WV', 'VA', 'TN', 'NC', 'MS', 'AL', 'GA', 'SC', 'FL', 'LA', 'AR') ~ 'southeast',
    State %in% c('ME', 'VT', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'PA', 'DE', 'MD', 'DC') ~ 'northeast',
    State %in% c('HI', 'AK', 'VI') ~ 'noncontiguous'
  )) %>% 
  distinct() %>% 
  dplyr::select(Dataset, SiteCode, Date, year, month, doy, season, region, State, Latitude, Longitude, 
                epsg, Elevation, grid_id_10km, smokePM_pred, AL, AS, BR, CA, EC1, EC2, EC3,
                EC, OC1,OC2, OC3, OC4, OC, OP, CL, CR, CU, FE, PB, MG, MN, MF, MO, NI,
                RCFM, NO3, P, K, RB, `NA`, SE, SI, S, SR, SOIL, SO4, TI, V, ZN, ZR, 
                ammNO3, ammSO4, TC, CHL, MT, RCTM, N2, CM_calculated, spec_units)

return(pm_speciation_df)
} # end function
