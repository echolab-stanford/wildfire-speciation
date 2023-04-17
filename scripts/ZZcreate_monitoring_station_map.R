# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: create a map of the different monitoring stations

# loadd(clean_pm_spec_df, cache = drake_cache)

# create a map of CSN vs improve sites

create_map_of_monitoring_stations <- function(clean_pm_spec_df) {
  
  # create sf object of sites
  monitor_sf <- clean_pm_spec_df %>% 
    distinct(Dataset, SiteCode, Latitude, Longitude, epsg) %>% 
    st_as_sf(coords = c('Longitude', 'Latitude'), 
             crs = unique(clean_pm_spec_df$epsg))
  
  
  map <- mapview(monitor_sf, zcol = 'Dataset')
  map
  
  
  map <- leaflet(options = leafletOptions(attributionControl=FALSE)) %>%
    addProviderTiles(providers$CartoDB.PositronNoLabels, options(opacity = .5), group = 'base') %>% 
    leaflet.extras::setMapWidgetStyle(list(background = "white")) %>% 
    addPolygons(data = monitor_sf,
                weight = .2, # control the weight of the borders
                smoothFactor = 0.5,
                opacity = .9,
                fillOpacity = 0.9,
                fillColor = ~var_pal(bin_label),
                group = "Precipitation") 

  map
  
  
  # output as png in Box folder
  mapshot(
    map,
    file = file.path(
      local_box_path, file.path("results", 'bangladesh_wbb_precip_white_bg.png'))
  )
  
  
  
  
}
  
