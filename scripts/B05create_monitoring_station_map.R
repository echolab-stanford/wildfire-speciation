# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: create a map of the different monitoring stations

# loadd(clean_pm_spec_df, cache = drake_cache)

# create a map of CSN vs improve sites

create_map_of_monitoring_stations <- function(clean_pm_spec_df) {
  
  # create sf object of sites
  monitor_sf <- clean_pm_spec_df %>% 
    distinct(Dataset, SiteCode, Latitude, Longitude, epsg) %>% 
    # st_as_sf(coords = c('Longitude', 'Latitude'), 
    #          crs = 4326) %>% 
    mutate(Dataset = ifelse(Dataset == 'EPACSN', "CSN", "IMPROVE"))
  
  # read in USA shapefile
  us_sf <- readRDS(file.path(data_gdrive_fp, 'boundaries/all_national_states.rds')) %>% 
    st_as_sf() %>% 
    st_transform(crs = 4326) %>% 
    filter(!STUSPS %in% c('AK', 'HI', "PR", 'VI', 'AS', 'GU', 'MP'))

  # map <- mapview(monitor_sf, zcol = 'Dataset') # get middle of country coordinates for zoom on leaflet
  # map
  
  # ste up palatte
  mon_pal <- colorFactor(c('steelblue', 'forestgreen'), monitor_sf$Dataset) 
  
  mon_map <- leaflet(monitor_sf) %>%
    #addProviderTiles(providers$CartoDB.PositronNoLabels, options(opacity = .7), group = 'base') %>% 
    # leaflet(options = leafletOptions(attributionControl=FALSE)) %>% 
    leaflet.extras::setMapWidgetStyle(list(background = "white")) %>% 
    setView(-96,40, zoom = 4) %>% 
    addPolygons(data = us_sf,
                color = 'black',
                weight = .7,
                smoothFactor = 0.5,
                opacity = .9,
                stroke = T,
                fillColor = "grey100",
                group = "us states border") %>%
    addCircleMarkers(lng = ~Longitude, 
                     lat = ~Latitude, 
                     radius = 4.5, 
                     fillColor = ~mon_pal(Dataset),
                     stroke=FALSE,
                     fillOpacity = .8) %>%
    addLegend("topright", 
              pal = mon_pal, 
              values = ~Dataset, 
              labels = "Monitor Type ", 
              title = "Monitor Type")
  mon_map

  
  # output as png in Box folder
  mapshot(
    mon_map,
    file = file.path(wip_gdrive_fp, 'figures/species_monitor_map.png'))
  
  
  return(monitor_sf)
}
  
