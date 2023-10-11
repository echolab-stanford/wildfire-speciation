# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: create a map of the different monitoring stations

# loadd(clean_pm_spec_df, cache = drake_cache)

# create a map of CSN vs improve sites

create_map_of_monitoring_stations <- function(clean_pm_spec_df) {
  
  # get duration of each site:
  monitor_duration <- clean_pm_spec_df %>% 
    mutate(Dataset = ifelse(Dataset == 'EPACSN', "CSN", "IMPROVE")) %>% 
    distinct(Dataset, SiteCode, year) %>% 
    group_by(SiteCode) %>% 
    # get how many years each site is online
    dplyr::summarise(yrs_online = n_distinct(year)) %>% 
    ungroup() %>% 
    # create categorical var for duration
    mutate(duration_cat = case_when(
      yrs_online < 5 ~ 'less than 5 years',
      yrs_online >= 5 & yrs_online < 10 ~ '5-10 years', 
      yrs_online >= 10 ~ '10-15 years')) %>% 
    mutate(duration_cat = factor(duration_cat, levels = c('less than 5 years',
                                                          '5-10 years',
                                                          '10-15 years')))
    #hist(monitor_duration$yrs_online)
  
  # create sf object of sites
  monitor_sf <- clean_pm_spec_df %>% 
    distinct(Dataset, SiteCode, Latitude, Longitude, epsg) %>% 
    # st_as_sf(coords = c('Longitude', 'Latitude'), 
    #          crs = 4326) %>% 
    mutate(Dataset = ifelse(Dataset == 'EPACSN', "CSN", "IMPROVE")) %>% 
    left_join(monitor_duration, by = 'SiteCode') %>% 
    mutate(icon = case_when(
      Dataset == 'CSN' ~ 'circle',
      Dataset == 'IMPROVE' ~ 'triangle')) %>% 
    mutate(color = case_when(
      duration_cat == 'less than 5 years' ~ 'goldenrod',
      duration_cat == '5-10 years' ~ 'firebrick',
      duration_cat == '10-15 years' ~ 'darkorchid4'))
  
  # read in USA shapefile
  us_sf <- readRDS(file.path(data_gdrive_fp, 'boundaries/all_national_states.rds')) %>% 
    st_as_sf() %>% 
    st_transform(crs = 4326) %>% 
    filter(!STUSPS %in% c('AK', 'HI', "PR", 'VI', 'AS', 'GU', 'MP'))

  # create awesome icons
  # my_icons <- awesomeIcons(icon = monitor_sf$icon,
  #                          markerColor = monitor_sf$color,
  #                          library = "glyphicon")
  
  # map <- mapview(monitor_sf, zcol = 'Dataset') # get middle of country coordinates for zoom on leaflet
  # map

  
  # ste up palatte
  # mon_pal <- colorFactor(c('steelblue', 'navy'), monitor_sf$Dataset) # color by dataset
   mon_pal <- colorFactor(c('goldenrod', 'firebrick', 'darkorchid4'), monitor_sf$duration_cat) # color by duration
  
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
                fillColor = "white",
                group = "us states border") %>%
    addCircleMarkers(lng = ~Longitude,
               lat = ~Latitude,
               radius = 4.5,
               fillColor = ~mon_pal(duration_cat),
               #icon = ~makeIcon(Dataset),
               stroke=FALSE,
               fillOpacity = .8) %>%
    addLegend("topright",
              pal = mon_pal,
              values = ~duration_cat,
              labels = "Monitor Type ",
              title = "Duration monitor online")
  mon_map

  
  # output as png in Box folder
  mapshot(
    mon_map,
    file = file.path(wip_gdrive_fp, 'figures/species_monitor_map.png'))
  
  
  return(monitor_sf)
}
  
