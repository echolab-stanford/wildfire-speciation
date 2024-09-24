# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: August 16, 2023
# Description: create a map to go along with regional coefficients

# us_states_fp = file.path(data_fp, 'raw/all_national_states.rds')
# loadd(region_pal, cache = drake_cache)

# function
create_us_region_map <- function(us_states_fp, region_pal) {
  
# read in USA shapefile
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


# Create a color palette for the map:
# pal <- colorFactor(wes_palette("Rushmore1", 6, type = "continuous"))
pal <- colorFactor(palette = region_pal, domain = us_regions$region)

region_map <- leaflet(us_regions) %>%
  #addProviderTiles(providers$CartoDB.PositronNoLabels, options(opacity = .7), group = 'base') %>% 
  # leaflet(options = leafletOptions(attributionControl=FALSE)) %>% 
  leaflet.extras::setMapWidgetStyle(list(background = "white")) %>% 
  setView(-96,40, zoom = 4) %>% 
  addPolygons(fillColor = ~pal(region),
              fillOpacity = 0.75, 
              color = 'white',
              weight = 1,
              smoothFactor = 0.5,
              opacity = .4,
              stroke = T,
              group = "us regions") #%>%

# region_map

# output as png in Box folder
mapview::mapshot(
  region_map,
  file = file.path(results_fp, 'Fig3/US_regions.pdf'))

return(us_regions)
}

