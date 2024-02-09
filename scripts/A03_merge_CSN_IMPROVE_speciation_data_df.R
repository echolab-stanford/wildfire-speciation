# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Sept 27, 2023
# Description: Binds together the CSN and IMPROVE speciation data - does some cleaning/flagging of data

# loadd(c(us_shp, raw_CSN_spec_df, raw_IMPROVE_spec_df), cache =  drake::drake_cache(".drake"))


# FUNCTION
merge_CSN_IMPROVE_speciation_data <- function(raw_CSN_spec_df, raw_IMPROVE_spec_df, us_shp) {
  
  all_raw_speciation_df <- bind_rows(raw_CSN_spec_df, 
                                     raw_IMPROVE_spec_df) %>% 
    # drop vars that are only in one monitoring network
    dplyr::select(-c(N2, ammNO3, ammSO4))
  
  # TRANSFORM SITE COORDS TO SAME EPSG (PROJECTION)
  # list all the unique EPSGs
  epsg_vec <- unique(all_raw_speciation_df$epsg)
  #  current_epsg <- epsg_vec[1] # test to see if this works
  
  # map over each crs in the list, assign the appropriate crs, then convert to one crs
  all_spec_transformed_df <- map_df(epsg_vec, function(current_epsg) {
    
    current_sf <- all_raw_speciation_df %>% 
      filter(epsg == current_epsg) %>% 
      # create a spatial dataframe with the given epsg
      st_as_sf(crs = current_epsg, coords = c('long', 'lat')) %>% 
      # transform the dataframe to be WGS84 (or EPSG 4326)
      st_transform(crs = 4326) %>% 
      # grab the newly projected coordinates
      mutate(Longitude = st_coordinates(.)[,1],
             Latitude = st_coordinates(.)[,2]) %>% 
      # drop spatial features of the dataframe
      st_drop_geometry() %>% 
      # reset the projection column to be the newly projected coordinates
      mutate(epsg = 4326)
    
  }) %>% 
    bind_rows() 
  
  # drop sites that are not in the contiguous US, find these sites by intersecting 
  # our datapoints with the US shapefile, any site that no longer has coordinates, drop
 CONUS_sites_df <- all_spec_transformed_df %>%
    distinct(site_id, Longitude, Latitude) %>%
    st_as_sf(coords = c(x = 'Longitude', y = 'Latitude'), crs = 4326) %>% 
   st_intersection(us_shp %>% dplyr::select(state_name = 'NAME_1'))
  
 # now join
 CONUS_spec_df <- all_spec_transformed_df %>% 
   left_join(CONUS_sites_df, by= 'site_id') %>% 
   filter(!is.na(state_name)) %>% 
   st_drop_geometry() %>% 
   dplyr::select(Dataset, state_name, long = 'Longitude', lat= 'Latitude', 
                 epsg, site_id, Date, units, AQSCode, MF, RCFM, AL, AS, BR,
                 CA, CHL, CL, CR, CU, EC, FE, K, MG, MN, `NA`, NI, NO3, OC,
                 P, PB, RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR) %>% 
   distinct()

  # quick check of the map to make sure things look okay
  # mapview(CONUS_spec_df %>%
  #           distinct(site_id, long, lat) %>%
  #           st_as_sf(coords = c(x = 'long', y = 'lat'), crs = 4326),
  #         legend = FALSE, cex = 3)
  #  
  
  return(CONUS_spec_df)
}