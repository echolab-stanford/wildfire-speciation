# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Sept 25, 2023
# Description: This script reads in downloaded speciation data on HAPS from
# https://aqs.epa.gov/aqsweb/airdata/download_files.html



# loadd(c(us_shp,spec_w_smoke_pm_df), cache = drake_cache)
# haps_spec_file_list = list.files(file.path(data_fp, "raw"),
#                                  pattern = "daily_HAPS", full.names= TRUE)
# haps_xwalk_fp = file.path(data_fp, 'intermediate/haps_xwalk_IUR_weight.csv')

# FUNCTION ------------------------------------------------------------------------------------
read_in_HAPS_speciation_data <- function(us_shp, haps_spec_file_list, haps_xwalk_fp, spec_w_smoke_pm_df) {
 
  # step 1 ----------------------------------------------
  # testing one file at a time
  # current_file <- haps_spec_file_list[1]
  
  # read each in file, clean up, and ensure that parameters match those in the improve data
  haps_spec_files = map_df(haps_spec_file_list, function(current_file) {
    
    # print current index
    print(paste0("current index is: ", match(current_file, unique(haps_spec_file_list)), ' of ', length(haps_spec_file_list)))
    
    # read in text files on EPA CSN chemical speciation
    temp_spec_df = read.csv(current_file) %>% 
      dplyr::select(State.Code, County.Code, Site.Number = 'Site.Num', 
                    Parameter = 'Parameter.Name',
                    AQSParamCode = 'Parameter.Code', 
                    lat = 'Latitude', long = 'Longitude', Datum, 
                    Date = 'Date.Local', sample_dur = 'Sample.Duration',
                    units = 'Units.of.Measure', 
                    conc = 'Arithmetic.Mean', 
                    st_name = 'State.Name', co_name = County.Name) %>% 
      filter(sample_dur == '24 HOUR') %>% 
      # recreate site code to merge with later steps
      mutate(State.Code = str_remove(State.Code, "^0")) %>% 
      mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
             Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
      mutate(SiteCode = paste0("0", State.Code, County.Code, Site.Number)) %>% 
      mutate(SiteCode = case_when(
        str_length(SiteCode) <= 9 ~ SiteCode,
        str_length(SiteCode) > 9 ~ str_remove(SiteCode, "^0")
      )) %>% 
      dplyr::select(-c(State.Code, County.Code, Site.Number, sample_dur)) %>% 
      # make sure we only have one obs per site-date
      group_by(SiteCode, Date, st_name, lat, long, Datum, units, Parameter) %>% 
      dplyr::summarise(conc = mean(conc, na.rm = TRUE), .groups = 'drop') 
       

  }) %>% 
    bind_rows()  
  
# wider
haps_spec_df <- haps_spec_files %>% 
  mutate(flag = ifelse(
    str_detect(Parameter, "PM10"), 'drop', 'keep'
  )) %>% 
  filter(flag == 'keep') %>% 
  filter(units != 'Nanograms/cubic meter (25 C)') %>% 
  pivot_wider(id_cols = c('SiteCode', 'Date', 'st_name',
                          'lat', 'long', 'Datum', 'units'),
              names_from = 'Parameter',
              values_from = 'conc') %>% 
  dplyr::select(SiteCode:units, 
                Cd = "Cadmium PM2.5 LC",
                Mn = "Manganese PM2.5 LC",
                Be = "Beryllium PM2.5 LC",
                Acetaldehyde:`trans-1,3-Dichloropropene`, 
                Acrolein = 'Acrolein - Verified') %>% 
  filter_at(vars(Cd:Acrolein),
            any_vars(!is.na(.))) %>% 
  rename(site_id = 'SiteCode') %>% 
  mutate(Dataset = 'CSN') %>% 
  mutate(epsg = case_when(
           Datum == 'WGS84' ~ 4326,
           Datum == 'NAD83' ~ 4629)) %>% 
  dplyr::select(-Datum)

# update the crs:
# TRANSFORM SITE COORDS TO SAME EPSG (PROJECTION)
# list all the unique EPSGs
epsg_vec <- unique(haps_spec_df$epsg)
#  current_epsg <- epsg_vec[1] # test to see if this works

# map over each crs in the list, assign the appropriate crs, then convert to one crs
all_spec_transformed_df <- map_df(epsg_vec, function(current_epsg) {
  
  current_sf <- haps_spec_df %>% 
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
  st_as_sf(coords = c(x = 'Longitude', 
                      y = 'Latitude'), crs = 4326) %>% 
  st_join(us_shp %>% 
            st_transform(4326) %>% 
            dplyr::select(state_name = 'NAME_1')) %>% 
  filter(!is.na(state_name)) %>% 
  mutate(Longitude = st_coordinates(.)[,1],
         Latitude = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  distinct()

# now join
CONUS_haps_spec_df <- CONUS_sites_df %>% 
  left_join(all_spec_transformed_df %>% 
              dplyr::select(-c(Longitude, Latitude)) %>% 
              distinct(), by= 'site_id') %>% 
  pivot_longer(cols = c(Cd:Acrolein), names_to = 'species', values_to = 'conc') %>% 
  left_join(read.csv(haps_xwalk_fp), by = 'species') %>% 
  filter(!is.na(weight_g_mol)) %>% 
   # Concentration (µg/m3) = molecular weight x concentration (ppb) ÷ 24.45
  mutate(conc_ug_m3 = case_when(
    units == 'Parts per billion Carbon' ~ (conc*weight_g_mol)/24.45,
    TRUE ~ conc
  )) %>% 
  mutate(units = 'ug/m3') %>% 
  dplyr::select(-conc, -CAS, -IUR, -weight_g_mol) %>% 
  distinct() %>% 
  filter(!is.na(conc_ug_m3)) %>% 
  filter(site_id %in% spec_w_smoke_pm_df$site_id) %>% 
  pivot_wider(names_from = 'species', values_from = 'conc_ug_m3')

  return(CONUS_haps_spec_df)

} # end function

