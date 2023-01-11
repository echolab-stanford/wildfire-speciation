# Description: This script reads in the location of AQS sites + monitors to creates a dataframe of CSN + IMPROVE sites 
# 
# aqs_sites_fp = file.path(wip_gdrive_fp, 'raw/aqs_sites.csv')
# monitors_fp = file.path(wip_gdrive_fp, 'raw/aqs_monitors.csv')
# csn_sites_fp = file.path(wip_gdrive_fp, 'raw/csn_sites.xlsx')
# improve_sites_fp = file.path(wip_gdrive_fp, 'raw/improve_sites.xlsx')
# improve_prelim_sites_fp = file.path(wip_gdrive_fp, 'raw/improve_prelim_sites.xlsx')
# loadd(raw_spec_df, cache = drake_cache)


# function
merge_projections_w_spec_df <- function(raw_spec_df, 
                                        aqs_sites_fp, 
                                        monitors_fp,
                                        csn_sites_fp, 
                                        improve_sites_fp, 
                                        improve_prelim_sites_fp) {
  
  
  # STEP 1) READ IN AQS INFORMATION WHICH CONTAINS DATUMS  
  # READ IN AQS SITE INFO
  aqs_sites <- read.csv(aqs_sites_fp) %>% 
    # select relevant columns
    dplyr::select(State.Name, State.Code, County.Code, Site.Number, Latitude, Longitude, Datum, Elevation) %>% 
    mutate(State.Code = str_remove(State.Code, "^0")) %>% 
    mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
           Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
    mutate(SiteCode = paste0("0", State.Code, County.Code, Site.Number))
  
  # read in AQS monitor locations
  aqs_monitors <- read.csv(monitors_fp) %>% 
    distinct(State.Code, County.Code, Site.Number, Latitude, Longitude, Datum, Networks) %>% 
    mutate(State.Code = str_remove(State.Code, "^0")) %>% 
    mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
           Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
    mutate(SiteCode = paste0("0", State.Code, County.Code, Site.Number))
  
  # merge monitors and sites together
  aqs_merged <- left_join(aqs_monitors, aqs_sites, 
                          by = c("Latitude", "Longitude", "Datum", 
                                 "SiteCode", "Site.Number",
                                 "State.Code", "County.Code")) %>% 
    mutate(epsg = case_when(
      Datum == 'WGS84' ~ 4326,
      Datum == 'NAD83' ~ 4629,
      Datum == 'NAD27' ~ 4267,
      is.na(Datum) ~ 4326)) %>% # assume 4326 projection for missing datums
    filter(!is.na(epsg)) %>% 
    # add in a column replicating sitecode because sometime site code is not the same as epa code
    mutate(SiteCode = case_when(
      str_length(SiteCode) <= 9 ~ SiteCode,
      str_length(SiteCode) > 9 ~ str_remove(SiteCode, "^0")
    )) %>% 
    distinct()
  
  # STEP 2) READ IN CSN SITE LIST OF ALL SITES AND MATCH WITH DATUM
  csn_sites_epsg <- read_excel(csn_sites_fp) %>%
    distinct(AQSCode, LandUseCode, Latitude, Longitude, Elevation) %>%
    rename(SiteCode = 'AQSCode') %>%
    left_join(aqs_merged %>% 
                distinct(epsg, SiteCode), 
              by = 'SiteCode') %>%
    rename(site_code = 'SiteCode') %>% 
    mutate(Dataset = 'EPACSN')
  
  # STEP 3) READ IN IMPROVE SITE LIST OF ALL SITES AND MATCH WITH DATUM
    improve_sites_epsg <- read_excel(improve_sites_fp) %>%
      mutate(Dataset = 'IMPROVE') %>% 
      bind_rows(read_excel(improve_prelim_sites_fp) %>% 
                  mutate(Dataset = 'IMPROVE')) %>%
      filter(Country =='US') %>%
      distinct(SiteCode, AQSCode, Dataset) %>%
      rename(SiteNameID = 'SiteCode',
             SiteCode = 'AQSCode') %>%
      left_join(aqs_merged %>% 
                  distinct(epsg, SiteCode), 
                by = 'SiteCode') %>%
      # create new var to identify sites uniquely
      # mutate(site_code = ifelse(is.na(SiteCode), SiteNameID, SiteCode)) %>% 
      mutate(epsg = ifelse(is.na(epsg), 4326, epsg))  # assume 4326 projecttion for missing datums
  
    #STEP 4) MERGE WITH OUR SPECIATION DATA:
    # improve sites
    improve_sites_proj <- raw_spec_df %>% 
      distinct(SiteCode, Dataset) %>% 
      filter(Dataset == 'IMPROVE') %>% 
      rename(site_code = SiteCode) %>% 
      left_join(improve_sites_epsg %>% 
                  rename(site_code = 'SiteNameID'), 
                by = c('Dataset', 'site_code')) %>% 
      # assume 4326 projection for missing datum
      mutate(epsg = ifelse(is.na(epsg), 4326, epsg)) %>% 
      distinct(Dataset, site_code,epsg)
      
    # csn sites
    csn_sites_proj <- raw_spec_df %>% 
      distinct(SiteCode, Dataset) %>% 
      filter(Dataset == 'EPACSN') %>% 
      rename(site_code = SiteCode) %>% 
      left_join(csn_sites_epsg, 
                by = c('Dataset', 'site_code')) %>% 
      filter(!is.na(Latitude)) %>% 
      bind_rows(raw_spec_df %>% 
                  distinct(SiteCode, Dataset) %>% 
                  filter(Dataset == 'EPACSN') %>% 
                  rename(site_code = SiteCode) %>% 
                  left_join(csn_sites_epsg, 
                            by = c('Dataset', 'site_code')) %>% 
                  filter(is.na(Latitude)) %>% 
                  distinct(site_code, Dataset) %>%
                  rename(SiteCode = site_code) %>% 
                  left_join(aqs_merged %>% 
                              distinct(SiteCode, Latitude, Longitude, Elevation, epsg), 
                            by = 'SiteCode') %>% 
                  rename(site_code = SiteCode)
                )
        
    # all speciation sites in our dataset with their epsg's
    all_monitor_sites <- improve_sites_proj %>% 
      rename(SiteCode = site_code) %>% 
      left_join(raw_spec_df %>% 
                  distinct(SiteCode,Latitude, Longitude, Dataset, Elevation),
                by = c('Dataset', 'SiteCode')) %>% 
      bind_rows(csn_sites_proj %>%
                  rename(SiteCode = site_code))
    
    # now merge in with speciation data to provide the correct crs for the speciation data:
    all_spec_w_proj <- raw_spec_df %>% 
      dplyr::select(-c(Latitude, Longitude, Elevation)) %>% 
      left_join(all_monitor_sites, by= c('SiteCode', 'Dataset')) %>% 
      filter(!is.na(Latitude))
    
    
    # NOW TRANSFORM TO SAME EPSG
    # list all the unique EPSGs
    epsg_vec <- unique(all_spec_w_proj$epsg)
    #  current_epsg <- epsg_vec[1] # test to see if this works
    
    # map over each crs in the list, assign the appropriate crs, then convert to one crs
    all_spec_proj_clean_df <- map_df(epsg_vec, function(current_epsg) {
      
      current_sf <- all_spec_w_proj %>% 
        filter(epsg == current_epsg) %>% 
        st_as_sf(crs = current_epsg, coords = c('Longitude', 'Latitude')) %>% 
        st_transform(crs = 4326) %>% 
        mutate(Longitude = st_coordinates(.)[,1],
               Latitude = st_coordinates(.)[,2]) %>% 
        st_drop_geometry() %>% 
        mutate(epsg = 4326)
        
      
    }) %>% 
      bind_rows() 
    
    cleaned_spec_df <- all_spec_proj_clean_df %>% 
      filter(!SiteCode %in% c('VIIS1', '720210010', 'BYIS1', 'BYISX'))
    
    # quick check of the map to make sure things look okay
    # mapview(cleaned_spec_df %>%
    #           distinct(SiteCode, Longitude, Latitude) %>% 
    #           st_as_sf(coords = c(x = 'Longitude', y = 'Latitude'), crs = 4326),
    #         legend = FALSE, cex = 3)
            
  return(cleaned_spec_df)

} # end function

