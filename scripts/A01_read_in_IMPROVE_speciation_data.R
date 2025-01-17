# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Sept 27, 2023
# Description: This script reads in downloaded IMPROVE speciation data from: https://views.cira.colostate.edu/fed/DataFiles/
# http://views.cira.colostate.edu/fed/Express/ImproveData.aspx


# testing function
# improve_spec_file_list = list.files(file.path(data_fp, 'raw'), pattern = 'IMPAER', full.names = TRUE)
# improve_site_fp = file.path(data_fp, 'raw/IMPROVE_sites.xlsx')
# aqs_sites_fp = file.path(data_fp, 'raw/aqs_sites.csv')
# aqs_monitors_fp = file.path(data_fp, 'raw/aqs_monitors.csv')

# FUNCTION
read_in_IMPROVE_specation_data <- function(improve_spec_file_list, 
                                           improve_site_fp,
                                           aqs_monitors_fp,
                                           aqs_sites_fp) {

# step 1 -----------------------------------------------------------------------
# read in improve speciation files and clean up
  
# testing one file at a time
# current_file <- improve_spec_file_list[3]

# read each in file, clean up, and merge all speciation data from IMPROVE together
improve_spec_files = map_df(improve_spec_file_list, function(current_file) {
  
  # print current index
  print(paste0("current index is: ", match(current_file, unique(improve_spec_file_list)), ' of ', length(improve_spec_file_list)))
  
  temp_spec_df <- read.csv(current_file, 
                             sep = "|", 
                             header= TRUE) %>% 
    dplyr::select(SiteCode, Date = 'FactDate', ParamCode, conc = 'FactValue', Units) %>% 
    mutate(ParamCode = str_remove(ParamCode, pattern = 'f$')) %>% 
    filter(Units == 'ug/m^3') %>% 
    filter(conc != '-999') %>% 
    mutate(conc = as.numeric(conc)) %>% 
    group_by(SiteCode, Date, ParamCode, Units) %>% 
    dplyr::summarise(conc = mean(conc, na.rm = TRUE), .groups = 'drop') %>% 
    # pivot wide again so that each monitor-day has one value for each measured species
    pivot_wider(id_cols = c('SiteCode', 'Date', 'Units'), 
                names_from = 'ParamCode', 
                values_from = 'conc') 
    
}) %>% 
  bind_rows()


# step 2 -----------------------------------------------------------------------
# get spatial data from other datasets on the AQS website

# read in AQS sites to get the spatial data associated with each site
aqs_monitors <- read.csv(aqs_monitors_fp, na.strings="") %>% 
  # recreate site code to merge with later steps
  mutate(State.Code = str_remove(State.Code, "^0")) %>% 
  mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
         Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
  mutate(AQSCode = paste0("0", State.Code, County.Code, Site.Number)) %>% 
  mutate(AQSCode = case_when(
    str_length(AQSCode) <= 9 ~ AQSCode,
    str_length(AQSCode) > 9 ~ str_remove(AQSCode, "^0")
  )) %>% 
  distinct() %>% 
  filter(!is.na(Datum)) %>% 
  dplyr::select(State.Name, State.Code, County.Code, AQSCode,Site.Number, Latitude, Longitude, Datum) 


# read in additional AQS Site data that is missing in the monitoring data
aqs_sites <- read.csv(aqs_sites_fp, na.strings="") %>% 
  distinct(State.Code, County.Code, Site.Number, Latitude, Longitude, Datum) %>% 
  # recreate site code to merge with later steps
  mutate(State.Code = str_remove(State.Code, "^0")) %>% 
  mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
         Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
  mutate(AQSCode = paste0("0", State.Code, County.Code, Site.Number)) %>% 
  mutate(AQSCode = case_when(
    str_length(AQSCode) <= 9 ~ AQSCode,
    str_length(AQSCode) > 9 ~ str_remove(AQSCode, "^0")
  )) %>% 
  distinct() 


# merge all missing spatial data together (monitors and sites together)
aqs_merged <- aqs_monitors %>% 
  bind_rows(aqs_sites %>% 
              filter(!AQSCode %in% aqs_monitors$AQSCode)) %>% 
  distinct(Latitude, Longitude, Datum, AQSCode, Site.Number, Site.Number) %>% 
  mutate(epsg = case_when(
    Datum == 'WGS84' ~ 4326,
    Datum == 'NAD83' ~ 4629,
    Datum == 'NAD27' ~ 4267,
    is.na(Datum) ~ 4326)) %>% # assume 4326 projection for missing datums
    filter(Latitude != 0) %>% 
  distinct() 


# Step 3 -----------------------------------------------------------------------
# merge with the improve sites so that each site has lat, long and datum
improve_sites <- read_xlsx(improve_site_fp) %>%
  filter(Country == 'US') %>% 
  dplyr::select(SiteCode, st_name = "State", co_name = "County", AQSCode) %>%
  distinct() %>% 
  left_join(aqs_merged, by = 'AQSCode') %>% 
  filter(!is.na(AQSCode))
  

# create AQS codes for the missing ones
missing_AQS <- read_xlsx(improve_site_fp) %>%
  filter(Country == 'US') %>%
  dplyr::select(SiteCode, st_name = "State", co_name = "County", AQSCode, Latitude, Longitude) %>%
  distinct() %>%
  filter(is.na(AQSCode))

# bind together and clean up
all_improve_site_spatial <- bind_rows(missing_AQS, improve_sites) %>% 
  # assume EPSG 4326 for missing datums
  mutate(epsg = ifelse(is.na(epsg), 4326, epsg)) %>% 
  mutate(site_id = SiteCode) %>% 
  dplyr::select(-c(Datum, Site.Number, SiteCode, co_name))


# step 4 -----------------------------------------------------------------------
# merge the aqs sites with the improve sites + speciation data
raw_IMPROVE_spec_df <- improve_spec_files %>% 
  mutate(site_id = SiteCode) %>% 
  left_join(all_improve_site_spatial, by = c('site_id')) %>% 
  filter(!is.na(Latitude)) %>% 
  mutate(Dataset = 'IMPROVE') %>% 
  dplyr::select(Dataset, st_name, AQSCode, site_id, 
                lat = 'Latitude', long = 'Longitude', epsg, Date, Units,
    MF, RCFM, AL, AS, BR, CA, CHL, CL, CR, CU, EC, 
    FE, K, MG, MN, N2, `NA`,NI, NO3, OC, P, PB, RB, 
    S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR, ammNO3, ammSO4
  ) %>% 
  filter_at(vars(MF, RCFM, AL, AS, BR, CA, CHL, CL, CR, CU, EC, 
                 FE, K, MG, MN, N2, `NA`,NI, NO3, OC, P, PB, RB, 
                 S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR, ammNO3, ammSO4),
          any_vars(!is.na(.)))

  return(raw_IMPROVE_spec_df)
} # end function
