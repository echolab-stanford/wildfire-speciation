# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Sept 25, 2023
# Description: This script reads in downloaded CSN speciation data from: 
# http://views.cira.colostate.edu/fed/Express/ImproveData.aspx
# PM + Speciation data here: https://aqs.epa.gov/aqsweb/airdata/download_files.html


# get list of raw EPA PM2.5 files in raw data folder
# csn_pm25_list = list.files(file.path(data_fp, 'raw'), pattern = 'daily_88101', full.names = TRUE)
# csn_spec_file_list = list.files(file.path(data_fp, 'raw'), pattern = 'daily_SPEC', full.names = TRUE)
# improve_params_fp = file.path(data_fp, 'raw/IMPROVE_parameters.csv')
# csn_param_fp <- file.path(data_fp, 'raw/CSN_parameters.csv')

# FUNCTION ------------------------------------------------------------------------------------
read_in_CSN_speciation_data <- function(csn_spec_file_list, csn_pm25_list, csn_param_fp) {

# step 1 ----------------------------------------------
  # read in CSN paramater file
  csn_params <- read.csv(csn_param_fp) %>% 
    mutate(flag = ifelse(
      str_detect(Parameter, "2.5"), 'keep', 'drop'
    )) %>% 
    filter(flag == 'keep') %>% 
    mutate(flag = ifelse(
      str_detect(Parameter, "10"), 'drop', 'keep'
    )) %>% 
    filter(flag == 'keep') %>% 
    mutate(flag = ifelse(
      str_detect(Parameter, "(STP|TSP)"), 'drop', 'keep'
    )) %>% 
    filter(flag == 'keep') %>% 
    mutate(flag = ifelse(
      str_detect(Parameter, "nm"), 'drop', 'keep'
    )) %>% 
    filter(flag == 'keep') %>% 
    filter(Standard.Units == 'Micrograms/cubic meter (LC)') %>% 
    dplyr::select(AQSParamCode = 'Parameter.Code', Parameter, CAS.Number) %>% 
    filter(Parameter %in% c("Aluminum PM2.5 LC", "Ammonium Nitrate PM2.5 LC", 
                            "Ammonium Sulfate PM2.5 LC","Arsenic PM2.5 LC","Bromine PM2.5 LC", 
                            "Calcium PM2.5 LC", "Chloride PM2.5 LC", "Chlorine PM2.5 LC", "Chromium PM2.5 LC", 
                            "Chromium VI PM2.5 LC", "Copper PM2.5 LC", "EC CSN PM2.5 LC TOT",
                            "EC PM2.5 LC TOR", "EC PM2.5 LC TOT","Iron PM2.5 LC", "Lead PM2.5 LC", "Magnesium PM2.5 LC", 
                            "Manganese PM2.5 LC", "Nickel PM2.5 LC","OC CSN Unadjusted PM2.5 LC TOT", 
                            "OC PM2.5 LC TOR", "OC PM2.5 LC TOT", "OCH PM2.5 LC TOT","Organic Carbon Mass PM2.5 LC", 
                            "Phosphorus PM2.5 LC", "PM2.5 - Local Conditions",
                            "Potassium PM2.5 LC", "Reconstructed Mass PM2.5 LC",
                            "Rubidium PM2.5 LC", "Selenium PM2.5 LC", "Silicon PM2.5 LC", 
                            "Sodium PM2.5 LC", "Soil PM2.5 LC", "Strontium PM2.5 LC", "Sulfate PM2.5 LC",
                            "Sulfur PM2.5 LC",  "Titanium PM2.5 LC", "Total Nitrate PM2.5 LC", 
                            "Vanadium PM2.5 LC","Zinc PM2.5 LC","Zirconium PM2.5 LC")) %>% 
    mutate(ParamCode = case_when(
      Parameter == "Aluminum PM2.5 LC" ~ 'AL',
      Parameter == "Ammonium Nitrate PM2.5 LC" ~ 'ammNO3',
      Parameter == "Ammonium Sulfate PM2.5 LC" ~ 'ammS04',
      Parameter == "Arsenic PM2.5 LC" ~ 'AS',
      Parameter == "Bromine PM2.5 LC" ~ 'BR', 
      Parameter == "Calcium PM2.5 LC" ~ 'CA', 
      Parameter == "Chloride PM2.5 LC" ~ 'CHL',
      Parameter == "Chlorine PM2.5 LC"~ 'CL',
      Parameter == "Chromium PM2.5 LC"~ 'CR', 
      Parameter == "Chromium VI PM2.5 LC" ~ 'CR VI', 
      Parameter == "Copper PM2.5 LC" ~ 'CU', 
      Parameter == "EC CSN PM2.5 LC TOT" ~ 'EC',
      Parameter == "Iron PM2.5 LC"~ 'FE',
      Parameter == "EC PM2.5 LC TOR"~ 'EC',
      Parameter == "EC PM2.5 LC TOT"~ 'EC',
      Parameter == "Lead PM2.5 LC"~ 'PB',
      Parameter == "Magnesium PM2.5 LC"~ 'MG',
      Parameter == "Manganese PM2.5 LC" ~ 'MN',
      Parameter == "Nickel PM2.5 LC" ~ 'NI',
      Parameter == "OCH PM2.5 LC TOT"~ 'OC', 
      Parameter == "Rubidium PM2.5 LC" ~ 'RB',
      Parameter == "Zirconium PM2.5 LC"~ 'ZR',
      Parameter == "Zinc PM2.5 LC"~ 'ZN', 
      Parameter == "OC PM2.5 LC TOR"~ 'OC', 
      Parameter == "Organic Carbon Mass PM2.5 LC"~ 'OC', 
      Parameter == "OC PM2.5 LC TOT"~ 'OC', 
      Parameter == "Phosphorus PM2.5 LC"~ 'P', 
      Parameter == "PM2.5 - Local Conditions"~ 'MF', 
      Parameter == "Potassium PM2.5 LC"~ 'K',
      Parameter == "Reconstructed Mass PM2.5 LC" ~ 'RCFM',
      Parameter == "Selenium PM2.5 LC" ~ 'SE',
      Parameter == "Silicon PM2.5 LC" ~ 'SI',
      Parameter == "Sodium PM2.5 LC"~ 'NA',  
      Parameter == "Soil PM2.5 LC" ~ 'SOIL', 
      Parameter == "Strontium PM2.5 LC"~ 'SR', 
      Parameter == "Sulfate PM2.5 LC" ~ 'SO4',
      Parameter == "Titanium PM2.5 LC"~ 'TI',
      Parameter == "Sulfur PM2.5 LC"~ 'S',
      Parameter == "Vanadium PM2.5 LC" ~ 'V', 
      Parameter == "Total Nitrate PM2.5 LC" ~ 'NO3'))
    

# get a list of the parameters that are common to both sets of monitoring networks
# read in IMPROVE parameter data, and use the improve parameter 
# codes to create crosswalk between codes + parameter names in CSN data
# improve_params <- read.csv(improve_params_fp) %>% 
#   dplyr::select(ParamCode, ParamName, AQSParamCode = 'AQSCode') %>% 
#   # get rid of parameters we are not interested in
#   filter(!is.na(AQSParamCode)) %>% 
#   # there are multiple names for certain code, identify where there are duplicates and drop
#   group_by(AQSParamCode) %>% 
#   mutate(num_codes = row_number()) %>% 
#   ungroup() %>% 
#   filter(num_codes == 1) %>% 
#   dplyr::select(-num_codes)
  

# step 2 ----------------------------------------------
# testing one file at a time
# current_file <- csn_spec_file_list[3]

# read each in file, clean up, and ensure that parameters match those in the improve data
csn_spec_files = map_df(csn_spec_file_list, function(current_file) {
  
  # print current index
  print(paste0("current index is: ", match(current_file, unique(csn_spec_file_list)), ' of ', length(csn_spec_file_list)))
 
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
    filter(units == 'Micrograms/cubic meter (LC)') %>% 
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
    left_join(csn_params, by = c('Parameter', 'AQSParamCode')) %>% 
    filter(!is.na(ParamCode)) %>% 
    # drop unnecessary vars
    dplyr::select(-c(State.Code, County.Code, Site.Number, sample_dur)) %>% 
    # make sure we only have one obs per site-date
    group_by(SiteCode, Date, st_name, lat, long, Datum, units, ParamCode) %>% 
    dplyr::summarise(conc = mean(conc, na.rm = TRUE), .groups = 'drop') %>% 
    # pivot wide again so that each monitor-day has one value for each measured species
    pivot_wider(id_cols = c('SiteCode', 'Date', 'st_name',
                             'lat', 'long', 'Datum', 'units'),
                names_from = 'ParamCode',
                values_from = 'conc')

}) %>% 
  bind_rows()  


# step 3 ----------------------------------------------
# read in EPA CSN Total PM2.5 data

# test one file at a time
# current_pm_file <- csn_pm25_list[1]

# read each in PM2.5 file (88101 - fine PM2.5 - this measure is captured at both CSN and IMPROVE SITES)
csn_pm25_df = map_df(csn_pm25_list, function(current_pm_file) {
  
  # print current index
  print(paste0("current index is: ", 
               match(current_pm_file, unique(csn_pm25_list)), ' of ', length(csn_pm25_list)))
  
  # read in text files on EPA CSN chemical speciation
  pm_temp_df = read.csv(current_pm_file) %>% 
    # select variables of interest
    dplyr::select(State.Code, County.Code, Site.Number = 'Site.Num', 
                  AQSParamCode = 'Parameter.Code',
                  lat = 'Latitude', long = 'Longitude', Datum, 
                  Date = 'Date.Local', Parameter = 'Parameter.Name', 
                  MF = 'Arithmetic.Mean', st_name = 'State.Name', co_name = County.Name) %>% 
    # recreate site code to merge with later steps
    mutate(State.Code = str_remove(State.Code, "^0")) %>% 
    mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
           Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
    mutate(SiteCode = paste0("0", State.Code, County.Code, Site.Number)) %>% 
    mutate(SiteCode = case_when(
      str_length(SiteCode) <= 9 ~ SiteCode,
      str_length(SiteCode) > 9 ~ str_remove(SiteCode, "^0")
    )) %>% 
    dplyr::select(SiteCode, lat, long, Datum, Date, MF)
  
}) %>% 
  bind_rows() 

# make sure there is only one PM value per date + monitor
unique_csn_pm <- csn_pm25_df %>% 
  rename(site_id = 'SiteCode') %>% 
  group_by(site_id, Date, lat, long) %>% 
  dplyr::summarise(MF = mean(MF, na.rm = TRUE), .groups = 'drop') 


# step 4 ----------------------------------------------
# merge the PM data with the CSN raw speciation data
raw_csn_spec_df <- csn_spec_files %>% 
  mutate(site_id = SiteCode,
         Dataset = 'CSN') %>% 
  mutate(epsg = case_when(
    Datum == 'WGS84' ~ 4326,
    Datum == 'NAD83' ~ 4629)) %>% 
  dplyr::select(Dataset, st_name, lat, long, epsg, site_id, AQSCode = 'SiteCode', Date, units, RCFM, 
                AL, AS, BR, CA, CHL, CL, CR, CU, EC, FE, K, MG, MN, 
                `NA`, NI, NO3, OC, P, PB, RB, S,
                SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR) %>% 
  # add in speciation data
  left_join(unique_csn_pm) %>% 
  # filter out any rows that are entirely missing
  filter_at(vars(RCFM, AL, AS, BR, CA, CHL, CL, CR, CU, EC, FE, K, MG, MN, 
                 `NA`, NI, NO3, OC, P, PB, RB, S,
                 SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
          any_vars(!is.na(.))) 
  
return(raw_csn_spec_df)

} # end function

