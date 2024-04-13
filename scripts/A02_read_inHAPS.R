# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Sept 25, 2023
# Description: This script reads in downloaded speciation data on HAPS from
# https://aqs.epa.gov/aqsweb/airdata/download_files.html


# get list of raw EPA PM2.5 files in raw data folder
# csn_pm25_list = list.files(file.path(data_fp, 'raw'), pattern = 'daily_88101', full.names = TRUE)
# csn_spec_file_list = list.files(file.path(data_fp, 'raw'), pattern = 'daily_SPEC', full.names = TRUE)
# improve_params_fp = file.path(data_fp, 'raw/IMPROVE_parameters.csv')
# csn_param_fp <- file.path(data_fp, 'raw/CSN_parameters.csv')

# FUNCTION ------------------------------------------------------------------------------------
read_in_HAPS_speciation_data <- function(csn_spec_file_list, csn_pm25_list, csn_param_fp) {
 
  # step 1 ----------------------------------------------
  # testing one file at a time
  haps_spec_file_list = list.files(file.path(data_fp, "raw"), pattern = "daily_HAPS", full.names= TRUE)
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

  return(haps_spec_df)
  
} # end function

