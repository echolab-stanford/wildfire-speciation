# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 16, 2023
# Description: read in total PM for CSN sites

# loadd(cleaned_spec_df, cache = drake_cache)

read_in_csn_totalPM_merge_w_spec <- function(cleaned_spec_df) {
  
  # list all the files speifically containing the measure of total PM for CSN sites
  csn_pm_file_list = list.files(file.path(wip_gdrive_fp, '/raw'), 
                                pattern = 'daily_88502', full.names = TRUE)
  
  # read in all files
  # current_file <- csn_pm_file_list[1] #test
  all_csn_pm_files <- purrr::map_df(csn_pm_file_list, function(current_file) {
    
    current_pm_file <- read.csv(current_file) %>% 
      # select variables of interest
      dplyr::select(State.Code, County.Code, Site.Number = Site.Num,
                    lat = Latitude, long =Longitude, Datum, Date = Date.Local, param = Parameter.Name, 
                    sampl_dur = Sample.Duration, obs_count = Observation.Count, 
                    obs_pct = Observation.Percent, PMconc = Arithmetic.Mean) %>% 
      # recreate site code to merge with later steps
      mutate(State.Code = str_remove(State.Code, "^0")) %>% 
      mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
             Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
      mutate(SiteCode = paste0("0", State.Code, County.Code, Site.Number)) %>% 
      mutate(SiteCode = case_when(
        str_length(SiteCode) <= 9 ~ SiteCode,
        str_length(SiteCode) > 9 ~ str_remove(SiteCode, "^0")
      )) %>% 
      distinct(SiteCode, PMconc, lat, long, Datum, Date, 
               param, sampl_dur, obs_count, obs_pct) 
    
    # split into three datasets based on the type of measurement taken, and use data in a hierarchical way
    # then split into the 24 hour avg
    df24hr <- current_pm_file %>% 
      filter(sampl_dur == '24 HOUR') %>% 
      group_by(SiteCode, lat, long, Datum, Date, param, sampl_dur) %>% 
      dplyr::summarise(PM = mean(PMconc, na.rm = TRUE)) %>% 
      ungroup()
    
    # split into the 24 hour block avg
    df24hr_block <- current_pm_file %>% 
      filter(sampl_dur == '24-HR BLK AVG') %>% 
      group_by(SiteCode, lat, long, Datum, Date, param, sampl_dur) %>% 
      dplyr::summarise(PM = mean(PMconc, na.rm = TRUE)) %>% 
      ungroup()
    
    # then split into the 1 hour avg
    df1hr <- current_pm_file %>% 
      filter(sampl_dur == '1 HOUR') %>% 
      group_by(SiteCode, lat, long, Datum, Date, param, sampl_dur) %>% 
      dplyr::summarise(PM = mean(PMconc, na.rm = TRUE)) %>% 
      ungroup()
    
    # check how many unique sites are in the 24 hour block that are not in the 24 hour avg
    check24 <- df24hr_block %>% 
      filter(!SiteCode %in% df24hr$SiteCode) %>% 
      distinct(SiteCode)
    print(paste0(nrow(check24), ' unique site(s) are in the 24 hour block avg that are not in the 24 hour avg'))
    
    
    # USE A HIERARCHY OF DATA TO KEEP,
    # first bind together the 24 hour data with 24 hour block
    df24hr_plus24hrblock <- df24hr %>% 
      bind_rows(df24hr_block %>% 
                  filter(!SiteCode %in% df24hr$SiteCode)) 
    
    # now determine what new sites are in the 1-hour dataframe 
    # that are not already included in the 24 hr + 24 hr block data
    
    # check how many unique sites are in the 1 hour block that are not in the 24 hour df
    check1 <- df1hr %>% 
      filter(!SiteCode %in% df24hr_plus24hrblock$SiteCode) %>% 
      distinct(SiteCode)
    print(paste0(nrow(check1), ' unique site(s) are in the 1 hour that is/are not in the 24 hour data'))
    
    # now bind together the 1 hour data with 24 hour
    all_current_pm_df <- df24hr_plus24hrblock %>% 
      bind_rows(df1hr %>% 
                  filter(!SiteCode %in% df24hr_plus24hrblock$SiteCode)) 
    
  }) %>% 
    bind_rows() %>% 
    dplyr::select(SiteCode, Date, sampl_dur, PM2.5 = 'PM') %>% 
    distinct()
  
  # merge csn pm with speciation data
  all_spec_df <- all_csn_pm_files %>% 
    left_join()
  
  
  return()
}