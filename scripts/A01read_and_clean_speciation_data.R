# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Dec 9, 2022
# Description: This script reads in downloaded speciation data from: 
# http://views.cira.colostate.edu/fed/Express/ImproveData.aspx


# improve_spec_fp = file.path(wip_gdrive_fp, 'intermediate/IMPROVE_Aerosol_combined.csv')
# loadd(c(csn_spec_files), cache = drake::new_cache("scripts/.drake"))


# function
read_and_clean_speciation_data <- function(improve_spec_fp, csn_spec_files) {
  
    # read in improve speciation combined file, clean up by dropping variables that are not of interest
    improve_spec_df <- read.csv(improve_spec_fp) %>%
      dplyr::select(-POC, -contains("Ref"), -contains("Trans"), -contains('FlowRate'),
                    -contains('Unit'),-contains('OP'),
                    -contains('SampDur'), -contains('SeaSalt'), -contains('fAbs')) %>% 
      rename_at(vars(matches("_Val")), ~str_remove(., "f_Val")) %>% 
      rename_at(vars(matches("_Val")), ~str_remove(., "_Val")) %>% 
      mutate(Date = as.Date(Date, format = '%m/%d/%Y'),
             EPACode = as.character(EPACode)) %>% 
      filter(!is.na(Date)) %>% 
      mutate_if(is.numeric, ~ifelse(. == -999, NA, .)) %>% 
      mutate(Dataset = 'IMPROVE') %>% 
      dplyr::select(Dataset, SiteCode, EPACode, Date, State, Latitude, Longitude, Elevation,
                    AL, ammNO3, ammSO4, AS, BR, EC, OC, TC, CA, CHL, 
                    CL, CR, CU, FE, NI, PB, MG, MN, PM2.5 = 'MF', RCFM, 
                    `NA`, NO3, N2, P, K, RB, SE, SI, SR, SOIL, 
                    S, SO4, TI, V, ZN, ZR) %>% 
      distinct() 
    

    # clean up CSN speciation file
    csn_spec_df <- csn_spec_files %>% 
      dplyr::select(-c(contains("881"), contains("882"), contains("883"),
                       contains("68"), contains("42"))) %>% 
      rename_at(vars(matches("_Val")), ~str_remove(., "f_Val")) %>% 
      rename_at(vars(matches("_Val")), ~str_remove(., "_Val")) %>% 
      mutate(Date = as.Date(Date, format = '%m/%d/%Y')) %>% 
      filter(!is.na(Date)) %>% 
      dplyr::select(Dataset, SiteCode, EPACode, Date, State, PM2.5 = `88502`,
                    AL, AS, BR, CA, EC, OC, OP, CL, CR, CU, FE, PB, MG, MN, MO, NI, 
                    RCFM, NO3, P, K, RB, `NA`, SE, SI, S, SR, SOIL, SO4, TI, V, ZN, ZR) %>% 
      distinct() %>% 
      # get the avg concentration of a species at a given monitor-day when there are more than one observation per day
      pivot_longer(cols = c('PM2.5', 'RCFM', 'AL', 'AS', 'BR', 
                            'CA', 'EC', 'OC', 'OP', 'CL',
                            'CR', 'CU', 'FE', 'K', 'MG', 'MN',
                            'MO', 'NA', 'NI', 'NO3', 'P', 
                            'SO4', 'SR', 'TI', 'V', 'ZN', 'ZR',
                            'PB', 'RB', 'S', 'SE', 'SI', 'SOIL'), 
                   names_to = 'lhs', values_to = 'conc') %>%   
      group_by(Dataset, SiteCode, EPACode, Date, State, lhs) %>% 
      # get avg conc of a species on a given date
      dplyr::summarise(conc = mean(conc, na.rm = TRUE)) %>% 
      ungroup() %>% 
      #mutate(conc = ifelse(is.NaN(conc), NA, conc))# %>% 
      # pivot wide again so that each monitor-day has one value for each measured species
      pivot_wider(id_cols = c('Dataset', 'SiteCode', 'EPACode', 'Date', 'State'), 
                  names_from = 'lhs', 
                  values_from = 'conc') 
      

    # now read in the additional measures of PM for CSN files 
    # (this measure should more closely match the PM measure in IMPROVE files)
    # list all the files specifically containing the measure of total PM for CSN sites
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
                      obs_pct = Observation.Percent, PM2.5 = Arithmetic.Mean) %>% 
        # recreate site code to merge with later steps
        mutate(State.Code = str_remove(State.Code, "^0")) %>% 
        mutate(County.Code = str_pad(County.Code, 3, pad = "0"), # #
               Site.Number = str_pad(Site.Number, 4, pad = "0")) %>% 
        mutate(SiteCode = paste0("0", State.Code, County.Code, Site.Number)) %>% 
        mutate(SiteCode = case_when(
          str_length(SiteCode) <= 9 ~ SiteCode,
          str_length(SiteCode) > 9 ~ str_remove(SiteCode, "^0")
        )) %>% 
        distinct(SiteCode, PM2.5, lat, long, Datum, Date, 
                 param, sampl_dur, obs_count, obs_pct) 
      
      # split into three datasets based on the type of measurement taken, and use data in a hierarchical way
      # then split into the 24 hour avg
      df24hr <- current_pm_file %>% 
        filter(sampl_dur == '24 HOUR') %>% 
        group_by(SiteCode, lat, long, Datum, Date, param, sampl_dur) %>% 
        dplyr::summarise(PM2.5 = mean(PM2.5, na.rm = TRUE)) %>% 
        ungroup()
      
      # split into the 24 hour block avg
      df24hr_block <- current_pm_file %>% 
        filter(sampl_dur == '24-HR BLK AVG') %>% 
        group_by(SiteCode, lat, long, Datum, Date, param, sampl_dur) %>% 
        dplyr::summarise(PM2.5 = mean(PM2.5, na.rm = TRUE)) %>% 
        ungroup()
      
      # then split into the 1 hour avg
      df1hr <- current_pm_file %>% 
        filter(sampl_dur == '1 HOUR') %>% 
        group_by(SiteCode, lat, long, Datum, Date, param, sampl_dur) %>% 
        dplyr::summarise(PM2.5 = mean(PM2.5, na.rm = TRUE)) %>% 
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
      dplyr::select(SiteCode, Date, PM2.5) %>% 
      mutate(Date = as.Date(Date, format = '%Y-%m-%d')) %>% 
      distinct()
    
    # merge csn pm with speciation data
    all_csn_spec_long <- csn_spec_df %>% 
      dplyr::select(SiteCode, Date, PM2.5) %>% 
      group_by(SiteCode, Date) %>% 
      dplyr::summarise(PM2.5 = mean(PM2.5, na.rm = TRUE), .groups = 'drop') %>% 
      mutate(source = 'csn_spec_files') %>% 
      # now add more rows for the additional PM data for CSN
      bind_rows(all_csn_pm_files %>% 
                  mutate(source = 'daily_files')) %>% 
      # pivot wide again
      pivot_wider(names_from = 'source', 
                  values_from = 'PM2.5',
                  values_fill = NA) %>% 
      # create flag that indicates if the PM values from the daily files match the PM
      # values in the speciation files
      mutate(diff_flag = case_when(
        csn_spec_files == daily_files ~ 'match',
        TRUE ~ 'mismatch'
      )) %>% 
      mutate(csn_spec_files = ifelse(is.nan(csn_spec_files), NA, csn_spec_files))
                  
    # separate into cases that have matching data and not matching data for values of PM
    # on a given monitoring date
    
    # matching values for PM (keep this as is)
    matching_csn_pm <- all_csn_spec_long %>% 
      filter(diff_flag == 'match') %>% 
      dplyr::select(SiteCode, Date, PM2.5 = daily_files)
    
    # mismatch numbers
    mismatch_csn_pm <- all_csn_spec_long %>% 
      filter(diff_flag == 'mismatch') %>% 
      mutate(keep_flag = case_when(
        is.na(daily_files) ~ 'keep_csn',
        is.na(csn_spec_files) ~ 'keep_daily',
        daily_files != csn_spec_files ~ 'drop'
      )) %>% 
      filter(keep_flag == 'keep_csn' |  keep_flag == 'keep_daily') %>% 
      # drop the flags
      dplyr::select(-c(diff_flag, keep_flag)) %>% 
      # pivot longer again, keeping only the values that are not NA on a given date
      pivot_longer(cols = c(csn_spec_files, daily_files),
                   names_to = 'file_type', values_to = 'PM2.5') %>% 
      filter(!is.na(PM2.5)) %>% 
      dplyr::select(-file_type)
     
    # bind together all PM data for CSN
    all_cleaned_pm_csn_df <- bind_rows(matching_csn_pm, mismatch_csn_pm) %>% 
      mutate(Dataset = 'EPACSN') %>% 
      # make sure each site only has one value for PM for a given monitor-day
      group_by(SiteCode, Date, Dataset) %>%
      dplyr::summarise(PM2.5 = mean(PM2.5, na.rm = TRUE)) %>%
      ungroup() %>% 
      # only keep PM data where we have other chemical speciation data on the same date
      filter(SiteCode %in% csn_spec_df$SiteCode)
    
    # join the rest of the speciation data by monitoring site and date but only for those sites that 
    # measure other speciation data
    all_joined_csn_spec_df <- csn_spec_df %>% 
      dplyr::select(-PM2.5) %>% 
      left_join(all_cleaned_pm_csn_df,
                  by = c('SiteCode', 'Date', 'Dataset')) %>% 
      # make sure NaNs are converted to NA
      mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .))
    
    # --------------------------------------------------------------------------------
    # merge speciation data together:
    all_spec_df <- bind_rows(all_joined_csn_spec_df, improve_spec_df) %>% 
      filter_at(vars(AL, AS, BR, CA, EC, OC, OP, TC, CHL, CL, CR, 
                     CU, FE, PB, MG, MN, MO, NI, 
                     RCFM, NO3, P, K, RB, `NA`, SE, SI, S, SR, SOIL, SO4, 
                     TI, N2, V, ZN, ZR, ammNO3, ammSO4, PM2.5),
                any_vars(!is.na(.))) %>% 
      distinct() 
    
  return(all_spec_df)
  
} # end function
