# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 16, 2023
# Description: set up speciation and smoke data for descriptive statistics

# loadd(spec_smoke_plumes_df, cache = drake_cache)


final_speciation_cleaning_pre_analysis <- function(spec_smoke_plumes_df) {
  
  # step 1: clean up the speciation data
  cleaned_spec_df <- spec_w_smoke_pm_df %>% 
    # change inf to NA + NaNs to NAs
    mutate_if(is.numeric, ~ifelse(. == Inf, NA, .)) %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    mutate_if(is.numeric, ~ifelse(. == -999, NA, .)) %>%
    mutate_if(is.numeric, ~ifelse(. == -999.00000, NA, .)) %>%
    # drop a row if concentrations for all chemicals are NA
    filter_at(vars(MF, RCFM, smokePM2.5, AL, AS, BR, CA, CHL, CL, CR, CU,
                    EC, FE, K, MG, MN, `NA`, NI, NO3, N2, OC, P, PB,
                    RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
              any_vars(!is.na(.))) %>%
    # change any negative value in the data to NA
    mutate_at(vars(MF, RCFM, smokePM2.5, AL, AS, BR, CA, CHL, CL, CR, CU,
                   EC, FE, K, MG, MN, `NA`, NI, NO3, N2, OC, P, PB,
                   RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
              ~ifelse(. < 0, NA, .)) %>%
    mutate(year = year(Date),
           month = month(Date),
           doy = yday(Date)) %>% 
    # drop anything prior to 2005 because we are limited by PM2.5 predictions (childs et al, 2022)
    filter(year > 2005 & year < 2021) %>%
    group_by(site_id) %>%
    mutate(st_date = year(min(Date)),
           end_date = year(max(Date))) %>%
    ungroup() %>%
    # get duration that a site is online
    mutate(duration = (end_date - st_date) + 1) %>%
    # create seasons 
    mutate(season = case_when(
      month %in% c(12,1,2) ~ 'winter',
      month %in% c(3,4,5) ~ 'spring',
      month %in% c(6,7,8) ~ 'summer',
      month %in% c(9,10,11) ~ 'autumn',
    )) %>% 
    # create regions
    mutate(region = case_when(
      state_name %in% c('Oregon', 'California', 'Washington') ~ 'pacific',
      state_name %in% c('Nevada', 'Utah', 'Colorado', 'Wyoming', 'Montana', 'Idaho') ~ 'rocky_mountain',
      state_name %in% c('Arizona', 'New Mexico', 'Texas', 'Oklahoma') ~ 'southwest',
      state_name %in% c('North Dakota', 'South Dakota', 'Nebraska', 'Kansas', 'Minnesota', 'Iowa', 
                        'Missouri', 'Wisconsin', 'Illinois', 'Indiana', 'Michigan', 'Ohio') ~ 'midwest',
      state_name %in% c('Kentucky', 'West Virginia', 'Virginia', 'Tennessee', 'North Carolina', 'Mississippi', 
                        'Alabama', 'Georgia', 'South Carolina', 'Florida', 'Louisiana', 'Arkansas') ~ 'southeast',
      state_name %in% c('Maine', 'Vermont', 'New Hampshire', 'Massachusetts', 'Rhode Island', 'Connecticut', 
                        'New York', 'New Jersey', 'Pennsylvania', 'Delaware', 'Maryland', 'District of Columbia') ~ 'northeast',
    )) %>% 
    # add in units column for each species, all species are same units
    mutate(units = 'ug_m3') %>%
    # ensure unique observations
    distinct() %>% 
    # do some cleaning of PM2.5 data, when measure for total PM is missing, set = to smoke PM
    mutate(MF_adj = ifelse(is.na(MF), smokePM2.5, MF),
           RCFM_adj = ifelse(is.na(RCFM), smokePM2.5, RCFM)) %>% 
    # now drop if all speciation vars are NA for a row
    # drop a row if concentrations for all chemicals are NA
    filter_at(vars(AL, AS, BR, CA, CHL, CL, CR, CU,
                   EC, FE, K, MG, MN, `NA`, NI, NO3, N2, OC, P, PB,
                   RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
              any_vars(!is.na(.))) %>% 
    # calculate station smoke and non smoke using total PM var (PM2.5)
    mutate(nonsmokePM_MF = MF_adj - smokePM2.5,
           nonsmokePM_RCFM = RCFM_adj - smokePM2.5) %>% 
    dplyr::select(Dataset, state_name, region, season, year, month, doy, 
                  duration, st_date, end_date, Date, site_id, 
                  MF_adj, RCFM_adj, smokePM = 'smokePM2.5', nonsmokePM_MF, nonsmokePM_RCFM,
                  AL, AS, BR, CA, CHL, CL, CR, CU,
                  EC, FE, K, MG, MN, `NA`, NI, NO3, N2, OC, P, PB,
                  RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR, units, 
                  long, lat, epsg, MF, RCFM) 
  
  return(cleaned_spec_df)       
} # end function


# OLD CODE
# station data was calculated using old MF col, need to update so commenting out code
# stationPM_fp = file.path(wip_gdrive_fp, 'intermediate/pm_plume_speciation_at_sites_monitorPM.csv') 

# # read in calculated station non smoke and smoke data from Ayako
# # this calculated station data will serve as a check
# station_calc_df <- read.csv(stationPM_fp) %>% 
#   dplyr::select(SiteCode, Date, smokePM_MF) %>% 
#   group_by(SiteCode, Date) %>% 
#   dplyr::summarise(station_calc_smokePM = mean(smokePM_MF, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   distinct() %>% 
#   mutate(Date = as.Date(Date, format = '%Y-%m-%d'))

# set up regression dataframe
# adj_pm_spec_df <- clean_pm_spec_df %>% 
#   left_join(station_calc_df, by = c('SiteCode', 'Date')) %>% 
#   
# station_non_smokePM = ifelse(!is.na(PM2.5), PM2.5 - station_calc_smokePM, NA))  

# mutate(NH4 = ifelse(!is.na(ammNO3) & !is.na(ammSO4),
#   ((18/(18+62)*ammNO3)) + (18/((18+96)*ammSO4)), 
#   NA)) %>%
# # calculate metals category
# mutate(tot_metals = AL + AS + CA + CR + CU + FE + PB + NI + MG + 
#          MN + `NA` + TI + K + RB + SR + V + ZN + ZR) %>% # dropping MO bc lots of mis

# calculate fractions
# mutate(OC_frac = ifelse(MF != 0, OC/MF, NA),
#        EC_frac = ifelse(MF != 0, EC/MF, NA), 
#        NO3_frac = ifelse(MF != 0, NO3/MF, NA), 
#        dust_frac = ifelse(MF != 0, SOIL/MF, NA), 
#        SO4_frac = ifelse(MF != 0, SO4/MF, NA), 
#        NH4_frac = ifelse(MF != 0, NH4/MF, NA), 
#        metals_frac = ifelse(MF != 0, tot_metals/MF, NA)) %>%

  