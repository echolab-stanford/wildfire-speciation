# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 16, 2023
# Description: set up speciation and smoke data for descriptive statistics

# loadd(pm_plume_speciation_df, cache = drake_cache)
# stationPM_fp = file.path(wip_gdrive_fp, 'intermediate/pm_plume_speciation_at_sites_monitorPM.csv')


calculate_and_clean_spec_categories <- function(pm_plume_speciation_df) {
  
  # step 1: clean up the speciation data
  cleaned_spec_df <- pm_plume_speciation_df %>% 
    # what does 0 concentration actually mean in these cases?
    # mutate_if(is.numeric, ~ifelse(. == 0, NA, .)) %>% 
    # drop a row if concentrations for all chemicals are NA
    filter(if_any(c(AL, AS, BR, CA, EC, OC, CL,PM2.5, 
                    CR, CU, FE, PB, MG, MN, OP,
                    MO, NI, RCFM, NO3, P, K, RB,
                    `NA`, SE, SI, S, SR, SOIL, SO4, TI, TC,
                    V, ZN, ZR, ammNO3, ammSO4, CHL, N2), 
                  ~!is.na(.))) %>%
    mutate(NH4 = ((18/(18+62)*ammNO3)) + (18/((18+96)*ammSO4))) %>%
    # calculate metals category
    mutate(tot_metals = AL + AS + CA + CR + CU + FE + PB + NI + MG + 
             MN + `NA` + TI + K + RB + SR + V + ZN + ZR) %>% # dropping MO bc lots of missing 
    mutate(year = year(Date),
           month = month(Date),
           doy = yday(Date)) %>% 
    # drop anything prior to 2005 because we are limited by PM2.5 predictions (childs et al, 2022)
    filter(year > 2005) %>% 
    # create seasons 
    mutate(season = case_when(
      month %in% c(12,1,2) ~ 'winter',
      month %in% c(3,4,5) ~ 'spring',
      month %in% c(6,7,8) ~ 'summer',
      month %in% c(9,10,11) ~ 'autumn',
    )) %>% 
    # create regions
    mutate(region = case_when(
      State %in% c('OR', 'CA', 'WA') ~ 'pacific',
      State %in% c('NV', 'UT', 'CO', 'WY', 'MT', 'ID') ~ 'rocky_mountain',
      State %in% c('AZ', 'NM', 'TX', 'OK') ~ 'southwest',
      State %in% c('ND', 'SD', 'NE', 'KS', 'MN', 'IA', 'MO', 'WI', 'IL', 'IN', 'MI', 'OH') ~ 'midwest',
      State %in% c('KY', 'WV', 'VA', 'TN', 'NC', 'MS', 'AL', 'GA', 'SC', 'FL', 'LA', 'AR') ~ 'southeast',
      State %in% c('ME', 'VT', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'PA', 'DE', 'MD', 'DC') ~ 'northeast',
      State %in% c('HI', 'AK', 'VI') ~ 'noncontiguous',
      is.na(State)  ~ 'outside US',
      State %in% c('AB', 'ON') ~ 'outside US'
    )) %>% 
    # Exclude regions that are  'noncontiguous' and 'outside US'
    filter(!region %in% c('noncontiguous', 'outside US')) %>% 
    # add in units column for each species, all species are same units
    mutate(units = 'ug_m3') %>%
    # ensure unique observations
    distinct() %>% 
    # change inf to NA
    mutate_if(is.numeric, ~ifelse(. == Inf, NA, .)) %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    # calculate station smoke and non smoke using total PM var (PM2.5)
    mutate(non_smokePM_cont = ifelse(!is.na(PM2.5), PM2.5 - smokePM_pred, NA)) %>% 
    dplyr::select(Dataset, SiteCode, Date, State, region, season, year, month, doy, 
                  km2fire, smoke_day, PM2.5, smokePM_pred,  non_smokePM_cont, RC_PM2 = 'RCFM', 
                  NH4, tot_metals, ammNO3, ammSO4, AL, 
                  AS, BR, CA, EC, OC, TC, OP, CHL, CL, CR, 
                  CU, FE, K, MG, MN,
                  MO, `NA`, NI, NO3, N2, P, PB, RB, 
                  S, SE, SI, SOIL, SO4, SR, V, ZN, TI,
                  ZR, units, Elevation, Latitude, Longitude, epsg) %>% 
    mutate(monitor_month = paste0(SiteCode,"_",month)) %>% 
    rename(totPM2.5 = 'PM2.5', 
           smokePM = 'smokePM_pred',
           nonsmokePM = 'non_smokePM_cont') %>% 
    mutate(pm_flag = case_when(
      totPM2.5 < smokePM ~'drop', 
      TRUE ~ 'keep')) %>% 
    filter(pm_flag == 'keep')
                  
  return(cleaned_spec_df)       
} # end function


# station data was calculated using old MF col, need to update so commenting out code
 
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
  



# OLD CODE
  # group_by(SiteCode) %>% 
  # mutate(st_date = year(min(Date)),
  #        end_date = year(max(Date))) %>% 
  # ungroup() %>% 
  # # get duration that a site is online
  # mutate(duration = (end_date - st_date) + 1) %>% 

# calculate fractions
# mutate(OC_frac = ifelse(MF != 0, OC/MF, NA),
#        EC_frac = ifelse(MF != 0, EC/MF, NA), 
#        NO3_frac = ifelse(MF != 0, NO3/MF, NA), 
#        dust_frac = ifelse(MF != 0, SOIL/MF, NA), 
#        SO4_frac = ifelse(MF != 0, SO4/MF, NA), 
#        NH4_frac = ifelse(MF != 0, NH4/MF, NA), 
#        metals_frac = ifelse(MF != 0, tot_metals/MF, NA)) %>%

  