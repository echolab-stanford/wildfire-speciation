# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 16, 2023
# Description: set up speciation and smoke data for analysis, final cleaning steps

# loadd(spec_w_smoke_pm_df, cache = drake_cache)

# speciation data
final_speciation_cleaning_pre_analysis <- function(spec_w_smoke_pm_df) {
  
  # step 1: clean up the speciation data
  cleaned_spec_df <- spec_w_smoke_pm_df %>% 
    # change inf to NA + NaNs to NAs
    mutate_if(is.numeric, ~ifelse(. == Inf, NA, .)) %>% 
    # change NaN to NA values
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    # if there are extreme negative values (e.g., 999 this indicates there was
    #a data entry error or the data was missing, so change to NA)
    # small negative values are fine (which can happen when a monitor needs to be calibrated)
    mutate_at(vars(MF, RCFM, smokePM2.5, AL, AS, BR, CA, CHL, CL, CR, CU,
                   EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
                   RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR
                   ), ~ifelse(. < -100, NA, .)) %>%
    # drop a row if concentrations for all chemicals are NA because we don't need that row
    filter_at(vars(MF, RCFM, smokePM2.5, AL, AS, BR, CA, CHL, CL, CR, CU,
                    EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
                    RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
              any_vars(!is.na(.))) %>%
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
    # do some cleaning of PM2.5 data, when measure for total PM is missing, 
    # use smoke PM value on that day
    mutate(MF_adj = ifelse(is.na(MF), smokePM2.5, MF),
           RCFM_adj = ifelse(is.na(RCFM), smokePM2.5, RCFM)) %>%
    # when smoke PM is greater than total PM, set smoke PM = total PM
    # if total PM is less than smoke PM, set smoke PM to be the same as total PM
    mutate(smokePM2.5 = ifelse(smokePM2.5 > MF_adj, MF_adj, smokePM2.5)) %>% # 2,088 entries where smoke was greater than measured total PM
    # calculate station smoke and non smoke using total PM var (PM2.5)
    mutate(nonsmokePM_MF = MF_adj - smokePM2.5,
           nonsmokePM_RCFM = RCFM_adj - smokePM2.5) %>% 
    # add an indicator for smoke day or nonsmoke day
    mutate(smoke_day = case_when(
      smokePM2.5 == 0 ~ 'nonsmoke day',
      smokePM2.5 != 0 ~ "smoke day")) %>% 
    # add a monitor month for FEs
    mutate(monitor_month = paste0(site_id,"_",month)) %>% 
    dplyr::select(Dataset, state_name, region, season, year, month, doy, 
                  duration, st_date, end_date, monitor_month, Date, site_id, smoke_day, 
                  MF_adj, RCFM_adj, smokePM = 'smokePM2.5', nonsmokePM_MF, nonsmokePM_RCFM,
                  AL, AS, BR, CA, CL, CR, CU,
                  EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
                  RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR, units, 
                  long, lat, epsg, MF, RCFM)
  
  return(cleaned_spec_df)       
} # end function


  