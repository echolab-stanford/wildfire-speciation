# Description: set up speciation and smoke data for descriptive statistics

# loadd(pm_plume_speciation_df, cache = drake_cache)


calculate_and_clean_spec_categories <- function(pm_plume_speciation_df) {
  
  # step 1: adjust negatives by adding the absolute value of the lowest conc of each compound to every value
  # get all the lowest negative values for each compound and add by that value to every value
  # min_vals <- pm_plume_speciation_df %>% 
  #   dplyr::select(-c(EC1, EC2, EC3, OC1, OC2, OC3, OC4,
  #      OP, RCTM, CM_calculated, MT, TC)) %>% 
  #   dplyr::select(AL:N2) %>% 
  #   mutate_all(funs(adj = . + (-1)*min(., na.rm = TRUE))) %>% 
  #   dplyr::select(ends_with("adj")) %>% 
  #   mutate(row = row_number()) 

  # now clean up the speciation data
  cleaned_spec_df <- pm_plume_speciation_df %>% 
    dplyr::select(-c(EC1, EC2, EC3, OC1, OC2, OC3, OC4,
                     OP, RCTM, CM_calculated, MT, TC)) %>% 
    # mutate(row = row_number()) %>% 
    # left_join(min_vals, by = 'row') %>%
    #dplyr::select(-c(AL:N2, row)) %>% 
    mutate_if(is.numeric, ~ifelse(. == 0, NA, .)) %>% # if concentration is 0,
    filter(if_any(c(AL, AS, BR, CA, EC, OC, CL,
                    CR, CU, FE, PB, MG, MN, 
                    MF, MO, NI, RCFM, NO3, P, K, RB,
                    `NA`, SE, SI, S, SR, SOIL, SO4, TI, 
                    V, ZN, ZR, ammNO3, ammSO4, CHL, N2), 
                  ~!is.na(.))) %>%
    mutate(NH4 = ((18/(18+62)*ammNO3)) + (18/((18+96)*ammSO4))) %>%
    # calculate metals category
    mutate(tot_metals = AL + AS + CA + CR + CU + FE + PB + NI + MG + 
             MN + `NA` + TI + K + RB + SR + V + ZN + ZR) %>% # dropping MO bc lots of missing 
    mutate(year = year(Date),
           month = month(Date),
           doy = yday(Date)) %>% 
    mutate(season = case_when(
      month %in% c(12,1,2) ~ 'winter',
      month %in% c(3,4,5) ~ 'spring',
      month %in% c(6,7,8) ~ 'summer',
      month %in% c(9,10,11) ~ 'autumn',
    )) %>% 
    mutate(region = case_when(
      State %in% c('OR', 'CA', 'WA') ~ 'pacific',
      State %in% c('NV', 'UT', 'CO', 'WY', 'MT', 'ID') ~ 'rocky_mountain',
      State %in% c('AZ', 'NM', 'TX', 'OK') ~ 'southwest',
      State %in% c('ND', 'SD', 'NE', 'KS', 'MN', 'IA', 'MO', 'WI', 'IL', 'IN', 'MI', 'OH') ~ 'midwest',
      State %in% c('KY', 'WV', 'VA', 'TN', 'NC', 'MS', 'AL', 'GA', 'SC', 'FL', 'LA', 'AR') ~ 'southeast',
      State %in% c('ME', 'VT', 'NH', 'MA', 'RI', 'CT', 'NY', 'NJ', 'PA', 'DE', 'MD', 'DC') ~ 'northeast',
      State %in% c('HI', 'AK', 'VI') ~ 'noncontiguous',
      is.na(State) | State %in% c('AB', 'ON') ~ 'outside US'
    )) %>% 
    # add in units column for each species, all species are same units
    mutate(units = 'ug_m3')  %>%
    # mutate(OC_frac = ifelse(MF != 0, OC/MF, NA),
    #        EC_frac = ifelse(MF != 0, EC/MF, NA), 
    #        NO3_frac = ifelse(MF != 0, NO3/MF, NA), 
    #        dust_frac = ifelse(MF != 0, SOIL/MF, NA), 
    #        SO4_frac = ifelse(MF != 0, SO4/MF, NA), 
    #        NH4_frac = ifelse(MF != 0, NH4/MF, NA), 
    #        metals_frac = ifelse(MF != 0, tot_metals/MF, NA)) %>%
    dplyr::select(Dataset, SiteCode, Date, State, region, season, year, month, doy, 
                  km2fire, smoke_day, low_count, med_count, high_count, density_missing, 
                  smokePM_pred, MF, RC_PM2 = 'RCFM', 
                  # OC_frac, EC_frac, NO3_frac, dust_frac,
                  # SO4_frac, NH4_frac, metals_frac,
                  NH4, tot_metals, ammNO3, ammSO4, AL, 
                  AS, BR, CA, EC, OC, CHL, CL, CR, 
                  CU, FE, K, MG, MN,
                  MO, `NA`, NI, NO3, N2, P, PB, RB, 
                  S, SE, SI, SOIL, SO4, SR, V, ZN, 
                  ZR, units, Elevation, Latitude, Longitude, epsg,
                  ammNO3, ammSO4, AL, 
                  AS, BR, CA, EC, OC, CHL, CL, CR, 
                  CU, FE, K, MG, MN,
                  MO, `NA`, NI, NO3, N2, P, PB, RB, 
                  S, SE, SI, SOIL, SO4, SR, V, ZN, ZR) %>% 
    distinct() %>% 
    filter(year > 2005) %>% 
    # Exclude region =  'noncontiguous' and 'outside US'
    filter(region != 'noncontiguous' & region != 'outside US')

  return(cleaned_spec_df)       
} # end function

 

  
  # group_by(SiteCode) %>% 
  # mutate(st_date = year(min(Date)),
  #        end_date = year(max(Date))) %>% 
  # ungroup() %>% 
  # # get duration that a site is online
  # mutate(duration = (end_date - st_date) + 1) %>% 

  