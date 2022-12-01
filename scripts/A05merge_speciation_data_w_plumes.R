# Description: merge speciation data with plume data

# loadd(c(cleaned_spec_df, improve_smoke_dens_fire_dist), cache = drake_cache)

merge_speciation_and_plumes <- function(cleaned_spec_df, improve_smoke_dens_fire_dist) {
  
  
  merged_spec_plumes <- improve_smoke_dens_fire_dist %>% 
    left_join(cleaned_spec_df,
              by = c('Dataset', 'SiteCode', 'Date', 'epsg', 'Elevation', 'Longitude', 'Latitude')) 
  
  
  # final cleaning + add relevant vars, e.g., regions, dates, etc
  all_plume_speciation_df <- merged_spec_plumes %>% 
    filter_at(vars(AL, AS, BR, CA, EC1, EC2, EC3, EC, OC1, OC2, OC3, 
                   OC4, OC, OP, TC, CHL, CL, CR, CU, FE, PB, MG, MN, MT, MF, MO, NI, 
                   RCTM, RCFM, NO3, P, K, RB, `NA`, SE, SI, S, SR, SOIL, SO4, 
                   TI, N2, V, ZN, ZR, CM_calculated, ammNO3, ammSO4),
              any_vars(!is.na(.))) %>% 
    group_by(SiteCode) %>% 
    mutate(st_date = year(min(Date)),
           end_date = year(max(Date))) %>% 
    ungroup() %>% 
    # get duration that a site is online
    mutate(duration = (end_date - st_date) + 1) %>% 
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
        State %in% c('HI', 'AK', 'VI') ~ 'noncontiguous'
      )) %>% 
    mutate(units = 'ug_m3') %>% 
    distinct() %>% 
    dplyr::select(Dataset, SiteCode, Date, st_date, end_date, duration, 
                  year, month, doy, season, region, State, Latitude, Longitude, 
                  epsg, Elevation, LandUseCode, smoke_day, low_count, med_count, 
                  high_count, density_missing, km2fire, AL, AS, BR, CA, EC1, EC2, EC3,
                  EC, OC1,OC2, OC3, OC4, OC, OP, CL, CR, CU, FE, PB, MG, MN, MF, MO, NI,
                  RCFM, NO3, P, K, RB, `NA`, SE, SI, S, SR, SOIL, SO4, TI, V, ZN, ZR, 
                  ammNO3, ammSO4, TC, CHL, MT, RCTM, N2, CM_calculated, units)
    
  
  return(all_plume_speciation_df)
  
}
