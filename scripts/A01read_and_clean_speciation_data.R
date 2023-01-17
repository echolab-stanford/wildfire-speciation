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
                    AL, ammNO3, ammSO4, AS, BR, EC1, EC2, EC3, EC, 
                    OC1, OC2, OC3, OC4, OC, TC, CA, CHL, 
                    CL, CR, CU, FE, NI, PB, MG, MN, MT, MF, RCTM, RCFM, 
                    `NA`, NO3, N2, P, K, RB, SE, SI, SR, SOIL, 
                    S, SO4, TI, V, ZN, ZR, CM_calculated) %>% 
      distinct() 
    

    # clean up CSN speciation file
    csn_spec_df <- csn_spec_files %>% 
      dplyr::select(-c(contains("88"), contains("68"), contains("42"))) %>% 
      rename_at(vars(matches("_Val")), ~str_remove(., "f_Val")) %>% 
      rename_at(vars(matches("_Val")), ~str_remove(., "_Val")) %>% 
      mutate(Date = as.Date(Date, format = '%m/%d/%Y')) %>% 
      filter(!is.na(Date)) %>% 
      dplyr::select(Dataset, SiteCode, EPACode, Date, State, 
                    AL, AS, BR, CA, EC1, EC2, EC3, EC, OC1, OC2, OC3, 
                    OC4, OC, OP, CL, CR, CU, FE, PB, MG, MN, MF, MO, NI, 
                    RCFM, NO3, P, K, RB, `NA`, SE, SI, S, SR, SOIL, SO4, TI, V, ZN, ZR) 
    
   
    # merge speciation data together:
    all_spec_df <- bind_rows(csn_spec_df, improve_spec_df) %>% 
      filter_at(vars(AL, AS, BR, CA, EC1, EC2, EC3, EC, OC1, OC2, OC3, 
                     OC4, OC, OP, TC, CHL, CL, CR, CU, FE, PB, MG, MN, MT, MF, MO, NI, 
                     RCTM, RCFM, NO3, P, K, RB, `NA`, SE, SI, S, SR, SOIL, SO4, 
                     TI, N2, V, ZN, ZR, CM_calculated, ammNO3, ammSO4),
                any_vars(!is.na(.))) %>% 
      distinct()
    
  return(all_spec_df)
  
} # end function
