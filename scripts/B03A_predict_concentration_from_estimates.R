# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species
# 
# grid_fp = file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')
# loadd(c(clean_PMspec_df,
#         regionalPMcoeffs_normalized, full_samplePM_df),
#       cache = drake::drake_cache(".drake"))

# function
predict_concentration_from_estimates <- function(grid_fp, clean_PMspec_df, 
                                                   full_samplePM_df, regionalPMcoeffs_normalized) {
  # read in grid
  grid_10km <- st_read(grid_fp) %>%
    st_transform(4326)
  
  # make points spatial data for monitors in df
  points <- clean_PMspec_df %>% 
    # get a unique list of sites
    distinct(Dataset, site_id, long, lat) %>% 
    # convert so sf object
    st_as_sf(coords = c(x ='long', y="lat"), crs = st_crs(grid_10km))
  
  # pull id cells that intersect with the monitoring sites
  points$grid_id_10km <- apply(
    st_intersects(grid_10km, points, sparse = FALSE), 2, 
    function(col) {
      grid_10km[which(col), ]$ID
    }) 
  
  sites_w_grid_cells <- points %>% 
    st_drop_geometry() %>% 
    left_join(grid_10km %>% 
                rename(grid_id_10km = 'ID') %>% 
                st_drop_geometry(), by = 'grid_id_10km') 
  
  # select the vars that are needed for the
  reg_df <- clean_PMspec_df %>%
    mutate(monitor_month = paste0(site_id,"_",month)) %>%
    dplyr::select(-ZR, -SOIL)
  
  rm(grid_10km, points, clean_PMspec_df) # drop to save memory
  
  # ------------------------------------------------------------------
  # REGIONAL ANALYSIS ----------------------------------------
  # ------------------------------------------------------------------
  coeffs <- regionalPMcoeffs_normalized %>% 
    dplyr::select(region, species, species_name, pm_type, Estimate) %>% 
    pivot_wider(names_from = 'pm_type', names_prefix = "est_", values_from = 'Estimate')
  
  rm(regionalPMcoeffs_normalized)
  
  # -----------------------------------------------------------------------------
  # map over each species and predict
  # -----------------------------------------------------------------------------
  # predict over our sample but with regional betas
  predicted_spec_smoke_conc <- reg_df %>%
    dplyr::select(Dataset, region, year, month, Date,
                  smoke_day, site_id, MF_adj, smokePM,
                  nonsmokePM_MF, AL:ZN) %>%
    pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>%
    filter(!is.na(conc_val)) %>%
    left_join(coeffs, by = c('region', 'species')) %>%
    # filter(nonsmokePM_MF >= 0 & smokePM >= 0) %>% 
    mutate(pred_spec_nonsmoke = est_nonsmokePM_MF*nonsmokePM_MF,
           pred_spec_smoke = est_smokePM*smokePM) %>%
    # filter(pred_spec_smoke >= 0 & pred_spec_nonsmoke >= 0) %>% 
    mutate(pred_spec_allPM = pred_spec_smoke + pred_spec_nonsmoke) %>% 
    # make sure we know which grid cell the monitors are in
    left_join(sites_w_grid_cells) 
  
  rm(coeffs, sites_w_grid_cells)
  
  # parellelize
  #future::plan(multisession)
  
  print(paste('merging predictions and concentration for all days'))
  # add to full dataset before averaging, need all days not just smoke days  
  pred_conc_all_days <- full_samplePM_df %>% 
    left_join(predicted_spec_smoke_conc %>%
                dplyr::select(Dataset, region, year, month, Date, site_id, grid_id_10km, 
                              species, species_name, pred_spec_nonsmoke, pred_spec_smoke, pred_spec_allPM),
              by = c("Date","grid_id_10km", "region")) 
  # pred_conc_all_days <- read_fst(file.path(data_fp, "intermediate/predicted_concentrations_all_days.fst"))
  
  rm(full_samplePM_df, predicted_spec_smoke_conc)
  
  # pivot reg df longer to prep for regression
  obs_species_conc_longer <- reg_df %>% 
    pivot_longer(cols = c(AL:ZN), names_to = 'species', values_to = 'obs_conc') %>% 
    dplyr::select(region, year, month, monitor_month, Date, site_id, species, obs_conc) 
  rm(reg_df)
  
  # join
  pred_obs_joined_df <- obs_species_conc_longer %>% 
    left_join(pred_conc_all_days, 
              by = c('species', 'Date', 'year', 'month', 'site_id', "region"))
  rm(pred_conc_all_days)
  
  return(pred_obs_joined_df) 
 }
#   
#   rm(obs_species_conc_longer)
#   
#   # run model to determine if the concentration in absolute terms (aka levels)
#   # is significantly increasing over time
#   pred_smoke_levels_df <-pred_obs_joined_df %>% 
#     # select only predicted concentrations in smoke
#     dplyr::select(Date, species, region, year, month, pred_spec_smoke, site_id, grid_id_10km, monitor_month) %>% 
#     pivot_wider(names_from = 'species', values_from = 'pred_spec_smoke') 
#   
#   # rm(pred_conc_all_days)
#   
#   print(paste('running model for smoke levels'))
#   smoke_levels_inc_model = feols(c(AL,AS,BR,CA,CL, CR,CU, EC, FE, K, MG, MN, 
#                                    `NA`, NI, NO3, OC, P, PB, RB, S, 
#                                    SE, SI, SO4, SR, TI, V, ZN)
#                                  ~  year | monitor_month,
#                                  data = pred_smoke_levels_df, cluster = 'site_id')
#   # etable(smoke_levels_inc_model)
#   # summary(smoke_levels_inc_model)
#   
#   ## calculate 95 CI%
#   CIsmoke <- confint(smoke_levels_inc_model) %>% 
#     rename(s_CI25 = `2.5 %`,
#            s_CI975 = `97.5 %`,
#            species = 'lhs') %>% 
#     dplyr::select(-id, -coefficient)
#   
#   # get coefficients and prepare for plotting
#   coeffs_smoke <- coeftable(smoke_levels_inc_model) %>%
#     rename(s_pval = 'Pr(>|t|)',
#            s_se = 'Std. Error',
#            species = 'lhs',
#            pred_smoke_coeff = 'Estimate') %>%
#     mutate(s_pval = round(s_pval, digits = 3)) %>%
#     dplyr::select(-id, -coefficient) %>%
#     left_join(CIsmoke, by = 'species') %>% 
#     #mutate(across(where(is.numeric), ~ round(.x, 9))) %>% 
#     mutate(s_CI = paste0("[", s_CI25, ", ", s_CI975, "]")) %>% 
#     dplyr::select(species, pred_smoke_coeff, s_CI, s_pval) 
#   
#   
#   # # RUN FOR NONSMOKE TO SEE IF THIS SIGNIFICANTLY DECREASING OVER TIME:
#   # # run model to determine if the concentration in absolute terms (aka levels)
#   # # is significantly increasing over time
#   pred_nonsmoke_levels_df <-pred_obs_joined_df %>%
#     # select only predicted concentrations in smoke
#     dplyr::select(Date, species, region, year, month, pred_spec_nonsmoke, site_id, grid_id_10km, monitor_month) %>%
#     pivot_wider(names_from = 'species', values_from = 'pred_spec_nonsmoke')
#   
#   print(paste('running model for nonsmoke levels'))
#   nonsmoke_levels_model = feols(c(AL,AS,BR,CA,CL, CR,CU, EC, FE, K, MG, MN,
#                                   `NA`, NI, NO3, OC, P, PB, RB, S,
#                                   SE, SI, SO4, SR, TI, V, ZN)
#                                 ~  year | monitor_month,
#                                 data = pred_nonsmoke_levels_df, cluster = 'site_id')
#   # etable(nonsmoke_levels_model)
#   # summary(nonsmoke_levels_model)
#   
#   ## calculate 95 CI%
#   CInonsmoke <- confint(nonsmoke_levels_model) %>%
#     rename(ns_CI25 = `2.5 %`,
#            ns_CI975 = `97.5 %`,
#            species = 'lhs') %>%
#     dplyr::select(-id, -coefficient)
#   
#   # get coefficients and prepare for plotting
#   coeffs_nonsmoke <- coeftable(nonsmoke_levels_model) %>%
#     rename(ns_pval = 'Pr(>|t|)',
#            ns_se = 'Std. Error',
#            species = 'lhs',
#            pred_nonsmoke_coeff = 'Estimate') %>%
#     mutate(ns_pval = round(ns_pval, digits = 3)) %>%
#     dplyr::select(-id, -coefficient) %>%
#     left_join(CInonsmoke, by = 'species') %>%
#     dplyr::select(species, pred_nonsmoke_coeff, ns_pval, ns_CI25, ns_CI975) %>%
#     mutate(ns_CI = paste0("[", ns_CI25, ", ", ns_CI975, "]")) %>%
#     dplyr::select(species, pred_nonsmoke_coeff,
#                   ns_CI,
#                   ns_pval)
#   
#   # combine both predicted levels
#   predicted_levels_df <- left_join(
#     coeffs_nonsmoke, coeffs_smoke) %>%
#     dplyr::select(-ns_CI, -s_CI) %>%
#     mutate(across(where(is.numeric), ~ signif(.x, 6)))
#   # write.csv(predicted_levels_df, file.path(data_fp, 'intermediate/SI_table4.csv'))
#   # future::plan(NULL)
#   
#   # do some adjusting to print the right amount of significant digits
#   # knitr::kable(predicted_levels_df, "latex")
#   rm(predicted_levels_df, coeffs_nonsmoke, coeffs_smoke,
#      nonsmoke_levels_model, CIsmoke, CInonsmoke, pred_nonsmoke_levels_df, pred_smoke_levels_df)
#   
#   return(pred_obs_joined_df) 
#   
# } # end function
# 
# 
# 
