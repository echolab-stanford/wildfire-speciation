# Emma Krasovich Southworth, emmars@stanford.edu
# Created: October 29, 2022 | Last Updated: Jan 3, 2024
# Description: Full workflow for reading in, cleaning, data for the wildfire speciation project.
# This script will load all relevant libraries + functions that are necessary for the workflow,
# and run regressions of interest.


# source the functions and libraries
source(here::here("0_functions.R"))
source(here::here("0_packages.R"))
data_fp = "~/BurkeLab Dropbox/projects/ws_project/data"
#dropbox_data_fp = '/Users/ekrasovich/Library/CloudStorage/GoogleDrive-emmars@stanford.edu/Shared drives/echolab:data'
gee_data_fp = "~/My Drive (emmars@stanford.edu)/Research/G_ECHOLab/GEE_Echolab_Export"
results_fp = "~/BurkeLab Dropbox/projects/ws_project/figures"
boundaries_fp = '~/BurkeLab Dropbox/data/boundaries'
# get rid of scientific notation
options(scipen = 999)


# SOURCE FUNCTIONS:
list.files(here::here(), full.names=T, pattern = '(^A|^B|^C)') %>%
  str_subset("\\.R$") %>%
  str_subset('workflow', negate = TRUE) %>%
  str_subset('combine_speciation_w_CASno_and_RSLs', negate = TRUE) %>%
  str_subset('C02_CSN_V_improve_sensitivty', negate = TRUE) %>%
  #str_subset('estimating_conc_response', negate = TRUE) %>%
 # str_subset('B04D_randomization', negate = TRUE) %>% 
walk(source)


wildfire_plan <- drake_plan(
  
  # read in files
  us_shp = target(
    st_read(file_in(!!file.path(data_fp, "raw/gadm36_USA_shp/gadm36_USA_1.shp"))) %>% 
      filter(!NAME_1 %in% c('Alaska', 'Hawaii'))
    ),
  
  # # --------------------------------------------------------------------------------
  # # Step 1) read in raw IMPROVE speciation files, read them in, and clean up 
  # # --------------------------------------------------------------------------------
  # 1a. get list of raw IMPROVE speciation files
  improve_spec_file_list = list.files(file.path(data_fp, 'raw'), 
                                      pattern = 'IMPAER', full.names = TRUE),
  
  # 1b. read in and clean IMPROVE speciation data
  raw_IMPROVE_spec_df = target(
    read_in_IMPROVE_specation_data(
      improve_spec_file_list,
      improve_site_fp = file_in(!!file.path(data_fp, "raw/IMPROVE_sites.xlsx")),
      aqs_sites_fp = file_in(!!file.path(data_fp, "raw/aqs_sites.csv")),
      aqs_monitors_fp = file_in(!!file.path(data_fp, "raw/aqs_monitors.csv"))
      )
    ),

  # # --------------------------------------------------------------------------------
  # # Step 2) read in raw CSN speciation files, read them in, and clean up 
  # # --------------------------------------------------------------------------------
  # 2a. get list of raw CSN speciation files
  csn_spec_file_list = list.files(file.path(data_fp, 'raw'), 
                                  pattern = 'daily_SPEC', full.names = TRUE), 
  
  # 2b. get list of raw PM2.5 files from EPA (88101)
  csn_pm25_list = list.files(file.path(data_fp, 'raw'), 
                             pattern = 'daily_88101', full.names = TRUE),


  # 2b. read in and clean CSN speciation data and merge with EPA PM2.5 data
  raw_CSN_spec_df = target(
    read_in_CSN_speciation_data(csn_spec_file_list, 
                                csn_pm25_list,
                                csn_param_fp = file_in(!!file.path(data_fp, 'raw/CSN_parameters.csv')))
  ),

  # # --------------------------------------------------------------------------------
  # # Step 3) merge CSN and IMPROVE speciation data together + transform to same coordinate system
  # # --------------------------------------------------------------------------------
  # 3a. create a combined dataset of CONUS speciation sites from IMPROVE and CSN data
  CONUS_spec_df = target(
      merge_CSN_IMPROVE_speciation_data(
        raw_CSN_spec_df, raw_IMPROVE_spec_df, us_shp)),
    
  
  # 3b. save out "raw" speciation sites, pre-dropping monitor-days for analysis
    CONUS_spec_sites_df = CONUS_spec_df %>% 
      distinct(Dataset, site_id, state_name, long, lat, epsg),
  
  # # --------------------------------------------------------------------------------
  # # Step 4) merge smoke PM2.5 data with speciation data
  # # --------------------------------------------------------------------------------
  # 4a. read in gridded Smoke PM2.5 data and merge in with monitor speciation data by matching grid cells
    spec_w_smoke_pm_df = target(
      join_speciation_w_gridded_smokePM(
        CONUS_spec_df,
        pm_fp = file_in(!!file.path(data_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')),
        grid_fp = file_in(!!file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp'))
        )),
  
  # --------------------------------------------------------------------------------
  # Step 5) final cleaning and prepping data for analysis
  # --------------------------------------------------------------------------------
  # 6a. final cleaning steps
  clean_PMspec_df = target(
      final_speciation_cleaning_pre_analysis(spec_w_smoke_pm_df)),

    # save out
    clean_PMspec_df_out = write.csv(clean_PMspec_df,
                                         file_out(!!file.path(data_fp, "clean/station_PMspeciation_2006_2020.csv")),
                                     row.names = FALSE),

  # 6b. save a list of final sites
  clean_PMspec_sites_df = clean_PMspec_df %>%
    distinct(Dataset, site_id, lat, long, epsg),
  clean_PMspec_sites_df_out = write.csv(clean_PMspec_sites_df,
                                   file_out(!!file.path(data_fp, "clean/final_PMspeciation_sites.csv")),
                                   row.names = FALSE),

  # 6c. read in constructed a parameter crosswalk for the figures
  parameter_categories = read.csv(file_in(!!file.path(data_fp, 'intermediate/xwalk_species_type.csv'))) %>% 
    mutate(species = str_to_upper(species)),
  

  # --------------------------------------------------------------------------------
  # Step 6) expand PM data so that every day in sample is accounted for
  # --------------------------------------------------------------------------------
  full_samplePM_df = target(
    expandPM_data_to_all_days(pm_fp = file_in(!!file.path(data_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')),
                              grid_fp = file_in(!!file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')),
                              us_states_fp = file_in(!!(file.path(data_fp, 'raw/all_national_states.rds'))))),

  # # --------------------------------------------------------------------------
  # # --------------------------------------------------------------------------
  # # PART B: FIGURES + ANALYSIS
  # # --------------------------------------------------------------------------
  # # --------------------------------------------------------------------------
  # set up palette for each species
      spec_pal = c("#35274A", "#B40F20", "#E58601", "#E2D200", "#46ACC8", "#0B775E","#CAB38C", "#899DA4"),
      pm_pal = c( "grey70", "#9B110E"),
      region_pal = c("#DC863B", "#FAD510","#649373", "#1B5656", "#5A283E", "#F2300F"),
  
  # # --------------------------------------------------------------------------------
  # # FIGURE 1
  # # --------------------------------------------------------------------------------
  # # a) Figure 1A - create map of monitoring stations 
  monitor_map = target(
    create_map_of_monitoring_stations(
      clean_PMspec_df,
      us_states_fp = file_in(!!(file.path(data_fp, 'raw/all_national_states.rds')))
      )),
  

  # b) Figure 1B - create time series of raw PM data
  avg_PM_monthly_mean = target(
    plot_monthly_mean_raw_PM_data(clean_PMspec_df, parameter_categories, pm_pal)),
  
  
  # c) Figure 1C - create time series of monthly chemical species averages
  avg_spec_monthly_mean = target(
    plot_monthly_mean_raw_all_specices(clean_PMspec_df, 
                                       parameter_categories, 
                                       spec_pal)),

  # # --------------------------------------------------------------------------------
  # # FIGURE 2
  # # --------------------------------------------------------------------------------
  # 2a) create plots of smoke coefficients - both versions (pct change + log scale of levels)
  full_samp_PMcoeffs_normalized = target(
    create_smoke_coeff_plots(clean_PMspec_df,
                             parameter_categories,
                             spec_pal)),
  full_samp_coeffs_out = write.csv(full_samp_PMcoeffs_normalized,
                                        file_out(!!file.path(data_fp, "clean/full_sample_coeffs.csv")),
                                        row.names = FALSE),
  
  # create table for coeffs to put in supplementary info
  tab_df = full_samp_PMcoeffs_normalized %>% 
    dplyr::select(species = 'species_name', est = 'Estimate', se, tval = `t value`, pval, 
                  baseline_NS_conc = 'avg_nonsmoke_spec_conc', 
                  norm_est, norm_CI25, norm_CI975) %>% 
    mutate(norm_est = 100*norm_est, 
           norm_CI25= 100*norm_CI25, norm_CI975 = 100*norm_CI975),
  #print(xtable(tab_df, digits = 4), include.rownames=FALSE)

  # 2b) create regional map for regional plot
  us_region_map = target(
    create_us_region_map(
      us_states_fp = file_in(!!(file.path(data_fp, 'raw/all_national_states.rds'))),
      region_pal)),
  
  # 2c) Plot the regional coefficients
  regionalPMcoeffs_normalized = target(
    create_regional_coeff_plots(clean_PMspec_df, parameter_categories, region_pal)),
  regional_coeffs_out = write.csv(regionalPMcoeffs_normalized,
                                  file_out(!!file.path(data_fp, "clean/regional_coeffs.csv")),
                                  row.names = FALSE),

  # create table for coeffs to put in supplementary info
  reg_tab_df = regionalPMcoeffs_normalized %>% 
    filter(pm_type == 'smokePM') %>% 
    dplyr::select(region, species = 'species_name', est = 'Estimate', se, pval, 
                  baseline_avg, norm_est, norm_CI25, norm_CI975) %>% 
    mutate(norm_est = 100*norm_est, 
           norm_CI25= 100*norm_CI25, norm_CI975 = 100*norm_CI975),
  
  
  # print(xtable(reg_tab_df, digits = 4), include.rownames=FALSE)

  # --------------------------------------------------------------------------------
  # FIGURE 3 - attributable fraction over time
  # --------------------------------------------------------------------------------
  # 3a) create data using regional model coeffs + predicting attributable fraction in smoke 
  pred_regional_attributable_preds = target(
    create_attributable_frac_w_regional_coeffs(
      grid_fp = file_in(!!file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')),
      clean_PMspec_df, parameter_categories, pm_pal, full_samplePM_df, regionalPMcoeffs_normalized)),

  # # 3b) create plot attributable fraction due to smoke plot
  attr_frac_df = target(
    plot_attributable_frac_trends(pred_regional_attributable_preds, spec_pal)),
  
  ## 3c) make table of starting and ending concentration for wildfire attributable concentrations
  change_in_wildfire_att_conc_table = pred_regional_attributable_preds %>% 
    group_by(year, species, label) %>% 
    dplyr::summarize(avg_ann_conc = mean(conc)) %>% 
    ungroup() %>% 
    filter(label == 'Attributable concentration to smoke PM2.5') %>% 
    filter(year %in% c(2006, 2020)) %>% 
    filter(species != 'CHL') %>% 
    dplyr::select(-label) %>% 
    pivot_wider(names_from = 'year', values_from = 'avg_ann_conc') %>% 
    mutate(increase_x = `2020`/`2006`),
  # print(xtable(change_in_wildfire_att_conc_table, digits = 7), include.rownames=FALSE)
  
  # --------------------------------------------------------------------------------
  # FIGURE 4 BURNED STRUCTURES + RANDOMIZATION INFERENCE
  # --------------------------------------------------------------------------------
  # 4a) Process fire + structures burned data
  globfire_structure_joined_df = target(
    processing_burned_structures_data(
      globfire_fp = file_in(!!(file.path(data_fp, 'intermediate/globfire/globfire_na_final_area_2006-2020.shp'))),
      mtbs_fp = file_in(!!(file.path(data_fp, 'intermediate/mtbs/mtbs_perims_DD.shp'))),
      nifc_fp = file_in(!!(file.path(data_fp, "intermediate/fire_locations/WFIGS_-_Wildland_Fire_Locations_Full_History.csv"))),
      damaged_struct_fp = file_in(!!(file.path(data_fp, 'intermediate/HE_Structures_Destroyed_2022.xlsx')))
    )),
  
  # 4b) merge burned structures with speciation data and clean
  burned_struc_smoke_spec_df = target(
    merge_burned_structures_w_speciation(
      globfire_structure_joined_df, clean_PMspec_sites_df, clean_PMspec_df,
      grid_fp = file_in(!!file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp'))
    )),
  
  # 4c) estimate how an additional burned structure impacts the concentration of certain chemicals
  burned_structure_coeffs = target(
    estimating_conc_response_2_burned_structures(
      burned_struc_smoke_spec_df,
      parameter_categories,
      spec_pal)),
  
  # grab coefficients
  bs_coeffs_tab = burned_structure_coeffs %>% 
    filter(species %in% c('MG','TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR')) %>% 
    filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
    mutate(Est_per_1000 = Estimate*1000) %>% 
    dplyr::select(species_long, avg_no_struc_spec_conc, Estimate, Est_per_1000, CI25, CI975,
                  pct_change_struct_vs_nostruct, pct_change_CI25, pct_change_CI975),
  # print(xtable(bs_coeffs_tab, digits = -2), include.rownames=FALSE)

  
  # 4d) run permutation test through randomization inference
  permutation_results_df = target(
    randomization_inference_permutation_test(
      burned_struc_smoke_spec_df, parameter_categories)),
  
  #4e) calculate the p-value associated with each species randomization-style inference
  # The randomization inference p-value is the proportion of times the placebo treatment 
  # effect was larger than the estimated treatment effect.
  permutation_results_w_pval = permutation_results_df %>% 
    left_join(bs_coeffs_tab %>% 
                dplyr::select(species_long, Estimate)) %>% 
    group_by(species) %>% 
    mutate(coeff_over_est = case_when(
          Estimate > 0 ~ sum(coef > Estimate),
          Estimate < 0 ~ sum(coef < Estimate))) %>%
        ungroup() %>%
    distinct(species, coeff_over_est) %>% 
    mutate(pval = coeff_over_est/1000) %>%
      left_join(coef_df %>%
                  dplyr::select(-pval)) %>%
      dplyr::select(species, obs_est = 'Estimate', pval) %>%
      arrange(pval),
  # print(xtable(permutation_results_w_pval, digits = -2, include.rownames=FALSE))
  
  perm_result_df = target(
    plot_burned_structure_coeffs_rand_plot(
      bs_coeffs_tab, parameter_categories, permutation_results_w_pval)),
    

  # --------------------------------------------------------------------------------
  # FIGURE 5 
  # -------------------------------------------------------------------------------
  #5a) run pixel regression + plot the gridded concentrations for each sample period
  # AS_avg_period_gridded_preds_df = target(
  #   run_pixel_regression_and_plot_5yr_maps(
  #     clean_PMspec_df,
  #     parameter_categories,
  #     grid_fp = file_in(!!file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')),
  #     full_samplePM_df,
  #     current_species = 'AS')),

  # NI_avg_period_gridded_preds_df = target(
  #   run_pixel_regression_and_plot_5yr_maps(
  #     clean_PMspec_df,
  #     parameter_categories,
  #     grid_fp = file_in(!!file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')),
  #     full_samplePM_df,
  #     current_species = 'NI')),
  
  # combine datasets
  comb_tox_preds_df = read_fst(file.path(data_fp, "clean/AS_gridded_preds.fst")) %>% 
    dplyr::select(grid_id_10km, Date, AS_pred_conc = 'pred_grid_conc', year) %>%
    left_join(read.fst(file.path(data_fp, "clean/NI_gridded_preds.fst")) %>% 
                dplyr::select(grid_id_10km, Date, NI_pred_conc = 'pred_grid_conc', year),
              by = c('grid_id_10km', 'Date', 'year')) %>% 
    left_join(read.fst(file.path(data_fp, "clean/PB_gridded_preds.fst")) %>% 
                dplyr::select(grid_id_10km, Date, PB_pred_conc = 'pred_grid_conc', year),
              by = c('grid_id_10km', 'Date', 'year')),
  # save out
  comb_tox_preds_df_out = write.fst(comb_tox_preds_df,
    file_out(!!file.path(data_fp, 'clean/predicted_gridded_AS_PB_NI_combined.fst')),
    compress = 100),

  #5b) calculate population excess cancer cases
  annual_excess_cancer_cases = target(
    calculate_risk_from_pop_exposure(csn_param_fp = file_in(!!file.path(data_fp, "raw/CSN_parameters.csv")), 
                           grid_fp = file_in(!!file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')), 
                           tox_fp = file_in(!!file.path(data_fp, "raw/RSL Table/oehha_ca_IURs.csv")), 
                           us_pop_fp = file_in(!!file.path(data_fp, "raw/united-states-population-2024-02-08.csv")), 
                           us_region_map, 
                           region_pal)),
  
  cancer_table = annual_excess_cancer_cases %>% 
    filter(scenario == 'lifetime') %>% 
    filter(samp_period != '2011-2015') %>% 
    dplyr::select(-scenario) %>% 
    pivot_wider(names_from = 'samp_period', values_from = 'canc_burden') %>% 
    mutate(increase_x = `2016-2020`/`2006-2010`) %>% 
    group_by(species) %>% 
    mutate(tot_early_cases = sum(`2006-2010`, na.rm = TRUE),
           tot_late_cases = sum(`2016-2020`, na.rm = TRUE)) %>% 
    arrange(species, region),
    #xtable(cancer_table, digits = 1, include.rownames=FALSE)
) # end drake plan

# make plan
if(file.exists(".drake")) {
  drake_cache = drake::drake_cache(".drake")
} else {
  drake_cache = drake::new_cache(".drake")
}
make(wildfire_plan, cache = drake_cache)
  
  # --------------------------------------------------------------------------------
  # SUPPLEMENTARY FIGURES
  # -------------------------------------------------------------------------------
  # coeffs_robustness = target(
  #   create_coeff_sensitivity_SI_figures(
  #     clean_PMspec_df,
  #     parameter_categories,
  #     spec_pal)),
  # 
  # SB_agg_df = target(
  #   make_SI_burned_structures_time_series_plot(
  #     globfire_structure_joined_df))
    
  
 


