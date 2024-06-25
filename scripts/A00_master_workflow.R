# Emma Krasovich Southworth, emmars@stanford.edu
# Created: October 29, 2022 | Last Updated: Jan 3, 2024
# Description: Full workflow for reading in, cleaning, data for the wildfire speciation project.
# This script will load all relevant libraries + functions that are necessary for the workflow,
# and run regressions of interest.


# source the functions and libraries
source(here::here("0_functions.R"))
source(here::here("0_packages.R"))
data_fp = '/Users/ekrasovich/Library/CloudStorage/Dropbox-BurkeLab/projects/ws_project/data'
gdrive_data_fp = '/Users/ekrasovich/Library/CloudStorage/GoogleDrive-emmars@stanford.edu/Shared drives/echolab:data'
gee_data_fp = '/Users/ekrasovich/My Drive/Research/G_ECHOLab/GEE_Echolab_Export'
results_fp = '/Users/ekrasovich/Library/CloudStorage/Dropbox-BurkeLab/projects/ws_project'
boundaries_fp = '/Users/ekrasovich/Library/CloudStorage/Dropbox-BurkeLab/data/boundaries/'
# get rid of scientific notation
options(scipen = 999)


# SOURCE FUNCTIONS:
list.files(here::here(), full.names=T, pattern = '(^A|^B|^C)') %>%
  str_subset("\\.R$") %>%
  str_subset('workflow', negate = TRUE) %>%
  str_subset('combine_speciation_w_CASno_and_RSLs', negate = TRUE) %>%
  str_subset('calculate_pop', negate = TRUE) %>%
  str_subset('C02_CSN_V_improve_sensitivty', negate = TRUE) %>%
  str_subset('CCcombine_gridded_predictions_w_toxicity_vals', negate = TRUE) %>%
  #str_subset('estimating_conc_response', negate = TRUE) %>%
 # str_subset('B04D_randomization', negate = TRUE) %>% 
walk(source)


wildfire_plan <- drake_plan(
  
  # read in files
  us_shp = target(
    st_read(file_in(!!file.path(gdrive_data_fp, 'boundaries/gadm/gadm41_USA_shp/gadm41_USA_1.shp'))) %>% 
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
      improve_site_fp = file_in(!!file.path(data_fp, 'raw/IMPROVE_sites.xlsx')),
      aqs_sites_fp = file_in(!!file.path(data_fp, 'raw/aqs_sites.csv')),
      aqs_monitors_fp = file_in(!!file.path(data_fp, 'raw/aqs_monitors.csv'))
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
    
  
  #3c. read in haps data and clean
  # haps_spec_file_list = list.files(file.path(data_fp, "raw"), 
  #                                  pattern = "daily_HAPS", full.names= TRUE),
  # 
   # haps_spec_df = target(
   #  read_in_HAPS_speciation_data(
   #    us_shp,
   #    haps_spec_file_list,
   #    haps_xwalk_fp = file_in(!!file.path(data_fp, 'intermediate/haps_xwalk_IUR_weight.csv')),
   #    spec_w_smoke_pm_df
   #    )),
  
  # # --------------------------------------------------------------------------------
  # # Step 4) merge smoke PM2.5 data with speciation data
  # # --------------------------------------------------------------------------------
  # 4a. read in gridded Smoke PM2.5 data and merge in with monitor speciation data by matching grid cells
    spec_w_smoke_pm_df = target(
      join_speciation_w_gridded_smokePM(
        CONUS_spec_df,
        pm_fp = file_in(!!file.path(data_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')),
        #pm_fp = file_in(!!file.path(data_fp, 'intermediate/smokePM_predictions_20060101-20230630.rds')),
        grid_fp = file_in(!!file.path(data_fp, 'intermediate/10km_grid_wgs84.shp'))
        )),
  
  # #4b. merge with haps
  # all_spec_df = target(
  #   merge_spec_smoke_w_HAPS(
  #     spec_w_smoke_pm_df, haps_spec_df)),
  
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
                              grid_fp = file_in(!!file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')),
                              us_states_fp = file_in(!!(file.path(boundaries_fp, 'all_national_states.rds'))))),
  
  
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
      us_states_fp = file_in(!!(file.path(boundaries_fp, 'all_national_states.rds')))
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
      us_states_fp = file_in(!!(file.path(boundaries_fp, 'all_national_states.rds'))),
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
  # pred_regional_attributable_preds = target(
  #   create_attributable_frac_w_regional_coeffs(
  #     grid_fp = file_in(!!file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')),
  #     clean_PMspec_df, parameter_categories, pm_pal, full_samplePM_df, regionalPMcoeffs_normalized)),
  # 
  # # 3b) create plot attributable fraction due to smoke plot
  # attr_frac_df = target(
  #   plot_attributable_frac_trends(pred_regional_attributable_preds, spec_pal)),
  

  # --------------------------------------------------------------------------------
  # FIGURE 4 BURNED STRUCTURES + RANDOMIZATION INFERENCE
  # --------------------------------------------------------------------------------
  # 4a) Process fire + structures burned data
  globfire_structure_joined_df = target(
    processing_burned_structures_data(
      globfire_fp = file_in(!!(file.path(data_fp, 'intermediate/globfire/globfire_na_final_area_2006-2020.shp'))),
      mtbs_fp = file_in(!!(file.path(data_fp, 'intermediate/mtbs/mtbs_perims_DD.shp'))),
      nifc_fp = file_in(!!(file.path(data_fp, "intermediate/fire_locations/WFIGS_-_Wildland_Fire_Locations_Full_History.csv"))),
      damaged_struct_fp = file_in(!!(file.path(data_fp, 'clean/HE_Structures_Destroyed_2022.xlsx')))
    )),
  
  # 4b) merge burned structures with speciation data and clean
  burned_struc_smoke_spec_df = target(
    merge_burned_structures_w_speciation(
      globfire_structure_joined_df, clean_PMspec_sites_df, clean_PMspec_df,
      grid_fp = file_in(!!(file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')))
    )),
  
  # 4c) estimate how an additional burned structure impacts the concentration of certain chemicals
  burned_structure_coeffs = target(
    estimating_conc_response_2_burned_structures(
      burned_struc_smoke_spec_df,
      parameter_categories,
      spec_pal)),
  
  bs_coeffs_tab = burned_structure_coeffs %>% 
    filter(species %in% c('MG','TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR')) %>% 
    filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
    mutate(pct_est = 100*no_struc_smoke_est,
           pct_CI25 = no_struc_smoke_CI25*100,
           pct_CI975 =no_struc_smoke_CI975*100) %>% 
    dplyr::select(species_long, avg_no_struc_spec_conc, no_struc_smoke_est, CI25, CI975,
                  pct_est, pct_CI25, pct_CI975),
  # print(xtable(bs_coeffs_tab, digits = -2), include.rownames=FALSE)
  
  # 4d) run permutation test through randomization inference
  # permutation_results_df = target(
  #   randomization_inference_permutation_test(
  #     burned_struc_smoke_spec_df, parameter_categories)),

  # --------------------------------------------------------------------------------
  # FIGURE 5 
  # -------------------------------------------------------------------------------
  # 5a) run pixel regression + plot the gridded concentrations for each sample period
  # PB_avg_period_gridded_preds_df = target(
  #   run_pixel_regression_and_plot_5yr_maps(
  #     clean_PMspec_df,
  #     parameter_categories,
  #     grid_fp = file_in(!!file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')),
  #     full_samplePM_df, current_species = 'PB')),

  # AS_avg_period_gridded_preds_df = target(
  #   run_pixel_regression_and_plot_5yr_maps(
  #     clean_PMspec_df,
  #     parameter_categories,
  # grid_fp = file_in(!!file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')),
  # full_samplePM_df, current_species = 'AS')),
  # 
  # NI_avg_period_gridded_preds_df = target(
  #   run_pixel_regression_and_plot_5yr_maps(
  #     clean_PMspec_df,
  #     parameter_categories,
  #     grid_fp = file_in(!!file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')),
  #     full_samplePM_df, current_species = 'NI')),

  
  #5b) calculate population excess cancer cases
  # annual_excess_cases = target(
  #   calculate_pop_exposure(
  #     csn_param_fp = file_in(!!file.path(data_fp, 'raw/CSN_parameters.csv')),
  #     grid_fp = file_in(!!file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')),
  #     tox_fp = file_in(!!file.path(data_fp, 'raw/RSL Table/oehha_ca_IURs.csv'))
  #     ))
  
  
  # --------------------------------------------------------------------------------
  # SUPPLEMENTARY FIGURES
  # -------------------------------------------------------------------------------
  coeffs_robustness = target(
    create_coeff_sensitivity_SI_figures(
      clean_PMspec_df,
      parameter_categories,
      spec_pal)),
  
  SB_agg_df = target(
    make_SI_burned_structures_time_series_plot(
      globfire_structure_joined_df))
    
  
 ) # end plan

# make plan
if(file.exists(".drake")) {
  drake_cache = drake::drake_cache(".drake")
} else {
  drake_cache = drake::new_cache(".drake")
}
make(wildfire_plan, cache = drake_cache)


