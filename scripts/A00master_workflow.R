# Emma Krasovich Southworth, emmars@stanford.edu
# Sam Heft-Neal, samhn@stanford.edu
# Ayako Kawano, akawano@stanford.edu
# Created: October 29, 2022 | Last Updated: Jan 31, 2023
# Description: Full workflow for reading in, cleaning, data for the wildfire speciation project.
# This script will load all relevant libraries + functions that are necessary for the workflow,
# and run regressions of interest.

# source the functions and libraries
source(here::here("0_functions.R"))
source(here::here("0_packages.R"))
data_gdrive_fp = '/Users/ekrasovich/Library/CloudStorage/GoogleDrive-emmars@stanford.edu/Shared drives/echolab:data/'
wip_gdrive_fp = paste0(data_gdrive_fp, 'wildfire_speciation')
# get rid of scientific notation
options(scipen = 999)


# SOURCE FUNCTIONS:
list.files(here::here(), full.names=T, pattern = 'A0') %>%
  str_subset("\\.R$") %>%
  str_subset('workflow', negate = TRUE) %>%
  # str_subset('ZZ', negate = TRUE) %>%
  # str_subset('compare_high_smoke_to_low_smoke_days', negate = TRUE) %>%
  walk(source)


wildfire_plan <- drake_plan(
  
  # read in files
  us_shp = target(
    st_read(file_in(!!file.path(data_gdrive_fp, 'boundaries/gadm/gadm41_USA_shp/gadm41_USA_1.shp'))) %>% 
      filter(!NAME_1 %in% c('Alaska', 'Hawaii'))
    ),
  
  # # --------------------------------------------------------------------------------
  # # Step 1) read in raw IMPROVE speciation files, read them in, and clean up 
  # # --------------------------------------------------------------------------------
  # 1a. get list of raw IMPROVE speciation files
  improve_spec_file_list = list.files(file.path(wip_gdrive_fp, 'raw'), 
                                      pattern = 'IMPAER', full.names = TRUE),
  
  # 1b. read in and clean IMPROVE speciation data
  raw_IMPROVE_spec_df = target(
    read_in_IMPROVE_specation_data(
      improve_spec_file_list,
      improve_site_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/IMPROVE_sites.xlsx')),
      aqs_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/aqs_sites.csv')),
      aqs_monitors_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/aqs_monitors.csv'))
      )
    ),
  
  
  # # --------------------------------------------------------------------------------
  # # Step 2) read in raw CSN speciation files, read them in, and clean up 
  # # --------------------------------------------------------------------------------
  # 2a. get list of raw CSN speciation files
  csn_spec_file_list = list.files(file.path(wip_gdrive_fp, 'raw'), 
                                  pattern = 'daily_SPEC', full.names = TRUE), 
  
  # 2b. get list of raw PM2.5 files from EPA (88101)
  csn_pm25_list = list.files(file.path(wip_gdrive_fp, 'raw'), 
                             pattern = 'daily_88101', full.names = TRUE),
  
 
  # 2b. read in and clean CSN speciation data and merge with EPA PM2.5 data
  raw_CSN_spec_df = target(
    read_in_CSN_speciation_data(csn_spec_file_list, 
                                csn_pm25_list,
                                improve_params_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/IMPROVE_parameters.csv')))
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
      distinct(site_id, state_name, long, lat, epsg),
    CONUS_spec_sites_df_out = write.csv(CONUS_spec_sites_df,
                                         file_out(!!file.path(
                                           wip_gdrive_fp, "intermediate/CONUS_combined_speciation_sites_2006_2020.csv"))),
  
  # # --------------------------------------------------------------------------------
  # # Step 4) merge smoke PM2.5 data with speciation data
  # # --------------------------------------------------------------------------------
  # 4a. read in gridded Smoke PM2.5 data and merge in with monitor speciation data by matching grid cells
    spec_w_smoke_pm_df = target(
      join_speciation_w_gridded_smokePM(
        CONUS_spec_df,
        pm_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')),
        #pm_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/smokePM_predictions_20060101-20230630.rds')),
        grid_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/10km_grid_wgs84.shp'))
        )),
  
  # # --------------------------------------------------------------------------------
  # # Step 5)  merge speciation sites with smoke plume data
  # # --------------------------------------------------------------------------------
  spec_smoke_plumes_df = target(
    merge_sites_w_smoke_plumes(
      smoke_plumes_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/hms_smoke_plumes.rds')),
      spec_w_smoke_pm_df,
      CONUS_spec_sites_df)),
  
  
  
  # # --------------------------------------------------------------------------------
  # # Step 6) final cleaning and prepping data for analysis
  # # --------------------------------------------------------------------------------
    clean_PMspec_df = target(
      final_speciation_cleaning_pre_analysis(spec_smoke_plumes_df)),
  
    # save out
    clean_PMspec_df_out = write.csv(clean_PMspec_df,
                                         file_out(!!file.path(wip_gdrive_fp, "clean/station_PMspeciation_2006_2020.csv")),
                                     row.names = FALSE),

  # save a list of final sites
  clean_PMspec_sites_df = clean_PMspec_df %>%
    distinct(Dataset, site_id, lat, long, epsg),
  clean_PMspec_sites_df_out = write.csv(clean_PMspec_sites_df,
                                   file_out(!!file.path(wip_gdrive_fp, "clean/final_PMspeciation_sites.csv")),
                                   row.names = FALSE),
 
) # end plan

# make plan
if(file.exists(".drake")){
  drake_cache = drake::drake_cache(".drake")
} else {
  drake_cache = drake::new_cache(".drake")
}
make(wildfire_plan, cache = drake_cache)

  # # --------------------------------------------------------------------------
  # # --------------------------------------------------------------------------
  # # FIGURES + ANALYSIS
  # # --------------------------------------------------------------------------
  # # --------------------------------------------------------------------------
  
  # # ----------------------------------------------------------------------------
  # # # Step 7)  create descriptive figures + figures for paper
  # # ---------------------------------------------------------------------------
  # 
  # # 7a) Figure 1
  # # create map of monitoring stations for Figure 1
  # monitor_sf = target(
  #   create_map_of_monitoring_stations(clean_pm_spec_df)),
  # 
  # # 7b) Figure 2
  # create plots of smoke coefficients - both versions (pct change + log scale of levels)
  #reg_df = target(create_smoke_coeff_plots(clean_PMspec_df)),
  # 
  # # 7c) Figure 3
  # # create regional map for regional plot
  # us_region_map = target(
  #   create_us_region_map(st_fp = file_in(!!file.path(data_gdrive_fp, 'boundaries/all_national_states.rds')))
  #   ),
  # 
  # # create plot of regional coefficients in logs 
  # regional_coeffs = target(
  #   create_regional_coeff_plots(clean_pm_spec_df)),
  # 
  # # 7d) Figure 4
  # # create time series trend plots to compare contribution of smoke vs nonsmoke
  # pred_spec_conc_long_df = target(
  #   create_attributable_frac_time_series(clean_pm_spec_df))
  
  
  
  # TO DO
  # . get list of raw CSN speciation files
  # create a site crosswalk
  # create a chemical and cas crosswalk
  
  
  

  

  
  
  
  
  
  
  #   # read each in
  #   csn_spec_files = map_df(csn_spec_file_list, function(current_file) {
  #     # print current index
  #     print(paste0("current index is: ", match(current_file, unique(csn_spec_file_list)), ' of ', length(csn_spec_file_list)))
  #     # read in
  #     current_df <- read_excel(current_file)
  #   }) %>% 
  #     bind_rows(),
  #   
  #   # save out
  #   csn_spec_files_out = write.csv(csn_spec_files, 
  #                                      file_out(!!file.path(wip_gdrive_fp, "intermediate/CSN_Aerosol_combined.csv"))),
  # 
  # 
  # # 1b. now read in all the speciation files from both CSN and Improve
  # raw_spec_df = target(
  #   read_and_clean_speciation_data(
  #     improve_spec_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/IMPROVE_Aerosol_combined.csv')),
  #     csn_spec_files)
  # ),
  # 
  # --------------------------------------------------------------------------------
  # STEP 2) read in monitors and sites from EPA's AQS and merge with location data (lat/long/epsg)
  # --------------------------------------------------------------------------------
  # 2a) retrieve EPSG for all sites + merge into speciation data so that every site has lat, long, epsg, elevation
  # cleaned_spec_df = target(
  #   merge_projections_w_spec_df(
  #     raw_spec_df,
  #     aqs_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/aqs_sites.csv')),
  #     monitors_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/aqs_monitors.csv')),
  #     csn_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/csn_sites.xlsm')),
  #     improve_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/improve_sites.xlsm')),
  #     improve_prelim_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/improve_prelim_sites.xlsx'))
  #   )),
  # 
  # 
  # # 2c) save out a dataframe of the speciation sites for CSN + IMPROVE
  # all_spec_sites = cleaned_spec_df %>%
  #   distinct(SiteCode, Latitude, Longitude, Dataset, Elevation, epsg) %>%
  #   filter(!is.na(Latitude)),
  # all_spec_sites_out = write.csv(all_spec_sites,
  #                                file_out(!!file.path(wip_gdrive_fp, 'intermediate/all_improve_csn_speciation_sites.csv'))),
  # 
  # 
 
  # 
  # # --------------------------------------------------------------------------------
  # # # Step 4)  merge improve smoke density data with fires, calculate distance to nearest monitor
  # # --------------------------------------------------------------------------------
  # improve_smoke_dens_fire_dist = target(
  #   calculate_distance_monitor_nearest_fire(
  #     fires_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/hms_fires.RDS')),
  #     improve_smoke_dens_df,
  #     all_spec_sites
  #     )),
  # 
  # 
  # # --------------------------------------------------------------------------------
  # # # Step 5)  merge speciation data with smoke and plume data
  # # --------------------------------------------------------------------------------
  #   pm_plume_speciation_df = target(
  #     merge_speciation_plumes_and_pm(
  #       cleaned_spec_df,
  #       improve_smoke_dens_fire_dist,
  #       pm_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')),
  #       grid_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/10km_grid_wgs84.shp'))
  #       )),
  # 
  # # --------------------------------------------------------------------------------
  # # # Step 6)  calculate different fractions + add vars relevant for data exploration, clean data
  # # --------------------------------------------------------------------------------
  # clean_pm_spec_df = target(
  #   calculate_and_clean_spec_categories(pm_plume_speciation_df)
  #   ),
  # 
  # # save out
  # clean_pm_spec_df_out = write.csv(clean_pm_spec_df,
  #                                      file_out(!!file.path(wip_gdrive_fp, "clean/pm_plume_speciation_at_sites.csv")),
  #                                  row.names = FALSE),
  # 
  # clean_pm_spec_sites_df = clean_pm_spec_df %>% 
  #   distinct(Dataset, SiteCode, Latitude, Longitude),
  # clean_pm_spec_sites_df_out = write.csv(clean_pm_spec_sites_df,
  #                                  file_out(!!file.path(wip_gdrive_fp, "clean/reg_spec_sites.csv")),
  #                                  row.names = FALSE),
  # 
  
    


