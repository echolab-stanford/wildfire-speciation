# Emma Krasovich Southworth, emmars@stanford.edu
# Sam Heft-Neal, samhn@stanford.edu
# Ayako Kawano, akawano@stanford.edu
# Created: October 29, 2022 | Last Updated: Jan 31, 2023
# Description: Full workflow for reading in, cleaning, data for the wildfire speciation project.
# This script will load all relevant libraries + functions that are necessary for the workflow,
# and run regressions of interest.

# set working directory:
setwd('/Users/ekrasovich/Desktop/ECHOLab Local/wildfire-speciation')

# source the functions and libraries
source("scripts/0_functions.R")
source("scripts/0_packages.R")
wip_gdrive_fp = '/Users/ekrasovich/Library/CloudStorage/GoogleDrive-emmars@stanford.edu/Shared drives/echolab:data/wildfire_speciation'


# SOURCE FUNCTIONS:
list.files('scripts', full.names=T, pattern = '(A|B)') %>%
  str_subset("\\.R$") %>%
  str_subset('workflow', negate = TRUE) %>%
  str_subset('attributable', negate = TRUE) %>%
  str_subset('regress', negate = TRUE) %>%
  walk(source)


wildfire_plan <- drake_plan(
  
  # --------------------------------------------------------------------------------
  # Step 1) read in raw speciation files, read them in, and clean up 
  # --------------------------------------------------------------------------------
    # 1a. get list of raw CSN files
    csn_spec_file_list = list.files(file.path(wip_gdrive_fp, 'raw'), 
                                    pattern = 'emmars', full.names = TRUE),
    # read each in
    csn_spec_files = map_df(csn_spec_file_list, function(current_file) {
      # print current index
      print(paste0("current index is: ", match(current_file, unique(csn_spec_file_list)), ' of ', length(csn_spec_file_list)))
      # read in
      current_df <- read_excel(current_file)
    }) %>% 
      bind_rows(),
    
    # save out
    csn_spec_files_out = write.csv(csn_spec_files, 
                                       file_out(!!file.path(wip_gdrive_fp, "intermediate/CSN_Aerosol_combined.csv"))),
  
  
  # 1b. now read in all the speciation files from both CSN and Improve
  raw_spec_df = target(
    read_and_clean_speciation_data(
      improve_spec_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/IMPROVE_Aerosol_combined.csv')),
      csn_spec_files)
  ),
  
  # --------------------------------------------------------------------------------
  # STEP 2) read in monitors and sites from EPA's AQS and merge with location data (lat/long/epsg)
  # --------------------------------------------------------------------------------
  # 2a) retrieve EPSG for all sites + merge into speciation data so that every site has lat, long, epsg, elevation
  cleaned_spec_df = target(
    merge_projections_w_spec_df(
      raw_spec_df,
      aqs_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/aqs_sites.csv')),
      monitors_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/aqs_monitors.csv')),
      csn_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/csn_sites.xlsm')),
      improve_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/improve_sites.xlsm')),
      improve_prelim_sites_fp = file_in(!!file.path(wip_gdrive_fp, 'raw/improve_prelim_sites.xlsx'))
    )),
  
  
  # 2c) save out a dataframe of the speciation sites for CSN + IMPROVE
  all_spec_sites = cleaned_spec_df %>%
    distinct(SiteCode, Latitude, Longitude, Dataset, Elevation, epsg) %>%
    filter(!is.na(Latitude)),
  all_spec_sites_out = write.csv(all_spec_sites,
                                 file_out(!!file.path(wip_gdrive_fp, 'intermediate/all_improve_csn_speciation_sites.csv'))),

  
  # --------------------------------------------------------------------------------
  # Step 3)  merge speciation sites with smoke plume data
  # --------------------------------------------------------------------------------
  improve_smoke_dens_df = target(
    merge_sites_w_smoke_plumes(
      smoke_plumes_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/hms_smoke_plumes.rds')),
      cleaned_spec_df,
      all_spec_sites)),

  # --------------------------------------------------------------------------------
  # # Step 4)  merge improve smoke density data with fires, calculate distance to nearest monitor
  # --------------------------------------------------------------------------------
  improve_smoke_dens_fire_dist = target(
    calculate_distance_monitor_nearest_fire(
      fires_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/hms_fires.RDS')),
      improve_smoke_dens_df,
      all_spec_sites
      )),


  # --------------------------------------------------------------------------------
  # # Step 5)  merge speciation data with smoke and plume data
  # --------------------------------------------------------------------------------
    pm_plume_speciation_df = target(
      merge_speciation_plumes_and_pm(
        cleaned_spec_df,
        improve_smoke_dens_fire_dist,
        pm_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')),
        grid_fp = file_in(!!file.path(wip_gdrive_fp, 'intermediate/10km_grid_wgs84.shp'))
        )),

  # --------------------------------------------------------------------------------
  # # Step 6)  calculate different fractions + add vars relevant for data exploration, clean data
  # --------------------------------------------------------------------------------
  clean_pm_spec_df = target(
    calculate_and_clean_spec_categories(pm_plume_speciation_df)
    ),

  # save out
  clean_pm_spec_df_out = write.csv(clean_pm_spec_df,
                                       file_out(!!file.path(wip_gdrive_fp, "clean/pm_plume_speciation_at_sites.csv")),
                                   row.names = FALSE),
  
  clean_pm_spec_sites_df = clean_pm_spec_df %>% 
    distinct(Dataset, SiteCode, Latitude, Longitude),
  clean_pm_spec_sites_df_out = write.csv(clean_pm_spec_sites_df,
                                   file_out(!!file.path(wip_gdrive_fp, "clean/reg_spec_sites.csv")),
                                   row.names = FALSE),


  # ----------------------------------------------------------------------------
  # # Step 7)  run regressions 
  # ---------------------------------------------------------------------------
  # reg_results_df = run_regressions(clean_pm_spec_df),
  # # save out
  # reg_results_df_out = write.csv(reg_results_df,
  #                                  file_out(!!file.path(wip_gdrive_fp, "results/regression_results.csv")),
  #                                  row.names = FALSE),
  
  # ----------------------------------------------------------------------------
  # # Step 8)  create descriptive figures + figures for paper
  # ---------------------------------------------------------------------------
  # create time series plots for each species
  # conc_rel2avg_df = target(
  #   create_indiv_species_time_series_plots(clean_pm_spec_df))
  
  # create plots demonstrating how much wildfire drives concentration of a given species
  # create_attributable_frac_plots(clean_pm_spec_sites_df)
    
    
  
) # end plan

# make plan
if(file.exists("scripts/.drake")){
  drake_cache = drake::drake_cache("scripts/.drake")
} else {
  drake_cache = drake::new_cache("scripts/.drake")
}
make(wildfire_plan, cache = drake_cache)

