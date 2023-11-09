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
data_fp = '/Users/ekrasovich/Library/CloudStorage/Dropbox-BurkeLab/data/wildfire_speciation'
gdrive_data_fp = '/Users/ekrasovich/Library/CloudStorage/GoogleDrive-emmars@stanford.edu/Shared drives/echolab:data'
results_fp = '/Users/ekrasovich/Library/CloudStorage/Dropbox-BurkeLab/projects/ws_project'
# get rid of scientific notation
options(scipen = 999)


# SOURCE FUNCTIONS:
list.files(here::here(), full.names=T, pattern = '(A|B)') %>%
  str_subset("\\.R$") %>%
  str_subset('workflow', negate = TRUE) %>%
  str_subset('B06compare_high_smoke_to_low_smoke_days', negate = TRUE) %>%
  str_subset('combine_speciation_w_CASno_and_RSLs', negate = TRUE) %>%
  str_subset('B04create_attributable_frac_timeseries_plots', negate = TRUE) %>%
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
                                improve_params_fp = file_in(!!file.path(data_fp, 'raw/IMPROVE_parameters.csv')))
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
    CONUS_spec_sites_df_out = write.csv(CONUS_spec_sites_df,
                                         file_out(!!file.path(
                                           data_fp, "intermediate/CONUS_combined_speciation_sites_2006_2020.csv"))),
  
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
  
  # # --------------------------------------------------------------------------------
  # # Step 5)  merge speciation sites with smoke plume data
  # # --------------------------------------------------------------------------------
  # IMPROVEspec_smoke_plumes_df = target(
  #   merge_sites_w_smoke_plumes(
  #     smoke_plumes_fp = file_in(!!file.path(data_fp, 'intermediate/hms_smoke_plumes.rds')),
  #     spec_w_smoke_pm_df,
  #     CONUS_spec_sites_df %>% filter(Dataset == 'IMPROVE'))),
  # 
  # CSNspec_smoke_plumes_df = target(
  #   merge_sites_w_smoke_plumes(
  #     smoke_plumes_fp = file_in(!!file.path(data_fp, 'intermediate/hms_smoke_plumes.rds')),
  #     spec_w_smoke_pm_df,
  #     CONUS_spec_sites_df %>% filter(Dataset == 'CSN'))),
  
  
  # # --------------------------------------------------------------------------------
  # # Step 6) final cleaning and prepping data for analysis
  # # --------------------------------------------------------------------------------
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
  parameter_categories = read.csv(file_in(!!file.path(data_fp, 'intermediate/xwalk_species_type.csv'))),
  
  
  # # --------------------------------------------------------------------------
  # # --------------------------------------------------------------------------
  # # PART B: FIGURES + ANALYSIS
  # # --------------------------------------------------------------------------
  # # --------------------------------------------------------------------------
  # set up palette for each species
      spec_pal = c("#35274A", "#B40F20", "#E58601", "#E2D200", "#46ACC8", "#0B775E","#CAB38C", "#899DA4"),
      pm_pal = c( "grey70", "#9B110E"),
      reg_pal = c("#DC863B", "#FAD510","#649373", "#1B5656", "#5A283E", "#F2300F")
  
  # # --------------------------------------------------------------------------------
  # # FIGURE 1
  # # --------------------------------------------------------------------------------
  # # a) Figure 1A - create map of monitoring stations 
  # monitor_sf = target(
  #   create_map_of_monitoring_stations(clean_PMspec_df)),
  
  # b) Figure 1B - create time series of raw PM data
  # raw_PM_time_series = target(
  #   plot_raw_PM_data(clean_PMspec_df, pm_pal, parameter_categories)),
  
  # c) Figure 1C - create time series of monthyl chemical species averages
  # avg_w_spec_type = target(
  #   plot_seasonality_w_monthly_average(clean_PMspec_df, 
  #                                      parameter_categories, 
  #                                      spec_pal)),
  
  # # --------------------------------------------------------------------------------
  # # FIGURE 2
  # # --------------------------------------------------------------------------------
  # 2a) create plots of smoke coefficients - both versions (pct change + log scale of levels)
  # reg_df = target(
  #   create_smoke_coeff_plots(clean_PMspec_df, 
  #                            parameter_categories, 
  #                            spec_pal)),
  # # 
  # # c) Figure 3
  # # create regional map for regional plot
  # us_region_map = target(
  #   create_us_region_map(st_fp = file_in(!!file.path(data_gdrive_fp, 'boundaries/all_national_states.rds')))
  #   ),
  
  
  
) # end plan

# make plan
if(file.exists(".drake")){
  drake_cache = drake::drake_cache(".drake")
} else {
  drake_cache = drake::new_cache(".drake")
}
make(wildfire_plan, cache = drake_cache)


  



  # # ----------------------------------------------------------------------------
  # # # Step 7)  create descriptive figures + figures for paper
  # # ---------------------------------------------------------------------------
  # 
  
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

# parameter_categories = tibble(
#   species = c("AL","AS","BR", "CA", "CL", "CHL","CR", "CU", "EC", "FE", 
#               "K", "MG","MN", "NA", "NI", "NO3", "N2", "OC", "P",  "PB", "RB",
#               "S",  "SE", "SI", "SO4", "SOIL", "SR", "TI", "V",  "ZN", "ZR")) %>% 
#   mutate(spec_type = case_when(
#     species == "AL" ~ 'toxic metal',
#     species == "AS" ~ 'heavy metal',
#     species == "BR" ~ 'toxic nonmetal',
#     species == 'CA' ~ 'non-toxic metal',
#     species == "CL" ~ 'toxic nonmetal',
#     species == 'CHL' ~ 'toxic nonmetal',
#     species == "CR" ~ 'heavy metal',
#     species == "CU" ~ 'toxic metal',
#     species == "EC" ~ 'secondary organic',
#     species == "FE" ~ 'toxic metal',
#     species == "K" ~ 'toxic metal',
#     species == "MG" ~ 'non-toxic metal',
#     species == "MN" ~ 'toxic metal',
#     species == "NA" ~ 'toxic metal',
#     species == "NI" ~ 'toxic metal',
#     species == "NO3" ~ 'secondary inorganic',
#     species == "N2" ~ 'non-toxic nonmetal',
#     species == "OC" ~ 'secondary organic',
#     species == "P" ~ 'toxic nonmetal',
#     species == "PB" ~ 'heavy metal',
#     species == 'RB' ~ 'non-toxic metal',
#     species == "S" ~ 'toxicity potentiator',
#     species == "SE" ~ 'toxic nonmetal',
#     species == "SI" ~ 'non-toxic metal',
#     species == 'SO4' ~ 'secondary inorganic',
#     species == "SOIL" ~ 'non-toxic nonmetal',
#     species == 'SR' ~ 'non-toxic metal',
#     species == "TI" ~ 'toxic metal',
#     species == "V" ~ 'toxic metal',
#     species == "ZN" ~ 'toxic metal',
#     species == "ZR" ~ 'toxic metal'
#   )) %>% 
#   mutate(spec_type = str_to_sentence(spec_type)) %>% 
#   mutate(species_long = case_when(
#     species == "AL" ~ 'Alumnium (Al)',
#     species == "AS" ~ 'Arsenic (As)',
#     species == "BR" ~ 'Bromine (Br)',
#     species == 'CA' ~ 'Calcium (Ca)',
#     species == "CL" ~ 'Chlorine (Cl)',
#     species == 'CHL' ~ 'Chloride (Chl)',
#     species == "CR" ~ 'Chromium (Cr)',
#     species == "CU" ~ 'Copper (Cu)',
#     species == "EC" ~ 'Elemental carbon (EC)',
#     species == "FE" ~ 'Iron (Fe)',
#     species == "K" ~ 'Potassium (K)',
#     species == "MG" ~ 'Magnesium (Mg)',
#     species == "MN" ~ 'Manganese (Mn)',
#     species == "NA" ~ 'Sodium (Na)',
#     species == "NI" ~ 'Nickel (Ni)',
#     species == "NO3" ~ 'Nitrate (NO3)',
#     species == "N2" ~ 'Nitrogen (N2)',
#     species == "OC" ~ 'Organic Carbon (OC)',
#     species == "P" ~ 'Phosphorus (P)',
#     species == "PB" ~ 'Lead (Pb)',
#     species == 'RB' ~ 'Rubidium (Rb)',
#     species == "S" ~ 'Sulfur (S)',
#     species == "SE" ~ 'Selenium (Se)',
#     species == "SI" ~ 'Silicon (Si)',
#     species == 'SO4' ~ 'Sulfate (SO4)',
#     species == "SOIL" ~ 'Soil',
#     species == 'SR' ~ 'Strontium (Sr)',
#     species == "TI" ~ 'Titanium (Ti)',
#     species == "V" ~ 'Vanadium (V)',
#     species == "ZN" ~ 'Zinc (Zn)',
#     species == "ZR" ~ 'Zirconium (Zr)'
#   ))


    


