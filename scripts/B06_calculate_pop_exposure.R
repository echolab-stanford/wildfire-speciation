# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 11, 2024
# Description: create population weighted exposure


# loadd(c(reg_pal, us_region_map), cache = drake_cache)
# csn_param_fp = file.path(data_fp, 'raw/CSN_parameters.csv')
# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# tox_fp = file.path(data_fp, 'raw/RSL Table/oehha_ca_IURs.csv')
# st_pop_fp <- file.path(data_fp, 'raw/NST-EST2023-POP.xlsx')
# us_states_fp = file.path(boundaries_fp, 'all_national_states.rds')
# us_pop_fp = file.path(data_fp, 'raw/united-states-population-2024-02-08.csv')


calculate_pop_exposure <- function(csn_param_fp, grid_fp, tox_fp, current_species) {

# Read in grid + transform
grid_10km <- st_read(grid_fp) %>% 
  dplyr::select(grid_id_10km = 'ID')

# see which grid cells are in which states
grid_in_states_regions <- st_join(us_region_map, grid_10km) %>% 
  st_drop_geometry() 

# -----------------------------------
# read in population data from GEE
# get list of files to read in:
file_list <- list.files(file.path(gee_data_fp), pattern = '_grid10km.csv', full.names = TRUE)
# current_fp <- file_list[1]

all_pop_df <- map_df(file_list, function(current_fp) {

  current_pop_df <- read.csv(current_fp) %>%
    # clean up data
    dplyr::select(ID, pop_ct = 'sum') %>%
    mutate(year = str_sub(current_fp, -15, -14)) %>%
    mutate(year = paste0('20', year)) %>% 
    distinct()

}) %>%
  bind_rows() %>%
  rename(grid_id_10km = 'ID') %>%
  distinct()


# read in annual US population data
us_annual_pop <- read.csv(us_pop_fp) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%y")) %>% 
  mutate(year = as.character(year(date)))


# GET TOXICITY INFORMATION ----------------------------------------------------
# read in the IURs from CALIFORNIA, EPA, ETC
tox_thresholds <- read_csv(tox_fp) %>%
  dplyr::select(CAS.Number = 'CAS_no', IUR_m3_ug, RfCi_mg_m3, SL_residential_cancer) %>% 
  mutate(RfCi_ug_m3 = 1000*RfCi_mg_m3)

# now adjust to only keep the three relevant parameters
IUR_vals <- read.csv(csn_param_fp) %>%
  filter(CAS.Number %in% tox_thresholds$CAS.Number) %>%
  mutate(flag = ifelse(
    str_detect(Parameter, "2.5"), 'keep', 'drop'
  )) %>%
  filter(flag == 'keep') %>%
  mutate(flag = ifelse(
    str_detect(Parameter, "10"), 'drop', 'keep'
  )) %>%
  filter(flag == 'keep') %>%
  mutate(flag = ifelse(
    str_detect(Parameter, "(STP|TSP)"), 'drop', 'keep'
  )) %>%
  filter(flag == 'keep') %>%
  mutate(flag = ifelse(
    str_detect(Parameter, "nm"), 'drop', 'keep'
  )) %>%
  filter(flag == 'keep') %>%
  dplyr::select(AQSParamCode = 'Parameter.Code', Parameter, CAS.Number) %>%
  mutate(ParamCode = case_when(
    Parameter == "Arsenic PM2.5 LC" ~ 'AS',
    Parameter == "Nickel PM2.5 LC" ~ 'NI',
    Parameter == "Lead PM2.5 LC"~ 'PB')) %>%
  # merge CAS numbers into tox thresholds
  left_join(tox_thresholds, by = 'CAS.Number') %>% 
  dplyr::select(-RfCi_mg_m3)


# parellelize
future::plan(multisession)

# list files in the folder path + read in
spec_fps <- list.files(file.path(data_fp, 'clean'),
                      pattern =  '(PB|AS|NI)_gridded_preds.fst',
                      full.name = TRUE)

#----------------------------------------------------------------------
# # read in predicted concentration for every grid cell for every day in sample
# pred_conc_df <- read.fst(spec_fps[1]) %>%
#   dplyr::select(grid_id_10km, Date, AS_pred_conc = 'pred_grid_conc', year) %>%
#   left_join(read.fst(spec_fps[2]) %>%
#               dplyr::select(grid_id_10km, Date, NI_pred_conc = 'pred_grid_conc', year),
#             by = c('grid_id_10km', 'Date', 'year')) %>% 
#   left_join(read.fst(spec_fps[3]) %>%
#               dplyr::select(grid_id_10km, Date, PB_pred_conc = 'pred_grid_conc', year),
#             by = c('grid_id_10km', 'Date', 'year'))
# write.fst(pred_conc_df, file.path(data_fp, 'clean/predicted_gridded_AS_PB_NI_combined.fst'))
#----------------------------------------------------------------------

# this file is huge so saved out and read in - instead of merging
pred_conc_df <- read.fst(file.path(data_fp, 'clean/predicted_gridded_AS_PB_NI_combined.fst')) %>% 
  mutate(across(where(is.numeric), ~replace(., . < 0, NA))) %>%
  mutate(month = month(Date))

# CANCER RISK ASSESSMENT ---------------------------------------------
# get annual avg and seasonal averages:
avg_ann_conc_df <- pred_conc_df %>%
  group_by(grid_id_10km, year) %>% 
  dplyr::summarise(avg_annual_PB = mean(PB_pred_conc, na.rm=TRUE),
                   avg_annual_NI = mean(NI_pred_conc, na.rm=TRUE),
                   avg_annual_AS = mean(AS_pred_conc, na.rm=TRUE)) %>% 
  ungroup() 

# avg_seasonal_conc_df <- pred_conc_df %>% 
#   mutate(wf_season = ifelse(month %in% c(6:11), "wf season", 'non wf season')) %>% 
#   group_by(grid_id_10km, year, wf_season) %>% 
#   dplyr::summarise(avg_seasonal_PB = mean(PB_pred_conc, na.rm=TRUE),
#                    avg_seasonal_NI = mean(NI_pred_conc, na.rm=TRUE),
#                    avg_seasonal_AS = mean(AS_pred_conc, na.rm=TRUE)) %>%
#   ungroup() 

# EXPOSURE ESTIMATES -------------------------------------------------------
# calculate excess cancer risk  
# (using guidance from: https://www.atsdr.cdc.gov/pha-guidance/resources/ATSDR-EDG-Inhalation-508.pdf)

# Population risk: Cancer burden = IUR * air concentration * EF
# - where IUR is taken from EPA or CA OEHHA
# - EF for cancer = (24 hrs/day * 7 days/week * 52.15 weeks/year) * 30 year / (24 hrs/day * 7 days/week * 52.15 weeks/year * 70 yrs)
# - air conc = avg conc in one year or one wildfire season
# - divide by lifetime to get annual attributable number of future cancers
# - for cancer at population level use 70 year duration
# lifetime cases = (1 - fraction of time indoors) * conc of outdoor air * IUR * population

# "Note that a 70-year exposure duration is required to estimate cancer burden or provide an estimate of population-wide risk"
# https://oehha.ca.gov/media/downloads/crnr/2015guidancemanual.pdf
# EF for cancer = (24 hrs/day * 7 days/week * 52.15 weeks/year) * 70 year / (24 hrs/day * 7 days/week * 52.15 weeks/year * 70 yrs)
EF_cancer =  (24 * 350 * 70) / (24 * 365 * 70)

# extract current species IUR
AS_IUR <- IUR_vals %>%
  filter(ParamCode == 'AS') %>%
  filter(!is.na(IUR_m3_ug)) %>% 
  distinct(IUR_m3_ug) %>%
  as.numeric()

PB_IUR <- IUR_vals %>%
  filter(ParamCode == 'PB') %>%
  filter(!is.na(IUR_m3_ug)) %>% 
  distinct(IUR_m3_ug) %>%
  as.numeric()

NI_IUR <- IUR_vals %>%
  filter(ParamCode == 'NI') %>%
  filter(!is.na(IUR_m3_ug)) %>% 
  distinct(IUR_m3_ug) %>%
  as.numeric()

# calculate the annual risk using the formula
gridded_ann_risk_df <- avg_ann_conc_df %>% 
  mutate(year = as.character(year)) %>% 
  left_join(all_pop_df, by = c('year', 'grid_id_10km')) %>% 
  mutate(AS_lifetime_cancer_risk = avg_annual_AS * AS_IUR * EF_cancer * pop_ct,
         NI_lifetime_cancer_risk = avg_annual_NI * NI_IUR * EF_cancer * pop_ct,
         PB_lifetime_cancer_risk = avg_annual_PB * PB_IUR * EF_cancer * pop_ct) %>% 
  group_by(year) %>% 
  dplyr::summarise(AS_canc_burden = sum(AS_lifetime_cancer_risk, na.rm = TRUE),
                   NI_canc_burden = sum(NI_lifetime_cancer_risk, na.rm = TRUE),
                   PB_canc_burden = sum(PB_lifetime_cancer_risk, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(mon_yr = as.Date(paste0(year,"-01-01"))) %>% 
  pivot_longer(cols = c('AS_canc_burden', 'NI_canc_burden', 'PB_canc_burden'), 
               names_to = 'species', values_to = 'canc_burden') %>% 
  mutate(species = str_remove(species, '_canc_burden')) 

datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")        
                                
# calculate the seasonal risk using the formula
# gridded_ssnl_risk_df <- avg_seasonal_conc_df %>% 
#   mutate(year = as.character(year)) %>% 
#   left_join(all_pop_df, by = c('year', 'grid_id_10km')) %>% 
#   mutate(AS_lifetime_cancer_risk = avg_seasonal_AS * AS_IUR * EF_cancer * pop_ct,
#          NI_lifetime_cancer_risk = avg_seasonal_NI * NI_IUR * EF_cancer * pop_ct,
#          PB_lifetime_cancer_risk = avg_seasonal_PB * PB_IUR * EF_cancer * pop_ct) %>% 
#   group_by(wf_season, year) %>% 
#   dplyr::summarise(AS_canc_burden_ssnl = sum(AS_lifetime_cancer_risk, na.rm = TRUE),
#                    NI_canc_burden_ssnl = sum(NI_lifetime_cancer_risk, na.rm = TRUE),
#                    PB_canc_burden_ssnl = sum(PB_lifetime_cancer_risk, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   mutate(mon_yr = as.Date(paste0(year,"-01-01"))) %>% 
#   pivot_longer(cols = c('AS_canc_burden_ssnl', 'PB_canc_burden_ssnl', 'NI_canc_burden_ssnl'), 
#                names_to = 'species', values_to = 'canc_burden') %>% 
#   mutate(species = str_remove(species, '_canc_burden_ssnl')) 
# 
#     
#     
#     # PLOT THE TIME SERIES OF EXCESS CANCER CASES BY WILDFIRE SEASON
#     time_series_cases <- ggplot(gridded_ssnl_risk_df %>% 
#                                   filter(wf_season == 'wf season'),
#                                 aes(x = mon_yr, 
#                                     y = canc_burden, color = species)) +
#       geom_line(linewidth = 2) +
#       scale_color_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) + 
#       scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#       labs(title = 'Population-level excess cancer burden',
#            x = "Year",
#            y = "Excess cancer burden") +
#       theme_minimal() + 
#       theme(panel.grid.major = element_blank(),  # Remove major grid lines
#             panel.grid.minor = element_blank(),  # Remove minor grid lines
#             axis.line = element_line(color = "black"),
#             legend.title = element_text(),
#             plot.title = element_text(face = "bold"),
#             legend.key.width = unit(.4, "cm"),  # Set the legend key width
#             legend.key.height = unit(.4, "cm")  # Set the legend key height
#       ) #+
#       #guides(color = 'none')
#     time_series_cases
    
    # ggsave(
    #   filename = paste0('Fig5B_excess_cancer_cases_time_series_wildfire_season_PB+AS+NI.pdf'),
    #   plot = time_series_cases,
    #   path = file.path(results_fp, 'figures/Fig5'),
    #   scale = 1,
    #   width = 7,
    #   height = 6,
    #   dpi = 320)

    
    # ANNUAL RISK
    time_series_ann_burden <- ggplot(gridded_ann_risk_df,
                                aes(x = mon_yr, 
                                    y = canc_burden, color = species)) +
      geom_line(linewidth = 1.5) +
      scale_color_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) +
      scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
      labs(title = 'Population-level excess cancer burden',
           x = "Year",
           y = "Excess cancer burden") +
      theme_minimal() + 
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),
        legend.title = element_text(),
        plot.title = element_text(face = "bold"),
        legend.key.width = unit(.4, "cm"),  # Set the legend key width
        legend.key.height = unit(.4, "cm")  # Set the legend key height
      ) #+
    #guides(color = 'none')
    time_series_ann_burden
    
    ggsave(
      filename = paste0('Fig5B_excess_cancer_cases_time_series_annual_PB+AS+NI.pdf'),
      plot = time_series_ann_burden,
      path = file.path(results_fp, 'figures/Fig5'),
      scale = 1,
      width = 7,
      height = 4,
      dpi = 320)
    
    
    # CANCER SCREENING LEVEL ASSESMENT ---------------------------------------------
    # get number of days above screening levels
    # extract current species IUR
    AS_SL <- IUR_vals %>%
      filter(ParamCode == 'AS') %>%
      distinct(SL_residential_cancer) %>%
      filter(!is.na(SL_residential_cancer)) %>% 
      as.numeric()
    
    PB_SL <- IUR_vals %>%
      filter(ParamCode == 'PB') %>%
      distinct(SL_residential_cancer) %>%
      filter(!is.na(SL_residential_cancer)) %>% 
      as.numeric()
    
    NI_SL <- IUR_vals %>%
      filter(ParamCode == 'NI') %>%
      distinct(SL_residential_cancer) %>%
      filter(!is.na(SL_residential_cancer)) %>% 
      as.numeric()
    
    # calculate if a grid cell day is above the SL
    grid_days_above_SL_df <- pred_conc_df %>%
      mutate(AS_above = ifelse(AS_pred_conc >= AS_SL, TRUE, FALSE),
             PB_above = ifelse(PB_pred_conc >= PB_SL, TRUE, FALSE),
             NI_above = ifelse(NI_pred_conc >= NI_SL, TRUE, FALSE)) %>% 
      group_by(grid_id_10km, year) %>% 
      dplyr::summarise(AS_grid_ann_above = sum(AS_above, na.rm = TRUE),
                       PB_grid_ann_above = sum(PB_above, na.rm = TRUE),
                       NI_grid_ann_above = sum(NI_above, na.rm = TRUE)) %>% 
      ungroup()
    
    # get annual summary
    ann_sum_grid_days <- grid_days_above_SL_df %>% 
      dplyr::select(-PB_grid_ann_above, -NI_grid_ann_above) %>% 
      mutate(year = as.character(year)) %>% 
      left_join(all_pop_df, by = c('year', 'grid_id_10km')) %>% 
      mutate(numerator = pop_ct*AS_grid_ann_above,
             denominator =pop_ct*365) %>% 
      group_by(year) %>% 
      dplyr::summarise(numerator_sum = sum(numerator, na.rm = TRUE),
                       denominator_smm = sum(denominator, na.rm = TRUE)) %>% 
      ungroup() %>% 
      mutate(mon_yr = as.Date(paste0(year,"-01-01"))) %>% 
      mutate(pct_pop_days_above = 100*(numerator_sum/denominator_smm))
    
    
    
    time_series_days_above <- ggplot(ann_sum_grid_days,
                                     aes(x = mon_yr, 
                                         y = pct_pop_days_above, color = species)) +
      geom_line(linewidth = 2) +
      scale_color_manual(values = c('AS' = '#E2D200')) +
      scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
      labs(title = ' cancer risk screening level',
           x = "Year",
           y = "% of population-days above threshold") +
      theme_minimal() + 
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),
        legend.title = element_text(),
        plot.title = element_text(face = "bold"),
        legend.key.width = unit(.4, "cm"),  # Set the legend key width
        legend.key.height = unit(.4, "cm")  # Set the legend key height
      ) #+
    #guides(color = 'none')
    time_series_days_above
    
    ggsave(
      filename = paste0('Fig5B_AS_days_above_screening_level.pdf'),
      plot = time_series_days_above,
      path = file.path(results_fp, 'figures/Fig5'),
      scale = 1,
      width = 7,
      height = 3,
      dpi = 320)
    

  # future::plan(NULL)
  
  return(annual_gridded_species_exposure_df)
}



# calculate for each year how many excess cases there could have been if exposed for a lifetime
# aggregated_annual_risk_df <- gridded_risk_df %>% 
#   group_by(year) %>% 
#   dplyr::summarise(AS_annual_risk = sum(AS_daily_cancer_risk, na.rm=TRUE),
#                    PB_annual_risk = sum(PB_daily_cancer_risk, na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   left_join(us_annual_pop, by = 'year') %>% 
#   dplyr::select(year, pop = 'Population', PB_annual_risk, AS_annual_risk) %>% 
#   pivot_longer(cols = c('PB_annual_risk', 'AS_annual_risk'), names_to = 'species', values_to = 'risk') %>% 
#   mutate(species = str_remove(species, '_annual_risk')) %>% 
#   mutate(us_excess_cases = risk*pop) %>% 
#   mutate(cdc_cancer_risk = case_when(
#     risk >= 1E-4 ~ 'A concern for increased cancer risk',
#     risk > 1E-6 & risk < 1E-4 ~ 'Likely a concern for increased risk',
#     risk <= 1E-6 ~ 'No concern for increased cancer risk')) %>% 
#   mutate(mon_yr = as.Date(paste0(year,"-01-01")))

# 
# 
# regional_cancer_risk <- gridded_risk_df %>% 
#   left_join(grid_in_states_regions, relationship = "many-to-many", by = 'grid_id_10km') %>% 
#   group_by(region, year) %>% 
#   dplyr::summarise(tot_regional_AS_risk = sum(AS_daily_cancer_risk, na.rm = TRUE),
#                    tot_regional_PB_risk = sum(PB_daily_cancer_risk, na.rm = TRUE)) %>% 
#   ungroup() 
# 
# 
# # get regional pop + merge
# regional_cancer_risk_wide <- regional_cancer_risk %>% 
#   pivot_longer(cols = c('tot_regional_AS_risk', 'tot_regional_PB_risk'), 
#                names_to = 'species', values_to = 'risk') %>% 
#   mutate(species = str_remove(species, 'tot_regional_')) %>% 
#   mutate(species = str_remove(species, '_risk')) 
# 
# datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "48 month")
# 
# # BAR PLOT -------------------------
# bar_plot_risk <- ggplot(regional_cancer_risk_wide %>% 
#                           mutate(mon_yr = as.Date(paste0(year,"-01-01"))),
#                         aes(x = mon_yr, 
#                             y = risk, 
#                             color = region,
#                             fill = region)) +
#   # geom_line() +
#   geom_bar(stat = "identity") +
#   scale_color_manual(values = reg_pal) +
#   scale_fill_manual(values = reg_pal) +
#   facet_wrap(~species) +
#   scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#   scale_y_continuous(breaks = c(0.000001, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008),
#                      labels = c(
#                        ".000001" = "1 in 1,000,000", 
#                        ".00001" = "1 in 100,000", 
#                        ".00002" = "2 in 100,000", 
#                        ".00003" = "3 in 100,000",
#                        ".00004" = "4 in 100,000",
#                        ".00005" = "5 in 100,000",
#                        ".00006" = "6 in 100,000",
#                        ".00007" = "7 in 100,000",
#                        ".00008" = "8 in 100,000"
#                      ))+ # Specify custom labels for y-axis
#   labs(x = "Year",
#        y = "Lifetime Cancer Risk") +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     axis.line = element_line(color = "black"),
#     legend.title = element_text(),
#     plot.title = element_text(face = "bold"),
#     legend.key.width = unit(.4, "cm"),  # Set the legend key width
#     legend.key.height = unit(.4, "cm")  # Set the legend key height
#   ) 
# bar_plot_risk
# 
# ggsave(
#   filename = paste0('Fig5B_excess_cancer_risk_bar_plot_AS+PB.pdf'),
#   plot = bar_plot_risk,
#   path = file.path(results_fp, 'figures/Fig5'),
#   scale = 1,
#   width = 12,
#   height = 6,
#   dpi = 320)
# 
# 
# 
# 



  # 
  # 
  # #1)  calculate weighted average by year + region to make time series
  # weighted_avg_region_yr <- full_pop_conc_df %>%
  #   mutate(month = month(Date)) %>%
  #   # for a given region with a grid cell in a year, region, and species, what is the pop weighted avg exposure
  #   group_by(year, region, species) %>%
  #   # if weighted avg conc is high -> conc was high where there are more people
  #   dplyr::summarise(region_yr_wgt_exposure = weighted.mean(pred_grid_conc, pop_ct, na.rm = TRUE)) %>%
  #   ungroup()
  # 
  
     # for each grid cell get avg risk and total risk across all years
     # aggregated_tot_gridded_risk_df <- aggregated_annual_gridded_risk_df %>%
     #   group_by(grid_id_10km) %>% 
     #   dplyr::summarise(AS_total_cancer_risk_in_grid = sum(AS_total_cancer_risk_in_grid, na.rm=TRUE),
     #                   AS_avg_cancer_risk_in_grid = mean(AS_avg_cancer_risk_in_grid, na.rm=TRUE),
     #                   PB_total_cancer_risk_in_grid = sum(PB_total_cancer_risk_in_grid, na.rm=TRUE),
     #                   PB_avg_cancer_risk_in_grid = mean(PB_avg_cancer_risk_in_grid, na.rm=TRUE)) %>% 
     #   ungroup() %>% 
     #   pivot_longer(cols = c('AS_total_cancer_risk_in_grid', 'PB_total_cancer_risk_in_grid'), 
     #                names_to = 'species', values_to = 'tot_risk') %>% 
     #   mutate(species = str_remove(species, "_total_cancer_risk_in_grid")) %>% 
     #   left_join(all_pop_df %>% 
     #               filter(year == 2020), by = 'grid_id_10km') %>% 
     #   mutate(cases = pop_ct *tot_risk)
     
     # MAKE A MAP!!
     # cumulative_cancer_risk_map <- ggplot() +
     #   geom_sf(data = aggregated_tot_gridded_risk_df %>% 
     #             filter(species == 'AS') %>% 
     #             left_join(grid_10km, by ='grid_id_10km') %>%
     #             st_as_sf(),
     #           aes(color = tot_risk,
     #               fill = tot_risk)) +
     #   scale_color_distiller(type = "seq",
     #                         palette = "Greys") +
     # scale_fill_distiller(type = "seq",
     #                       palette = "Greys") +
     #   # facet_wrap(~species) +
     #   theme_minimal() +  # Apply a minimal theme
     #   theme(
     #     panel.grid.major = element_blank(),  # Remove major grid lines
     #     panel.grid.minor = element_blank(),  # Remove minor grid lines
     #     axis.text = element_blank(),  # Remove axis labels
     #     axis.title = element_text(),  # Remove axis titles
     #     legend.title = element_text(),
     #     plot.title = element_text(face = "bold"),
     #     legend.key.width = unit(.4, "cm"),  # Set the legend key width
     #     legend.key.height = unit(.4, "cm")  # Set the legend key height
     #   ) +
     #   labs(title = "",
     #        color = paste0('')) +
     #   guides(fill = 'none')
     # cumulative_cancer_risk_map
  




