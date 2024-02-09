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
    mutate(year = paste0('20', year))

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
  dplyr::select(CAS.Number = 'CAS_no', IUR_m3_ug, RfCi_mg_m3) %>% 
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
    Parameter == "Lead PM2.5 LC"~ 'PB')) %>%
  # merge CAS numbers into tox thresholds
  left_join(tox_thresholds, by = 'CAS.Number') %>% 
  dplyr::select(-RfCi_mg_m3)


# parellelize
future::plan(multisession)

# list files in the folder path + read in
# spec_fps <- list.files(file.path(data_fp, 'clean'), 
#                       pattern =  '(PB|AS)_gridded_preds.fst',
#                       full.name = TRUE)
# 
# # read in predicted concentration for every grid cell for every day in sample
# pred_conc_df <- read.fst(spec_fps[1]) %>%
#   dplyr::select(grid_id_10km, Date, AS_pred_conc = 'pred_grid_conc', year) %>% 
#   left_join(read.fst(spec_fps[2]) %>% 
#               dplyr::select(grid_id_10km, Date, PB_pred_conc = 'pred_grid_conc', year),
#             by = c('grid_id_10km', 'Date', 'year'))
# this file is huge so saved out and read in - instead of merging
pred_conc_df <- read.fst(file.path(data_fp, 'clean/predicted_gridded_AS_PB_combined.fst')) 

# EXPOSURE ESTIMATES -------------------------------------------------------
# calculate excess cancer risk  
# (using guidance from: https://www.atsdr.cdc.gov/pha-guidance/resources/ATSDR-EDG-Inhalation-508.pdf)
# Cancer risk = IUR * air concentration * EF
# - where IUR is taken from EPA or CA OEHHA
# - EF for cancer = (hrs/day * days/week * years) / (24hrs/day * 7 days/week * 70 yrs)
# - air conc = avg over whatever time period of interest
EF = 8 / (24 * 365 * 70) # this days fraction of your life

# extract current species IUR
AS_IUR <- IUR_vals %>%
  filter(ParamCode == 'AS') %>%
  distinct(IUR_m3_ug) %>%
  as.numeric()

PB_IUR <- IUR_vals %>%
  filter(ParamCode == 'PB') %>%
  distinct(IUR_m3_ug) %>%
  as.numeric()

# calculate the risk using the formula
gridded_risk_df <- pred_conc_df %>% 
      mutate(year = as.character(year)) %>% # each file is one species predictions
  mutate(AS_daily_cancer_risk = AS_pred_conc * AS_IUR * EF,
         PB_daily_cancer_risk = PB_pred_conc * PB_IUR * EF) 

rm(pred_conc_df) # drop to save memory

    # aggregate risk to count up each day to get total risk in a grid cell
    # aggregated_annual_gridded_risk_df <- gridded_risk_df %>% 
    # # aggregate the risk
    #   #mutate(num_excess_canc_cases_in_grid = daily_excess_cancer_risk*pop_ct) %>% 
    #   group_by(grid_id_10km, year) %>% 
    #   dplyr::summarise(AS_total_cancer_risk_in_grid = sum(AS_daily_cancer_risk, na.rm=TRUE),
    #                    AS_avg_cancer_risk_in_grid = mean(AS_daily_cancer_risk, na.rm=TRUE),
    #                    PB_total_cancer_risk_in_grid = sum(PB_daily_cancer_risk, na.rm=TRUE),
    #                    PB_avg_cancer_risk_in_grid = mean(PB_daily_cancer_risk, na.rm=TRUE)) %>% 
    #   ungroup() 
    
    # calculate for each year how many excess cases there could have been if exposed for a lifetime
    aggregated_annual_risk_df <- gridded_risk_df %>% 
      group_by(year) %>% 
      dplyr::summarise(AS_annual_risk = sum(AS_daily_cancer_risk, na.rm=TRUE),
                       PB_annual_risk = sum(PB_daily_cancer_risk, na.rm=TRUE)) %>% 
      ungroup() %>% 
      left_join(us_annual_pop, by = 'year') %>% 
      dplyr::select(year, pop = 'Population', PB_annual_risk, AS_annual_risk) %>% 
      pivot_longer(cols = c('PB_annual_risk', 'AS_annual_risk'), names_to = 'species', values_to = 'risk') %>% 
      mutate(species = str_remove(species, '_annual_risk')) %>% 
      mutate(us_excess_cases = risk*pop) %>% 
      mutate(cdc_cancer_risk = case_when(
        risk >= 1E-4 ~ 'A concern for increased cancer risk',
        risk > 1E-6 & risk < 1E-4 ~ 'Likely a concern for increased risk',
        risk <= 1E-6 ~ 'No concern for increased cancer risk')) %>% 
      mutate(mon_yr = as.Date(paste0(year,"-01-01")))
    
   
    
    datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "48 month")
    
    # PLOT THE TIME SERIES OF EXCESS CANCER CASES
    time_series_cases <- ggplot(aggregated_annual_risk_df,
                                aes(x = mon_yr, 
                                    y = us_excess_cases, color = species)) +
      geom_line() +
      scale_color_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8')) + 
      scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 20000),
                    labels = c("0", "10", "100", "1,000", "10,000", "20,000")) +  
      labs(title = '',
           x = "Year",
           y = "Estimated excess cancer cases") +
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
    time_series_cases
    
    ggsave(
      filename = paste0('Fig5B_excess_cancer_cases_time_series_PB+AS.pdf'),
      plot = time_series_cases,
      path = file.path(results_fp, 'figures/Fig5'),
      scale = 1,
      width = 7,
      height = 6,
      dpi = 320)

    time_series_risk <- ggplot(aggregated_annual_risk_df,
                                aes(x = mon_yr, 
                                    y = risk, color = species)) +
      geom_line() +
      facet_wrap(~species, scales = 'free_y') +
      scale_color_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8')) + 
      scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
      # scale_y_continuous(breaks = c( 0.000001, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008),
      #                    labels = c(
      #                      ".000001" = "1 in 1,000,000", 
      #                      ".00001" = "1 in 100,000", 
      #                      ".00002" = "2 in 100,000", 
      #                      ".00003" = "3 in 100,000",
      #                      ".00004" = "4 in 100,000",
      #                      ".00005" = "5 in 100,000",
      #                      ".00006" = "6 in 100,000",
      #                      ".00007" = "7 in 100,000",
      #                      ".00008" = "8 in 100,000"
      #                    ))+ # Specify custom labels for y-axis
      labs(title = '',
           x = "Year",
           y = "Cumulative estimated cancer risk") +
      theme_minimal() + 
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),
        # axis.text = element_blank(),  # Remove axis labels
        # axis.title = element_text(),  # Remove axis titles
        legend.title = element_text(),
        plot.title = element_text(face = "bold"),
        legend.key.width = unit(.4, "cm"),  # Set the legend key width
        legend.key.height = unit(.4, "cm")  # Set the legend key height
      ) +
      guides(color = 'none')
    time_series_risk
    
    
     # merge with states
     regional_cancer_risk <- gridded_risk_df %>% 
       left_join(grid_in_states_regions, relationship = "many-to-many", by = 'grid_id_10km') %>% 
       group_by(region, year) %>% 
       dplyr::summarise(tot_regional_AS_risk = sum(AS_daily_cancer_risk, na.rm = TRUE),
                        tot_regional_PB_risk = sum(PB_daily_cancer_risk, na.rm = TRUE)) %>% 
       ungroup() 
     
     
     # get regional pop + merge
     regional_cancer_risk_wide <- regional_cancer_risk %>% 
       pivot_longer(cols = c('tot_regional_AS_risk', 'tot_regional_PB_risk'), 
                    names_to = 'species', values_to = 'risk') %>% 
       mutate(species = str_remove(species, 'tot_regional_')) %>% 
       mutate(species = str_remove(species, '_risk')) 
     
     datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "48 month")
     
     # BAR PLOT -------------------------
     bar_plot_risk <- ggplot(regional_cancer_risk_wide %>% 
                               mutate(mon_yr = as.Date(paste0(year,"-01-01"))),
                             aes(x = mon_yr, 
                                 y = risk, 
                                 color = region,
                                 fill = region)) +
       # geom_line() +
       geom_bar(stat = "identity") +
       scale_color_manual(values = reg_pal) +
       scale_fill_manual(values = reg_pal) +
       facet_wrap(~species) +
       scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
       scale_y_continuous(breaks = c(0.000001, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008),
                          labels = c(
                            ".000001" = "1 in 1,000,000", 
                            ".00001" = "1 in 100,000", 
                            ".00002" = "2 in 100,000", 
                            ".00003" = "3 in 100,000",
                            ".00004" = "4 in 100,000",
                            ".00005" = "5 in 100,000",
                            ".00006" = "6 in 100,000",
                            ".00007" = "7 in 100,000",
                            ".00008" = "8 in 100,000"
                          ))+ # Specify custom labels for y-axis
       labs(x = "Year",
            y = "Lifetime Cancer Risk") +
       theme_minimal() +
       theme(
         panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank(),  # Remove minor grid lines
         axis.line = element_line(color = "black"),
         legend.title = element_text(),
         plot.title = element_text(face = "bold"),
         legend.key.width = unit(.4, "cm"),  # Set the legend key width
         legend.key.height = unit(.4, "cm")  # Set the legend key height
       ) 
     bar_plot_risk
     
     ggsave(
       filename = paste0('Fig5B_excess_cancer_risk_bar_plot_AS+PB.pdf'),
       plot = bar_plot_risk,
       path = file.path(results_fp, 'figures/Fig5'),
       scale = 1,
       width = 12,
       height = 6,
       dpi = 320)
     
     
     
     
     avg_risk <- gridded_risk_df %>% 
       group_by(year, grid_id_10km) %>% 
       dplyr::summarise(AS_grid_avg_conc = mean(AS_pred_conc, na.rm= TRUE), 
                        PB_grid_avg_conc = mean(PB_pred_conc, na.rm= TRUE)) %>% 
       ungroup()
     
     EF_yr = 8 * 365 / (24 * 365 * 70)
     
     summed_risk <- avg_risk %>% 
       left_join(all_pop_df) %>% 
       mutate(AS_cancer_cases = EF_yr*AS_grid_avg_conc*AS_IUR*pop_ct) %>% 
       mutate(PB_cancer_cases = EF_yr*AS_grid_avg_conc*PB_IUR*pop_ct) %>% 
       group_by(year) %>% 
       mutate(cases_AS = sum(AS_cancer_cases),
              cases_PB = sum(PB_cancer_cases)) %>% 
       ungroup()
     
     year_df <- summed_risk %>% 
       distinct(year, cases_AS, cases_PB)
     

  # future::plan(NULL)
  
  return(annual_gridded_species_exposure_df)
}



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
  




