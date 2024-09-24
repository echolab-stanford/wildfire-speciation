# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 11, 2024
# Description: create population weighted exposure

# 
# loadd(c(region_pal, us_region_map), cache = drake_cache)
# csn_param_fp = file.path(data_fp, 'raw/CSN_parameters.csv')
# grid_fp = file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')
# tox_fp = file.path(data_fp, 'raw/RSL Table/oehha_ca_IURs.csv')
# us_pop_fp = file.path(data_fp, 'raw/united-states-population-2024-02-08.csv')

calculate_risk_from_pop_exposure <- function(csn_param_fp, grid_fp, tox_fp, us_pop_fp, us_region_map, region_pal) {

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

# check annual pop
test_pop <- all_pop_df %>% 
  group_by(year) %>% 
  dplyr::summarise(ann_pop = sum(pop_ct, na.rm =T)) %>% 
  ungroup()

# read in annual US population data
us_annual_pop <- read.csv(us_pop_fp) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%y")) %>% 
  mutate(year = as.character(year(date)))

# join population data from the two sources and calculate a reweighting factor 
# since the grid cells provide too high of a population count
reweight_pop <- test_pop %>% 
  left_join(us_annual_pop %>% 
              dplyr::select(year, census_pop = 'Population'), by = 'year') %>% 
  mutate(reweight_factor = census_pop/ann_pop)

# now reweight the gridded population data
all_pop_reweighted_df <- all_pop_df %>% 
  left_join(reweight_pop %>% 
              dplyr::select(year, reweight_factor), by = 'year') %>% 
  mutate(pop_adj = pop_ct*reweight_factor)

#check pop counts again:
# test_pop2 <- all_pop_rewieghted_df %>% 
#   group_by(year) %>% 
#   dplyr::summarise(ann_pop = sum(pop_adj, na.rm =T)) %>% 
#   ungroup()


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
    Parameter == "Lead PM2.5 LC"~ 'PB'
    )) %>%
  # merge CAS numbers into tox thresholds
  left_join(tox_thresholds, by = 'CAS.Number') %>% 
  dplyr::select(-RfCi_mg_m3)

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

# parellelize
future::plan(multisession)

# this file is huge so saved out and read in - instead of merging
pred_conc_df <- read.fst(file.path(data_fp, 'clean/predicted_gridded_AS_PB_NI_combined.fst')) %>% 
  mutate(across(where(is.numeric), ~replace(., . < 0, NA))) %>%
  mutate(month = month(Date))

# EXPOSURE ASSESSMENT FROM DICKINSON ET AL 2022 ----------------------------------------------------------------------
# EC = CA * ET * EF * ED) / AT
# CANCER RISK = IUR * EC

# CA = ambient concentration
# ET = 24 hours, fire was occurring 24 hours per day
# EF = 365 days (fire event) 
# ED = 30 years (residential scenario) or 70 years (lifetime scenario)
# AT = 24 hours * 365 days * 70 years

# get tot 5 yr conc
tot_5yr_conc_df <- pred_conc_df %>%
  mutate(samp_period = case_when(
    year > 2005 & year < 2011 ~ '2006-2010',
    year > 2010 & year < 2016 ~ '2011-2015',
    year > 2015 & year < 2021 ~ '2016-2020'
  )) %>% 
  group_by(grid_id_10km, samp_period) %>%
  dplyr::summarise(tot_5yr_PB = sum(PB_pred_conc, na.rm=TRUE),
                   tot_5yr_NI = sum(NI_pred_conc, na.rm=TRUE),
                   tot_5yr_AS = sum(AS_pred_conc, na.rm=TRUE)) %>%
  ungroup() %>% 
  left_join(grid_in_states_regions, by = 'grid_id_10km')


# RESIDENTIAL SCENARIO: 26 years
gridded_5yr_cancer_risk_df <- tot_5yr_conc_df %>% 
  left_join(all_pop_reweighted_df %>%
              mutate(samp_period = case_when(
                year > 2005 & year < 2011 ~ '2006-2010',
                year > 2010 & year < 2016 ~ '2011-2015',
                year > 2015 & year < 2021 ~ '2016-2020'
              )) %>% 
              group_by(grid_id_10km, samp_period) %>%
              # take max, assume this is the last year in sample
              dplyr::summarise(pop_pd = max(pop_adj, na.rm = TRUE)) %>% 
              ungroup(), 
            by = c('samp_period', 'grid_id_10km')) %>% 
  # calculate lifetime cancer risk due to exposure to 5 years of wildfire attributable concentrations
  mutate(grid_AS_lifetime_cancer_risk = AS_IUR * tot_5yr_AS * 70 * (1/(70*365)) * pop_pd,
         grid_NI_lifetime_cancer_risk = NI_IUR * tot_5yr_NI *  70 * (1/(70*365)) * pop_pd,
         grid_PB_lifetime_cancer_risk = PB_IUR * tot_5yr_PB * 70 * (1/(70*365)) * pop_pd,
         grid_AS_residential_cancer_risk = AS_IUR * tot_5yr_AS * 30 * (1/(70*365)) * pop_pd,
         grid_NI_residential_cancer_risk = NI_IUR * tot_5yr_NI *  30 * (1/(70*365)) * pop_pd,
         grid_PB_residential_cancer_risk = PB_IUR * tot_5yr_PB * 30 * (1/(70*365)) * pop_pd)

# regional
regional_5yr_cancer_risk <- gridded_5yr_cancer_risk_df %>% 
  group_by(samp_period, region) %>% 
  dplyr::summarise(AS_lifetime_canc_burden = sum(grid_AS_lifetime_cancer_risk, na.rm = TRUE),
                   NI_lifetime_canc_burden = sum(grid_NI_lifetime_cancer_risk, na.rm = TRUE),
                   PB_lifetime_canc_burden = sum(grid_PB_lifetime_cancer_risk, na.rm = TRUE),
                   AS_res_canc_burden = sum(grid_AS_residential_cancer_risk, na.rm = TRUE),
                   NI_res_canc_burden = sum(grid_NI_residential_cancer_risk, na.rm = TRUE),
                   PB_res_canc_burden = sum(grid_PB_residential_cancer_risk, na.rm = TRUE)
                   ) %>% 
  ungroup()

# pivot longer for easier plotting
risk_5yr_long_df <- regional_5yr_cancer_risk %>% 
  pivot_longer(cols = c('AS_lifetime_canc_burden', 'NI_lifetime_canc_burden', 'PB_lifetime_canc_burden'), 
               names_to = 'species', values_to = 'canc_burden') %>% 
  mutate(species = str_remove(species, '_lifetime_canc_burden')) %>% 
  mutate(scenario = 'lifetime') %>% 
  dplyr::select(species, samp_period, region, canc_burden, scenario) %>% 
  bind_rows(regional_5yr_cancer_risk %>% 
              pivot_longer(cols = c('AS_res_canc_burden', 'NI_res_canc_burden', 'PB_res_canc_burden'), 
               names_to = 'species', values_to = 'canc_burden') %>% 
              mutate(species = str_remove(species, '_res_canc_burden')) %>% 
              mutate(scenario = 'residential') %>% 
              dplyr::select(species, samp_period, region, canc_burden, scenario)) %>% 
  mutate(samp_period = fct_relevel(samp_period,
                                   c("2006-2010",
                                   "2011-2015",
                                   "2016-2020")))

# 5yr RISK
risk_scenario_burden <- ggplot(risk_5yr_long_df %>% 
                                 filter(samp_period != "2011-2015") %>% 
                                 filter(scenario == 'lifetime'),
                                 aes(x = samp_period, 
                                     y = canc_burden, color = region, 
                                     group = region)) +
  geom_line(linewidth = .5 ) +
  geom_point(size = 3) +
  scale_color_manual(values=region_pal) +
  facet_wrap(~species, ncol = 1, scales = 'free_y') +
  theme_minimal() + 
theme(
  axis.line = element_line(color = "black"),
  legend.title = element_text(),
  plot.title = element_text(face = "bold"),
  legend.position = "bottom",
  legend.key.width = unit(.4, "cm"),  # Set the legend key width
  legend.key.height = unit(.4, "cm")  # Set the legend key height
) 
risk_scenario_burden

ggsave(
  filename = paste0('Fig5B_excess_cancer_cases_lifetime_PB+AS+NI.pdf'),
  plot = risk_scenario_burden,
  path = file.path(results_fp, 'Fig5'),
  scale = 1,
  width = 5.5,
  height = 7.5,
  dpi = 320)

ggsave(
  filename = paste0('Fig5B_excess_cancer_cases_lifetime_PB+AS+NI.png'),
  plot = risk_scenario_burden,
  path = file.path(results_fp, 'Fig5'),
  scale = 1,
  width = 4.5,
  height = 7.5,
  dpi = 320)


future::plan(NULL)

return(risk_5yr_long_df)
} # end function


# ----------------------------------------------------------------
# ANNUAL RISK 
# ----------------------------------------------------------------
# get annual avg and seasonal averages:
# tot_ann_conc_df <- pred_conc_df %>%
#   group_by(grid_id_10km, year) %>%
#   dplyr::summarise(tot_annual_PB = sum(PB_pred_conc, na.rm=TRUE),
#                    tot_annual_NI = sum(NI_pred_conc, na.rm=TRUE),
#                    tot_annual_AS = sum(AS_pred_conc, na.rm=TRUE)) %>%
#   ungroup() %>% 
#   left_join(grid_in_states_regions, by = 'grid_id_10km')
# 
# # RESIDENTIAL SCENARIO: 26 years
# gridded_ann_cancer_risk_df <- tot_ann_conc_df %>% 
#   mutate(year = as.character(year)) %>% 
#   left_join(all_pop_reweighted_df, by = c('year', 'grid_id_10km')) %>% 
#   mutate(grid_AS_lifetime_cancer_risk = AS_IUR * tot_annual_AS * 24 * 1 * 26 * (1/(24*70*365)) * pop_adj,
#          grid_NI_lifetime_cancer_risk = NI_IUR * tot_annual_NI * 24 * 1 * 26 *(1/(24*70*365)) * pop_adj,
#          grid_PB_lifetime_cancer_risk = PB_IUR * tot_annual_PB * 24 * 1 * 26 * (1/(24*70*365)) * pop_adj) 
# 
# # regional
# regional_cancer_risk <- gridded_ann_cancer_risk_df %>% 
#   group_by(year, region) %>% 
#   dplyr::summarise(AS_canc_burden = sum(grid_AS_lifetime_cancer_risk, na.rm = TRUE),
#                    NI_canc_burden = sum(grid_NI_lifetime_cancer_risk, na.rm = TRUE),
#                    PB_canc_burden = sum(grid_PB_lifetime_cancer_risk, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   mutate(mon_yr = as.Date(paste0(year,"-01-01"))) #
# 
# # pivot longer for easier plotting
# risk_long_df <- regional_cancer_risk %>% 
#   pivot_longer(cols = c('AS_canc_burden', 'NI_canc_burden', 'PB_canc_burden'), 
#                names_to = 'species', values_to = 'canc_burden') %>% 
#   mutate(species = str_remove(species, '_canc_burden')) %>% 
#   dplyr::select(species, year, region,mon_yr, canc_burden) 
# 
# # set plotting date breaks
# datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "12 month")     
# 
# 
# time_series_ann_burden <- ggplot(risk_long_df,# %>% 
#                                  #filter(year %in% c(2006:2010, 2016:2020)),
#                                  aes(x = mon_yr, 
#                                      y = canc_burden, color = region, group = region)) +
#   geom_line(linewidth = .5 ) +
#   geom_point(size = 3) +
#   # geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = species), alpha = 0.3) +
#   # scale_color_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) +
#   # scale_fill_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) +
#   scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#   facet_wrap(~species, ncol = 1, scales = 'free_y') +
#   labs(title = 'Population-level excess cancer burden from annual exposure',
#        x = "Year",
#        y = "Excess cancer burden") +
#   theme_minimal() + 
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     axis.line = element_line(color = "black"),
#     legend.title = element_text(),
#     plot.title = element_text(face = "bold"),
#     legend.key.width = unit(.4, "cm"),  # Set the legend key width
#     legend.key.height = unit(.4, "cm")  # Set the legend key height
#   ) #+
# #guides(color = 'none')
# time_series_ann_burden
# 
# ggsave(
#   filename = paste0('Fig5B_excess_cancer_cases_time_series_annual_PB+AS+NI.pdf'),
#   plot = time_series_ann_burden,
#   path = file.path(results_fp, 'figures/Fig5'),
#   scale = 1,
#   width = 6,
#   height = 8,
#   dpi = 320)
# ggsave(
#   filename = paste0('Fig5B_excess_cancer_cases_time_series_annual_PB+AS+NI.png'),
#   plot = time_series_ann_burden,
#   path = file.path(results_fp, 'figures/Fig5'),
#   scale = 1,
#   width = 6,
#   height = 8,
#   dpi = 320)



# # EXPOSURE ESTIMATES -------------------------------------------------------
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
# EF_cancer =  (24 * 350 * 70) / (24 * 365 * 70)
# gridded_ann_risk_df <- tot_ann_conc_df %>% 
#   mutate(year = as.character(year)) %>% 
#   left_join(all_pop_reweighted_df, by = c('year', 'grid_id_10km')) %>% 
#   mutate(m_AS_lifetime_cancer_risk = tot_annual_AS * AS_IUR * 1/(70*365) * pop_adj,
#          m_NI_lifetime_cancer_risk = tot_annual_NI * NI_IUR * 1/(70*365) * pop_adj,
#          m_PB_lifetime_cancer_risk = tot_annual_PB * PB_IUR * 1/(70*365) * pop_adj) %>% 
#   group_by(year, region) %>% 
#   dplyr::summarise(m_AS_canc_burden = sum(m_AS_lifetime_cancer_risk, na.rm = TRUE),
#                    m_NI_canc_burden = sum(m_NI_lifetime_cancer_risk, na.rm = TRUE),
#                    m_PB_canc_burden = sum(m_PB_lifetime_cancer_risk, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   mutate(mon_yr = as.Date(paste0(year,"-01-01"))) #
# 
# 
# # # calculate the annual risk using the formula
# # gridded_ann_risk_df <- avg_ann_conc_df %>% 
# #   mutate(year = as.character(year)) %>% 
# #   left_join(all_pop_df, by = c('year', 'grid_id_10km')) %>% 
# #   mutate(m_AS_lifetime_cancer_risk = avg_annual_AS * AS_IUR * EF_cancer * pop_ct,
# #          m_NI_lifetime_cancer_risk = avg_annual_NI * NI_IUR * EF_cancer * pop_ct,
# #          m_PB_lifetime_cancer_risk = avg_annual_PB * PB_IUR * EF_cancer * pop_ct,
# #          # upper and lower quantiles
# #          pct95_AS_lifetime_cancer_risk = pct95_annual_AS * AS_IUR * EF_cancer * pop_ct,
# #          pct95_NI_lifetime_cancer_risk = pct95_annual_NI * NI_IUR * EF_cancer * pop_ct,
# #          pct95_PB_lifetime_cancer_risk = pct95_annual_PB * PB_IUR * EF_cancer * pop_ct,
# #          # lower
# #          pct5_AS_lifetime_cancer_risk = pct5_annual_AS * AS_IUR * EF_cancer * pop_ct,
# #          pct5_NI_lifetime_cancer_risk = pct5_annual_NI * NI_IUR * EF_cancer * pop_ct,
# #          pct5_PB_lifetime_cancer_risk = pct5_annual_PB * PB_IUR * EF_cancer * pop_ct) %>% 
# #   group_by(year) %>% 
# #   dplyr::summarise(m_AS_canc_burden = sum(m_AS_lifetime_cancer_risk, na.rm = TRUE),
# #                    m_NI_canc_burden = sum(m_NI_lifetime_cancer_risk, na.rm = TRUE),
# #                    m_PB_canc_burden = sum(m_PB_lifetime_cancer_risk, na.rm = TRUE),
# #                    pct95_AS_canc_burden = sum(pct95_AS_lifetime_cancer_risk, na.rm = TRUE),
# #                    pct95_NI_canc_burden = sum(pct95_NI_lifetime_cancer_risk, na.rm = TRUE),
# #                    pct95_PB_canc_burden = sum(pct95_PB_lifetime_cancer_risk, na.rm = TRUE),
# #                    pct5_AS_canc_burden = sum(pct5_AS_lifetime_cancer_risk, na.rm = TRUE),
# #                    pct5_NI_canc_burden = sum(pct5_NI_lifetime_cancer_risk, na.rm = TRUE),
# #                    pct5_PB_canc_burden = sum(pct5_PB_lifetime_cancer_risk, na.rm = TRUE)) %>% 
# #   ungroup() %>% 
# #   mutate(mon_yr = as.Date(paste0(year,"-01-01"))) #
# 
# ann_risk_long_df <- gridded_ann_risk_df %>% 
#   pivot_longer(cols = c('m_AS_canc_burden', 'm_NI_canc_burden', 'm_PB_canc_burden'), 
#                names_to = 'species', values_to = 'canc_burden') %>% 
#   mutate(species = str_remove(species, '_canc_burden'),
#          species = str_remove(species, 'm_')) %>% 
#   dplyr::select(species, year, region,mon_yr, canc_burden) 
#   # add in UCL
#   # left_join(gridded_ann_risk_df %>% 
#   #             pivot_longer(cols = c('pct95_AS_canc_burden', 'pct95_NI_canc_burden', 'pct95_PB_canc_burden'), 
#   #                          names_to = 'species', values_to = 'pct95_canc_burden') %>% 
#   #             mutate(species = str_remove(species, 'pct95_'), 
#   #                      species = str_remove(species, '_canc_burden')) %>% 
#   #             dplyr::select(species, year, mon_yr, UCL = 'pct95_canc_burden')) %>% 
#   
# 
# datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "168 month")     
# 
# # ANNUAL RISK
# time_series_ann_burden <- ggplot(ann_risk_long_df %>% 
#                                    filter(year %in% c(2006, 2020)),
#                                  aes(x = mon_yr, 
#                                      y = canc_burden, color = region, group = region)) +
#   geom_line(linewidth = .5 ) +
#   geom_point(size = 3) +
#   # geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = species), alpha = 0.3) +
#   # scale_color_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) +
#   # scale_fill_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) +
#   scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#   facet_wrap(~species, ncol = 1, scales = 'free_y') +
#   labs(title = 'Population-level excess cancer burden from annual exposure',
#        x = "Year",
#        y = "Excess cancer burden") +
#   theme_minimal() + 
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     axis.line = element_line(color = "black"),
#     legend.title = element_text(),
#     plot.title = element_text(face = "bold"),
#     legend.key.width = unit(.4, "cm"),  # Set the legend key width
#     legend.key.height = unit(.4, "cm")  # Set the legend key height
#   ) #+
# #guides(color = 'none')
# time_series_ann_burden
# 
# ggsave(
#   filename = paste0('Fig5B_excess_cancer_cases_time_series_annual_PB+AS+NI.pdf'),
#   plot = time_series_ann_burden,
#   path = file.path(results_fp, 'figures/Fig5'),
#   scale = 1,
#   width = 6,
#   height = 8,
#   dpi = 320)
# ggsave(
#   filename = paste0('Fig5B_excess_cancer_cases_time_series_annual_PB+AS+NI.png'),
#   plot = time_series_ann_burden,
#   path = file.path(results_fp, 'figures/Fig5'),
#   scale = 1,
#   width = 6,
#   height = 8,
#   dpi = 320)
# 
# 
# # WILDFIRE SEASON RISK ONLY            
# # alternatively could conduct risk assessment for a wildfire season instead 
# avg_wf_ssn_conc_df <- pred_conc_df %>%
#   mutate(wf_season = ifelse(month %in% c(6:10), "wf season", 'non wf season')) %>%
#   group_by(grid_id_10km, year, wf_season) %>%
#   dplyr::summarise(avg_seasonal_PB = mean(PB_pred_conc, na.rm=TRUE),
#                    avg_seasonal_NI = mean(NI_pred_conc, na.rm=TRUE),
#                    avg_seasonal_AS = mean(AS_pred_conc, na.rm=TRUE),
#                    pct95_seasonal_PB = quantile(PB_pred_conc, .95, na.rm=TRUE),
#                    pct95_seasonal_NI = quantile(NI_pred_conc, .95, na.rm=TRUE),
#                    pct95_seasonal_AS = quantile(AS_pred_conc, .95, na.rm=TRUE),
#                    pct5_seasonal_PB = quantile(PB_pred_conc, .05, na.rm=TRUE),
#                    pct5_seasonal_NI = quantile(NI_pred_conc, .05, na.rm=TRUE),
#                    pct5_seasonal_AS = quantile(AS_pred_conc, .05, na.rm=TRUE)) %>% 
#   ungroup() 
# 
# write.fst(avg_wf_ssn_conc_df, file.path(data_fp, 'clean/avg_wfseason_AS_PB_NI_conc.fst')) 
# 
# 
# # calculate the annual risk using the formula
# gridded_ssnl_risk_df <- avg_wf_ssn_conc_df %>% 
#   mutate(year = as.character(year)) %>% 
#   left_join(all_pop_df, by = c('year', 'grid_id_10km')) %>% 
#   mutate(m_AS_lifetime_cancer_risk = avg_seasonal_AS * AS_IUR * EF_cancer * pop_ct,
#          m_NI_lifetime_cancer_risk = avg_seasonal_NI * NI_IUR * EF_cancer * pop_ct,
#          m_PB_lifetime_cancer_risk = avg_seasonal_PB * PB_IUR * EF_cancer * pop_ct,
#          # upper and lower quantiles
#          pct95_AS_lifetime_cancer_risk = pct95_seasonal_AS * AS_IUR * EF_cancer * pop_ct,
#          pct95_NI_lifetime_cancer_risk = pct95_seasonal_NI * NI_IUR * EF_cancer * pop_ct,
#          pct95_PB_lifetime_cancer_risk = pct95_seasonal_PB * PB_IUR * EF_cancer * pop_ct,
#          # lower
#          pct5_AS_lifetime_cancer_risk = pct5_seasonal_AS * AS_IUR * EF_cancer * pop_ct,
#          pct5_NI_lifetime_cancer_risk = pct5_seasonal_NI * NI_IUR * EF_cancer * pop_ct,
#          pct5_PB_lifetime_cancer_risk = pct5_seasonal_PB * PB_IUR * EF_cancer * pop_ct) %>% 
#   group_by(year, wf_season) %>% 
#   dplyr::summarise(m_AS_canc_burden = sum(m_AS_lifetime_cancer_risk, na.rm = TRUE),
#                    m_NI_canc_burden = sum(m_NI_lifetime_cancer_risk, na.rm = TRUE),
#                    m_PB_canc_burden = sum(m_PB_lifetime_cancer_risk, na.rm = TRUE),
#                    pct95_AS_canc_burden = sum(pct95_AS_lifetime_cancer_risk, na.rm = TRUE),
#                    pct95_NI_canc_burden = sum(pct95_NI_lifetime_cancer_risk, na.rm = TRUE),
#                    pct95_PB_canc_burden = sum(pct95_PB_lifetime_cancer_risk, na.rm = TRUE),
#                    pct5_AS_canc_burden = sum(pct5_AS_lifetime_cancer_risk, na.rm = TRUE),
#                    pct5_NI_canc_burden = sum(pct5_NI_lifetime_cancer_risk, na.rm = TRUE),
#                    pct5_PB_canc_burden = sum(pct5_PB_lifetime_cancer_risk, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   mutate(mon_yr = as.Date(paste0(year,"-01-01"))) #
# 
# ssnl_risk_long_df <- gridded_ssnl_risk_df %>% 
#   pivot_longer(cols = c('m_AS_canc_burden', 'm_NI_canc_burden', 'm_PB_canc_burden'), 
#                names_to = 'species', values_to = 'canc_burden') %>% 
#   mutate(species = str_remove(species, '_canc_burden'),
#          species = str_remove(species, 'm_')) %>% 
#   dplyr::select(species, year, wf_season, mon_yr, canc_burden) %>% 
#   # add in UCL
#   left_join(gridded_ssnl_risk_df %>% 
#               pivot_longer(cols = c('pct95_AS_canc_burden', 'pct95_NI_canc_burden', 'pct95_PB_canc_burden'), 
#                            names_to = 'species', values_to = 'pct95_canc_burden') %>% 
#               mutate(species = str_remove(species, 'pct95_'), 
#                      species = str_remove(species, '_canc_burden')) %>% 
#               dplyr::select(species, year, wf_season,  mon_yr, UCL = 'pct95_canc_burden')) %>% 
#   left_join(gridded_ssnl_risk_df %>% 
#               pivot_longer(cols = c('pct5_AS_canc_burden', 'pct5_NI_canc_burden', 'pct5_PB_canc_burden'), 
#                            names_to = 'species', values_to = 'pct5_canc_burden') %>% 
#               mutate(species = str_remove(species, 'pct5_'),
#                      species = str_remove(species, '_canc_burden')) %>% 
#               dplyr::select(species, year,wf_season, mon_yr, LCL = 'pct5_canc_burden'))
# 
# 
#     # PLOT THE TIME SERIES OF EXCESS CANCER CASES BY WILDFIRE SEASON
#     ssnl_time_series_cases <- ggplot(ssnl_risk_long_df %>% 
#                                        filter(wf_season == 'wf season'),
#                                 aes(x = mon_yr,
#                                     y = canc_burden, color = species)) +
#       geom_line(linewidth = 2) +
#       geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = species), alpha = 0.3) +
#       scale_color_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) +
#       scale_fill_manual(values = c('AS' = '#E2D200', 'PB' = '#46ACC8', 'NI' = '#E58601')) +
#       scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#       facet_wrap(~species, ncol = 1, scales = 'free_y') +
#       scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#       labs(title = 'Population-level excess cancer burden from wildfire season exposure',
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
#     ssnl_time_series_cases
#     
#     ggsave(
#       filename = paste0('SIFig_excess_cancer_cases_time_series_wildfire_season_PB+AS+NI.pdf'),
#       plot = ssnl_time_series_cases,
#       path = file.path(results_fp, 'figures/SI Figs'),
#       scale = 1,
#       width = 7,
#       height = 8,
#       dpi = 320)
#     ggsave(
#       filename = paste0('SIFIG_excess_cancer_cases_time_series_wildfire_season_PB+AS+NI.png'),
#       plot = ssnl_time_series_cases,
#       path = file.path(results_fp, 'figures/SI Figs'),
#       scale = 1,
#       width = 7,
#       height = 8,
#       dpi = 320)
# 
#     
#     #  ---------------------------------------------
#     # CANCER SCREENING LEVEL ASSESMENT 
#     # ---------------------------------------------
#     # get number of days above screening levels
#     AS_SL <- IUR_vals %>%
#       filter(ParamCode == 'AS') %>%
#       distinct(SL_residential_cancer) %>%
#       filter(!is.na(SL_residential_cancer)) %>% 
#       as.numeric()
#     
#     PB_SL <- IUR_vals %>%
#       filter(ParamCode == 'PB') %>%
#       distinct(SL_residential_cancer) %>%
#       filter(!is.na(SL_residential_cancer)) %>% 
#       as.numeric()
#     
#     NI_SL <- IUR_vals %>%
#       filter(ParamCode == 'NI') %>%
#       distinct(SL_residential_cancer) %>%
#       filter(!is.na(SL_residential_cancer)) %>% 
#       as.numeric()
#     
#     # calculate if a grid cell day is above the SL
#     grid_days_above_SL_df <- pred_conc_df %>%
#       mutate(AS_above = ifelse(AS_pred_conc >= AS_SL, TRUE, FALSE),
#              PB_above = ifelse(PB_pred_conc >= PB_SL, TRUE, FALSE),
#              NI_above = ifelse(NI_pred_conc >= NI_SL, TRUE, FALSE)) %>% 
#       group_by(grid_id_10km, year) %>% 
#       dplyr::summarise(AS_grid_ann_above = sum(AS_above, na.rm = TRUE),
#                        PB_grid_ann_above = sum(PB_above, na.rm = TRUE),
#                        NI_grid_ann_above = sum(NI_above, na.rm = TRUE)) %>% 
#       ungroup()
#     
#     # get annual summary
#     ann_sum_grid_days <- grid_days_above_SL_df %>% 
#       dplyr::select(-PB_grid_ann_above, -NI_grid_ann_above) %>% 
#       mutate(year = as.character(year)) %>% 
#       left_join(all_pop_df, by = c('year', 'grid_id_10km')) %>% 
#       mutate(numerator = pop_ct*AS_grid_ann_above,
#              denominator =pop_ct*365) %>% 
#       group_by(year) %>% 
#       dplyr::summarise(numerator_sum = sum(numerator, na.rm = TRUE),
#                        denominator_smm = sum(denominator, na.rm = TRUE)) %>% 
#       ungroup() %>% 
#       mutate(mon_yr = as.Date(paste0(year,"-01-01"))) %>% 
#       mutate(pct_pop_days_above = 100*(numerator_sum/denominator_smm))
#     
#     
#     
#     time_series_days_above <- ggplot(ann_sum_grid_days,
#                                      aes(x = mon_yr, 
#                                          y = pct_pop_days_above, color = species)) +
#       geom_line(linewidth = 2) +
#       scale_color_manual(values = c('AS' = '#E2D200')) +
#       scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#       labs(title = ' cancer risk screening level',
#            x = "Year",
#            y = "% of population-days above threshold") +
#       theme_minimal() + 
#       theme(
#         panel.grid.major = element_blank(),  # Remove major grid lines
#         panel.grid.minor = element_blank(),  # Remove minor grid lines
#         axis.line = element_line(color = "black"),
#         legend.title = element_text(),
#         plot.title = element_text(face = "bold"),
#         legend.key.width = unit(.4, "cm"),  # Set the legend key width
#         legend.key.height = unit(.4, "cm")  # Set the legend key height
#       ) #+
#     #guides(color = 'none')
#     time_series_days_above
#     
#     ggsave(
#       filename = paste0('SIFig_AS_days_above_risk_screening_level.pdf'),
#       plot = time_series_days_above,
#       path = file.path(results_fp, 'figures/Fig5'),
#       scale = 1,
#       width = 7,
#       height = 3,
#       dpi = 320)
#     



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
  




