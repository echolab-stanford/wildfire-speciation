# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 11, 2024
# Description: create population weighted exposure

# 
# loadd(c(region_pal, us_region_map), cache = drake_cache)
# csn_param_fp = file.path(data_fp, 'raw/CSN_parameters.csv')
# grid_fp = file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')
# tox_fp = file.path(data_fp, 'raw/RSL Table/oehha_ca_IURs.csv')
# us_pop_fp = file.path(data_fp, 'raw/united-states-population-2024-02-08.csv')

calculate_risk_from_pop_exposure <- function(comb_tox_preds_df, 
                                             csn_param_fp, 
                                             grid_fp, 
                                             tox_fp, 
                                             us_pop_fp, 
                                             us_region_map, 
                                             region_pal) {

# Read in grid + transform
grid_10km <- st_read(grid_fp) %>% 
  dplyr::select(grid_id_10km = 'ID')

# see which grid cells are in which states
grid_in_states_regions <- st_join(us_region_map, grid_10km) %>% 
  st_drop_geometry() 

# -----------------------------------
# read in population data from GEE
# get list of files to read in:
file_list <- list.files(file.path(gee_data_fp), pattern = '^pop', full.names = TRUE)
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

# check pop counts again:
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
pred_conc_df <- comb_tox_preds_df %>% 
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

