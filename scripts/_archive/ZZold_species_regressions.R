# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 16, 2023
# Description: run regressions for each species

# loadd(clean_pm_spec_df, cache = drake_cache)


# FUNCTION ----------- ----------- ----------- ----------- -----------
# run_regressions <- function(clean_pm_spec_df) {

reg_df <- clean_pm_spec_df %>% 
  # drop vars that are not needed for regs
  dplyr::select(-c(units, Elevation, Latitude, Longitude, epsg, km2fire, smoke_day)) %>% 
  # classify date vars as factors for regression in later steps
  mutate(SiteCode = as.factor(SiteCode),
         year = as.factor(year),
         month = as.factor(month),
         moy = as.factor(paste0(month,'-', year)),
         monitor_month = paste0(SiteCode, "_", month)) %>% 
  rename(nonsmokePM = 'non_smokePM_cont',
         smokePM = 'smokePM_pred',
         totPM2.5 = 'PM2.5') %>% 
  pivot_longer(cols = totPM2.5:ZR,
               names_to = 'species_name',
               values_to = 'conc_val') %>% 
  group_by(Dataset, SiteCode, Date, State, season, region, year, 
           month, doy, species_name, smoke_day) %>% 
  dplyr::summarise(conc_val = mean(conc_val, na.rm=TRUE)) %>% 
  ungroup()

  # now get the regional concentration average (aka baseline) for each species and region
  # pivot_longer(cols = totPM2.5:ZR,
  #              names_to = 'lhs', 
  #              values_to = 'conc') %>% 
  # group_by(region, lhs) %>% 
  # mutate(region_avg_conc = mean(conc, na.rm = TRUE)) %>% 
  # ungroup() %>% 
  # mutate(region_conc_ratio = conc/region_avg_conc) %>% 
  # pivot_wider(id_cols = c('Dataset', 'SiteCode', 'Date', 'State', 
  #                         'season', 'region', 'year', 
  #                         'month', 'doy', 'moy', 'monitor_month'), 
  #             names_from = 'lhs', 
  #             values_from = c('conc', 'region_avg_conc', 'region_conc_ratio')) %>% 
  distinct() 


#######################################################################################
# RUN DIFFERENT MODEL VERSIONS 
#######################################################################################
# -------------------------------------------------------------------------------
# 1) Monitor-day Species conc = smokePMconc + nonsmokePMconc + monitor-monthFE + yearFE 
# -------------------------------------------------------------------------------
predPM_reg = feols(c(conc_MF, conc_RC_PM2, conc_NH4, conc_tot_metals, 
                         conc_ammNO3, conc_ammSO4, conc_AL,conc_AS, conc_BR,conc_CA,
                         conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR,
                         conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
                         conc_NA,conc_NI,conc_NO3,conc_N2,
                         conc_P,conc_PB,conc_RB,conc_S,conc_SE,
                         conc_SI,conc_SOIL,conc_SO4,conc_SR,
                         conc_V,conc_ZN,
                     conc_ZR) ~ smokePM_pred + non_smokePM_cont | monitor_month + year, reg_df)

# Save coefficients to table
predPM_coeffs <- coeftable(predPM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(species = str_remove(lhs, 'conc_')) %>% 
  mutate(pm_type = case_when(
    coefficient == 'smokePM_pred' ~ 'smoke PM',
    coefficient == 'non_smokePM_cont' ~ 'non-smoke PM')) %>% 
  dplyr::select(-lhs, -id)
    
# plot coefficients for the speciation at the avg monitor ---------------------------
pred_reg_plot <- ggplot(predPM_coeffs, 
                        aes(x = reorder(species, Estimate), 
                            y = Estimate, color=pm_type,
                            ymin = pmax(0, Estimate - se), 
                            ymax = Estimate + se)) +
  geom_point(size=3, aes(shape = sig), alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(width = 0.1), stat = "identity") +
  # LEAVE THE BELOW COMMENTED IF PLOTTING IN LEVELS, UNCOMMENT IF WANT PLOT TO BE IN LOG SCALE
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(scale = 100)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'PM Type',
       # title = "Speciation of smoke vs nonsmoke PM2.5 (Marissa data), levels") +
       title = "Speciation of smoke vs nonsmoke PM2.5 (Marissa data), logs") +
       #subtitle = 'Monitor-day species conc = smokePMconc + nonsmokePMconc + monitor-monthFE + yearFE') +
  coord_flip() +
  theme_light() 
pred_reg_plot


# -------------------------------------------------------------------------------
# 2) Monitor-day Species conc = smokePMconc + nonsmokePMconc + monitor-monthFE + yearFE 
# -------------------------------------------------------------------------------
stationPM_reg = feols(c(conc_MF, conc_RC_PM2, conc_NH4, conc_tot_metals, 
                         conc_ammNO3, conc_ammSO4, conc_AL,conc_AS,conc_BR,conc_CA,
                         conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR,
                         conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
                         conc_NA,conc_NI,conc_NO3,conc_N2,
                         conc_P,conc_PB,conc_RB,conc_S,conc_SE,
                         conc_SI,conc_SOIL,conc_SO4,conc_SR,
                         conc_V,conc_ZN,conc_ZR) ~ station_calc_smokePM + station_non_smokePM | monitor_month + year, reg_df)

# Save coefficients to table
station_PM_coeffs <- coeftable(stationPM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(species = str_remove(lhs, 'conc_')) %>% 
  mutate(pm_type = case_when(
    coefficient == 'station_calc_smokePM' ~ 'station smoke PM',
    coefficient == 'station_non_smokePM' ~ 'station non-smoke PM')) %>% 
  dplyr::select(-lhs, -id)

# plot coefficients for the speciation at the avg monitor ---------------------------
station_PM_plot <- ggplot(station_PM_coeffs, 
                        aes(x = reorder(species, Estimate), 
                            y = Estimate, color=pm_type,
                            ymin = pmax(0, Estimate - se), 
                            ymax = Estimate + se)) +
  geom_point(size=3, aes(shape = sig), alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(width = 0.1), stat = "identity") +
  # LEAVE THE BELOW COMMENTED IF PLOTTING IN LEVELS, UNCOMMENT IF WANT PLOT TO BE IN LOG SCALE
  # scale_y_continuous(trans = log10_trans(),
  #                    breaks = trans_breaks("log10", function(x) 10^x),
  #                    labels = scales::percent_format(scale = 100)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Significance based on p < .05',
       title = "Speciation of smoke vs nonsmoke PM2.5 (station data), levels") +
       #title = "Speciation of smoke vs nonsmoke PM2.5 (station data), logs") +
  coord_flip() +
  theme_light() 
station_PM_plot



#######################################################################################
# RUN THE MAIN MODEL -----------------------------------------------------------------
#######################################################################################
# set up main regression for each species:
# Monitor-day Species conc = smokePMconc + monitorFE + yearFE + monthFE
# run the main model
main_model_reg = feols(c(conc_MF, conc_RC_PM2, conc_NH4, conc_tot_metals, 
                         conc_ammNO3, conc_ammSO4, conc_AL,conc_AS,conc_BR,conc_CA,
                         conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR,
                         conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
                         conc_MO,conc_NA,conc_NI,conc_NO3,conc_N2,
                         conc_P,conc_PB,conc_RB,conc_S,conc_SE,
                         conc_SI,conc_SOIL,conc_SO4,conc_SR,conc_V,conc_ZN,conc_ZR) ~ smokePM_pred | SiteCode + year + month, reg_df)

# Save coefficients to table
main_model_reg <- summary(main_model_reg, se = "hetero")
main_coeffs <- coefficients(main_model_reg) %>% 
  rename(coef = 'smokePM_pred') %>% 
  mutate(lhs = str_remove(lhs, 'conc_'))

# get standard errors:
main_std_errors <- fixest::se(main_model_reg) %>% 
  rename(se = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

# get pvalues
main_pval <- pvalue(main_model_reg) %>%
  rename(pval = 'smokePM_pred') %>%
  dplyr::select(-id, -lhs) %>% 
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))

main_reg_res_df <- cbind(main_coeffs, main_std_errors, main_pval) %>% # main_se, main_conf, main_pval) %>% 
  mutate(model= 'Main Model w/ Monitor, Year, Month FE') 

# plot coefficients for the speciation at the avg monitor ---------------------------
main_reg_plot <- ggplot(main_reg_res_df, 
                        aes(x = reorder(lhs, coef), 
                            y = coef, color=sig,
                            ymin = pmax(0, coef - se), 
                            ymax = coef + se)) +
  geom_point(size=3, shape = 15, alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(width = 0.1), stat = "identity") +
  # LEAVE THE BELOW COMMENTED IF PLOTTING IN LEVELS, UNCOMMENT IF WANT PLOT TO BE IN LOG SCALE
  # scale_y_continuous(trans = log10_trans(),
  #                    breaks = trans_breaks("log10", function(x) 10^x),
  #                    labels = scales::percent_format(scale = 100)) + 
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Significance based on p < .05',
       title = "Speciation of smoke PM2.5 (contiuous measure of smoke), levels",
       #title = "Speciation of smoke PM2.5 (contiuous measure of smoke), logs",
       subtitle = 'Monitor-day species conc = smokePM2.5 + monitor FE + year FE + month FE') +
  coord_flip() +
  theme_light() 
main_reg_plot

rm(main_coeffs, main_pval, main_std_errors, main_reg_plot) # save memory
#######################################################################################
# (1) RUN MAIN MODEL BUT DIVIDE EACH SPECIES CONCENTRATION BY REGIONAL MEANS
# (2) RUN MAIN MODEL, THEN DIVIDE COEFFICIENTS BY REGIONAL MEANS
#######################################################################################
# (1) set up df where each concentration has been divided
regional1_main_model_reg = feols(c(region_conc_ratio_MF, region_conc_ratio_RC_PM2, 
                                   region_conc_ratio_NH4, region_conc_ratio_tot_metals, 
                                   region_conc_ratio_ammNO3, region_conc_ratio_ammSO4, 
                                   region_conc_ratio_AL,region_conc_ratio_AS,
                                   region_conc_ratio_BR,region_conc_ratio_CA,
                                   region_conc_ratio_EC,region_conc_ratio_OC,region_conc_ratio_CHL,
                                   region_conc_ratio_CL,region_conc_ratio_CR,
                                   region_conc_ratio_CU,region_conc_ratio_FE,region_conc_ratio_K,
                                   region_conc_ratio_MG, region_conc_ratio_MN,
                                   region_conc_ratio_MO,region_conc_ratio_NA,region_conc_ratio_NI,
                                   region_conc_ratio_NO3,region_conc_ratio_N2, region_conc_ratio_P,
                                   region_conc_ratio_PB,region_conc_ratio_RB,region_conc_ratio_S,
                                   region_conc_ratio_SE, region_conc_ratio_SI,
                                   region_conc_ratio_SOIL,region_conc_ratio_SO4,region_conc_ratio_SR,
                                   region_conc_ratio_V,region_conc_ratio_ZN,
                                   region_conc_ratio_ZR) ~ smokePM_pred | SiteCode + year + month, reg_df, fsplit = ~region)

# Save coefficients to table
regional1_main_model_reg <- summary(regional1_main_model_reg, se = "hetero")
regional1_main_coeffs <- coefficients(regional1_main_model_reg) %>% 
  rename(coef = 'smokePM_pred') 

# get standard errors:
regional1_main_std_errors <- fixest::se(regional1_main_model_reg) %>% 
  rename(se = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs, -sample.var, -sample)

# get pvalues
regional1_main_pval <- pvalue(regional1_main_model_reg) %>%
  rename(pval = 'smokePM_pred') %>%
  dplyr::select(-id, -lhs, -sample.var, -sample) %>% 
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))

regional1_main_reg_res_df <- cbind(regional1_main_coeffs, regional1_main_std_errors, regional1_main_pval) %>% 
  dplyr::select(-id) %>% 
  mutate(lhs = str_remove(lhs, 'region_conc_ratio_')) %>% 
  mutate(model= 'Main Model w/ FE for year, month, site | lhs: ratio - conc/reg avg conc') 

# (2) run regression with normal conc and then divide by regional means
regional2_main_model_reg = feols(c(conc_MF, conc_RC_PM2, conc_NH4, conc_tot_metals, 
                                   conc_ammNO3, conc_ammSO4, conc_AL,conc_AS,conc_BR,conc_CA,
                                   conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR,
                                   conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
                                   conc_MO,conc_NA,conc_NI,conc_NO3,conc_N2,
                                   conc_P,conc_PB,conc_RB,conc_S,conc_SE,
                                   conc_SI,conc_SOIL,conc_SO4,conc_SR,conc_V,
                                   conc_ZN,conc_ZR) ~ smokePM_pred | SiteCode + year + month, reg_df, fsplit = ~region)

# Save coefficients to table
regional2_main_model_reg <- summary(regional2_main_model_reg, se = "hetero")
regional2_main_coeffs <- coefficients(regional2_main_model_reg) %>% 
  rename(coef = 'smokePM_pred') 

# get standard errors:
regional2_main_std_errors <- fixest::se(regional2_main_model_reg) %>% 
  rename(se = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs, -sample.var, -sample)

# get pvalues
regional2_main_pval <- pvalue(regional2_main_model_reg) %>%
  rename(pval = 'smokePM_pred') %>%
  dplyr::select(-id, -lhs, -sample.var, -sample) %>% 
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))

# GET REGIONAL AVG CONC TO DIVIDE ESTIMATES BY
reg_avgs <- reg_df %>% 
  dplyr::select(sample = 'region', region_avg_conc_MF, region_avg_conc_RC_PM2, region_avg_conc_NH4, region_avg_conc_tot_metals, 
                region_avg_conc_ammNO3, region_avg_conc_ammSO4, region_avg_conc_AL,region_avg_conc_AS,
                region_avg_conc_BR,region_avg_conc_CA, region_avg_conc_EC,region_avg_conc_OC,region_avg_conc_CHL,
                region_avg_conc_CL,region_avg_conc_CR, region_avg_conc_CU,region_avg_conc_FE,
                region_avg_conc_K,region_avg_conc_MG, region_avg_conc_MN, region_avg_conc_MO,
                region_avg_conc_NA,region_avg_conc_NI,region_avg_conc_NO3,region_avg_conc_N2,
                region_avg_conc_P,region_avg_conc_PB,region_avg_conc_RB,region_avg_conc_S,region_avg_conc_SE,
                region_avg_conc_SI,region_avg_conc_SOIL,region_avg_conc_SO4,region_avg_conc_SR,region_avg_conc_V,
                region_avg_conc_ZN,region_avg_conc_ZR) %>% 
  distinct() %>% 
  pivot_longer(cols = c('region_avg_conc_MF', 'region_avg_conc_RC_PM2', 'region_avg_conc_NH4', 'region_avg_conc_tot_metals', 
                        'region_avg_conc_ammNO3', 'region_avg_conc_ammSO4', 'region_avg_conc_AL','region_avg_conc_AS',
                        'region_avg_conc_BR','region_avg_conc_CA', 'region_avg_conc_EC','region_avg_conc_OC','region_avg_conc_CHL',
                        'region_avg_conc_CL','region_avg_conc_CR', 'region_avg_conc_CU','region_avg_conc_FE',
                        'region_avg_conc_K','region_avg_conc_MG', 'region_avg_conc_MN', 'region_avg_conc_MO',
                        'region_avg_conc_NA','region_avg_conc_NI','region_avg_conc_NO3','region_avg_conc_N2',
                        'region_avg_conc_P','region_avg_conc_PB','region_avg_conc_RB','region_avg_conc_S','region_avg_conc_SE',
                        'region_avg_conc_SI','region_avg_conc_SOIL','region_avg_conc_SO4','region_avg_conc_SR','region_avg_conc_V',
                        'region_avg_conc_ZN','region_avg_conc_ZR'),
               names_to = 'lhs', 
               values_to = 'avg_reg_conc') %>% 
  mutate(lhs = str_remove(lhs, 'region_avg_conc_')) %>% 
  distinct() 

# now create a coefficient that has been divided by regional averages
regional2_main_reg_res_df <- cbind(regional2_main_coeffs, regional2_main_std_errors, regional2_main_pval) %>% 
  dplyr::select(-id) %>% 
  mutate(lhs = str_remove(lhs, 'conc_')) %>% 
  mutate(model= 'Main Model w/ FE for year, month, site | lhs: coeff/reg avg') %>% 
  left_join(reg_avgs, by = c('sample', 'lhs')) %>% 
  mutate(coeff_div_reg_avg = coef/avg_reg_conc)

# BIND BOTH VERSION OF THE REGIONAL REGRESSION TOGETHER + PLOT
both_regional_reg_df <- bind_rows(regional1_main_reg_res_df, regional2_main_reg_res_df)

# plot coefficients for each version of handling the regional baseline ---------------------------
regional_reg_plot <- ggplot(both_regional_reg_df, 
                        aes(x = reorder(lhs, coef), 
                            y = coef,
                            ymin = coef - se, 
                            ymax = coef + se)) +
  geom_point(size=3, alpha = 0.6, aes(shape = sig, color = sample), stat = "identity") +
  #geom_errorbar(aes(width = 0.1), stat = "identity") +
  # LEAVE THE BELOW COMMENTED IF PLOTTING IN LEVELS, UNCOMMENT IF WANT PLOT TO BE IN LOG SCALE
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(scale = 100)) +
  scale_color_manual(values=c('black', "coral","steelblue", "forestgreen", 'plum', 'goldenrod', 'mediumpurple')) +
  scale_shape_manual(name="Significance", values=c(18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  facet_wrap(~model) +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Significance based on p < .05',
       # title = "Speciation of smoke PM2.5 (contiuous measure of smoke), levels",
       title = "Speciation of smoke PM2.5 (contiuous measure of smoke), logs",
       subtitle = 'Comparing regional coefficients, one model per region') +
  coord_flip() +
  theme_light() 
regional_reg_plot

# drop extra objects to save memory
rm(regional1_main_model_reg, regional1_main_pval, regional1_main_coeffs, regional1_main_std_errors,
   regional2_main_model_reg, regional2_main_pval, regional2_main_coeffs, regional2_main_std_errors)

#######################################################################################
# SMOKE PM VS TOTAL PM2.5 ANALYSIS
#######################################################################################
# calculate non-wildfire PM at the station level by just subtracting wildfire PM from station-measured total PM on that day
pm_types_df <- reg_df %>% 
  filter(!is.na(non_smokePM))

pm_types_reg = feols(c(conc_MF, conc_RC_PM2, conc_NH4, conc_tot_metals, 
                       conc_ammNO3, conc_ammSO4, conc_AL,conc_AS,conc_BR,conc_CA,
                       conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR,
                       conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
                       conc_NA,conc_NI,conc_NO3,conc_N2,
                       conc_P,conc_PB,conc_RB,conc_S,conc_SE,
                       conc_SI,conc_SOIL,conc_SO4,conc_SR,conc_V,
                       conc_ZN,conc_ZR) ~ non_smokePM | SiteCode + year + month, pm_types_df)

pm_types_reg <- summary(pm_types_reg, se = "hetero")
# get coefficients
pm_coeffs <- coefficients(pm_types_reg) %>% 
  dplyr::select(-id) %>% 
  rename(coef = 'non_smokePM') %>% 
  mutate(lhs = str_remove(lhs, 'conc_'))

# get standard errors:
pm_std_errors <- fixest::se(pm_types_reg) %>% 
  rename(se = 'non_smokePM') %>% 
  dplyr::select(-id, -lhs)

# get pvalues
pm_pval <- pvalue(pm_types_reg) %>%
  rename(pval = 'non_smokePM') %>%
  dplyr::select(-id, -lhs) %>% 
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))

# bind smoke, nonsmoke, and total pm together
pm_reg_res_df <- cbind(pm_coeffs, pm_std_errors, pm_pval) %>% # main_se, main_conf, main_pval) %>% 
  mutate(model= '*NONSMOKE* | Main Model w/ FE for year, month, site')


# total PM -------------------
total_pm_types_reg = feols(c(conc_RC_PM2, conc_NH4, conc_tot_metals, 
                       conc_ammNO3, conc_ammSO4, conc_AL,conc_AS,conc_BR,conc_CA,
                       conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR,
                       conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
                       conc_NA,conc_NI,conc_NO3,conc_N2,
                       conc_P,conc_PB,conc_RB,conc_S,conc_SE,
                       conc_SI,conc_SOIL,conc_SO4,conc_SR,conc_V,
                       conc_ZN,conc_ZR) ~ conc_MF | SiteCode + year + month, pm_types_df)

total_pm_types_reg <- summary(total_pm_types_reg, se = "hetero")
# get coefficients
total_pm_coeffs <- coefficients(total_pm_types_reg) %>% 
  dplyr::select(-id) %>% 
  rename(coef = 'conc_MF') %>% 
  mutate(lhs = str_remove(lhs, 'conc_'))

# get standard errors:
total_pm_std_errors <- fixest::se(total_pm_types_reg) %>% 
  rename(se = 'conc_MF') %>% 
  dplyr::select(-id, -lhs)

# get pvalues
total_pm_pval <- pvalue(total_pm_types_reg) %>%
  rename(pval = 'conc_MF') %>%
  dplyr::select(-id, -lhs) %>% 
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))

# bind together
total_pm_reg_res_df <- cbind(total_pm_coeffs, total_pm_std_errors, total_pm_pval) %>%  
  mutate(model= '*TOTAL PM* | Main Model w/ FE for year, month, site')

##############
all_smoke_nonsmoke <- pm_reg_res_df %>% 
  bind_rows(main_reg_res_df %>% 
              dplyr::select(-id) %>% 
              mutate(model = '*SMOKE* | Main Model w/ FE for year, month, site')) #%>% 
  #bind_rows(total_pm_reg_res_df)


# PLOT A COMPARISON
pm_types_reg_plot <- ggplot(all_smoke_nonsmoke, 
                            aes(x = reorder(lhs, coef), 
                                y = coef,
                                ymin = coef - se, 
                                ymax = coef + se)) +
  geom_point(size=3, alpha = 0.6, aes(shape = sig, color = model), stat = "identity") +
  geom_errorbar(aes(width = 0.1), stat = "identity") +
  scale_color_manual(name="Model",values=c('dodgerblue', "forestgreen", 'plum')) +
  scale_shape_manual(name="Significance", values=c(18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  # LEAVE THE BELOW COMMENTED IF PLOTTING IN LEVELS, UNCOMMENT IF WANT PLOT TO BE IN LOG SCALE
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(scale = 100)) +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Significance based on p < .05',
       title = "Speciation of PM2.5 (logs)",
       # title = "Speciation of PM2.5 (levels)",
       subtitle = 'Comparing smoke PM, nonsmoke PM, and total PM') +
  coord_flip() +
  theme_light() 
pm_types_reg_plot

#######################################################################################
# Sensitivity analysis using different FE:
#######################################################################################
# 1. Monitor + year + month of year ----------------------------------------------------
# moy_FE = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
#                              AL,AS,BR,CA,EC,OC,CHL,CL,CR,
#                              CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
#                              P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + year + moy, reg_df)
# 
# # Save coefficients to table (this is necessary anyways when we make a plot)
# moy_coeffs <- coefficients(moy_FE) %>% 
#   rename(coef = 'smokePM_pred')
# 
# # get standard errors
# moy_se <- fixest::se(moy_FE) %>% 
#   rename(se = 'smokePM_pred') %>% 
#   dplyr::select(-id, -lhs)
# 
# # get pvalues
# moy_pval <- fixest::pvalue(moy_FE) %>% 
#   rename(pval = 'smokePM_pred') %>% 
#   dplyr::select(-id, -lhs)
# 
# moyFE_res_df <- cbind(moy_coeffs, moy_se, moy_pval) %>% 
#   mutate(model= 'Monitor + year + month-of-year FE')
# 
# # 2. Monitor + month-of-sample ----------------------------------------------------
# month_of_sample_FE = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
#                  AL,AS,BR,CA,EC,OC,CHL,CL,CR,
#                  CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
#                  P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + month, reg_df)
# 
# # Save coefficients to table (this is necessary anyways when we make a plot)
# month_of_sample_coeffs <- coefficients(month_of_sample_FE) %>% 
#   rename(coef = 'smokePM_pred')
# 
# # get standard errors
# month_of_sample_se <- fixest::se(month_of_sample_FE) %>% 
#   rename(se = 'smokePM_pred') %>% 
#   dplyr::select(-id, -lhs)
# 
# # get pvalues
# month_of_sample_pval <- fixest::pvalue(month_of_sample_FE) %>% 
#   rename(pval = 'smokePM_pred') %>% 
#   dplyr::select(-id, -lhs)
# 
# month_of_sample_res_df <- cbind(month_of_sample_coeffs, month_of_sample_se, month_of_sample_pval) %>% 
#   mutate(model= 'Monitor + month-of-sample FE')
# 
# # 3. Monitor * month of year + year----------------------------------------------------
# # 4. Drop negatives -----------------------------------------------------------------
# no_negs <- reg_df %>% 
#   mutate_if(is.numeric, ~ifelse(. < 0, NA, .)) 
# 
# no_negs_reg = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
#                              AL,AS,BR,CA,EC,OC,CHL,CL,CR,
#                              CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
#                              P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + year + month, reg_df)
# 
# # Save coefficients to table (this is necessary anyways when we make a plot)
# no_negs_coeffs <- coefficients(no_negs_reg) %>% 
#   rename(coef = 'smokePM_pred')
# 
# # get standard errors
# no_negs_se <- fixest::se(no_negs_reg) %>% 
#   rename(se = 'smokePM_pred') %>% 
#   dplyr::select(-id, -lhs)
# 
# # get pvalues
# no_negs_pval <- fixest::pvalue(no_negs_reg) %>% 
#   rename(pval = 'smokePM_pred') %>% 
#   dplyr::select(-id, -lhs)
# 
# no_negs_res_df <- cbind(no_negs_coeffs, no_negs_se, no_negs_pval) %>% 
#   mutate(model= 'Main, no negatives')
# 
# ##################################################################
# # 5. Change negatives to detection limit values ---------
# ##################################################################
# # bind all results together and plot
# all_results <- bind_rows(main_reg_res_df, moyFE_res_df, 
#                          month_of_sample_res_df, no_negs_res_df) # %>% 
#   # drop the measures of carbon
#   # filter(!lhs %in% c('RC_PM2', 'MF', 'OC'))
# 
# all_results$model <- fct_rev(factor(all_results$model, 
#                                      levels = c("Main w/ Monitor, Year, Month FE",
#                                                 "Monitor + year + month-of-year FE",
#                                                 "Monitor + month-of-sample FE",
#                                                 "Main, no negatives")))
# 
# col <- c('dodgerblue', 'forestgreen', 'salmon', 'plum')
# 
# # EXAMPLE PLOT
# # plot avg monitor
# all_reg_plot <- ggplot(all_results, 
#                         aes(x = lhs, y = coef, color=model)) +
#   geom_point(aes(shape=model),size=3, alpha = 0.6) +
#   scale_color_manual(name="Model",values=c("coral","steelblue", "forestgreen", 'plum')) +
#   scale_shape_manual(name="Model",values=c(15, 17,18, 19)) +
#   geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
#   #geom_linerange(aes(ymin = coef-se, ymax = coef+se)) +
#   scale_y_continuous(
#     trans = log10_trans(),
#     breaks = trans_breaks("log10", function(x) 10^x),
#     labels = trans_format("log10", math_format(10^.x))
#   ) +
#   theme(axis.text.x=element_text(angle=-45)) +
#   ggtitle("Regression Results") +
#   labs(y = 'Coefficient',
#        x = 'Chemical Species',
#        color = 'Model',
#        title = "Sensitivity Analysis") +
#   coord_flip() +
#   theme_light() 
# all_reg_plot


##################################################################
# 5. divide each species by total PM (MF), use smoke v non-smoke day as the RHS 
##################################################################
# - coefficient is relative to non-smoke day
# - regional plots, smoke vs non-smoke day, overall smoke
# - make one plot per region using that map

# species_frac <- clean_pm_spec_df %>% 
#   mutate_at(c('NH4', 'tot_metals', 'ammNO3', 'ammSO4',
#               'AL','AS','BR','CA','EC','OC','CHL','CL','CR',
#               'CU','FE','K','MG', 'MN','MO','NA','NI',
#               'NO3','N2','P','PB','RB','S','SE','SI','SOIL',
#               'SO4','SR','V','ZN','ZR'),  list(frac = ~./MF)) %>% 
#   filter(!is.na(MF))
#   
# 
# # Run a regression for each fraction
# smoke_day_model_reg = feols(c(NH4_frac,tot_metals_frac,
#                               ammNO3_frac, ammSO4_frac,
#                               AL_frac, AS_frac, BR_frac, CA_frac,
#                               EC_frac, OC_frac, CHL_frac, CL_frac, CR_frac,
#                               CU_frac,FE_frac, K_frac, MG_frac, MN_frac,
#                               NA_frac, NI_frac, NO3_frac,N2_frac,
#                               P_frac,PB_frac,RB_frac,S_frac,SE_frac,
#                               SI_frac,SOIL_frac,SO4_frac,SR_frac,
#                               V_frac,ZN_frac,ZR_frac) ~ as.factor(smoke_day) | SiteCode + year + month, 
#                             species_frac, fsplit = ~region)
#                              
# smoke_day_model_reg <- summary(smoke_day_model_reg, se = "hetero")
# region_coeffs <- coefficients(smoke_day_model_reg) %>% 
#   rename(coef = 'as.factor(smoke_day)1') %>% 
#   mutate(model= 'Binary Smoke Day regressed on Species Fraction of Total PM2.5') 
# 
# # # get pvalues
# region_pval <- pvalue(smoke_day_model_reg) %>%
#   rename(pval = 'as.factor(smoke_day)1') %>%
#   dplyr::select(-id, -lhs, -sample, -sample.var) %>% 
#   mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))
# 
# region_model_df <- cbind(region_coeffs, region_pval) %>% 
#   dplyr::select(lhs, coef, pval, sig, region = 'sample')
# 
# # PLOT
# all_region_plot <- ggplot(region_model_df, 
#                           aes(x = reorder(lhs, coef), y = coef)) +
#   geom_point(aes(shape=sig, color = region), size=3, alpha = 0.6) +
#   scale_color_manual(name="Region",
#                      values=c("coral","steelblue", "forestgreen", 'plum', 'mediumpurple', 'goldenrod2', 'darkblue')) +
#   scale_shape_manual(name="Significance", values=c(18, 19)) +
#   geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
#   ylim(-.5, .5) +
#   # scale_y_continuous(
#   #   trans = log10_trans(),
#   #   breaks = trans_breaks("log10", function(x) 10^x),
#   #   labels = scales::percent_format(scale = 100, accuracy = 1)) +
#   facet_wrap(~region, nrow = 2) +
#   theme(axis.text.x=element_text(angle=-45)) +
#   labs(y = 'Coefficient',
#        x = 'Chemical Species',
#        color = 'Region',
#        title = "Smoke day regressed on fraction of species concentration") +
#   coord_flip() +
#   theme_light() 
# all_region_plot
# 
# 
# 
# # ------------------------------------------------------------------------------


# return(all_results)
# 
# } # end function


# pd <- position_dodge(0.1) # move them .05 to the left and right
# plot <- ggplot(all_results, aes(x = lhs, y = log(coef), color=model)) +
#   geom_point(aes(shape=model),size=3, alpha = 0.6, position=position_jitter(width=.5, seed = 50)) +
#   scale_color_manual(name="Model",values=c("coral","steelblue", "forestgreen", 'plum')) +
#   scale_shape_manual(name="Model",values=c(15, 17,18, 19)) +
#   # scale_x_continuous("model", breaks=1:length(model), labels=model) +
#   # scale_y_continuous("Speciation") +
#   geom_linerange(aes(ymin = coef-se, ymax = coef+se),
#                  position = position_jitter(width = 0.5, seed = 50)) +
#   #geom_errorbar(aes(ymin=coef-se,ymax=coef+se),width=0.1,position=pd) +
#   theme(axis.text.x=element_text(angle=-45)) +
#   ggtitle("Regression Results") +
#   labs(y = 'Coefficient',
#        x = 'Chemical Species',
#        color = 'Model',
#        title = "Speciation of Smoke PM2.5") +
#   coord_flip() +
#   theme_bw() 
# plot

# create scatter of estimates w/standard error bars
# coeff_plot <- ggplot(all_results, 
#                        aes(x = lhs, y = coef)) +
#   geom_point(aes(x = lhs, y = coef, color = model), size = 1, position = "jitter", alpha = .4) +
#   scale_colour_manual(values= col) +
#   
#   # add poulation density base
#   geom_hline(yintercept=0, linetype="dashed", lwd=.5, color='grey') +
#   ggtitle("Regression Results") +
#   labs(y = 'Coefficient',
#        x = 'Chemical Species',
#        color = 'Model',
#        title = "Speciation of Smoke PM2.5") +
#   coord_flip() +
#   theme_bw() 
#   # theme(plot.title =element_text(size=14, face='bold')) +
#   # # theme(plot.subtitle = element_text(size=12)) +
#   # theme(axis.text=element_text(size=11),
#   #       axis.title=element_text(size=12,face="bold")) +
#   # + theme(legend.text=element_text(size=14))
# 
# coeff_plot

# calculate smokePM using plumes and totalPM by monitors data
# to compare the regression results with those using Marissa's)
# df <- reg_df %>%
#   dplyr::select(SiteCode, Date, smoke_day, MF = 'conc_MF', RC_PM2 = 'conc_RC_PM2') %>% 
#   dplyr::filter(smoke_day == 0) 
# df$yr_month_chr <- format(df$Date, "%Y-%m")
# 
# # Calculate 3-year median total PM2.5 using MF and RC_PM2 respectively ----
# pm_pre <- df %>% mutate(Date_pre = Date - years(1)) %>% 
#   dplyr::select(SiteCode, Date, Date_pre,smoke_day) %>% 
#   dplyr::filter(smoke_day == 0) 
# 
# pm_pre <- dplyr::left_join(pm_pre, df, by=c('SiteCode'='SiteCode', 'Date_pre'='Date', 'smoke_day'='smoke_day'))%>% 
#   rename(MF_pre = MF, RC_PM2_pre = RC_PM2)
# 
# pm_post <- df %>% mutate(Date_post = Date + years(1)) %>% 
#   dplyr::select(SiteCode, Date, Date_post,smoke_day) %>% 
#   dplyr::filter(smoke_day == 0)
# 
# pm_post <- dplyr::left_join(pm_post, df, by=c('SiteCode'='SiteCode', 'Date_post'='Date', 'smoke_day'='smoke_day')) %>% 
#   rename(MF_post = 'MF', RC_PM2_post = 'RC_PM2')
# 
# # Join altogether
# df <- dplyr::left_join(df, pm_pre, by=c('SiteCode', 'Date'))
# df <- dplyr::left_join(df, pm_post, by=c('SiteCode', 'Date'))
# df <- df %>% dplyr::select(SiteCode, Date, yr_month_chr.x, MF, RC_PM2, Date_pre,
#                            MF_pre,RC_PM2_pre, Date_post, MF_post,RC_PM2_post) %>% 
#   mutate(yr_month_pre = format(Date_pre, "%Y-%m"),
#          yr_month_post = format(Date_post, "%Y-%m")) %>% 
#   rename(yr_month = yr_month_chr.x) %>% 
#   dplyr::select(SiteCode, Date, yr_month, MF, RC_PM2, yr_month_pre,
#                 MF_pre,RC_PM2_pre, yr_month_post, MF_post,RC_PM2_post)
# 
# df_gb <- df %>% group_by(SiteCode, yr_month) %>% summarize(MF = sum(MF),
#                                                            RC_PM2 = sum(RC_PM2),
#                                                            MF_pre = sum(MF_pre),
#                                                            RC_PM2_pre = sum(RC_PM2_pre),
#                                                            MF_post = sum(MF_post),
#                                                            RC_PM2_post = sum(RC_PM2_post))
# df_gb <- df_gb %>% rowwise() %>% 
#   mutate(MF_med = median(c(MF, MF_pre, MF_post), na.rm = TRUE),
#          RC_PM2_med = median(c(RC_PM2, RC_PM2_pre, RC_PM2_post, na.rm = TRUE)))
# 
# # Add the calculated 3-year median to the original dataset ----
# monitor <- reg_df
# monitor$yr_month <- format(monitor$Date, "%Y-%m")
# 
# df_gb <- df_gb %>% dplyr::select(SiteCode, yr_month, MF_med, RC_PM2_med)
# 
# fin <- dplyr::left_join(monitor, df_gb, by=c('SiteCode', 'yr_month'))
# 
# # calculated smoke 
# smoke_comp_df <- fin %>% 
#   mutate(MF_anom = conc_MF - MF_med,
#          RC_PM2_anom = conc_RC_PM2 - RC_PM2_med,
#          smokePM_MF = ifelse(smoke_day, pmax(MF_anom, 0), 0),
#          smokePM_RC = ifelse(smoke_day, pmax(RC_PM2_anom, 0), 0))
# 
# rm(df_gb, df, fin, monitor, pm_post, pm_pre)        
# 
# # Compare the results for these two different measures of smoke ----------------
# smoke_comp_reg = feols(c(conc_NH4, conc_tot_metals, 
#                          conc_ammNO3, conc_ammSO4, conc_AL,conc_AS,conc_BR,conc_CA,
#                          conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR,
#                          conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
#                          conc_MO, conc_NA,conc_NI,conc_NO3,conc_N2,
#                          conc_P,conc_PB,conc_RB,conc_S,conc_SE,
#                          conc_SI,conc_SOIL,conc_SO4,conc_SR,conc_V,
#                          conc_ZN,conc_ZR) ~ smokePM_MF | SiteCode + year + month, smoke_comp_df)
