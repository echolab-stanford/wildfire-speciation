# # Emma Krasovich Southworth, emmars@stanford.edu
# # Last Updated: Sept 27, 2023
# # Description: This script matches the CAS numbers for the parameters in our data to identify potential health risks
# # create a crosswalk of parameters to CAS numbers to identify health risks
# # useful infor: https://www.epa.gov/risk/conducting-human-health-risk-assessment
# # https://www.epa.gov/expobox/exposure-assessment-tools-routes-inhalation
# #  
# improve_params_fp = file.path(data_fp, 'raw/IMPROVE_parameters.csv')
# csn_param_fp = file.path(data_fp, 'raw/CSN_parameters.csv')
# RSL_table_fp = file.path(data_fp, 'raw/RSL Table/chronic_RSL.xlsx')
# loadd(c(full_samp_PMcoeffs_normalized), cache = drake::drake_cache(".drake"))
# 
# 
# # create a crosswalk between parameters and CAS numbers and toxicity thresholds
# #combine_speciation_w_CASno_and_RSL <- function(improve_params_fp, CAS_nos_fp, RSL_table_fp, clean_PMspec_df)
# CAS_nos <- read.csv(csn_param_fp) %>%
#   mutate(flag = ifelse(
#     str_detect(Parameter, "2.5"), 'keep', 'drop'
#   )) %>%
#   filter(flag == 'keep') %>%
#   mutate(flag = ifelse(
#     str_detect(Parameter, "10"), 'drop', 'keep'
#   )) %>%
#   filter(flag == 'keep') %>%
#   mutate(flag = ifelse(
#     str_detect(Parameter, "(STP|TSP)"), 'drop', 'keep'
#   )) %>%
#   filter(flag == 'keep') %>%
#   mutate(flag = ifelse(
#     str_detect(Parameter, "nm"), 'drop', 'keep'
#   )) %>%
#   filter(flag == 'keep') %>%
#   filter(Standard.Units == 'Micrograms/cubic meter (LC)') %>%
#   dplyr::select(AQSParamCode = 'Parameter.Code', Parameter, CAS.Number) %>%
#   filter(Parameter %in% c("Aluminum PM2.5 LC", "Ammonium Nitrate PM2.5 LC",
#                           "Ammonium Sulfate PM2.5 LC","Arsenic PM2.5 LC","Bromine PM2.5 LC",
#                           "Calcium PM2.5 LC", "Chloride PM2.5 LC", "Chlorine PM2.5 LC", "Chromium PM2.5 LC",
#                           "Chromium VI PM2.5 LC", "Copper PM2.5 LC", "EC CSN PM2.5 LC TOT",
#                           "EC PM2.5 LC TOR", "EC PM2.5 LC TOT","Iron PM2.5 LC", "Lead PM2.5 LC", "Magnesium PM2.5 LC",
#                           "Manganese PM2.5 LC", "Nickel PM2.5 LC","OC CSN Unadjusted PM2.5 LC TOT",
#                           "OC PM2.5 LC TOR", "OC PM2.5 LC TOT", "OCH PM2.5 LC TOT","Organic Carbon Mass PM2.5 LC",
#                           "Phosphorus PM2.5 LC", "PM2.5 - Local Conditions",
#                           "Potassium PM2.5 LC", "Reconstructed Mass PM2.5 LC",
#                           "Rubidium PM2.5 LC", "Selenium PM2.5 LC", "Silicon PM2.5 LC",
#                           "Sodium PM2.5 LC", "Strontium PM2.5 LC", "Sulfate PM2.5 LC",
#                           "Sulfur PM2.5 LC",  "Titanium PM2.5 LC", "Total Nitrate PM2.5 LC",
#                           "Vanadium PM2.5 LC","Zinc PM2.5 LC")) %>%
#   mutate(ParamCode = case_when(
#     Parameter == "Aluminum PM2.5 LC" ~ 'AL',
#     Parameter == "Ammonium Nitrate PM2.5 LC" ~ 'ammNO3',
#     Parameter == "Ammonium Sulfate PM2.5 LC" ~ 'ammS04',
#     Parameter == "Arsenic PM2.5 LC" ~ 'AS',
#     Parameter == "Bromine PM2.5 LC" ~ 'BR',
#     Parameter == "Calcium PM2.5 LC" ~ 'CA',
#     Parameter == "Chloride PM2.5 LC" ~ 'CHL',
#     Parameter == "Chlorine PM2.5 LC"~ 'CL',
#     Parameter == "Chromium PM2.5 LC"~ 'CR',
#     Parameter == "Chromium VI PM2.5 LC" ~ 'CR VI',
#     Parameter == "Copper PM2.5 LC" ~ 'CU',
#     Parameter == "EC CSN PM2.5 LC TOT" ~ 'EC',
#     Parameter == "Iron PM2.5 LC"~ 'FE',
#     Parameter == "EC PM2.5 LC TOR"~ 'EC',
#     Parameter == "EC PM2.5 LC TOT"~ 'EC',
#     Parameter == "Lead PM2.5 LC"~ 'PB',
#     Parameter == "Magnesium PM2.5 LC"~ 'MG',
#     Parameter == "Manganese PM2.5 LC" ~ 'MN',
#     Parameter == "Nickel PM2.5 LC" ~ 'NI',
#     Parameter == "OCH PM2.5 LC TOT"~ 'OC',
#     Parameter == "Rubidium PM2.5 LC" ~ 'RB',
#     Parameter == "Zirconium PM2.5 LC"~ 'ZR',
#     Parameter == "Zinc PM2.5 LC"~ 'ZN',
#     Parameter == "OC PM2.5 LC TOR"~ 'OC',
#     Parameter == "Organic Carbon Mass PM2.5 LC"~ 'OC',
#     Parameter == "OC PM2.5 LC TOT"~ 'OC',
#     Parameter == "Phosphorus PM2.5 LC"~ 'P',
#     Parameter == "PM2.5 - Local Conditions"~ 'MF',
#     Parameter == "Potassium PM2.5 LC"~ 'K',
#     Parameter == "Reconstructed Mass PM2.5 LC" ~ 'RCFM',
#     Parameter == "Selenium PM2.5 LC" ~ 'SE',
#     Parameter == "Silicon PM2.5 LC" ~ 'SI',
#     Parameter == "Sodium PM2.5 LC"~ 'NA',
#     Parameter == "Strontium PM2.5 LC"~ 'SR',
#     Parameter == "Sulfate PM2.5 LC" ~ 'SO4',
#     Parameter == "Titanium PM2.5 LC"~ 'TI',
#     Parameter == "Sulfur PM2.5 LC"~ 'S',
#     Parameter == "Vanadium PM2.5 LC" ~ 'V',
#     Parameter == "Total Nitrate PM2.5 LC" ~ 'NO3'))
# 
# # merge CAS numbers into species list
# chem_df <- tibble(ParamCode = c('AL', 'AS','BR', 'CA', 'CL', 'CHL','CR', 'CU',
#                                 'EC', 'FE', 'K', 'MG','MN', 'NA', 'NI', 'NO3',
#                                 'OC', 'P',  'PB', 'RB', 'S',  'SE', 'SI', 'SO4',
#                                 'SR', 'TI', 'V', 'ZN')) %>%
#   left_join(CAS_nos %>% distinct(AQSParamCode, ParamCode, CAS.Number), by = 'ParamCode') %>%
#   mutate(CAS.Number = ifelse(CAS.Number == "", NA, CAS.Number)) %>%
#   filter(!is.na(CAS.Number))
# 
# # read in RSL table information
# RSL_chronic <- read_xlsx(file.path(RSL_table_fp), sheet = 4) %>%
#   filter(`CAS No.` %in% chem_df$CAS.Number) %>%
#   dplyr::select(Analyte, CAS_num = `CAS No.`, IUR_ug_m3 = `IUR\r\n(ug/m3)-1`, RfCi_mg_m3 = `RfCi(mg/m3)`,
#                 carcSL_ug_m3 = `Carcinogenic SL\r\nTR=1E-06\r\n(ug/m3)`,
#                 noncarcSL_ug_m3 = `Noncarcinogenic SL\r\nTHI=1\r\n(ug or fibers/m3)`) %>%
#   # convert to ug/m3
#   mutate(RfCi_ug_m3 = RfCi_mg_m3*1000) %>%
#   dplyr::select(-RfCi_mg_m3) %>%
#   # drop a row if reference doses for all chemicals are NA
#   filter_at(vars(IUR_ug_m3, carcSL_ug_m3, noncarcSL_ug_m3, RfCi_ug_m3),
#             any_vars(!is.na(.))) %>%
#   mutate(species = case_when(
#     Analyte == "Aluminum" ~ 'AL',
#     Analyte == "Arsenic, Inorganic" ~ 'AS',
#     Analyte == "Chlorine" ~ 'CL',
#     Analyte == "Chromium, Total" ~ 'CR',
#     Analyte == "~Lead and Compounds" ~ 'PB',
#     Analyte == "Manganese (Diet)" ~ 'MN',
#     Analyte == "Copper" ~ 'CU',
#     Analyte == "Iron" ~ 'FE',
#     Analyte == "Manganese (Non-diet)" ~ 'MN',
#     Analyte == "Nickel Soluble Salts" ~ 'NI',
#     Analyte == "Nitrate (measured as nitrogen)" ~ 'NO3',
#     Analyte == "Phosphorus, White" ~ 'P',
#     Analyte == "Selenium" ~ 'SE',
#     Analyte == "Strontium, Stable" ~ 'SR',
#     Analyte == "Zinc and Compounds" ~ 'ZN',
#     Analyte == "Vanadium and Compounds" ~ 'V',
#     Analyte == "Zirconium" ~ 'ZR'))
# 
# 
# 
# # # use the gridded predictions to see if any species are over 
# # 
# # # DETERMINE COUNTERFACTUAL
# # # how many days absent wildfire were we over the threshold?
# # # how many days did wildfire contribute?
# # counterfactual <- predicted_spec_smoke_conc %>% 
# #   dplyr::select(Dataset, year, month, region, Date, site_id, MF_adj, smokePM, nonsmokePM_MF, 
# #                 species_type, species_long, species_name, species, conc_val, pred_spec_smoke,
# #                 IUR_ug_m3, carcSL_ug_m3, noncarcSL_ug_m3, RfCi_ug_m3) %>% 
# #   mutate(absent_wildfire_conc = conc_val - pred_spec_smoke) %>% 
# #   mutate(above_IUR_nosmoke_flag = ifelse(
# #     absent_wildfire_conc >= IUR_ug_m3, TRUE, FALSE)) %>% 
# #   mutate(above_IUR_smoke_flag = ifelse(
# #     pred_spec_smoke >= IUR_ug_m3, TRUE, FALSE)) %>% 
# #   mutate(above_IUR_tot_conc = ifelse(
# #     conc_val >= IUR_ug_m3, TRUE, FALSE)) %>% 
# #   distinct()
# # 
# # 
# # days_above_IUR <- counterfactual %>%
# #   group_by(species) %>%
# #   dplyr::summarise(tot_num_days = n_distinct(Date, site_id)) %>%
# #   ungroup() %>%
# #   group_by(species) %>%
# #   left_join(
# #     test <- counterfactual %>%
# #       group_by(species) %>%
# #       dplyr::summarise(smoke_days = sum(above_IUR_smoke_flag, na.rm = TRUE),
# #                        nonsmoke_days = sum(above_IUR_nosmoke_flag, na.rm = TRUE),
# #                        totconc_days = sum(above_IUR_tot_conc, na.rm = TRUE)) %>%
# #       ungroup(),
# #     by = 'species')
# # 
# # 
# # 
# # 
# # 
# # 
# # # 
# # # 
# # # 
# # # 
# # # #  -----------------------------------
# # # # plot each species daily plot -----------------------------------
# # # species_list_w_RSL <- unique(predicted_spec_smoke_conc$species)
# # # # current_species <-species_list_w_RSL[3]
# # # 
# # # all_toxic_exposures_df <- map_df(species_list_w_RSL, function(current_species) {
# # #   
# # #   current_spec_df <- predicted_spec_smoke_conc %>% 
# # #     filter(species == current_species) %>% 
# # #     mutate(color_id = case_when(
# # #       species == 'AL' ~ "#B40F20",
# # #       species == 'AS' ~ "#E2D200",
# # #       species == 'CL' ~ "#0B775E",
# # #       species == 'MN' ~ "#B40F20", 
# # #       species == 'NI' ~ "#B40F20",
# # #       species == 'PB' ~ "#E2D200",
# # #       species == 'SE' ~ "#35274A",
# # #       species == 'V' ~ "#B40F20"))
# # # 
# # #   max_lim = max(current_spec_df$pred_spec_smoke)
# # #   
# # #   # plot time series 
# # #   current_threshold_timeseries <- ggplot(current_spec_df) +
# # #     geom_line(aes(x = Date,
# # #                   y = conc_val),
# # #               color = 'grey60', alpha = 0.6, stat = "identity") +
# # #     geom_line(aes(x = Date,
# # #                   y = pred_spec_smoke), color = "goldenrod", alpha = 0.6, stat = "identity", show.legend = TRUE) +
# # #     scale_y_continuous(limits =c(0, max_lim)) +
# # #     geom_hline(aes(yintercept = unique(IUR_ug_m3), color = 'IUR'), linetype = 'solid') + 
# # #     geom_hline(aes(yintercept = unique(RfCi_ug_m3), color = "RFCi"), linetype = 'solid') + 
# # #     geom_hline(aes(yintercept = unique(carcSL_ug_m3), color = "Carc_SL"), linetype = 'solid') + 
# # #     geom_hline(aes(yintercept = unique(carcSL_ug_m3), color = "NonCarc_SL"), linetype = 'solid') + 
# # #     geom_hline(aes(yintercept = mean(pred_spec_smoke, na.rm= TRUE), color = "avg pred. mean"), linetype = 'dotted') +
# # #     scale_color_manual(values = c('IUR' = 'brown',
# # #                                   'RFCi' = 'darkorange',
# # #                                   'Carc_SL' = 'red',
# # #                                   'NonCarc_SL' = 'tomato2',
# # #                                   'avg pred. mean' = 'black'),
# # #                        labels = c('IUR',
# # #                                   'RFCi',
# # #                                   'Carc_SL',
# # #                                   'NonCarc_SL',
# # #                                   'avg pred. mean'),
# # #                        guide = guide_legend(title = "Threshold Type")) +
# # #     labs(y = 'Predicted Daily Species Concentration (ug/m3)',
# # #          color = 'Species Category', 
# # #          title = paste("Toxicity Thresholds for", unique(current_spec_df$species_long))) +
# # #     # facet_wrap(~region) +
# # #     theme_minimal() +
# # #     guides(color = guide_legend(title = "Threshold Type",
# # #                                 override.aes = list(color = c("brown", 
# # #                                                               "darkorange", 
# # #                                                               'red','tomato2', 'black'))))
# # #   current_threshold_timeseries
# # #   
# # #   
# # #   ggsave(
# # #     filename = paste0('Fig6_',unique(current_spec_df$species), '_toxicity_plot.pdf'),
# # #     plot = current_threshold_timeseries,
# # #     path = file.path(results_fp, 'figures/Fig6'),
# # #     scale = 1,
# # #     width = 10,
# # #     height = 6,
# # #     dpi = 320) 
# # #   
# # #   current_spec_df
# # #   
# # # 
# # # }) %>% 
# # #   bind_rows
# # # 
# # # 
# # # 
# # # # Calculating the Potential Dose for intake processes:
# # # 
# # # # Average concentration in air 
# # # # C = 0.0005620278 # this it the predicted avg concentration in totPM2.5 over sample period for Arsenic
# # # # ED = 15*365*24 # 15 years -> hours
# # # # 
# # # # Dpot = C*IR*ED # average potential dose across US and across sample
# # # # BW = 70 # avg bodyweight
# # # # AT = 15*365
# # # # 
# # # # ADDpot = Dpot/ (BW * AT)
# # # # 
# # # # Risk  = 0.000015 * C 
# # # 
# # # 
# # #  
# # # avg_daily_pred_conc_w_RSL <- clean_PMspec_df %>% 
# # #   left_join(RSL_chronic, by = 'species')
# # # 
# # #  
# # 
# # 
# # # 
# # # 
# # # mean_AS_totPM <- mean(predicted_concs$tot_pred_conc, na.rm = TRUE)
# # # 
# # # # count <- predicted_concs %>% group_by(over_under_flag) %>% tally()
# # # 
# # # 
# # # 
# # # # AS_threshold_timeseries <- ggplot(avg_conc_daily_CONUS) +
# # # #   geom_line(aes(x = Date,
# # # #                 y = avg_tot_conc),
# # # #             color = 'grey30', alpha = 0.6, stat = "identity") +
# # # #   geom_line(aes(x = Date,
# # # #                 y = avg_smoke_conc),
# # # #             color = 'goldenrod', alpha = 0.6, stat = "identity") +
# # # #   geom_hline(yintercept = 0.000015, col = "red", linetype = 'dashed') + # cancer screening level
# # # #   geom_hline(yintercept = mean_AS_totPM, col = "grey10", linetype = 'solid') + # non cancer screening level
# # # #   facet_wrap(~region, scales = 'free_y') +
# # # #   theme_light()
# # # #   AS_threshold_timeseries
# # # 
# # # # # select the vars that are needed for the
# # # # reg_df <- clean_PMspec_df %>%
# # # #   mutate(monitor_month = paste0(site_id,"_", month)) %>% 
# # # #   dplyr::select(Dataset, Date, year, month, monitor_month, site_id, state_name, region,
# # # #                 MF_adj, smokePM, nonsmokePM_MF, AS, MN, NI, V, CL, PB, AL, SE) 
# # # # 
# # # # 
# # # # # # run the regressions ----------------------------------------------------------
# # # # AS_reg = feols(AS ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # MN_reg = feols(MN ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # NI_reg = feols(AS ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # V_reg = feols(V ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # CL_reg = feols(CL ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # PB_reg = feols(PB ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # AL_reg = feols(AL ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # SE_reg = feols(SE ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# # # # 
# # # # 
# # # # # NOW USE THE RESULTS OF EACH REGRESSION TO PREDICT
# # # # nonsmoke0_df <- reg_df %>%
# # # #   # 0 out nonsmoke
# # # #   mutate(nonsmokePM_MF = 0)
# # # # 
# # # # smoke0_df <- reg_df %>%
# # # #   # 0 out nonsmoke
# # # #   mutate(smokePM = 0)
# # # # 
# # # # # now predict with new data
# # # # # predict the smoke concentration
# # # # predicted_smoke_conc_df <- reg_df %>%
# # # #   # this gets the fraction of a species attributable to smoke
# # # #   mutate(AS_smoke_pred = predict(AS_reg, newdata = nonsmoke0_df),
# # # #          MN_smoke_pred = predict(MN_reg, newdata = nonsmoke0_df),
# # # #          NI_smoke_pred = predict(NI_reg, newdata = nonsmoke0_df),
# # # #          V_smoke_pred = predict(V_reg, newdata = nonsmoke0_df),
# # # #          CL_smoke_pred = predict(CL_reg, newdata = nonsmoke0_df),
# # # #          PB_smoke_pred = predict(PB_reg, newdata = nonsmoke0_df),
# # # #          AL_smoke_pred = predict(AL_reg, newdata = nonsmoke0_df),
# # # #          SE_smoke_pred = predict(SE_reg, newdata = nonsmoke0_df)) %>% 
# # # #   dplyr::select(Dataset:site_id, month, region, monitor_month, MF_adj, smokePM, nonsmokePM_MF, 
# # # #                 AS_smoke_pred, MN_smoke_pred, NI_smoke_pred, V_smoke_pred,
# # # #                 CL_smoke_pred,PB_smoke_pred, AL_smoke_pred, SE_smoke_pred) %>% 
# # # #   pivot_longer(cols = c(AS_smoke_pred, MN_smoke_pred, NI_smoke_pred, V_smoke_pred,
# # # #                         CL_smoke_pred,PB_smoke_pred, AL_smoke_pred, SE_smoke_pred),
# # # #                names_to = 'species', values_to = "smokePM_pred_conc") %>% 
# # # #   filter(smokePM_pred_conc > 0) %>% 
# # # #   mutate(species= str_remove(species, "_smoke_pred"))
# # # # 
# # # # 
# # # # # now predict with new data
# # # # predicted_nonsmoke_conc_df <- reg_df %>%
# # # #   # this gets the fraction of a species attributable to nonsmoke
# # # #   mutate(AS_nonsmoke_pred = predict(AS_reg, newdata = smoke0_df),
# # # #          MN_nonsmoke_pred = predict(MN_reg, newdata = smoke0_df),
# # # #          NI_nonsmoke_pred = predict(NI_reg, newdata = smoke0_df),
# # # #          V_nonsmoke_pred = predict(V_reg, newdata = smoke0_df),
# # # #          CL_nonsmoke_pred = predict(CL_reg, newdata = smoke0_df),
# # # #          PB_nonsmoke_pred = predict(PB_reg, newdata = smoke0_df),
# # # #          AL_nonsmoke_pred = predict(AL_reg, newdata = smoke0_df),
# # # #          SE_nonsmoke_pred = predict(SE_reg, newdata = smoke0_df)) %>% 
# # # #   dplyr::select(Dataset:site_id, month, region, monitor_month, MF_adj, smokePM, nonsmokePM_MF, 
# # # #                 AS_nonsmoke_pred, MN_nonsmoke_pred, NI_nonsmoke_pred, V_nonsmoke_pred,
# # # #                 CL_nonsmoke_pred,PB_nonsmoke_pred, AL_nonsmoke_pred, SE_nonsmoke_pred) %>% 
# # # #   pivot_longer(cols = c(AS_nonsmoke_pred, MN_nonsmoke_pred, NI_nonsmoke_pred, V_nonsmoke_pred,
# # # #                         CL_nonsmoke_pred,PB_nonsmoke_pred, AL_nonsmoke_pred, SE_nonsmoke_pred),
# # # #                names_to = 'species', values_to = "nonsmokePM_pred_conc") %>% 
# # # #   filter(nonsmokePM_pred_conc > 0) %>% 
# # # #   mutate(species= str_remove(species, "_nonsmoke_pred"))
# # # # 
# # # # # rm(MN_reg, NI_reg, PB_reg, CL_reg, AS_reg, AL_reg, nonsmoke0_df, smoke0_df, reg_df, SE_reg, V_reg) # drop to save memory
# # # # 
# # # # # calculate total concentration in PM2.5
# # # # predicted_concs <- predicted_smoke_conc_df %>%
# # # #   left_join(predicted_nonsmoke_conc_df) %>% 
# # # #   mutate(totPM_conc = smokePM_pred_conc + nonsmokePM_pred_conc)
# # # # 
# # # # 
# # # # # calculate avg attributable fraction of species concentration to total concentration
# # # # avg_conc_daily_CONUS <- predicted_concs %>%
# # # #   # first average at each site then average across all days
# # # #   group_by(year, month, Date, region, site_id, species) %>%
# # # #   dplyr::summarise(avg_smoke_conc = mean(smokePM_pred_conc, na.rm = TRUE),
# # # #                    avg_nonsmoke_conc = mean(nonsmokePM_pred_conc, na.rm = TRUE),
# # # #                    avg_tot_conc = mean(totPM_conc, na.rm = TRUE)) %>%
# # # #   ungroup() %>%
# # # #   # now average across all sites for a given date
# # # #   group_by(year, month, Date, species) %>% # can add in regions
# # # #   dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
# # # #                    avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
# # # #                    avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>%
# # # #   ungroup() 
# # #   
# # 
