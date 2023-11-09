# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Sept 27, 2023
# Description: This script matches the CAS numbers for the parameters in our data to identify potential health risks
# create a crosswalk of parameters to CAS numbers to identify health risks
# useful infor: https://www.epa.gov/risk/conducting-human-health-risk-assessment
# https://www.epa.gov/expobox/exposure-assessment-tools-routes-inhalation
#  
improve_params_fp = file.path(wip_gdrive_fp, 'raw/IMPROVE_parameters.csv')
CAS_nos_fp = file.path(wip_gdrive_fp, 'raw/CSN_parameters.csv')
RSL_table_fp = file.path(wip_gdrive_fp, 'raw/RSL Table/chronic_RSL.xlsx')
loadd(clean_PMspec_df, cache = drake::drake_cache(".drake"))
avg_mon_pred_species <- read.csv(file.path(wip_gdrive_fp, 'intermediate/avg_monthly_attributable_conc.csv'))

# create a crosswalk between parameters and CAS numbers and toxicity thresholds
#combine_speciation_w_CASno_and_RSL <- function(improve_params_fp, CAS_nos_fp, RSL_table_fp, clean_PMspec_df)

# read in all the possible parameters in the improve data (which we use to determine what parameters to keep in CSN data)
improve_params <- read.csv(improve_params_fp) %>% 
  dplyr::select(AQSCode, ParamName, ParamCode)
  
# read in CAS numbers
CAS_nos <- read.csv(file.path(CAS_nos_fp), na.strings ="") %>% 
  filter(!is.na(CAS.Number)) %>% 
  filter(Parameter.Code %in% improve_params$AQSCode) %>% 
  dplyr::select(AQSCode = 'Parameter.Code', CAS.Number) %>% 
  left_join(improve_params, by = 'AQSCode')
 

# Calculating the Potential Dose for intake processes:

# Average concentration in air 
# C = 0.0005620278 # this it the predicted avg concentration in totPM2.5 over sample period for Arsenic
# ED = 15*365*24 # 15 years -> hours
# 
# Dpot = C*IR*ED # average potential dose across US and across sample
# BW = 70 # avg bodyweight
# AT = 15*365
# 
# ADDpot = Dpot/ (BW * AT)
# 
# Risk  = 0.000015 * C 

RSL_chronic <- read_xlsx(file.path(RSL_table_fp), sheet = 4) %>% 
  filter(`CAS No.` %in% CAS_nos$CAS.Number) %>% 
  dplyr::select(Analyte, CAS_num = `CAS No.`, IUR_ug_m3 = `IUR\r\n(ug/m3)-1`, RfCi_mg_m3 = `RfCi(mg/m3)`, 
                carcSL_ug_m3 = `Carcinogenic SL\r\nTR=1E-06\r\n(ug/m3)`, 
                noncarcSL_ug_m3 = `Noncarcinogenic SL\r\nTHI=1\r\n(ug or fibers/m3)`) %>% 
  # convert to ug/m3
  mutate(RfCi_ug_m3 = RfCi_mg_m3*1000) %>% 
  dplyr::select(-RfCi_mg_m3) %>% 
  mutate(species = case_when(
    Analyte == "Aluminum" ~ 'AL',
    Analyte == "Arsenic, Inorganic" ~ 'AS',
    Analyte == "Chlorine" ~ 'CL',
    Analyte == "Chromium, Total" ~ 'CR',
    Analyte == "~Lead and Compounds" ~ 'PB',
    Analyte == "Manganese (Diet)" ~ 'MN',
    Analyte == "Copper" ~ 'CU',
    Analyte == "Iron" ~ 'FE',
    Analyte == "Manganese (Non-diet)" ~ 'MN',
    Analyte == "Nickel Soluble Salts" ~ 'NI',
    Analyte == "Nitrate (measured as nitrogen)" ~ 'NO3',
    Analyte == "Phosphorus, White" ~ 'P',
    Analyte == "Selenium" ~ 'SE',
    Analyte == "Strontium, Stable" ~ 'SR',
    Analyte == "Zinc and Compounds" ~ 'ZN',
    Analyte == "Vanadium and Compounds" ~ 'V',
    Analyte == "Zirconium" ~ 'ZR')) 
 
avg_daily_pred_conc_w_RSL <- avg_conc_daily_CONUS %>% 
  left_join(RSL_chronic, by = 'species')

 
species_list_w_RSL <- unique(avg_daily_pred_conc_w_RSL$species)
#current_species <-species_list_w_RSL[6]

all_toxic_exposures_df <- map_df(species_list_w_RSL, function(current_species) {
  
  current_spec_df <- avg_daily_pred_conc_w_RSL %>% 
    filter(species == current_species) %>% 
    mutate(color_id = case_when(
      species == 'AL' ~ "#B40F20",
      species == 'AS' ~ "#E2D200",
      species == 'CL' ~ "#0B775E",
      species == 'MN' ~ "#B40F20", 
      species == 'NI' ~ "#B40F20",
      species == 'PB' ~ "#E2D200",
      species == 'SE' ~ "#35274A",
      species == 'V' ~ "#B40F20"))

  
  # plot time series 
  current_threshold_timeseries <- ggplot(current_spec_df) +
    geom_line(aes(x = Date,
                  y = avg_tot_conc),
              color = 'grey60', alpha = 0.6, stat = "identity") +
    geom_line(aes(x = Date,
                  y = avg_smoke_conc), color = "goldenrod", alpha = 0.6, stat = "identity", show.legend = TRUE) +
    #scale_color_manual(values=c(unique(current_spec_df$color_id))) +
    geom_hline(aes(yintercept = unique(IUR_ug_m3), color = 'IUR'), linetype = 'solid') + 
    geom_hline(aes(yintercept = unique(RfCi_ug_m3), color = "RFCi"), linetype = 'solid') + 
    geom_hline(aes(yintercept = unique(carcSL_ug_m3), color = "Carc_SL"), linetype = 'solid') + 
    geom_hline(aes(yintercept = unique(carcSL_ug_m3), color = "NonCarc_SL"), linetype = 'solid') + 
    geom_hline(aes(yintercept = mean(avg_tot_conc, na.rm= TRUE), color = "avg pred. mean"), linetype = 'dashed') +
    scale_color_manual(values = c('IUR' = 'brown',
                                 'RFCi' = 'darkorange',
                                 'Carc_SL' = 'red',
                                 'NonCarc_SL' = 'tomato2',
                                 'avg pred. mean' = 'black'),
                      labels = c('IUR',
                                 'RFCi',
                                 'Carc_SL',
                                 'NonCarc_SL',
                                 'avg pred. mean'),
                      guide = guide_legend(title = "Threshold Type")) +
    labs(y = 'Avg Predicted Daily Concentration (ug/m3)',
         color = 'Species Category', 
         title = paste("Toxicity Thresholds for", unique(current_spec_df$Analyte))) +
    # facet_wrap(~region) +
    theme_minimal() +
    guides(color = guide_legend(title = "Threshold Type",
                                override.aes = list(color = c("brown", "darkorange", 'red','tomato2', 'black'))))
  current_threshold_timeseries
  
  ggsave(
    filename = paste0('Fig6_',unique(current_spec_df$species), '_toxicity_plot.pdf'),
    plot = current_threshold_timeseries,
    path = file.path(wip_gdrive_fp, 'figures/Fig6'),
    scale = 1,
    width = 10,
    height = 6,
    dpi = 320) 
  
  current_spec_df
  
  # Compare how many days above threshold with smoke and without smoke
  
  
}) %>% 
  bind_rows




mean_AS_totPM <- mean(predicted_concs$tot_pred_conc, na.rm = TRUE)

# count <- predicted_concs %>% group_by(over_under_flag) %>% tally()



# AS_threshold_timeseries <- ggplot(avg_conc_daily_CONUS) +
#   geom_line(aes(x = Date,
#                 y = avg_tot_conc),
#             color = 'grey30', alpha = 0.6, stat = "identity") +
#   geom_line(aes(x = Date,
#                 y = avg_smoke_conc),
#             color = 'goldenrod', alpha = 0.6, stat = "identity") +
#   geom_hline(yintercept = 0.000015, col = "red", linetype = 'dashed') + # cancer screening level
#   geom_hline(yintercept = mean_AS_totPM, col = "grey10", linetype = 'solid') + # non cancer screening level
#   facet_wrap(~region, scales = 'free_y') +
#   theme_light()
#   AS_threshold_timeseries

# # select the vars that are needed for the
# reg_df <- clean_PMspec_df %>%
#   mutate(monitor_month = paste0(site_id,"_", month)) %>% 
#   dplyr::select(Dataset, Date, year, month, monitor_month, site_id, state_name, region,
#                 MF_adj, smokePM, nonsmokePM_MF, AS, MN, NI, V, CL, PB, AL, SE) 
# 
# 
# # # run the regressions ----------------------------------------------------------
# AS_reg = feols(AS ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# MN_reg = feols(MN ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# NI_reg = feols(AS ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# V_reg = feols(V ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# CL_reg = feols(CL ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# PB_reg = feols(PB ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# AL_reg = feols(AL ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# SE_reg = feols(SE ~ smokePM + nonsmokePM_MF | monitor_month + year, reg_df, cluster = 'site_id')
# 
# 
# # NOW USE THE RESULTS OF EACH REGRESSION TO PREDICT
# nonsmoke0_df <- reg_df %>%
#   # 0 out nonsmoke
#   mutate(nonsmokePM_MF = 0)
# 
# smoke0_df <- reg_df %>%
#   # 0 out nonsmoke
#   mutate(smokePM = 0)
# 
# # now predict with new data
# # predict the smoke concentration
# predicted_smoke_conc_df <- reg_df %>%
#   # this gets the fraction of a species attributable to smoke
#   mutate(AS_smoke_pred = predict(AS_reg, newdata = nonsmoke0_df),
#          MN_smoke_pred = predict(MN_reg, newdata = nonsmoke0_df),
#          NI_smoke_pred = predict(NI_reg, newdata = nonsmoke0_df),
#          V_smoke_pred = predict(V_reg, newdata = nonsmoke0_df),
#          CL_smoke_pred = predict(CL_reg, newdata = nonsmoke0_df),
#          PB_smoke_pred = predict(PB_reg, newdata = nonsmoke0_df),
#          AL_smoke_pred = predict(AL_reg, newdata = nonsmoke0_df),
#          SE_smoke_pred = predict(SE_reg, newdata = nonsmoke0_df)) %>% 
#   dplyr::select(Dataset:site_id, month, region, monitor_month, MF_adj, smokePM, nonsmokePM_MF, 
#                 AS_smoke_pred, MN_smoke_pred, NI_smoke_pred, V_smoke_pred,
#                 CL_smoke_pred,PB_smoke_pred, AL_smoke_pred, SE_smoke_pred) %>% 
#   pivot_longer(cols = c(AS_smoke_pred, MN_smoke_pred, NI_smoke_pred, V_smoke_pred,
#                         CL_smoke_pred,PB_smoke_pred, AL_smoke_pred, SE_smoke_pred),
#                names_to = 'species', values_to = "smokePM_pred_conc") %>% 
#   filter(smokePM_pred_conc > 0) %>% 
#   mutate(species= str_remove(species, "_smoke_pred"))
# 
# 
# # now predict with new data
# predicted_nonsmoke_conc_df <- reg_df %>%
#   # this gets the fraction of a species attributable to nonsmoke
#   mutate(AS_nonsmoke_pred = predict(AS_reg, newdata = smoke0_df),
#          MN_nonsmoke_pred = predict(MN_reg, newdata = smoke0_df),
#          NI_nonsmoke_pred = predict(NI_reg, newdata = smoke0_df),
#          V_nonsmoke_pred = predict(V_reg, newdata = smoke0_df),
#          CL_nonsmoke_pred = predict(CL_reg, newdata = smoke0_df),
#          PB_nonsmoke_pred = predict(PB_reg, newdata = smoke0_df),
#          AL_nonsmoke_pred = predict(AL_reg, newdata = smoke0_df),
#          SE_nonsmoke_pred = predict(SE_reg, newdata = smoke0_df)) %>% 
#   dplyr::select(Dataset:site_id, month, region, monitor_month, MF_adj, smokePM, nonsmokePM_MF, 
#                 AS_nonsmoke_pred, MN_nonsmoke_pred, NI_nonsmoke_pred, V_nonsmoke_pred,
#                 CL_nonsmoke_pred,PB_nonsmoke_pred, AL_nonsmoke_pred, SE_nonsmoke_pred) %>% 
#   pivot_longer(cols = c(AS_nonsmoke_pred, MN_nonsmoke_pred, NI_nonsmoke_pred, V_nonsmoke_pred,
#                         CL_nonsmoke_pred,PB_nonsmoke_pred, AL_nonsmoke_pred, SE_nonsmoke_pred),
#                names_to = 'species', values_to = "nonsmokePM_pred_conc") %>% 
#   filter(nonsmokePM_pred_conc > 0) %>% 
#   mutate(species= str_remove(species, "_nonsmoke_pred"))
# 
# # rm(MN_reg, NI_reg, PB_reg, CL_reg, AS_reg, AL_reg, nonsmoke0_df, smoke0_df, reg_df, SE_reg, V_reg) # drop to save memory
# 
# # calculate total concentration in PM2.5
# predicted_concs <- predicted_smoke_conc_df %>%
#   left_join(predicted_nonsmoke_conc_df) %>% 
#   mutate(totPM_conc = smokePM_pred_conc + nonsmokePM_pred_conc)
# 
# 
# # calculate avg attributable fraction of species concentration to total concentration
# avg_conc_daily_CONUS <- predicted_concs %>%
#   # first average at each site then average across all days
#   group_by(year, month, Date, region, site_id, species) %>%
#   dplyr::summarise(avg_smoke_conc = mean(smokePM_pred_conc, na.rm = TRUE),
#                    avg_nonsmoke_conc = mean(nonsmokePM_pred_conc, na.rm = TRUE),
#                    avg_tot_conc = mean(totPM_conc, na.rm = TRUE)) %>%
#   ungroup() %>%
#   # now average across all sites for a given date
#   group_by(year, month, Date, species) %>% # can add in regions
#   dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
#                    avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
#                    avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>%
#   ungroup() 
#   

