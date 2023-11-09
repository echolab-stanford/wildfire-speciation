# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species


# loadd(c(parameter_categories, pm_pal, clean_PMspec_df), cache = drake::drake_cache(".drake"))

# function
create_attributable_frac_time_series <- function(clean_PMspec_df) {

  # select the vars that are needed for the regression
  reg_df <- clean_PMspec_df %>% 
    dplyr::select(Dataset, state_name, region, year, month, Date, 
                  monitor_month, smoke_day, site_id, MF_adj, smokePM, 
                  nonsmokePM_MF, AL:ZR) 
  
################################################################################
# How much is concentration due to smoke changing over time?
################################################################################
  # -----------------------------------------------------------------------------
  # RUN REGRESSION FOR SPECIES
  #  A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
  # -----------------------------------------------------------------------------
  # the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
  full_sampPM_regMF = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                              K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                              S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                            ~ smokePM + nonsmokePM_MF | 
                              monitor_month + year, reg_df, cluster = 'site_id')
  
  ## calculate 95 CI%
  CIs <- confint(
    full_sampPM_regMF) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'MF') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  betas <- coeftable(full_sampPM_regMF) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    dplyr::select(-id, -pval, -`t value`)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(spec_type = fct_relevel(spec_type,
                                   c("Heavy metal", "Toxic metal", "Toxic nonmetal", 
                                     "Toxicity potentiator", "Non-toxic metal", 
                                     "Non-toxic nonmetal", "Secondary inorganic", "Secondary organic"
                                   ))) %>% 
    mutate(measure = 'MF') %>% 
    filter(pm_type == 'smokePM')
  
# predict over our sample
predicted_spec_smoke_conc <- reg_df %>% 
  dplyr::select(Dataset, state_name, year, month, Date, 
                smoke_day, site_id, MF_adj, smokePM, 
                nonsmokePM_MF, AL:ZR) %>% 
  pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>% 
  left_join(betas, by = 'species') %>% 
  mutate(pred_spec_nonsmoke = Estimate*nonsmokePM_MF,
         pred_spec_smoke = Estimate*smokePM) %>% 
  mutate(pred_spec_allPM = pred_spec_smoke+pred_spec_nonsmoke) %>% 
  distinct() 

  # calculate avg attributable fraction of species concentration to total concentration
  avg_conc_site_month_yr <- predicted_spec_smoke_conc %>%
    # get avg concentration across all sites for a given month and year
    group_by(species, year, month, site_id) %>% 
    dplyr::summarise(avg_smoke_conc = mean(pred_spec_smoke, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(pred_spec_nonsmoke, na.rm = TRUE),
                     avg_tot_conc = mean(pred_spec_allPM, na.rm = TRUE)) %>% 
    ungroup() %>% 
    # get avg concentration across all sites for a given month and year
    group_by(species, year, month) %>% 
    dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                     avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>% 
    ungroup() 
  
  
  # now get the average concentration across all sites, for a given month and year
  avg_attributable_fraction_vs_obs <- avg_conc_site_month_yr %>% 
    # add in the real observed "total concentration"
    left_join(reg_df %>% 
                dplyr::select(Dataset, year, month, Date, 
                              smoke_day, site_id, MF_adj, smokePM, 
                              nonsmokePM_MF, AL:ZR) %>% 
                pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>% 
                group_by(species, year, month, site_id) %>% 
                dplyr::summarise(avg_obs_tot_conc = mean(conc_val, na.rm = TRUE)) %>% 
                ungroup() %>% 
                group_by(species, year, month) %>% 
                dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>% 
                ungroup(), 
              by = c('species', 'year', 'month')) %>% 
    filter(species != 'N2')
  
  # calculate the number of monitoring stations per month
  num_monthly_stations_by_species <- predicted_spec_smoke_conc %>%
    # get avg concentration across all sites for a given month and year
    group_by(species, year, month, site_id) %>% 
    dplyr::summarise(avg_smoke_conc = mean(pred_spec_smoke, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(pred_spec_nonsmoke, na.rm = TRUE),
                     avg_tot_conc = mean(pred_spec_allPM, na.rm = TRUE)) %>% 
    ungroup() %>% 
    group_by(species, year,  month) %>% 
    dplyr::summarise(num_stations = n_distinct(site_id)) %>% 
    ungroup() %>% 
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) 
    
  # test plotting
    # plot <- ggplot(reg_df %>%
    #                  # filter(species == 'RB') %>% 
    #                  mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d"))) +
    #   geom_line(aes(x = mon_yr,
    #                    y = nonsmokePM_MF), alpha = 0.8) +
    #   labs(title = expression("Predicted vs observed species concentrations"),
    #        y = expression("Concentration (\u00b5g/m"^3*")"),
    #        x = "",
    #   )
    # plot

  
    # prep df for plotting
    plot_df <- avg_attributable_fraction_vs_obs %>% 
      dplyr::select(-avg_obs_tot_conc, -avg_nonsmoke_conc, -avg_tot_conc) %>% 
      mutate(label = 'Attributable concentration to smoke PM2.5') %>% 
      bind_rows(avg_attributable_fraction_vs_obs %>% 
                  dplyr::select(avg_smoke_conc = 'avg_obs_tot_conc', year, month, species) %>% 
                  mutate(label = 'Observed Concetration in Total PM2.5')) %>% 
      mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
      rename(conc = 'avg_smoke_conc') %>% 
      left_join(parameter_categories, by = 'species') %>% 
      mutate(label = fct_relevel(label, c( "Observed Concetration in Total PM2.5", 
                                           "Attributable concentration to smoke PM2.5"))) 
    
    limits <- c("Organic Carbon (OC)", "Elemental carbon (EC)", # secondary organics
                #secondary inorganic
                "Sulfate (SO4)", "Nitrate (NO3)",
                # nontoxic metals
                "Strontium (Sr)",  "Calcium (Ca)", "Silicon (Si)", "Magnesium (Mg)", "Rubidium (Rb)",
                # nontoxic nonmetals
                "Soil", 
                # toxicity potentiator
                "Sulfur (S)",
                # toxic nonmetals
                "Phosphorus (P)", "Bromine (Br)", "Selenium (Se)", "Chlorine (Cl)", "Chloride (Chl)",
                # toxic metals
                "Potassium (K)",  "Zinc (Zn)", "Manganese (Mn)", "Titanium (Ti)",
                "Alumnium (Al)", "Iron (Fe)", "Copper (Cu)",
                "Vanadium (V)", "Nickel (Ni)", "Sodium (Na)", "Zirconium (Zr)",
                # heavy metals
                "Arsenic (As)", "Lead (Pb)", "Chromium (Cr)")
    
    datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")
    
    
# PLOTTING
    avg_mon_pred_spec_plot <- ggplot(plot_df) +
      geom_ribbon(aes(x = mon_yr,
                       y = conc,
                       # fill = label,
                       #color = label), 
                      stat = "identity", alpha = 0.8) +
      labs(title = expression("Predicted vs observed species concentrations"),
           y = expression("Concentration (\u00b5g/m"^3*")"),
           x = "",
      ) +
      scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
      geom_hline(aes(yintercept = 0)) +
      theme_minimal() +
      scale_color_manual(values = pm_pal) +
      scale_fill_manual(values = pm_pal) +
      facet_wrap(~factor(species_long, levels = rev(limits)), scales = 'free_y', ncol = 4) +
      theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
            axis.text = element_text(),
            axis.title = element_text(),
            legend.text = element_text(),
            legend.title = element_blank(),
            axis.ticks.x = element_line(linetype = 1),
            axis.ticks.length=unit(-0.2, "cm"),
            panel.grid = element_blank(),
            legend.position = "bottom",
            axis.line.y = element_line(color = "grey10"))  # Add y-axis line   # Add y-axis ticks
    avg_mon_pred_spec_plot
    
 
  # save file
  ggsave(
    filename = 'Fig4_smoke_attributable_fraction_trend_monthly.pdf',
    plot = avg_mon_pred_spec_plot,
    path = file.path(wip_gdrive_fp, 'figures/Fig4'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)   
  
  
  # PLOT THE NUMBER OF STATIONS ONLINE OVER THE TIME SERIES
 num_stations_recording <- ggplot(num_monthly_stations_by_species %>% 
                                    left_join(parameter_categories, by = 'species') %>% 
                                    filter(species != 'N2')) +
    geom_line(aes(x = mon_yr,
                  y = num_stations
                  ), color = 'black') +
    labs(title = expression("Average number of stations per month contributing to predictions"),
         y = "Number of stations",
         x = "",
    ) +
    scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
    geom_hline(aes(yintercept = 0)) +
    theme_minimal() +
    scale_color_manual(values = spec_pal) +
    #facet_wrap(~factor(species_long, levels = rev(limits)), scales = 'free_y', ncol = 4) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
          axis.text = element_text(),
          axis.title = element_text(),
          legend.text = element_text(),
          legend.title = element_blank(),
          axis.ticks.x = element_line(linetype = 1),
          axis.ticks.length=unit(-0.2, "cm"),
          panel.grid = element_blank(),
          #legend.position = "bottom",
          axis.line.y = element_line(color = "grey10"))  # Add y-axis line   # Add y-axis ticks
 num_stations_recording
  
 # save file
 ggsave(
   filename = 'Fig4b_avg_num_stations_recording_time_series.pdf',
   plot = num_stations_recording,
   path = file.path(wip_gdrive_fp, 'figures/Fig4'),
   height = 3,
   width =4,
   scale = 1,
   dpi = 320)   
 

  return(avg_attributable_fraction_vs_obs)
}
  
# selected_spec_df <- clean_PMspec_df %>% 
#   mutate(monitor_month = paste0(site_id,"_",month)) %>% 
#   dplyr::select(Dataset:month, monitor_month,Date, site_id, MF_adj, smokePM, nonsmokePM_MF, AL:ammSO4) #%>% 
#   # only keep rows that have all values so that its the same selection of sites for the analysis
#   #na.omit()  
# 
# # run the regressions ----------------------------------------------------------
# V_reg = feols(V ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# PB_reg = feols(PB ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# FE_reg = feols(FE ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# CU_reg = feols(CU ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# CR_reg = feols(CR ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# MN_reg = feols(MN ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# ZN_reg = feols(ZN ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# NI_reg = feols(NI ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# S_reg = feols(S ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# EC_reg = feols(EC ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# OC_reg = feols(OC ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
# 
# 
# # NOW USE THE RESULTS OF EACH REGRESSION TO PREDICT
# nonsmoke0_df <- selected_spec_df %>%
#   # 0 out nonsmoke
#   mutate(nonsmokePM = 0) 
# 
# smoke0_df <- selected_spec_df %>%
#   # 0 out nonsmoke
#   mutate(smokePM = 0)
# 
# # now predict with new data
# # predict the smoke concentration
# predicted_smoke_conc_df <- selected_spec_df %>% 
#   # this gets the fraction of a species attributable to smoke
#   mutate(V = predict(V_reg, newdata = nonsmoke0_df)) %>% 
#   mutate(PB = predict(PB_reg, newdata = nonsmoke0_df)) %>%
#   mutate(FE = predict(FE_reg, newdata = nonsmoke0_df)) %>%
#   mutate(CU = predict(CU_reg, newdata = nonsmoke0_df)) %>%
#   mutate(CR = predict(CR_reg, newdata = nonsmoke0_df)) %>%
#   mutate(MN = predict(MN_reg, newdata = nonsmoke0_df)) %>%
#   mutate(ZN = predict(ZN_reg, newdata = nonsmoke0_df)) %>%
#   mutate(NI = predict(NI_reg, newdata = nonsmoke0_df)) %>%
#   mutate(S = predict(S_reg, newdata = nonsmoke0_df)) %>%
#   mutate(EC = predict(EC_reg, newdata = nonsmoke0_df)) %>%
#   mutate(OC = predict(OC_reg, newdata = nonsmoke0_df)) %>% 
#   pivot_longer(cols = AL:ammSO4, 
#                names_to = 'species', 
#                values_to = 'smoke_conc') 
#  
# # now predict with new data
# predicted_nonsmoke_conc_df <- selected_spec_df %>% 
#   # this gets the fraction of a species attributable to nonsmoke
#   mutate(V = predict(V_reg, newdata = smoke0_df)) %>% 
#   mutate(PB = predict(PB_reg, newdata = smoke0_df)) %>%
#   mutate(FE = predict(FE_reg, newdata = smoke0_df)) %>%
#   mutate(CU = predict(CU_reg, newdata = smoke0_df)) %>%
#   mutate(CR = predict(CR_reg, newdata = smoke0_df)) %>%
#   mutate(MN = predict(MN_reg, newdata = smoke0_df)) %>%
#   mutate(ZN = predict(ZN_reg, newdata = smoke0_df)) %>%
#   mutate(NI = predict(NI_reg, newdata = smoke0_df)) %>%
#   mutate(S = predict(S_reg, newdata = smoke0_df)) %>%
#   mutate(EC = predict(EC_reg, newdata = smoke0_df)) %>%
#   mutate(OC = predict(OC_reg, newdata = smoke0_df)) %>% 
#   # pivot longer
#   pivot_longer(cols = AL:ammSO4, 
#                names_to = 'species', 
#                values_to = 'nonsmoke_conc')
# 
# # calculate total concentration in PM2.5 
# predicted_concs <- predicted_smoke_conc_df %>% 
#   left_join(predicted_nonsmoke_conc_df) %>% 
#   filter(nonsmoke_conc > 0) %>% 
#   filter(smoke_conc > 0) %>% 
#   mutate(tot_pred_conc = nonsmoke_conc + smoke_conc)


# rm(predicted_smoke_conc_df, predicted_nonsmoke_conc_df, nonsmoke0_df, smoke0_df) # drop to save memory




