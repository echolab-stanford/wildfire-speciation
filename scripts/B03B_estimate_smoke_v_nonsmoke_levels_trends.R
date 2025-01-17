# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species
# 
# loadd(pred_obs_joined_df,cache = drake::drake_cache(".drake"))

# function
estimate_smoke_v_nonsmoke_levels_trends <- function(pred_obs_joined_df) {
  
  # run model to determine if the concentration in absolute terms (aka levels)
  # is significantly increasing over time
  pred_smoke_levels_df <-pred_obs_joined_df %>% 
  # select only predicted concentrations in smoke
    dplyr::select(Date, species, region, year, month, pred_spec_smoke, site_id, grid_id_10km, monitor_month) %>% 
    pivot_wider(names_from = 'species', values_from = 'pred_spec_smoke') 
  
  
  print(paste('running model for smoke levels'))
  smoke_levels_inc_model = feols(c(AL,AS,BR,CA,CL, CR,CU, EC, FE, K, MG, MN, 
                           `NA`, NI, NO3, OC, P, PB, RB, S, 
                           SE, SI, SO4, SR, TI, V, ZN)
                         ~  year | monitor_month,
                         data = pred_smoke_levels_df, cluster = 'site_id')
  # etable(smoke_levels_inc_model)
  # summary(smoke_levels_inc_model)
  
  ## calculate 95 CI%
  CIsmoke <- confint(smoke_levels_inc_model) %>% 
    rename(s_CI25 = `2.5 %`,
           s_CI975 = `97.5 %`,
           species = 'lhs') %>% 
    dplyr::select(-id, -coefficient)
  
  # get coefficients and prepare for plotting
  coeffs_smoke <- coeftable(smoke_levels_inc_model) %>%
    rename(s_pval = 'Pr(>|t|)',
           s_se = 'Std. Error',
           species = 'lhs',
           pred_smoke_coeff = 'Estimate') %>%
    mutate(s_pval = round(s_pval, digits = 3)) %>%
    dplyr::select(-id, -coefficient) %>%
    left_join(CIsmoke, by = 'species') %>% 
    #mutate(across(where(is.numeric), ~ round(.x, 9))) %>% 
    mutate(s_CI = paste0("[", s_CI25, ", ", s_CI975, "]")) %>% 
    dplyr::select(species, pred_smoke_coeff, s_CI, s_pval) 
    
  
  # # RUN FOR NONSMOKE TO SEE IF THIS SIGNIFICANTLY DECREASING OVER TIME:
  # # run model to determine if the concentration in absolute terms (aka levels)
  # # is significantly increasing over time
  pred_nonsmoke_levels_df <-pred_obs_joined_df %>%
    # select only predicted concentrations in smoke
    dplyr::select(Date, species, region, year, month, pred_spec_nonsmoke, site_id, grid_id_10km, monitor_month) %>%
    pivot_wider(names_from = 'species', values_from = 'pred_spec_nonsmoke')

  print(paste('running model for nonsmoke levels'))
  nonsmoke_levels_model = feols(c(AL,AS,BR,CA,CL, CR,CU, EC, FE, K, MG, MN,
                                   `NA`, NI, NO3, OC, P, PB, RB, S,
                                   SE, SI, SO4, SR, TI, V, ZN)
                                 ~  year | monitor_month,
                                 data = pred_nonsmoke_levels_df, cluster = 'site_id')
  # etable(nonsmoke_levels_model)
  # summary(nonsmoke_levels_model)

  ## calculate 95 CI%
  CInonsmoke <- confint(nonsmoke_levels_model) %>%
    rename(ns_CI25 = `2.5 %`,
           ns_CI975 = `97.5 %`,
           species = 'lhs') %>%
    dplyr::select(-id, -coefficient)

  # get coefficients and prepare for plotting
  coeffs_nonsmoke <- coeftable(nonsmoke_levels_model) %>%
    rename(ns_pval = 'Pr(>|t|)',
           ns_se = 'Std. Error',
           species = 'lhs',
           pred_nonsmoke_coeff = 'Estimate') %>%
    mutate(ns_pval = round(ns_pval, digits = 3)) %>%
    dplyr::select(-id, -coefficient) %>%
    left_join(CInonsmoke, by = 'species') %>%
    dplyr::select(species, pred_nonsmoke_coeff, ns_pval, ns_CI25, ns_CI975) %>%
    mutate(ns_CI = paste0("[", ns_CI25, ", ", ns_CI975, "]")) %>%
    dplyr::select(species, pred_nonsmoke_coeff,
                  ns_CI,
                  ns_pval)

  # combine both predicted levels
  predicted_levels_df <- left_join(
    coeffs_nonsmoke, coeffs_smoke) %>%
    dplyr::select(-ns_CI, -s_CI) %>%
    mutate(across(where(is.numeric), ~ signif(.x, 6)))

  rm(coeffs_nonsmoke, coeffs_smoke, nonsmoke_levels_model, 
     CIsmoke, CInonsmoke, pred_nonsmoke_levels_df, pred_smoke_levels_df)
  

  return(predicted_levels_df)
}

#   # TEST IF THE ATTRIBUTABLE FRACTION IS INCREASING OVER TIME
#   print(paste('running model for fraction'))
#   att_frac_model = feols(c(AL,AS,BR,CA,CL, CR,CU, EC, FE, K, MG, MN, 
#                            `NA`, NI, NO3, OC, P, PB, RB, S, 
#                            SE, SI, SO4, SR, TI, V, ZN)
#                          ~  year | monitor_month,
#                          data = calc_frac, cluster = 'site_id')
#   # etable(att_frac_model)
#   
#   
#   ## calculate 95 CI%
#   CIs <- confint(att_frac_model) %>% 
#     rename(CI25 = `2.5 %`,
#            CI975 = `97.5 %`,
#            species = 'lhs') %>% 
#     dplyr::select(-id)
#   
#   # get coefficients and prepare for plotting
#   coeffs <- coeftable(att_frac_model) %>%
#     rename(pval = 'Pr(>|t|)',
#            se = 'Std. Error',
#            species = 'lhs',
#            coeff = 'Estimate') %>%
#     mutate(pval = round(pval, digits = 3)) %>%
#     dplyr::select(-id) %>%
#     left_join(CIs, by = c('species', 'coefficient')) %>% 
#     mutate(sig = ifelse(pval < .05, "sig increasing", "cannot detect increase"))
#   
# rm(att_frac_model)
# 
# print(paste('averaging preds for each region'))
# 
#   # combine all days with predicted dataset
#   avg_conc_site_month_yr <- pred_conc_all_days %>%
#   # calculate avg attributable fraction of species concentration to total concentration
#     # get avg concentration across all sites for a given month and year
#     group_by(region, species, year, month, site_id) %>%
#     dplyr::summarise(avg_smoke_conc = mean(pred_spec_smoke, na.rm = TRUE),
#                      avg_nonsmoke_conc = mean(pred_spec_nonsmoke, na.rm = TRUE),
#                      avg_tot_conc = mean(pred_spec_allPM, na.rm = TRUE)) %>%
#     ungroup() %>%
#     group_by(species, year, month) %>%
#     dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
#                      avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
#                      avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>%
#     ungroup()
# 
#   # get monthly avg conc for the observed "total concentration" of a species
#   obs_improve_monthly_mean <- reg_df %>%
#     dplyr::select(Dataset, region, year, month, Date,
#                   smoke_day, site_id, MF_adj, smokePM,
#                   nonsmokePM_MF, AL:ZN) %>%
#     pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>%
#     group_by(region, species, year, month, site_id) %>%
#     dplyr::summarise(avg_obs_tot_conc = mean(conc_val, na.rm = TRUE)) %>%
#     ungroup() %>%
#     group_by(region,species, year, month) %>%
#     dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>%
#     ungroup() %>%
#     group_by(species, year, month) %>%
#     dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>%
#     ungroup()
# 
#   # now merge back with predicted data
#   print(paste('merging predictions with observed data'))
#   avg_attributable_fraction_vs_obs <- avg_conc_site_month_yr %>%
#     left_join(obs_improve_monthly_mean, 
#               by = c('species', 'year', 'month'))
# 
#   # # prep df for plotting
#   plot_df <- avg_attributable_fraction_vs_obs %>%
#     dplyr::select(-avg_obs_tot_conc, -avg_nonsmoke_conc, -avg_tot_conc) %>%
#     mutate(label = 'Attributable concentration to smoke PM2.5') %>%
#     rename(conc = 'avg_smoke_conc') %>%
#     bind_rows(avg_attributable_fraction_vs_obs %>%
#                 dplyr::select(conc = 'avg_obs_tot_conc', year, month, species) %>%
#                 mutate(label = 'Observed Concentration in Total PM2.5')) %>%
#     mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>%
#     left_join(parameter_categories, by = 'species') %>%
#     mutate(label = fct_relevel(label, c("Observed Concentration in Total PM2.5",
#                                         "Attributable concentration to smoke PM2.5"))) %>%
#     filter(!is.na(species)) %>% 
#     left_join(coeffs %>% rename(year_trend = 'coeff'), by = 'species')
#   
#   datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")
# 
#   limits = c(
#     # "Organics"
#     "Organic Carbon (OC)", "Elemental Carbon (EC)",
#     # "Halogens"
#     "Bromine (Br)", "Chlorine (Cl)", 
#     #  "Nonmetals"
#     "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
#     # "Other metals"
#     "Aluminum (Al)", "Lead (Pb)",
#     # "Metalloids"
#     "Silicon (Si)", "Arsenic (As)",
#     # "Transition metals"
#     "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
#     # "Alkali metals"
#     "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
#     # "Alkaline-earth metals"
#     "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")
# 
# 
#   # PLOTTING
#   avg_mon_pred_spec_plot <- ggplot(plot_df) +
#     geom_density(aes(x = mon_yr,
#                      y = conc,
#                      fill = label,
#                      color = label),
#                  stat = "identity", alpha = 0.8) +
#     labs(title = expression("Chemical species concentration attributable to wildfire smoke PM2.5"),
#          y = expression("Concentration (\u00b5g/m"^3*")"),
#          x = "") +
#     scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#     geom_hline(aes(yintercept = 0)) +
#     theme_minimal() +
#     scale_color_manual(values = pm_pal) +
#     scale_fill_manual(values = pm_pal) +
#     facet_wrap(~factor(species_long, levels = rev(limits)), scales = 'free_y', ncol = 4) +
#     theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
#           axis.text = element_text(),
#           axis.title = element_text(),
#           legend.text = element_text(),
#           legend.title = element_blank(),
#           axis.ticks.x = element_line(linetype = 1),
#           axis.ticks.length=unit(-0.2, "cm"),
#           panel.grid = element_blank(),
#           legend.position = "bottom",
#           axis.line.y = element_line(color = "grey10"))
#   avg_mon_pred_spec_plot
# 
#   # save file
#   ggsave(
#     filename = 'SIFig6_ALLsmoke_attributable_fraction_trend_regional_model.png',
#     plot = avg_mon_pred_spec_plot,
#     path = file.path(results_fp, 'SI Figs'),
#     scale = 1,
#     width = 12,
#     height = 8,
#     dpi = 320)
#   
#   # save file
#   ggsave(
#     filename = 'SIFig6_ALLsmoke_attributable_fraction_trend_regional_model.pdf',
#     plot = avg_mon_pred_spec_plot,
#     path = file.path(results_fp, 'SI Figs'),
#     scale = 1,
#     width = 12,
#     height = 8,
#     dpi = 320)
#   
#   
#   future::plan(NULL)
#   
#   return(plot_df)#
# }
# 
#   
# 
