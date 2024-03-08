# Emma Krasovich Southworth, emmars@stanford.edu
# 
# # ---------------------------------------------------------------------------
# # TEMPORAL HETEROGENEITY
# # ---------------------------------------------------------------------------
# temporal_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
#                         ~ smokePM + nonsmokePM | 
#                           monitor_month + year, reg_df, fsplit = ~time_period) 
# etable(temporal_reg) # look at regression results
# 
# # get coefficients and prepare for plotting
# temporalPM_coeffs <- coeftable(temporal_reg) %>% 
#   rename(pval = 'Pr(>|t|)',
#          se = 'Std. Error',
#          species = 'lhs',
#          pm_type = 'coefficient',
#          time_per = 'sample') %>% 
#   # get pvalues
#   mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
#   dplyr::select(-id, -sample.var)  %>% 
#   mutate(spec_type = case_when(
#     species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
#     species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
#     species=="EC" | species=="OC" ~ "Secondary organic",
#     species=="S" ~ "Toxicity potentiator",
#     TRUE ~ NA)) %>% 
#   mutate(spec_type = factor(spec_type, 
#                             levels = c("Toxic metal", "Transition metal",
#                                        "Secondary organic",
#                                        "Toxicity potentiator"))) %>% 
#   mutate(spec_f = factor(species, 
#                          levels = c("V", "CR", "PB", 
#                                     "CU", "MN", "NI", "ZN","FE", 
#                                     "S", 
#                                     "EC", "OC"))) 
# 
# # marginal plot
# temporal_reg_plot <- ggplot(temporalPM_coeffs, 
#                             aes(x = time_per,
#                                 y = Estimate, 
#                                 color=spec_type, shape = pm_type)) +
#   geom_point(size=2, alpha = 0.7, stat = "identity") + # shape = pm_type
#   geom_linerange(aes(ymin = (Estimate - se), 
#                      ymax = (Estimate + se)), stat = "identity") +
#   scale_x_discrete(limits = c("Full sample", "2006-2010","2010-2015", "2015-2020")) +
#   scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
#   scale_shape_manual(values=c(15, 17,18, 19)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#   facet_wrap(~spec_f, ncol = 3, scales = 'free_y') +
#   labs(y = 'Concentration (ug/m3)',
#        x = 'Time Period in Sample',
#        color = 'Species Category',
#        shape = 'PM Type',
#        title = "Comparing chemical composition of nonsmoke and smoke PM over time",
#        subtitle = 'nonsmoke + smoke on RHS, monitor-month and year FE') +
#   theme_light() 
# temporal_reg_plot
# 
# 
# 
# 
# # ---------------------------------------------------------------------------
# # 3.CSN VS IMPROVE INTERACTION
# # ---------------------------------------------------------------------------
# monitor_regional_PM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
#                         ~ smokePM + nonsmokePM + Dataset + smokePM*Dataset + nonsmokePM*Dataset | 
#                           monitor_month + year, reg_df, fsplit = ~region) 
# etable(monitor_regional_PM_reg) # look at regression results
# 
# # get coefficients and prepare for plotting
# monitor_regionalPM_coeffs <- coeftable(monitor_regional_PM_reg) %>% 
#   rename(pval = 'Pr(>|t|)',
#          se = 'Std. Error',
#          species = 'lhs',
#          pm_type = 'coefficient',
#          region = 'sample') %>% 
#   # get pvalues
#   mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
#   dplyr::select(-id, -sample.var)  %>% 
#   mutate(spec_type = case_when(
#     species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
#     species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
#     species=="EC" | species=="OC" ~ "Secondary organic",
#     species=="S" ~ "Toxicity potentiator",
#     TRUE ~ NA)) %>% 
#   mutate(spec_type = factor(spec_type, 
#                             levels = c("Toxic metal", "Transition metal",
#                                        "Secondary organic",
#                                        "Toxicity potentiator"))) 
# 
# 
# # merge sample avg for each species and divide each species' betas by full sample avg for each species
# regionalPMcoeffs_normalized <- regionalPM_coeffs %>% 
#   left_join(spec_regional_avgs_df, by = c('region', 'species')) %>% 
#   mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
#          norm_se = se/baseline_avg) 
# 
# 
# # plot coefficients for speciation at the avg monitor which tells us how much 
# # of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
# regional_reg_plot <- ggplot(regionalPMcoeffs_normalized, 
#                             aes(x = pm_type,
#                                 y = 100*norm_est, 
#                                 color=region)) +
#   geom_point(size=2, alpha = 0.6, stat = "identity") + # shape = pm_type
#   geom_errorbar(aes(ymin = (norm_est - norm_se)*100, 
#                     ymax = (norm_est + norm_se)*100, 
#                     width =.2), stat = "identity") +
#   #scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
#   #scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
#   scale_shape_manual(values=c(15, 17,18, 19)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#   facet_wrap(~species, ncol = 6, scales = 'free_y') +
#   labs(y = '% Change (ug/m3) relative to regional baseline',
#        x = 'Chemical Species (ug/m3)',
#        color = 'Species Category',
#        shape = 'PM Type',
#        title = "Percent change in species concentration relative to regional baseline for a 1 ug/m3 increase in smoke and nonsmoke",
#        subtitle = 'nonsmoke + smoke on RHS, monitor-month and year FE') +
#   theme_light() 
# regional_reg_plot
# 
# 
# 
# 
# 
# 
# 
# # csn_only <- reg_df %>% 
# #   filter(Dataset == 'EPACSN')
# # 
# # csn_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN)
# #                      ~ smokePM + nonsmokePM | 
# #                        monitor_month + year, csn_only)
# # 
# # # IMPROVE SITES ONLY
# # improve_only <- reg_df %>% 
# #   filter(Dataset == 'IMPROVE')
# # 
# # improve_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN)
# #                 ~ smokePM + nonsmokePM | 
# #                   monitor_month + year, improve_only)
# #                                
# # # get coeffs for both improve and csn
# # site_coeffs <- coeftable(improve_reg) %>%
# #   mutate(site = 'IMPROVE') %>% 
# #   bind_rows(coeftable(csn_reg) %>% 
# #               mutate(site = 'EPACSN')) %>% 
# #   rename(pval = 'Pr(>|t|)',
# #          se = 'Std. Error',
# #          species = 'lhs') %>% 
# #   # get pvalues
# #   mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
# #   mutate(spec_type = case_when(
# #     species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
# #     species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
# #     species=="EC" | species=="OC" ~ "Secondary organic",
# #     species=="S" ~ "Toxicity potentiator",
# #     TRUE ~ NA)) %>% 
# #   mutate(spec_type = factor(spec_type, 
# #                             levels = c("Toxic metal", "Transition metal",
# #                                        "Secondary organic",
# #                                        "Toxicity potentiator"))) 
# #  
# # # plot the comparison between CSN vs Improve
# # site_reg_plot <- ggplot(site_coeffs, 
# #                         aes(x = species,
# #                             y = Estimate, 
# #                             color=spec_type)) +
# #   geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
# #   scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
# #   scale_y_continuous(trans = log10_trans(),
# #                      breaks = trans_breaks("log10", function(x) 10^x),
# #                      labels = scales::percent_format(accuracy = .01)) +
# #   scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
# #   scale_shape_manual(values=c(15, 17,18, 19)) +
# #   geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
# #   facet_wrap(~site)+
# #   labs(x = 'Chemical Species (ug/m3)',
# #        y = 'Estimate',
# #        color = 'Site Type',
# #        shape = 'PM Type',
# #        title = "Comparing CSN vs Improve Speciation") +
# #   theme_light() 
# # site_reg_plot
# # 
# # # calculate the sample averages for each chemical species:
# # csn_baseline_avgs_df <- tibble(
# #   species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
# #   baseline_avg = c(mean(csn_only$V, na.rm = TRUE), 
# #                    mean(csn_only$PB, na.rm = TRUE), 
# #                    mean(csn_only$FE, na.rm = TRUE), 
# #                    mean(csn_only$CU, na.rm = TRUE),  
# #                    mean(csn_only$CR, na.rm = TRUE),  
# #                    mean(csn_only$MN, na.rm = TRUE), 
# #                    mean(csn_only$ZN, na.rm = TRUE),  
# #                    mean(csn_only$NI, na.rm = TRUE),  
# #                    mean(csn_only$EC, na.rm = TRUE),  
# #                    mean(csn_only$OC, na.rm = TRUE),  
# #                    mean(csn_only$S, na.rm = TRUE)),
# #   site ='EPACSN'
# # ) %>% 
# #   bind_rows(tibble(
# #     species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
# #     baseline_avg = c(mean(improve_only$V, na.rm = TRUE), 
# #                      mean(improve_only$PB, na.rm = TRUE), 
# #                      mean(improve_only$FE, na.rm = TRUE), 
# #                      mean(improve_only$CU, na.rm = TRUE),  
# #                      mean(improve_only$CR, na.rm = TRUE),  
# #                      mean(improve_only$MN, na.rm = TRUE), 
# #                      mean(improve_only$ZN, na.rm = TRUE),  
# #                      mean(improve_only$NI, na.rm = TRUE),  
# #                      mean(improve_only$EC, na.rm = TRUE),  
# #                      mean(improve_only$OC, na.rm = TRUE),  
# #                      mean(improve_only$S, na.rm = TRUE)),
# #     site ='IMRPOVE'
# #   ))
#   
# # merge with coeffs:
# norm_site_coeffs <- site_coeffs %>% 
#   left_join(csn_baseline_avgs_df, by = c('site', 'species')) %>% 
#   mutate(norm_coeff = Estimate/baseline_avg)
# 
# # PLOT NORMALIZED CSN VS IMPROVE
# # plot the comparison between CSN vs Improve
# site_reg_plot <- ggplot(norm_site_coeffs, 
#                         aes(x = species,
#                             y = Estimate, 
#                             color=spec_type)) +
#   geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
#   scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
#   scale_y_continuous(trans = log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = scales::percent_format(accuracy = .01)) +
#   scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
#   scale_shape_manual(values=c(15, 17,18, 19)) +
#   geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
#   facet_wrap(~site)+
#   labs(x = 'Chemical Species (ug/m3)',
#        y = 'Estimate',
#        color = 'Site Type',
#        shape = 'PM Type',
#        title = "Comparing CSN vs Improve Speciation, normalized") +
#   theme_light() 
# site_reg_plot
# 
# 
# # ---------------------------------------------------------------------------
# # RUN A MODEL WHERE EACH PM TYPE IS PLOTTED ON SAME
# # ---------------------------------------------------------------------------
# smokePM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ smokePM | monitor_month + year, reg_df)
# nonsmokePM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ nonsmokePM | monitor_month + year, reg_df)
# totPM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ totPM2.5 | monitor_month + year, reg_df)
# 
# # COMBINED COEFFS FOR PLOTTING
# pm_type_coeffs <- coeftable(smokePM_reg) %>%
#   mutate(pm_type = 'smoke PM2.5') %>% 
#   bind_rows(coeftable(nonsmokePM_reg) %>%
#               mutate(pm_type = 'nonsmoke PM2.5')) %>% 
#   bind_rows(coeftable(totPM_reg) %>%
#               mutate(pm_type = 'total PM2.5')) %>% 
#   rename(pval = 'Pr(>|t|)',
#          se = 'Std. Error',
#          species = 'lhs') %>% 
#   # get pvalues
#   mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
#   mutate(spec_type = case_when(
#     species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
#     species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
#     species=="EC" | species=="OC" ~ "Secondary organic",
#     species=="S" ~ "Toxicity potentiator",
#     TRUE ~ NA)) %>% 
#   mutate(spec_type = factor(spec_type, 
#                             levels = c("Toxic metal", "Transition metal",
#                                        "Secondary organic",
#                                        "Toxicity potentiator"))) 
# 
# # join with species baseline avgs + divide to normalize betas
# norm_pm_type_coeffs <- pm_type_coeffs %>% 
#   left_join(spec_baseline_avgs_df, by = 'species') %>%
#   mutate(norm_est = Estimate/baseline_avg)
# 
# # plot the normalized values for each smoke/nonsmoke/total pm
# # plot the comparison between CSN vs Improve
# pmtype_reg_plot <- ggplot(norm_pm_type_coeffs, 
#                         aes(x = species,
#                             y = Estimate, 
#                             color=spec_type)) +
#   geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
#   scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
#   scale_y_continuous(trans = log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = scales::percent_format(accuracy = .01)) +
#   scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
#   scale_shape_manual(values=c(15, 17,18, 19)) +
#   geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
#   labs(x = 'Chemical Species (ug/m3)',
#        y = 'Estimate',
#        color = 'Species Category',
#        shape = 'PM Type',
#        title = "Comparing PM Types, normalized") +
#   theme_light() 
# pmtype_reg_plot
# 
# 
# #  -----------------------------------------------------------------------------
# # DO SOME REGRESSION CHECKS 
# # ------------------------------------------------------------------------------
# # CHECK 1) ARE BETAS THE SAME FOR WHEN WE RUN TOTAL PM ON LHS WITH AND W/O SMOKE
# nonsmoke_control_reg = feols(conc_val_totPM2.5 ~ conc_val_smokePM + conc_val_nonsmokePM | monitor_month + year, conc_rel2avg_df) 
# summary(nonsmoke_control_reg)
# no_control_reg = feols(conc_val_totPM2.5 ~ conc_val_smokePM | monitor_month + year, conc_rel2avg_df) 
# summary(no_control_reg)
# nonosmoke_control_reg = feols(conc_val_totPM2.5 ~ conc_val_nonsmokePM | monitor_month + year, conc_rel2avg_df) 
# summary(nonosmoke_control_reg)
# 
# 
# 




# make time series attributable fraction (af) plot
# PLOT ONE AT A TIME:
# species_frac_plot <- ggplot(all_spec_attrib_frac_df,
#                             aes(group = spec_f)) + 
#   geom_line(aes(x = mon_yr,
#                 y = avg_smoke_attrib_lvl*100, # multiply by 100 to get in percent
#                 color = 'Smoke PM'), size = .75) +
#   geom_line(aes(x = mon_yr,
#                 y = avg_nonsmoke_attrib_lvl*100, # multiply by 100 to get in percent
#                 color = 'Nonsmoke PM'), size = .75) +
#   geom_line(aes(x = mon_yr,
#                 y = avg_tot_conc*100, # multiply by 100 to get in percent
#                 color = 'Total PM2.5'), size = .5, alpha =.7) +
#   scale_color_manual(name = "PM Type", values = c("Smoke PM" = "coral", "Nonsmoke PM" = "steelblue", "Total PM2.5" = 'black')) +
#   scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
#   # scale_color_manual(name = "PM Type", values = c("Y1" = "darkblue", "Y2" = "red", 'Y3' = 'black')) +
#   # scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',"red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
#   facet_wrap(~spec_f, scales = 'free', ncol = 2) +
#   labs(x = 'Date',
#        y = paste('% of Species Concentration Attributable to Each PM2.5 Type'),
#        color = 'PM Type',
#        title = 'Trend in wildfire smoke contribution to species concentration in total PM2.5') +
#   theme(axis.text.x = element_text(angle=45, hjust = 1)) +
#   theme_minimal()
# species_frac_plot
# 
# # plot one at a time
# each_species_frac_plot <- ggplot(all_spec_attrib_frac_df,
#                             aes(x = mon_yr,
#                                 y = avg_mon_smoke_attrib_frac*100, # multiply by 100 to get in percent 
#                                 color = spec_type, 
#                                 group = spec_type)) + 
#   geom_line(aes(color = spec_type), size = 1) +
#   geom_smooth(method=lm, se=FALSE, col='black', size=.5) +
#   scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
#   scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
#   # scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',"red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
#   facet_wrap(~spec_f, ncol = 3, scales = 'free_y') +
#   labs(x = 'Date',
#        y = paste('% of Species Concentration Attributable to Wildfire Smoke PM2.5'),
#        title = 'Trend in wildfire smoke contribution to species concentration in total PM2.5') +
#   theme(axis.text.x = element_text(angle=45, hjust = 1)) +
#   theme_minimal()
# each_species_frac_plot
# 
# 
# # SEASONAL PLOT
# seasonal_frac_plot <- ggplot(all_spec_attrib_frac_df,
#                        aes(x = mon_yr,
#                            y = avg_ssn_smoke_attrib_frac*100, # multiply by 100 to get in percent 
#                            color = spec_type, 
#                            group = spec_type)) + 
#   geom_line(aes(color = spec_type), size = 1) +
#   geom_smooth(method=lm, se=FALSE, col='black', size=.5) +
#   scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
#   scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
#   # scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',"red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
#   facet_wrap(~spec_f, scales = 'free', ncol = 4) +
#   labs(x = 'Date',
#        y = paste('% of Species Concentration Attributable to Wildfire Smoke PM2.5'),
#        title = 'Seasonal trend in wildfire smoke contribution to species concentrations in total PM2.5') +
#   theme(axis.text.x = element_text(angle=45, hjust = 1)) +
#   theme_minimal()
# seasonal_frac_plot


# ---------------------------------------------------------------------------
# CHECK STATES IN EACH REGION
# ---------------------------------------------------------------------------
# get betas to then weight smoke and nonsmoke concentrations
# betas <- pred_base_coeffs %>% 
#   mutate(smoke_day = case_when(
#     coefficient == 'smokePM' ~ 1,
#     coefficient == 'nonsmokePM' ~ 0
#   )) %>% 
#   dplyr::select(smoke_day, Estimate, species)


# test
# current_species <- species_list[1]
# map over each species
# all_spec_attrib_fracSTATE_df <- purrr::map_df(species_list, function(current_species) {
#   
#   # step one: 
#   current_weighted_conc_df <- selected_spec_df %>% 
#     mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
#     dplyr::select(-c(season, monitor_month)) %>% 
#     dplyr::select(Dataset, SiteCode, State, Date, year, month, region, 
#                   smoke_day, smokePM, nonsmokePM, !!sym(current_species)) %>% 
#     left_join(betas %>% 
#                 filter(species == current_species), 
#               by = 'smoke_day') %>% 
#     filter(!!sym(current_species) >= 0) %>% 
#     mutate(pred_nonsmoke_spec_conc = Estimate* nonsmokePM,
#            pred_smoke_spec_conc = Estimate*smokePM) %>% 
#     mutate(pred_allPM_spec_conc = pred_nonsmoke_spec_conc + pred_smoke_spec_conc) %>% 
#     mutate(smoke_attrib_frac = pred_smoke_spec_conc/pred_allPM_spec_conc) %>% 
#     mutate(nonsmoke_attrib_frac = pred_nonsmoke_spec_conc/pred_allPM_spec_conc) %>%
#     distinct() %>% 
#     mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) 
#   
#   # get avg attributable fraction of species' concentration to total concentration
#   avg_spec <- current_weighted_conc_df %>%
#     # get avg concentration across all sites for a given month-year within a state
#     group_by(year, month, SiteCode, State,) %>% 
#     dplyr::summarise(avg_smoke_attrib_lvl = mean(pred_smoke_spec_conc, na.rm = TRUE),
#                      avg_nonsmoke_attrib_lvl = mean(pred_nonsmoke_spec_conc, na.rm = TRUE),
#                      avg_tot_conc = avg_smoke_attrib_lvl + avg_nonsmoke_attrib_lvl
#     ) %>% 
#     ungroup() %>% 
#     # now avg the site averages across full sample
#     group_by(year, month, State) %>% 
#     dplyr::summarise(avg_smoke_attrib_lvl = mean(avg_smoke_attrib_lvl, na.rm = TRUE),
#                      avg_nonsmoke_attrib_lvl = mean(avg_nonsmoke_attrib_lvl, na.rm = TRUE),
#                      avg_tot_conc = avg_smoke_attrib_lvl + avg_nonsmoke_attrib_lvl
#     ) %>% 
#     ungroup() %>% 
#     mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
#     mutate(species = current_species) %>% 
#     distinct(year, month, State, mon_yr, species, 
#              avg_smoke_attrib_lvl, avg_nonsmoke_attrib_lvl, avg_tot_conc) 
#   
#   
# }) %>% 
#   bind_rows() %>% 
#   mutate(spec_type = case_when(
#     species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
#     species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
#     species=="EC" | species=="OC" ~ "Secondary organic",
#     species=="S" ~ "Toxicity potentiator",
#     TRUE ~ NA)) %>% 
#   mutate(spec_f = factor(species, levels = c('CR', "V", "PB", 'NI', "CU", "MN", "ZN","FE", 'S', "EC", "OC")))
# 
# 
# # adjust dataframe so that it can plot stacked area plot
# all_state_spec_long_df <- all_spec_attrib_fracSTATE_df %>% 
#   pivot_longer(cols = c(avg_nonsmoke_attrib_lvl, avg_smoke_attrib_lvl), 
#                names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
#   mutate(pm_type = factor(pm_type, levels = c("avg_smoke_attrib_lvl",
#                                               "avg_nonsmoke_attrib_lvl"))) %>% 
#   filter(species == 'OC') %>% 
#   # select states from each region
#   filter(State %in% c('CA', 'CO', 'TX', 'NC', 'MN', 'NJ'))
# 
# # MAKE AN AREA PLOT OF AVG PM OVER TIME
# state_OC_area_plot <- ggplot(all_state_spec_long_df,
#                              aes(x = mon_yr,
#                                  y = avg_spec_conc,
#                                  fill = pm_type
#                              )) + 
#   geom_area(alpha =.6) + # position = 'fill'
#   scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
#   scale_fill_manual(name = 'PM2.5 Type',
#                     values=c("firebrick3", "navy"),
#                     labels = c('Smoke attributable concentration (ug/m3)',
#                                'Nonmoke attributable concentration (ug/m3)')) +
#   facet_wrap(~State, ncol = 3) +
#   labs(x = 'Date',
#        y = paste('Organic Carbon Concentration (ug/m3)'),
#        title = 'Trend in Smoke vs Nonsmoke PM2.5 contribution to Organic Carbon concentration') +
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "grey30", size = .5),
#         plot.title = element_text(size=14, face='bold'),
#         axis.text.x = element_text(face="plain", size = 12),
#         axis.text.y = element_text(face="plain", size = 12),
#         legend.position = 'bottom') 
# state_OC_area_plot
# 
# # save file
# ggsave(
#   filename = 'Fig4alt_state_conc_attributable_fracTS.png',
#   plot = state_OC_area_plot,
#   path = file.path(wip_gdrive_fp, 'figures/'),
#   scale = 1,
#   width = 10,
#   height = 6,
#   dpi = 320) 

