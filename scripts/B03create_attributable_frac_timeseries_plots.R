# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species

# loadd(clean_pm_spec_df, cache = drake_cache)
species_list <-c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")

################################################################################
# How much of total exposure to a given species, 
# or total exposure above some concentration threshold, 
# is driven by wildfire smoke vs other PM sources?
################################################################################
selected_spec_df <- clean_pm_spec_df %>% 
  dplyr::select(Dataset:month, smoke_day, V, PB, FE, 
                CU, CR, MN, ZN, NI,
                EC, OC, S, smokePM, 
                nonsmokePM, totPM2.5, monitor_month) %>% 
  # only keep rows that have all values so that its the same selection of sites for the analysis
  na.omit() # %>% 
  # group_by(SiteCode, year) %>% 
  # mutate(num_obs_per_yr_site = n_distinct(Date)) %>% 
  # ungroup()

# run regression
base_reg = feols(c(V, PB, FE, 
                   CU, CR, MN, ZN, NI,
                   EC, OC, S) ~ smokePM + nonsmokePM | 
                   monitor_month + year, selected_spec_df) # add in region,run a model for each region

base_reg = feols(S ~ smokePM + nonsmokePM | monitor_month + year, selected_spec_df)

# prediction
selected_spec_df <- selected_spec_df %>% 
  mutate(predS = predict(base_reg))

# predict using full model and then plot!
# new df
smoke_species_df <- selected_spec_df %>% 
  mutate(nonsmokePM = 0) 

selected_spec_df <- selected_spec_df %>% 
  mutate(smokeS = predict(base_reg, newdata = smoke_species_df))

 

# Save coefficients to table
pred_base_coeffs <- coeftable(base_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  rename(species = 'lhs') %>% 
  dplyr::select(-id) 


# get betas to then weight smoke and nonsmoke concentrations
betas <- pred_base_coeffs %>% 
  mutate(smoke_day = case_when(
    coefficient == 'smokePM' ~ 1,
    coefficient == 'nonsmokePM' ~ 0
  )) %>% 
  dplyr::select(smoke_day, Estimate, species)


# test
# current_species <- species_list[1]
# map over each species
all_spec_attrib_frac_df <- purrr::map_df(species_list, function(current_species) {
  
  # step one: 
  current_weighted_conc_df <- selected_spec_df %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    dplyr::select(-c(State, season, monitor_month)) %>% 
    dplyr::select(Dataset, SiteCode, Date, year, month, region, 
                  smoke_day, smokePM, nonsmokePM, !!sym(current_species)) %>% 
    left_join(betas %>% 
                filter(species == current_species), 
              by = 'smoke_day') %>% 
    filter(!!sym(current_species) >= 0) %>% 
    mutate(pred_nonsmoke_spec_conc = Estimate* nonsmokePM,
           pred_smoke_spec_conc = Estimate*smokePM) %>% 
    mutate(pred_allPM_spec_conc = pred_nonsmoke_spec_conc + pred_smoke_spec_conc) %>% 
    mutate(smoke_attrib_frac = pred_smoke_spec_conc/pred_allPM_spec_conc) %>% 
    mutate(nonsmoke_attrib_frac = pred_nonsmoke_spec_conc/pred_allPM_spec_conc) %>%
    distinct() %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) 
  
  # get avg attributable fraction of species' concentration to total concentration
  avg_spec <- current_weighted_conc_df %>%
    # get avg concentration across all sites for a given month-year
    group_by(year, month, SiteCode) %>% 
    dplyr::summarise(avg_smoke_attrib_lvl = mean(pred_smoke_spec_conc, na.rm = TRUE),
                     avg_nonsmoke_attrib_lvl = mean(pred_nonsmoke_spec_conc, na.rm = TRUE),
                     avg_tot_conc = avg_smoke_attrib_lvl + avg_nonsmoke_attrib_lvl
                     ) %>% 
    ungroup() %>% 
    # now avg the site averages across full sample
      group_by(year, month) %>% 
      dplyr::summarise(avg_smoke_attrib_lvl = mean(avg_smoke_attrib_lvl, na.rm = TRUE),
                       avg_nonsmoke_attrib_lvl = mean(avg_nonsmoke_attrib_lvl, na.rm = TRUE),
                       avg_tot_conc = avg_smoke_attrib_lvl + avg_nonsmoke_attrib_lvl
      ) %>% 
      ungroup() %>% 
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
    mutate(species = current_species) %>% 
    distinct(year, month, mon_yr, species, 
             avg_smoke_attrib_lvl, avg_nonsmoke_attrib_lvl, avg_tot_conc) 

  
}) %>% 
  bind_rows() %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_f = factor(species, levels = c('CR', "V", "PB", 'NI', "CU", "MN", "ZN","FE", 'S', "EC", "OC")))


# adjust dataframe so that it can plot stacked area plot
ns_spec_long_df <- all_spec_attrib_frac_df %>% 
  pivot_longer(cols = c(avg_nonsmoke_attrib_lvl, avg_tot_conc), 
               names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
  mutate(pm_type = factor(pm_type, levels = c("avg_tot_conc",
                                              "avg_nonsmoke_attrib_lvl"))) 
  

# MAKE AN AREA PLOT OF AVG PM OVER TIME ------------------------------------------
ns_stacked_species_area_plot <- ggplot(ns_spec_long_df,
                                    aes(x = mon_yr,
                                        y = avg_spec_conc,
                                        fill = pm_type
                                        )) + 
  geom_area(alpha =.6) + # position = 'fill'
  scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c("grey80", "navy"), 
                    labels = c('Total PM2.5 species concentration (ug/m3)',
                               'Nonsmoke attributable species concentration (ug/m3)')) +
  facet_wrap(~spec_f, scales = 'free', ncol = 3) +
  labs(x = 'Date',
       y = paste('Species Concentration (ug/m3)'),
       title = 'Contribution of PM2.5 type to species concentration over time') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30", size = .5),
        plot.title = element_text(size=12, face='bold'),
        axis.text.x = element_text(face="plain", size = 10),
        axis.text.y = element_text(face="plain", size = 10),
        legend.position = c(0.85, 0.1)) 
ns_stacked_species_area_plot

# save file
ggsave(
  filename = 'Fig4_nonsmoke_conc_attributable_fracTS.png',
  plot = ns_stacked_species_area_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 12,
  height = 6,
  dpi = 320) 



smoke_spec_long_df <- all_spec_attrib_frac_df %>% 
  pivot_longer(cols = c(avg_smoke_attrib_lvl, avg_tot_conc), 
               names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
  mutate(pm_type = factor(pm_type, levels = c("avg_tot_conc",
                                              "avg_smoke_attrib_lvl"))) 


# MAKE AN AREA PLOT OF AVG PM OVER TIME ------------------------------------------
s_stacked_species_area_plot <- ggplot(smoke_spec_long_df,
                                    aes(x = mon_yr,
                                        y = avg_spec_conc,
                                        fill = pm_type
                                    )) + 
  geom_area(alpha =.6) + # position = 'fill'
  scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c("grey80", "firebrick3"), 
                    labels = c('Total PM2.5 species concentration (ug/m3)',
                               'Smoke attributable species concentration (ug/m3)')) +
  facet_wrap(~spec_f, scales = 'free', ncol = 3) +
  labs(x = 'Date',
       y = paste('Species Concentration (ug/m3)'),
       title = 'Contribution of PM2.5 type to species concentration over time') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30", size = .5),
        plot.title = element_text(size=12, face='bold'),
        axis.text.x = element_text(face="plain", size = 10),
        axis.text.y = element_text(face="plain", size = 10),
        legend.position = c(0.85, 0.1)) 
s_stacked_species_area_plot

# plot <- ggarrange(s_stacked_species_area_plot, ns_stacked_species_area_plot)


# save file
ggsave(
  filename = 'Fig4_smoke_conc_attributable_fracTS.png',
  plot = s_stacked_species_area_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 12,
  height = 6,
  dpi = 320)   


# NOW DO THIS WITH SMOKE AND TOTAL PM
# adjust dataframe so that it can plot stacked area plot
all_spec_long_df <- all_spec_attrib_frac_df %>% 
  pivot_longer(cols = c(avg_nonsmoke_attrib_lvl, avg_smoke_attrib_lvl), 
               names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
  mutate(pm_type = factor(pm_type, levels = c("avg_nonsmoke_attrib_lvl",
                                              "avg_smoke_attrib_lvl"))) 


# MAKE AN AREA PLOT OF AVG PM OVER TIME
sns_stacked_species_area_plot <- ggplot(all_spec_long_df,
                                    aes(x = mon_yr,
                                        y = avg_spec_conc,
                                        fill = pm_type
                                    )) + 
  geom_area(alpha =.6) + # position = 'fill'
  scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c("navy", "firebrick3"), 
                    labels = c('Nonsmoke attributable species concentration (ug/m3)',
                               'Smoke attributable species concentration (ug/m3)')) +
  facet_wrap(~spec_f, scales = 'free', ncol = 3) +
  labs(x = 'Date',
       y = paste('Species Concentration (ug/m3)'),
       title = 'Trend in Smoke vs Nonsmoke PM2.5 contribution to species concentration') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30", size = .5),
        plot.title = element_text(size=12, face='bold'),
        axis.text.x = element_text(face="plain", size = 10),
        axis.text.y = element_text(face="plain", size = 10),
        legend.position = c(0.85, 0.1)) 
sns_stacked_species_area_plot


# save file
ggsave(
  filename = 'Fig4alt_both_conc_attributable_fracTS.png',
  plot = sns_stacked_species_area_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 12,
  height = 6,
  dpi = 320)   

# ---------------------------------------------------------------------------
# CHECK STATES IN EACH REGION
# ---------------------------------------------------------------------------
# get betas to then weight smoke and nonsmoke concentrations
betas <- pred_base_coeffs %>% 
  mutate(smoke_day = case_when(
    coefficient == 'smokePM' ~ 1,
    coefficient == 'nonsmokePM' ~ 0
  )) %>% 
  dplyr::select(smoke_day, Estimate, species)


# test
# current_species <- species_list[1]
# map over each species
all_spec_attrib_fracSTATE_df <- purrr::map_df(species_list, function(current_species) {
  
  # step one: 
  current_weighted_conc_df <- selected_spec_df %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    dplyr::select(-c(season, monitor_month)) %>% 
    dplyr::select(Dataset, SiteCode, State, Date, year, month, region, 
                  smoke_day, smokePM, nonsmokePM, !!sym(current_species)) %>% 
    left_join(betas %>% 
                filter(species == current_species), 
              by = 'smoke_day') %>% 
    filter(!!sym(current_species) >= 0) %>% 
    mutate(pred_nonsmoke_spec_conc = Estimate* nonsmokePM,
           pred_smoke_spec_conc = Estimate*smokePM) %>% 
    mutate(pred_allPM_spec_conc = pred_nonsmoke_spec_conc + pred_smoke_spec_conc) %>% 
    mutate(smoke_attrib_frac = pred_smoke_spec_conc/pred_allPM_spec_conc) %>% 
    mutate(nonsmoke_attrib_frac = pred_nonsmoke_spec_conc/pred_allPM_spec_conc) %>%
    distinct() %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) 
  
  # get avg attributable fraction of species' concentration to total concentration
  avg_spec <- current_weighted_conc_df %>%
    # get avg concentration across all sites for a given month-year within a state
    group_by(year, month, SiteCode, State,) %>% 
    dplyr::summarise(avg_smoke_attrib_lvl = mean(pred_smoke_spec_conc, na.rm = TRUE),
                     avg_nonsmoke_attrib_lvl = mean(pred_nonsmoke_spec_conc, na.rm = TRUE),
                     avg_tot_conc = avg_smoke_attrib_lvl + avg_nonsmoke_attrib_lvl
    ) %>% 
    ungroup() %>% 
    # now avg the site averages across full sample
    group_by(year, month, State) %>% 
    dplyr::summarise(avg_smoke_attrib_lvl = mean(avg_smoke_attrib_lvl, na.rm = TRUE),
                     avg_nonsmoke_attrib_lvl = mean(avg_nonsmoke_attrib_lvl, na.rm = TRUE),
                     avg_tot_conc = avg_smoke_attrib_lvl + avg_nonsmoke_attrib_lvl
    ) %>% 
    ungroup() %>% 
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
    mutate(species = current_species) %>% 
    distinct(year, month, State, mon_yr, species, 
             avg_smoke_attrib_lvl, avg_nonsmoke_attrib_lvl, avg_tot_conc) 
  
  
}) %>% 
  bind_rows() %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_f = factor(species, levels = c('CR', "V", "PB", 'NI', "CU", "MN", "ZN","FE", 'S', "EC", "OC")))


# adjust dataframe so that it can plot stacked area plot
all_state_spec_long_df <- all_spec_attrib_fracSTATE_df %>% 
  pivot_longer(cols = c(avg_nonsmoke_attrib_lvl, avg_smoke_attrib_lvl), 
               names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
  mutate(pm_type = factor(pm_type, levels = c("avg_smoke_attrib_lvl",
                                              "avg_nonsmoke_attrib_lvl"))) %>% 
  filter(species == 'OC') %>% 
  # select states from each region
  filter(State %in% c('CA', 'CO', 'TX', 'NC', 'MN', 'NJ'))

# MAKE AN AREA PLOT OF AVG PM OVER TIME
state_OC_area_plot <- ggplot(all_state_spec_long_df,
                                        aes(x = mon_yr,
                                            y = avg_spec_conc,
                                            fill = pm_type
                                        )) + 
  geom_area(alpha =.6) + # position = 'fill'
  scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c("firebrick3", "navy"),
                    labels = c('Smoke attributable concentration (ug/m3)',
                               'Nonmoke attributable concentration (ug/m3)')) +
  facet_wrap(~State, ncol = 3) +
  labs(x = 'Date',
       y = paste('Organic Carbon Concentration (ug/m3)'),
       title = 'Trend in Smoke vs Nonsmoke PM2.5 contribution to Organic Carbon concentration') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30", size = .5),
        plot.title = element_text(size=14, face='bold'),
        axis.text.x = element_text(face="plain", size = 12),
        axis.text.y = element_text(face="plain", size = 12),
        legend.position = 'bottom') 
state_OC_area_plot

# save file
ggsave(
  filename = 'Fig4alt_state_conc_attributable_fracTS.png',
  plot = state_OC_area_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 


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

