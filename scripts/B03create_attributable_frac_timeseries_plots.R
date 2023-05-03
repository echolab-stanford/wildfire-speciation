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
  dplyr::select(Dataset:monitor_month, V, PB, FE, 
                CU, CR, MN, ZN, NI,
                EC, OC, S, smokePM, 
                nonsmokePM, totPM2.5) 

# get coeffs
base_reg = feols(c(V, PB, FE, 
                   CU, CR, MN, ZN, NI,
                   EC, OC, S) ~ smokePM + nonsmokePM | 
                   monitor_month + year, selected_spec_df) # add in region,run a model for each region

res <- etable(base_reg)

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
    dplyr::select(-c(State, season, doy, monitor_month)) %>% 
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
    distinct() 
  
  # get avg attributable fraction of species' concentration to total concentation
  avg_spec <- current_weighted_conc_df %>% 
    mutate(season = case_when(
      month %in% c(12, 1, 2) ~ 'winter',
      month %in% c(3, 4, 5) ~ 'spring',
      month %in% c(6, 7, 8) ~ 'summer',
      month %in% c(9, 10, 11) ~ 'fall'
    )) %>% 
    # get avg concentration across all sites for a given month-year
    group_by(year, month) %>% 
    dplyr::summarise(avg_smoke_attrib_lvl = mean(pred_smoke_spec_conc, na.rm = TRUE),
                     avg_nonsmoke_attrib_lvl = mean(pred_nonsmoke_spec_conc, na.rm = TRUE)
                     ) %>% 
    ungroup() %>% 
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
    mutate(species = current_species) %>% 
    # group_by(year, season) %>% 
    # dplyr::mutate(avg_ssn_smoke_attrib_frac = mean(smoke_attrib_frac, na.rm = TRUE)) %>% 
    # ungroup() %>% 
    distinct(year, month, mon_yr, species, 
             avg_smoke_attrib_lvl, avg_nonsmoke_attrib_lvl)

  
}) %>% 
  bind_rows() %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_f = factor(species, levels = c("V", "PB", "CR", "FE","CU", "MN", "ZN", "NI", "EC", "OC", "S")))


# adjust dataframe so that it can plot stacked area plot
all_spec_long_df <- all_spec_attrib_frac_df %>% 
  pivot_longer(cols = c(avg_smoke_attrib_lvl, avg_nonsmoke_attrib_lvl), 
               names_to = 'pm_type', values_to = 'avg_spec_conc')

# MAKE THIS IN AREA PLOTS
stacked_species_area_plot <- ggplot(all_spec_long_df,
                                    aes(x = mon_yr,
                                        y = avg_spec_conc,
                                        fill = pm_type
                                        )) + 
  geom_area(alpha =.6) +
scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c("steelblue", "coral"),
                    labels = c('Nonsmoke Species Conc (ug/m3)',
                               'Smoke Species Conc (ug/m3)')) +
  facet_wrap(~spec_f, scales = 'free', ncol = 2) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(x = 'Date',
       y = paste('Species Concentration (ug/m3)'),
       title = 'Trend in smoke vs nonsmoke PM2.5 contribution to species concentration') +
  theme_minimal()
stacked_species_area_plot



# make time series attributable fraction (af) plot
# PLOT ONE AT A TIME:
species_frac_plot <- ggplot(all_spec_attrib_frac_df,
                            aes(color = spec_type,
                                x = mon_yr,
                                y = avg_smoke_attrib_lvl,
                                group = spec_type)) + 
  # geom_line(aes(x = mon_yr, 
  #               y = avg_smoke_attrib_frac*100, # multiply by 100 to get in percent 
  #   color = spec_type), size = 1) +
  geom_line(aes(color = spec_type), size = 1) +
  geom_smooth(method=lm, se=FALSE, col='black', size=.5) +
  scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  # scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',"red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
  facet_wrap(~spec_f, scales = 'free', ncol = 2) +
  labs(x = 'Date',
       y = paste('% of Species Concentration Attributable to Nonsmoke PM2.5'),
       title = 'Trend in wildfire smoke contribution to species concentration in total PM2.5') +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  theme_minimal()
species_frac_plot

# plot one at a time
each_species_frac_plot <- ggplot(all_spec_attrib_frac_df,
                            aes(x = mon_yr,
                                y = avg_mon_smoke_attrib_frac*100, # multiply by 100 to get in percent 
                                color = spec_type, 
                                group = spec_type)) + 
  geom_line(aes(color = spec_type), size = 1) +
  geom_smooth(method=lm, se=FALSE, col='black', size=.5) +
  scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  # scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',"red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
  facet_wrap(~spec_f, ncol = 3, scales = 'free_y') +
  labs(x = 'Date',
       y = paste('% of Species Concentration Attributable to Wildfire Smoke PM2.5'),
       title = 'Trend in wildfire smoke contribution to species concentration in total PM2.5') +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  theme_minimal()
each_species_frac_plot


# SEASONAL PLOT
seasonal_frac_plot <- ggplot(all_spec_attrib_frac_df,
                       aes(x = mon_yr,
                           y = avg_ssn_smoke_attrib_frac*100, # multiply by 100 to get in percent 
                           color = spec_type, 
                           group = spec_type)) + 
  geom_line(aes(color = spec_type), size = 1) +
  geom_smooth(method=lm, se=FALSE, col='black', size=.5) +
  scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
  scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  # scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',"red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
  facet_wrap(~spec_f, scales = 'free', ncol = 4) +
  labs(x = 'Date',
       y = paste('% of Species Concentration Attributable to Wildfire Smoke PM2.5'),
       title = 'Seasonal trend in wildfire smoke contribution to species concentrations in total PM2.5') +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  theme_minimal()
seasonal_frac_plot

