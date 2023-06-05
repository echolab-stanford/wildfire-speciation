# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species

# loadd(clean_pm_spec_df, cache = drake_cache)
# species_list <-c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")

thresholds <- read_xlsx(file.path(wip_gdrive_fp, 'species_toxicity_RSL.xlsx')) %>% 
  dplyr::select(species, IUR_ug_m3) %>% # can change which threshold we are working with RfC v IUR
  filter(!is.na(IUR_ug_m3)) 
################################################################################
# How much of total exposure to a given species, 
# or total exposure above some concentration threshold, 
# is driven by wildfire smoke vs other PM sources?
################################################################################
selected_spec_df <- clean_pm_spec_df %>% 
  dplyr::select(Dataset:monitor_month, 
                smokePM, nonsmokePM, totPM2.5,
                AL, AS, BR, CA, EC, OC, CHL, CL, CR, CU,
                FE, K, MG, MN, `NA`, NI, NO3, N2, P, PB, RB,
                S, SE, SI, SOIL, SO4, SR, V, ZN, TI, ZR) 

# get coeffs
base_reg = feols(c(AL, AS, BR, CA, EC, OC, CHL, CL, CR, CU,
                   FE, K, MG, MN, `NA`, NI, NO3, N2, P, PB, RB,
                   S, SE, SI, SOIL, SO4, SR, V, ZN, TI, ZR) ~ smokePM + nonsmokePM | 
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

species_list <- c('AL', 'AS', 'BR', 'CA', 'EC', 'OC', 'CHL', 'CL', 'CR', 'CU',
                  'FE', 'K', 'MG', 'MN', 'NA', 'NI', 'NO3', 'N2', 'P', 'PB', 'RB',
                  'S', 'SE', 'SI', 'SO4', 'SR', 'V', 'ZN', 'TI', 'ZR')
# test
# current_species <- species_list[20]
# map over each species
daily_conc_above_threshold <- purrr::map_df(species_list, function(current_species) {
  
  print(current_species)
  
  # step one: 
  current_weighted_conc_df <- selected_spec_df %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    dplyr::select(-c(State, season, doy, monitor_month)) %>% 
    dplyr::select(Dataset, SiteCode, Date, year, month, region, 
                  smoke_day, totPM2.5, smokePM, nonsmokePM, !!sym(current_species)) %>% 
    filter(!is.na(totPM2.5)) %>% 
    left_join(betas %>% 
                filter(species == current_species), 
              by = 'smoke_day') %>% 
    filter(!!sym(current_species) >= 0) %>% 
    mutate(pred_nonsmoke_spec_conc = Estimate* nonsmokePM,
           pred_smoke_spec_conc = Estimate*smokePM) %>% 
    mutate(totPM_pred_spec_conc = pred_nonsmoke_spec_conc + pred_smoke_spec_conc) %>% 
    distinct() 
  
  # get current human health threshold
    current_threshold <- thresholds %>% 
      filter(species == current_species) %>% 
      distinct(IUR_ug_m3) %>% 
      as.numeric()
    
  # create a flag to determine if a daily concentration is above or below this threshold
    flagged_current_spec_df <- current_weighted_conc_df %>% 
      mutate(s_threshold_flag = case_when(
        pred_smoke_spec_conc >= current_threshold ~ 'HH risk from smoke',
        pred_smoke_spec_conc < current_threshold ~ 'no immediate risk from smoke',
        TRUE ~ NA
      )) %>% 
      mutate(ns_threshold_flag = case_when(
        pred_nonsmoke_spec_conc >= current_threshold ~ 'HH risk from nonsmoke',
        pred_nonsmoke_spec_conc < current_threshold ~ 'no immediate risk from nonsmoke',
        TRUE ~ NA
      )) %>% 
      mutate(totPM_threshold_flag = case_when(
        totPM_pred_spec_conc >= current_threshold ~ 'HH risk from totPM2.5',
        totPM_pred_spec_conc < current_threshold ~ 'no immediate risk from totPM2.5',
        TRUE ~ NA
      )) 
  
   summed_days_above <- tibble(
     species = current_species,
     threshold = current_threshold,
     total_days = nrow(flagged_current_spec_df),
     nonsmoke_days_above = sum(flagged_current_spec_df$ns_threshold_flag == 'HH risk from nonsmoke'),
     smoke_days_above = sum(flagged_current_spec_df$s_threshold_flag == 'HH risk from smoke'),
     totPM2.5_days_above = sum(flagged_current_spec_df$totPM_threshold_flag == 'HH risk from totPM2.5'),
   )
  
  
}) %>% 
  bind_rows() 




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

