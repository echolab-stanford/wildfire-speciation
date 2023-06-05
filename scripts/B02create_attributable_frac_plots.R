# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species

# loadd(clean_pm_spec_df, cache = drake_cache)

create_attributable_frac_plots <- function(clean_pm_spec_df)

# set up list of species we care about
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
# current_species <- species_list[3]
# map over each species
all_spec_quantiles_df <- purrr::map_df(species_list, function(current_species) {
  
  # step one: 
  current_weighted_conc_df <- selected_spec_df %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    rename_at(vars(matches("^conc_val_")), ~str_remove(., "^conc_val_")) %>% 
    dplyr::select(-c(starts_with('conc_val_'), State, season, doy, moy, monitor_month)) %>% 
    dplyr::select(Dataset, SiteCode, Date, region, 
                  smoke_day, smokePM, nonsmokePM, !!sym(current_species)) %>% 
    left_join(betas %>% 
                filter(species == current_species), 
              by = 'smoke_day') %>% 
    filter(!!sym(current_species) >= 0) %>% 
    mutate(pred_nonsmoke = Estimate* nonsmokePM,
           pred_smoke = Estimate*smokePM) %>% 
    mutate(pred_allPM = pred_nonsmoke+pred_smoke) %>% 
    distinct() 
  
  
  # get quantiles for summed weighted total PM
  current_quantiles <- quantile(current_weighted_conc_df$pred_allPM, 
                                seq(0, 1, .1), na.rm = TRUE)
  
  # create dataframe where all concentrations of a species have been identified to be above or below
  quant_df <- current_weighted_conc_df %>%
    mutate(q0 = ifelse(pred_allPM >= current_quantiles[1], '0%', NA),
           q10 = ifelse(pred_allPM >= current_quantiles[2], '10%', NA),
           q20 = ifelse(pred_allPM >= current_quantiles[3], '20%', NA),
           q30 = ifelse(pred_allPM >= current_quantiles[4], '30%', NA),
           q40 = ifelse(pred_allPM >= current_quantiles[5], '40%', NA),
           q50 = ifelse(pred_allPM >= current_quantiles[6], '50%', NA),
           q60 = ifelse(pred_allPM >= current_quantiles[7], '60%', NA),
           q70 = ifelse(pred_allPM >= current_quantiles[8], '70%', NA),
           q80 = ifelse(pred_allPM >= current_quantiles[9], '80%', NA),
           q90 = ifelse(pred_allPM >= current_quantiles[10], '90%', NA),
           q100 = ifelse(pred_allPM >= current_quantiles[11], '100%', NA))
           
  
  # get list of quantiles
  q_list <- c('q0', 'q10', 'q20', 'q30', 'q40', 'q50', 'q60', 'q70', 'q80', 'q90', 'q100')
  # current_quant <- q_list[1] # test
  
  # map over the each quantile and create a dataframe for the specific species
  # that should help us plot how much of total exposure to a given species is driven by wildfire smoke
  frac_species_smoke_df <- purrr::map_df(q_list, function(current_quant) {
    
    # create a dateframe for each quantile and bind together:
    current_q_df <- quant_df %>% 
      # filter to the first quantile by converting the string to a symbol so R recognizes it as a column name
      filter(!is.na(!!sym(current_quant))) %>% 
      mutate(q_val = min(pred_allPM)) %>% 
      mutate(tot_pred_nonsmoke = sum(pred_nonsmoke, na.rm=TRUE),
             tot_pred_smoke = sum(pred_smoke, na.rm= TRUE)) %>% 
      mutate(num_days = max(row_number())) %>% 
      mutate(fraction_smoke = 100*(tot_pred_smoke/(tot_pred_smoke + tot_pred_nonsmoke))) %>% 
      distinct(species, num_days, !!sym(current_quant), q_val, fraction_smoke, tot_pred_smoke, tot_pred_nonsmoke) %>% 
      rename(quantile = !!sym(current_quant)) 
    
  }) %>% 
    bind_rows() 
  
  # now create an order to the quantiles for plotting
  frac_species_smoke_df <- frac_species_smoke_df %>% 
    mutate(quantile = factor(quantile, 
                             levels = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"))) %>% 
    mutate(spec_type = case_when(
      species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
      species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
      species=="EC" | species=="OC" ~ "Secondary organic",
      species=="S" ~ "Toxicity potentiator",
      TRUE ~ NA)) 
  
  
  # PLOT ONE AT A TIME:
  # current_species_frac_plot <- ggplot(frac_species_smoke_df,
  #                                     aes(x = quantile,
  #                                         y = fraction_smoke, color = spec_type, group = spec_type)) + # to get in percent
  #   geom_line(size=.5) +
  #   geom_point(size=1.5, shape = 16, alpha = 0.6) +
  #   geom_text(aes(label = num_days), vjust = 0) +
  #   scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',
  #                                      "red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
  #                                        facet_wrap(spec_type ~ species, scales = 'free', ncol = 3) +
  #   ylim(0, 50) +
  #   labs(x = 'Daily Concentration >= quantile',
  #        y = '% of contribution from wildfire smoke',
  #        color = 'Species',
  #        title = "Contribution of wildfire smoke to overall and extreme daily species concentrations") +
  #   theme_minimal()
  # species_frac_plot
  
}) %>% 
  bind_rows() %>%  
  distinct()

classified_pred_spec <- all_spec_quantiles_df %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(q_val = prettyNum(round(q_val, digits = 5)))


# PLOT FIG 3 EXAMPLE:
species_frac_plot <- ggplot(classified_pred_spec %>% 
                              filter(spec_type == 'Secondary organic'),
                            aes(x = quantile,
                                y = fraction_smoke, color = spec_type, group = spec_type)) + # to get in percent
  geom_line(size=.5) +
  geom_point(size=1.5, shape = 16, alpha = 0.6) +
  geom_text(aes(label = prettyNum(num_days, big.mark = ",")), vjust = -2, color = 'black', size = 2) +
  # geom_text(aes(label = prettyNum(q_val)), vjust = -6, color = 'black', size = 2) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod3', 'navy',
                                     "red", "aquamarine", 'violetred', 'yellow', 'darkorange')) +
                                       facet_wrap(spec_type ~ species, scales = 'free', ncol = 3) +
  ylim(0, 50) +
  labs(x = 'Daily Concentration >= quantile',
       y = '% of contribution from wildfire smoke',
       color = 'Species',
       title = "Contribution of wildfire smoke to overall and extreme daily species concentrations") +
  theme_minimal()
species_frac_plot

# -------------------------------------------------------------------------------------------
# CREATE ATTRIBUTABLE FRACTION PLOTS WITH BOTTOM 25, MIDDLE 50, and TOP 25
# -------------------------------------------------------------------------------------------
# current_species <- species_list[1]
summ_stat_quantiles_df <- purrr::map_df(species_list, function(current_species) {
  
  # step one: 
  current_weighted_conc_df <- selected_spec_df %>% 
    mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
    dplyr::select(-c(State, season, doy, monitor_month)) %>% 
    dplyr::select(Dataset, SiteCode, Date, region, totPM2.5,
                  smoke_day, smokePM, nonsmokePM, !!sym(current_species)) %>% 
    left_join(betas %>% 
                filter(species == current_species), 
              by = 'smoke_day') %>% 
    filter(!!sym(current_species) >= 0) %>% 
    mutate(pred_nonsmoke = Estimate* nonsmokePM,
           pred_smoke = Estimate*smokePM) %>% 
    mutate(pred_allPM = pred_nonsmoke+pred_smoke) %>% 
    distinct() %>% 
    filter(!is.na(totPM2.5))
  
  
  # get quantiles for summed weighted total PM
  current_quantiles <- quantile(current_weighted_conc_df$pred_allPM, 
                                seq(0, 1, .25), na.rm = TRUE)
  
  # create dataframe where all concentrations of a species have been identified to be above or below a quantile
  quant_df <- current_weighted_conc_df %>%
    mutate(quant_flag = case_when(
          pred_allPM < current_quantiles[2] ~ 'bottom 25%',
           pred_allPM >= current_quantiles[2] & pred_allPM <= current_quantiles[4] ~ 'middle 50%',
           pred_allPM > current_quantiles[4] ~ 'top 25%')) %>% 
    mutate(q_conc = case_when(
      pred_allPM < current_quantiles[2] ~ paste('below', current_quantiles[2]),
      pred_allPM >= current_quantiles[2] & pred_allPM <= current_quantiles[4] ~ paste0('between ', current_quantiles[2], ' - ', current_quantiles[4]),
      pred_allPM > current_quantiles[4] ~ paste('above', current_quantiles[4]))
    )

  # get list of quantiles
  q_list <- c('bottom 25%', 'middle 50%', 'top 25%')
  # current_quant <- q_list[1] # test
  
  # map over the each quantile and create a dataframe for the specific species
  # that should help us plot how much of total exposure to a given species is driven by wildfire smoke
  frac_species_smoke_df <- purrr::map_df(q_list, function(current_quant) {
    
    # create a dateframe for each quantile and bind together:
    current_q_df <- quant_df %>% 
      # filter to the first quantile by converting the string to a symbol so R recognizes it as a column name
      filter(quant_flag == current_quant) %>% 
      # sum all smoke concentrations and nonsmoke concentrations within a given quantile cutoff
      mutate(tot_pred_nonsmoke = sum(pred_nonsmoke, na.rm=TRUE),
             tot_pred_smoke = sum(pred_smoke, na.rm= TRUE),
             tot_pred_conc = sum(pred_allPM, na.rm = TRUE)) %>% 
      # how many days are we summing over in this quantile
      mutate(num_days = max(row_number())) %>% 
      mutate(conc_frac_attrib2smoke = 100*(tot_pred_smoke/tot_pred_conc)) %>% 
      distinct(species, num_days, quant_flag, q_conc, conc_frac_attrib2smoke, tot_pred_smoke, tot_pred_nonsmoke, tot_pred_conc) 

  }) %>% 
    bind_rows() # for all species, bind together
  
  # now create an order to the quantiles for plotting
  frac_species_smoke_df <- frac_species_smoke_df %>% 
    mutate(quant_flag = factor(quant_flag, 
                             levels = c("bottom 25%", "middle 50%", "top 25%"))) %>% 
    mutate(spec_type = case_when(
      species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
      species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
      species=="EC" | species=="OC" ~ "Secondary organic",
      species=="S" ~ "Toxicity potentiator",
      TRUE ~ NA)) 
  
}) %>% 
  bind_rows() %>%  
  distinct()#  %>% 
  # mutate(q_val = prettyNum(round(q_val, digits = 5)))

# adjust dataframe so that it can plot stacked area plot
quantiles_long_df <- summ_stat_quantiles_df %>% 
  dplyr::select(spec_type, species, quant_flag, tot_pred_smoke, tot_pred_nonsmoke, tot_pred_conc) %>% 
  pivot_longer(cols = c(tot_pred_smoke, tot_pred_nonsmoke, tot_pred_conc), 
               names_to = 'pm_type', values_to = 'tot_pred_conc_in_quantile') %>% 
  mutate(pm_type = factor(pm_type, levels = c("tot_pred_conc", 
                                              "tot_pred_nonsmoke",
                                              "tot_pred_smoke")))  


# PLOT FIG 3: MAKE STACKED AREA PLOTS
stacked_species_quantile_bar_plot <- ggplot(quantiles_long_df,
                                             aes(x = quant_flag,
                                                 y = tot_pred_conc_in_quantile,
                                                 fill = pm_type
                                             )) + 
  geom_bar(position="fill", stat="identity", alpha =.6) +
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c('grey60',"steelblue", "coral"),
                    labels = c('Total Species Conc in PM2.5 (ug/m3)',
                               'Nonsmoke PM2.5 Species Conc (ug/m3)',
                               'Smoke PM2.5 Species Conc (ug/m3)'
                               )) +
  facet_wrap(~species, scales = 'free', ncol = 2) +
  scale_x_discrete(limits = c("bottom 25%", "middle 50%", "top 25%")) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(x = 'Quantiles of Predicted Species Concentration (ug/m3)',
       y = paste('Percent Attributable to PM2.5 Type'),
       title = 'Contribution of smoke and nonsmoke PM2.5 to species concentration') +
  theme_minimal()
stacked_species_quantile_bar_plot


species_frac_plot <- ggplot(summ_stat_quantiles_df,
                            aes(x = quant_flag,
                                y = conc_frac_attrib2smoke, 
                                color = spec_type, group = spec_type)) + # to get in percent
  geom_line(size=.5) +
  geom_point(size=1.5, shape = 16, alpha = 0.6) +
  # geom_text(aes(label = species), vjust = -2, color = 'black', size = 2) +
  # geom_text(aes(label = prettyNum(num_days, big.mark = ",")), vjust = -2, color = 'black', size = 2) +
  facet_wrap(~species, scales = 'free', ncol = 2) +
  # geom_text(aes(label = prettyNum(q_val)), vjust = -6, color = 'black', size = 2) +
  scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
                                       # facet_wrap(spec_type ~ species, scales = 'free', ncol = 3) +
  labs(x = 'Daily Species Concentration Quantiles',
       y = 'Percent Attributable to Smoke PM2.5',
       color = 'Species',
       title = "Average contribution of smoke and nonsmoke PM2.5 to species concentration at different concentration quantiles") +
  theme_minimal()

species_frac_plot


# ------------------------------------------------------------------------------------
# RUN SPECIES BY REGION
# ------------------------------------------------------------------------------------
region_reg = feols(c(conc_val_V, conc_val_PB, conc_val_FE, 
                     conc_val_CU, conc_val_CR, conc_val_MN, conc_val_ZN, conc_val_NI,
                     conc_val_EC, conc_val_OC, conc_val_S)
                   ~ conc_val_smokePM + conc_val_nonsmokePM | 
                     monitor_month + year, fsplit = ~region, data = selected_spec_df) # add in region,run a model for each region


# Save coefficients to table
pred_region_coeffs <- coeftable(region_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         region = 'sample') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(species = str_remove(lhs, 'conc_')) %>% 
  mutate(pm_type = case_when(
    coefficient == 'conc_val_smokePM' ~ 'smokePM',
    coefficient == 'conc_val_nonsmokePM' ~ 'nonsmokePM')) %>% 
  dplyr::select(-lhs, -id, - sample.var) 


##################
# How much of total exposure to a given species, 
# or total exposure above some concentration threshold, 
# is driven by wildfire smoke vs other PM sources?
##################

# get betas to then weight smoke and nonsmoke concentrations
region_betas <- pred_region_coeffs %>% 
  mutate(smoke_day = case_when(
    pm_type == 'smokePM' ~ 1,
    pm_type == 'nonsmokePM' ~ 0
  )) %>% 
  dplyr::select(region, smoke_day, Estimate, species)



# get list of regions:
region_list <- unique(pred_region_coeffs$region)
# test
# current_region <- region_list[4]

all_regions_df <- purrr::map_df(region_list, function(current_region) {
  
  current_region_betas <- betas %>% 
    filter(region == current_region)
  
  # test
  # current_species <- species_list[4]
  # map over each species
  all_spec_quantiles_df <- purrr::map_df(species_list, function(current_species) {
    
    # step one: lead
    current_weighted_conc_df <- selected_spec_df %>% 
      rename_at(vars(matches("^conc_")), ~ str_remove(., "^conc_")) %>% 
      dplyr::select(-c(starts_with('region_conc'), State, season, doy, moy, monitor_month)) %>% 
      dplyr::select(region, Dataset, SiteCode, Date, region, smoke_day, conc_val_smokePM, conc_val_nonsmokePM, !!sym(current_species)) %>% 
      left_join(current_region_betas %>% 
                  filter(species == current_species), 
                by = c('smoke_day', 'region')) %>% 
      mutate(pred_nonsmoke = Estimate* non_smokePM_cont,
             pred_smoke = Estimate*smokePM_pred) %>% 
      mutate(pred_allPM = pred_nonsmoke+pred_smoke) %>% 
      distinct()
    
    # get quantiles for summed weighted total PM
    current_quantiles <- quantile(current_weighted_conc_df$pred_allPM, seq(0, 1, .1), na.rm = TRUE)
    
    # create dataframe where all concentrations of a species have been identified to be above or below
    quant_df <- current_weighted_conc_df %>%
      mutate(q0 = ifelse(pred_allPM >= current_quantiles[1], '0%', NA),
             q10 = ifelse(pred_allPM >= current_quantiles[2], '10%', NA),
             q20 = ifelse(pred_allPM >= current_quantiles[3], '20%', NA),
             q30 = ifelse(pred_allPM >= current_quantiles[4], '30%', NA),
             q40 = ifelse(pred_allPM >= current_quantiles[5], '40%', NA),
             q50 = ifelse(pred_allPM >= current_quantiles[6], '50%', NA),
             q60 = ifelse(pred_allPM >= current_quantiles[7], '60%', NA),
             q70 = ifelse(pred_allPM >= current_quantiles[8], '70%', NA),
             q80 = ifelse(pred_allPM >= current_quantiles[9], '80%', NA),
             q90 = ifelse(pred_allPM >= current_quantiles[10], '90%', NA),
             q100 = ifelse(pred_allPM >= current_quantiles[11], '100%', NA)) 
    
    
    # get list of quantiles
    q_list <- c('q0', 'q10', 'q20', 'q30', 'q40', 'q50', 'q60', 'q70', 'q80', 'q90', 'q100')
    # current_quant <- q_list[1] # test
    
    # map over the each quantile and create a dataframe for the specific species
    # that should help us plot how much of total exposure to a given species is driven by wildfire smoke
    frac_species_smoke_df <- purrr::map_df(q_list, function(current_quant) {
      
      # create a dateframe for each quantile and bind together:
      current_q_df <- quant_df %>% 
        # filter to the first quantile by converting the string to a symbol so R recognizes it as a column name
        filter(!is.na(!!sym(current_quant))) %>% 
        mutate(q_val = min(pred_allPM)) %>% 
        mutate(tot_pred_nonsmoke = sum(pred_nonsmoke, na.rm=TRUE),
               tot_pred_smoke = sum(pred_smoke, na.rm= TRUE)) %>% 
        mutate(num_days = max(row_number())) %>% 
        mutate(fraction_smoke = 100*(tot_pred_smoke/(tot_pred_smoke + tot_pred_nonsmoke))) %>% 
        distinct(region, species, num_days, !!sym(current_quant), q_val, fraction_smoke, tot_pred_smoke, tot_pred_nonsmoke) %>% 
        rename(quantile = !!sym(current_quant)) 
      
    }) %>% 
      bind_rows() 
    
    # now create an order to the quantiles for plotting
    frac_species_smoke_df <- frac_species_smoke_df %>% 
      mutate(quantile = factor(quantile, 
                               levels = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")))
    
    
  }) %>% 
    bind_rows() %>% # end binding for all species 
    distinct()
  
}) %>% 
  bind_rows() # end for all regions

# PLOT FIG 3 EXAMPLE:
regional_species_frac_plot <- ggplot(all_regions_df, 
                                     aes(x = quantile, 
                                         y = fraction_smoke, color = species, group = species)) + # to get in percent
  geom_line(size=.5) +
  geom_point(size=3.5, shape = 16, alpha = 0.6) + 
  #scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod', 'navy')) +
  # geom_hline(yintercept = frac_species_smoke_dfquantile, linetype = "dashed", color = 'grey') +
  facet_wrap(~region) +
  ylim(0, 50) +
  labs(x = 'Daily Concentration >= quantile',
       y = '% of contribution from wildfire smoke',
       color = 'Species',
       title = "Contribution of wildfire smoke to overall and extreme daily species concentrations by region") +
  theme_light() 
regional_species_frac_plot 


# ------------------------------------------------------------------------------------
# RUN SPECIES BY MONITOR
# ------------------------------------------------------------------------------------
test <- reg_df %>% 
  dplyr::select(year, Dataset, SiteCode, Date, conc_MF, conc_RC_PM2, smokePM_pred, non_smokePM_cont, conc_OC)

# run a separate regression for CSN vs IMPROVE
predPM_reg = feols(c(conc_MF, conc_RC_PM2, conc_NH4, conc_tot_metals, 
                     conc_ammNO3, conc_ammSO4, conc_AL,conc_AS,conc_BR,conc_CA,
                     conc_EC,conc_OC,conc_CHL,conc_CL,conc_CR, conc_CU,conc_FE,conc_K,conc_MG, conc_MN,
                     conc_NA,conc_NI,conc_NO3,conc_N2,
                     conc_P,conc_PB,conc_RB,conc_S,conc_SE,
                     conc_SI,conc_SOIL,conc_SO4,conc_SR, conc_V,conc_ZN,
                     conc_ZR) ~ smokePM_pred + non_smokePM_cont | monitor_month + year, reg_df, fsplit = ~Dataset) # add in region,run a model for each region

# Save coefficients to table
predPM_coeffs <- coeftable(predPM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         Dataset = 'sample') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(species = str_remove(lhs, 'conc_')) %>% 
  mutate(pm_type = case_when(
    coefficient == 'smokePM_pred' ~ 'smoke PM',
    coefficient == 'non_smokePM_cont' ~ 'non-smoke PM')) %>% 
  dplyr::select(-lhs, -id, - sample.var) 


# GET BETAS FOR SPECIES TO START
# set up list of species we care about
#species_list <- unique(predPM_coeffs$species) # c('PB', 'V', 'CR', 'SO4', 'NO3', 'AS')
# species_list <- c('OC', 'EC', 'RC_PM2', 'MF', )
species_list <- c('tot_metals', 'NH4', 'ammNO3', 'ammSO4', 'AL', 'AS', 'BR', 'CA', 'CHL', 'CL',
                  'CR', 'CU', 'FE', 'K', 'MG', 'MN', 'NA', 'NI', 'NO3', 'N2', 'P', 'RB',
                  'RB', 'S', 'SE', 'SI', 'SOIL', 'SO4', 'SR', 'V', 'ZN', 'ZR')
# species_list <- c('tot_metals', 'NH4', 'ammNO3', 'ammSO4', 
#                   'PB', 'V', 'CR', 'AS', 'FE')


# get betas to then weight smoke and nonsmoke concentrations
betas <- predPM_coeffs %>% 
  filter(species %in% species_list) %>% 
  mutate(smoke_day = case_when(
    coefficient == 'smokePM_pred' ~ 1,
    coefficient == 'non_smokePM_cont' ~ 0
  )) %>% 
  dplyr::select(Dataset, smoke_day, Estimate, species)


# get list of regions:
monitor_list <- unique(predPM_coeffs$Dataset)
# test
# current_monitor <- monitor_list[2]

all_monitors_df <- purrr::map_df(monitor_list, function(current_monitor) {
  
  current_region_betas <- betas %>% 
    filter(Dataset == current_monitor)
  
  # test
  # current_species <- species_list[4]
  # map over each species
  all_spec_quantiles_df <- purrr::map_df(species_list, function(current_species) {
    
    # step one: lead
    current_weighted_conc_df <- reg_df %>% 
      rename_at(vars(matches("^conc_")), ~ str_remove(., "^conc_")) %>% 
      dplyr::select(-c(starts_with('region_conc'), State, season, doy, moy, 
                       station_calc_smokePM, station_non_smokePM, monitor_month)) %>% 
      dplyr::select(Dataset, SiteCode, Date, region, smoke_day, smokePM_pred, non_smokePM_cont, !!sym(current_species)) %>% 
      left_join(current_region_betas %>% 
                  filter(species == current_species), 
                by = c('smoke_day', 'Dataset')) %>% 
      mutate(pred_nonsmoke = Estimate* non_smokePM_cont,
             pred_smoke = Estimate*smokePM_pred) %>% 
      mutate(pred_allPM = pred_nonsmoke+pred_smoke) %>% 
      distinct()
    
    # get quantiles for summed weighted total PM
    current_quantiles <- quantile(current_weighted_conc_df$pred_allPM, seq(0, 1, .1), na.rm = TRUE)
    
    # create dataframe where all concentrations of a species have been identified to be above or below
    quant_df <- current_weighted_conc_df %>%
      mutate(q0 = ifelse(pred_allPM >= current_quantiles[1], '0%', NA),
             q10 = ifelse(pred_allPM >= current_quantiles[2], '10%', NA),
             q20 = ifelse(pred_allPM >= current_quantiles[3], '20%', NA),
             q30 = ifelse(pred_allPM >= current_quantiles[4], '30%', NA),
             q40 = ifelse(pred_allPM >= current_quantiles[5], '40%', NA),
             q50 = ifelse(pred_allPM >= current_quantiles[6], '50%', NA),
             q60 = ifelse(pred_allPM >= current_quantiles[7], '60%', NA),
             q70 = ifelse(pred_allPM >= current_quantiles[8], '70%', NA),
             q80 = ifelse(pred_allPM >= current_quantiles[9], '80%', NA),
             q90 = ifelse(pred_allPM >= current_quantiles[10], '90%', NA),
             q100 = ifelse(pred_allPM >= current_quantiles[11], '100%', NA)) 
    
    
    # get list of quantiles
    q_list <- c('q0', 'q10', 'q20', 'q30', 'q40', 'q50', 'q60', 'q70', 'q80', 'q90', 'q100')
    # current_quant <- q_list[1] # test
    
    # map over the each quantile and create a dataframe for the specific species
    # that should help us plot how much of total exposure to a given species is driven by wildfire smoke
    frac_species_smoke_df <- purrr::map_df(q_list, function(current_quant) {
      
      # create a dateframe for each quantile and bind together:
      current_q_df <- quant_df %>% 
        # filter to the first quantile by converting the string to a symbol so R recognizes it as a column name
        filter(!is.na(!!sym(current_quant))) %>% 
        mutate(q_val = min(pred_allPM)) %>% 
        mutate(tot_pred_nonsmoke = sum(pred_nonsmoke, na.rm=TRUE),
               tot_pred_smoke = sum(pred_smoke, na.rm= TRUE)) %>% 
        mutate(num_days = max(row_number())) %>% 
        mutate(fraction_smoke = 100*(tot_pred_smoke/(tot_pred_smoke + tot_pred_nonsmoke))) %>% 
        distinct(Dataset, species, num_days, !!sym(current_quant), q_val, fraction_smoke, tot_pred_smoke, tot_pred_nonsmoke) %>% 
        rename(quantile = !!sym(current_quant)) 
      
    }) %>% 
      bind_rows() 
    
    # now create an order to the quantiles for plotting
    frac_species_smoke_df <- frac_species_smoke_df %>% 
      mutate(quantile = factor(quantile, 
                               levels = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")))
    
    
  }) %>% 
    bind_rows() %>% # end binding for all species 
    distinct()
  
}) %>% 
  bind_rows() # end for all regions

# PLOT FIG 3 EXAMPLE:
montior_species_frac_plot <- ggplot(all_monitors_df, 
                                    aes(x = quantile, 
                                        y = fraction_smoke, color = species, group = species)) + # to get in percent
  geom_line(size=.5) +
  geom_point(size=3.5, shape = 16, alpha = 0.6) + 
  #scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum', 'goldenrod', 'navy')) +
  # geom_hline(yintercept = frac_species_smoke_dfquantile, linetype = "dashed", color = 'grey') +
  facet_wrap(~Dataset) +
  ylim(0, 50) +
  labs(x = 'Daily Concentration >= quantile',
       y = '% of contribution from wildfire smoke',
       color = 'Species',
       title = "Contribution of wildfire smoke to overall and extreme daily species concentrations by Monitor") +
  theme_light() 
montior_species_frac_plot 

