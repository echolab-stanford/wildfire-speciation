# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: make a time series plot for each species

# loadd(clean_pm_spec_df, cache = drake_cache)

# set up list of species we care about
species_list <-c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")

# select the vars that are needed for the
reg_df <- clean_pm_spec_df %>% 
  dplyr::select(Dataset:smoke_day, monitor_month,
                totPM2.5, 
                smokePM,
                nonsmokePM,
                V, PB, FE, CU, CR, MN, ZN, NI, EC, OC, S)

# calculate the sample averages for each chemical species:
spec_baseline_avgs_df <- tibble(
  species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
  baseline_avg = c(mean(reg_df$V, na.rm = TRUE), 
                   mean(reg_df$PB, na.rm = TRUE), 
                   mean(reg_df$FE, na.rm = TRUE), 
                   mean(reg_df$CU, na.rm = TRUE),  
                   mean(reg_df$CR, na.rm = TRUE),  
                   mean(reg_df$MN, na.rm = TRUE), 
                   mean(reg_df$ZN, na.rm = TRUE),  
                   mean(reg_df$NI, na.rm = TRUE),  
                   mean(reg_df$EC, na.rm = TRUE),  
                   mean(reg_df$OC, na.rm = TRUE),  
                   mean(reg_df$S, na.rm = TRUE))
)

# -----------------------------------------------------------------------------
# RUN REGRESSION FOR SPECIES IN LEVELS + PLOT IN LOGS
# -----------------------------------------------------------------------------
levelsPM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                   ~ smokePM + nonsmokePM | 
                     monitor_month + year, reg_df) 

etable(levelsPM_reg)

# get coefficients and prepare for plotting
levelsPM_coeffs <- coeftable(levelsPM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  rename(species = 'lhs',
         pm_type = 'coefficient') %>% 
  dplyr::select(-id)  %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) 


# merge sample avg for each species and divide through betas by sample avg
PMcoeffs_normalized <- levelsPM_coeffs %>% 
  left_join(spec_baseline_avgs_df, by = 'species') %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_se = se/baseline_avg)

# plot coefficients for the speciation at the avg monitor ---------------------------
norm_levels_reg_plot <- ggplot(PMcoeffs_normalized, 
                        aes(x = species,
                            y = 100*norm_est, 
                            color=spec_type)) +
                            # ymin = pmax(0, Estimate - se), 
                            # ymax = Estimate + se)) +
  geom_point(size=3, aes(shape = pm_type), alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(ymin = (norm_est - norm_se)*100, 
                    ymax = (norm_est + norm_se)*100, 
                    width =.2), stat = "identity") +
  scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  # scale_y_continuous(trans = log10_trans(),
  #                    breaks = trans_breaks("log10", function(x) 10^x)) +
  #                    #labels = scales::percent_format(accuracy = .01)) +
  scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = '% Change (ug/m3)',
       x = 'Chemical Species (ug/m3)',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Percent change in species concentration relative to baseline for a 1 ug/m3 increase in smoke and nonsmoke",
       subtitle = 'nonsmoke + smoke on RHS, monitor-month and year FE') +
  theme_light() 
norm_levels_reg_plot


# CHECK 1) ARE BETAS THE SAME FOR WHEN WE RUN TOTAL PM ON LHS WITH AND W/O SMOKE
nonsmoke_control_reg = feols(conc_val_totPM2.5 ~ conc_val_smokePM + conc_val_nonsmokePM | monitor_month + year, conc_rel2avg_df) 
summary(nonsmoke_control_reg)
no_control_reg = feols(conc_val_totPM2.5 ~ conc_val_smokePM | monitor_month + year, conc_rel2avg_df) 
summary(no_control_reg)
nonosmoke_control_reg = feols(conc_val_totPM2.5 ~ conc_val_nonsmokePM | monitor_month + year, conc_rel2avg_df) 
summary(nonosmoke_control_reg)

# ---------------------------------------------------------------------------
# CHECK 2) SPLIT REGRESSION BY CSN VS IMPROVE
# ---------------------------------------------------------------------------
# CSN SITES ONLY
csn_only <- reg_df %>% 
  filter(Dataset == 'EPACSN')

csn_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN)
                     ~ smokePM + nonsmokePM | 
                       monitor_month + year, csn_only)

# IMPROVE SITES ONLY
improve_only <- reg_df %>% 
  filter(Dataset == 'IMPROVE')

improve_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN)
                ~ smokePM + nonsmokePM | 
                  monitor_month + year, improve_only)
                               
# get coeffs for both improve and csn
site_coeffs <- coeftable(improve_reg) %>%
  mutate(site = 'IMPROVE') %>% 
  bind_rows(coeftable(csn_reg) %>% 
              mutate(site = 'EPACSN')) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) 
 
# plot the comparison between CSN vs Improve
site_reg_plot <- ggplot(site_coeffs, 
                        aes(x = species,
                            y = Estimate, 
                            color=spec_type)) +
  geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
  scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(accuracy = .01)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  facet_wrap(~site)+
  labs(x = 'Chemical Species (ug/m3)',
       y = 'Estimate',
       color = 'Site Type',
       shape = 'PM Type',
       title = "Comparing CSN vs Improve Speciation") +
  theme_light() 
site_reg_plot

# calculate the sample averages for each chemical species:
csn_baseline_avgs_df <- tibble(
  species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
  baseline_avg = c(mean(csn_only$V, na.rm = TRUE), 
                   mean(csn_only$PB, na.rm = TRUE), 
                   mean(csn_only$FE, na.rm = TRUE), 
                   mean(csn_only$CU, na.rm = TRUE),  
                   mean(csn_only$CR, na.rm = TRUE),  
                   mean(csn_only$MN, na.rm = TRUE), 
                   mean(csn_only$ZN, na.rm = TRUE),  
                   mean(csn_only$NI, na.rm = TRUE),  
                   mean(csn_only$EC, na.rm = TRUE),  
                   mean(csn_only$OC, na.rm = TRUE),  
                   mean(csn_only$S, na.rm = TRUE)),
  site ='EPACSN'
) %>% 
  bind_rows(tibble(
    species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
    baseline_avg = c(mean(improve_only$V, na.rm = TRUE), 
                     mean(improve_only$PB, na.rm = TRUE), 
                     mean(improve_only$FE, na.rm = TRUE), 
                     mean(improve_only$CU, na.rm = TRUE),  
                     mean(improve_only$CR, na.rm = TRUE),  
                     mean(improve_only$MN, na.rm = TRUE), 
                     mean(improve_only$ZN, na.rm = TRUE),  
                     mean(improve_only$NI, na.rm = TRUE),  
                     mean(improve_only$EC, na.rm = TRUE),  
                     mean(improve_only$OC, na.rm = TRUE),  
                     mean(improve_only$S, na.rm = TRUE)),
    site ='IMRPOVE'
  ))
  
# merge with coeffs:
norm_site_coeffs <- site_coeffs %>% 
  left_join(csn_baseline_avgs_df, by = c('site', 'species')) %>% 
  mutate(norm_coeff = Estimate/baseline_avg)

# PLOT NORMALIZED CSN VS IMPROVE
# plot the comparison between CSN vs Improve
site_reg_plot <- ggplot(norm_site_coeffs, 
                        aes(x = species,
                            y = Estimate, 
                            color=spec_type)) +
  geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
  scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(accuracy = .01)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  facet_wrap(~site)+
  labs(x = 'Chemical Species (ug/m3)',
       y = 'Estimate',
       color = 'Site Type',
       shape = 'PM Type',
       title = "Comparing CSN vs Improve Speciation, normalized") +
  theme_light() 
site_reg_plot


# ---------------------------------------------------------------------------
# RUN A MODEL WHERE EACH PM TYPE IS PLOTTED ON SAME
# ---------------------------------------------------------------------------
smokePM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ smokePM | monitor_month + year, reg_df)
nonsmokePM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ nonsmokePM | monitor_month + year, reg_df)
totPM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ totPM2.5 | monitor_month + year, reg_df)

# COMBINED COEFFS FOR PLOTTING
pm_type_coeffs <- coeftable(smokePM_reg) %>%
  mutate(pm_type = 'smoke PM2.5') %>% 
  bind_rows(coeftable(nonsmokePM_reg) %>%
              mutate(pm_type = 'nonsmoke PM2.5')) %>% 
  bind_rows(coeftable(totPM_reg) %>%
              mutate(pm_type = 'total PM2.5')) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) 

# join with species baseline avgs + divide to normalize betas
norm_pm_type_coeffs <- pm_type_coeffs %>% 
  left_join(spec_baseline_avgs_df, by = 'species') %>%
  mutate(norm_est = Estimate/baseline_avg)

# plot the normalized values for each smoke/nonsmoke/total pm
# plot the comparison between CSN vs Improve
pmtype_reg_plot <- ggplot(norm_pm_type_coeffs, 
                        aes(x = species,
                            y = Estimate, 
                            color=spec_type)) +
  geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
  scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(accuracy = .01)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  labs(x = 'Chemical Species (ug/m3)',
       y = 'Estimate',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Comparing PM Types, normalized") +
  theme_light() 
pmtype_reg_plot


return(reg_df)

}



