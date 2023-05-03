# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 30, 2023
# Description: testing harmonized data  vs our data

 loadd(clean_pm_spec_df, cache = drake_cache)
# set up list of species we care about
species_list <-c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")

# read in kara's data
harm_df <- readRDS(file.path(wip_gdrive_fp, 'intermediate/Kara_IMPROVE_CSN_smokeday_smokePM_data.rds'))

harm_spec_df <- harm_df %>% 
  filter(year != '0') %>% 
  filter(!is.na(smokePM_pred)) %>% 
  dplyr::select(Network, Site, State, year, month, date, day, 
                V, Pb, Fe, 
                Cu, Cr, Zn,
                EC, OC, smokePM = 'smokePM_pred', 
                totPM2.5 = 'PM2.5') %>% 
  mutate(nonsmokePM = totPM2.5- smokePM) %>% 
  mutate(pm_flag = case_when(
    totPM2.5 < smokePM ~'drop', 
    TRUE ~ 'keep')) %>% 
  filter(pm_flag == 'keep') %>% 
  mutate(monitor_month = paste0(Site,"_", month))

spec_baseline_avgs_df <- tibble(
  species = c("V", "Pb", "Fe", "Cu","Cr", "Mn", "Zn", "Ni", "EC", "OC", "S"),
  baseline_avg = c(mean(harm_spec_df$V, na.rm = TRUE), 
                   mean(harm_spec_df$Pb, na.rm = TRUE), 
                   mean(harm_spec_df$Fe, na.rm = TRUE), 
                   mean(harm_spec_df$Cu, na.rm = TRUE),  
                   mean(harm_spec_df$Cr, na.rm = TRUE),  
                   mean(harm_spec_df$Mn, na.rm = TRUE), 
                   mean(harm_spec_df$Zn, na.rm = TRUE),  
                   mean(harm_spec_df$Ni, na.rm = TRUE),  
                   mean(harm_spec_df$EC, na.rm = TRUE),  
                   mean(harm_spec_df$OC, na.rm = TRUE),  
                   mean(harm_spec_df$S, na.rm = TRUE))
)


# RUN MAIN MODEL REGRESSION
# get coeffs
kara_base_reg = feols(c(V, Pb, Fe, 
                   Cu, Cr, Zn,
                   EC, OC) ~ smokePM + nonsmokePM | 
                   monitor_month + year, harm_spec_df) # add in region,run a model for each region

res <- etable(kara_base_reg)

# get coefficients and prepare for plotting
kara_levelsPM_coeffs <- coeftable(kara_base_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  rename(species = 'lhs',
         pm_type = 'coefficient') %>% 
  dplyr::select(-id)  %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="Pb" | species=="Cr" ~ "Toxic metal",
    species=="Fe" | species=="Cu" | species=="Mn" | species=="Zn" | species=="Ni" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) 


# merge sample avg for each species and divide through betas by sample avg
karas_PMcoeffs_normalized <- kara_levelsPM_coeffs %>% 
  left_join(spec_baseline_avgs_df, by = 'species') %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_se = se/baseline_avg)

# plot coefficients for the speciation at the avg monitor ---------------------------
kara_norm_levels_reg_plot <- ggplot(karas_PMcoeffs_normalized, 
                               aes(x = species,
                                   y = 100*norm_est, 
                                   color=spec_type)) +
  geom_point(size=3, aes(shape = pm_type), alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(ymin = (norm_est - norm_se)*100,
                    ymax = (norm_est + norm_se)*100,
                    width =.2), stat = "identity") +
  scale_x_discrete(limits = c("V", "Cr", "Pb", "Cu", "Mn", "Ni", "Zn","Fe", "S", "EC", "OC")) +
  scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = '% Change (ug/m3)',
       x = 'Chemical Species (ug/m3)',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Percent change in species concentration relative to baseline for a 1 ug/m3 increase in smoke and nonsmoke",
       subtitle = 'Karas Harmonized Data') +
  theme_light() 
kara_norm_levels_reg_plot


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### NOW PLOT THE NORMAL REG ---------------------------------------------------
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

# MERGE ALL COEFFS TOGETHER:
all_coeffs <- bind_rows(karas_PMcoeffs_normalized %>% 
                          mutate(data = 'kara'), PMcoeffs_normalized %>% 
                          mutate(data = 'ours')) %>% 
  mutate(species = toupper(species))




# plot coefficients for the speciation at the avg monitor ---------------------------
all_levels_reg_plot <- ggplot(all_coeffs, 
                               aes(x = species,
                                   y = 100*norm_est, 
                                   color=data)) +
  geom_point(size=3, aes(shape = pm_type), alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(ymin = (norm_est - norm_se)*100, 
                    ymax = (norm_est + norm_se)*100, 
                    width =.2), stat = "identity") +
  scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c("red","navyblue")) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = '% Change (ug/m3)',
       x = 'Chemical Species (ug/m3)',
       color = 'Data',
       shape = 'PM Type',
       title = "Percent change in species concentration relative to baseline for a 1 ug/m3 increase in smoke and nonsmoke",
       subtitle = 'comparing coefficients between our data') +
  theme_light() 
all_levels_reg_plot



