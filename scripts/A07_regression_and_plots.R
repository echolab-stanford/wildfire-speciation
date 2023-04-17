# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: make a time series plot for each species

# loadd(conc_rel2avg_df, cache = drake_cache)

# set up list of species we care about
species_list <-c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")
# ------------------------------------------------------------------------------------
# RUN SPECIES BY SPECIES REGRESSION - SPECIES DIVIDED BY SAMPLE BASELINE
# ------------------------------------------------------------------------------------
predPM_CR_reg = feols(c(CR_val_AL, CR_val_AS, CR_val_BR, CR_val_CA, 
                        CR_val_CHL, CR_val_CL, CR_val_CR, CR_val_CU, 
                        CR_val_EC, CR_val_FE, CR_val_K, CR_val_MG, 
                        CR_val_MN, CR_val_N2, CR_val_NA,
                        CR_val_NH4, CR_val_NI, CR_val_NO3, CR_val_OC, 
                        CR_val_OP, CR_val_P, CR_val_PB, CR_val_RB, CR_val_S,
                        CR_val_SE, CR_val_SI, CR_val_SO4, CR_val_SOIL, 
                        CR_val_SR, CR_val_TC, CR_val_TI,CR_val_V, CR_val_ZN,
                        CR_val_ZR, CR_val_ammNO3, CR_val_ammSO4, CR_val_tot_metals) ~ CR_val_smokePM + CR_val_nonsmokePM | 
                        monitor_month + year, conc_rel2avg_df) # add in region,run a model for each region

# saveRDS(predPM_CR_reg, 
#         file.path(wip_gdrive_fp, 'intermediate/base_speciation_model_baseline_adj.rds'))

# Save coefficients to table
predPM_CR_coeffs <- coeftable(predPM_CR_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(species = str_remove(lhs, 'CR_val_'),
         coefficient = str_remove(coefficient, 'CR_val_')) %>% 
  filter(species %in% species_list) %>% 
  dplyr::select(-lhs, -id) %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) 

saveRDS(predPM_CR_coeffs,
        file.path(wip_gdrive_fp, 'results/conc_ratio_model_reg_coeffs.rds'))

# plot coefficients for the speciation at the avg monitor ---------------------------
pred_reg_plot <- ggplot(predPM_CR_coeffs, 
                        aes(x = reorder(species, Estimate), 
                            y = Estimate, color=spec_type,
                            ymin = pmax(0, Estimate - se), 
                            ymax = Estimate + se)) +
  geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(width = 0.1), stat = "identity") +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  scale_x_discrete(limits = c("CR", "PB", "V",  "CU", "FE", "MN", "NI", "ZN", "EC", "OC", "S")) +
  # geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Species Type',
       shape = 'PM Type',
       title = "Speciation of smoke vs nonsmoke PM2.5",
       subtitle = 'concentrations divided by sample average') +
  theme_light() 
pred_reg_plot


# ------------------------------------------------------------------------------------
# RUN SPECIES BY SPECIES REGRESSION - SPECIES IN LOG
# ------------------------------------------------------------------------------------
predPM_reg = feols(c(log(conc_val_AL), log(conc_val_AS), log(conc_val_BR), log(conc_val_CA), 
                     log(conc_val_CHL), log(conc_val_CL), log(conc_val_CR), log(conc_val_CU), 
                     log(conc_val_EC), log(conc_val_FE), log(conc_val_K), log(conc_val_MG), 
                     log(conc_val_MN), log(conc_val_N2), log(conc_val_NA), log(conc_val_NH4), 
                     log(conc_val_NI), log(conc_val_NO3), log(conc_val_OC), log(conc_val_OP), 
                     log(conc_val_P), log(conc_val_PB), log(conc_val_RB), log(conc_val_S),
                     log(conc_val_SE), log(conc_val_SI), log(conc_val_SO4), log(conc_val_SOIL), 
                     log(conc_val_SR), log(conc_val_TC), log(conc_val_TI), log(conc_val_V), 
                     log(conc_val_ZN), log(conc_val_ZR), log(conc_val_ammNO3), log(conc_val_ammSO4)) 
                   ~ conc_val_smokePM + conc_val_nonsmokePM | monitor_month + year, conc_rel2avg_df) 

# saveRDS(predPM_reg,
#         file.path(wip_gdrive_fp, 'results/base_model_speciation.rds'))

# Save coefficients to table
predPM_coeffs <- coeftable(predPM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(species = str_remove(lhs, 'log[(]conc_val_'),
         coefficient = str_remove(coefficient, 'conc_val_')) %>% 
  mutate(species = str_remove(species, '[)]')) %>% 
  filter(species %in% species_list) %>% 
  dplyr::select(-lhs, -id)  %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) %>% 
  arrange(spec_type)


# saveRDS(predPM_coeffs,
#         file.path(wip_gdrive_fp, 'results/base_model_reg_coeffs.rds'))

# plot coefficients for the speciation at the avg monitor ---------------------------
pred_reg_plot <- ggplot(predPM_coeffs, 
                        aes(x = species,
                            y = Estimate*100, 
                            color=spec_type, 
                            ymin = pmax(0, Estimate*100 - se*100), 
                            ymax = Estimate*100 + se*100)) +
  geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
  geom_errorbar(aes(width = 0.1), stat = "identity") +
  scale_x_discrete(limits = c("CR", "PB", "V",  "CU", "FE", "MN", "NI", "ZN", "EC", "OC", "S")) +
  # scale_y_continuous(trans = log10_trans(),
  #                    breaks = trans_breaks("log10", function(x) 10^x),
  #                    labels = scales::percent_format(scale = 10)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  labs(y = '% Change',
       x = 'Chemical Species (ug/m3)',
       color = 'Species Type',
       shape = 'PM Type',
       title = "Speciation of smoke vs nonsmoke PM2.5, base model") +
  theme_light() 
pred_reg_plot