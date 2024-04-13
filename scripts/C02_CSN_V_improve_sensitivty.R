# -----------------------------------------------------------------------------
# SEPARATE BY IMPROVE VS CSN
# -----------------------------------------------------------------------------
csn <- reg_df %>% 
  filter(Dataset == 'CSN')

improve <- reg_df %>% 
  filter(Dataset == 'IMPROVE')

# the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
CSN_model = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                    K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                    S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                  ~ smokePM  | 
                    monitor_month + year, data = csn, cluster = 'site_id') 

## calculate 95 CI%
CIs <- confint(
  CSN_model) %>% 
  rename(species = 'lhs',
         pm_type = 'coefficient',
         CI25 = '2.5 %',
         CI975 = '97.5 %') %>% 
  mutate(measure = 'MF') %>% 
  dplyr::select(-id)

# get coefficients and prepare for plotting
csn_model_coefs <- coeftable(CSN_model) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient') %>% 
  mutate(pval = round(pval, digits = 3)) %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id)  %>% 
  left_join(parameter_categories, by = 'species') %>% 
  mutate(spec_type = fct_relevel(spec_type,
                                 c("Heavy metal", "Toxic metal", "Toxic nonmetal", 
                                   "Toxicity potentiator", "Non-toxic metal", 
                                   "Non-toxic nonmetal", "Secondary inorganic", "Secondary organic"
                                 ))) %>% 
  mutate(measure = 'MF') %>% 
  left_join(CIs, by = c('species', 'pm_type', 'measure')) %>% 
  mutate(model = 'csn only')


# the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
improve_model = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                        K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                        S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                      ~ smokePM  | 
                        monitor_month + year, data = improve, cluster = 'site_id') 

## calculate 95 CI%
CIs <- confint(
  improve_model) %>% 
  rename(species = 'lhs',
         pm_type = 'coefficient',
         CI25 = '2.5 %',
         CI975 = '97.5 %') %>% 
  mutate(measure = 'MF') %>% 
  dplyr::select(-id)

# get coefficients and prepare for plotting
improve_model_coefs <- coeftable(improve_model) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient') %>% 
  mutate(pval = round(pval, digits = 3)) %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id)  %>% 
  left_join(parameter_categories, by = 'species') %>% 
  mutate(spec_type = fct_relevel(spec_type,
                                 c("Heavy metal", "Toxic metal", "Toxic nonmetal", 
                                   "Toxicity potentiator", "Non-toxic metal", 
                                   "Non-toxic nonmetal", "Secondary inorganic", "Secondary organic"
                                 ))) %>% 
  mutate(measure = 'MF') %>% 
  left_join(CIs, by = c('species', 'pm_type', 'measure')) %>% 
  mutate(model = 'improve only')


all_coefs <- main_model_coefs %>% 
  bind_rows(no_nonsmoke_model_coefs) %>% 
  bind_rows(csn_model_coefs) %>% 
  bind_rows(improve_model_coefs) %>% 
  filter(pm_type == 'smokePM') %>% 
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
  mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
         norm_CI25 = CI25/avg_nonsmoke_spec_conc,
         norm_CI975 = CI975/avg_nonsmoke_spec_conc)
