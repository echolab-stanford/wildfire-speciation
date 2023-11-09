# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: model sensitivity checks
  #1: check that CSN and improve don't have really different coefficients
  #2: check that taking out nonsmoke PM does not change our coefficients
  #3: check these things too
# main + no negatives
# main + dropping things near MDLs
# main + monitor + year + month-of-year FE
# main + month of sample FE
# main + monitor + year + month FE

# loadd(c(parameter_categories, clean_PMspec_df), cache = drake::drake_cache(".drake"))

exploring_model_sensitivity <- function(parameter_categories, clean_PMspec_df)
  
  # select the vars that are needed for the
  reg_df <- clean_PMspec_df %>% 
    dplyr::select(Dataset, state_name, region, year, month, Date, 
                  monitor_month, smoke_day, site_id, MF_adj, smokePM, 
                  nonsmokePM_MF, AL:ZR) 
  

# -----------------------------------------------------------------------------
# 1. GET BASELINE AVERAGES FOR EACH SPECIES (FULL SAMPLE + REGIONAL SAMPLE)
# -----------------------------------------------------------------------------
# calculate the baseline for each chemical species across full sample by filtering to nonsmoke days
# rather than all days, bc some places may be really smoky 
# relative to others 
spec_ns_samp_avgs_df <- reg_df %>% 
  dplyr::select(site_id, Date, AL:ZR, smoke_day) %>% 
  pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>% 
  filter(smoke_day == 'nonsmoke day') %>% 
  group_by(species) %>% 
  dplyr::summarise(avg_nonsmoke_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')

  # -----------------------------------------------------------------------------
  # MAIN SPECIFICATION REGRESSION
  # -----------------------------------------------------------------------------
  # the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
  main_model = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                              K, MG,MN, `NA`, NI, NO3, N2, OC, P,  PB, RB,
                              S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                            ~ smokePM + nonsmokePM_MF | 
                              monitor_month + year, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs <- confint(
    main_model) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'MF') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  main_model_coefs <- coeftable(main_model) %>% 
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
    mutate(model = 'main')
  
  # -----------------------------------------------------------------------------
  # TAKE OUT NONSMOKE PM
  # -----------------------------------------------------------------------------
  # the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
  no_nonsmoke_model = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                       K, MG,MN, `NA`, NI, NO3, N2, OC, P,  PB, RB,
                       S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                     ~ smokePM  | 
                       monitor_month + year, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs <- confint(
    no_nonsmoke_model) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'MF') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  no_nonsmoke_model_coefs <- coeftable(no_nonsmoke_model) %>% 
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
    mutate(model = 'drop nonsmoke in model')
  
  
  # -----------------------------------------------------------------------------
  # TAKE OUT SMOKE PM
  # -----------------------------------------------------------------------------
  # the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
  smoke_model = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                              K, MG,MN, `NA`, NI, NO3, N2, OC, P,  PB, RB,
                              S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                            ~ nonsmokePM_MF  | 
                              monitor_month + year, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs <- confint(
    smoke_model) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'MF') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  smoke_model_coefs <- coeftable(smoke_model) %>% 
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
    mutate(model = 'drop smoke in model')
  
  
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
   
  
  # -----------------------------------------------------------------------------
  # CHANGE FEs in MAIN MODEL
  # -----------------------------------------------------------------------------
  # the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
  changeFEs_model = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                       K, MG,MN, `NA`, NI, NO3, N2, OC, P,  PB, RB,
                       S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                     ~ smokePM + nonsmokePM_MF | 
                       monitor_month + year, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs <- confint(
    main_model) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'MF') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  main_model_coefs <- coeftable(main_model) %>% 
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
    mutate(model = 'main')
  
  
  # -------------------------------------------------------------------------------
  # PLOT MAIN VS MAIN MODEL WITHOUT NONSMOKE COEFF
  # -------------------------------------------------------------------------------
 main_plot <- ggplot(all_coefs %>%
                       filter(model %in% c('main', 
                                           'drop nonsmoke in model')), # 'drop nonsmoke in model'
                     aes(x = species_long,
                         y = 100*norm_est, 
                         color= model)) +
    geom_point(size=3, alpha = 0.6, stat = "identity", 
               position = position_dodge(width = .5)) + #
    geom_linerange(aes(ymin = 100*norm_CI25, 
                       ymax = 100*norm_CI975), stat = "identity",
                   position = position_dodge(width = .5)) +
    # scale_shape_manual(values = c(15,16,17,19))+
    scale_color_manual(values = c('forestgreen', 'black'))+
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') + # 0.0000045184477
    # scale_y_continuous(trans='log10') + #, 
    #facet_wrap(~species_long) +
    labs(y = '% change relative to concentration baseline on nonsmoke days',
         x = 'Species',
         color = 'model type', 
         title = "Comparing main model to nonsmoke term") + 
    theme_light() +
    coord_flip() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') #+
    # theme(axis.text.x = element_text(angle = 90, hjust = 1))
  main_plot
  
  comp_monitor_plot <- ggplot(all_coefs %>%
                        filter(model %in% c('main', 'csn only', 
                                            'improve only')), # 'drop nonsmoke in model'
                      aes(x = species_long,
                          y = norm_est*100, 
                          color= model)) +
    geom_point(size=3, alpha = 0.6, stat = "identity", 
               position = position_dodge(width = .5)) + #
    geom_linerange(aes(ymin = norm_CI25*100, 
                       ymax = norm_CI975*100), stat = "identity",
                   position = position_dodge(width = .5)) +
    # scale_shape_manual(values = c(15,16,17,19))+
    scale_color_manual(values = c('steelblue', 'firebrick', 'black'))+
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') + # 0.0000045184477
    #scale_y_continuous(trans='log10') + #, 
    labs(y = '% change relative to concentration baseline on nonsmoke days',
         x = 'Species',
         color = 'model type', 
         title = "Comparing main model coeffs w/ CSN vs Improve only") + 
    theme_light() +
    coord_flip() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') #+
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))
  comp_monitor_plot
  
  
  
 

