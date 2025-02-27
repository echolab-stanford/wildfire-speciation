# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 1, 2024
# Description: create plot demonstrating that the results are similar regardless of PM measure used

# loadd(c(parameter_categories, spec_pal, clean_PMspec_df), cache = drake::drake_cache(".drake"))

create_coeff_sensitivity_SI_figures <- function(clean_PMspec_df, parameter_categories, spec_pal) {
  
  # select the vars that are needed for the regression
  reg_df <- clean_PMspec_df %>% 
    dplyr::select(Dataset, state_name, region, year, month, Date, 
                  monitor_month, smoke_day, site_id, MF_adj, RCFM_adj, smokePM, 
                  nonsmokePM_MF, nonsmokePM_RCFM, AL:ZR) %>% 
    # drop Soil and zirconium, not needed for this analysis
    dplyr::select(-SOIL, -ZR)
  
  # -----------------------------------------------------------------------------
  # 1. GET BASELINE AVERAGES FOR EACH SPECIES (FULL SAMPLE + REGIONAL SAMPLE)
  # -----------------------------------------------------------------------------
  # calculate the baseline for each chemical species across full sample by filtering to nonsmoke days
  # rather than all days, bc some places may be really smoky 
  # relative to others 
  spec_ns_samp_avgs_df <- reg_df %>% 
    dplyr::select(site_id, Date, AL:ZN, smoke_day) %>% 
    pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>% 
    filter(smoke_day == 'nonsmoke day') %>% 
    group_by(species) %>% 
    dplyr::summarise(avg_nonsmoke_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')
  
  
  # -----------------------------------------------------------------------------
  # 3. RUN REGRESSION FOR SPECIES
  #   A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
  # -----------------------------------------------------------------------------
  # run original regression
  full_sampPM_regMF = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE, 
                              K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                              S,  SE, SI, SO4, SR, TI, V,  ZN)
                            ~ smokePM + nonsmokePM_MF | 
                              monitor_month + year, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs_MF <- confint(
    full_sampPM_regMF) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'MF') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  full_sampPM_coeffsMF <- coeftable(full_sampPM_regMF) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))
    ) %>% 
    mutate(measure = 'MF') %>% 
    left_join(CIs_MF, by = c('species', 'pm_type', 'measure'))
  
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  full_samp_PMcoeffs_normalized <- full_sampPM_coeffsMF %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM')
  
  # run regression for alternative PM measure ------------------------------
  sens_analysis = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                          K, MG,MN, `NA`, NI, NO3, OC, P, PB, RB,
                          S,  SE, SI, SO4, SR, TI, V,  ZN)
                        ~ smokePM + nonsmokePM_RCFM | 
                          monitor_month + year, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs_RCFM <- confint(
    sens_analysis) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'RCFM') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs_RCFM <- coeftable(sens_analysis) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))
    ) %>% 
    mutate(measure = 'RCFM') %>% 
    left_join(CIs_RCFM, by = c('species', 'pm_type', 'measure'))
  
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  coeffs_RCFM_normalized <- coeffs_RCFM %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM')
  
  
  # bind together:
  all_coeffs <- bind_rows(coeffs_RCFM_normalized, full_samp_PMcoeffs_normalized)
  
  # plot coefficients for speciation at the avg monitor which tells us how much 
  # of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
  # --------------------------------------------------------------------------------
  # plot percent change for all species
  # --------------------------------------------------------------------------------
  PM_robustness <- ggplot(all_coeffs, 
                          aes(x = species_long,
                              y = 100*norm_est, 
                              color=measure, 
                          )) +
    geom_point(size=3, alpha = 0.6, stat = "identity") +
    geom_linerange(aes(ymin = (100*norm_CI25), 
                       ymax = (100*norm_CI975)), stat = "identity") +
    scale_x_discrete(limits = c(
      # "Organics"
      "Organic Carbon (OC)", "Elemental Carbon (EC)",
      # "Halogens"
      "Bromine (Br)", "Chlorine (Cl)", 
      #  "Nonmetals"
      "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
      # "Other metals"
      "Aluminum (Al)", "Lead (Pb)",
      # "Metalloids"
      "Silicon (Si)", "Arsenic (As)",
      # "Transition metals"
      "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
      # "Alkali metals"
      "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
      # "Alkaline-earth metals"
      "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)"
    )) +
    scale_color_manual(values= c("steelblue", 'navy')) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(y = '% change in concentration relative to nonsmoke day baseline',
         x = 'Species',
         color = 'PM2.5 Measure', 
         title = "Results robust to different measures of PM2.5") + 
    theme_light() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50') +
    theme(legend.position = c(0.85, 0.90))
    
  PM_robustness
  
  # save file
  ggsave(
    filename = 'SIFig2A_PM_Robustness_checks.pdf',
    plot = PM_robustness,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 5,
    height = 7,
    dpi = 320) 
  
  ggsave(
    filename = 'SIFig2A_PM_Robustness_checks.png',
    plot = PM_robustness,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 5,
    height = 7,
    dpi = 320) 
  
  

  # -----------------------------------------------------------------------------
  # RUN REGRESSION WITH AND WITHOUT NONSMOKE + CHANGE MODEL SPECIFICATION
  # -----------------------------------------------------------------------------
  # run regression for alternative PM measure ------------------------------
  sens_analysis2 = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                          K, MG,MN, `NA`, NI, NO3, OC, P, PB, RB,
                          S,  SE, SI, SO4, SR, TI, V,  ZN)
                        ~ smokePM | 
                          monitor_month + year, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs_no_NS <- confint(
    sens_analysis2) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'no covariate') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs_no_NS <- coeftable(sens_analysis2) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))
    ) %>% 
    mutate(measure = 'no covariate') %>% 
    left_join(CIs_no_NS, by = c('species', 'pm_type', 'measure'))
  
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  coeffs_no_NS_normalized <- coeffs_no_NS %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM') 
   
  
  # bind together:
  all_coeffs2 <- bind_rows(coeffs_no_NS_normalized, full_samp_PMcoeffs_normalized) %>% 
    mutate(cov_flag = ifelse(measure == 'MF', 'nonsmoke as covariate', measure))
  
  # plot coefficients for speciation at the avg monitor which tells us how much 
  # of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
  # --------------------------------------------------------------------------------
  # plot percent change for all species
  # --------------------------------------------------------------------------------
  covariate_robustness <- ggplot(all_coeffs2, 
                          aes(x = species_long,
                              y = 100*norm_est, 
                              color=cov_flag, 
                          )) +
    geom_point(size=3, alpha = 0.7, stat = "identity") +
    geom_linerange(aes(ymin = (100*norm_CI25), 
                       ymax = (100*norm_CI975)), stat = "identity") +
    scale_x_discrete(limits = c(
      # "Organics"
      "Organic Carbon (OC)", "Elemental Carbon (EC)",
      # "Halogens"
      "Bromine (Br)", "Chlorine (Cl)", 
      #  "Nonmetals"
      "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
      # "Other metals"
       "Aluminum (Al)", "Lead (Pb)",
      # "Metalloids"
      "Silicon (Si)", "Arsenic (As)",
      # "Transition metals"
      "Manganese (Mn)", "Zinc (Zn)","Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
      # "Alkali metals"
      "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
      # "Alkaline-earth metals"
      "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)"
    )) +
    scale_color_manual(values= c("goldenrod1", 'coral'),
                       labels = c("No", "Yes"),
                       name = "Covariate Inclusion") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(y = '% change in concentration relative to nonsmoke day baseline',
         x = 'Species',
         title = "Consistent treatment effects regardless of covariate inclusion") + 
    theme_light() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    theme(legend.position = c(0.85, 0.90)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50')
  covariate_robustness
  
  
  # save file
  ggsave(
    filename = 'SIFig2C_Covariate_Robustness_checks.pdf',
    plot = covariate_robustness,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 5,
    height = 7,
    dpi = 320) 
  ggsave(
    filename = 'SIFig2C_Covariate_Robustness_checks.png',
    plot = covariate_robustness,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 5,
    height = 7,
    dpi = 320) 
  
  # -----------------------------------------------------------------------------
  # RUN DIFFERENT MODEL SPECIFICATIONS
  # -----------------------------------------------------------------------------
  # run regression for alternative PM measure ------------------------------
  mod1 = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                           K, MG,MN, `NA`, NI, NO3, OC, P, PB, RB,
                           S,  SE, SI, SO4, SR, TI, V,  ZN)
                         ~ smokePM + nonsmokePM_MF | 
                           site_id + year, reg_df, cluster = 'site_id') 
  ## calculate 95 CI%
  CIs1 <- confint(mod1) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(FEs = 'monitor + year') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs1 <- coeftable(mod1) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))) %>% 
    left_join(CIs1, by = c('species', 'pm_type'))
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  coeffs1_normalized <- coeffs1 %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM') 
  
  # run regression for alternative PM measure ------------------------------
  mod2 = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                 K, MG,MN, `NA`, NI, NO3, OC, P, PB, RB,
                 S,  SE, SI, SO4, SR, TI, V,  ZN)
               ~ smokePM + nonsmokePM_MF | 
                 site_id + month, reg_df, cluster = 'site_id') 
  
  ## calculate 95 CI%
  CIs2 <- confint(mod2) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(FEs = 'monitor + month') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs2 <- coeftable(mod2) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))
    ) %>% 
    left_join(CIs2, by = c('species', 'pm_type'))
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  coeffs2_normalized <- coeffs2 %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM') 
  
  
  # run regression for alternative PM measure ------------------------------
  mod3 = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                 K, MG,MN, `NA`, NI, NO3, OC, P, PB, RB,
                 S,  SE, SI, SO4, SR, TI, V,  ZN)
               ~ smokePM + nonsmokePM_MF | 
                 region + year, reg_df, cluster = 'region') 
  ## calculate 95 CI%
  CIs3 <- confint(mod3) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(FEs = 'region + year') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs3 <- coeftable(mod3) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))
    ) %>% 
    left_join(CIs3, by = c('species', 'pm_type'))
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  coeffs3_normalized <- coeffs3 %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM') 
  
  # run regression for alternative PM measure ------------------------------
  mod4 = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                 K, MG,MN, `NA`, NI, NO3, OC, P, PB, RB,
                 S,  SE, SI, SO4, SR, TI, V,  ZN)
               ~ smokePM + nonsmokePM_MF | 
                 region + month, reg_df, cluster = 'region') 
  ## calculate 95 CI%
  CIs4 <- confint(mod4) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(FEs = 'region + month') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs4 <- coeftable(mod4) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))
    ) %>% 
    left_join(CIs4, by = c('species', 'pm_type'))
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  coeffs4_normalized <- coeffs4 %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM') 
  
  
  mod5 = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                 K, MG,MN, `NA`, NI, NO3, OC, P, PB, RB,
                 S,  SE, SI, SO4, SR, TI, V,  ZN)
               ~ smokePM + nonsmokePM_MF , reg_df, cluster = 'site_id') 
  ## calculate 95 CI%
  CIs5 <- confint(mod5) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(FEs = 'none') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs5 <- coeftable(mod5) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 3)) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))) %>% 
    left_join(CIs5, by = c('species', 'pm_type'))
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  coeffs5_normalized <- coeffs5 %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM') 
  
  
  # bind together:
  all_coeffs_mods <- bind_rows(coeffs1_normalized, 
                               coeffs3_normalized, 
                               coeffs2_normalized, 
                               coeffs4_normalized,
                               coeffs5_normalized,
                               full_samp_PMcoeffs_normalized %>% 
                                 mutate(FEs = 'monitor-month + year (main)')) 
  all_coeffs_mods <- all_coeffs_mods %>% 
    mutate(FEs = fct_relevel(FEs, 
                             c("none",
                               "region + month", 
                               "region + year",
                               "monitor + month",
                               "monitor + year",
                               "monitor-month + year (main)")))
  
  
  # plot coefficients for speciation at the avg monitor which tells us how much 
  # of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
  # --------------------------------------------------------------------------------
  # plot percent change for all species
  # --------------------------------------------------------------------------------
  specification_robustness <- ggplot(all_coeffs_mods, 
                                 aes(x = species_long,
                                     y = 100*norm_est, 
                                     color=FEs, 
                                 )) +
    geom_point(size=3, alpha = 0.6, stat = "identity") +
    geom_linerange(aes(ymin = (100*norm_CI25), 
                       ymax = (100*norm_CI975)), stat = "identity") +
    scale_x_discrete(limits = c(
      # "Organics"
      "Organic Carbon (OC)", "Elemental Carbon (EC)",
      # "Halogens"
      "Bromine (Br)", "Chlorine (Cl)",  
      #  "Nonmetals"
      "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
      # "Other metals"
       "Aluminum (Al)", "Lead (Pb)",
      # "Metalloids"
      "Silicon (Si)", "Arsenic (As)",
      # "Transition metals"
      "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
      # "Alkali metals"
      "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
      # "Alkaline-earth metals"
      "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)" 
    )) +
    scale_color_manual(values= c("lightpink", 'palevioletred','deeppink', 'purple', 'slateblue4', 'grey15'),
                       name = "Model specification") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(y = '% change in concentration relative to nonsmoke day baseline',
         x = 'Species',
         title = "Consistent treatment effects regardless of model specification") + 
    theme_light() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    theme(legend.position = c(0.85, 0.80)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50')
  specification_robustness
  
  
  # save file
  ggsave(
    filename = 'SIFig2B_ModelSpec_Robustness_checks.pdf',
    plot = specification_robustness,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 5,
    height = 7,
    dpi = 320) 
  
  ggsave(
    filename = 'SIFig2B_ModelSpec_Robustness_checks.png',
    plot = specification_robustness,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 5,
    height = 7,
    dpi = 320) 
  
  # ----------------------------------------------------------------------------
  # R1, Q2: VALIDATE OUR MODEL BY INCLUDING ADDITIONAL ENVIRONMENTAL COVARIATES
  # ----------------------------------------------------------------------------
  # read in precipation data from CHIRPS
  ppt_files <- list.files(file.path(gee_data_fp), pattern = 'precip', full.names = TRUE)
  # current_fp <- ppt_files[1]
  
  ppt_df <- map_df(ppt_files, function(current_fp) {
    
    current_df <- read.csv(current_fp) %>% 
      dplyr::select(grid_id_10km = 'ID', ppt_mm_per_day = 'sum') %>% 
      mutate(year = gsubfn::strapply(current_fp, "\\d{4}", as.numeric, simplify = TRUE))
    
  })
  
  # read in temperature data from ERA5
  temp_files <- list.files(file.path(gee_data_fp), pattern = 'temp', full.names = TRUE)
  # current_fp <- temp_files[1]
  
  temp_df <- map_df(temp_files, function(current_fp) {
    
    current_df <- read.csv(current_fp) %>% 
      dplyr::select(grid_id_10km = 'ID', temp = 'mean_2m_air_temperature') %>% 
      mutate(year = gsubfn::strapply(current_fp, "\\d{4}", as.numeric, simplify = TRUE)) %>% 
      mutate(tempF = weathermetrics::kelvin.to.fahrenheit(temp, round = 2))
    
  })
  
  covar_df <-  left_join(temp_df, ppt_df) %>% 
    dplyr::select(-temp)
  
  # load these datasets
  loadd(c(clean_PMspec_df,clean_PMspec_sites_df, parameter_categories), cache = drake_cache)
  
  # Load 10 km grid
  grid_10km <- st_read(file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')) %>% 
    st_transform(4326)
  
  # Match monitor locations in speciation data with grid 10km ID ####
  # pull id cells that intersect with the monitoring sites
  pts_in_grid_df <- clean_PMspec_sites_df %>% 
    st_as_sf(coords = c('long', 'lat'), crs = 4326) %>% 
    st_join(grid_10km) %>% 
    st_drop_geometry() %>% 
    rename(grid_id_10km = 'ID') %>% 
    dplyr::select(-epsg, -COORDX, -COORDY)
  
  # Add grid_id_10km to the speciation df, join the smoke PM data by grid cell
  station_gridded_df <- clean_PMspec_df %>% 
    left_join(pts_in_grid_df, 
              by = c("site_id", 'Dataset')) %>% 
    left_join(covar_df) %>% 
    dplyr::select(Dataset:region, year, month, Date, 
                  monitor_month, grid_id_10km, smoke_day, site_id, 
                  tempF, ppt_mm_per_day, MF_adj, smokePM, 
                  nonsmokePM_MF, AL:ZN) %>% 
    # drop Soil and zirconium, not needed for this analysis
    dplyr::select(-SOIL)
  
  # GET BASELINE AVERAGES FOR EACH SPECIES --------------------------------
  # calculate the baseline for each chemical species across full sample by filtering to nonsmoke days
  # rather than all days, bc some places may be really smoky 
  # relative to others 
  spec_ns_samp_avgs_df <- station_gridded_df %>% 
    dplyr::select(site_id, Date, AL:ZN, smoke_day) %>% 
    pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>% 
    filter(smoke_day == 'nonsmoke day') %>% 
    group_by(species) %>% 
    dplyr::summarise(avg_nonsmoke_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')
  
  
  # RUN REGRESSION FOR SPECIES -------------------------------------------------
  # RUN IN LEVELS ACROSS FULL SAMPLE, 
  # DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
  
  covar_reg = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE, 
                      K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                      S,  SE, SI, SO4, SR, TI, V,  ZN)
                    ~ smokePM + nonsmokePM_MF + tempF + ppt_mm_per_day | 
                      monitor_month + year, station_gridded_df, cluster = 'site_id')
  # summary(covar_reg) # look at regression results
  ## calculate 95 CI%
  CIs <- confint(
    covar_reg) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'MF') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  covar_coeffs <- coeftable(covar_reg) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    mutate(pval = round(pval, digits = 4)) %>% 
    # get pvalues
    mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
    dplyr::select(-id)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals", 
                                        "Transition metals", "Metalloids", "Other metals", 
                                        "Nonmetals",  "Halogens", "Organics"))
    ) %>% 
    mutate(measure = 'MF') %>% 
    left_join(CIs, by = c('species', 'pm_type', 'measure'))
  
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  covar_coeffs_normalized <- covar_coeffs %>% 
    left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
    mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
           norm_CI25 = CI25/avg_nonsmoke_spec_conc,
           norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
    filter(pm_type == 'smokePM') %>% 
    mutate(model_type = 'environmental covariates included')
  
  # load the original coefficients
  loadd(spec_pal, full_samp_PMcoeffs_normalized)
  
  og_coeffs <- full_samp_PMcoeffs_normalized %>% 
    mutate(model_type = 'no environmental covariates')
  
  # add together:
  env_coeffs_coeffs <- bind_rows(covar_coeffs_normalized, og_coeffs)
  
  # plot percent change for all species ------------------------------------------
  pct_change_env_cov_plot <- ggplot(env_coeffs_coeffs, 
                                     aes(x = species_long,
                                         y = 100*norm_est, 
                                         color=model_type,
                                     )) +
    geom_point(size=3, alpha = 0.6, stat = "identity") +
    geom_linerange(aes(ymin = (100*norm_CI25), 
                       ymax = (100*norm_CI975)), stat = "identity") +
    scale_x_discrete(limits = c(
      # "Organics"
      "Organic Carbon (OC)", "Elemental Carbon (EC)",
      # "Halogens"
      "Bromine (Br)", "Chlorine (Cl)", 
      #  "Nonmetals"
      "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
      # "Other metals"
      "Aluminum (Al)", "Lead (Pb)",
      # "Metalloids"
      "Silicon (Si)", "Arsenic (As)",
      # "Transition metals"
      "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)","Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
      # "Alkali metals"
      "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
      # "Alkaline-earth metals"
      "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)" 
    )) +
    scale_color_manual(values= c('lightgreen', 'forestgreen')) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(y = expression(paste('% change relative to average nonsmoke day concentration')),
         x = 'Species',
         shape = 'Model inclusions',
         title = expression(paste("% change in concentration for 1 ug/", m^3, " increase in smoke ", "PM"["2.5"]))) + 
    theme_light() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50')
  pct_change_env_cov_plot
  
  
  # save file
  ggsave(
    filename = 'SIFig2D_EnvCov_Robustness_checks.png',
    plot = pct_change_env_cov_plot,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 7,
    height = 7,
    dpi = 320)  
  
  ggsave(
    filename = 'SIFig2D_EnvCov_Robustness_checks.pdf',
    plot = pct_change_env_cov_plot,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 7,
    height = 7,
    dpi = 320) 
  
  
  sensitivity_checks_coeffs <- bind_rows(all_coeffs, all_coeffs2, all_coeffs_mods, env_coeffs_coeffs)

  
  return(sensitivity_checks_coeffs)
  
}
  
  
  
  
  
  
  
  
  
