# # # Emma Krasovich Southworth, emmars@stanford.edu, Jeff Wen
# # # Last Updated: April 19, 2024
# # # Description: run regression on burned structures and smoke speciation + plot response curves at different levels of smoke
#
# loadd(c(burned_struc_smoke_spec_df, parameter_categories, permutation_results_df), cache = drake_cache)
#
# ----------------------------------------------------------------
#  RUN PERMUTATION TEST TO TEST SHARP NULL
# ----------------------------------------------------------------
plot_burned_structure_coeffs_rand_plot <- function(burned_struc_smoke_spec_df, parameter_categories, permutation_results_df) {
  
  TREATMENT_VAR <- "contrib_daily_structures_destroyed"
  INTEREST_VAR <- "smokePM:contrib_daily_structures_destroyed"
  LINE_COLOR <- "firebrick"
  DEFAULT_COLOR <- "gray60"
  
  # set up regression dataframe
  reg_df <- burned_struc_smoke_spec_df %>%
    filter(contrib_daily_structures_destroyed > 0) %>%
    mutate(monitor_month = paste0(site_id, "-", month))
  
  # Run a regression using your sample
  sample_mod = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE,
                       K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                       S,  SE, SI, SO4, SR, TI, V,  ZN)
                     ~ smokePM*contrib_daily_structures_destroyed |
                       monitor_month + year,
                     reg_df, cluster = 'site_id')
  # get coeffs
  sample_coeffs <- coeftable(sample_mod) %>%
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs') %>%
    mutate(pval = round(pval, digits = 3)) %>%
    dplyr::select(-id)  %>%
    left_join(confint(
      sample_mod) %>%
        rename(species = 'lhs',
               CI25 = '2.5 %',
               CI975 = '97.5 %') %>%
        dplyr::select(-id))
  
 avg_bs <- mean(reg_df$contrib_daily_structures_destroyed, na.rm=TRUE)
  
  coef_df <- sample_coeffs %>% 
    filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
    left_join(parameter_categories) %>% 
    filter(species %in% c('MG','TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR'))  
    #mutate(Estimate = Estimate*avg_bs)
  
  coef_df$species_long <- factor(coef_df$species_long, 
                                 levels = c("Magnesium (Mg)", "Chromium (Cr)", "Nickel (Ni)", 
                                            "Copper (Cu)", "Zinc (Zn)", "Manganese (Mn)",
                                            "Arsenic (As)", "Lead (Pb)", "Titanium (Ti)"))
  
  permutation_results_df$species_long <- factor(permutation_results_df$species_long, 
                                                levels = c("Magnesium (Mg)", "Chromium (Cr)", "Nickel (Ni)", 
                                                           "Copper (Cu)", "Zinc (Zn)", "Manganese (Mn)",
                                                           "Arsenic (As)", "Lead (Pb)", "Titanium (Ti)"))
  
  limits <- c("Magnesium (Mg)", "Chromium (Cr)", "Nickel (Ni)", 
              "Copper (Cu)", "Zinc (Zn)", "Manganese (Mn)",
              "Arsenic (As)", "Lead (Pb)", "Titanium (Ti)")
  
  # get percentiles of rand tests
  merged_coeffs <- permutation_results_df %>% 
    left_join(coef_df) %>%    
    group_by(species) %>% 
    mutate(coeff_over_est = case_when(
      Estimate > 0 ~ sum(coef > Estimate),
      Estimate < 0 ~ sum(coef < Estimate))) %>% 
    ungroup() %>%
    distinct(species, coeff_over_est) %>% 
    mutate(pval = coeff_over_est/1000) %>% 
    left_join(coef_df %>% 
                dplyr::select(-pval)) %>% 
    dplyr::select(species, obs_est = 'Estimate', pval) %>% 
    arrange(pval) 
    
  
  print(xtable(merged_coeffs, digits =-2), include.rownames=FALSE)
  
  
  # # plot permuation test
  p_perm_test <- ggplot(data=permutation_results_df %>% 
                          mutate(coef = coef*avg_bs), aes(color=species_long)) +
    theme_minimal() +
    theme(
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_markdown()
    ) +
    geom_histogram(aes(x=coef, fill=species_long), bins=100, alpha=0.3, 
                   position="identity", size=0.25) +
    geom_vline(data=coef_df, 
               aes(xintercept=Estimate), lty = 2, color='black', size=0.25) +
    #coord_cartesian(xlim=c(-1e-6,1e-6), expand=T) +
    facet_wrap(species_long~., scales="free", ncol = 1) +
    labs(x="Coefficient",
         y="Frequency<br>", fill="", color="") +
    scale_color_manual(values= c("#35274A", "#E58601", "#E58601", "#E58601","#E58601", "#E58601",
                                 "#E2D200", "#46ACC8", "#46ACC8")) + 
    scale_fill_manual(values= c("#35274A", "#E58601", "#E58601", "#E58601","#E58601", "#E58601",
                                "#E2D200", "#46ACC8", "#46ACC8"))

  p_perm_test
  
  # save file
  ggsave(filename = 'Fig4B_rand_inference_plot.pdf',
         plot = p_perm_test,
         path = file.path(results_fp, 'figures/Fig4'),
         scale = 1,
         width = 8,
         height = 10,
         dpi = 320)
  
  ggsave(filename = 'Fig4B_rand_inference_plot.png',
         plot = p_perm_test,
         path = file.path(results_fp, 'figures/Fig4'),
         scale = 1,
         width = 8,
         height = 10,
         dpi = 320)
  
  
  return(perm_result_df)
}

