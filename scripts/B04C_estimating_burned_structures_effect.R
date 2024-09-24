# # Emma Krasovich Southworth, emmars@stanford.edu
# # Last Updated: April 2, 2024
# # Description: run regression on burned structures and smoke speciation + plot response curves at different levels of smoke

# loadd(c(burned_struc_smoke_spec_df, parameter_categories, spec_pal), cache = drake_cache)

# function
estimating_conc_response_2_burned_structures <- function(burned_struc_smoke_spec_df, 
                                                         parameter_categories, spec_pal) {

  # ----------------------------------------------------------------
  # STEP 1. RUN INTERACTION MODEL THAT LOOKS AT ASSOCIATED BW STRUCTURES + SPECIES CONC
  # ----------------------------------------------------------------
  # set up regression dataframe
  reg_df <- burned_struc_smoke_spec_df %>% 
    filter(!is.na(contrib_daily_structures_destroyed)) %>% 
    dplyr::select(ID, Dataset:site_id, contrib_daily_structures_destroyed,
                  smokePM:ZR, nonsmokePM_MF) %>% 
    filter_at(vars(AL, AS, BR, CA, CL, CR, CU,
                   EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
                   RB, S, SE, SI, SO4, SR, TI, V, ZN, ZR),
              any_vars(!is.na(.))) %>% 
   filter(contrib_daily_structures_destroyed > 0) %>% 
    mutate(monitor_month = paste0(site_id, "-", month)) %>% 
    mutate(region = case_when(
      state_name %in% c('Oregon', 'California', 'Washington') ~ 'pacific',
      state_name %in% c('Nevada', 'Utah', 'Colorado', 'Wyoming', 'Montana', 'Idaho') ~ 'rocky_mountain',
      state_name %in% c('Arizona', 'New Mexico', 'Texas', 'Oklahoma') ~ 'southwest',
      state_name %in% c('North Dakota', 'South Dakota', 'Nebraska', 'Kansas', 'Minnesota', 'Iowa', 
                        'Missouri', 'Wisconsin', 'Illinois', 'Indiana', 'Michigan', 'Ohio') ~ 'midwest',
      state_name %in% c('Kentucky', 'West Virginia', 'Virginia', 'Tennessee', 'North Carolina', 'Mississippi', 
                        'Alabama', 'Georgia', 'South Carolina', 'Florida', 'Louisiana', 'Arkansas') ~ 'southeast',
      state_name %in% c('Maine', 'Vermont', 'New Hampshire', 'Massachusetts', 'Rhode Island', 'Connecticut', 
                        'New York', 'New Jersey', 'Pennsylvania', 'Delaware', 'Maryland', 'District of Columbia') ~ 'northeast',
    )) 
  
# calculate nonsmoke day baseline
  no_structures_burned_smoke_conc_avg <- burned_struc_smoke_spec_df %>%
    mutate(region = case_when(
      state_name %in% c('Oregon', 'California', 'Washington') ~ 'pacific',
      state_name %in% c('Nevada', 'Utah', 'Colorado', 'Wyoming', 'Montana', 'Idaho') ~ 'rocky_mountain',
      state_name %in% c('Arizona', 'New Mexico', 'Texas', 'Oklahoma') ~ 'southwest',
      state_name %in% c('North Dakota', 'South Dakota', 'Nebraska', 'Kansas', 'Minnesota', 'Iowa', 
                        'Missouri', 'Wisconsin', 'Illinois', 'Indiana', 'Michigan', 'Ohio') ~ 'midwest',
      state_name %in% c('Kentucky', 'West Virginia', 'Virginia', 'Tennessee', 'North Carolina', 'Mississippi', 
                        'Alabama', 'Georgia', 'South Carolina', 'Florida', 'Louisiana', 'Arkansas') ~ 'southeast',
      state_name %in% c('Maine', 'Vermont', 'New Hampshire', 'Massachusetts', 'Rhode Island', 'Connecticut', 
                        'New York', 'New Jersey', 'Pennsylvania', 'Delaware', 'Maryland', 'District of Columbia') ~ 'northeast',
    )) %>% 
    filter(contrib_daily_structures_destroyed == 0 | is.na(contrib_daily_structures_destroyed)) %>% 
    filter(smokePM > 0) %>% 
    pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>%
    group_by(species) %>%
    dplyr::summarise(avg_no_struc_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')
  
  # nonsmoke_day_avg <- burned_struc_smoke_spec_df %>%
  #   filter(smokePM == 0) %>% 
  #   pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>%
  #   group_by(species) %>%
  #   dplyr::summarise(avg_nonsmoke_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')
  # 
  # all_avg_scenarios <- no_structures_burned_avg %>% 
  #   left_join(nonsmoke_day_avg)
    
    # Run a regression using your sample
    # clustered model
    sample_mod = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE,
                         K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                         S,  SE, SI, SO4, SR, TI, V,  ZN)
                       ~ smokePM*contrib_daily_structures_destroyed |
                         monitor_month + year, 
                       reg_df, cluster = 'site_id')
    # get coeffs
    sample_coeffs <- coeftable(sample_mod) %>% 
      #filter(sample != 'Full sample') %>% 
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
      dplyr::select(-id)) %>% 
      left_join(parameter_categories, by = 'species') %>%
      mutate(species_type = fct_relevel(species_type,
                                        c("Alkaline-earth metals", "Alkali metals",
                                          "Transition metals", "Metalloids", "Other metals",
                                          "Nonmetals",  "Halogens", "Organics"))) 
    
    # merge sample avg for each species and divide each species' betas by full sample avg for each species
    burned_struc_coeffs_normalized <- sample_coeffs %>% 
      left_join(no_structures_burned_smoke_conc_avg, by = 'species') %>% 
      mutate(pct_change_struct_vs_nostruct = 100*(Estimate/avg_no_struc_spec_conc), # how much a species has changed relative to its baseline
             pct_change_CI25 = 100*(CI25/avg_no_struc_spec_conc),
             pct_change_CI975 = 100*(CI975/avg_no_struc_spec_conc)) %>% 
      dplyr::select(species, coefficient, species_name, species_long, species_type, Estimate, CI25, CI975,
                    avg_no_struc_spec_conc, pct_change_struct_vs_nostruct, pct_change_CI25, pct_change_CI975)
    
    limits <- c("Magnesium (Mg)", "Chromium (Cr)", "Nickel (Ni)", 
                "Copper (Cu)", "Titanium (Ti)", "Zinc (Zn)", "Manganese (Mn)",
                "Arsenic (As)", "Lead (Pb)")
          

    
    
    plot_chems <- burned_struc_coeffs_normalized %>% 
      filter(species %in% c('MG','TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR')) %>% 
      filter(coefficient == 'smokePM:contrib_daily_structures_destroyed')
    
    # PLOT ALL TOGETHER NO FACETS
    burned_structures_reg_plot <- ggplot(plot_chems,
                                aes(x = species_long,
                                    y = pct_change_struct_vs_nostruct,
                                    color = species_type)) +
      geom_point(size=4, alpha = 0.6, stat = "identity", position = position_dodge(width = .8)) +
      geom_linerange(aes(ymin = (pct_change_CI25),
                         ymax = (pct_change_CI975)), stat = "identity", position = position_dodge(width = .8)) +
     scale_x_discrete(limits = rev(limits)) +
      scale_color_manual(values= c("#35274A", "#E58601", "#E2D200", "#46ACC8")) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
      labs(title = 'Effect of an additional structure burned \nper fire per day on species concentrations',
        y = '% change in concentration relative to no structures burning in smoke PM2.5',
           x = 'Species',
           color = 'Species Type',
      ) +
      theme_light() +
      coord_flip()+
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            title= element_text(size=12, face='bold'),
            legend.position = "top",
            axis.title.x = element_text(size=11, face = 'plain'),
            axis.title.y = element_text(size=11, face = 'plain')) +
      theme(axis.text.y = element_text(angle = 0, hjust = 1)) 
    burned_structures_reg_plot

    
    # save file
    ggsave(filename = 'Fig4A_selected_species_smoke_coeff_burned_figs.pdf',
           plot = burned_structures_reg_plot,
           path = file.path(results_fp, 'Fig4'),
           scale = 1,
           width = 8,
           height = 10,
           dpi = 320)
    ggsave(filename = 'Fig4A_selected_species_smoke_coeff_burned_figs.png',
           plot = burned_structures_reg_plot,
           path = file.path(results_fp, 'Fig4'),
           scale = 1,
           width = 8,
           height = 10,
           dpi = 320)
    
    
    
    # PLOT ALL TOGETHER NO FACETS
    all_coeffs_burned_structures_reg_plot <- ggplot(burned_struc_coeffs_normalized %>% 
                                           filter(coefficient == 'smokePM:contrib_daily_structures_destroyed'),
                                         aes(x = species_long,
                                             y = pct_change_struct_vs_nostruct,
                                             color = species_type)) +
      geom_point(size=4, alpha = 0.6, stat = "identity", position = position_dodge(width = .8)) +
      geom_linerange(aes(ymin = (pct_change_CI25),
                         ymax = (pct_change_CI975)), stat = "identity", position = position_dodge(width = .8)) +
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
      scale_color_manual(values= spec_pal) +
      #scale_color_manual(values= c("#35274A", "#E58601", "#E2D200", "#46ACC8")) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
      labs(title = 'Effect of an additional structure burned \nper fire per day on species concentrations',
           y = '% change in concentration relative to no structures burning in smoke PM2.5',
           x = 'Species',
           color = 'Species Type',
      ) +
      theme_light() +
      coord_flip()+
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            title= element_text(size=12, face='bold'),
            legend.position = "bottom",
            axis.title.x = element_text(size=11, face = 'plain'),
            axis.title.y = element_text(size=11, face = 'plain')) +
      theme(axis.text.y = element_text(angle = 0, hjust = 1)) 
    all_coeffs_burned_structures_reg_plot
    
    
    # save file
    ggsave(filename = 'SIFig5_all_species_smoke_coeff_burned_figs.pdf',
           plot = all_coeffs_burned_structures_reg_plot,
           path = file.path(results_fp, 'SI Figs'),
           scale = 1,
           width = 8,
           height = 10,
           dpi = 320)
    ggsave(filename = 'SIFig5_all_species_smoke_coeff_burned_figs.png',
           plot = all_coeffs_burned_structures_reg_plot,
           path = file.path(results_fp, 'Fig4'),
           scale = 1,
           width = 8,
           height = 10,
           dpi = 320)
    
    return(burned_struc_coeffs_normalized)
} 
    
# no_structures_burned_avg <- burned_struc_smoke_spec_df %>%
#   mutate(period = case_when(
#     year > 2006 & year <= 2010 ~ 'early sample',
#     year > 2010 & year <= 2015 ~ 'mid sample',
#     year > 2015 & year <= 2020 ~ 'late sample',
#   )) %>% 
#   filter(contrib_daily_structures_destroyed == 0 | is.na(contrib_daily_structures_destroyed)) %>% 
#   filter(smokePM > 0) %>% 
#   pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>%
#   group_by(species, period) %>%
#   mutate(avg_per_no_struc_spec_conc = mean(conc_val, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   group_by(species) %>%
#   mutate(full_samp_avg_no_struc_spec_conc = mean(conc_val, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   distinct(species, period, avg_per_no_struc_spec_conc, full_samp_avg_no_struc_spec_conc) %>% 
#   mutate(period = ifelse(is.na(period), 'Full sample', period))
# 
# time_periods <- reg_df %>% 
#   mutate(period = case_when(
#     year > 2006 & year <= 2010 ~ 'early sample',
#     year > 2010 & year <= 2015 ~ 'mid sample',
#     year > 2015 & year <= 2020 ~ 'late sample',
#   ))
# 
# period_mod = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE,
#                      K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
#                      S,  SE, SI, SO4, SR, TI, V,  ZN)
#                    ~ smokePM*contrib_daily_structures_destroyed |
#                      monitor_month + year, fsplit = ~period,
#                    time_periods, cluster = 'site_id')
# 
# # get coeffs
# period_coeffs <- coeftable(period_mod) %>%
#   filter(lhs %in% c('MG','TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR')) %>% 
#   filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
#   rename(pval = 'Pr(>|t|)',
#          se = 'Std. Error',
#          species = 'lhs', 
#          period = 'sample') %>%
#   mutate(pval = round(pval, digits = 3)) %>%
#   dplyr::select(-id, -sample.var) 
#   
#   CIs <- confint(period_mod) %>% 
#                    filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
#   dplyr::select(-id, -sample.var)  %>%
#       rename(species = 'lhs',
#              CI25 = '2.5 %',
#              CI975 = '97.5 %',
#              period = 'sample') %>% 
#   left_join(parameter_categories, by = 'species') %>%
#   mutate(species_type = fct_relevel(species_type,
#                                     c("Alkaline-earth metals", "Alkali metals",
#                                       "Transition metals", "Metalloids", "Other metals",
#                                       "Nonmetals",  "Halogens", "Organics"))) %>% 
#     filter(species %in% c('MG','TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR')) %>%
#     left_join(period_coeffs)
# 
# # merge sample avg for each species and divide each species' betas by full sample avg for each species
# burned_struc_coeffs_normalized <- CIs %>% 
#   left_join(no_structures_burned_avg) %>% 
#   mutate(per_est = Estimate/avg_per_no_struc_spec_conc, # how much a species has changed relative to its baseline
#          per_CI25 = CI25/avg_per_no_struc_spec_conc,
#          per_CI975 = CI975/avg_per_no_struc_spec_conc) %>% 
#   mutate(fs_est = Estimate/full_samp_avg_no_struc_spec_conc, # how much a species has changed relative to its baseline
#          fs_CI25 = CI25/full_samp_avg_no_struc_spec_conc,
#          fs_CI975 = CI975/full_samp_avg_no_struc_spec_conc) %>% 
#   dplyr::select(species, coefficient, species_name, species_long, species_type, period,Estimate, CI25, CI975,
#                 period_baseline = 'avg_per_no_struc_spec_conc', 
#                 fullsamp_baseline ='full_samp_avg_no_struc_spec_conc', per_est, per_CI25, per_CI975,
#                 fs_est, fs_CI25, fs_CI975) 
# 
# norm_coeffs <- burned_struc_coeffs_normalized %>% 
#   filter(period == 'Full sample') %>% 
#   dplyr::select(species, coefficient, species_name, species_long, species_type, period, 
#                 baseline = 'fullsamp_baseline', norm_est ='fs_est', CI25 ='fs_CI25', CI975 = 'fs_CI975') %>% 
#   bind_rows(burned_struc_coeffs_normalized %>% 
#             filter(period != 'Full sample') %>% 
#               dplyr::select(species, coefficient, species_name, species_long, species_type, period,
#                             baseline ='period_baseline', norm_est = 'per_est', CI25 ='per_CI25', CI975 ='per_CI975'))
#   
# 
# # limits <- c("Magnesium (Mg)", "Arsenic (As)", "Titanium (Ti)", "Lead (Pb)", 
# #             "Chromium (Cr)", "Nickel (Ni)", "Manganese (Mn)", 
# #             "Copper (Cu)", "Zinc (Zn)"
# # )
# 
# limits <- c("early sample", "mid sample", "late sample", "Full sample")
# 
# 
# avg_burned_structures <- mean(reg_df$contrib_daily_structures_destroyed, na.rm = TRUE)
# 
# # PLOT ALL TOGETHER NO FACETS
# burned_structures_reg_plot <- ggplot(norm_coeffs %>% 
#                                        filter(period != 'Full sample'),
#                                      aes(x = species_long,
#                                          y = 100*norm_est*avg_burned_structures,
#                                          color = period)) +
#   geom_point(size=4, alpha = 0.6, 
#              stat = "identity",
#              position = position_dodge(width = .8)) +
#   geom_linerange(aes(ymin = (CI25*100*avg_burned_structures),
#                      ymax = (CI975*100*avg_burned_structures)), stat = "identity", 
#                  position = position_dodge(width = .8)) +
#   #scale_x_discrete(limits = c('early sample', 'mid sample', 'late sample')) +
#   scale_color_manual(values= c('goldenrod', 'salmon', 'firebrick')) + #"#DC863B" "#FAD510" "#649373" "#1B5656" "#5A283E" "#F2300F"
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#   labs(title = 'Effect of average structures burned in fire/day on species concentrations',
#        y = '% change in concentration relative to no structures burning in smoke PM2.5',
#        x = 'Species',
#        color = 'Period') +
#   theme_light() +
#   facet_grid(.~species_long, ncol = 1) +
#   coord_flip()+
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         title= element_text(size=12, face='bold'),
#         legend.position = "top",
#         axis.title.x = element_text(size=11, face = 'plain'),
#         axis.title.y = element_text(size=11, face = 'plain')) +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1))
# burned_structures_reg_plot
# 
# 
# # save file
# ggsave(filename = 'Fig4A_selected_species_smoke_coeff_burned_figs.pdf',
#        plot = burned_structures_reg_plot,
#        path = file.path(results_fp, 'figures/Fig4'),
#        scale = 1,
#        width = 8,
#        height = 8,
#        dpi = 320)
# ggsave(filename = 'Fig4A_selected_species_smoke_coeff_burned_figs.png',
#        plot = burned_structures_reg_plot,
#        path = file.path(results_fp, 'figures/Fig4'),
#        scale = 1,
#        width = 8,
#        height = 8,
#        dpi = 320)



    # ------------------------------------------------------------------------------------
    # get response for CU/K and PB/K as per Boaggio et al
    # BOOTSTRAPPING RESPONSE CURVE
    # ------------------------------------------------------------------------------------
    # ratio_df <- reg_df %>%
    #   dplyr::select(ID:smokePM, PB, CU, K, monitor_month) %>%
    #   mutate(PB_CU_K_ratio = (PB+CU)/K) %>%
    #   mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>%
    #   filter(PB_CU_K_ratio>=0)
    # 
    # 
    # # PLOT RATIO OVER TIME
    # plot_ratio <- ggplot(ratio_df,
    #                              aes(x = log(contrib_daily_structures_destroyed),
    #                                  y = PB_CU_K_ratio)) +
    #   geom_point(size = 0.6) +
    #   geom_smooth() +
    #   theme_minimal() +
    #   theme(plot.title = element_text(size = 14, face = "bold"),
    #         axis.text = element_text(size = 9),
    #         axis.title = element_text(size = 12),
    #         axis.text.x = element_text(vjust = 0),
    #         axis.ticks.x = element_line(linetype = 1),
    #         axis.ticks.length=unit(-0.1, "cm"),
    #         panel.grid = element_blank(),
    #         panel.grid.major.y = element_line(color = "gray95")) +
    #   geom_hline(aes(yintercept = 0))
    # plot_ratio
    # 
    # # RUN REGRESSIOn
    # ratio_mod = feols(PB_CU_K_ratio ~ log(contrib_daily_structures_destroyed),
    #                   ratio_df)
    # # # get coeffs
    # ratio_coeffs <- broom::tidy(ratio_mod) %>%
    #   rename(pval = 'p.value',
    #          se = 'std.error',
    #          coefficient = 'term')

    # # BOOTSTRAPPING CONFIDENCE INTERVALS
    # # -----------------------
    # # bootstrap at the level of clustering
    # # stratified bootstrap function ----
    # set.seed(123) #set for replicable results
    # 
    # # parellelize
    # future::plan(multisession)
    # 
    # # set parameters
    # units_uq <- unique(ratio_df$monitor_month)
    # num_simulations <- 2
    # sample_units <- length(units_uq)
    # 
    # # set up empty dataframe to store
    # pred_results_df <- c(NA)
    # 
    # start_time <- Sys.time() # see how long it takes
    # print(start_time)
    # # Simulate the sampling distributions
    # for (i in 1:num_simulations) {
    #   # 1. Sample units from the population 
    #   sample_df <- ratio_df %>% 
    #     group_by(monitor_month) %>% # group by the clusters
    #     sample_n(sample_units, replace = TRUE) %>% 
    #     ungroup()
    #   
    #   # 2. Run a regression using your sample
    #   # sample_mod <- fixest::feols(c(PB_K_ratio, CU_K_ratio) ~
    #   #                             smokePM*contrib_daily_structures_destroyed | monitor_month + year,
    #   #                             sample_df, cluster = 'monitor_month')
    #   # # get coefficients frmo regression
    #   # sample_coeffs <- coeftable(sample_mod) %>% 
    #   #   rename(ratio = 'lhs')
    #   
    #   #ratio_list <- unique(sample_coeffs$ratio)
    #   # current_ratio <- ratio_list[1]
    #   
    #   # map over each species ratio
    #   # all_species_marg_slopes <- map_df(ratio_list, function(current_ratio) {
    #   ratio_mod = feols(PB_CU_K_ratio ~ smokePM*contrib_daily_structures_destroyed |
    #                       monitor_month + year,
    #                     ratio_df, cluster = 'monitor_month')
    #   # get coeffs
    #   ratio_coeffs <- broom::tidy(ratio_mod) %>% 
    #     rename(pval = 'p.value',
    #            se = 'std.error',
    #            coefficient = 'term') 
    # 
    #     # grab coefficients
    #     B1 <- ratio_coeffs %>% 
    #       filter(coefficient == 'smokePM') %>% 
    #       pull(estimate) 
    #     
    #     B2 <- ratio_coeffs %>% 
    #       filter(coefficient == 'contrib_daily_structures_destroyed') %>% 
    #       pull(estimate) 
    #     
    #     B3 <- ratio_coeffs %>% 
    #       filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
    #       pull(estimate) 
    #     
    #     # calculate the marginal effects of a structure burned over a given range
    #     # pred conc = sb_val(B2 + B3 * smoke_avg)
    #     # set to average smoke level
    #     #smoke_avg <- mean(ratio_df$smokePM, na.rm=T) # this is just the average smoke PM2.5 in each bin
    #     sb_vals <- seq(0, 2000, 1) # projected structures burned 
    #     results_df <- tibble()
    #     
    #     for (sb_val in sb_vals) {
    #       
    #       pred_ratio = B1 + B3*sb_val
    #       
    #       y_results <- tibble(pred_ratio = pred_ratio,
    #                           sb = sb_val)
    #       
    #       # Combine results for the current species
    #       results_df <- results_df %>%
    #         bind_rows(y_results) 
    #         
    #     } # end loop over each structure burned to predict ratio
    #     
    #       current_spec_results <- results_df %>%
    #         mutate(ratio = current_ratio) %>% 
    #         mutate(index = i)
    #       
    #   # save this bootstrap iteration results in a dataframe
    #   pred_results_df <- rbind(pred_results_df, current_spec_results) 
    #   
    #   } # end loop over simulation
    # 
    # # now calculate std errors
    # bootstrapped_df <- pred_results_df %>% 
    #   filter(!is.na(index)) %>% 
    #   group_by(sb) %>% 
    #   dplyr::summarise(est = mean(pred_ratio, na.rm = TRUE),
    #                    lb = quantile(pred_ratio, 0.025)[[1]],
    #                    ub = quantile(pred_ratio, 0.975)[[1]]) %>% 
    #   ungroup()
    # 
    # plot <- ggplot(bootstrapped_df,
    #                aes(x = sb, 
    #                    y = est), 
    #                    color = 'goldenrod', fill = 'goldenrod') +
    #   geom_line(alpha = .8, linewidth = .5)  +
    #   geom_ribbon(aes(x = sb, y = est, ymin= est - lb, ymax=est +ub), 
    #               alpha = 0.2, show.legend = FALSE, colour = NA) +
    #   labs(x = '# structures burned/day',
    #        y = 'Marginal effect of an additional structure burned') +
    #   theme_minimal() + 
    #   theme(panel.border = element_blank(),
    #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #         axis.line = element_line(colour = "black"),
    #         title= element_text(size=12, face='bold'),
    #         axis.title.x = element_text(size=11, face = 'plain'),
    #         axis.title.y = element_text(size=11, face = 'plain')) +
    #   geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5, color = 'grey60') 
    # plot
    # 
    #   
    # 
    


# ----------------------------------------------------------------
# STEP 1. PLOT HISTOGRAM OF BURNED STRUCTURES
# ----------------------------------------------------------------
# plot burned structures histogram
# bs_histogram <- ggplot(burned_struc_smoke_spec_df, aes(x = contrib_daily_structures_destroyed)) +
#   geom_histogram(bins = 100, color = 'darkred', fill = 'violetred4', alpha =.6) +
#   xlim(0, 2000) +
#   ylim(0, 100) +
#   xlab("# structures burned per day per fire") +
#   ylab("\n") +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.background = element_blank())
# bs_histogram



   

    # 
    # ## percentile of observed values
    # ecdf(result_df %>% filter(species=="Copper (CU)") %>% .$coef)(coef_df %>% filter(species=="Copper (CU)") %>% .$est)
    # ecdf(result_df %>% filter(species=="Lead (PB)") %>% .$coef)(coef_df %>% filter(species=="Lead (PB)") %>% .$est)
    # 
    # ## figures
    # p_final <- plot_grid(p_coef, p_perm_test, ncol=1)
    # save_plot("etc/coef_rand_test.pdf", plot=p_final, base_width=7, base_height=4, device=cairo_pdf)
    # 
    # 
    # # Other ----
    # 
    # 
    # ## NOT USED alternative method for merging
    # # globfire_spatial_joined_df <- pblapply(year_list, function(main_year_){
    # #     
    # #     ## spatial join doesnt also join on year so filter then spatial join
    # #     ## subset to relevant year
    # #     damaged_structure_temp_df <- damaged_structure_df %>% 
    # #         filter(year==main_year_)
    # #     
    # #     globfire_temp_df <- globfire_overlap_df %>% 
    # #         filter(year==main_year_) 
    # #     
    # #     ## spatial join and filter out unlikely matches starting with very tight bounds on date match
    # #     temp_globfire_structure_df <- globfire_temp_df %>% 
    # #         st_join(damaged_structure_temp_df %>% dplyr::select(-year), left=T) %>% 
    # #         filter(((IDate-3) <= start_date) & (start_date<=IDate+2)) %>% 
    # #         dplyr::select(Id, year, IDate, FDate, 
    # #                       Event_ID, irwinID, Incid_Name, Incid_Type, Ig_Date, start_date, 
    # #                       BurnBndAc, BurnBndLat, BurnBndLon,Pre_ID, Post_ID, Perim_ID,
    # #                       incident_number, incident_name, burn_area, 
    # #                       intersection_area, coverage_perc, structures_destroyed) %>% 
    # #         st_drop_geometry()
    # #     
    # #     return(temp_globfire_structure_df)        
    # #     
    # # }) %>% bind_rows()
    # # 
    # # ## merge unmatched records by name, state, and year (after removing spatial joined obs)
    # # globfire_text_joined_df <- globfire_overlap_df %>% 
    # #     filter(Id %not_in% c(NA, unique(globfire_spatial_joined_df$Id))) %>% 
    # #     mutate(incid_name_lower = str_to_lower(Incid_Name),
    # #            state=str_sub(Event_ID, start=1, end=2)) %>% 
    # #     filter(incid_name_lower %not_in% c(NA, "unnamed")) %>% 
    # #     left_join(damaged_structure_df %>% 
    # #                   st_drop_geometry() %>% 
    # #                   filter(incident_number %not_in% globfire_spatial_joined_df$incident_number), 
    # #               by=c("year", "state"="st", "incid_name_lower"="incident_name_lower")) %>% 
    # #     filter(!is.na(structures_destroyed)) %>% 
    # #     dplyr::select(Id, year, IDate, FDate, 
    # #                   Event_ID, irwinID, Incid_Name, Incid_Type, Ig_Date, start_date, 
    # #                   BurnBndAc, BurnBndLat, BurnBndLon,Pre_ID, Post_ID, Perim_ID,
    # #                   incident_number, incident_name, burn_area, 
    # #                   intersection_area, coverage_perc, structures_destroyed) %>% 
    # #     st_drop_geometry()
    # # 
    # # globfire_structure_joined_df <- bind_rows(globfire_spatial_joined_df, globfire_text_joined_df)
    # # 
    # # globfire_structure_joined_df$Id[duplicated(globfire_structure_joined_df$Id)]
    # 
    # # # Soil geochemistry exploration ----
    # # soil_geochem_df <- readr::read_csv("data/soils/soil_geochemistry_extract.csv") %>% 
    # #     mutate(lon = Longitude,
    # #            lat = Latitude) %>% 
    # #     st_as_sf(coords=c("Longitude","Latitude"), crs=4326, agr="constant") %>% 
    # #     st_transform("epsg:5070") %>% 
    # #     dplyr::select(Top5_Cr) %>% 
    # #     mutate(
    # #            # geometry = st_geometry(.) %>% 
    # #            #     st_combine(.) %>% 
    # #            #     st_voronoi(.) %>% 
    # #            #     st_collection_extract(.),
    # #            Top5_Cr=replace(Top5_Cr, Top5_Cr == "N.S.", NA),
    # #            Top5_Cr=replace(Top5_Cr, Top5_Cr == "<1", 0),
    # #            Top5_Cr=as.integer(Top5_Cr),
    # #            cr_capped=ifelse(Top5_Cr>220.67, 220,Top5_Cr)) %>% 
    # #     st_intersection(us_shape) # clip to US shape
    # # 
    # # US_BBOX <- st_bbox(us_shape)
    # # 
    # # ggplot() +
    # #     # map_theme + 
    # #     # geom_sf(data=us_shape) + 
    # #     # geom_sf(data=soil_geochem_df, size=0.3)
    # #     geom_sf(data=soil_geochem_df,aes(color=cr_capped)) + 
    # #     # geom_sf(data=soil_geochem_df,aes(fill=cr_capped, color=cr_capped)) + 
    # #     coord_sf(xlim = c(st_bbox(US_BBOX)$xmin, st_bbox(US_BBOX)$xmax), 
    # #              ylim = c(st_bbox(US_BBOX)$ymin, st_bbox(US_BBOX)$ymax))
    # 
    # # custom_trans <- function(){
    # #     trans <- function(x){ifelse(x<0,0,x^(1/4))}
    # #     inv <- function(x){x^4}
    # #     scales::trans_new('custom', trans, inv)
    # # }
    # # 
    # # ggplot(data=globfire_structure_joined_df,
    # #        aes(x=BurnBndAc, y=structures_destroyed)) + 
    # #     theme_classic() +
    # #     theme(
    # #         legend.position = c(0.2, 0.9),
    # #         legend.title=element_text(color="black", size=8),
    # #         text = element_text(color = "#22211d"),
    # #         plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    # #         panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    # #         legend.background = element_rect(fill = "#f5f5f2", color = NA)
    # #     ) +
    # #     geom_hex(bins=40) + 
    # #     geom_smooth(method='lm', se=FALSE, lwd=0.5, linetype=2, fullrange=T, color='red') + 
    # #     scale_x_continuous(trans='custom', labels=scales::comma,
    # #                        breaks=c(500,100000,250000,500000,1000000)) +
    # #     scale_y_continuous(trans='custom', labels=scales::comma) +
    # #     labs(y="# structures destroyed", x="Acres burned") + 
    # #     scale_fill_viridis_c(name="# observations", 
    # #                          breaks=c(1,5,10,15,30),
    # #                          guide = guide_legend(keyheight = unit(2.5, units = "mm"),
    # #                                               keywidth=unit(10, units = "mm"), 
    # #                                               label.position = "bottom", 
    # #                                               title.position = 'top', nrow=1)) + 
    # #     ggpubr::stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
    # #                      r.accuracy = 0.01,
    # #                      label.x = 0, label.y = 10, size = 4)
    # 
    
    
    
    
#     # bootstrap at the level of clustering
#     # stratified bootstrap function ----
#     set.seed(123) #set for replicable results
#     
#     # parellelize
#     # future::plan(multisession)
#     
#     # set parameters
#     units_uq <- unique(reg_df$monitor_month)
#     num_simulations <- 1
#     sample_units <- length(units_uq)
#     
#     # set up empty dataframe to store
#     pred_results_df <- NA
#     
#     start_time <- Sys.time() # see how long it takes
#     print(start_time)
#     
#     # Simulate the sampling distributions
#     for (i in 1:num_simulations) {
#       
#       # 1. Sample units from the population 
#       sample_df <- reg_df %>% 
#         mutate(burned_bins = case_when(
#           contrib_daily_structures_destroyed > 0 & contrib_daily_structures_destroyed <= 10 ~ 'low structures burned',
#           contrib_daily_structures_destroyed > 10 & contrib_daily_structures_destroyed <= 100 ~ 'medium structures burned',
#           contrib_daily_structures_destroyed > 100 ~ 'high structure burned'
#         )) %>% 
#         group_by(monitor_month) %>% # group by the clusters
#         sample_n(sample_units, replace = TRUE) %>% 
#         ungroup()
#     
#     
#     # 2. loop over each species to calculate the predicted concentration with this iteration of sample data
#     species <- unique(plot_chems$species)
#     # current_spec <- species[1]
#   
#     all_species_marg_slopes <- map_df(species, function(current_spec) {
#       
#       formula <- as.formula(
#         paste0(current_spec, " ~ smokePM*burned_bins |
#                          monitor_month + year"))
#       
#       sample_mod <- feols(fml = formula, sample_df, cluster = 'monitor_month')
#       # marg_slopes <- avg_slopes(sample_mod, by = TRUE, vcov = TRUE) %>% 
#       #   mutate(species = current_spec)
#       current_spec_coeffs <- coeftable(sample_mod) %>% 
#         as_data_frame() %>% 
#         mutate(coeff_id = row_number())
#   
#       # calculate the marginal effects:
#       #struc_burned_vec <- seq(0, 100, 1) # projected number of structures burned per day
#       smoke_vals <- seq(0, 100, 1) # projected smoke concentration
#       
#       B1 <- current_spec_coeffs %>% 
#         filter(coeff_id == 1) %>% 
#         pull(Estimate) 
#       
#       # grab coefficients
#       B2 <- current_spec_coeffs %>% 
#         filter(coeff_id == 2) %>% 
#         pull(Estimate) 
#       
#       B3 <- current_spec_coeffs %>% 
#         filter(coeff_id == 3) %>% 
#         pull(Estimate) 
#       
#       smoke_avg <- mean(reg_df$smokePM, na.rm= TRUE) # this is just the average smoke PM2.5 in each bin
#       #struct_burned_vals <- unique(reg_df$SB_avg) # this is the average structures burned in each bin
#       
#       results_df <- tibble()
#       
#       # Calculate y for each combination of x at average smoke level
#       
#           # calculate the marginal effects of a structure burning over a given range
#           # pred_conc = sb_val(B2 + B3*smoke_bin_avg)
#           # B1 = smokePM
#           # B2 = contrib_daily_structures_destroyed
#           # B3 = smokePM:contrib_daily_structures_destroyed
#           # sb_val will be the structures burned in a fire on a given day, range from 0-2100, 
#           # but support from 0-25, so predict over 0-25
#           # here "smoke_avg" represents the average smoke in each bin (low or high scenario)
#       
#       # Iterate over x_vec and smoke_vals to calculate y
#       #for (sb in struc_burned_vec) {
#          for (smoke_val in smoke_vals) {
#         # calculate marginal effect
#         marg_eff_structures_burned <- (B2 + B3 * smoke_avg)
#         # calculate predicted concentration at each level of structures burned
#           pred_conc <- sb * marg_eff_structures_burned
#           
#           y_results <- tibble(SB_val = sb,
#                               avg_smoke = smoke_avg,
#                               spec_conc = pred_conc,
#                               marg_eff = marg_eff_structures_burned,
#                               B1 =B1,
#                               B2 = B2,
#                               B3 =B3)
#                           
#           
#           # Combine results for the current species
#           results_df <- results_df %>%
#             bind_rows(y_results)
#         #}
#       }
#       
#       current_spec_results <- results_df %>% 
#         mutate(species = current_spec) #%>% 
#         # mutate(smoke_bins = case_when(
#         #   avg_smoke_level_in_bin <= 5 ~ 'low smoke',
#         #   avg_smoke_level_in_bin >= 25 ~ 'high smoke',
#         # ))
#     }) %>%  # end map over each species
#       mutate(index = i) # add index for iteration of bootstrap
#     
#     # 5. save this bootstrap iteration results in a dataframe
#     pred_results_df <- rbind(pred_results_df, all_species_conc_results)
# 
#   } # end bootstrap calculations loop
# 
#   
#   # After obtaining your bootstrapped estimates, do the following:
#   # 6. calculate the mean of the estimate and confidence intervals
#   boostrapped_CIs_df <- pred_results_df %>% 
#     filter(!is.na(SB_val)) %>% 
#     group_by(species, SB_val) %>% 
#     dplyr::summarise(est_conc = mean(spec_conc, na.rm = TRUE),
#                      lb = quantile(spec_conc, 0.025)[[1]],
#                      ub = quantile(spec_conc, 0.975)[[1]]) %>% 
#     ungroup() %>% 
#     left_join(parameter_categories, by = 'species')
#   
#   
#   # PLOT
#   # set bin order
#   # bin_order <- c("0-5", "5-10", "10-25","25-50", ">50" )
#   # bin_order <- c("low smoke", "high smoke")
#   
# 
#   plot <- ggplot(boostrapped_CIs_df,
#                  aes(x = SB_val, y = est_conc, color = species_type)) +
#     geom_line(alpha = .8, linewidth = .5)  +
#     geom_ribbon(aes(x = SB_val, y = est_conc, ymin=lb, ymax=ub, 
#                     fill = species_type), 
#                 alpha = 0.2, show.legend = FALSE, colour = NA) +
#     #scale_x_discrete(limits = rev(limits)) +
#     scale_color_manual(values = spec_pal) +
#     labs(x = '# structures burned/day',
#          y = 'Estimated species concentration in smoke') +
#     facet_wrap(~factor(species_long, levels = rev(limits)), scales = 'free_y', ncol = 4) +
#     theme_minimal() + 
#     theme(panel.border = element_blank(),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           title= element_text(size=12, face='bold'),
#           axis.title.x = element_text(size=11, face = 'plain'),
#           axis.title.y = element_text(size=11, face = 'plain')) +
#     geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5, color = 'grey60') 
#   plot
#   
#   
#   # plot precip histogram
#   ppt_histogram <- ggplot(reg_df, aes(x = contrib_daily_structures_destroyed)) +
#     geom_histogram(bins = 100, color = 'darkred', fill = 'violetred4', alpha =.6) +
#     xlim(0, 100) +
#     ylim(0, 1000) +
#     xlab("# structures burned per day per fire") +
#     ylab("\n") +
#     theme_minimal() +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.background = element_blank())
#   
#   
#   # save file
#   ggsave(filename = 'Fig3A_conc_smoke_response_2_burned_figs.pdf',
#          plot = plot,
#          path = file.path(results_fp, 'figures/Fig3'),
#          scale = 1,
#          width = 12,
#          height = 8,
#          dpi = 320)
#   
#   ggsave(filename = 'Fig3A_conc_smoke_response_2_burned_figs.png',
#          plot = plot,
#          path = file.path(results_fp, 'figures/Fig3'),
#          scale = 1,
#          width = 12,
#          height = 8,
#          dpi = 320)
#   
#   ggsave(
#     filename = 'Fig3B_conc_smoke_response_2_burned_figs.pdf',
#     plot = ppt_histogram,
#     path = file.path(results_fp, 'figures/Fig3'),
#     scale = 1,
#     width = 12,
#     height = 8,
#     dpi = 320)
#   
#   ggsave(
#     filename = 'Fig3B_conc_smoke_response_2_burned_figs.png',
#     plot = ppt_histogram,
#     path = file.path(results_fp, 'figures/Fig3'),
#     scale = 1,
#     width = 12,
#     height = 8,
#     dpi = 320)
#   
#   end_time <- Sys.time()
#   print(end_time)
#   future::plan(NULL)
#   
#   return(boostrapped_CIs_df)
# }
# 
# 
# # old code before bootstrapping:
# # run model at the grid cell level
# burned_struc_mod = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE,
#                            K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
#                            S,  SE, SI, SO4, SR, TI, V,  ZN)
#                          ~  smokePM + smokePM*contrib_daily_structures_destroyed |
#                            ID + month + year, reg_df, cluster = 'ID')
# 
# # Compute marginal effects for interaction term
# m <- margins(burned_struc_mod)
# margins_summary(burned_struc_mod)
# 
# # 
# 
# 
# 
# # get coeffs"
# coeffs2 <- coeftable(burned_struc_mod) %>%
#   rename(pval = 'Pr(>|t|)',
#          se = 'Std. Error',
#          species = 'lhs') %>%
#   mutate(pval = round(pval, digits = 3)) %>%
#   dplyr::select(-id)  %>%
#   left_join(parameter_categories, by = 'species') %>%
#   mutate(species_type = fct_relevel(species_type,
#                                     c("Alkaline-earth metals", "Alkali metals",
#                                       "Transition metals", "Metalloids", "Other metals",
#                                       "Nonmetals",  "Halogens", "Organics"))
#   ) 
# #
# species <- unique(coeffs$species)
# # current_spec <- species[1]
# 
# # loop over each species
# all_species_conc_results <- map_df(species, function(current_spec) {
#   
#   # create a loop to go through each species
#   current_coeff_df <- coeffs2 %>%
#     filter(species == current_spec)
#   
#   # calculate the marginal effects of a structure burning over a given range
#   # pred_conc = sb_val(B3 + B4*smoke_bin_avg)
#   # B1 = nonsmokePM
#   # B2 = smokePM
#   # B3 = contrib_daily_structures_destroyed
#   # B4 = smokePM:contrib_daily_structures_destroyed
#   # sb_val will be the structures burned in a fire on a given day, range from 0-2100, 
#   # but support from 0-25, so predict over 0-25
#   # here "smoke_bin_avg" represents the average smoke in each bin (low or high scenario)
#   
#   # implementation steps
#   smoke_vec <- seq(0, 100, 1)
#   
#   B2 <- current_coeff_df %>%
#     filter(coefficient == 'contrib_daily_structures_destroyed') %>%
#     pull(Estimate)
#   
#   B3 <- current_coeff_df %>%
#     filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>%
#     pull(Estimate)
#   
#   SB_mean <- mean(reg_df$contrib_daily_structures_destroyed,  na.rm = TRUE)
#   
#   results_df <- tibble()
#   # goal is to calculate the following:
#   #  y = x*(B2 + B3*smoke_level_avg_in_category)
#   # Calculate y for each combination of x and smoke_vals
#   # Iterate over x_vec and smoke_vals to calculate y
#   for (smoke_val in smoke_vec) {
#     # for (SB_val in SB_vals) {
#     partial_SB <- (B2 + B3 * SB_mean)
#     calc_conc = partial_SB*smoke_val
#     
#     y_results <- tibble(partial = partial_SB,
#                         spec_conc = calc_conc,
#                         smokePM_conc = smoke_val,
#                         SB = SB_mean
#     )
#     
#     # Combine results for the current species
#     results_df <- results_df %>%
#       bind_rows(y_results)
#   }
#   
#   
#   current_spec_results <- results_df %>%
#     mutate(species = current_spec) 
#   
# })
# 
# 
# # # --------------------------------------------------------------------------------
# # bin_order = c("low structures burned","high structures burned")
# 
# burned_structures_coeff_plot <- ggplot(all_species_conc_results,
#                                        aes(x = smokePM_conc, y = spec_conc)) +
#   geom_line(alpha = .8, linewidth = .5)  +
#   # geom_ribbon(aes(x = smokePM_conc, y = spec_conc,  
#   #                 fill = smoke_bins), 
#   #             alpha = 0.2, show.legend = FALSE, colour = NA) +
#   # scale_color_manual(limits = bin_order,
#   #                    values = c('goldenrod', 'chocolate', 'firebrick', 'darkorchid4')) +
#   # scale_fill_manual(limits = bin_order,
#   #                    values = c('goldenrod', 'chocolate', 'firebrick', 'darkorchid4')) +
#   # scale_color_manual(limits = bin_order,
#   #                    values = c('goldenrod', 'darkorchid4')) +
#   # scale_fill_manual(limits = bin_order,
#   #                   values = c('goldenrod', 'darkorchid4')) +
# labs(x = '# structures burned/day',
#      y = 'Estimated species concentration in smoke') +
#   facet_wrap(~species, scales = 'free_y', ncol = 4) +
#   theme_minimal() + 
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         title= element_text(size=12, face='bold'),
#         axis.title.x = element_text(size=11, face = 'plain'),
#         axis.title.y = element_text(size=11, face = 'plain')) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5, color = 'grey60') 
# burned_structures_coeff_plot
# # 
