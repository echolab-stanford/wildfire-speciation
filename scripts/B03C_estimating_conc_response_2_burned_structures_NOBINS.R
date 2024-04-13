# # Emma Krasovich Southworth, emmars@stanford.edu
# # Last Updated: April 2, 2024
# # Description: run regression on burned structures and smoke speciation + plot response curves at different levels of smoke

# loadd(c(burned_struc_smoke_spec_df, parameter_categories), cache = drake_cache)

# function
estimating_conc_response_2_burned_structures <- function(burned_struc_smoke_spec_df, parameter_categories) {
  
  # set up regression dataframe
  reg_df <- burned_struc_smoke_spec_df %>% 
    mutate(smoke_bins = case_when(
      smokePM >= 0 & smokePM < 10 ~ 'low smoke',
      smokePM >= 10 ~ 'high smoke',
    )) %>%
    group_by(smoke_bins) %>% 
    mutate(smoke_bin_avg = mean(smokePM, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(SB_bins = case_when(
      contrib_daily_structures_destroyed >= 0 & contrib_daily_structures_destroyed < 2 ~ '< 2 structures burned/day',
      contrib_daily_structures_destroyed >= 2 ~ '>2 structures burned/day',
    )) %>%
    group_by(SB_bins) %>% 
    mutate(SB_avg = mean(contrib_daily_structures_destroyed, na.rm = TRUE)) %>% 
    ungroup()
  
  # BOOTSTRAPPING CONFIDENCE INTERVALS
  # -----------------------
  # bootstrap at the level of clustering
  # stratified bootstrap function ----
  set.seed(123) #set for replicable results
  
  # parellelize
  future::plan(multisession)
  
  # set parameters
  units_uq <- unique(reg_df$ID)
  num_simulations <- 4
  sample_units <- length(units_uq)
  
  # set up empty dataframe to store
  ames_df <- c(NA)
  
  start_time <- Sys.time() # see how long it takes
  print(start_time)
  # Simulate the sampling distributions
  for (i in 1:num_simulations) {
    # 1. Sample units from the population 
    sample_df <- reg_df %>% 
      group_by(ID) %>% # group by the clusters
      sample_n(sample_units, replace = TRUE) %>% 
      ungroup()
    
    colnames(sample_df) <- tolower(colnames(sample_df))
    
    # 2. Run a regression using your sample
    species <- tolower(unique(sample_coeffs$species))
    # current_spec <- species[2]
    
    all_species_marg_slopes <- map_df(species, function(current_spec) {
      
      formula <- as.formula(
        paste0(current_spec, "~ smokepm + contrib_daily_structures_destroyed + contrib_daily_structures_destroyed*smokepm | id + month + year"))
      
      sample_mod <- feols(fml = formula, sample_df, cluster = 'id')
      marg_slopes <- avg_slopes(sample_mod, by = TRUE, vcov = TRUE) %>% 
        mutate(species = current_spec)
      
    }) %>% 
      bind_rows() %>% 
      mutate(index = i)
 
    # 5. save this bootstrap iteration results in a dataframe
    ames_df <- rbind(ames_df, all_species_marg_slopes) 
    
  }   
    
  # now calculate std errors
  bootstrapped_ames_df <- ames_df %>% 
    filter(!is.na(index)) %>% 
    group_by(species, term) %>% 
    dplyr::summarise(est = mean(estimate, na.rm = TRUE),
                     lb = quantile(estimate, 0.025)[[1]],
                     ub = quantile(estimate, 0.975)[[1]],
                     se = mean(std.error, na.rm=TRUE),
                     pval = mean(p.value, na.rm= TRUE)) %>% 
    ungroup() %>% 
    mutate(species = toupper(species)) %>% 
    left_join(parameter_categories, by = 'species') %>% 
    filter(term == 'contrib_daily_structures_destroyed') %>% 
    mutate(sig = case_when(
      pval > .05 ~ 'not significant',
      pval <= .05 ~ 'significant',
    ))
  
  # Create the plot
 plot<-  ggplot(bootstrapped_ames_df %>% 
                  filter(term == 'contrib_daily_structures_destroyed'), 
                aes(x = species, y = est)) +
    geom_point(size = 3) +
   geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.1) +
    labs(title = "Avg Marginal Effect of structures burned on specices concentration", x = "species", y = "marginal effect") + 
   theme_minimal() + 
   # facet_wrap(~species, scales = 'free_y') +
   geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5, color = 'grey60') 
 plot

  
    # 3. get coeffs
    sample_coeffs <- coeftable(sample_mod) %>%
      rename(pval = 'Pr(>|t|)',
             se = 'Std. Error',
             species = 'lhs') %>%
      mutate(pval = round(pval, digits = 3)) %>%
      dplyr::select(-id)  %>%
      left_join(parameter_categories, by = 'species') %>%
      mutate(species_type = fct_relevel(species_type,
                                        c("Alkaline-earth metals", "Alkali metals",
                                          "Transition metals", "Metalloids", "Other metals",
                                          "Nonmetals",  "Halogens", "Organics"))
      ) %>% 
      mutate(index = i) 
    
    # 4. loop over each species to calculate the predicted concentration with this iteration of sample data
    species <- unique(sample_coeffs$species)
    # current_spec <- species[2]
    
    # loop over each species
    all_species_conc_results <- map_df(species, function(current_spec) {
      
      # create a loop to go through each species
      current_coeff_df <- sample_coeffs %>% 
        filter(species == current_spec)
      
      # calculate the marginal effects:
      smoke_vals <- seq(0, 100, 1) # projected smoke
      #x_vec <- seq(0, 100, 1) # projected smoke concentration
      
      # grab coefficients
      B1 <- current_coeff_df %>% 
        filter(coefficient == 'smokePM') %>% 
        pull(Estimate) 
      
      B2 <- current_coeff_df %>% 
        filter(coefficient == 'contrib_daily_structures_destroyed') %>% 
        pull(Estimate) 
      
      B3 <- current_coeff_df %>% 
        filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
        pull(Estimate) 
      
      smoke_avg <- mean(reg_df$smokePM, na.rm=T) # this is just the average smoke PM2.5 in each bin
      struct_burned_avg <- mean(reg_df$contrib_daily_structures_destroyed, na.rm=T) # this is the average structures burned in each bin
   
      partial_eff_smoke_on_conc <- (B1 + B3 * struct_burned_avg)
      partial_eff_sb_on_conc <- (B2 + B3 * smoke_avg)
          
    
     
 
  
  # After obtaining your bootstrapped estimates, do the following:
  # 6. calculate the mean of the estimate and confidence intervals
  boostrapped_CIs_df <- pred_results_df %>% 
    filter(!is.na(index)) %>% 
    group_by(species, smoke_bins, SB_val) %>% 
    dplyr::summarise(est_conc = mean(spec_conc, na.rm = TRUE),
                     lb = quantile(spec_conc, 0.025)[[1]],
                     ub = quantile(spec_conc, 0.975)[[1]]) %>% 
    ungroup() %>% 
    left_join(parameter_categories, by = 'species')
  
  
  # PLOT
  # set bin order
  # bin_order <- c("0-5", "5-10", "10-25","25-50", ">50" )
  bin_order <- c("low smoke", "high smoke")
  
  # set order of plots:
  limits = c(
    # "Organics"
    "Organic Carbon (OC)", "Elemental Carbon (EC)",
    # "Halogens"
    "Bromine (Br)", "Chlorine (Cl)", "Chloride (Chl)",
    #  "Nonmetals"
    "Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Phosphorus (P)", "Selenium (Se)",
    # "Other metals"
    "Titanium (Ti)", "Aluminum (Al)", "Lead (Pb)",
    # "Metalloids"
    "Silicon (Si)", "Arsenic (As)",
    # "Transition metals"
    "Zinc (Zn)", "Manganese (Mn)",  "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
    # "Alkali metals"
    "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
    # "Alkaline-earth metals"
    "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)"
  )
  
  plot <- ggplot(boostrapped_CIs_df,
                 aes(x = SB_val, y = est_conc, group = smoke_bins, color = smoke_bins)) +
    geom_line(alpha = .8, linewidth = .5)  +
    geom_ribbon(aes(x = SB_val, y = est_conc, ymin=lb, ymax=ub, 
                    fill = smoke_bins), 
                alpha = 0.2, show.legend = FALSE, colour = NA) +
    scale_color_manual(limits = bin_order,
                       values = c('goldenrod', 'darkorchid4')) +
    scale_fill_manual(limits = bin_order,
                      values = c('goldenrod', 'darkorchid4')) +
    labs(x = '# structures burned/day',
         y = 'Estimated species concentration in smoke') +
    facet_wrap(~factor(species_long, levels = rev(limits)), scales = 'free_y', ncol = 4) +
    theme_minimal() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5, color = 'grey60') 
  plot
  
  
  # plot precip histogram
  ppt_histogram <- ggplot(reg_df, aes(x = contrib_daily_structures_destroyed)) +
    geom_histogram(bins = 100, color = 'darkred', fill = 'violetred4', alpha =.6) +
    xlim(0, 100) +
    ylim(0, 1000) +
    xlab("# structures burned per day per fire") +
    ylab("\n") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank())
  
  
  # save file
  ggsave(filename = 'Fig3A_conc_smoke_response_2_burned_figs.pdf',
         plot = plot,
         path = file.path(results_fp, 'figures/Fig3'),
         scale = 1,
         width = 12,
         height = 8,
         dpi = 320)
  
  ggsave(filename = 'Fig3A_conc_smoke_response_2_burned_figs.png',
         plot = plot,
         path = file.path(results_fp, 'figures/Fig3'),
         scale = 1,
         width = 12,
         height = 8,
         dpi = 320)
  
  ggsave(
    filename = 'Fig3B_conc_smoke_response_2_burned_figs.pdf',
    plot = ppt_histogram,
    path = file.path(results_fp, 'figures/Fig3'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
  ggsave(
    filename = 'Fig3B_conc_smoke_response_2_burned_figs.png',
    plot = ppt_histogram,
    path = file.path(results_fp, 'figures/Fig3'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
  end_time <- Sys.time()
  print(end_time)
  future::plan(NULL)
  
  return(boostrapped_CIs_df)
}

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
#   # need to calculate the marginal effects of a structure burning
#   # y = x(B3 + B2*smoke)
#   # B1 = smokePM
#   # B2 = contrib_daily_structures_destroyed
#   # B3 = smokePM:contrib_daily_structures_destroyed
#   # x will be the structures burned in a fire on a given day, range from 0-2100, but support from 0-25
#   # here "smoke" represents the average smoke in each bin
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
# 
# # running regressions ----
# # fixest::feols(conc_PB~smokePM_pred_coded + smokePM_pred_coded*contrib_daily_structures_destroyed | ID + year + month, data=final_df_subset, cluster="State")
# # fixest::feols(conc_CU~smokePM_pred_coded + smokePM_pred_coded*contrib_daily_structures_destroyed | ID + year + month, data=final_df_subset, cluster="State")
# 
# 
# # permutation/ randomization inference ----
# ## sampling within state, across years 
# withinAcrossResampler <- function(data, within_condition, across_condition, interest_var){
#   
#   data_dt <- data %>% as.data.table()
#   
#   # iterate through within condition
#   outer_dt <- parallel::mclapply(data[[within_condition]] %>% unique(), function(within_unit){
#     
#     # get list of possible values of across condition
#     across_condition_list <-  data_dt[get(within_condition)==within_unit][[across_condition]] %>% 
#       unique()
#     
#     across_condition_dt <- parallel::mclapply(across_condition_list, function(across_unit){
#       
#       # number to sample needs to match how many obs we are trying to replace
#       inner_dt <- data_dt[get(within_condition)==within_unit & 
#                             get(across_condition)==across_unit] 
#       
#       num_samples <- nrow(inner_dt)
#       
#       # filter for within unit but not across unit, then sample without replacement
#       # from set of possible options
#       possible_new_values <- data_dt[get(within_condition)==within_unit & 
#                                        get(across_condition)!=across_unit
#       ][[interest_var]] 
#       
#       # get num_samples from possible_new_values
#       new_values <- possible_new_values[sample.int(length(possible_new_values), size=num_samples)]
#       
#       inner_dt[[stringr::str_glue('{interest_var}_new')]] <- new_values
#       
#       return(inner_dt)
#     }) %>% bind_rows()
#     
#     return(across_condition_dt)
#   }) %>% bind_rows()
#   
#   return(outer_dt)
# }
# 
# 
# permTest <- function(i, outcome, rhs, data, treatment, interest, cluster_err, across_condition, sample_group=NULL){
#   
#   ## sample and create new data with shuffled treatment
#   
#   if(!is.null(sample_group)){
#     new_data <- withinAcrossResampler(data, within_condition=sample_group, across_condition=across_condition, 
#                                       interest_var=treatment)
#     
#   } else{
#     new_data <- data %>% 
#       mutate(!!paste0(treatment,"_new") := sample(get(c(treatment)), size=n())) %>% 
#       ungroup()
#   }
#   
#   
#   fixest_df <- lapply(outcome, function(x){
#     fmla <- paste0(x,"~",rhs)
#     fixe_est <- fixest::feols(as.formula(stringr::str_replace_all(fmla, treatment, paste0(treatment,"_new"))),
#                               data=new_data)
#     
#     ## no clustering SEs here because we are just interested in the coef estimates
#     data.frame("coef"=c(summary(fixe_est)$coefficients[[paste0(interest,"_new")]]),
#                "species"=c(as.character(fixe_est$fml)[2]))
#   }) %>%
#     bind_rows()
#   
#   return(fixest_df)
# }
# 
# TREATMENT_VAR <- "contrib_daily_structures_destroyed"
# INTEREST_VAR <- "smokePM_pred_coded:contrib_daily_structures_destroyed"
# LINE_COLOR <- "blue"
# DEFAULT_COLOR <- "gray60"
# 
# mod_pb <- fixest::feols(conc_PB~smokePM_pred_coded + smokePM_pred_coded*contrib_daily_structures_destroyed | ID + year + month, data=final_df_subset, cluster="State")
# mod_cu <- fixest::feols(conc_CU~smokePM_pred_coded + smokePM_pred_coded*contrib_daily_structures_destroyed | ID + year + month, data=final_df_subset, cluster="State")
# 
# coef_df <- data.frame("species"=c("Lead (PB)", "Copper (CU)"),
#                       "est"=c(mod_pb$coeftable[INTEREST_VAR,'Estimate'],
#                               mod_cu$coeftable[INTEREST_VAR,'Estimate']),
#                       "se"=c(mod_pb$coeftable[INTEREST_VAR,'Std. Error'],
#                              mod_cu$coeftable[INTEREST_VAR,'Std. Error']))
# 
# p_coef <- ggplot(data=coef_df, aes(color=species)) + 
#   theme_minimal() + 
#   theme(
#     legend.position="none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_blank(),
#     axis.title.y = ggtext::element_markdown()
#   ) + 
#   geom_hline(yintercept = 0, colour = DEFAULT_COLOR, lty = 2, size=0.25) + 
#   geom_linerange(aes(x=species, ymin=est-qnorm(0.975)*se, ymax=est+qnorm(0.975)*se),lwd=0.7,position=position_dodge2(width=.25)) + 
#   geom_point(aes(x=species, y=est), position=position_dodge2(width=.25), shape = 21, fill = "WHITE", lwd=2) + 
#   labs(y="&Delta; &mu;g/m<sup>3</sup> per additional<br>structure destroyed", x="", color="Species") + 
#   facet_grid(~species, scales="free_x", space="free_x") + 
#   scale_color_manual(values=c("#d9a574", "#858583"))
# 
# 
# start_time <- Sys.time()
# result_df <- pbmclapply(1:1000, permTest, 
#                         outcome=c("conc_PB", "conc_CU"),
#                         rhs="smokePM_pred_coded + smokePM_pred_coded*contrib_daily_structures_destroyed | ID + year + month", 
#                         data=final_df_subset, 
#                         treatment=TREATMENT_VAR,
#                         interest=INTEREST_VAR,
#                         across_condition="year",
#                         sample_group="State") %>%
#   bind_rows()
# 
# result_df <- result_df %>% 
#   mutate(species=fifelse(species == "conc_PB","Lead (PB)", "Copper (CU)"))
# 
# end_time <- Sys.time()
# 
# p_perm_test <- ggplot(data=result_df, aes(color=species)) + 
#   theme_minimal() + 
#   theme(
#     legend.position="none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title.y = element_markdown()
#   ) + 
#   geom_histogram(aes(x=coef, fill=species), bins=100, alpha=0.3, position="identity", size=0.25) + 
#   geom_vline(data=coef_df, aes(xintercept=est), lty = 2, color=LINE_COLOR, size=0.25) +
#   coord_cartesian(xlim=c(-1e-6,1e-6), expand=F) +
#   facet_grid(species~., scales="free") + 
#   labs(x="Estimated coefficient per additional structure destroyed", 
#        y="Frequency<br>", fill="", color="") +
#   scale_fill_manual(values=c("#d9a574", "#858583")) + 
#   scale_color_manual(values=c("#d9a574", "#858583"))
# 
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
