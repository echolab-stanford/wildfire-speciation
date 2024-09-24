# # # Emma Krasovich Southworth, emmars@stanford.edu
# # # Last Updated: April 19, 2024
# # # Description: make time series of burned structures in data for SI

# loadd(c(globfire_structure_joined_df, burned_struc_smoke_spec_df, parameter_categories), cache = drake_cache)

make_SI_burned_structures_time_series_plot <- function(globfire_structure_joined_df, burned_struc_smoke_spec_df) {
  
  agg_df <- globfire_structure_joined_df %>% 
    ungroup() %>% 
    st_drop_geometry() %>% 
    distinct(year, incident_name, structures_destroyed) %>% 
    group_by(year) %>% 
    tally(structures_destroyed)
  
  # Plotting
  ts_burned_structures <- ggplot(agg_df, aes(x = year, y = n)) +
    geom_line(color = 'deeppink4') +
    labs(title = "Structures destroyed in wildfires per year",
         x = "Year",
         y = "# of structures destroyed") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"))
  ts_burned_structures
  
  ggsave(filename = 'SIFig_burned_structures_time_series.pdf',
         plot = ts_burned_structures,
         path = file.path(results_fp, 'SI Figs'),
         scale = 1,
         width = 8,
         height = 4,
         dpi = 320)
  
  ggsave(filename = 'SIFig_burned_structures_time_series.png',
         plot = ts_burned_structures,
         path = file.path(results_fp, 'SI Figs'),
         scale = 1,
         width = 8,
         height = 4,
         dpi = 320)
  
  
  
  # ----------------------------------------------------------------
  # STEP 1. RUN INTERACTION MODEL THAT LOOKS AT ASSOCIATED BW STRUCTURES + SPECIES CONC
  # ----------------------------------------------------------------
  # set up regression dataframe
  reg_df <- burned_struc_smoke_spec_df %>% 
    filter(!is.na(contrib_daily_structures_destroyed)) %>% 
    dplyr::select(ID, Dataset:site_id, contrib_daily_structures_destroyed,
                  smokePM:ZR, nonsmokePM_MF) %>% 
    filter_at(vars(AL, AS, BR, CA, CHL, CL, CR, CU,
                   EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
                   RB, S, SE, SI, SO4, SR, TI, V, ZN, ZR),
              any_vars(!is.na(.))) %>% 
    filter(contrib_daily_structures_destroyed > 0) %>% 
    mutate(monitor_month = paste0(site_id, "-", month)) #%>% 
  
  # calculate nonsmoke day baseline
  no_structures_burned_avg <- burned_struc_smoke_spec_df %>%
    filter(contrib_daily_structures_destroyed == 0 | is.na(contrib_daily_structures_destroyed)) %>% 
    filter(smokePM > 0) %>% 
    pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>%
    group_by(species) %>%
    dplyr::summarise(avg_no_struc_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')
  
  # Run a regression using your sample
  # clustered model
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
        dplyr::select(-id)) %>% 
    left_join(parameter_categories, by = 'species') %>%
    mutate(species_type = fct_relevel(species_type,
                                      c("Alkaline-earth metals", "Alkali metals",
                                        "Transition metals", "Metalloids", "Other metals",
                                        "Nonmetals",  "Halogens", "Organics"))) 
  
  # merge sample avg for each species and divide each species' betas by full sample avg for each species
  burned_struc_coeffs_normalized <- sample_coeffs %>% 
    left_join(no_structures_burned_avg, by = 'species') %>% 
    mutate(no_struc_smoke_est = Estimate/avg_no_struc_spec_conc, # how much a species has changed relative to its baseline
           no_struc_smoke_CI25 = CI25/avg_no_struc_spec_conc,
           no_struc_smoke_CI975 = CI975/avg_no_struc_spec_conc) %>% 
    dplyr::select(species, coefficient, species_name, species_long, species_type, Estimate, CI25, CI975,
                  avg_no_struc_spec_conc, no_struc_smoke_est, no_struc_smoke_CI25, no_struc_smoke_CI975)
  
  limits = c(
    # "Organics"
    "Organic Carbon (OC)", "Elemental Carbon (EC)",
    # "Halogens"
    "Bromine (Br)", "Chlorine (Cl)", "Chloride (Chl)", 
    #  "Nonmetals"
    "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
    # "Other metals"
    "Titanium (Ti)", "Aluminum (Al)", "Lead (Pb)",
    # "Metalloids"
    "Silicon (Si)", "Arsenic (As)",
    # "Transition metals"
    "Manganese (Mn)", "Zinc (Zn)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
    # "Alkali metals"
    "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
    # "Alkaline-earth metals"
    "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")
  
  plot_chems <- burned_struc_coeffs_normalized %>% 
    #filter(species %in% c('MG','TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR')) %>% 
    filter(coefficient == 'smokePM:contrib_daily_structures_destroyed')
  
  avg_burned_structures <- mean(reg_df$contrib_daily_structures_destroyed, na.rm = TRUE)
  
  # PLOT ALL TOGETHER NO FACETS
  all_chems_bs_reg_plot <- ggplot(plot_chems,
                                       aes(x = species_long,
                                           y = 100*no_struc_smoke_est,
                                           color = species_type)) +
    geom_point(size=4, alpha = 0.6, stat = "identity", position = position_dodge(width = .8)) +
    geom_linerange(aes(ymin = (no_struc_smoke_CI25*100),
                       ymax = (no_struc_smoke_CI975*100)), stat = "identity", position = position_dodge(width = .8)) +
    scale_x_discrete(limits = limits) +
    scale_color_manual(values= spec_pal) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(title = 'Effect of an additional structure burned per day on species concentrations',
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
  all_chems_bs_reg_plot
  
  
  # save file
  ggsave(filename = 'SIFig_all_species_smoke_coeff_burned_figs.pdf',
         plot = all_chems_bs_reg_plot,
         path = file.path(results_fp, 'figures/SI Figs'),
         scale = 1,
         width = 8,
         height = 10,
         dpi = 320)
  ggsave(filename = 'SIFig_all_species_smoke_coeff_burned_figs.png',
         plot = all_chems_bs_reg_plot,
         path = file.path(results_fp, 'figures/SI Figs'),
         scale = 1,
         width = 8,
         height = 10,
         dpi = 320)
  
  
  
  return(agg_df)
  
}


