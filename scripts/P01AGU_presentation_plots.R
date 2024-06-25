# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: December 7, 2023
# Description: make figures for AGU presentation


# FULL MODEL COEFF + SUBSET TO LARGEST EFFECT SIZES -----------------------------------------------------
loadd(c(parameter_categories, spec_pal, 
        full_samp_PMcoeffs_normalized, regionalPMcoeffs_normalized, regoin_pal), 
      cache = drake::drake_cache(".drake"))

limits <- c(
  # "Lead (Pb)", 
  "Phosphorus (P)", 
  "Potassium (K)",
  "Elemental Carbon (EC)", "Organic Carbon (OC)")
 
  # plot coefficients for speciation at the avg monitor which tells us how much 
  # of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
  # --------------------------------------------------------------------------------
  # plot percent change for all species
  # --------------------------------------------------------------------------------
  pct_change_samp_reg_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                       filter(species %in% c('EC', 'OC', 'K', 'P')), 
                                     aes(x = species_long,
                                         y = 100*norm_est, 
                                         color=species_type, 
                                     )) +
    geom_point(size=4, alpha = 0.6, stat = "identity") +
    geom_linerange(aes(ymin = (100*norm_CI25), 
                       ymax = (100*norm_CI975)), stat = "identity") +
    scale_x_discrete(limits = rev(limits)) +
     #spec_pal = c("#35274A", , "#E58601", "#E2D200", , "#0B775E","#CAB38C", )
    scale_color_manual(values= c("#B40F20","#E58601","#46ACC8", "#899DA4"),
                       name = "Species type") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(y = expression(paste('% change relative to average nonsmoke day concentration')))+
         #x = 'Species') +
         #color = 'Species Category',
         #title = expression(paste("% change in concentration for 1 ug/", m^3, " increase in smoke ", "PM"["2.5"]))) + 
    theme_light() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=16, face='bold'),
          axis.title.x = element_text(size=16, face = 'bold'),
          axis.title.y = element_text(size=16, face = 'bold'),
          axis.text.x = element_text(size = 16),  # Adjust the size of x-axis tick labels
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.position = c(0.02, 0.98),  # Adjust the position of the legend
          legend.justification = c(0, 1)) +  # Adjust the justification of the legend) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis text by 90 degrees
    guides(color = 'none')
  pct_change_samp_reg_plot
  
  # save file
  ggsave(
    filename = 'Pct_change_PM_P_K_OC_EC.png',
    plot = pct_change_samp_reg_plot,
    path = file.path(results_fp, 'figures/Presentation Figs'),
    scale = 1,
    width = 8,
    height = 8,
    dpi = 320) 
  # --------------------------------------------------------------------------------
  
  pct_change_samp_full_plot <- ggplot(full_samp_PMcoeffs_normalized, 
                                     aes(x = species_long,
                                         y = 100*norm_est, 
                                         color=species_type, 
                                     )) +
    geom_point(size=4, alpha = 0.6, stat = "identity") +
    geom_linerange(aes(ymin = (100*norm_CI25), 
                       ymax = (100*norm_CI975)), stat = "identity") +
    scale_x_discrete(limits = rev(c(
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
    ))) +
    scale_color_manual(values= spec_pal, name = 'Species Type') +
    # coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(y = expression(paste('% change relative to average nonsmoke day concentration')),
         x = 'Species',
         #color = 'Species Category',
         title = expression(paste("% change in concentration for 1 ug/", m^3, " increase in smoke ", "PM"["2.5"]))) + 
    theme_light() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=16, face='bold'),
          axis.title.x = element_text(size=16, face = 'bold'),
          axis.title.y = element_text(size=16, face = 'bold'),
          axis.text.x = element_text(size = 16),  # Adjust the size of x-axis tick labels
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.position = c(0.02, 0.98),  # Adjust the position of the legend
          legend.justification = c(0, 1)) +  # Adjust the justification of the legend) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  pct_change_samp_full_plot
  
  # save file
  ggsave(
    filename = 'Pct_change_PM_spec_smokeMF_ALL_chems.png',
    plot = pct_change_samp_full_plot,
    path = file.path(results_fp, 'figures/presentations'),
    scale = 1,
    width = 18,
    height = 11,
    dpi = 320) 
  

  # --------------------------------------------------------------------------------
# REGIONAL MODEL COEFF + SUBSET TO LARGEST EFFECT SIZES
  # --------------------------------------------------------------------------------
  
  # PLOT ALL TOGETHER NO FACETS
  regional_reg_plot <- ggplot(regionalPMcoeffs_normalized %>% 
                                filter(species %in% c('EC', 'OC', 'K', 'CU', 'PB')),
                              aes(x = species_long,
                                  y = 100*norm_est,
                                  color = region)) +
    geom_point(size=3, alpha = 0.7, stat = "identity", position = position_dodge(width = .4)) +
    geom_linerange(aes(ymin = (100*norm_CI25),
                       ymax = (100*norm_CI975)), stat = "identity", position = position_dodge(width = .4)) +
    scale_x_discrete(limits = limits) +
    scale_color_manual(values=region_pal) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    labs(y = '% change relative to regional baseline on nonsmoke days',
         x = 'Chemical Species',
         color = 'Region',
    ) +
    theme_light() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          title= element_text(size=12, face='bold'),
          axis.title.x = element_text(size=11, face = 'plain'),
          axis.title.y = element_text(size=11, face = 'plain')) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') #+
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis text by 90 degrees
  regional_reg_plot
  
  # save file
  ggsave(
    filename = 'REGIONAL_pct_change_PM_spec_smokeMF_3_chems.png',
    plot = regional_reg_plot,
    path = file.path(results_fp, 'figures/presentations'),
    scale = 1,
    width = 8,
    height = 8,
    dpi = 320) 



# TIME SERIES + SUBSET TO LARGEST EFFECT SIZES
  
  