# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species


# loadd(c(parameter_categories, pm_pal, clean_PMspec_df), cache = drake::drake_cache(".drake"))

# function
create_attributable_frac_time_series <- function(clean_PMspec_df, pm_pal, parameter_categories) {

  # select the vars that are needed for the regression
  reg_df <- clean_PMspec_df %>% 
    dplyr::select(Dataset, state_name, region, year, month, Date, 
                  monitor_month, smoke_day, site_id, MF_adj, smokePM, 
                  nonsmokePM_MF, AL:ZN) 
  
################################################################################
# How much is concentration due to smoke changing over time?
################################################################################
  # -----------------------------------------------------------------------------
  # RUN REGRESSION FOR SPECIES
  #  A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
  # -----------------------------------------------------------------------------
  full_sampPM_reg_MF = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                              K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                              S,  SE, SI, SO4, SR, TI, V, ZN)
                            ~ smokePM + nonsmokePM_MF | 
                              monitor_month + year, reg_df, cluster = 'site_id')
  
  ## calculate 95 CI%
  CIs <- confint(
    full_sampPM_reg_MF) %>% 
    rename(species = 'lhs',
           pm_type = 'coefficient',
           CI25 = '2.5 %',
           CI975 = '97.5 %') %>% 
    mutate(measure = 'RCFM') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  betas <- coeftable(full_sampPM_reg_MF) %>% 
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           pm_type = 'coefficient') %>% 
    dplyr::select(-id, -pval, -`t value`)  %>% 
    left_join(parameter_categories, by = 'species') %>% 
    mutate(measure = 'MF') %>% 
    filter(pm_type == 'smokePM')
  
# predict over our sample
predicted_spec_smoke_conc <- reg_df %>% 
  dplyr::select(Dataset, state_name, year, month, Date, 
                smoke_day, site_id, MF_adj, smokePM, 
                nonsmokePM_MF, AL:ZN) %>% 
  dplyr::select(-SOIL) %>% 
  pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>% 
  filter(!is.na(conc_val)) %>% 
  left_join(betas, by = 'species') %>% 
  mutate(pred_spec_nonsmoke = Estimate*nonsmokePM_MF,
         pred_spec_smoke = Estimate*smokePM) %>% 
  mutate(pred_spec_allPM = pred_spec_smoke+pred_spec_nonsmoke) %>% 
  distinct() 

  # calculate avg attributable fraction of species concentration to total concentration
  avg_conc_site_month_yr <- predicted_spec_smoke_conc %>%
    # get avg concentration across all sites for a given month and year
    group_by(species, year, month, site_id) %>% 
    dplyr::summarise(avg_smoke_conc = mean(pred_spec_smoke, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(pred_spec_nonsmoke, na.rm = TRUE),
                     avg_tot_conc = mean(pred_spec_allPM, na.rm = TRUE)) %>% 
    ungroup() %>% 
    # get avg concentration across all sites for a given month and year
    group_by(species, year, month) %>% 
    dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                     avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>% 
    ungroup() 
  
  
    # get monthly avg conc for the observed "total concentration" of a species
   obs_improve_monthly_mean <- reg_df %>% 
     dplyr::select(Dataset, year, month, Date, 
                  smoke_day, site_id, MF_adj, smokePM, 
                  nonsmokePM_MF, AL:ZN) %>% 
     dplyr::select(-SOIL) %>% 
     pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>% 
     group_by(Dataset, species, year, month, site_id) %>% 
     dplyr::summarise(avg_obs_tot_conc = mean(conc_val, na.rm = TRUE)) %>% 
     ungroup() %>% 
     group_by(Dataset, species, year, month) %>% 
     dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>% 
     ungroup() 
     
   # now merge back with predicted ata
   avg_attributable_fraction_vs_obs <- avg_conc_site_month_yr %>% 
     left_join(obs_improve_monthly_mean %>% 
                 pivot_wider(names_from = 'Dataset', values_from = 'avg_obs_tot_conc') %>% 
                 rename(avg_obs_tot_conc = 'IMPROVE'), 
               by = c('species', 'year', 'month'))
   
    # prep df for plotting
    plot_df <- avg_attributable_fraction_vs_obs %>% 
      dplyr::select(-avg_obs_tot_conc, -avg_nonsmoke_conc, -avg_tot_conc, -CSN) %>% 
      mutate(label = 'Attributable concentration to smoke PM2.5') %>% 
      rename(conc = 'avg_smoke_conc') %>% 
      bind_rows(avg_attributable_fraction_vs_obs %>% 
                  dplyr::select(conc = 'avg_obs_tot_conc', year, month, species) %>% 
                  mutate(label = 'Observed Concetration in Total PM2.5')) %>% 
      mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
      left_join(parameter_categories, by = 'species') %>% 
      mutate(label = fct_relevel(label, c("Observed Concetration in Total PM2.5", 
                                          "Attributable concentration to smoke PM2.5"))) 
    
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
    datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")
    
    
# PLOTTING
    avg_mon_pred_spec_plot <- ggplot(plot_df) +
      geom_density(aes(x = mon_yr,
                       y = conc,
                       fill = label,
                       color = label), 
                      stat = "identity", alpha = 0.8) +
      labs(title = expression("Predicted vs IMPROVE observed species concentrations"),
           y = expression("Concentration (\u00b5g/m"^3*")"),
           x = "",
      ) +
      scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
      geom_hline(aes(yintercept = 0)) +
      theme_minimal() +
      scale_color_manual(values = pm_pal) +
      scale_fill_manual(values = pm_pal) +
      # facet_wrap(~species_long, scales = 'free_y', ncol = 4) +
      facet_wrap(~factor(species_long, levels = rev(limits)), scales = 'free_y', ncol = 4) +
      theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
            axis.text = element_text(),
            axis.title = element_text(),
            legend.text = element_text(),
            legend.title = element_blank(),
            axis.ticks.x = element_line(linetype = 1),
            axis.ticks.length=unit(-0.2, "cm"),
            panel.grid = element_blank(),
            legend.position = "bottom",
            axis.line.y = element_line(color = "grey10"))  # Add y-axis line   # Add y-axis ticks
    avg_mon_pred_spec_plot
    
 
  # save file
  ggsave(
    filename = 'Fig4_smoke_attributable_fraction_monthly_trend_full_model.pdf',
    plot = avg_mon_pred_spec_plot,
    path = file.path(results_fp, 'figures/Fig4'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
    return(predicted_spec_smoke_conc)
}
