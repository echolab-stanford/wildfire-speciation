# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: create a bar plot for each region's breakdown for 1 ug/m3 of pm2.5

# loadd(clean_pm_spec_df, cache = drake_cache)

# function
create_regional_bar_plot_w_map <- function(clean_pm_spec_df) {
  
  # set up 
  selected_spec_df <- clean_pm_spec_df %>% 
    dplyr::select(Dataset:month, V, PB, FE, 
                  CU, CR, MN, ZN, NI,
                  EC, OC, S, smokePM, 
                  nonsmokePM, totPM2.5, monitor_month) %>% 
    # only keep rows that have all values so that its the same selection of sites for the analysis
    na.omit() %>% 
    pivot_longer(cols = c(V, PB, FE, CU, CR, MN, ZN, NI, EC, OC, S), names_to = 'species', values_to = 'conc') 
  
  # set up list of species we care about
  species_list <- c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")
  
  # current_species <- species_list[1] # test
  # map over each species 
  reg_spec_df <- map_df(species_list, function(current_species) {
    
    current_species_df <- selected_spec_df %>% 
      filter(species == current_species)
    
    # get list of regions
    region_list <- unique(current_species_df$region)
    
    # current_region <- region_list[1]
    # map over each region
    region_df <- map_df(region_list, function(current_region) {
      
      # region 
      current_region_df <- current_species_df %>% 
        filter(region == current_region)
      
      # run the regressions ----------------------------------------------------------
      current_species_region_reg = feols(conc ~ smokePM + nonsmokePM | monitor_month + year, current_region_df)
      
      # NOW USE THE REGRESSION RESULTS TO PREDICT FRACTION OF CHEMICAL IN PM2.5
      nonsmoke0_df <- current_region_df %>%
        # 0 out nonsmoke
        mutate(nonsmokePM = 0) 
      
      # now predict with new data
      # predict the smoke concentration
      predicted_smoke_conc_df <- nonsmoke0_df %>% 
        # this gets the fraction of a species attributable to smoke
        mutate(pred_Sconc = predict(current_species_region_reg)) %>% 
        dplyr::select(Dataset:region, year, month, monitor_month, 
                      species, pred_Sconc) 
      
      # smoke0_df <- current_region_df %>%
      #   # 0 out nonsmoke
      #   mutate(smokePM = 0)
      
      # predict the nonsmoke concentration
    #   predicted_nonsmoke_conc_df <- smoke0_df %>% 
    #     # this gets the fraction of a species attributable to nonsmoke
    #     mutate(pred_NSconc = predict(current_species_region_reg)) %>% 
    #     dplyr::select(Dataset:region, year, month, monitor_month, 
    #                   nonsmokePM, totPM2.5, species, pred_NSconc)
    #   
    #   
    # pred_pm_concs <- predicted_smoke_conc_df %>% 
    #   left_join(predicted_nonsmoke_conc_df) 
        
    }) %>%
      bind_rows() # end map over the regions

  }) %>%
    bind_rows() # end map over the species
  

  reg_spec_wide <- reg_spec_df %>% 
    filter(pred_Sconc > 0) %>% 
    group_by(Dataset, SiteCode, Date, State, region, year, month, monitor_month) %>%
    mutate(pred_allS_conc = sum(pred_Sconc, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(frac_in1ug_smoke = pred_Sconc/pred_allS_conc) %>% 
    group_by(species, year, month, SiteCode) %>% 
    dplyr::summarise(avg_frac = mean(frac_in1ug_smoke, na.rm = TRUE)) %>% 
    ungroup()
  
      
  # calculate total concentration in PM2.5 
  predicted_concs <- predicted_smoke_conc_df %>% 
    left_join(predicted_nonsmoke_conc_df) %>% 
    filter(nonsmoke_conc > 0) %>% 
    filter(smoke_conc > 0) %>% 
    mutate(tot_pred_conc = nonsmoke_conc + smoke_conc)
  
  # calculate avg attributable fraction of species concentration to total concentration
  avg_conc_site_month_yr <- predicted_concs %>%
    # get avg concentration across all sites for a given month and year
    group_by(species, year, month, SiteCode) %>% 
    dplyr::summarise(avg_smoke_conc = mean(smoke_conc, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(nonsmoke_conc, na.rm = TRUE),
                     avg_tot_conc = mean(tot_pred_conc, na.rm = TRUE)) %>% 
    ungroup()
  
  # now get the average concentration across all sites, for a given month and year
  avg_conc_month_yr <- avg_conc_site_month_yr %>%
    # get avg concentration across all sites for a given month and year
    group_by(species, year, month) %>% 
    dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                     avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
    mutate(spec_type = case_when(
      species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
      species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
      species=="EC" | species=="OC" ~ "Secondary organic",
      species=="S" ~ "Toxicity potentiator",
      TRUE ~ NA)) %>% 
    mutate(spec_f = factor(species, levels = c('CR', "V", "PB", 'NI', "CU", "MN", "ZN","FE", 'S', "EC", "OC"))) 
  
  
  # adjust dataframe so that it can plot stacked area plot
  pred_spec_conc_long_df <- avg_conc_month_yr %>% 
    pivot_longer(cols = c(avg_smoke_conc, avg_nonsmoke_conc, avg_tot_conc), 
                 names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
    mutate(pm_type_f = factor(pm_type, 
                              levels = c('avg_smoke_conc', 'avg_nonsmoke_conc', 'avg_tot_conc')))
  
  datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-01-01"), by = "24 month")
  
  
  # MAKE AN AREA PLOT OF AVG PM OVER TIME ------------------------------------------
  stacked_species_area_plot <- ggplot(pred_spec_conc_long_df %>% 
                                        filter(pm_type != 'avg_tot_conc'),
                                      aes(x = mon_yr, y = avg_spec_conc, 
                                          fill = pm_type_f)) + 
    geom_area(alpha =.6) +
    # add in line for total predicted concentration
    # add in scale
    scale_x_date(labels = date_format("%Y"), breaks = datebreaks) +
    # add in fill manually
    scale_fill_manual(name = 'PM2.5 Type',
                      values=c("firebrick3", "grey70", 'black'), 
                      labels = c('Smoke attributable concentration (ug/m3)',
                                 'Nonmoke attributable concentration (ug/m3)',
                                 'Total predicted species concentration in PM2.5'
                      )) +
    facet_wrap(~spec_f, scales = 'free', ncol = 3) +
    labs(x = 'Year',
         y = paste('Species Concentration (ug/m3)'),
         title = 'Contribution of smoke and nonsmoke PM2.5 to species concentration over time') +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "grey30", size = .5),
          plot.title = element_text(size=12, face='bold'),
          axis.text.x = element_text(face="plain", size = 10),
          axis.text.y = element_text(face="plain", size = 10),
          legend.position = c(0.85, 0.1)) 
  stacked_species_area_plot
  
  # save file
  ggsave(
    filename = 'Fig4_smoke_attributable_fraction_trendAREA.pdf',
    plot = stacked_species_area_plot,
    path = file.path(wip_gdrive_fp, 'figures/'),
    scale = 1,
    width = 12,
    height = 6,
    dpi = 320)   
  
  # MAKE AN AREA PLOT OF AVG PM OVER TIME BY YEAR ------------------------------------------
  avg_conc_yr <- avg_conc_site_month_yr %>%
    # get avg concentration across all sites for a given month and year
    group_by(species, year) %>% 
    dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                     avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(yr = as.Date(paste0(year, "-", "01", "-", "01"), format = "%Y-%m-%d")) %>% 
    mutate(spec_type = case_when(
      species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
      species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
      species=="EC" | species=="OC" ~ "Secondary organic",
      species=="S" ~ "Toxicity potentiator",
      TRUE ~ NA)) %>% 
    mutate(spec_f = factor(species, levels = c('CR', "V", "PB", 'NI', "CU", "MN", "ZN","FE", 'S', "EC", "OC"))) 
  
  
  # adjust dataframe so that it can plot stacked area plot
  pred_spec_conc_long_df <- avg_conc_yr %>% 
    pivot_longer(cols = c(avg_smoke_conc, avg_nonsmoke_conc, avg_tot_conc), 
                 names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
    mutate(pm_type_f = factor(pm_type, 
                              levels = c('avg_smoke_conc', 'avg_nonsmoke_conc', 'avg_tot_conc')))
  
  datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-01-01"), by = "24 month")
  
  
  # MAKE AN AREA PLOT OF AVG PM OVER TIME ------------------------------------------
  YRstacked_species_area_plot <- ggplot(pred_spec_conc_long_df %>% 
                                          filter(pm_type != 'avg_tot_conc'),
                                        aes(x = yr, y = avg_spec_conc, 
                                            fill = pm_type_f)) + 
    geom_area(alpha =.6) +
    # add in scale
    scale_x_date(labels = date_format("%Y"), breaks = datebreaks) +
    # add in fill manually
    scale_fill_manual(name = 'PM2.5 Type',
                      values=c("firebrick3", "grey70"), 
                      labels = c('Smoke attributable concentration (ug/m3)',
                                 'Nonmoke attributable concentration (ug/m3)')) +
    facet_wrap(~spec_f, scales = 'free', ncol = 3) +
    labs(x = 'Year',
         y = paste('Species Concentration (ug/m3)'),
         title = 'Contribution of smoke and nonsmoke PM2.5 to species concentration over time') +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "grey30", size = .5),
          plot.title = element_text(size=12, face='bold'),
          axis.text.x = element_text(face="plain", size = 10),
          axis.text.y = element_text(face="plain", size = 10),
          legend.position = c(0.85, 0.1)) 
  YRstacked_species_area_plot
  
  # save file
  ggsave(
    filename = 'Fig4_YRsmoke_attributable_fraction_trendAREA.png',
    plot = YRstacked_species_area_plot,
    path = file.path(wip_gdrive_fp, 'figures/'),
    scale = 1,
    width = 12,
    height = 6,
    dpi = 320)   
  
  
  return(pred_spec_conc_long_df)
}

