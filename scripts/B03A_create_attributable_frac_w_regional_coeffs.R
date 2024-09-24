# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species

# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(c(parameter_categories, pm_pal, clean_PMspec_df,
#         regionalPMcoeffs_normalized, full_samplePM_df),
#       cache = drake::drake_cache(".drake"))

# function
create_attributable_frac_w_regional_coeffs <- function(grid_fp, clean_PMspec_df, parameter_categories, 
                                                       pm_pal, full_samplePM_df, regionalPMcoeffs_normalized) {
  # read in grid
  grid_10km <- st_read(grid_fp) %>%
    st_transform(4326)
  
  # make points spatial data for monitors in df
  points <- clean_PMspec_df %>% 
    # get a unique list of sites
    distinct(Dataset, site_id, long, lat) %>% 
    # convert so sf object
    st_as_sf(coords = c(x ='long', y="lat"), crs = st_crs(grid_10km))
  
  # pull id cells that intersect with the monitoring sites
  points$grid_id_10km <- apply(
    st_intersects(grid_10km, points, sparse = FALSE), 2, 
    function(col) {
      grid_10km[which(col), ]$ID
    }) 
  
  sites_w_grid_cells <- points %>% 
    st_drop_geometry() %>% 
    left_join(grid_10km %>% 
                rename(grid_id_10km = 'ID') %>% 
                st_drop_geometry(), by = 'grid_id_10km') 
  
  # select the vars that are needed for the
  reg_df <- clean_PMspec_df %>%
    mutate(monitor_month = paste0(site_id,"_",month)) %>%
    dplyr::select(-ZR, -SOIL)
  
  # parellelize
  future::plan(multisession)
  
  # ------------------------------------------------------------------
  # REGIONAL ANALYSIS ----------------------------------------
  # ------------------------------------------------------------------
  coeffs <- regionalPMcoeffs_normalized %>% 
    dplyr::select(region, species, species_name, pm_type, Estimate) %>% 
  pivot_wider(names_from = 'pm_type', names_prefix = "est_", values_from = 'Estimate')
  
  rm(regionalPMcoeffs_normalized)
  
  # -----------------------------------------------------------------------------
  # map over each species and predict
  # -----------------------------------------------------------------------------
  # predict over our sample but with regional betas
  predicted_spec_smoke_conc <- reg_df %>%
    dplyr::select(Dataset, region, year, month, Date,
                  smoke_day, site_id, MF_adj, smokePM,
                  nonsmokePM_MF, AL:ZN) %>%
    pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>%
    filter(!is.na(conc_val)) %>%
    left_join(coeffs, by = c('region', 'species')) %>%
    # filter(nonsmokePM_MF >= 0 & smokePM >= 0) %>% 
    mutate(pred_spec_nonsmoke = est_nonsmokePM_MF*nonsmokePM_MF,
           pred_spec_smoke = est_smokePM*smokePM) %>%
    # filter(pred_spec_smoke >= 0 & pred_spec_nonsmoke >= 0) %>% 
    mutate(pred_spec_allPM = pred_spec_smoke + pred_spec_nonsmoke) %>% 
    # make sure we know which grid cell the monitors are in
    left_join(sites_w_grid_cells) 
    
  rm(grid_10km, points, coeffs)
  print(paste('merging predictions and concentration for all days'))
  # add to full dataset before averaging, need all days not just smoke days  
  pred_conc_all_days <- full_samplePM_df %>% 
    left_join(predicted_spec_smoke_conc %>%
                dplyr::select(Dataset, region, year, month, Date, site_id, grid_id_10km, 
                              species, species_name, pred_spec_nonsmoke, pred_spec_smoke, pred_spec_allPM),
              by = c("Date","grid_id_10km", "region")) 
  
  rm(full_samplePM_df, predicted_spec_smoke_conc)
  
  # pivot reg df longer
  obs_species_conc_longer <- reg_df %>% 
    pivot_longer(cols = c(AL:ZN), names_to = 'species', values_to = 'obs_conc') %>% 
    dplyr::select(region, year, month, monitor_month, Date, site_id, species, obs_conc) 

  
  # join
  increasing_frac_df <- obs_species_conc_longer %>% 
    left_join(pred_conc_all_days, 
              by = c('species', 'Date', 'year', 'month', 'site_id', "region"))
  
  # calc fraction
  calc_frac <- increasing_frac_df %>% 
    mutate(att_frac = pred_spec_smoke/obs_conc) %>% 
    dplyr::select(Date, species, region, year, month, att_frac, site_id, grid_id_10km, monitor_month) %>% 
    pivot_wider(names_from = 'species', values_from = 'att_frac') 
  
  
  # TEST IF THE ATTRIBUTABLE FRACTION IS INCREASING OVER TIME
  print(paste('running model'))
  att_frac_model = feols(c(AL,AS,BR,CA,CL, CR,CU, EC, FE, K, MG, MN, 
                           `NA`, NI, NO3, OC, P, PB, RB, S, 
                           SE, SI, SO4, SR, TI, V, ZN)
                         ~  year | monitor_month,
                         data = calc_frac, cluster = 'site_id')
  # etable(att_frac_model)
  
  ## calculate 95 CI%
  CIs <- confint(att_frac_model) %>% 
    rename(CI25 = `2.5 %`,
           CI975 = `97.5 %`,
           species = 'lhs') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs <- coeftable(att_frac_model) %>%
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           coeff = 'Estimate') %>%
    mutate(pval = round(pval, digits = 3)) %>%
    dplyr::select(-id) %>%
    left_join(CIs, by = c('species', 'coefficient')) %>% 
    mutate(sig = ifelse(pval < .05, "sig increasing", "cannot detect increase"))
  
rm(att_frac_model)

print(paste('averaging preds for each region'))

  # combine all days with predicted dataset
  avg_conc_site_month_yr <- pred_conc_all_days %>%
  # calculate avg attributable fraction of species concentration to total concentration
    # get avg concentration across all sites for a given month and year
    group_by(region, species, year, month, site_id) %>%
    dplyr::summarise(avg_smoke_conc = mean(pred_spec_smoke, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(pred_spec_nonsmoke, na.rm = TRUE),
                     avg_tot_conc = mean(pred_spec_allPM, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(species, year, month) %>%
    dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                     avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>%
    ungroup()

  # get monthly avg conc for the observed "total concentration" of a species
  obs_improve_monthly_mean <- reg_df %>%
    dplyr::select(Dataset, region, year, month, Date,
                  smoke_day, site_id, MF_adj, smokePM,
                  nonsmokePM_MF, AL:ZN) %>%
    pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>%
    group_by(region, species, year, month, site_id) %>%
    dplyr::summarise(avg_obs_tot_conc = mean(conc_val, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(region,species, year, month) %>%
    dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(species, year, month) %>%
    dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>%
    ungroup()

  # now merge back with predicted data
  print(paste('merging predictions with observed data'))
  avg_attributable_fraction_vs_obs <- avg_conc_site_month_yr %>%
    left_join(obs_improve_monthly_mean, 
              by = c('species', 'year', 'month'))

  # # prep df for plotting
  plot_df <- avg_attributable_fraction_vs_obs %>%
    dplyr::select(-avg_obs_tot_conc, -avg_nonsmoke_conc, -avg_tot_conc) %>%
    mutate(label = 'Attributable concentration to smoke PM2.5') %>%
    rename(conc = 'avg_smoke_conc') %>%
    bind_rows(avg_attributable_fraction_vs_obs %>%
                dplyr::select(conc = 'avg_obs_tot_conc', year, month, species) %>%
                mutate(label = 'Observed Concentration in Total PM2.5')) %>%
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>%
    left_join(parameter_categories, by = 'species') %>%
    mutate(label = fct_relevel(label, c("Observed Concentration in Total PM2.5",
                                        "Attributable concentration to smoke PM2.5"))) %>%
    filter(!is.na(species)) %>% 
    left_join(coeffs %>% rename(year_trend = 'coeff'), by = 'species')
  
  datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")

  limits = c(
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
    "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")


  # # PLOTTING
  avg_mon_pred_spec_plot <- ggplot(plot_df) +
    geom_density(aes(x = mon_yr,
                     y = conc,
                     fill = label,
                     color = label),
                 stat = "identity", alpha = 0.8) +
    labs(title = expression("Chemical species concentration attributable to wildfire smoke PM2.5"),
         y = expression("Concentration (\u00b5g/m"^3*")"),
         x = "") +
    scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
    geom_hline(aes(yintercept = 0)) +
    theme_minimal() +
    scale_color_manual(values = pm_pal) +
    scale_fill_manual(values = pm_pal) +
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
          axis.line.y = element_line(color = "grey10"))
  avg_mon_pred_spec_plot

  # save file
  ggsave(
    filename = 'SIFig_ALLsmoke_attributable_fraction_trend_regional_model.png',
    plot = avg_mon_pred_spec_plot,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
  # save file
  ggsave(
    filename = 'SIFig_ALLsmoke_attributable_fraction_trend_regional_model.pdf',
    plot = avg_mon_pred_spec_plot,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
  
  future::plan(NULL)
  
  return(plot_df)#
}

  

