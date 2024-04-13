# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: run full sample new regression model for gridded data for each species

# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(c(clean_PMspec_df, parameter_categories, full_samplePM_df), cache = drake_cache)


run_pixel_regression_and_plot_5yr_maps <- function(clean_PMspec_df, parameter_categories, 
                                                   grid_fp, full_samplePM_df, current_species) {
  
  # parellelize
  future::plan(multisession)
  
  # read in grids + get centroid of each grid cell
  grid_10km <- st_read(grid_fp) %>%
      st_transform(4326) %>% 
      # calculate centroids of each grid cell to be used in predictions
      mutate(grid_center = st_centroid(geometry)) %>% 
      mutate(long_grid =  st_coordinates(grid_center)[,1],
             lat_grid =  st_coordinates(grid_center)[,2]) %>% 
      dplyr::select(-COORDX, -COORDY, -grid_center)
  
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
    dplyr::select(Dataset, year, month, monitor_month, Date,
                    site_id, MF_adj, smokePM, nonsmokePM_MF,
                    AS, PB, NI, long, lat) %>% 
    left_join(sites_w_grid_cells,
              by = c("Dataset", "site_id")) 
  
  # merge lat long with smoke data for all days (not just smoke days)
  gridded_PM_all_days <- full_samplePM_df %>% 
    left_join(grid_10km %>% 
                rename(grid_id_10km = 'ID') %>% 
                st_drop_geometry(), 
              by = 'grid_id_10km') %>% 
    dplyr::select(grid_id_10km, Date, smokePM, long_grid, lat_grid) 
  
  rm(full_samplePM_df) # to save memory
  
  # -----------------------------------------------------------------------------
  # RUN REGRESSION FOR SPECIES
  # -----------------------------------------------------------------------------
  # lat long regression model-
  latlong_model = feols(c(AS, PB, NI) ~ smokePM + smokePM*lat + smokePM*long +  smokePM*long*lat |
                            monitor_month + year,
                          data = reg_df,
                          cluster = 'site_id')

  ## calculate 95 CI%
  CIs <- confint(latlong_model) %>% 
    rename(CI25 = `2.5 %`,
           CI975 = `97.5 %`,
           species = 'lhs') %>% 
    dplyr::select(-id)
    
  # get coefficients and prepare for plotting
  coeffs <- coeftable(latlong_model) %>%
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           B = 'Estimate') %>%
    mutate(pval = round(pval, digits = 3)) %>%
    dplyr::select(-id) %>%
    left_join(CIs, by = c('species', 'coefficient'))
  
  rm(latlong_model) # drop to save memory
  
  # -----------------------------------------------------------------------------
  # map over each species and predict
  # -----------------------------------------------------------------------------
    # grab coefficients for current species
    current_coeffs <- coeffs %>% 
      filter(species == current_species) %>% 
      dplyr::select(species, coefficient, B) %>% 
      pivot_wider(names_prefix = "B_", 
                  names_from = 'coefficient', 
                  values_from = 'B') %>% 
      left_join(parameter_categories)
      
    # using current coefficients, predict for each grid cell using that grid cell's lat long
    current_preds_df <- gridded_PM_all_days %>% 
      mutate(B_smokePM = current_coeffs$B_smokePM,
             B_smokePM_lat = current_coeffs$`B_smokePM:lat`,
             B_smokePM_long = current_coeffs$`B_smokePM:long`,
             B_smokePM_lat_long = current_coeffs$`B_smokePM:lat:long`) %>%
      # calculate marginal effect (BsmokePM + BsmokePM*lat + BsmokePM*long + BsmokePM*lat*long)*smokePMconc
      mutate(pred_grid_conc = (B_smokePM + B_smokePM_lat*lat_grid + B_smokePM_long*long_grid + B_smokePM_lat_long*long_grid*lat_grid)*smokePM) %>%
      mutate(ratio = (B_smokePM + B_smokePM_lat*lat_grid + B_smokePM_long*long_grid + B_smokePM_lat_long*long_grid*lat_grid)) %>% 
      mutate(species = current_species) %>% 
    # prep for plotting + exposure analysis -> want to see how things have changed spatially over this time period
    mutate(year = year(Date)) %>% 
      # categorize into 3 periods of 5 years
    mutate(samp_period = case_when(
      year > 2005 & year < 2011 ~ '2006-2010',
      year > 2010 & year < 2016 ~ '2011-2015',
      year > 2015 & year < 2021 ~ '2016-2020'
    ))
    
    # save current species predictions
    # write.fst(current_preds_df, file.path(data_fp, paste0('clean/', current_species, '_gridded_preds.fst')))
    rm(gridded_PM_all_days)
    # ------------------------------------------------------------------------
    # DATA CHECK - make sure the betas added up seem reasonable when compared to regional model
    # test <- current_preds_df %>% 
    #   filter(Date == '2020-01-10') %>% 
    #   mutate(ratio = (B_smokePM + B_smokePM_lat*lat_grid + B_smokePM_long*long_grid + B_smokePM_lat_long*long_grid*lat_grid))
    # ------------------------------------------------------------------------
      # summarize across sample periods to get avg exposure for five year periods in each grid cell
    sample_avg_predictions <- current_preds_df %>%
      filter(samp_period != '2011-2015') %>% 
      mutate(across(where(is.numeric), ~replace(., . < 0, NA))) %>% 
        group_by(grid_id_10km, species, long_grid, lat_grid, samp_period) %>%
        dplyr::summarise(avg_grid_conc = mean(pred_grid_conc, na.rm = TRUE)) %>%
        ungroup() 
      
    
    rm(current_preds_df) # drop to save memory
    
    # PLOT!!!!!!!!!!!!!!!!
    if (current_species == 'PB') {
      current_pred_map <- ggplot() +
        geom_sf(data = sample_avg_predictions %>%
                  left_join(grid_10km %>% 
                              dplyr::select(grid_id_10km = 'ID', geometry), 
                            by ='grid_id_10km') %>% 
                  st_as_sf(), 
                aes(color = avg_grid_conc, 
                    fill = avg_grid_conc)) +  
        scale_colour_viridis_c(option = "mako", direction = -1) +
        scale_fill_viridis_c(option = "mako", direction = -1) +
        facet_wrap(~samp_period) +
        theme_minimal() +  # Apply a minimal theme
        theme(
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.text = element_blank(),  # Remove axis labels
          axis.title = element_text(),  # Remove axis titles
          legend.title = element_text(),
          plot.title = element_text(face = "bold"),
          legend.key.width = unit(.4, "cm"),  # Set the legend key width
          legend.key.height = unit(.4, "cm")  # Set the legend key height
        ) +
        labs(title = paste0(current_species),
             color = paste0('Concentration (ug/m3)')) +
        guides(fill = 'none')
      # Specify the order for the legend
      
      # current_pred_map
      ggsave(
        filename = paste0('Fig5A_continuous_gridded_predictions', current_species, '_conc.png'),
        plot = current_pred_map,
        path = file.path(results_fp, 'figures/Fig5'),
        scale = 1,
        width = 12,
        height = 8,
        dpi = 320)
      ggsave(
        filename = paste0('Fig5A_continuous_gridded_predictions', current_species, '_conc.pdf'),
        plot = current_pred_map,
        path = file.path(results_fp, 'figures/Fig5'),
        scale = 1,
        width = 12,
        height = 8,
        dpi = 320)
      
    } 
    
    # change color scale based on 
    if (current_species == 'AS') {
      current_pred_map <- ggplot() +
        geom_sf(data = sample_avg_predictions %>%
                  left_join(grid_10km %>% 
                              dplyr::select(grid_id_10km = 'ID', geometry), 
                            by ='grid_id_10km') %>% 
                  st_as_sf(), 
                aes(color = avg_grid_conc, 
                    fill = avg_grid_conc)) +  
        scale_color_continuous_sequential(palette = "ag_GrnYl", rev = T) +
        scale_fill_continuous_sequential(palette = "ag_GrnYl", rev = T) +
        facet_wrap(~samp_period) +
        theme_minimal() +  # Apply a minimal theme
        theme(
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.text = element_blank(),  # Remove axis labels
          axis.title = element_text(),  # Remove axis titles
          legend.title = element_text(),
          plot.title = element_text(face = "bold"),
          legend.key.width = unit(.4, "cm"),  # Set the legend key width
          legend.key.height = unit(.4, "cm")  # Set the legend key height
        ) +
        labs(title = paste0(current_species),
             color = paste0('Concentration (ug/m3)')) +
        guides(fill = 'none')
      current_pred_map
      
      ggsave(
        filename = paste0('Fig5A_continuous_gridded_predictions', current_species, '_conc.png'),
        plot = current_pred_map,
        path = file.path(results_fp, 'figures/Fig5'),
        scale = 1,
        width = 12,
        height = 8,
        dpi = 320)
      ggsave(
        filename = paste0('Fig5A_continuous_gridded_predictions', current_species, '_conc.pdf'),
        plot = current_pred_map,
        path = file.path(results_fp, 'figures/Fig5'),
        scale = 1,
        width = 12,
        height = 8,
        dpi = 320)
      
    }
    
    if (current_species == 'NI') {
      current_pred_map <- ggplot() +
        geom_sf(data = sample_avg_predictions %>%
                  left_join(grid_10km %>% 
                              dplyr::select(grid_id_10km = 'ID', geometry), 
                            by ='grid_id_10km') %>% 
                  st_as_sf(), 
                aes(color = avg_grid_conc, 
                    fill = avg_grid_conc)) +  
        scale_colour_viridis_c(option = "magma", direction = -1) +
        scale_fill_viridis_c(option = "magma", direction = -1) +
        facet_wrap(~samp_period) +
        theme_minimal() +  # Apply a minimal theme
        theme(
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.text = element_blank(),  # Remove axis labels
          axis.title = element_text(),  # Remove axis titles
          legend.title = element_text(),
          plot.title = element_text(face = "bold"),
          legend.key.width = unit(.4, "cm"),  # Set the legend key width
          legend.key.height = unit(.4, "cm")  # Set the legend key height
        ) +
        labs(title = paste0(current_species),
             color = paste0('Concentration (ug/m3)')) +
        guides(fill = 'none')
      current_pred_map
      
      ggsave(
        filename = paste0('Fig5A_continuous_gridded_predictions', current_species, '_conc.png'),
        plot = current_pred_map,
        path = file.path(results_fp, 'figures/Fig5'),
        scale = 1,
        width = 12,
        height = 8,
        dpi = 320)
      ggsave(
        filename = paste0('Fig5A_continuous_gridded_predictions', current_species, '_conc.pdf'),
        plot = current_pred_map,
        path = file.path(results_fp, 'figures/Fig5'),
        scale = 1,
        width = 12,
        height = 8,
        dpi = 320)
    }
      
      sample_avg_predictions
      # End the parallel session
      future::plan(NULL)

    return(sample_avg_predictions)
  }


