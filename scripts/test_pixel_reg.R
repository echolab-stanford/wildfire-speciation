# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: run full sample new regression model for gridded data for each species
# 
# pm_fp = file.path(data_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')
# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(c(clean_PMspec_df, clean_PMspec_sites_df, parameter_categories), cache = drake_cache)

run_pixel_regression_and_plot <- function(clean_PMspec_df, clean_PMspec_sites_df, 
                                          parameter_categories, pm_fp, grid_fp) {
  
  # smoke pm 
  smoke_pm <- readRDS(pm_fp) %>% 
    rename(Date = 'date') %>% 
    distinct()
  
  # read in gridded smoke data
  grid_10km <- st_read(grid_fp) %>%
    st_transform(4326) %>% 
    # calculate centroids of each grid cell to be used in predictions
    mutate(grid_center = st_centroid(geometry)) %>% 
    mutate(long_grid =  st_coordinates(grid_center)[,1],
           lat_grid =  st_coordinates(grid_center)[,2]) %>% 
    dplyr::select(-COORDX, -COORDY, -grid_center)
  
  # adjust the smoke PM
  gridded_PM <- smoke_pm %>% 
    left_join(grid_10km %>% 
                rename(grid_id_10km = 'ID') %>% 
                st_drop_geometry(), 
              by = 'grid_id_10km') %>% 
    dplyr::select(grid_id_10km, Date, smokePM = 'smokePM_pred', long_grid, lat_grid) 
  
  # make points spatial data for monitors in df
  points <- clean_PMspec_sites_df %>% 
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
  
  # merge grid cell coordinates with monitoring station coordinates
  sites_w_grid_cells <- points %>% 
    st_drop_geometry() %>% 
    left_join(grid_10km %>% 
                rename(grid_id_10km = 'ID') %>% 
                st_drop_geometry(), by = 'grid_id_10km') 
  
  # Add grid_id_10km to the speciation df, join the smoke PM data by grid cell
  reg_df <- clean_PMspec_df %>% 
    mutate(Date = as.Date(Date, format = '%Y-%m-%d')) %>% 
    left_join(sites_w_grid_cells,
              by = c("Dataset", "site_id")) %>% 
    dplyr::select(Date, site_id, month, year, monitor_month, 
                  AS, OC, PB, long, lat, smokePM, grid_id_10km, long_grid, lat_grid)
 
  # # get all the smoke PM data for sites in grid cells 
  # sites_w_PM_df <- reg_df %>% 
  #   left_join(gridded_PM, 
  #             by = c('Date', 'grid_id_10km')) %>% 
  #   mutate(month = month(Date),
  #          year = year(Date)) %>% 
  #   mutate(monitor_month = paste0(site_id, "_", month))
  
  # drop dfs to save memory
  rm(grid_10km, smoke_pm)
  
  # -----------------------------------------------------------------------------
  # 3. RUN REGRESSION FOR SPECIES
  #   A) RUN IN LEVELS ACROSS FULL SAMPLE
  # -----------------------------------------------------------------------------
  # the original list:
  latlong_modelOC = feols(OC ~ smokePM + smokePM*lat + smokePM*long +  smokePM*long*lat |
                            monitor_month + year,
                          data = reg_df,
                          cluster = 'site_id')
  
  # create data to predict conc
  newdata <- gridded_PM %>% 
    # adjust to match the other data, when NA, we assume smoke PM is 0
    mutate(smokePM = ifelse(is.na(smokePM), 0, smokePM)) %>% 
    left_join(reg_df %>% 
                dplyr::select(Date, month, year, monitor_month, grid_id_10km),
              by = c('Date','grid_id_10km')) %>% 
    rename(lat = 'lat_grid', long = 'long_grid')
  
  # predict using FE
  predOC_full_model_df <- newdata %>% 
    mutate(predOC = predict(latlong_modelOC, newdata = newdata))
  
  hist(predOC_full_model_df$predOC)
  
  # predict using the non FE method
  latlong_model = feols(c(OC, PB, AS) ~ smokePM + smokePM*lat + smokePM*long +  smokePM*long*lat |
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
  
  # # now get the lat long associated with each grid cell's centroid to use in predictions
  # gridded_PM <- smoke_pm %>%
  #   left_join(grid_10km %>%
  #               rename(grid_id_10km = 'ID') %>%
  #               st_drop_geometry(),
  #             by = 'grid_id_10km') %>%
  #   dplyr::select(grid_id_10km, Date, smokePM_pred, long_grid, lat_grid)
  
  # -----------------------------------------------------------------------------
  # map over each species and predict
  # -----------------------------------------------------------------------------
  species_list <- unique(coeffs$species)
  current_species <- species_list[1] # test
  
  # PIVOT WIDER AND JOIN TO DATAFRAME
  # all_pred_grid_cells_sf = map_df(species_list, function(current_species) {
    
    current_coeffs <- coeffs %>% 
      filter(species == current_species) %>% 
      dplyr::select(species, coefficient, B) %>% 
      pivot_wider(names_prefix = "B_", 
                  names_from = 'coefficient', 
                  values_from = 'B') %>% 
      left_join(parameter_categories)
    
    # add in column for predicted Beta
    current_preds_df <- gridded_PM %>% 
      mutate(B_smokePM = current_coeffs$B_smokePM,
             B_smokePM_lat = current_coeffs$`B_smokePM:lat`,
             B_smokePM_long = current_coeffs$`B_smokePM:long`,
             B_smokePM_lat_long = current_coeffs$`B_smokePM:lat:long`) %>%
      # calculate marginal effect (BsmokePM + BsmokePM*lat + BsmokePM*long + BsmokePM*lat*long)*smokePMconc
      mutate(pred_grid_conc = (B_smokePM + B_smokePM_lat*lat_grid + B_smokePM_long*long_grid + B_smokePM_lat_long*long_grid*lat_grid)*smokePM) %>% 
      mutate(species = current_species)
    
    # check  
    hist(current_preds_df$pred_grid_conc)
    
    
    # COMPARE PREDICTIONS AT SAME LOCATION + DATE ------------------------------
    compared_preds <- current_preds_df %>% 
      left_join(predOC_full_model_df %>% 
                  dplyr::select(Date, grid_id_10km, predOC), 
                by =c('Date', 'grid_id_10km')) %>% 
      filter(!is.na(predOC))
    
    ratio <- compared_preds %>% 
      mutate(`preds_wFE:noFE` = pred_grid_conc/predOC) %>% 
      dplyr::select(grid_id_10km, Date, predOC, predOC_noFE = 'pred_grid_conc') %>% 
      pivot_longer(cols = c('predOC', 'predOC_noFE'), names_to = 'df', values_to = 'conc') %>% 
      mutate(year = year(Date))
    
   # hist(ratio$`preds_wFE:noFE`)
    
    # Plot the time series 
    datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")
    
    # PLOTTING
    test_density <- ggplot(ratio) +
      geom_density(aes(x = Date,
                       y = conc,
                       fill = df,
                       color = df), 
                   stat = "identity",
                   alpha = 0.6) +
      scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
      geom_hline(aes(yintercept = 0)) +
      theme_minimal() 
    test_density

    
    
    
    
    
    
    # LOOK AT THE PLOTS ----------------------------------------
    avg_noFE <- current_preds_df %>% 
      mutate(month = month(Date),
             year = year(Date)) %>% 
      # group_by(grid_id_10km, month, year) %>% 
      # dplyr::summarise(avg_predOC = mean(pred_grid_conc, na.rm = TRUE)) %>% 
      # ungroup() %>% 
      group_by(month, year) %>% 
      dplyr::summarise(avg_predOC = mean(pred_grid_conc, na.rm = TRUE)) %>% 
      ungroup() %>%
      mutate(df = 'noFE')
    
    # get average by month for each data set and plot
    avg_wFE <- predOC_full_model_df %>% 
      mutate(month = month(Date),
             year = year(Date)) %>% 
      # group_by(grid_id_10km, month, year) %>% 
      # dplyr::summarise(avg_predOC = mean(predOC, na.rm = TRUE)) %>% 
      # ungroup() %>% 
      group_by(month, year) %>% 
      dplyr::summarise(avg_predOC = mean(predOC, na.rm = TRUE)) %>% 
      ungroup() %>% 
      mutate(df = 'wFE') 
    
    #merge together and plot
    all_data <- avg_wFE %>% 
      bind_rows(avg_noFE)
    
    # Plot the time series 
    plot <- ggplot(all_data, aes(x = as.factor(month), y = avg_predOC)) +
      geom_line(aes(color = df, group = df)) +
      facet_wrap(~year) +
      labs(title = "comparing predicted concentrations over time",
           x = "Date",
           y = "conc") +theme_bw()
    plot
    
    
    # plot by grid
    avg_noFE_grid <- current_preds_df %>% 
      mutate(month = month(Date),
             year = year(Date)) %>% 
      group_by(grid_id_10km, lat_grid, long_grid) %>%
      dplyr::summarise(avg_predOC = mean(pred_grid_conc, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(df = 'noFE') %>% 
      left_join(grid_10km %>% dplyr::select(grid_id_10km = 'ID', geometry), by = 'grid_id_10km')
    
    # get average by month for each data set and plot
    avg_wFE_grid <- predOC_full_model_df %>% 
      mutate(month = month(Date),
             year = year(Date)) %>% 
      group_by(grid_id_10km, lat, long) %>%
      dplyr::summarise(avg_predOC = mean(predOC, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(df = 'wFE') %>% 
      left_join(grid_10km %>% dplyr::select(grid_id_10km = 'ID', geometry), by = 'grid_id_10km')
    
    # all_data <- avg_wFE_grid %>% 
    #   bind_rows(avg_noFE_grid)
    # # summarize across full sample to get avg exposure over 15 yrs
    # sample_avg_predictions <- current_preds_df %>%
    #   group_by(grid_id_10km, species, long_grid, lat_grid) %>%
    #   dplyr::summarise(avg_pred_spec_conc = mean(pred_grid_conc, na.rm = TRUE)) %>%
    #   ungroup() %>% 
    #   left_join(grid_10km %>% dplyr::select(grid_id_10km = 'ID', geometry), by = 'grid_id_10km') 
    
    
    # PLOT QUANTILES FOR EACH SPECIES AND REGIONS COMBINED
    current_pred_map <- ggplot() +
      geom_sf(data = avg_wFE_grid %>%
                st_as_sf(), aes(color = avg_predOC), alpha =.70) +  # Add point data
      scale_colour_viridis(option = "A", direction = -1) + 
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
      labs(color = paste0('Average annual\n', current_coeffs$species_long,  ' \nConcentration (ug/m3)'))
    # Specify the order for the legend
     current_pred_map
    
    # save the map
    ggsave(filename = paste0("Fig5_continuous_gridded_predictions", current_species, ".pdf"),
           plot = current_pred_map,
           path = file.path(results_fp, 'figures/Fig5'),
           scale = 1,
           width = 7,
           height = 6,
           dpi = 320)
    
    sample_avg_predictions
    
  }) %>%
    bind_rows()
  
  
  # convert the sf object to a df by dropping geometry
  all_pred_grid_cells_df <- all_pred_grid_cells_sf %>%
    st_drop_geometry()
  
  
  return(all_pred_grid_cells_df)
  
}


# -----------------------------------------------------------------------
# get quantiles
# binned_predictions <- tibble(
#   bins = cut(sample_avg_predictions$avg_pred_spec_conc,
#              breaks = quantile(sample_avg_predictions$avg_pred_spec_conc,
#                                probs = seq(0,1,.20),
#                                na.rm = TRUE),
#              include.lowest = T, dig.lab = 5)
# ) %>%
#   distinct() %>%
#   arrange(bins) %>%
#   bind_cols(
#     tibble(
#       percentile = names(quantile(sample_avg_predictions$avg_pred_spec_conc,
#                                   probs = seq(0,1,.20),
#                                   na.rm = TRUE))[-1]
#     )
#   ) %>%
#   # format bin labels for legend
#   separate(bins, into = c("start", "end"), sep = ",", remove = FALSE) %>%
#   mutate(across(c("start","end"), ~str_remove(.,"\\(|\\)|\\[|\\]") %>% as.numeric(.))) %>%
#   mutate(across(c("start","end"), ~format(round(.,0), trim = TRUE, big.mark = ","), .names = "{.col}_label")) %>%
#   mutate(start_label = if_else(start_label=="0", as.character(signif(start, 1)), start_label),
#          end_label = if_else(end_label=="0", as.character(signif(end, 1)), end_label)) %>%
#   mutate(bin_label = paste0(percentile," ",str_extract(bins,"^\\(|\\["),start_label," to ",end_label,str_extract(bins,"\\)|\\]"))) %>%
#   # remove outer limits
#   mutate(bin_label = case_when(
#     row_number()==1 ~ paste0(percentile," (\u2264",end_label,"]"),
#     row_number()==n() ~ paste0(percentile," (>",start_label,")"),
#     TRUE ~ bin_label)) %>%
#   dplyr::select(bins, start, end, bin_label) %>%
#   mutate(order = row_number())
# 
# order_bins <- binned_predictions %>%
#   distinct(bin_label, order) %>%
#   arrange(order) %>%
#   distinct(bin_label) %>%
#   as.vector()
# 
# binned_predictions_ordered <- binned_predictions %>%
#   mutate(bin_label = fct_relevel(bin_label, order_bins))
# 
# # add bins back to current chemical species gridded predictions
# gridded_preds_w_bins <- sample_avg_predictions %>%
#   mutate(bins = cut(sample_avg_predictions$avg_pred_spec_conc,
#                     breaks = quantile(sample_avg_predictions$avg_pred_spec_conc,
#                                       probs = seq(0,1,.20),
#                                       na.rm = TRUE),
#                     include.lowest = T, dig.lab = 5)) %>%
#   left_join(binned_predictions_ordered, by = "bins") %>%
#   dplyr::select(-bins, -start, -end)

# # PLOT QUANTILES FOR EACH SPECIES AND REGIONS COMBINED
#   current_pred_map <- ggplot() +
#     geom_sf(data = gridded_preds_w_bins %>%
#               st_as_sf(), aes(color = bin_label), alpha =.70) +  # Add point data
#     scale_colour_viridis_d(option = "A", direction = -1)+
#     theme_minimal() +  # Apply a minimal theme
#     theme(
#       panel.grid.major = element_blank(),  # Remove major grid lines
#       panel.grid.minor = element_blank(),  # Remove minor grid lines
#       axis.text = element_blank(),  # Remove axis labels
#       axis.title = element_text(),  # Remove axis titles
#       legend.title = element_text(),
#       plot.title = element_text(face = "bold"),
#       legend.key.width = unit(.4, "cm"),  # Set the legend key width
#       legend.key.height = unit(.4, "cm")  # Set the legend key height
#     ) #+
#     #labs(color = paste0('Average annual\n', current_name,  ' \nConcentration (ug/m3)'))
#   # Specify the order for the legend
#   current_pred_map
