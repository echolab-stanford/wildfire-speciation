# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: November 14, 2023
# Description: This script uses our estimated coefficients to predict the concentration of species across a grid


# pm_fp = file.path(data_fp, 'intermediate/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds')
# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(c(regionalPMcoeffs_normalized, us_region_map), cache = drake_cache)

# to do - get grid cells in each region and then predict
create_regional_species_conc_predictions <- function(pm_fp, grid_fp, regionalPMcoeffs_normalized, us_region_map) {
  
  # COMBINE GRID CELLS + SMOKE PREDICTIONS:
  # smoke pm 
  smoke_pm <- readRDS(pm_fp) %>% 
    rename(Date = 'date') %>% 
    distinct()
  
  #read in gridded smoke data
  grid_10km <- st_read(grid_fp) %>%
    st_transform(4326) 
  
  # group geometry into regions + remake polygons
  regions_sf <- us_region_map %>% 
    group_by(region) %>% 
    dplyr::summarise(reg_geom = st_union(geometry)) %>% 
    ungroup() %>% 
    st_make_valid()
  
  # determine which grid cells are in which regions by doing a spatial join
  grid_10km_reg <- grid_10km %>% 
    st_join(regions_sf) %>% 
    filter(!is.na(region)) %>% 
    rename(grid_id_10km = 'ID') %>% 
    dplyr::select(-COORDX, -COORDY) 
  
  
  # DETERMINE WHICH GRID CELLS HAVE MONITORS IN THEM
  # CREATE A BUFFER AROUND MONITORS
  # ONLY ESTIMATE CONCENTRATION IN PIXELS AROUND MONITORS
  # MERGE WITH GRIDDED POPULATION + COUNT HOW MANY DAYS ABOVE THREHSOLD
  
  
  
    # create species list
    species_list <- c('AS', 'EC', 'OC', 'PB') # unique(regionalPMcoeffs_normalized$species)
    # current_species <- species_list[4]
    
    # map over each species
    all_species_preds_sf <- map_df(species_list, function(current_species) {
      
      print(current_species)
      
      # filter dataframe of betas to current species + grab estimate
      current_species_betaS <- regionalPMcoeffs_normalized %>%
        filter(species == current_species) %>% 
        filter(pm_type == 'smokePM') %>% 
        distinct(Estimate, region) 
      
      # grab long name
      current_species_name <- regionalPMcoeffs_normalized %>%
        filter(species == current_species) %>% 
        distinct(species_long) %>% 
        as.character()
      
      # predict species concentration due to smoke in each grid cell
      predicted_gridded_species_conc_df <- grid_10km_reg %>%
        # drop geometry to save
        st_drop_geometry() %>% 
        left_join(current_species_betaS, by = 'region') %>% 
        left_join(smoke_pm, by = 'grid_id_10km') %>% 
        mutate(pred_grid_conc = Estimate*smokePM_pred,
               species_long = paste0(current_species_name),
               species = current_species)
      

      # summarize across full sample to get avg exposure over 15 yrs
      sample_avg_predictions <- predicted_gridded_species_conc_df %>%
        mutate(year = year(Date)) %>% 
        mutate(month = month(Date)) %>% 
        group_by(grid_id_10km, region, species_long, month, species, year) %>%
        dplyr::summarise(sample_avg_pred_conc = mean(pred_grid_conc, na.rm = TRUE)) %>%
        ungroup() %>% 
        left_join(grid_10km_reg, by = c('region','grid_id_10km')) 
      
      long_diff <- sample_avg_predictions %>% 
        mutate(mon_yr = as.Date(paste0(year, "-", "01-01"))) %>% 
        group_by(species_long, month, species, year) %>%
        dplyr::summarise(sample_avg_pred_conc = mean(sample_avg_pred_conc, na.rm = TRUE)) %>%
        ungroup() 
  
        # pivot_wider(names_from = 'year', values_from = 'sample_avg_pred_conc') %>% 
        # mutate(change_in_conc = `2020`-`2006`) %>% 
        # mutate(pct_change = change_in_conc/`2006`)
      
      
      # first year of sample minus last year in sample
      # long_diff <- sample_avg_predictions %>% 
      #   filter(year %in% c(2006, 2020)) %>% 
      #   pivot_wider(names_from = 'year', values_from = 'sample_avg_pred_conc') %>% 
      #   mutate(change_in_conc = `2020`-`2006`) %>% 
      #   mutate(pct_change = change_in_conc/`2006`)
      
      
      # Plot the time series
      test <- ggplot(long_diff, aes(x = as.factor(month), y = sample_avg_pred_conc, group = year)) +
        geom_line(aes(color = as.factor(year))) +
        # geom_hline(yintercept = 0.000012, color = 'red') +
        labs(title = "Attributable lead concentration due to wildfire smoke by year",
             x = "Date",
             y = "Conc") +
        facet_wrap(~year, ncol = 3) + 
        theme_minimal()
      test
      
      # region_list <- unique(grid_10km_reg$region)
      
      
      # PLOT QUANTILES FOR EACH SPECIES AND REGIONS COMBINED
      current_pred_map <- ggplot() +
        geom_sf(data = long_diff %>%
                  st_as_sf(), aes(fill = `2020`, color = `2020`), alpha =.9) +  # Add point data
        scale_fill_viridis(option = "A", direction = -1) + 
        scale_color_viridis(option = "A", direction = -1) + 
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
        labs(color = paste0('Average ', current_species_name,  ' \nConcentration (ug/m3) \n in 2020')) +
        guides(fill = 'none')
      current_pred_map
      
      # save the map
      # ggsave(filename = paste0("Fig5_regional_gridded_preds_", current_species_code, "_", current_region, ".pdf"),
      #        plot = current_pred_map,
      #        path = file.path(results_fp, 'figures/Fig5'),
      #        scale = 1,
      #        width = 7,
      #        height = 6,
      #        dpi = 320)
      
      predicted_gridded_species_conc_df # include so that this is the dataframe that is added to in each map
      
    }) %>%
      bind_rows() # bind for all species within a given region
    
    
  }) %>%
    bind_rows() # bind all predictions together for all regions 
  
  return(all_reg_species_preds_df)
  
}

# PLOT QUANTILES FOR EACH SPECIES AND REGIONS COMBINED
# current_pred_map <- ggplot() +
#   # geom_sf(data = us_sf, fill = "grey95") +  # Add the base map
#   geom_sf(data = predicted_gridded_species_conc_df %>%
#             st_as_sf(), 
#           aes(color = pred_grid_conc), alpha =.70) +  # Add point data
#   scale_colour_viridis(option = "A", direction = -1)+
#   theme_minimal() +  # Apply a minimal theme
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     axis.text = element_blank(),  # Remove axis labels
#     axis.title = element_text(),  # Remove axis titles
#     legend.title = element_text(),
#     plot.title = element_text(face = "bold"),
#     legend.key.width = unit(.4, "cm"),  # Set the legend key width
#     legend.key.height = unit(.4, "cm")  # Set the legend key height
#   ) +
#   labs(color = paste0('Average annual\n', current_name,  ' \nConcentration (ug/m3)'))
# # Specify the order for the legend
# current_pred_map

# 
#   # save the predictions
# ggsave(filename = paste0("Fig5_gridded_preds_by20", current_name, ".pdf"),
#          plot = current_pred_map,
#          path = file.path(results_fp, 'figures/Fig5'),
#          scale = 1,
#          width = 7,
#          height = 6,
#          dpi = 320)



# get geometry for predictions within a region
# geom_sf <- matched_grid_smoke %>%
#   filter(region == current_region) %>%
#   dplyr::select(region, grid_id_10km) %>%
#   group_by(region, grid_id_10km) %>%
#   mutate(row = row_number()) %>%
#   ungroup() %>%
#   filter(row == 1)


# # current_pred_map <- ggplot() +
# #   # geom_sf(data = us_sf, fill = "grey95") +  # Add the base map
# #   geom_sf(data = sample_avg_predictions, aes(color = avg_pred_conc), alpha =.70) +  # Add point data
# #   scale_color_viridis(option = "A", direction = -1)+
# #   theme_minimal() +  # Apply a minimal theme
# #   theme(
# #     panel.grid.major = element_blank(),  # Remove major grid lines
# #     panel.grid.minor = element_blank(),  # Remove minor grid lines
# #     axis.text = element_blank(),  # Remove axis labels
# #     axis.title = element_text(),  # Remove axis titles
# #     legend.title = element_text(),
# #     plot.title = element_text(face = "bold"),
# #     legend.key.width = unit(.4, "cm"),  # Set the legend key width
# #     legend.key.height = unit(.4, "cm")  # Set the legend key height
# #   ) +
# #   labs(color = paste0('Average annual\n', current_name,  ' \nConcentration (ug/m3)'))
# 
# #pred_map
# # ggsave(filename = paste0("Fig5_gridded_preds_", unique(current_species_df$species), ".pdf"),
# #        plot = current_pred_map,
# #        path = file.path(results_fp, 'figures/Fig5'),
# #        scale = 1,
# #        width = 7,
# #        height = 6,
# #        dpi = 320)
# #
# # sample_avg_predictions
# #
# #
# # # set limits to order the plot
# # limits <- c(
# #   # "Organics"
# #   "Organic Carbon (OC)", "Elemental Carbon (EC)",
# #   # "Halogens"
# #   "Bromine (Br)", "Chlorine (Cl)", "Chloride (Chl)",
# #   #  "Nonmetals"
# #   "Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Phosphorus (P)", "Selenium (Se)",
# #   # "Other metals"
# #   "Titanium (Ti)", "Aluminum (Al)", "Lead (Pb)",
# #   # "Metalloids"
# #   "Silicon (Si)", "Arsenic (As)",
# #   # "Transition metals"
# #   "Zinc (Zn)", "Manganese (Mn)",  "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
# #   # "Alkali metals"
# #   "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
# #   # "Alkaline-earth metals"
# #   "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")
# 
# # plot all predictions in one map big map faceted by species
# # all_plots <- all_gridded_species_conc_pred_df %>%
# #   split(.$species) %>%
# #   map(~ ggplot(., aes(fill = avg_pred_conc), alpha =.70) +
# #         geom_sf() +
# #         facet_wrap(~species) +
# #         scale_fill_viridis(option = "A", direction = -1) +
# #         theme(
# #           panel.grid.major = element_blank(),  # Remove major grid lines
# #           panel.grid.minor = element_blank(),  # Remove minor grid lines
# #           axis.text = element_blank(),  # Remove axis labels
# #           axis.title = element_text(),  # Remove axis titles
# #           legend.title = element_text(),
# #           plot.title = element_text(face = "bold"),
# #           legend.key.width = unit(.4, "cm"),  # Set the legend key width
# #           legend.key.height = unit(.4, "cm")  # Set the legend key height
# #         ) +
# #         theme_minimal()) %>%
# #   cowplot::plot_grid(plotlist = .)
# #
# # all_plots
# # all_pred_map <- ggplot() +
# #   # geom_sf(data = us_sf, fill = "grey95") +  # Add the base map
# #   geom_sf(data = all_gridded_species_conc_pred_df, aes(color = avg_pred_conc), alpha =.70) +  # Add point data
# #   scale_color_viridis(option = "A", direction = -1) +
# #   facet_wrap(~species, scales = 'free', ncol = 6) +
# #   # facet_wrap(~factor(species, levels = rev(limits)), ncol = 4) +
# #   theme_minimal() +  # Apply a minimal theme
# #   theme(
# #     panel.grid.major = element_blank(),  # Remove major grid lines
# #     panel.grid.minor = element_blank(),  # Remove minor grid lines
# #     axis.text = element_blank(),  # Remove axis labels
# #     axis.title = element_text(),  # Remove axis titles
# #     legend.title = element_text(),
# #     plot.title = element_text(face = "bold"),
# #     legend.key.width = unit(.4, "cm"),  # Set the legend key width
# #     legend.key.height = unit(.4, "cm")  # Set the legend key height
# #   ) +
# #   labs(color = 'Concentration (ug/m3)')
# 
# 
# 
# # save
# # ggsave(filename = paste0("Fig5_gridded_preds_faceted.pdf"),
# #        plot = all_plots,
# #        path = file.path(results_fp, 'figures/Fig5'),
# #        scale = 1,
# #        width = 12,
# #        height = 14,
# #        dpi = 320)
# 
# 
# # 
# rm(predicted_gridded_species_conc_df)# drop to save memory
# # get quantiles
# binned_predictions <- tibble(
#   bins = cut(all_regions_pred_sf$avg_pred_conc,
#              breaks = quantile(all_regions_pred_sf$avg_pred_conc,
#                                probs = seq(0,1,.20),
#                                na.rm = TRUE),
#              include.lowest = T, dig.lab = 5)
# ) %>%
#   distinct() %>%
#   arrange(bins) %>%
#   bind_cols(
#     tibble(
#       percentile = names(quantile(all_regions_pred_sf$avg_pred_conc,
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
# gridded_preds_w_bins <- all_regions_pred_sf %>%
#   mutate(bins = cut(all_regions_pred_sf$avg_pred_conc,
#                     breaks = quantile(all_regions_pred_sf$avg_pred_conc,
#                                       probs = seq(0,1,.20),
#                                       na.rm = TRUE),
#                     include.lowest = T, dig.lab = 5)) %>%
#   left_join(binned_predictions_ordered, by = "bins") %>%
#   dplyr::select(-bins, -start, -end)
