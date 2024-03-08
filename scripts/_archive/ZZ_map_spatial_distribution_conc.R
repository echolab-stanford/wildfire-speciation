# # Emma Krasovich Southworth, emmars@stanford.edu
# # Last Updated: January 25, 2024
# # Description: Here wetake our average predictions and plot

# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(reg_pal, cache = drake_cache)

map_spatial_distribution_conc <- function(grid_fp, reg_pal) {

  # Plan for parallel processing
  plan(multiprocess, workers = 4)

  # Read in grid + transform
  grid_10km <- st_read(grid_fp) %>% 
    dplyr::select(grid_id_10km = 'ID')
  
  species <- c('AS', 'NI', 'PB')
  # current_species <- species[3]
  # current_species <- 'PB'
   
  # # for each species:
  all_spec_avg_regional_period_df <- map_df(species, function(current_species) {
    
    # list files in the folder path + read in
    files <- list.files(file.path(data_fp, 'intermediate'), 
                        pattern = 'regional_gridded(PB|NI|AS)*',
                        full.name = TRUE)
    
    species_fps <- files[str_detect(files, current_species)]
    # current_conc_fp <- species_fps[3]
    
    # map over each file and read it in
    avg_regional_period_df <- map_df(species_fps, function(current_conc_fp) {
      
      print(current_conc_fp)
      # read in predicted concentration for every grid cell for every day in sample
      current_region_spec_df <- read.fst(current_conc_fp) %>%  # each file is one species and one region
      # calculate avg conc in grid cell for each time period, the first five, middle 5 and last 5 years
        mutate(year = year(Date)) %>% 
        mutate(samp_period = case_when(
          year > 2005 & year < 2011 ~ '2006-2010',
          year > 2010 & year < 2016 ~ '2011-2015',
          year > 2015 & year < 2021 ~ '2016-2020'
        )) %>% 
        # group by grid cell, region, species, and sample_period
        group_by(region, grid_id_10km, species, species_long, samp_period) %>% 
        summarise(avg_grid_conc_in_period = mean(pred_grid_conc, na.rm = TRUE),
                  max_conc_in_grid_cell = max(pred_grid_conc, na.rm = TRUE)) %>% 
        ungroup() 
      
  
    }) %>%  
      bind_rows()
    
    
    # MAKE A REGIONAL PLOT FOR EACH SAMPLE PERIOD -------------------------
    regional_time_period_map <- ggplot() +
      geom_sf(data = avg_regional_period_df %>%
                left_join(grid_10km, by ='grid_id_10km') %>% 
                st_as_sf(), 
              aes(color = avg_grid_conc_in_period, 
                  fill = avg_grid_conc_in_period)) +  
      scale_colour_viridis(option = "A", direction = -1) +
      scale_fill_viridis(option = "A", direction = -1) +
      facet_grid(region~samp_period) +
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
           color = paste0('')) +
      guides(fill = 'none')
    regional_time_period_map
    
    ggsave(
      filename = paste0('Fig5A_regional_avg_', current_species, '_conc.png'),
      plot = regional_time_period_map,
      path = file.path(results_fp, 'figures/Fig5'),
      scale = 1,
      width = 12,
      height = 16,
      dpi = 320)
    
  }) %>% 
    bind_rows() # for all species
      
      
      return(all_spec_avg_regional_period_df)
      
}
      

