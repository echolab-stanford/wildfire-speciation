# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: create a map of the different monitoring stations

# loadd(clean_PMspec_df, cache = drake_cache)
# us_states_fp = file.path(data_fp, 'raw/all_national_states.rds')

# create a map of CSN vs improve sites

create_map_of_monitoring_stations <- function(clean_PMspec_df, us_states_fp) {
  
  # get duration of each site:
  monitor_duration <- clean_PMspec_df %>% 
    mutate(Dataset = ifelse(Dataset == 'EPACSN', "CSN", "IMPROVE")) %>% 
    distinct(Dataset, site_id, year) %>% 
    group_by(site_id) %>% 
    # get how many years each site is online
    dplyr::summarise(yrs_online = n_distinct(year)) %>% 
    ungroup() %>% 
    # create categorical var for duration
    mutate(duration_cat = case_when(
      yrs_online < 5 ~ 'less than 5 years',
      yrs_online >= 5 & yrs_online < 10 ~ '5-9 years', 
      yrs_online >= 10 ~ '10-15 years')) %>% 
    mutate(duration_cat = factor(duration_cat, levels = c('less than 5 years',
                                                          '5-9 years',
                                                          '10-15 years')))
    #hist(monitor_duration$yrs_online)
  
  # create sf object of sites
  monitor_sf <- clean_PMspec_df %>% 
    distinct(Dataset, site_id, lat, long, epsg) %>% 
    mutate(Dataset = ifelse(Dataset == 'CSN', "CSN", "IMPROVE")) %>% 
    left_join(monitor_duration, by = 'site_id')

  # drop time frame and only include points
  tally <- monitor_sf %>% 
    st_drop_geometry() %>% 
    group_by(duration_cat) %>% 
    tally()
  
  # read in USA shapefile
  us_sf <- readRDS(us_states_fp) %>% 
    st_as_sf() %>% 
    st_transform(crs = 4326) %>% 
    filter(!STUSPS %in% c('AK', 'HI', "PR", 'VI', 'AS', 'GU', 'MP'))

  # create sf object for plotting
  monitor_cat_sf <- monitor_sf %>% 
    st_as_sf(coords = c('long', 'lat'), crs= 4326) %>% 
    left_join(tally) %>% 
    mutate(mon_dur_cat = paste0(duration_cat, " (n=",n, ")" ))
  
  # create map
  mon_map <- ggplot() +
    geom_sf(data = us_sf, fill = "grey99") +  # Add the base map
    geom_sf(data = monitor_cat_sf, aes(color = mon_dur_cat, 
                                   shape =Dataset), size = 2, alpha =.70) +  # Add point data
    scale_color_manual(values=c('plum4', 'thistle3', 'lavenderblush3')) + 
    scale_shape_manual(values=c(16, 17)) + 
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
    labs(title = "A    Monitoring station distribution", 
         color = 'Duration online in sample', 
         shape = 'Monitoring network')
  
  mon_map
  
  ggsave(file.path(results_fp, "Fig1/Fig1A_monitor_stn_map_raw.pdf"), 
         plot = mon_map, dpi = 320, width = 7, height = 6)
  
  
  # drop time frame and only include points
  mon_facet <- monitor_cat_sf %>% 
    dplyr::select(-mon_dur_cat) %>% 
    left_join(tally) %>% 
    mutate(duration_cat_ct = paste0(duration_cat, " (n=", n, ")")) %>% 
    mutate(duration_cat_ct = factor(duration_cat_ct, levels = c('less than 5 years (n=142)',
                                                          '5-9 years (n=104)',
                                                          '10-15 years (n=454)')))
  
  # create supplementary figure for revisions that reflect 
  mon_map_no_duration <- ggplot() +
    geom_sf(data = us_sf, fill = "grey95") +  # Add the base map
    geom_sf(data = mon_facet, 
            aes(color = Dataset), size = 2.5, alpha =.70) +  # Add point data
    scale_color_manual(values=c('black', 'firebrick')) + 
    facet_wrap(~duration_cat_ct) +
    theme_minimal() +  # Apply a minimal theme
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.text = element_blank(),  # Remove axis labels
      axis.title = element_text(),  # Remove axis titles
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = 'bottom',
      legend.key.width = unit(.4, "cm"),  # Set the legend key width
      legend.key.height = unit(.4, "cm")  # Set the legend key height
    ) +
    labs(color = 'Monitoring network')
  
  mon_map_no_duration
  
  ggsave(file.path(results_fp, "SI Figs/SIFig2_monitor_stn_map_facet.png"), 
         plot = mon_map_no_duration, dpi = 320, width = 7, height = 6)
  ggsave(file.path(results_fp, "SI Figs/SIFig2_monitor_stn_map_facet.pdf"), 
         plot = mon_map_no_duration, dpi = 320, width = 7, height = 6)
  

  
  return(mon_map)
}
  
