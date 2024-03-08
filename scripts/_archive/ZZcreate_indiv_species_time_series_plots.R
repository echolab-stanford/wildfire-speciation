# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species

# loadd(clean_pm_spec_df, cache = drake_cache)

create_indiv_species_time_series_plots <- function(clean_pm_spec_df) {

# set up dataframe 
species_conc_long_df <- clean_pm_spec_df %>% 
  dplyr::select(-c(units, Elevation, Latitude, Longitude, epsg, km2fire)) %>% 
  rename(nonsmokePM = 'non_smokePM_cont',
         smokePM = 'smokePM_pred',
         totPM2.5 = 'PM2.5') %>% 
  pivot_longer(cols = totPM2.5:ZR,
               names_to = 'species_name',
               values_to = 'conc_val') %>% 
  group_by(Dataset, SiteCode, Date, State, season, region, year, 
           month, doy, species_name, smoke_day) %>% 
  dplyr::summarise(conc_val = mean(conc_val, na.rm=TRUE)) %>% 
  ungroup()


# map over each species and generate a time series for that species
# first get list of all species:
spec_list <- unique(species_conc_long_df$species_name)
# current_species <- spec_list[7] # testing

# for each species, generate a plot 
all_month_yr_avg_df <- purrr::map_df(spec_list, function(current_species) {
  
  # print current iteration
  print(paste('current species:', current_species, ' | current index:', 
              match(current_species, spec_list), 
              'out of', length(spec_list)))
  
  current_spec_df <- species_conc_long_df %>% 
    filter(species_name == current_species) %>% 
    filter(!is.na(conc_val)) %>% 
    # get avg concentration across all sites for a given month-year
    group_by(year, month, species_name) %>% 
    dplyr::summarise(avg_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop') %>% 
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d"))
    
  # make time series plot
  current_plot <- ggplot(data = current_spec_df, 
                         aes(x = mon_yr, y = avg_conc)) +
    geom_line(color = "grey60", size = .5) +
    geom_smooth(method=lm, se=FALSE, col='maroon', size=1) +
    scale_x_date(labels = date_format("%b-%Y"), minor_breaks = "year") +
    labs(x = 'Date',
         y = 'Concentration (ug/m3)',
         title = paste(current_species, "time series from 2006-2020")
         ) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    theme_light()

  # current_plot
  
  # save plot
  ggsave(file.path(wip_gdrive_fp, 'species_time_series/figures', paste(current_species, 'time_series.png')), 
         width=6, height=4, dpi=300)
 
  current_spec_df
}) %>% # create a df of all the summarized concentrations
  bind_rows()


# DIVIDE THROUGH BY THE AVERAGE:
conc_rel2avg_df <- species_conc_long_df %>%
  mutate(SiteCode = as.factor(SiteCode),
         year = as.factor(year),
         month = as.factor(month),
         moy = as.factor(paste0(month,'-', year)),
         monitor_month = paste0(SiteCode, "_", month)) %>% 
  group_by(species_name) %>% 
  mutate(baseline_avg = mean(conc_val, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(CR_val = ifelse(!is.na(conc_val), conc_val/baseline_avg, NA)) %>% 
  pivot_wider(id_cols = c('Dataset', 'SiteCode',  
                          'Date', 'State', 'smoke_day',
                          'season', 'region', 'year',
                          'month', 'doy', 'moy', 'monitor_month'),
    names_from = 'species_name', 
    values_from = c('conc_val', 'CR_val'))

return(conc_rel2avg_df)

} # end function




