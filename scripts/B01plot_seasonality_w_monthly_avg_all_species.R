# Ayako Kawano, Emma Krasovich Southworth
# Description: This function plots the raw data

#parameter_categories <- read.csv('/Users/ekrasovich/Desktop/xwalk_species.csv')
# loadd(c(clean_PMspec_df, spec_pal), cache = drake::drake_cache(".drake"))

# create function to plot raw averages for each species
plot_seasonality_w_monthly_average <- function(clean_PMspec_df, parameter_categories, spec_pal) {

# set up dataframe with species+date
avg_df <- clean_PMspec_df %>% 
  dplyr::select(site_id, year, Date, month, AL:NO3, OC:SO4, SR:ZR) %>% 
  pivot_longer(cols =c(AL:ZR), 
               names_to = 'species', values_to = 'conc') %>% 
  # get monthly averages for a given site + species within a year
  group_by(site_id, year, month, species) %>% 
  # count dates
  dplyr::summarise(distinct_dates_count = n_distinct(Date)) %>%
  ungroup() %>%
  group_by(site_id, year, month, species) %>% 
  dplyr::summarise(avg_site_mon_conc = mean(conc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # get monthly average for a given species across all sites
  group_by(year, month, species, distinct_dates_count) %>%
  dplyr::summarise(num_sites = n_distinct(site_id)) %>% 
  ungroup() %>% 
  group_by(year, month, species, distinct_dates_count) %>%
  dplyr::summarise(avg_mon_conc = mean(avg_site_mon_conc, na.rm = TRUE)) %>% 
  ungroup() 

# clean up dataframe
avg_w_spec_type <- avg_df %>% 
  left_join(parameter_categories %>% 
              mutate(species = str_to_upper(species)) %>% 
              dplyr::select(species, species_long, 
                            species_type), by = 'species') %>% 
  filter(!is.na(species_type)) %>% 
  mutate(species_type = fct_relevel(species_type,
                                 c("Alkaline-earth metals", "Alkali metals", 
                                   "Transition metals", "Metalloids", "Other metals", 
                                   "Nonmetals",  "Halogens", "Organics"))
                                 ) %>% 
  mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) 

data_ends <- avg_w_spec_type %>% 
  filter(mon_yr == max(mon_yr))

datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")

# -------------------------------------------------
# NOW PLOT IN LOGS 
# -------------------------------------------------
plot_monthly_means <- ggplot(avg_w_spec_type, 
                             aes(x = as.Date(mon_yr), y = avg_mon_conc, group = species_long)) +
  geom_line(size = 0.8, aes(color = species_type)) +
  ylab("Concentration (ug/m3)") +
  xlab("") +
  ggtitle("C.   Monthly average concentrations of chemical species, from 2006 to 2020") +
  scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
  scale_y_continuous(trans='log10') +
                     # limits = c(10^-5,4.5), 
                     # breaks=c(10^-6, 10^-5, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 4.5),
                     # labels=c("0", '0.00001', "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1.0", "2.0", "4.5")) +
  theme_minimal() +
  geom_point(data = data_ends, size =1.5,
             aes(x = as.Date(mon_yr), y = avg_mon_conc, color = species_type), show.legend = FALSE)+
  geom_text(data = data_ends,
            aes(x =as.Date("2021-01-01", format = "%Y"), 
                y = avg_mon_conc,  
                label = species), 
            nudge_x = 0.2, show.legend = FALSE, 
            hjust = 0, size = 3)+
  facet_wrap(~species_type, ncol = 2, scales = 'free_y') +
  scale_color_manual(values = spec_pal, na.translate = F) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(vjust = 0),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray95")) +
  geom_hline(aes(yintercept = 10^-6)) + 
  guides(color = 'none')
plot_monthly_means

ggsave(
  filename = 'Fig1C_LOG_spec_monthly_means.pdf',
  plot = plot_monthly_means,
  path = file.path(wip_gdrive_fp, 'figures/Fig1'),
  scale = 1,
  width = 10,
  height = 14,
  dpi = 320)


return(avg_w_spec_type)
} # end function

