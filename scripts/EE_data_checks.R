# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: November 6, 2023
# Description: Check individual plots of monitors for weird chemicals

loadd(c(clean_PMspec_df, parameter_categories, spec_pal), cache = drake::drake_cache(".drake"))

check_weird_vars <- clean_PMspec_df %>% 
  filter(duration == 15) %>% 
  #filter(year >2010 & year < 2019) %>% 
  dplyr::select(year, month, Date, site_id, RB, PB, ZR, SE) %>% #SE, PB
  pivot_longer(cols =c(PB, RB, ZR, SE), 
               names_to = 'species', values_to = 'conc') 


# calculate how many observations go into a monthly average at a given site
# obs_in_mon_avg <- check_weird_vars %>% 
#   group_by(species, year, month, site_id) %>% 
#   tally() %>% 
#   rename(num_obs_at_site_per_month = 'n')

#plot histogram by species
# num_obs_per_month_plot <- ggplot(obs_in_mon_avg, aes(x = num_obs_at_site_per_month)) +
#   geom_histogram(fill = 'grey80', color = "black") +
#   facet_wrap(~ species, ncol = 2) +  # Facet by the "group" variable with 2 columns
#   labs(title = "# of Obs in Monthly Site Avg", y = "Frequency") +
#   guides(color = 'none') + theme_minimal()
# num_obs_per_month_plot


# calculate how mnay sites are in each month for each species
monthly_mean_species_conc <- check_weird_vars %>% 
  group_by(species, year, month, site_id) %>% 
  dplyr::summarise(avg_conc = mean(conc, na.rm=T)) %>% 
  ungroup() %>% 
  group_by(species, year, month) %>%
  dplyr::summarise(avg_conc = mean(avg_conc, na.rm=T)) %>% 
  ungroup()  %>% 
  mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) 

data_ends <- monthly_mean_species_conc %>% 
  filter(mon_yr == max(mon_yr))

datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")


plot_monthly_means <- ggplot(monthly_mean_species_conc, 
                             aes(x = as.Date(mon_yr), y = avg_conc, group = species)) +
  geom_line(size = 0.8, aes(color = species)) +
  ylab("Concentration (ug/m3)") +
  xlab("") +
  ggtitle("C.   Monthly average concentrations of chemical species, from 2006 to 2020") +
  scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
  #scale_y_continuous(trans='log10') +
  # limits = c(10^-5,4.5), 
  # breaks=c(10^-6, 10^-5, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 4.5),
  # labels=c("0", '0.00001', "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1.0", "2.0", "4.5")) +
  theme_minimal() +
  geom_point(data = data_ends, size =1.5,
             aes(x = as.Date(mon_yr), y = avg_conc, color = species), show.legend = FALSE)+
  geom_text(data = data_ends,
            aes(x =as.Date("2021-01-01", format = "%Y"), 
                y = avg_conc,  
                label = species), 
            nudge_x = 0.2, show.legend = FALSE, 
            hjust = 0, size = 3)+
  facet_wrap(~species, ncol = 2, scales = 'free_y') +
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




#plot histogram by species
# sites_in_mon_avg_plot <- ggplot(sites_in_mon_avg, aes(x = num_sites_in_monthly_avg)) +
#   geom_histogram(fill = 'grey80', color = "black") +
#   facet_wrap(~species, ncol = 2) +  # Facet by the "group" variable with 2 columns
#   labs(title = "# of Sites in Monthly Site Avg", y = "Frequency") +
#   guides(color = 'none') + theme_minimal()
# sites_in_mon_avg_plot
# 
# sites_in_mon_avg <- check_weird_vars %>% 
#   group_by(species, year, month, site_id) %>% 
#   dplyr::summarise(avg_conc = mean(conc, na.rm=T)) %>% 
#   ungroup() 

# plot_boxplots <- ggplot(sites_in_mon_avg, 
#                            aes(x = as.factor(year), y = avg_conc, fill = species)) +
#   geom_boxplot() +
#   ylab("Concentration (ug/m3)") +
#   xlab("") +
#   facet_wrap(~species, scales = 'free_y') +
#   theme_minimal()
# plot_boxplots



  
# get random sample of sites
sample_sites <- sample(check_weird_vars$site_id, 30)

check_weird_vars_sample <- check_weird_vars %>% 
  filter(site_id %in% sample_sites)


plot_individual_monitors <- ggplot(check_weird_vars_sample, 
                             aes(x = Date, y = conc, group = species)) +
  geom_line(size = 0.8, aes(color = species)) +
  ylab("Concentration (ug/m3)") +
  xlab("") +
  theme_minimal() +
  facet_wrap(~site_id, scales = 'free_y') +
  scale_color_manual(values = spec_pal, na.translate = F) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(vjust = 0),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray95"))
plot_individual_monitors
