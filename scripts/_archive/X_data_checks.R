# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: November 6, 2023
# Description: Check individual plots of monitors for weird chemicals

loadd(c(clean_PMspec_df, parameter_categories, spec_pal), cache = drake::drake_cache(".drake"))


# 
# plot <- ggplot(clean_PMspec_df, aes(x = MF, y = RCFM)) +
#   geom_point() +
#   geom_abline(slope = 1, size = 0.5, color = 'red')
# plot

check_weird_vars <- clean_PMspec_df %>% 
  filter(duration == 15) %>% 
  #filter(year >2010 & year < 2019) %>% 
  dplyr::select(Dataset, year, month, Date, site_id, AL:ZR) %>% #SE, PB
  pivot_longer(cols =c(AL:ZR), 
               names_to = 'species', values_to = 'conc') #%>% 
  #filter(species %in% c('RB', 'ZR', 'PB', 'SE', 'OC'))



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
  mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
  left_join(parameter_categories %>% 
              mutate(species = str_to_upper(species)) %>% 
              dplyr::select(-species_long, -X), by = 'species')

# get quantiles:
# q10 <- quantile(monthly_mean_species_conc$avg_conc, .10, na.rm=TRUE)
# q90 <- quantile(monthly_mean_species_conc$avg_conc, .90, na.rm=TRUE)
# 
# quantiles <- monthly_mean_species_conc %>% 
#   # create flag for quantile:
#   mutate(quantiles = case_when(
#    avg_conc <= q10 ~ 'bottom 10',
#    avg_conc >= q90 ~ 'top 10', 
#    TRUE ~ 'middle 80'))
#    
# plot_df <- monthly_mean_species_conc %>% 
#   filter(species_type %in% c('Alkali metals', 'Alkaline-earth metals', 'Organics', 'Metalloid'))
         # 
         # [1] "Other metals"          "Metalloid"             "Halogens"              "Alkaline-earth metals"
         # [5] "Transition metals"     "Organics"              "Alkali metals"         ""                     
         # [9] "Nonmetals"             "Nonmetal" 


datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")


plot_monthly_means <- ggplot(monthly_mean_species_conc,
                             aes(x = as.Date(mon_yr), y = avg_conc)) + #group = site_id
  geom_line(size = 0.8) +#, aes(color = Dataset)) +
  ylab("Concentration (ug/m3)") +
  xlab("") +
  ggtitle("C.   Monthly average concentrations of chemical species, from 2006 to 2020") +
  scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
  theme_minimal() +
  facet_wrap(~species, scales = 'free_y') +
  #scale_color_manual(values = spec_pal, na.translate = F) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(vjust = 0, angle = 90),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray95")) +
  geom_hline(aes(yintercept = 10^-6))
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
