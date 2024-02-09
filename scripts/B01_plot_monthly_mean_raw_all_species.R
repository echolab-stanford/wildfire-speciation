# Ayako Kawano, Emma Krasovich Southworth
# Description: This function plots the raw data

# loadd(c(clean_PMspec_df, parameter_categories, spec_pal), cache = drake::drake_cache(".drake"))

# create function to plot raw averages for each species
plot_monthly_mean_raw_all_specices <- function(clean_PMspec_df, parameter_categories, spec_pal) {
  
# set up dataframe with species+date
avg_df <- clean_PMspec_df %>% 
  dplyr::select(-SOIL, -ZR) %>% 
  dplyr::select(Dataset, site_id, year, Date, month, AL:ZN) %>% 
  pivot_longer(cols =c(AL:ZN), 
               names_to = 'species', values_to = 'conc') %>% 
  # get monthly averages for a given site + species within a year
  group_by(site_id, year, month, species) %>% 
  dplyr::summarise(avg_site_mon_conc = mean(conc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # get monthly average for a given species across all sites
  group_by(year, month, species) %>%
  dplyr::summarise(avg_mon_conc = mean(avg_site_mon_conc, na.rm = TRUE)) %>% 
  ungroup() 

# clean up dataframe
avg_w_spec_type <- avg_df %>% 
  # filter(Dataset == 'CSN') %>% 
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
  mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
  group_by(species_type) %>% 
  mutate(max_lim = max(avg_mon_conc),
         min_lim = min(avg_mon_conc)) %>% 
  ungroup() %>% 
  mutate(color = case_when(
    species_type == 'Alkaline-earth metals' ~ "#35274A",
    species_type =='Alkali metals' ~ "#B40F20",
    species_type =='Transition metals'~ "#E58601",
    species_type =='Metalloids' ~ "#E2D200",
    species_type =='Other metals'~ "#46ACC8",
    species_type =='Nonmetals' ~ "#0B775E",
    species_type =='Halogens' ~ "#CAB38C",
    species_type =='Organics' ~ "#899DA4"
  ))


datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")
type_list <- unique(avg_w_spec_type$species_type)

# map over eeach species type and create fig for each
all_indiv_species_monthly_plots <- purrr:::map(type_list, function(current_type) {
# current_type <- type_list[5]

current_plot_df <- avg_w_spec_type %>% 
  filter(species_type == current_type)

data_ends <- current_plot_df %>% 
  filter(mon_yr == max(mon_yr))

min <- unique(current_plot_df$min_lim)
max <- unique(current_plot_df$max_lim)
current_color <- unique(current_plot_df$color)

# -------------------------------------------------
# NOW PLOT IN LOGS 
# -------------------------------------------------
plot_monthly_means <- ggplot(current_plot_df, 
                             aes(x = as.Date(mon_yr), y = avg_mon_conc, group = species_long)) +
  geom_line(size = 0.6, color = current_color) + #aes(color = current_color)) +
  ylab(expression("Concentration (\u00b5g/m"^3*")")) +
  xlab("") +
  ggtitle(paste0(unique(current_plot_df$species_type))) +
  scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
  scale_y_continuous(trans='log10',
                     limits = c(min, max)) +
  geom_point(data = data_ends, size = 1,
             aes(x = as.Date(mon_yr), y = avg_mon_conc), color = current_color, show.legend = FALSE) +
  geom_text(data = data_ends,
            aes(x =as.Date("2021-01-01", format = "%Y"), 
                y = avg_mon_conc,  
                label = species), 
            nudge_x = 0.2, show.legend = FALSE, 
            hjust = 1, size = 3)+
  geom_vline(xintercept = as.Date("2011-01-01"), linetype = 'dotted', color = "grey50") +  # method change 2011
  geom_vline(xintercept = as.Date("2015-12-01"), linetype = 'dotted', color = "grey20") +  # method change nov 2015
  geom_vline(xintercept = as.Date("2018-10-01"), linetype = 'dotted', color = "grey20") + 
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(vjust = 0),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray95")) +
  geom_hline(aes(yintercept = 0)) + 
  guides(color = 'none')
plot_monthly_means


ggsave(
  filename = paste0('Fig1C_LOG_spec_monthly_means_', unique(current_plot_df$species_type), '.pdf'),
  plot = plot_monthly_means,
  path = file.path(results_fp, 'figures/Fig1'),
  scale = 1,
  width = 5,
  height = 4,
  dpi = 320)

current_plot_df
})  # end loop



# -------------------------------------------------
# ALL DATA
# -------------------------------------------------
data_ends <- avg_w_spec_type %>% 
  filter(mon_yr == max(mon_yr))

plot_all_monthly_means <- ggplot(avg_w_spec_type, 
                             aes(x = as.Date(mon_yr), y = avg_mon_conc, group = species_long)) +
  geom_line(size = 0.6, aes(color = species_type)) +
  ylab(expression("Concentration (\u00b5g/m"^3*")")) +
  xlab("") +
  ggtitle("C.   Monthly average concentrations of chemical species, from 2006 to 2020") +
  scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
  scale_y_continuous(trans='log10',
                     limits = c(0.0001, 4.1), 
                     breaks=c(10^-5, 0.0001, 0.001, 0.01, 0.1, 1.0, 2.0, 4.1), 
                     labels=c("0", '0.0001', "0.001", "0.01", "0.1", "1.0", "2.0", "4")) +
  geom_vline(xintercept = as.Date("2011-01-01"), linetype = 'dotted', color = "grey50") +  # method change 2011
  geom_vline(xintercept = as.Date("2015-12-01"), linetype = 'dotted', color = "grey20") +  # method change nov 2015
  geom_vline(xintercept = as.Date("2018-10-01"), linetype = 'dotted', color = "grey20") + 
  geom_point(data = data_ends, size = 1,
             aes(x = as.Date(mon_yr), y = avg_mon_conc, color = species_type), show.legend = FALSE)+
  geom_text(data = data_ends,
            aes(x =as.Date("2021-01-01", format = "%Y"), 
                y = avg_mon_conc,  
                label = species), 
            nudge_x = 0.2, show.legend = FALSE, 
            hjust = 0, size = 3)+
  facet_wrap(~species_type, ncol = 4, scales = 'free_y') +
  scale_color_manual(values = spec_pal, na.translate = F) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(vjust = 0),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = "gray95")) +
  geom_hline(aes(yintercept = 0)) + 
  guides(color = 'none')
plot_all_monthly_means


ggsave(
  filename = 'Fig1C_LOG_all_spec_monthly_means_raw.pdf',
  plot = plot_all_monthly_means,
  path = file.path(results_fp, 'figures/Fig1'),
  scale = 1,
  width = 16,
  height = 8,
  dpi = 320)

# -------------------------------------------------------------
# PLOT MONTHLY MEANS NORMALIZED BY ANNUAL MEAN
# -------------------------------------------------------------
# get monthly means
# monthly_avg <- clean_PMspec_df %>% 
#   dplyr::select(-SOIL, -ZR) %>% 
#   dplyr::select(Dataset, site_id, year, Date, month, AL:ZN) %>% 
#   pivot_longer(cols =c(AL:ZN), 
#                names_to = 'species', values_to = 'conc') %>% 
#   # get monthly averages for a given site + species within a year
#   group_by(site_id, year, month, species) %>% 
#   dplyr::summarise(avg_site_mon_conc = mean(conc, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   # get monthly average for a given species across all sites
#   group_by(year, month, species) %>%
#   dplyr::summarise(avg_mon_conc = mean(avg_site_mon_conc, na.rm = TRUE)) %>% 
#   ungroup() 
# 
# # get annual means
# annual_avg <- clean_PMspec_df %>% 
#   dplyr::select(-SOIL, -ZR) %>% 
#   dplyr::select(Dataset, site_id, year, Date, month, AL:ZN) %>% 
#   pivot_longer(cols =c(AL:ZN), 
#                names_to = 'species', values_to = 'conc') %>% 
#   # get monthly averages for a given site + species within a year
#   group_by(site_id, year, species) %>% 
#   dplyr::summarise(avg_site_yr_conc = mean(conc, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   # get monthly average for a given species across all sites
#   group_by(year, species) %>%
#   dplyr::summarise(avg_yr_conc = mean(avg_site_yr_conc, na.rm = TRUE)) %>% 
#   ungroup() 
# 
# # join monthly and annual means:
# norm_mon_means <- monthly_avg %>% 
#   left_join(annual_avg, by = c('year', 'species')) %>% 
#   mutate(norm_conc = avg_mon_conc/avg_yr_conc) %>% 
#   
#   left_join(parameter_categories %>% 
#               mutate(species = str_to_upper(species)) %>% 
#               dplyr::select(species, species_long, 
#                             species_type), by = 'species') %>% 
#   filter(!is.na(species_type)) %>% 
#   mutate(species_type = fct_relevel(species_type,
#                                     c("Alkaline-earth metals", "Alkali metals", 
#                                       "Transition metals", "Metalloids", "Other metals", 
#                                       "Nonmetals",  "Halogens", "Organics"))
#   ) %>% 
#   mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
#   mutate(color = case_when(
#     species_type == 'Alkaline-earth metals' ~ "#35274A",
#     species_type =='Alkali metals' ~ "#B40F20",
#     species_type =='Transition metals'~ "#E58601",
#     species_type =='Metalloids' ~ "#E2D200",
#     species_type =='Other metals'~ "#46ACC8",
#     species_type =='Nonmetals' ~ "#0B775E",
#     species_type =='Halogens' ~ "#CAB38C",
#     species_type =='Organics' ~ "#899DA4"
#   ))
# 
# data_ends <- norm_mon_means %>% 
#   filter(mon_yr == max(mon_yr))
# 
# plot_monthly_means <- ggplot(norm_mon_means, 
#                              aes(x = as.Date(mon_yr), y = norm_conc, group = species_long)) +
#   geom_line(size = 0.6, aes(color = species_type)) +
#   ylab("Normalized Concentration (ug/m3)") +
#   xlab("") +
#   ggtitle("C.   Monthly mean species concentration normalized by annual mean, from 2006 to 2020") +
#   scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
#   geom_point(data = data_ends, size = 1,
#              aes(x = as.Date(mon_yr), y = norm_conc, color = species_type), show.legend = FALSE)+
#   geom_text(data = data_ends,
#             aes(x =as.Date("2021-01-01", format = "%Y"), 
#                 y = norm_conc,  
#                 label = species), 
#             nudge_x = 0.2, show.legend = FALSE, 
#             hjust = 0, size = 3)+
#   facet_wrap(~species_type, ncol = 2, scales = 'free_y') +
#   scale_color_manual(values = spec_pal, na.translate = F) +
#   theme_minimal() +
#   theme(plot.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 9),
#         axis.title = element_text(size = 12),
#         axis.text.x = element_text(vjust = 0),
#         axis.ticks.x = element_line(linetype = 1),
#         axis.ticks.length=unit(-0.1, "cm"),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major.y = element_line(color = "gray95")) +
#   geom_hline(aes(yintercept = 10^-6)) + 
#   guides(color = 'none')
# plot_monthly_means

return(avg_w_spec_type)
} # end function

