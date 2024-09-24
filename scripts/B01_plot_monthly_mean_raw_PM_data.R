# Ayako Kawano, Emma Krasovich Southworth
# Created: October 24, 2023
# Description: This function plots the raw smoke and total PM data for figure 1

# loadd(c(clean_PMspec_df, parameter_categories, pm_pal), cache = drake::drake_cache(".drake"))

# function to plot PM data
plot_monthly_mean_raw_PM_data <- function(clean_PMspec_df, parameter_categories, pm_pal) {
  
  # set up dataframe with species+date
  df <- clean_PMspec_df %>% 
    mutate(month = month(Date)) %>% 
    dplyr::select(site_id, year, month, totPM2.5 = 'MF_adj', smokePM, nonsmokePM = 'nonsmokePM_MF') %>% 
    pivot_longer(cols =c(totPM2.5, smokePM, nonsmokePM), 
                 names_to = 'species', values_to = 'conc')
  
# -------------------------------------------------
# NOW PLOT TOTAL PM and SMOKE PM
# -------------------------------------------------
#plot all region, all years, monthly means for PM2.5 ----
pm25 <- df %>% 
  group_by(site_id, year, month, species) %>% 
  dplyr::summarise(avg_site_mon_conc = mean(conc, na.rm = TRUE)) %>% 
  ungroup()  %>%
  # get monthly average for a given species across all sites
  group_by(year, month, species) %>%
  dplyr::summarise(avg_mon_conc = mean(avg_site_mon_conc, na.rm = TRUE)) %>%
  ungroup()  %>% 
  group_by(year, month, species) %>%
  dplyr::summarise(avg_mon_conc = mean(avg_mon_conc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(species = ifelse(species == 'smokePM', "Wildfire smoke PM2.5", "Total PM2.5")) %>% 
  mutate(species = fct_relevel(species, c( "Total PM2.5", "Wildfire smoke PM2.5"))) %>% 
  mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) 


data_ends <- pm25 %>% filter(month == 12)
datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-01-01"), by = "24 month")


avg_mon_PM_plot <- ggplot(pm25) +
  geom_density(aes(x = mon_yr, 
                   y = avg_mon_conc, 
                   fill =species, 
                   color = species), stat = "identity", alpha = 0.8) +  
  labs(title = expression("B    Total PM"[2.5]*" and wildfire smoke PM"[2.5]*" in sample"),
                          y = expression("Concentration (\u00b5g/m"^3*")"),
                          x = "",
                          ) +
  scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
  scale_y_continuous(limits = c(0,12) , breaks=c(0, 2, 4, 6, 8, 10, 12),
                     labels=c("0","2.0", "4.0", "6.0", "8.0", "10.0", "12.0")) +
  geom_hline(aes(yintercept = 0)) +
  theme_minimal() +
  scale_color_manual(values = pm_pal) +
  scale_fill_manual(values = pm_pal) +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0),
        axis.text = element_text(),
        axis.title = element_text(),
        # axis.text.x = element_text(vjust = 1),
        legend.text = element_text(),
        legend.title = element_blank(),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.grid = element_blank(),
        legend.position = c(0, 1),  # Adjust the legend position (top-left corner)
        legend.justification = c(0, 1))  # Justify the legend to the top-left) +

avg_mon_PM_plot

ggsave(
  filename = 'Fig1B_avg_monthlyPM_raw.pdf',
  plot = avg_mon_PM_plot,
  path = file.path(results_fp, 'Fig1'),
  scale = 1,
  width = 8,
  height = 6,
  dpi = 320)

return(pm25)
} # end function
