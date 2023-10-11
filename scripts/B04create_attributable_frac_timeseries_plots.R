# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species


# loadd(clean_PMspec_df, cache = drake::drake_cache("scripts/.drake"))

# function
create_attributable_frac_time_series <- function(clean_PMspec_df) {

################################################################################
# How much is concentration due to smoke changing over time?
################################################################################
selected_spec_df <- clean_PMspec_df %>% 
  mutate(monitor_month = paste0(site_id,"_",month)) %>% 
  dplyr::select(Dataset:month, monitor_month,Date, site_id, MF_adj, smokePM, nonsmokePM_MF, AL:ammSO4) #%>% 
  # only keep rows that have all values so that its the same selection of sites for the analysis
  #na.omit()  

# run the regressions ----------------------------------------------------------
V_reg = feols(V ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
PB_reg = feols(PB ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
FE_reg = feols(FE ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
CU_reg = feols(CU ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
CR_reg = feols(CR ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
MN_reg = feols(MN ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
ZN_reg = feols(ZN ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
NI_reg = feols(NI ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
S_reg = feols(S ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
EC_reg = feols(EC ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')
OC_reg = feols(OC ~ smokePM + nonsmokePM_MF | monitor_month + year, selected_spec_df, cluster = 'site_id')


# NOW USE THE RESULTS OF EACH REGRESSION TO PREDICT
nonsmoke0_df <- selected_spec_df %>%
  # 0 out nonsmoke
  mutate(nonsmokePM = 0) 

smoke0_df <- selected_spec_df %>%
  # 0 out nonsmoke
  mutate(smokePM = 0)

# now predict with new data
# predict the smoke concentration
predicted_smoke_conc_df <- selected_spec_df %>% 
  # this gets the fraction of a species attributable to smoke
  mutate(V = predict(V_reg, newdata = nonsmoke0_df)) %>% 
  mutate(PB = predict(PB_reg, newdata = nonsmoke0_df)) %>%
  mutate(FE = predict(FE_reg, newdata = nonsmoke0_df)) %>%
  mutate(CU = predict(CU_reg, newdata = nonsmoke0_df)) %>%
  mutate(CR = predict(CR_reg, newdata = nonsmoke0_df)) %>%
  mutate(MN = predict(MN_reg, newdata = nonsmoke0_df)) %>%
  mutate(ZN = predict(ZN_reg, newdata = nonsmoke0_df)) %>%
  mutate(NI = predict(NI_reg, newdata = nonsmoke0_df)) %>%
  mutate(S = predict(S_reg, newdata = nonsmoke0_df)) %>%
  mutate(EC = predict(EC_reg, newdata = nonsmoke0_df)) %>%
  mutate(OC = predict(OC_reg, newdata = nonsmoke0_df)) %>% 
  pivot_longer(cols = AL:ammSO4, 
               names_to = 'species', 
               values_to = 'smoke_conc') 
 
# now predict with new data
predicted_nonsmoke_conc_df <- selected_spec_df %>% 
  # this gets the fraction of a species attributable to nonsmoke
  mutate(V = predict(V_reg, newdata = smoke0_df)) %>% 
  mutate(PB = predict(PB_reg, newdata = smoke0_df)) %>%
  mutate(FE = predict(FE_reg, newdata = smoke0_df)) %>%
  mutate(CU = predict(CU_reg, newdata = smoke0_df)) %>%
  mutate(CR = predict(CR_reg, newdata = smoke0_df)) %>%
  mutate(MN = predict(MN_reg, newdata = smoke0_df)) %>%
  mutate(ZN = predict(ZN_reg, newdata = smoke0_df)) %>%
  mutate(NI = predict(NI_reg, newdata = smoke0_df)) %>%
  mutate(S = predict(S_reg, newdata = smoke0_df)) %>%
  mutate(EC = predict(EC_reg, newdata = smoke0_df)) %>%
  mutate(OC = predict(OC_reg, newdata = smoke0_df)) %>% 
  # pivot longer
  pivot_longer(cols = AL:ammSO4, 
               names_to = 'species', 
               values_to = 'nonsmoke_conc')

# calculate total concentration in PM2.5 
predicted_concs <- predicted_smoke_conc_df %>% 
  left_join(predicted_nonsmoke_conc_df) %>% 
  filter(nonsmoke_conc > 0) %>% 
  filter(smoke_conc > 0) %>% 
  mutate(tot_pred_conc = nonsmoke_conc + smoke_conc)


rm(predicted_smoke_conc_df, predicted_nonsmoke_conc_df, nonsmoke0_df, smoke0_df) # drop to save memory

# calculate avg attributable fraction of species concentration to total concentration
avg_conc_site_month_yr <- predicted_concs %>%
  # get avg concentration across all sites for a given month and year
  group_by(species, year, month, site_id) %>% 
  dplyr::summarise(avg_smoke_conc = mean(smoke_conc, na.rm = TRUE),
                   avg_nonsmoke_conc = mean(nonsmoke_conc, na.rm = TRUE),
                   avg_tot_conc = mean(tot_pred_conc, na.rm = TRUE)) %>% 
  ungroup()


# now get the average concentration across all sites, for a given month and year
avg_conc_month_yr <- avg_conc_site_month_yr %>%
  # get avg concentration across all sites for a given month and year
  group_by(species, year, month) %>% 
  dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                   avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                   avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_f = factor(species, levels = c('CR', "V", "PB", 'NI', "CU", "MN", "ZN","FE", 'S', "EC", "OC"))) 

rm(predicted_concs) #drop for memory

# adjust dataframe so that it can plot stacked area plot
pred_spec_conc_long_df <- avg_conc_month_yr %>% 
  pivot_longer(cols = c(avg_smoke_conc, avg_nonsmoke_conc, avg_tot_conc), 
               names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
  mutate(pm_type_f = factor(pm_type, 
                            levels = c('avg_smoke_conc', 'avg_nonsmoke_conc', 'avg_tot_conc'))) %>% 
  filter(!is.na(spec_f)) %>% 
  filter(pm_type != 'avg_tot_conc')

datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-01-01"), by = "24 month")


# MAKE AN AREA PLOT OF AVG PM OVER TIME ------------------------------------------
stacked_species_area_plot <- ggplot(pred_spec_conc_long_df, 
                                    aes(x = mon_yr, y = avg_spec_conc, 
                                        fill = pm_type_f)) + 
  geom_area(alpha =.6) +
  # add in line for total predicted concentration
  # add in scale
  scale_x_date(labels = date_format("%Y"), breaks = datebreaks) +
  # add in fill manually
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c("firebrick3", "grey70", 'black'), 
                    labels = c('Smoke attributable concentration (ug/m3)',
                               'Nonmoke attributable concentration (ug/m3)',
                               'Total predicted species concentration in PM2.5'
                               )) +
  facet_wrap(~spec_f, scales = 'free', ncol = 3) +
  labs(x = 'Year',
       y = paste('Species Concentration (ug/m3)'),
       title = 'Contribution of wildfire smoke and nonsmoke PM2.5 to species concentration over time') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30", size = .5),
        plot.title = element_text(size=12, face='bold'),
        axis.text.x = element_text(face="plain", size = 10),
        axis.text.y = element_text(face="plain", size = 10),
        legend.position = c(0.85, 0.1)) 
stacked_species_area_plot

# save file
ggsave(
  filename = 'Fig4_smoke_attributable_fraction_trend_monthly.pdf',
  plot = stacked_species_area_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 12,
  height = 6,
  dpi = 320)   

# MAKE AN AREA PLOT OF AVG PM OVER TIME BY YEAR ------------------------------------------
avg_conc_yr <- avg_conc_site_month_yr %>%
  # get avg concentration across all sites for a given month and year
  group_by(species, year) %>% 
  dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                   avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                   avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(yr = as.Date(paste0(year, "-", "01", "-", "01"), format = "%Y-%m-%d")) %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_f = factor(species, levels = c('CR', "V", "PB", 'NI', "CU", "MN", "ZN","FE", 'S', "EC", "OC"))) 


# adjust dataframe so that it can plot stacked area plot
pred_spec_conc_long_df <- avg_conc_yr %>% 
  pivot_longer(cols = c(avg_smoke_conc, avg_nonsmoke_conc, avg_tot_conc), 
               names_to = 'pm_type', values_to = 'avg_spec_conc') %>% 
  mutate(pm_type_f = factor(pm_type, 
                            levels = c('avg_smoke_conc', 'avg_nonsmoke_conc', 'avg_tot_conc')))%>% 
  filter(pm_type != 'avg_tot_conc') %>% 
  filter(!is.na(spec_f))

datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-01-01"), by = "24 month")


# MAKE AN AREA PLOT OF AVG PM OVER TIME ------------------------------------------
YRstacked_species_area_plot <- ggplot(pred_spec_conc_long_df,
                                    aes(x = yr, y = avg_spec_conc, 
                                        fill = pm_type_f)) + 
  geom_area(alpha =.6) +
  # add in scale
  scale_x_date(labels = date_format("%Y"), breaks = datebreaks) +
  # add in fill manually
  scale_fill_manual(name = 'PM2.5 Type',
                    values=c("firebrick3", "grey70"), 
                    labels = c('Smoke attributable concentration (ug/m3)',
                               'Nonmoke attributable concentration (ug/m3)')) +
  facet_wrap(~spec_f, scales = 'free', ncol = 3) +
  labs(x = 'Year',
       y = paste('Species Concentration (ug/m3)'),
       title = 'Contribution of smoke and nonsmoke PM2.5 to species concentration over time') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30", size = .5),
        plot.title = element_text(size=12, face='bold'),
        axis.text.x = element_text(face="plain", size = 10),
        axis.text.y = element_text(face="plain", size = 10),
        legend.position = c(0.85, 0.1)) 
YRstacked_species_area_plot

# save file
ggsave(
  filename = 'Fig4_smoke_attributable_fraction_trend_YEAR.png',
  plot = YRstacked_species_area_plot,
  path = file.path(wip_gdrive_fp, 'figures/Fig4'),
  scale = 1,
  width = 12,
  height = 6,
  dpi = 320)   


return(pred_spec_conc_long_df)
}


