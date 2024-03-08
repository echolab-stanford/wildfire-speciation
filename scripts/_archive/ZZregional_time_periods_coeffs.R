# Emma Krasovich Southworth, emmars@stanford.edu
# Description: run a regression for each time period and region and then normalize
# by the baseline average for the region and time period to get the percent
# change from the baseline. then plot!
# Created: 9/6/23 | Last Updated: 9/12/23

# loadd(clean_PMspec_df, cache = drake::drake_cache("scripts/.drake"))

# function to create regional plot of coefficient by time period
create_regional_time_period_plot <- function(clean_PMspec_df) {

# select the vars that are needed for the
reg_df <- clean_PMspec_df %>% 
  dplyr::select(Dataset:smoke_day, monitor_month,
                totPM2.5, 
                smokePM,
                nonsmokePM,
                V, PB, FE, CU, CR, 
                MN, ZN, NI, EC, OC, S) %>% 
  mutate(time_period = case_when(
    year >= 2006 & year <= 2010 ~ '2006-2010',
    year > 2010 & year <= 2015 ~ '2010-2015',
    year > 2015 ~ '2015-2020',
  )) %>% # only keep rows that have all values so that its the same selection of sites for the analysis
na.omit()  

# calculate species regional averages for each chemical species:
# calculate the sample averages for each chemical species across full sample:
spec_sample_avgs_df <- tibble(
  species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
  baseline_avg = c(mean(reg_df$V, na.rm = TRUE), 
                   mean(reg_df$PB, na.rm = TRUE), 
                   mean(reg_df$FE, na.rm = TRUE), 
                   mean(reg_df$CU, na.rm = TRUE),  
                   mean(reg_df$CR, na.rm = TRUE),  
                   mean(reg_df$MN, na.rm = TRUE), 
                   mean(reg_df$ZN, na.rm = TRUE),  
                   mean(reg_df$NI, na.rm = TRUE),  
                   mean(reg_df$EC, na.rm = TRUE),  
                   mean(reg_df$OC, na.rm = TRUE),  
                   mean(reg_df$S, na.rm = TRUE))) 

spec_regional_avgs_df <- reg_df %>% 
  dplyr::select(Dataset:region, V:S) %>% 
  pivot_longer(cols =c(V:S), names_to = 'species', values_to = 'conc') %>% 
  group_by(region, species) %>% 
  dplyr::summarise(baseline_avg = mean(conc, na.rm = TRUE), .groups = 'drop') %>% 
  bind_rows(spec_sample_avgs_df %>% mutate(region ='Full sample'))

#  ------------------------------------------------------------------
# REGIONAL ANALYSIS ----------------------------------------
#  ------------------------------------------------------------------
#early period
early_regional <- reg_df %>% 
  filter(time_period == '2006-2010')

early_regional_avgs_df <- early_regional %>% 
  dplyr::select(Dataset:region, V:S) %>% 
  pivot_longer(cols =c(V:S), names_to = 'species', values_to = 'conc') %>% 
  group_by(region, species) %>% 
  dplyr::summarise(baseline_avg = mean(conc, na.rm = TRUE), .groups = 'drop') %>% 
  mutate(time_period = '2006-2010')


early_regional_PM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                        ~ smokePM + nonsmokePM | 
                          monitor_month + year, early_regional, fsplit = ~region)

# get coefficients and prepare for plotting
early_regionalPM_coeffs <- coeftable(early_regional_PM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         region = 'sample') %>% 
  dplyr::select(-id, -sample.var)  %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) %>% 
  mutate(time_period = '2006-2010')

# middle period
mid_regional <- reg_df %>% 
  filter(time_period == '2010-2015')

mid_regional_avgs_df <- mid_regional %>% 
  dplyr::select(Dataset:region, V:S) %>% 
  pivot_longer(cols =c(V:S), names_to = 'species', values_to = 'conc') %>% 
  group_by(region, species) %>% 
  dplyr::summarise(baseline_avg = mean(conc, na.rm = TRUE), .groups = 'drop') %>% 
  mutate(time_period = '2010-2015')

mid_regional_PM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                        ~ smokePM + nonsmokePM | 
                          monitor_month + year, mid_regional, fsplit = ~region)

mid_regionalPM_coeffs <- coeftable(mid_regional_PM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         region = 'sample') %>% 
  dplyr::select(-id, -sample.var)  %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) %>% 
  mutate(time_period = '2010-2015')

# late period
late_regional <- reg_df %>% 
  filter(time_period == '2015-2020')

late_regional_avgs_df <- late_regional %>% 
  dplyr::select(Dataset:region, V:S) %>% 
  pivot_longer(cols =c(V:S), names_to = 'species', values_to = 'conc') %>% 
  group_by(region, species) %>% 
  dplyr::summarise(baseline_avg = mean(conc, na.rm = TRUE), .groups = 'drop') %>% 
  mutate(time_period = '2015-2020')

late_regional_PM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                        ~ smokePM + nonsmokePM | 
                          monitor_month + year, late_regional, fsplit = ~region)
 
late_regionalPM_coeffs <- coeftable(late_regional_PM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         region = 'sample') %>% 
  dplyr::select(-id, -sample.var)  %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) %>% 
  mutate(time_period = '2015-2020')




# all periods for all regions baselines
all_baseline_avg_concs <- bind_rows(late_regional_avgs_df, 
                                    mid_regional_avgs_df) %>% 
  bind_rows(early_regional_avgs_df)


# combine all time periods and then normalize coefficents
all_time_periods_regs_coeffs_normalized <- bind_rows(late_regionalPM_coeffs, 
                                          mid_regionalPM_coeffs) %>% 
  bind_rows(early_regionalPM_coeffs) %>% 
  dplyr::select(time_period, region, species, spec_type, pm_type, Estimate, se) %>% 
  filter(region != 'Full sample') %>% 
  left_join(all_baseline_avg_concs, by = c('region', 'species', 'time_period')) %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_se = se/baseline_avg) 


# marginal plots
regional_time_reg_plot <- ggplot(all_time_periods_regs_coeffs_normalized %>% 
                              filter(pm_type == 'smokePM'), # %>% filter(region %in% c('pacific', 'Full sample', 'southeast')), 
                            aes(x = time_period,
                                y = 100*norm_est, 
                                color = region,
                                #shape = time_period
                                )) +
  #shape = spec_type)) +
  geom_point(size=3, alpha = 0.7, stat = "identity") + # shape = pm_type
  geom_linerange(aes(ymin = (100*norm_est - 100*norm_se), 
                     ymax = (100*norm_est + 100*norm_se), 
                     width =.2), stat = "identity") +
  #scale_x_discrete(limits = c("CR", "V", "PB", "NI", "CU", "MN", "ZN","FE", "S", "EC", "OC")) +
  scale_x_discrete(limits = c("2006-2010", "2010-2015", "2015-2020")) +
  #scale_y_log10() +
  scale_color_manual(values=c("#DC863B", "#FAD510",
                              "#649373", "#1B5656", 
                              "#5A283E", "#F2300F")) +
  scale_shape_manual(values=c(15, 16, 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  facet_grid(rows =vars(region), cols = vars(species)) +
  labs(y = '% change (ug/m3) ',
       x = 'Chemical Species',
       color = 'Region',
       shape = 'Time Period',
       title = "Chemical composition of smoke PM2.5 relative to regional and period baseline"
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey')
regional_time_reg_plot


return(all_time_periods_regs_coeffs_normalized)

}

# etable(regional_PM_reg) # look at regression results

