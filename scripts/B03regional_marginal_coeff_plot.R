# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: June 20, 2023
# Description: create a regional marginal coefficient plot

# loadd(clean_PMspec_df, cache = drake::drake_cache("scripts/.drake"))

# function
create_regional_coeff_plots <- function(clean_PMspec_df) {
  
# set up list of species we care about
species_list <- c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")

# select the vars that are needed for the
reg_df <- clean_PMspec_df %>% 
  mutate(monitor_month = paste0(site_id,"_",month)) 

# calculate species regional averages for each chemical species:
spec_regional_avgs_df <- reg_df %>% 
  dplyr::select(Dataset:region,monitor_month, AL:ammSO4) %>% 
  pivot_longer(cols =c(AL:ammSO4), names_to = 'species', values_to = 'conc') %>% 
  group_by(region, species) %>% 
  dplyr::summarise(baseline_avg = mean(conc, na.rm = TRUE), .groups = 'drop') 

#  ------------------------------------------------------------------
# REGIONAL ANALYSIS ----------------------------------------
#  ------------------------------------------------------------------
regional_PM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                        ~ smokePM + nonsmokePM_MF | 
                          monitor_month + year, reg_df, fsplit = ~region, cluster = 'site_id') 
# etable(regional_PM_reg) # look at regression results

## calculate 95 CI%
CIs <- confint(
  regional_PM_reg) %>% 
  rename(species = 'lhs',
         region = 'sample',
         pm_type = 'coefficient',
         CI25 = '2.5 %',
         CI975 = '97.5 %') %>% 
  mutate(measure = 'MF') %>% 
  dplyr::select(-id)


# get coefficients and prepare for plotting
regionalPM_coeffs <- coeftable(regional_PM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         region = 'sample') %>% 
  filter(pm_type == 'smokePM') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id, -sample.var, -`t value`)  %>% 
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
  left_join(CIs, by = c('species', 'pm_type',  'region'))


# merge sample avg for each species and divide each species' betas by full sample avg for each species
regionalPMcoeffs_normalized <- regionalPM_coeffs %>% 
  left_join(spec_regional_avgs_df, by = c('region', 'species')) %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_CI25 = CI25/baseline_avg,
         norm_CI975 = CI975/baseline_avg) %>% 
  filter(region != 'Full sample') 


# marginal plots
regional_reg_plot <- ggplot(regionalPMcoeffs_normalized, 
                            aes(x = species,
                                y = 100*norm_est, 
                                color = region)) +
                                #shape = spec_type)) +
  geom_point(size=3, alpha = 0.7, stat = "identity") + # shape = pm_type
  geom_linerange(aes(ymin = (100*norm_CI25), 
                     ymax = (100*norm_CI975)), stat = "identity") +
  scale_x_discrete(limits = c("CR", "V", "PB", "NI", "CU", "MN", "ZN","FE", "S", "EC", "OC")) +
  #scale_y_log10() +
  scale_color_manual(values=c("#DC863B", "#FAD510",
                              "#649373", "#1B5656", 
                              "#5A283E", "#F2300F")) +
  #scale_shape_manual(values=c(0, 1, 2, 5, 8, 6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  #facet_wrap(~region, ncol = 6) +
  labs(y = '% Change relative to regional baseline',
       x = 'Chemical Species',
       color = 'Region',
       #shape = 'Species Category',
       title = "Regional differences in the chemical composition of wildfire smoke PM2.5"
       ) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey')
regional_reg_plot

# save file
ggsave(
  filename = 'Fig3_regional_smoke_coeffs.png',
  plot = regional_reg_plot,
  path = file.path(wip_gdrive_fp, 'figures/Fig3'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 



# PLOT BY REGION AS WELL -------------------------------------------------------
region_list <- unique(regionalPMcoeffs_normalized$region)

# current_region <- region_list[4]

plot_each_region <- map_df(region_list, function(current_region) {
  
  current_regionalPMcoeffs <- regionalPMcoeffs_normalized %>% 
    mutate(color_id = case_when(
      region == 'midwest' ~ "#DC863B",
      region == 'northeast' ~ "#FAD510", 
      region == 'pacific' ~ "#649373", 
      region == 'rocky_mountain' ~ "#1B5656", 
      region == 'southeast' ~ "#5A283E", 
      region == 'southwest' ~ "#F2300F"
    )) %>% 
    filter(region == current_region) 

# marginal plots
one_region_reg_plot <- ggplot(current_regionalPMcoeffs, 
                            aes(x = species,
                                y = 100*norm_est, 
                                color = region)) +
  #shape = spec_type)) +
  geom_point(size=3, alpha = 0.7, stat = "identity") + # shape = pm_type
  geom_linerange(aes(ymin = (100*norm_CI25), 
                     ymax = (100*norm_CI975)), stat = "identity") +
  scale_x_discrete(limits = c("CR", "V", "PB", "NI", "CU", "MN", "ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c(unique(current_regionalPMcoeffs$color_id))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = '% Change relative to regional baseline',
       x = 'Chemical Species',
       #color = 'Region',
       #shape = 'Species Category',
       title = paste0(str_to_title(current_region))
  ) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  guides(color = FALSE)
one_region_reg_plot

# save file
ggsave(
  filename = paste('Fig3_', current_region, '_smoke_coeffs.png'),
  plot = one_region_reg_plot,
  path = file.path(wip_gdrive_fp, 'figures/Fig3'),
  scale = 1,
  width = 6,
  height = 4,
  dpi = 320) 

current_regionalPMcoeffs
}) 


return(regionalPMcoeffs_normalized)
}
