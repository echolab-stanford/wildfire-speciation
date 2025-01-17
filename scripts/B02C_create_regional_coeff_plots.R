# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: June 20, 2023
# Description: create a regional marginal coefficient plot

# loadd(c(clean_PMspec_df,parameter_categories, region_pal), cache = drake::drake_cache(".drake"))

# function
create_regional_coeff_plots <- function(clean_PMspec_df, parameter_categories, region_pal) {
  
# select the vars that are needed for the
reg_df <- clean_PMspec_df %>% 
  mutate(monitor_month = paste0(site_id,"_",month))  %>% 
  # drop Soil and zirconium, not needed for this analysis
  dplyr::select(-SOIL, -ZR)

# calculate species regional averages for each chemical species:
spec_regional_avgs_df <- reg_df %>% 
  dplyr::select(Dataset:region,monitor_month, AL:ZN, smoke_day) %>% 
  pivot_longer(cols =c(AL:ZN), names_to = 'species', values_to = 'conc') %>% 
  # filter to nonsmoke days to get relative to normal concentration 
  filter(smoke_day == 'nonsmoke day') %>% 
  group_by(region, species) %>% 
  dplyr::summarise(baseline_avg = mean(conc, na.rm = TRUE), .groups = 'drop') 

#  ------------------------------------------------------------------
# REGIONAL ANALYSIS ----------------------------------------
#  ------------------------------------------------------------------
regional_PM_reg = feols(c(AL, AS, BR, CA, CL,CR, CU, EC, FE, 
                          K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                          S,  SE, SI, SO4, SR, TI, V,  ZN)
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
  dplyr::select(-id, -sample.var) %>% 
  filter(region != 'Full sample')


# get coefficients and prepare for plotting
regionalPM_coeffs <- coeftable(regional_PM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         region = 'sample') %>% 
  #filter(pm_type == 'smokePM') %>% 
  filter(region != 'Full sample') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id, -sample.var)  %>% 
  left_join(parameter_categories, by = 'species') %>% 
  mutate(measure = 'MF') %>% 
  left_join(CIs,
            by = c('region','species', 'pm_type', 'measure'))


# merge sample avg for each species and divide each species' betas by full sample avg for each species
regionalPMcoeffs_normalized <- regionalPM_coeffs %>% 
  left_join(spec_regional_avgs_df, by = c('region', 'species')) %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_CI25 = CI25/baseline_avg,
         norm_CI975 = CI975/baseline_avg) 

limits <- c(" ",
  # "Organics"
  "Organic Carbon (OC)", "Elemental Carbon (EC)", " ", " ", 
  # "Halogens"
  "Bromine (Br)", "Chlorine (Cl)", " "," ", 
  #  "Nonmetals"
  "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)", "Nitrate (NO3)", "Selenium (Se)", " "," ",
  # "Other metals"
  "Aluminum (Al)", "Lead (Pb)", " "," ", 
  # "Metalloids"
  "Silicon (Si)", "Arsenic (As)", " "," ", 
  # "Transition metals"
  "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)", " "," ", 
  # "Alkali metals"
  "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)", " "," ", 
  # "Alkaline-earth metals"
  "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)", " ")


# PLOT ALL TOGETHER NO FACETS
regional_reg_plot <- ggplot(regionalPMcoeffs_normalized %>% 
                              filter(pm_type == 'smokePM'),
                            aes(x = species_long,
                                y = 100*norm_est,
                                color = region)) +
  geom_point(size=3, alpha = 0.6, stat = "identity", position = position_dodge(width = .8)) +
  geom_linerange(aes(ymin = (100*norm_CI25),
                     ymax = (100*norm_CI975)), stat = "identity", position = position_dodge(width = .8)) +
  scale_x_discrete(limits = rev(limits)) +
  scale_color_manual(values=region_pal) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = '% change relative to regional baseline on nonsmoke days',
       x = 'Chemical Species',
       color = 'Region',
  ) +
  theme_light() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        legend.position = "top",
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + # Rotate x-axis text by 90 degrees
  geom_vline(xintercept = 7.5, linetype = "dotted", color = 'grey85') +
  geom_vline(xintercept = 14.5, linetype = "dotted", color = 'grey85') + 
  geom_vline(xintercept = 29.5, linetype = "dotted", color = 'grey85') +
  geom_vline(xintercept = 34.5, linetype = "dotted", color = 'grey85') +
  geom_vline(xintercept = 41.5, linetype = "dotted", color = 'grey85') +
  geom_vline(xintercept = 52.5, linetype = "dotted", color = 'grey85') +
  geom_vline(xintercept = 59.5, linetype = "dotted", color = 'grey85') 
regional_reg_plot

# save file
ggsave(
  filename = 'SIFig5_regional_smoke_coeffs_NOFACET_wide_raw.pdf',
  plot = regional_reg_plot,
  path = file.path(results_fp, 'SI Figs'), 
  scale = 1,
  width = 11,
  height = 8,
  dpi = 320) 

ggsave(
  filename = 'SIFig5_regional_smoke_coeffs_NOFACET_wide_raw.png',
  plot = regional_reg_plot,
  path = file.path(results_fp, 'SI Figs'),
  scale = 1,
  width = 11,
  height = 8,
  dpi = 320) 

return(regionalPMcoeffs_normalized)
}

