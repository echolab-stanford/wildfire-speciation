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

# long conc
# conc_long_df <- reg_df %>% 
#   dplyr::select(Dataset:region,monitor_month, AL:ZN, smoke_day) %>% 
#   pivot_longer(cols =c(AL:ZN), names_to = 'species', values_to = 'conc') %>% 
#   filter(!is.na(conc))

# Create a histogram of each regions number of stations
# hist_plot <- ggplot(conc_long_df, 
#                     aes(x = conc, fill = region)) +
#   geom_histogram(position = "identity", color = "black", alpha = 0.7, bins = 500) +
#   facet_wrap(~species, scales = "free") +
#   labs(title = "Distribution of concetration by species and region",
#        x = "Concentration (ug/m3)",
#        y = "Frequency") +
#   theme_minimal()
# hist_plot

#  ------------------------------------------------------------------
# REGIONAL ANALYSIS ----------------------------------------
#  ------------------------------------------------------------------
regional_PM_reg = feols(c(AL, AS, BR, CA, CL, CHL,CR, CU, EC, FE, 
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
  #filter(pm_type == 'smokePM') %>% 
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
  dplyr::select(-id, -sample.var, -`t value`)  %>% 
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

limits <- c(
  # "Organics"
  "Organic Carbon (OC)", "Elemental Carbon (EC)", " ", " ", 
  # "Halogens"
  "Bromine (Br)", "Chlorine (Cl)", "Chloride (Chl)", " "," ", 
  #  "Nonmetals"
  "Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Phosphorus (P)", "Selenium (Se)",  " "," ", 
  # "Other metals"
  "Titanium (Ti)", "Aluminum (Al)", "Lead (Pb)", " "," ", 
  # "Metalloids"
  "Silicon (Si)", "Arsenic (As)", " "," ", 
  # "Transition metals"
  "Zinc (Zn)", "Manganese (Mn)",  "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)", " "," ", 
  # "Alkali metals"
  "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)", " "," ", 
  # "Alkaline-earth metals"
  "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")

# PLOT ALL TOGETHER NO FACETS
regional_reg_plot <- ggplot(regionalPMcoeffs_normalized %>% 
                              filter(pm_type == 'smokePM'),
                            aes(x = species_long,
                                y = 100*norm_est,
                                color = region)) +
  geom_point(size=3, alpha = 0.7, stat = "identity", position = position_dodge(width = .8)) +
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
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis text by 90 degrees
regional_reg_plot

# save file
ggsave(
  filename = 'Fig3_regional_smoke_coeffs_NOFACET_wide_raw.pdf',
  plot = regional_reg_plot,
  path = file.path(results_fp, 'figures/Fig3'),
  scale = 1,
  width = 11,
  height = 8,
  dpi = 320) 

return(regionalPMcoeffs_normalized)
}

# 
# # PLOT BY REGION AS WELL -------------------------------------------------------
# region_list <- unique(regionalPMcoeffs_normalized$region)
# 
# # current_region <- region_list[4]
# 
# plot_each_region <- map_df(region_list, function(current_region) {
#   
#   current_regionalPMcoeffs <- regionalPMcoeffs_normalized %>% 
#     mutate(color_id = case_when(
#       region == 'midwest' ~ "#DC863B",
#       region == 'northeast' ~ "#FAD510", 
#       region == 'pacific' ~ "#649373", 
#       region == 'rocky_mountain' ~ "#1B5656", 
#       region == 'southeast' ~ "#5A283E", 
#       region == 'southwest' ~ "#F2300F"
#     )) %>% 
#     filter(region == current_region) 
# 
# # marginal plots
# one_region_reg_plot <- ggplot(current_regionalPMcoeffs, 
#                             aes(x = species,
#                                 y = 100*norm_est, 
#                                 color = region)) +
#   #shape = spec_type)) +
#   geom_point(size=3, alpha = 0.7, stat = "identity") + # shape = pm_type
#   geom_linerange(aes(ymin = (100*norm_CI25), 
#                      ymax = (100*norm_CI975)), stat = "identity") +
#   scale_x_discrete(limits = c("CR", "V", "PB", "NI", "CU", "MN", "ZN","FE", "S", "EC", "OC")) +
#   scale_color_manual(values=c(unique(current_regionalPMcoeffs$color_id))) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#   labs(y = '% Change relative to regional baseline',
#        x = 'Chemical Species',
#        #color = 'Region',
#        #shape = 'Species Category',
#        title = paste0(str_to_title(current_region))
#   ) +
#   theme_minimal() +
#   theme(panel.border = element_blank(), 
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         title= element_text(size=12, face='bold'),
#         axis.title.x = element_text(size=11, face = 'plain'),
#         axis.title.y = element_text(size=11, face = 'plain')) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#   guides(color = FALSE)
# one_region_reg_plot
# 
# # save file
# ggsave(
#   filename = paste('Fig3_', current_region, '_smoke_coeffs.png'),
#   plot = one_region_reg_plot,
#   path = file.path(wip_gdrive_fp, 'figures/Fig3'),
#   scale = 1,
#   width = 6,
#   height = 4,
#   dpi = 320) 
# 
# current_regionalPMcoeffs
# }) 



# axis_names <- regionalPMcoeffs_normalized %>% 
#   distinct(spec_type, species_long)
# # set up empty list to store plots
# plot_list = list()
# #Loop
# for(i in 1:length(unique(axis_names$spec_type))) {
# 
#   current_df <- regionalPMcoeffs_normalized %>% 
#     filter(spec_type == unique(axis_names$spec_type)[i])
#   
#   current_limits <- limits[limits %in% current_df$species_long]
#   
#   # marginal plots
#   regional_reg_plot <- ggplot(current_df, 
#                               aes(x = species_long,
#                                   y = 100*norm_est, 
#                                   color = region)) +
#     geom_point(size=2, alpha = 0.7, stat = "identity", position = position_dodge(width = .5)) +
#     geom_linerange(aes(ymin = (100*norm_CI25), 
#                        ymax = (100*norm_CI975)), stat = "identity", position = position_dodge(width = .5)) +
#     scale_x_discrete(limits = current_limits) +
#     scale_color_manual(values=c("#DC863B", "#FAD510",
#                                 "#649373", "#1B5656", 
#                                 "#5A283E", "#F2300F")) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#     labs(#y = '% Change relative to regional baseline on nonsmoke days',
#          #x = 'Chemical Species',
#          #color = 'Region',
#          title = paste0(unique(current_df$spec_type))
#     ) +
#     coord_flip() +
#     theme_light() +
#     theme(panel.border = element_blank(), 
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           title= element_text(size=12, face='bold'),
#           axis.title.x = element_text(size=12, face = 'plain'),
#           axis.title.y = element_text(size=12, face = 'plain')) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#     guides(color = 'none')
# 
#   plot_list[[i]] = regional_reg_plot
# 
# }
# #Plot
# all_regional_plots <- wrap_plots(plot_list, ncol = 1) + 
#   plot_annotation(title = "Regional differences in the chemical composition of wildfire smoke PM2.5") 
# all_regional_plots
# 
# 
# # save file
# ggsave(
#   filename = 'Fig3_regional_smoke_coeffs_facet.pdf',
#   plot = all_regional_plots,
#   path = file.path(wip_gdrive_fp, 'figures/Fig3'),
#   scale = 1,
#   width = 8,
#   height = 10,
#   dpi = 320) 


