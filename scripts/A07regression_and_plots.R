# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: make a time series plot for each species

# loadd(clean_pm_spec_df, cache = drake::drake_cache("scripts/.drake"))

# set up list of species we care about
species_list <- c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S")

# select the vars that are needed for the
reg_df <- clean_pm_spec_df %>% 
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
  ))

# -----------------------------------------------------------------------------
# 1. GET BASELINE AVERAGES FOR EACH SPECIES (FULL SAMPLE + REGIONAL SAMPLE)
# -----------------------------------------------------------------------------
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

# calculate species regional averages for each chemical species:
spec_regional_avgs_df <- reg_df %>% 
  dplyr::select(Dataset:region, V:S) %>% 
  pivot_longer(cols =c(V:S), names_to = 'species', values_to = 'conc') %>% 
  group_by(region, species) %>% 
  dplyr::summarise(baseline_avg = mean(conc, na.rm = TRUE), .groups = 'drop') %>% 
  bind_rows(spec_sample_avgs_df %>% mutate(region ='Full sample'))
  
# -----------------------------------------------------------------------------
# 3. RUN REGRESSION FOR SPECIES
#   A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
#   B) RUN IN LEVELS ACROSS REGIONS, DIVIDE BETAS BY REGIONAL AVG TO GET % CHANGE RELATIVE TO REGIONAL BASELINE
# -----------------------------------------------------------------------------
full_sampPM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                   ~ smokePM + nonsmokePM | 
                     monitor_month + year, reg_df) 
# etable(full_sampPM_reg) # look at regression results

# get coefficients and prepare for plotting
full_sampPM_coeffs <- coeftable(full_sampPM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id)  %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", 
                                       "Transition metal",
                                       "Toxicity potentiator",
                                       "Secondary organic"
                                       ))) %>% 
  mutate(spec_f = factor(species, 
                            levels = c("V", "CR", "PB", 
                                       "CU", "MN", "NI", "ZN","FE", 
                                       "S", 
                                       "EC", "OC"
                                       ))) 


# merge sample avg for each species and divide each species' betas by full sample avg for each species
full_samp_PMcoeffs_normalized <- full_sampPM_coeffs %>% 
  left_join(spec_sample_avgs_df, by = 'species') %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_se = se/baseline_avg) 
  

# plot coefficients for speciation at the avg monitor which tells us how much 
# of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
# SMOKE ONLY BY SPECIES CATEGORY ------------------------------------------------------------------
tox_metal_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                               filter(pm_type == 'smokePM') %>% 
                               filter(spec_type == 'Toxic metal'), 
                        aes(x = species,
                            y = Estimate, 
                            color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                    ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c('CR', "V", "PB")) +
  scale_color_manual(values=c("firebrick3")) +
  labs(title = "Toxic Metal") +
  theme_light() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30"),
        title= element_text(size=10, face='plain', hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 

tox_potentiator_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                 filter(pm_type == 'smokePM') %>% 
                                 filter(spec_type == 'Toxicity potentiator'), 
                               aes(x = species,
                                   y = Estimate, 
                                   color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c('S')) +
  scale_color_manual(values=c("deepskyblue1")) +
  labs(#y = 'Species Concentration (ug/m3)',
       #x = 'Wildfire Smoke PM2.5 (ug/m3)',
       #color = 'Species Category',
       title = "Toxicity Potentiator") +
  theme_light() +
  theme(panel.border = element_blank(), 
       # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30"),
        title= element_text(size=10, face='plain', hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 

trans_metal_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                       filter(pm_type == 'smokePM') %>% 
                                       filter(spec_type == 'Transition metal'), 
                                     aes(x = species,
                                         y = Estimate, 
                                         color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c('NI', "CU", "MN", 'ZN', 'FE')) +
  scale_color_manual(values=c("orange1")) +
  labs(#y = 'Species Concentration (ug/m3)',
       #x = 'Wildfire Smoke PM2.5 (ug/m3)',
       #color = 'Species Category',
       title = "Transition Metal") +
  theme_light() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30"),
        title= element_text(size=10, face='plain', hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 

sec_org_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                   filter(pm_type == 'smokePM') %>% 
                                   filter(spec_type == 'Secondary organic'), 
                                 aes(x = species,
                                     y = Estimate, 
                                     color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c('EC', "OC")) +
  scale_color_manual(values=c('grey80')) +
  labs(#y = 'Species Concentration (ug/m3)',
       #x = 'Wildfire Smoke PM2.5 (ug/m3)',
       #color = 'Species Category',
       title = "Secondary Organic") +
  theme_light() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey30"),
        title= element_text(size=10, face='plain', hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 

# COMBINE PLOT
comb_plot<-ggarrange(tox_metal_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"),
          trans_metal_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"),
          tox_potentiator_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"),
          sec_org_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"))
   
# combine and annotate
comb_plot_ann<-annotate_figure(comb_plot, 
                top = text_grob("Marginal effect of wildfire smoke PM2.5 on species concentration", 
                                color = "black", face = "bold", size = 14),
                bottom = text_grob("Species", color = "black", size = 12),
                left = text_grob("Predicted species concentration (ug/m3)", size = 12, color = "black", rot = 90))

# save file
ggsave(
  filename = 'Fig2alt_smoke_coeffs_by_spectype.png',
  plot = comb_plot_ann,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320)   

# --------------------------------------------------------------------------------
# plot percent change for all species
# --------------------------------------------------------------------------------
pct_change_samp_reg_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                     filter(pm_type == 'smokePM'), 
                             aes(x = species,
                                 y = 100*norm_est, 
                                 color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = 100*(norm_est - norm_se), 
                    ymax = 100*(norm_est + norm_se)), stat = "identity") +
  scale_x_discrete(limits = c( "CR", "V", "PB", "NI", "CU", "MN","ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = '% Change (ug/m3) relative to sample baseline',
       x = 'Species',
       color = 'Species Category',
       title = "Percent Change of a given species concentation for 1 ug/m3 increase of smoke PM2.5") +
  theme_light() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
pct_change_samp_reg_plot

# save file
ggsave(
  filename = 'Fig2alt_pct_change_all_spec_smoke.png',
  plot = pct_change_samp_reg_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 

# -------------------------------------------------------------------------------
# SMOKE ONLY ALL SPECIES ---------------------------------------------------------
# -------------------------------------------------------------------------------
all_species_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                 filter(pm_type == 'smokePM'), 
                               aes(x = species,
                                   y = Estimate, 
                                   color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c( "CR", 'V', "PB", "NI", "CU", "MN", "ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  labs(y = 'Predicted Species Concentration (ug/m3)',
       x = 'Species',
       color = 'Species Category',
       title = "Marginal effect of wildfire smoke PM2.5 on chemical species concentation") +
  theme_light() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
all_species_smoke_plot

# save file
ggsave(
  filename = 'Fig2alt_all_spec_coeffs_smoke.png',
  plot = all_species_smoke_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 

# LOG SCALE
LOGall_species_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                   filter(pm_type == 'smokePM'), 
                                 aes(x = species,
                                     y = Estimate, 
                                     color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_y_log10() +
  scale_x_discrete(limits = c( "CR", 'V', "PB", "NI", "CU","MN", "ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  # facet_wrap(~spec_type, ncol = 2, scales = 'free_y') +
  labs(y = 'Predicted Species Concentration (ug/m3)',
       x = 'Species',
       color = 'Species Category',
       title = "Marginal effect of wildfire smoke PM2.5 on chemical species concentation") +
  theme_light() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
LOGall_species_smoke_plot

# save file
ggsave(
  filename = 'Fig2alt_all_spec_coeffs_smokeLOG.png',
  plot = LOGall_species_smoke_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 

# ------------------------------------------------------------------------------
# NONSMOKE VERSUS SMOKE
# ------------------------------------------------------------------------------
all_species_smoke_nonsmoke_plot <- ggplot(full_samp_PMcoeffs_normalized, 
                                 aes(x = pm_type,
                                     y = Estimate, 
                                     color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c("smokePM", "nonsmokePM")) +
  scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  facet_wrap(~spec_f, ncol = 3, scales = 'free_y') +
  labs(y = 'Predicted Species Concentration (ug/m3)',
       x = 'PM2.5 Type',
       color = 'Species Category',
       title = "Marginal effect of smoke PM2.5 and nonsmoke PM2.5 on chemical species concentation") +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
  #theme(strip.background = element_rect(color="black", fill="white", size=1))
all_species_smoke_nonsmoke_plot

# save file
ggsave(
  filename = 'Fig3_allPM_by_species.png',
  plot = all_species_smoke_nonsmoke_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 


all_species_smoke_nonsmoke_stacked_plot <- ggplot(full_samp_PMcoeffs_normalized, 
                                          aes(x = species,
                                              y = Estimate, 
                                              color=spec_type, shape = pm_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c( "CR", 'V', "PB", "NI", "CU","MN", "ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  scale_y_log10() +
  labs(y = 'Predicted Species Concentration (ug/m3)',
       x = 'Species',
       color = 'Species Category',
       title = "Marginal effect of smoke PM2.5 and nonsmoke PM2.5 on chemical species concentation") +
  theme_light() +
  theme(panel.border = element_blank(), 
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
all_species_smoke_nonsmoke_stacked_plot

# save file
ggsave(
  filename = 'Fig3_allPM_by_speciesSTACKEDLOG.png',
  plot = all_species_smoke_nonsmoke_stacked_plot,
  path = file.path(wip_gdrive_fp, 'figures/'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 






#  ------------------------------------------------------------------
# REGIONAL ANALYSIS ----------------------------------------
#  ------------------------------------------------------------------
regional_PM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                     ~ smokePM + nonsmokePM | 
                       monitor_month + year, reg_df, fsplit = ~region) 
etable(regional_PM_reg) # look at regression results

# get coefficients and prepare for plotting
regionalPM_coeffs <- coeftable(regional_PM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         region = 'sample') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
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
                                       "Toxicity potentiator"))) 


# merge sample avg for each species and divide each species' betas by full sample avg for each species
regionalPMcoeffs_normalized <- regionalPM_coeffs %>% 
  left_join(spec_regional_avgs_df, by = c('region', 'species')) %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_se = se/baseline_avg) 
  

# plot coefficients for speciation at the avg monitor which tells us how much 
# of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
regional_reg_plot <- ggplot(regionalPMcoeffs_normalized, 
                               aes(x = pm_type,
                                   y = 100*norm_est, 
                                   color=region)) +
  geom_point(size=2, alpha = 0.6, stat = "identity") + # shape = pm_type
  geom_errorbar(aes(ymin = (norm_est - norm_se)*100, 
                    ymax = (norm_est + norm_se)*100, 
                    width =.2), stat = "identity") +
  #scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  #scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  facet_wrap(~species, ncol = 6, scales = 'free_y') +
  labs(y = '% Change (ug/m3) relative to regional baseline',
       x = 'Chemical Species (ug/m3)',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Percent change in species concentration relative to regional baseline for a 1 ug/m3 increase in smoke and nonsmoke",
       subtitle = 'nonsmoke + smoke on RHS, monitor-month and year FE') +
  theme_light() 
regional_reg_plot

# marginal plots
regional_reg_plot <- ggplot(regionalPMcoeffs_normalized %>% filter(region %in% c('pacific', 'Full sample', 'southeast')), 
                            aes(x = pm_type,
                                y = Estimate, 
                                color=region)) +
  geom_point(size=2, alpha = 0.7, stat = "identity") + # shape = pm_type
  geom_linerange(aes(ymin = (Estimate - se), 
                    ymax = (Estimate + se), 
                    width =.2), stat = "identity") +
  #scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c('black','steelblue', 'plum', "red","goldenrod1", "deepskyblue", 'forestgreen')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  facet_wrap(~species, ncol = 4, scales = 'free_y') +
  labs(y = 'Concentration (ug/m3)',
       x = 'Chemical Species (ug/m3)',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Comparing chemical composition of nonsmoke and smoke PM by region",
       subtitle = 'nonsmoke + smoke on RHS, monitor-month and year FE') +
  theme_light() 
regional_reg_plot

# ---------------------------------------------------------------------------
# TEMPORAL HETEROGENEITY
# ---------------------------------------------------------------------------
temporal_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                        ~ smokePM + nonsmokePM | 
                          monitor_month + year, reg_df, fsplit = ~time_period) 
etable(temporal_reg) # look at regression results

# get coefficients and prepare for plotting
temporalPM_coeffs <- coeftable(temporal_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         time_per = 'sample') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
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
  mutate(spec_f = factor(species, 
                         levels = c("V", "CR", "PB", 
                                    "CU", "MN", "NI", "ZN","FE", 
                                    "S", 
                                    "EC", "OC"))) 

# marginal plot
temporal_reg_plot <- ggplot(temporalPM_coeffs, 
                            aes(x = time_per,
                                y = Estimate, 
                                color=spec_type, shape = pm_type)) +
  geom_point(size=2, alpha = 0.7, stat = "identity") + # shape = pm_type
  geom_linerange(aes(ymin = (Estimate - se), 
                     ymax = (Estimate + se)), stat = "identity") +
  scale_x_discrete(limits = c("Full sample", "2006-2010","2010-2015", "2015-2020")) +
  scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  facet_wrap(~spec_f, ncol = 3, scales = 'free_y') +
  labs(y = 'Concentration (ug/m3)',
       x = 'Time Period in Sample',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Comparing chemical composition of nonsmoke and smoke PM over time",
       subtitle = 'nonsmoke + smoke on RHS, monitor-month and year FE') +
  theme_light() 
temporal_reg_plot




# ---------------------------------------------------------------------------
# 3.CSN VS IMPROVE INTERACTION
# ---------------------------------------------------------------------------
monitor_regional_PM_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                        ~ smokePM + nonsmokePM + Dataset + smokePM*Dataset + nonsmokePM*Dataset | 
                          monitor_month + year, reg_df, fsplit = ~region) 
etable(monitor_regional_PM_reg) # look at regression results

# get coefficients and prepare for plotting
monitor_regionalPM_coeffs <- coeftable(monitor_regional_PM_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient',
         region = 'sample') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
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
                                       "Toxicity potentiator"))) 


# merge sample avg for each species and divide each species' betas by full sample avg for each species
regionalPMcoeffs_normalized <- regionalPM_coeffs %>% 
  left_join(spec_regional_avgs_df, by = c('region', 'species')) %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_se = se/baseline_avg) 


# plot coefficients for speciation at the avg monitor which tells us how much 
# of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
regional_reg_plot <- ggplot(regionalPMcoeffs_normalized, 
                            aes(x = pm_type,
                                y = 100*norm_est, 
                                color=region)) +
  geom_point(size=2, alpha = 0.6, stat = "identity") + # shape = pm_type
  geom_errorbar(aes(ymin = (norm_est - norm_se)*100, 
                    ymax = (norm_est + norm_se)*100, 
                    width =.2), stat = "identity") +
  #scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  #scale_color_manual(values=c("red","goldenrod1", "deepskyblue", 'grey')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  facet_wrap(~species, ncol = 6, scales = 'free_y') +
  labs(y = '% Change (ug/m3) relative to regional baseline',
       x = 'Chemical Species (ug/m3)',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Percent change in species concentration relative to regional baseline for a 1 ug/m3 increase in smoke and nonsmoke",
       subtitle = 'nonsmoke + smoke on RHS, monitor-month and year FE') +
  theme_light() 
regional_reg_plot







# csn_only <- reg_df %>% 
#   filter(Dataset == 'EPACSN')
# 
# csn_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN)
#                      ~ smokePM + nonsmokePM | 
#                        monitor_month + year, csn_only)
# 
# # IMPROVE SITES ONLY
# improve_only <- reg_df %>% 
#   filter(Dataset == 'IMPROVE')
# 
# improve_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN)
#                 ~ smokePM + nonsmokePM | 
#                   monitor_month + year, improve_only)
#                                
# # get coeffs for both improve and csn
# site_coeffs <- coeftable(improve_reg) %>%
#   mutate(site = 'IMPROVE') %>% 
#   bind_rows(coeftable(csn_reg) %>% 
#               mutate(site = 'EPACSN')) %>% 
#   rename(pval = 'Pr(>|t|)',
#          se = 'Std. Error',
#          species = 'lhs') %>% 
#   # get pvalues
#   mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
#   mutate(spec_type = case_when(
#     species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
#     species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
#     species=="EC" | species=="OC" ~ "Secondary organic",
#     species=="S" ~ "Toxicity potentiator",
#     TRUE ~ NA)) %>% 
#   mutate(spec_type = factor(spec_type, 
#                             levels = c("Toxic metal", "Transition metal",
#                                        "Secondary organic",
#                                        "Toxicity potentiator"))) 
#  
# # plot the comparison between CSN vs Improve
# site_reg_plot <- ggplot(site_coeffs, 
#                         aes(x = species,
#                             y = Estimate, 
#                             color=spec_type)) +
#   geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
#   scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
#   scale_y_continuous(trans = log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = scales::percent_format(accuracy = .01)) +
#   scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
#   scale_shape_manual(values=c(15, 17,18, 19)) +
#   geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
#   facet_wrap(~site)+
#   labs(x = 'Chemical Species (ug/m3)',
#        y = 'Estimate',
#        color = 'Site Type',
#        shape = 'PM Type',
#        title = "Comparing CSN vs Improve Speciation") +
#   theme_light() 
# site_reg_plot
# 
# # calculate the sample averages for each chemical species:
# csn_baseline_avgs_df <- tibble(
#   species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
#   baseline_avg = c(mean(csn_only$V, na.rm = TRUE), 
#                    mean(csn_only$PB, na.rm = TRUE), 
#                    mean(csn_only$FE, na.rm = TRUE), 
#                    mean(csn_only$CU, na.rm = TRUE),  
#                    mean(csn_only$CR, na.rm = TRUE),  
#                    mean(csn_only$MN, na.rm = TRUE), 
#                    mean(csn_only$ZN, na.rm = TRUE),  
#                    mean(csn_only$NI, na.rm = TRUE),  
#                    mean(csn_only$EC, na.rm = TRUE),  
#                    mean(csn_only$OC, na.rm = TRUE),  
#                    mean(csn_only$S, na.rm = TRUE)),
#   site ='EPACSN'
# ) %>% 
#   bind_rows(tibble(
#     species = c("V", "PB", "FE", "CU","CR", "MN", "ZN", "NI", "EC", "OC", "S"),
#     baseline_avg = c(mean(improve_only$V, na.rm = TRUE), 
#                      mean(improve_only$PB, na.rm = TRUE), 
#                      mean(improve_only$FE, na.rm = TRUE), 
#                      mean(improve_only$CU, na.rm = TRUE),  
#                      mean(improve_only$CR, na.rm = TRUE),  
#                      mean(improve_only$MN, na.rm = TRUE), 
#                      mean(improve_only$ZN, na.rm = TRUE),  
#                      mean(improve_only$NI, na.rm = TRUE),  
#                      mean(improve_only$EC, na.rm = TRUE),  
#                      mean(improve_only$OC, na.rm = TRUE),  
#                      mean(improve_only$S, na.rm = TRUE)),
#     site ='IMRPOVE'
#   ))
  
# merge with coeffs:
norm_site_coeffs <- site_coeffs %>% 
  left_join(csn_baseline_avgs_df, by = c('site', 'species')) %>% 
  mutate(norm_coeff = Estimate/baseline_avg)

# PLOT NORMALIZED CSN VS IMPROVE
# plot the comparison between CSN vs Improve
site_reg_plot <- ggplot(norm_site_coeffs, 
                        aes(x = species,
                            y = Estimate, 
                            color=spec_type)) +
  geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
  scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(accuracy = .01)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  facet_wrap(~site)+
  labs(x = 'Chemical Species (ug/m3)',
       y = 'Estimate',
       color = 'Site Type',
       shape = 'PM Type',
       title = "Comparing CSN vs Improve Speciation, normalized") +
  theme_light() 
site_reg_plot


# ---------------------------------------------------------------------------
# RUN A MODEL WHERE EACH PM TYPE IS PLOTTED ON SAME
# ---------------------------------------------------------------------------
smokePM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ smokePM | monitor_month + year, reg_df)
nonsmokePM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ nonsmokePM | monitor_month + year, reg_df)
totPM_reg = feols(c(CR, CU, EC, FE, MN, OC, PB, S, V, ZN) ~ totPM2.5 | monitor_month + year, reg_df)

# COMBINED COEFFS FOR PLOTTING
pm_type_coeffs <- coeftable(smokePM_reg) %>%
  mutate(pm_type = 'smoke PM2.5') %>% 
  bind_rows(coeftable(nonsmokePM_reg) %>%
              mutate(pm_type = 'nonsmoke PM2.5')) %>% 
  bind_rows(coeftable(totPM_reg) %>%
              mutate(pm_type = 'total PM2.5')) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  mutate(spec_type = case_when(
    species=="V" | species=="PB" | species=="CR" ~ "Toxic metal",
    species=="FE" | species=="CU" | species=="MN" | species=="ZN" | species=="NI" ~ "Transition metal",
    species=="EC" | species=="OC" ~ "Secondary organic",
    species=="S" ~ "Toxicity potentiator",
    TRUE ~ NA)) %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic metal", "Transition metal",
                                       "Secondary organic",
                                       "Toxicity potentiator"))) 

# join with species baseline avgs + divide to normalize betas
norm_pm_type_coeffs <- pm_type_coeffs %>% 
  left_join(spec_baseline_avgs_df, by = 'species') %>%
  mutate(norm_est = Estimate/baseline_avg)

# plot the normalized values for each smoke/nonsmoke/total pm
# plot the comparison between CSN vs Improve
pmtype_reg_plot <- ggplot(norm_pm_type_coeffs, 
                        aes(x = species,
                            y = Estimate, 
                            color=spec_type)) +
  geom_point(size=3, aes(shape = coefficient), alpha = 0.6, stat = "identity") +
  scale_x_discrete(limits = c("V", "CR", "PB", "CU", "MN", "NI", "ZN","FE", "S", "EC", "OC")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = scales::percent_format(accuracy = .01)) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  labs(x = 'Chemical Species (ug/m3)',
       y = 'Estimate',
       color = 'Species Category',
       shape = 'PM Type',
       title = "Comparing PM Types, normalized") +
  theme_light() 
pmtype_reg_plot


#  -----------------------------------------------------------------------------
# DO SOME REGRESSION CHECKS 
# ------------------------------------------------------------------------------
# CHECK 1) ARE BETAS THE SAME FOR WHEN WE RUN TOTAL PM ON LHS WITH AND W/O SMOKE
nonsmoke_control_reg = feols(conc_val_totPM2.5 ~ conc_val_smokePM + conc_val_nonsmokePM | monitor_month + year, conc_rel2avg_df) 
summary(nonsmoke_control_reg)
no_control_reg = feols(conc_val_totPM2.5 ~ conc_val_smokePM | monitor_month + year, conc_rel2avg_df) 
summary(no_control_reg)
nonosmoke_control_reg = feols(conc_val_totPM2.5 ~ conc_val_nonsmokePM | monitor_month + year, conc_rel2avg_df) 
summary(nonosmoke_control_reg)



return(reg_df)

}



