# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: make a time series plot for each species

# loadd(clean_pm_spec_df, cache = drake::drake_cache("scripts/.drake"))

create_smoke_coeff_plots <- function(clean_pm_spec_df) {

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

  
# -----------------------------------------------------------------------------
# 3. RUN REGRESSION FOR SPECIES
#   A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
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
# --------------------------------------------------------------------------------
# plot percent change for all species
# --------------------------------------------------------------------------------
pct_change_samp_reg_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                     filter(pm_type == 'smokePM'), 
                             aes(x = species,
                                 y = 100*norm_est, 
                                 color=spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (100*norm_est - 100*norm_se), 
                    ymax = (100*norm_est + 100*norm_se)), stat = "identity") +
  scale_x_discrete(limits = c( "CR", "V", "PB", "NI", "CU", "MN","ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = expression(paste('% Change (ug/', m^3, ') relative to sample baseline')),
       x = 'Species',
       color = 'Species Category',
       title = expression(paste("% change in concentration for 1 ug/", m^3, " increase in smoke ", "PM"["2.5"]))) + 
  theme_light() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
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
  path = file.path(wip_gdrive_fp, 'figures/Fig2'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 

# -------------------------------------------------------------------------------
# SMOKE ONLY LOG SCALE
# -------------------------------------------------------------------------------
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
  labs(y = 
    expression(
      paste('Marginal effect of wildfire smoke', " PM"["2.5"])),
       x = 'Species',
       color = 'Species Category', 
       title = expression(paste("Marginal effect of wildfire smoke ", "PM"["2.5"], " on species concentrations (ug/", m^3, ")"))) + 
  theme_light() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
LOGall_species_smoke_plot

# save file

ggsave(
  filename = 'Fig2_all_spec_log_coeffs.png',
  plot = LOGall_species_smoke_plot,
  path = file.path(wip_gdrive_fp, 'figures/Fig2'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 

return(reg_df)

}

# SMOKE ONLY BY SPECIES CATEGORY ------------------------------------------------------------------
# tox_metal_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
#                                filter(pm_type == 'smokePM') %>% 
#                                filter(spec_type == 'Toxic metal'), 
#                         aes(x = species,
#                             y = Estimate, 
#                             color=spec_type)) +
#   geom_point(size=3, alpha = 0.6, stat = "identity") +
#   geom_linerange(aes(ymin = (Estimate - se), 
#                     ymax = (Estimate + se)), stat = "identity") +
#   scale_x_discrete(limits = c('CR', "V", "PB")) +
#   scale_color_manual(values=c("firebrick3")) +
#   labs(title = "Toxic Metal") +
#   theme_light() +
#   theme(panel.border = element_blank(), 
#         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "grey30"),
#         title= element_text(size=10, face='plain', hjust = 0.5)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
# 
# tox_potentiator_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
#                                  filter(pm_type == 'smokePM') %>% 
#                                  filter(spec_type == 'Toxicity potentiator'), 
#                                aes(x = species,
#                                    y = Estimate, 
#                                    color=spec_type)) +
#   geom_point(size=3, alpha = 0.6, stat = "identity") +
#   geom_linerange(aes(ymin = (Estimate - se), 
#                      ymax = (Estimate + se)), stat = "identity") +
#   scale_x_discrete(limits = c('S')) +
#   scale_color_manual(values=c("deepskyblue1")) +
#   labs(#y = 'Species Concentration (ug/m3)',
#        #x = 'Wildfire Smoke PM2.5 (ug/m3)',
#        #color = 'Species Category',
#        title = "Toxicity Potentiator") +
#   theme_light() +
#   theme(panel.border = element_blank(), 
#        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "grey30"),
#         title= element_text(size=10, face='plain', hjust = 0.5)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
# 
# trans_metal_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
#                                        filter(pm_type == 'smokePM') %>% 
#                                        filter(spec_type == 'Transition metal'), 
#                                      aes(x = species,
#                                          y = Estimate, 
#                                          color=spec_type)) +
#   geom_point(size=3, alpha = 0.6, stat = "identity") +
#   geom_linerange(aes(ymin = (Estimate - se), 
#                      ymax = (Estimate + se)), stat = "identity") +
#   scale_x_discrete(limits = c('NI', "CU", "MN", 'ZN', 'FE')) +
#   scale_color_manual(values=c("orange1")) +
#   labs(#y = 'Species Concentration (ug/m3)',
#        #x = 'Wildfire Smoke PM2.5 (ug/m3)',
#        #color = 'Species Category',
#        title = "Transition Metal") +
#   theme_light() +
#   theme(panel.border = element_blank(), 
#         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "grey30"),
#         title= element_text(size=10, face='plain', hjust = 0.5)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
# 
# sec_org_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
#                                    filter(pm_type == 'smokePM') %>% 
#                                    filter(spec_type == 'Secondary organic'), 
#                                  aes(x = species,
#                                      y = Estimate, 
#                                      color=spec_type)) +
#   geom_point(size=3, alpha = 0.6, stat = "identity") +
#   geom_linerange(aes(ymin = (Estimate - se), 
#                      ymax = (Estimate + se)), stat = "identity") +
#   scale_x_discrete(limits = c('EC', "OC")) +
#   scale_color_manual(values=c('grey80')) +
#   labs(#y = 'Species Concentration (ug/m3)',
#        #x = 'Wildfire Smoke PM2.5 (ug/m3)',
#        #color = 'Species Category',
#        title = "Secondary Organic") +
#   theme_light() +
#   theme(panel.border = element_blank(), 
#         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "grey30"),
#         title= element_text(size=10, face='plain', hjust = 0.5)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
# 
# # COMBINE PLOT
# library(ggpubr)
# comb_plot<-ggarrange(tox_metal_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"),
#           trans_metal_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"),
#           tox_potentiator_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"),
#           sec_org_smoke_plot + rremove("legend") + rremove("ylab") + rremove("xlab"))
#    
# # combine and annotate
# comb_plot_ann<-annotate_figure(comb_plot, 
#                 top = text_grob("Marginal effect of wildfire smoke PM2.5 on species concentration", 
#                                 color = "black", face = "bold", size = 14),
#                 bottom = text_grob("Species", color = "black", size = 12),
#                 left = text_grob("Predicted species concentration (ug/m3)", size = 12, color = "black", rot = 90))
# 
# # save file
# ggsave(
#   filename = 'Fig2alt_smoke_coeffs_by_spectype.png',
#   plot = comb_plot_ann,
#   path = file.path(wip_gdrive_fp, 'figures/'),
#   scale = 1,
#   width = 10,
#   height = 6,
#   dpi = 320) 

# -------------------------------------------------------------------------------
# SMOKE ONLY ALL SPECIES LEVELS ---------------------------------------------------------
# -------------------------------------------------------------------------------
# all_species_smoke_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
#                                  filter(pm_type == 'smokePM'), 
#                                aes(x = species,
#                                    y = Estimate, 
#                                    color=spec_type)) +
#   geom_point(size=3, alpha = 0.6, stat = "identity") +
#   geom_linerange(aes(ymin = (Estimate - se), 
#                      ymax = (Estimate + se)), stat = "identity") +
#   scale_x_discrete(limits = c( "CR", 'V', "PB", "NI", "CU", "MN", "ZN","FE", "S", "EC", "OC")) +
#   scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
#   scale_shape_manual(values=c(15, 17,18, 19)) +
#   labs(y = 'Predicted Species Concentration (ug/m3)',
#        x = 'Species',
#        color = 'Species Category',
#        title = "Marginal effect of wildfire smoke PM2.5 on chemical species concentation") +
#   theme_light() +
#   theme(panel.border = element_blank(), 
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         title= element_text(size=12, face='bold'),
#         axis.title.x = element_text(size=11, face = 'plain'),
#         axis.title.y = element_text(size=11, face = 'plain')) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
# all_species_smoke_plot
# 
# # save file
# ggsave(
#   filename = 'Fig2alt_all_spec_coeffs_smoke.png',
#   plot = all_species_smoke_plot,
#   path = file.path(wip_gdrive_fp, 'figures/'),
#   scale = 1,
#   width = 10,
#   height = 6,
#   dpi = 320) 
