
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





