# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: make a time series plot for each species

# loadd(clean_PMspec_df, cache = drake::drake_cache(".drake"))

create_smoke_coeff_plots <- function(clean_pm_spec_df) {
  
  # use this: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4144270/
  # set up parameters of interest
  params <- tibble(
    species = c("AL","AS","BR", "CA", "CL", "CHL","CR", "CU", "EC", "FE", 
                "K", "MG","MN", "NA", "NI", "NO3", "N2", "OC", "P",  "PB", "RB",
                "S",  "SE", "SI", "SO4", "SOIL", "SR", "TI", "V",  "ZN", "ZR")) %>% 
    mutate(spec_type = case_when(
      species == "AL" ~ 'toxic metal',
      species == "AS" ~ 'heavy metal',
      species == "BR" ~ 'toxic nonmetal',
      species == 'CA' ~ 'non-toxic metal',
      species == "CL" ~ 'toxic nonmetal',
      species == 'CHL' ~ 'toxic nonmetal',
      species == "CR" ~ 'heavy metal',
      species == "CU" ~ 'heavy metal',
      species == "EC" ~ 'secondary organic',
      species == "FE" ~ 'toxic metal',
      species == "K" ~ 'toxic metal',
      species == "MG" ~ 'non-toxic metal',
      species == "MN" ~ 'toxic metal',
      species == "NA" ~ 'toxic metal',
      species == "NI" ~ 'heavy metal',
      species == "NO3" ~ 'secondary inorganic',
      species == "N2" ~ 'non-toxic nonmetal',
      species == "OC" ~ 'secondary organic',
      species == "P" ~ 'toxic nonmetal',
      species == "PB" ~ 'heavy metal',
      species == 'RB' ~ 'non-toxic metal',
      species == "S" ~ 'toxicity potentiator',
      species == "SE" ~ 'toxic nonmetal',
      species == "SI" ~ 'non-toxic metal',
      species == 'SO4' ~ 'secondary inorganic',
      species == "SOIL" ~ 'non-toxic nonmetal',
      species == 'SR' ~ 'non-toxic metal',
      species == "TI" ~ 'heavy metal',
      species == "V" ~ 'heavy metal',
      species == "ZN" ~ 'heavy metal',
      species == "ZR" ~ 'toxic metal'
    )) %>% 
    mutate(spec_type = str_to_sentence(spec_type))


# select the vars that are needed for the
reg_df <- clean_PMspec_df %>% 
  mutate(monitor_month = paste0(site_id,"_",month))  %>% 
  dplyr::select(Dataset, state_name, region, year, month, Date, monitor_month, 
                site_id, MF_adj, smokePM, nonsmokePM_MF, AL:ZR) 

  
# -----------------------------------------------------------------------------
# 1. GET BASELINE AVERAGES FOR EACH SPECIES (FULL SAMPLE + REGIONAL SAMPLE)
# -----------------------------------------------------------------------------
# calculate the sample averages for each chemical species across full sample 
# on nonsmoke days - rather than all days, bc some places may be really smoky 
# relative to others and have more of
spec_sample_avgs_df <- reg_df %>% 
  dplyr::select(site_id, Date, AL:ZR, smoke_day) %>% 
  pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val')




  
# -----------------------------------------------------------------------------
# 3. RUN REGRESSION FOR SPECIES
#   A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
# -----------------------------------------------------------------------------
# the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
full_sampPM_regMF = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                            K, MG,MN, `NA`, NI, NO3, N2, OC, P,  PB, RB,
                            S,  SE, SI, SO4, SOIL, SR, TI, V,  ZN, ZR)
                   ~ smokePM + nonsmokePM_MF | 
                     monitor_month + year, reg_df, cluster = 'site_id') 

# etable(full_sampPM_regMF) # look at regression results

## calculate 95 CI%
CIs <- confint(
    full_sampPM_regMF) %>% 
      rename(species = 'lhs',
             pm_type = 'coefficient',
             CI25 = '2.5 %',
             CI975 = '97.5 %') %>% 
      mutate(measure = 'MF') %>% 
  dplyr::select(-id)

# get coefficients and prepare for plotting
full_sampPM_coeffsMF <- coeftable(full_sampPM_regMF) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient') %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id)  %>% 
  left_join(params, by = 'species') %>% 
  mutate(spec_type = factor(spec_type, 
                            levels = c("Toxic Metal", 
                                       "Heavy Metal",
                                       "Toxicity Potentiator",
                                       "Secondary Organic",
                                       "Toxic Nonmetal",
                                       "Toxic Semi-Metal"
                                       ))) %>% 
  mutate(measure = 'MF') %>% 
  left_join(CIs, by = c('species', 'pm_type', 'measure'))


# merge sample avg for each species and divide each species' betas by full sample avg for each species
full_samp_PMcoeffs_normalized <- full_sampPM_coeffsMF %>% 
  left_join(spec_sample_avgs_df, by = 'species') %>% 
  mutate(norm_est = Estimate/baseline_avg, # how much a species has changed relative to its baseline
         norm_CI25 = CI25/baseline_avg,
         norm_CI975 = CI975/baseline_avg) %>% 
  filter(pm_type == 'smokePM')
  

# plot coefficients for speciation at the avg monitor which tells us how much 
# of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
# --------------------------------------------------------------------------------
# plot percent change for all species
# --------------------------------------------------------------------------------
spec_type_pal <- wes_palette('FantasticFox1', 5, type = c("discrete"))[2:5]
spec_type_pal2 <- wes_palette('Rushmore1', 5, type = c("discrete"))[2:4]
spec_type_pal3 <- wes_palette('Royal1', 4, type = c("discrete"))[1]

colors <- c(spec_type_pal, spec_type_pal2, spec_type_pal3)
# "Toxic metal"          "Heavy metal"    "Non-toxic metal"  "Toxic nonmetal"  "Non-toxic nonmetal"
#  "Toxicity potentiator "Secondary organic"    "Secondary inorganic"    

# c("#B40F20", "#E58601", "#E2D200" "#46ACC8", "#0B775E", "#35274A","#CAB38C" "#899DA4")


pct_change_samp_reg_plot <- ggplot(full_samp_PMcoeffs_normalized %>% 
                                     filter(measure == 'MF'), 
                             aes(x = species,
                                 y = 100*norm_est, 
                                 color=spec_type, 
                                 #shape = measure
                                 )) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (100*norm_CI25), 
                    ymax = (100*norm_CI975)), stat = "identity") +
  scale_x_discrete(limits = c("SE", "CR", "PB","AS","S", "CL", "P", 
                               "NI", "V", "CU", "AL", "TI", "FE", 
                              "MN","ZN", "K",  "EC", "OC")) +
  scale_color_manual(values=c("#B40F20", "#E2D200", "#46ACC8", "#E58601", "#0B775E", "#35274A", "#899DA4"))+
  #scale_shape_manual(values = c(16,17)) +
 # facet_wrap(~spec_type)+
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = expression(paste('% Change relative to sample baseline')),
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
  filename = 'Fig2_pct_change_PM_spec_smokeMF_more_chems.png',
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
  geom_linerange(aes(ymin = CI25, 
                     ymax = CI975), stat = "identity") +
  scale_y_log10() +
  scale_x_discrete(limits = c("SE", "CR", "PB","AS","S", "CL", "P", 
                               "NI", "V", "CU", "AL", "TI", "FE", 
                              "MN","ZN", "K",  "EC", "OC")) +
  scale_color_manual(values=c("#B40F20", "#E2D200", "#46ACC8", "#E58601", "#0B775E", "#35274A", "#899DA4"))+
  # scale_x_discrete(limits = c( "CR", 'V', "PB", "NI", "CU", "FE", "MN", "ZN", "S", "EC", "OC")) +
  # scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = 
    expression(
      paste('Effect')),
       x = 'Species',
       color = 'Species Category', 
       title = expression(paste("Variation in how wildfire smoke ", "PM"["2.5"], " effects species' concentrations"))) + 
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
  filename = 'Fig2_PM_MF_spec_log_coeffs_more_chems.png',
  plot = LOGall_species_smoke_plot,
  path = file.path(wip_gdrive_fp, 'figures/Fig2'),
  scale = 1,
  width = 10,
  height = 6,
  dpi = 320) 

return(reg_df)

}


 
# 







# spec_sample_avgs_df <- tibble(
#   species = c("AL", "AS", "BR", "CA", "CHL","CL", "CR", "CU", "EC", "FE", "K",  
#               "MG", "MN", "NA", "NI", "NO3" ,"N2", "OC", "P", "PB", "RB", "S",
#               "SE", "SI", "SO4", "SOIL","SR", "TI", "V",  "ZN", "ZR"),
#   baseline_avg = c(mean(reg_df$AL, na.rm = TRUE), 
#                    mean(reg_df$AS, na.rm = TRUE), 
#                    mean(reg_df$BR, na.rm = TRUE), 
#                    mean(reg_df$CL, na.rm = TRUE),
#                    mean(reg_df$CR, na.rm = TRUE), 
#                    mean(reg_df$CU, na.rm = TRUE), 
#                    mean(reg_df$EC, na.rm = TRUE),  
#                    mean(reg_df$FE, na.rm = TRUE), 
#                    mean(reg_df$K, na.rm = TRUE),
#                    mean(reg_df$MN, na.rm = TRUE), 
#                    mean(reg_df$`NA`, na.rm = TRUE),  
#                    mean(reg_df$NI, na.rm = TRUE),
#                    mean(reg_df$OC, na.rm = TRUE), 
#                    mean(reg_df$P, na.rm = TRUE),  
#                    mean(reg_df$PB, na.rm = TRUE), 
#                    mean(reg_df$S, na.rm = TRUE),
#                    mean(reg_df$SE, na.rm = TRUE),  
#                    mean(reg_df$TI, na.rm = TRUE), 
#                    mean(reg_df$V, na.rm = TRUE), 
#                    mean(reg_df$ZN, na.rm = TRUE),  
#                    mean(reg_df$ZR, na.rm = TRUE))
# )
