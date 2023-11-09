# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 17, 2023
# Description: make a time series plot for each species

# loadd(c(parameter_categories, spec_pal, clean_PMspec_df), cache = drake::drake_cache(".drake"))

create_smoke_coeff_plots <- function(clean_PMspec_df, parameter_categories, spec_pal) {
  
# select the vars that are needed for the regression
reg_df <- clean_PMspec_df %>% 
  dplyr::select(Dataset, state_name, region, year, month, Date, 
                monitor_month, smoke_day, site_id, MF_adj, smokePM, 
                nonsmokePM_MF, AL:ZR) 

  
# -----------------------------------------------------------------------------
# 1. GET BASELINE AVERAGES FOR EACH SPECIES (FULL SAMPLE + REGIONAL SAMPLE)
# -----------------------------------------------------------------------------
# calculate the baseline for each chemical species across full sample by filtering to nonsmoke days
# rather than all days, bc some places may be really smoky 
# relative to others 
spec_ns_samp_avgs_df <- reg_df %>% 
  dplyr::select(site_id, Date, AL:ZR, smoke_day) %>% 
  pivot_longer(cols = c(AL:ZR), names_to ='species', values_to = 'conc_val') %>% 
  filter(smoke_day == 'nonsmoke day') %>% 
  group_by(species) %>% 
  dplyr::summarise(avg_nonsmoke_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')


# -----------------------------------------------------------------------------
# 3. RUN REGRESSION FOR SPECIES
#   A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
# -----------------------------------------------------------------------------
# the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
full_sampPM_regMF = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
                            K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
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
  mutate(pval = round(pval, digits = 3)) %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id)  %>% 
  left_join(parameter_categories, by = 'species') %>% 
  mutate(spec_type = fct_relevel(spec_type,
                                 c("Heavy metal", "Toxic metal", "Toxic nonmetal", 
                                   "Toxicity potentiator", "Non-toxic metal", 
                                   "Non-toxic nonmetal", "Secondary inorganic", "Secondary organic"
                                 ))) %>% 
  mutate(measure = 'MF') %>% 
  left_join(CIs, by = c('species', 'pm_type', 'measure'))


# merge sample avg for each species and divide each species' betas by full sample avg for each species
full_samp_PMcoeffs_normalized <- full_sampPM_coeffsMF %>% 
  left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
  mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
         norm_CI25 = CI25/avg_nonsmoke_spec_conc,
         norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
  filter(pm_type == 'smokePM')
  

# plot coefficients for speciation at the avg monitor which tells us how much 
# of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
# --------------------------------------------------------------------------------
# plot percent change for all species
# --------------------------------------------------------------------------------
    pct_change_samp_reg_plot <- ggplot(full_samp_PMcoeffs_normalized, 
                             aes(x = species_long,
                                 y = 100*norm_est, 
                                 color=spec_type, 
                                 )) +
  geom_point(size=4, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (100*norm_CI25), 
                    ymax = (100*norm_CI975)), stat = "identity") +
  scale_x_discrete(limits = c("Organic Carbon (OC)", "Elemental carbon (EC)", # secondary organics
                              #secondary inorganic
                              "Sulfate (SO4)", "Nitrate (NO3)", 
                              # nontoxic metals
                              "Strontium (Sr)",  "Calcium (Ca)", "Silicon (Si)", "Magnesium (Mg)", "Rubidium (Rb)",
                              # nontoxic nonmetals
                              "Soil", "Nitrogen (N2)",
                              # toxicity potentiator
                              "Sulfur (S)",
                              # toxic nonmetals
                              "Phosphorus (P)", "Bromine (Br)", "Selenium (Se)", "Chlorine (Cl)", "Chloride (Chl)",
                              # toxic metals
                              "Potassium (K)",  "Zinc (Zn)", "Manganese (Mn)", "Titanium (Ti)",
                              "Alumnium (Al)", "Iron (Fe)", "Copper (Cu)",
                              "Vanadium (V)", "Nickel (Ni)", "Sodium (Na)", "Zirconium (Zr)",  
                              # heavy metals
                              "Arsenic (As)", "Lead (Pb)", "Chromium (Cr)")) +
  scale_color_manual(values= spec_pal) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = expression(paste('% Change relative to average nonsmoke day concentration')),
       x = 'Species',
       #color = 'Species Category',
       title = expression(paste("% change in concentration for 1 ug/", m^3, " increase in smoke ", "PM"["2.5"]))) + 
  theme_light() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  guides(color = 'none')
pct_change_samp_reg_plot


# save file
ggsave(
  filename = 'Fig2_pct_change_PM_spec_smokeMF_all_chems.pdf',
  plot = pct_change_samp_reg_plot,
  path = file.path(wip_gdrive_fp, 'figures/Fig2'),
  scale = 1,
  width = 7,
  height = 10,
  dpi = 320) 

# -------------------------------------------------------------------------------
# SMOKE ONLY LOG SCALE
# -------------------------------------------------------------------------------
LOGall_species_smoke_plot <- ggplot(full_sampPM_coeffsMF %>% 
                                   filter(pm_type == 'smokePM'), 
                                 aes(x = species_long,
                                     y = Estimate, 
                                     color= spec_type)) +
  geom_point(size=4, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = CI25, 
                     ymax = CI975), stat = "identity") +
 # scale_y_log10() +
  scale_x_discrete(limits = c("Organic Carbon (OC)", "Elemental carbon (EC)", # secondary organics
                              #secondary inorganic
                              "Sulfate (SO4)", "Nitrate (NO3)", 
                              # nontoxic metals
                              "Strontium (Sr)",  "Calcium (Ca)", "Silicon (Si)", "Magnesium (Mg)", "Rubidium (Rb)",
                              # nontoxic nonmetals
                              "Soil", "Nitrogen (N2)",
                              # toxicity potentiator
                              "Sulfur (S)",
                              # toxic nonmetals
                              "Phosphorus (P)", "Bromine (Br)", "Selenium (Se)", "Chlorine (Cl)", "Chloride (Chl)",
                              # toxic metals
                              "Potassium (K)",  "Zinc (Zn)", "Manganese (Mn)", "Titanium (Ti)",
                              "Alumnium (Al)", "Iron (Fe)", "Copper (Cu)",
                              "Vanadium (V)", "Nickel (Ni)", "Sodium (Na)", "Zirconium (Zr)",  
                              # heavy metals
                              "Arsenic (As)", "Lead (Pb)", "Chromium (Cr)")) +
  scale_color_manual(values = spec_pal)+
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') + # 0.0000045184477
  scale_y_continuous(trans='log10', limits = c(0.000001,.5) , breaks=c(0.000001, 0.00005, 0.00001, .0005, .001, .005,
                                                                       0.01, 0.05, 0.1, 0.5),
                     labels=c('0.000001', '0.00005', '0.00001', '.0005', '.001', '.005',
                              '0.01', '0.05', '0.1', '0.5')) +
  labs(y = 'Effect on species concentration (ug/m3)',
       x = 'Species',
       color = 'Species Category', 
       title = expression(paste("Variation in how wildfire smoke ", "PM"["2.5"], " effects species' concentrations"))) + 
  theme_light() +
 coord_flip() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) 
  
LOGall_species_smoke_plot

# save file

ggsave(
  filename = 'Fig2_PM_MF_spec_log_coeffs.pdf',
  plot = LOGall_species_smoke_plot,
  path = file.path(wip_gdrive_fp, 'figures/Fig2'),
  scale = 1,
  width = 8,
  height = 10,
  dpi = 320) 

return(reg_df)

}

