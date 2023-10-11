# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Sept 12, 2023
# Description: create a map to go along with regional coefficients

# loadd(clean_PMspec_df, cache = drake::drake_cache("scripts/.drake"))

# set up df
reg_df <- clean_PMspec_df %>% 
  mutate(monitor_month = paste0(site_id,"_",month)) %>% 
  dplyr::select(Dataset:month, monitor_month, Date, site_id, MF_adj, smokePM, nonsmokePM_MF, 
                V, PB, FE, CU, CR, 
                MN, ZN, NI, EC, OC, S) #%>% 
 # filter(region == 'pacific')

# separate by high and low smoke
high_smoke_df <- reg_df %>% 
  filter(smokePM >= 25) 

low_smoke_df <- reg_df %>% 
  filter(smokePM < 25)

# -----------------------------------------------------------------------------
# run regressions for each 
lowsmoke_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN)
                        ~ smokePM + nonsmokePM_MF | 
                          monitor_month + year, low_smoke_df, cluster = 'site_id') 

# run regressions for each 
highsmoke_reg = feols(c(CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN) 
                      ~ smokePM + nonsmokePM_MF |
                        monitor_month + year, high_smoke_df, cluster = 'site_id') 


# get coefficients and prepare for plotting
coeffs <- coeftable(lowsmoke_reg) %>% 
  mutate(smoke_lvl = 'low smoke PM2.5') %>% 
  bind_rows(coeftable(highsmoke_reg) %>% 
              mutate(smoke_lvl = 'high smoke PM2.5')) %>% 
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
                         ))) %>% 
  filter(pm_type == 'smokePM')

## calculate 95 CI%
CIs <- confint(
  lowsmoke_reg) %>% 
  mutate(smoke_lvl = 'low smoke PM2.5') %>% 
  bind_rows(confint(highsmoke_reg) %>% 
              mutate(smoke_lvl = 'high smoke PM2.5')) %>% 
  rename(species = 'lhs',
         pm_type = 'coefficient',
         CI25 = '2.5 %',
         CI975 = '97.5 %') %>% 
  mutate(measure = 'MF') %>% 
  dplyr::select(-id)


low_high_smoke_df <- coeffs %>% 
  left_join(CIs)


PMlow2high_plot <- ggplot(low_high_smoke_df, aes(x = smoke_lvl,
                                      y = Estimate, 
                                      #shape = pm_type,
                                      color = spec_type)) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = CI25, 
                      ymax = CI975), stat = "identity") +
  #scale_x_discrete(limits = c('smokePM', 'nonsmokePM_MF')) +
  scale_x_discrete(limits = c('low smoke PM2.5', 'high smoke PM2.5')) +
  # scale_x_discrete(limits = c( "CR", "V", "PB", "NI", "CU", "MN","ZN","FE", "S", "EC", "OC")) +
  scale_color_manual(values=c("firebrick3","orange1", "deepskyblue1", 'grey80')) +
  #scale_shape_manual(values=c(16, 17)) +
  # scale_y_log10() +
  facet_wrap(~spec_f, scales= 'free_y', ncol = 3)+
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = 'Effect on concentration (ug/m3)',
       x = 'Smoke Level',
       color = 'Species Category',
       title = 'Species concentration on high smoke to low smoke days') + 
  theme_minimal() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') 
PMlow2high_plot
                                       


  