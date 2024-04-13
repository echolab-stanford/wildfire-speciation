# # Emma Krasovich Southworth, emmars@stanford.edu
# # Last Updated: April 8, 2024
# # Description: sensitivity around changing negative values
# 
# # loadd(c(spec_w_smoke_pm_df,spec_pal), cache = drake_cache)
# 
# 
# # step 1: clean up the speciation data
# cleaned_spec_df <- spec_w_smoke_pm_df %>% 
#   # change inf to NA + NaNs to NAs
#   mutate_if(is.numeric, ~ifelse(. == Inf, NA, .)) %>% 
#   mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .)) %>% 
#   mutate_if(is.numeric, ~ifelse(. < -100, NA, .)) %>%
#   # mutate_if(is.numeric, ~ifelse(. == -999.00000, NA, .)) %>%
#   # mutate_if(is.numeric, ~ifelse(. == -999.00000, NA, .)) %>%
#   # drop a row if concentrations for all chemicals are NA
#   filter_at(vars(MF, RCFM, smokePM2.5, AL, AS, BR, CA, CHL, CL, CR, CU,
#                  EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
#                  RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
#             any_vars(!is.na(.))) %>%
#   # change any negative value in the data to NA
#   # mutate_at(vars(MF, RCFM, smokePM2.5, AL, AS, BR, CA, CHL, CL, CR, CU,
#   #                EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
#   #                RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
#   #           ~ifelse(. < 0, NA, .)) %>%
#   # change any 0 value in the data to NA, except for smoke data
#   mutate_at(vars(AL, AS, BR, CA, CHL, CL, CR, CU,
#                  EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
#                  RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
#             ~ifelse(. == 0, NA, .)) %>%
#   mutate(year = year(Date),
#          month = month(Date),
#          doy = yday(Date)) %>% 
#   # drop anything prior to 2005 because we are limited by PM2.5 predictions (childs et al, 2022)
#   filter(year > 2005 & year < 2021) %>%
#   group_by(site_id) %>%
#   mutate(st_date = year(min(Date)),
#          end_date = year(max(Date))) %>%
#   ungroup() %>%
#   # get duration that a site is online
#   mutate(duration = (end_date - st_date) + 1) %>%
#   # create seasons 
#   mutate(season = case_when(
#     month %in% c(12,1,2) ~ 'winter',
#     month %in% c(3,4,5) ~ 'spring',
#     month %in% c(6,7,8) ~ 'summer',
#     month %in% c(9,10,11) ~ 'autumn',
#   )) %>% 
#   # create regions
#   mutate(region = case_when(
#     state_name %in% c('Oregon', 'California', 'Washington') ~ 'pacific',
#     state_name %in% c('Nevada', 'Utah', 'Colorado', 'Wyoming', 'Montana', 'Idaho') ~ 'rocky_mountain',
#     state_name %in% c('Arizona', 'New Mexico', 'Texas', 'Oklahoma') ~ 'southwest',
#     state_name %in% c('North Dakota', 'South Dakota', 'Nebraska', 'Kansas', 'Minnesota', 'Iowa', 
#                       'Missouri', 'Wisconsin', 'Illinois', 'Indiana', 'Michigan', 'Ohio') ~ 'midwest',
#     state_name %in% c('Kentucky', 'West Virginia', 'Virginia', 'Tennessee', 'North Carolina', 'Mississippi', 
#                       'Alabama', 'Georgia', 'South Carolina', 'Florida', 'Louisiana', 'Arkansas') ~ 'southeast',
#     state_name %in% c('Maine', 'Vermont', 'New Hampshire', 'Massachusetts', 'Rhode Island', 'Connecticut', 
#                       'New York', 'New Jersey', 'Pennsylvania', 'Delaware', 'Maryland', 'District of Columbia') ~ 'northeast',
#   )) %>% 
#   # add in units column for each species, all species are same units
#   mutate(units = 'ug_m3') %>%
#   # ensure unique observations
#   distinct() %>% 
#   # do some cleaning of PM2.5 data, when measure for total PM is missing, set = to smoke PM
#   mutate(MF_adj = ifelse(is.na(MF), smokePM2.5, MF),
#          RCFM_adj = ifelse(is.na(RCFM), smokePM2.5, RCFM)) %>% 
#   # now drop if all speciation vars are NA for a row
#   # drop a row if concentrations for all chemicals are NA
#   filter_at(vars(AL, AS, BR, CA, CHL, CL, CR, CU,
#                  EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
#                  RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR),
#             any_vars(!is.na(.))) %>% 
#   # calculate station smoke and non smoke using total PM var (PM2.5)
#   mutate(nonsmokePM_MF = MF_adj - smokePM2.5,
#          nonsmokePM_RCFM = RCFM_adj - smokePM2.5) %>% 
#   # add an indicator for smoke day or nonsmoke day
#   mutate(smoke_day = case_when(
#     smokePM2.5 == 0 ~ 'nonsmoke day',
#     smokePM2.5 != 0 ~ "smoke day")) %>% 
#   # add a monitor month for FEs
#   mutate(monitor_month = paste0(site_id,"_",month)) %>% 
#   dplyr::select(Dataset, state_name, region, season, year, month, doy, 
#                 duration, st_date, end_date, monitor_month, Date, site_id, smoke_day, 
#                 MF_adj, RCFM_adj, smokePM = 'smokePM2.5', nonsmokePM_MF, nonsmokePM_RCFM,
#                 AL, AS, BR, CA, CHL, CL, CR, CU,
#                 EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
#                 RB, S, SE, SI, SO4, SOIL, SR, TI, V, ZN, ZR, units, 
#                 long, lat, epsg, MF, RCFM)
# 
# 
# # select the vars that are needed for the regression
# reg_df <- cleaned_spec_df %>% 
#   dplyr::select(Dataset, state_name, region, year, month, Date, 
#                 monitor_month, smoke_day, site_id, MF_adj, smokePM, 
#                 nonsmokePM_MF, AL:ZR) %>% 
#   # drop Soil and zirconium, not needed for this analysis
#   dplyr::select(-SOIL, -ZR)
# 
# # -----------------------------------------------------------------------------
# # 1. GET BASELINE AVERAGES FOR EACH SPECIES (FULL SAMPLE + REGIONAL SAMPLE)
# # -----------------------------------------------------------------------------
# # calculate the baseline for each chemical species across full sample by filtering to nonsmoke days
# # rather than all days, bc some places may be really smoky 
# # relative to others 
# spec_ns_samp_avgs_df <- reg_df %>% 
#   dplyr::select(site_id, Date, AL:ZN, smoke_day) %>% 
#   pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>% 
#   filter(smoke_day == 'nonsmoke day') %>% 
#   group_by(species) %>% 
#   dplyr::summarise(avg_nonsmoke_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')
# 
# 
# # -----------------------------------------------------------------------------
# # 3. RUN REGRESSION FOR SPECIES
# #   A) RUN IN LEVELS ACROSS FULL SAMPLE, DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE
# # -----------------------------------------------------------------------------
# # the original list: CR, CU, EC, FE, MN, NI, OC, PB, S, V, ZN
# full_sampPM_regMF = feols(c(AL,AS,BR, CA, CL, CHL,CR, CU, EC, FE, 
#                             K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
#                             S,  SE, SI, SO4, SR, TI, V,  ZN)
#                           ~ smokePM + nonsmokePM_MF | 
#                             monitor_month + year, reg_df, cluster = 'site_id') 
# 
# # etable(full_sampPM_regMF) # look at regression results
# 
# ## calculate 95 CI%
# CIs <- confint(
#   full_sampPM_regMF) %>% 
#   rename(species = 'lhs',
#          pm_type = 'coefficient',
#          CI25 = '2.5 %',
#          CI975 = '97.5 %') %>% 
#   mutate(measure = 'MF') %>% 
#   dplyr::select(-id)
# 
# # get coefficients and prepare for plotting
# full_sampPM_coeffsMF <- coeftable(full_sampPM_regMF) %>% 
#   rename(pval = 'Pr(>|t|)',
#          se = 'Std. Error',
#          species = 'lhs',
#          pm_type = 'coefficient') %>% 
#   mutate(pval = round(pval, digits = 3)) %>% 
#   # get pvalues
#   mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
#   dplyr::select(-id)  %>% 
#   left_join(parameter_categories, by = 'species') %>% 
#   mutate(species_type = fct_relevel(species_type,
#                                     c("Alkaline-earth metals", "Alkali metals", 
#                                       "Transition metals", "Metalloids", "Other metals", 
#                                       "Nonmetals",  "Halogens", "Organics"))
#   ) %>% 
#   mutate(measure = 'MF') %>% 
#   left_join(CIs, by = c('species', 'pm_type', 'measure'))
# 
# 
# # merge sample avg for each species and divide each species' betas by full sample avg for each species
# full_samp_PMcoeffs_normalized <- full_sampPM_coeffsMF %>% 
#   left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
#   mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
#          norm_CI25 = CI25/avg_nonsmoke_spec_conc,
#          norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
#   filter(pm_type == 'smokePM')
# 
# 
# # plot coefficients for speciation at the avg monitor which tells us how much 
# # of a chemical species is in 1 ug/m3 of smoke PM2.5 and nonsmoke PM2.5
# # --------------------------------------------------------------------------------
# # plot percent change for all species
# # --------------------------------------------------------------------------------
# pct_change_samp_reg_plot <- ggplot(full_samp_PMcoeffs_normalized, 
#                                    aes(x = species_long,
#                                        y = 100*norm_est, 
#                                        color=species_type, 
#                                    )) +
#   geom_point(size=4, alpha = 0.6, stat = "identity") +
#   geom_linerange(aes(ymin = (100*norm_CI25), 
#                      ymax = (100*norm_CI975)), stat = "identity") +
#   scale_x_discrete(limits = c(
#     # "Organics"
#     "Organic Carbon (OC)", "Elemental Carbon (EC)",
#     # "Halogens"
#     "Bromine (Br)", "Chlorine (Cl)", "Chloride (Chl)", 
#     #  "Nonmetals"
#     "Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Phosphorus (P)", "Selenium (Se)", 
#     # "Other metals"
#     "Titanium (Ti)", "Aluminum (Al)", "Lead (Pb)",
#     # "Metalloids"
#     "Silicon (Si)", "Arsenic (As)",
#     # "Transition metals"
#     "Zinc (Zn)", "Manganese (Mn)",  "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
#     # "Alkali metals"
#     "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
#     # "Alkaline-earth metals"
#     "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)" 
#   )) +
#   scale_color_manual(values= spec_pal) +
#   coord_flip() +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
#   labs(y = expression(paste('% change relative to average nonsmoke day concentration')),
#        x = 'Species',
#        #color = 'Species Category',
#        title = expression(paste("% change in concentration for 1 ug/", m^3, " increase in smoke ", "PM"["2.5"]))) + 
#   theme_light() + 
#   theme(panel.border = element_blank(), 
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         title= element_text(size=12, face='bold'),
#         axis.title.x = element_text(size=11, face = 'plain'),
#         axis.title.y = element_text(size=11, face = 'plain')) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50') +
#   guides(color = 'none')
# pct_change_samp_reg_plot
