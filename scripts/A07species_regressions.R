# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Jan 16, 2023
# Description: run regressions for each species

# loadd(clean_pm_spec_df, cache = drake_cache)

run_regressions <- function(clean_pm_spec_df) {

# set up dataframe for regressions
reg_df <- clean_pm_spec_df %>% 
  mutate_if(is.numeric, ~ifelse(. == Inf, NA, .)) %>% 
  mutate(SiteCode = as.factor(SiteCode),
         year =as.factor(year),
         month = as.factor(month),
         moy = as.factor(paste0(month,'-', year))) 



#######################################################################################
# RUN THE MAIN MODEL -----------------------------------------------------------------
#######################################################################################
# set up main regression for each species:
# Monitor-day Species conc = smokePMconc + monitorFE + yearFE + monthFE

# run a fixest regression
main_model_reg = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
                    AL,AS,BR,CA,EC,OC,CHL,CL,CR,
                    CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
                    P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + year + month, reg_df)

# Save coefficients to table (this is necessary anyways when we make a plot)
main_model_reg <- summary(main_model_reg, se = "hetero")
main_coeffs <- coefficients(main_model_reg) %>% 
  rename(coef = 'smokePM_pred')

# get standard errors
# main_se <- fixest::se(main_model_reg) %>% 
#   rename(se = 'smokePM_pred') %>% 
#   dplyr::select(-id, -lhs)
# 
# # get pvalues
main_pval <- pvalue(main_model_reg) %>%
  rename(pval = 'smokePM_pred') %>%
  dplyr::select(-id, -lhs) %>% 
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))
# 
# main_conf <- confint(main_model_reg) %>% 
#   rename(low_ci = `2.5 %`,
#          upper_ci = `97.5 %`) %>%
#   dplyr::select(-id, -lhs)

main_reg_res_df <- cbind(main_coeffs, main_pval) %>% # main_se, main_conf, main_pval) %>% 
  mutate(model= 'Main w/ Monitor, Year, Month FE') 

# plot avg monitor
main_reg_plot <- ggplot(main_reg_res_df, 
                        aes(x = reorder(lhs, coef), y = coef, color=sig)) +
  geom_point(size=3, shape = 15, alpha = 0.6) +
  scale_color_manual(values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  # geom_hline(data= mean, aes(yintercept = mean_val, col='grey')) +
  #geom_linerange(aes(ymin = coef-se, ymax = coef+se)) +
  # scale_y_continuous(
  #   trans =  'log10',
  #   breaks = trans_breaks("log10", function(x) 10^x),
  #   labels = trans_format("log10", math_format(10^.x))
  # ) +
  #coord_trans(y ='log10', x='log10') +
  theme(axis.text.x=element_text(angle=-45)) +
  ggtitle("Regression Results") +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Significance based on p < .05',
       title = "Speciation of Smoke PM2.5",
       subtitle = 'Monitor-day species conc = smokePM2.5 + monitor FE + year FE + month FE') +
  coord_flip() +
  theme_light() 
main_reg_plot

#######################################################################################
# REGIONAL ANALYSIS
#######################################################################################
# run a fixest regression
region_model_reg = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
                         AL,AS,BR,CA,EC,OC,CHL,CL,CR,
                         CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
                         P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + year + month, reg_df, fsplit = ~region)

region_model_reg <- summary(region_model_reg, se = "hetero")
region_coeffs <- coefficients(region_model_reg) %>% 
  rename(coef = 'smokePM_pred') %>% 
  mutate(model= 'Main model w/ Monitor, Year, Month FE *BY REGION*') 

# # get pvalues
region_pval <- pvalue(region_model_reg) %>%
  rename(pval = 'smokePM_pred') %>%
  dplyr::select(-id, -lhs, -sample, -sample.var) %>% 
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant'))

region_model_df <- cbind(region_coeffs, region_pval) %>% 
  dplyr::select(lhs, coef, pval, sig, sample)

# PLOT
all_region_plot <- ggplot(region_model_df, 
                       aes(x = reorder(lhs, coef), y = coef)) +
  geom_point(aes(shape=sig, color = sample), size=3, alpha = 0.6) +
  scale_color_manual(name="Region",
                     values=c("coral","steelblue", "forestgreen", 'plum', 'mediumpurple', 'goldenrod2', 'darkblue')) +
  scale_shape_manual(name="Significance", values=c(18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme(axis.text.x=element_text(angle=-45)) +
  ggtitle("Regression Results") +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Region',
       title = "Speciation of Smoke PM2.5 by Region") +
  coord_flip() +
  theme_light() 
all_region_plot


#######################################################################################
# SMOKE VS PM2.5 ANALYSIS
#######################################################################################
# calculate non-wildfire PM at the station level by just subtracting wildfire PM from station-measured total PM on that day
pm_types_df <- reg_df %>% 
  mutate(non_smokePM = ifelse(!is.na(MF), MF - smokePM_pred, NA)) %>% 
  filter(!is.na(non_smokePM))
  
pm_types_reg = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
                       AL,AS,BR,CA,EC,OC,CHL,CL,CR,
                       CU,FE,K,MG, MN,`NA`,NI,NO3,N2,
                       P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ non_smokePM | SiteCode + year + month, pm_types_df)

pm_types_reg <- summary(pm_types_reg, se = "hetero")
pm_coeffs <- coefficients(pm_types_reg) %>% 
  rename(coef = 'non_smokePM') %>% 
  mutate(model= 'Main model w/ Monitor, Year, Month FE *NONSMOKE*') 

# plot avg monitor
# non_smoke_reg_plot <- ggplot(pm_coeffs, 
#                         aes(x = lhs, y = coef, color=model)) +
#   geom_point(aes(shape=model),size=3, alpha = 0.6) +
#   scale_color_manual(name="Model",values=c("coral","steelblue", "forestgreen", 'plum')) +
#   scale_shape_manual(name="Model",values=c(15, 17,18, 19)) +
#   geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
#   scale_y_continuous(
#     trans = log10_trans(),
#     breaks = trans_breaks("log10", function(x) 10^x),
#     labels = trans_format("log10", math_format(10^.x))
#   ) +
#   theme(axis.text.x=element_text(angle=-45)) +
#   ggtitle("Regression Results") +
#   labs(y = 'Coefficient',
#        x = 'Chemical Species',
#        #color = 'Model',
#        title = "Speciation of Non-Smoke PM2.5") +
#   coord_flip() +
#   theme_light() 
# non_smoke_reg_plot

##############
both_smoke_nonsmoke <- pm_coeffs %>% 
  bind_rows(main_coeffs %>% 
              mutate(model = 'Main model w/ Monitor, Year, Month FE *SMOKE*'))

# plot avg monitor
pm_types_reg_plot <- ggplot(both_smoke_nonsmoke, 
                            aes(x = lhs, y = coef, color=model)) +
  geom_point(aes(shape=model),size=3, alpha = 0.6) +
  scale_color_manual(name="Model",values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(name="Model",values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme(axis.text.x=element_text(angle=-45)) +
  ggtitle("Regression Results") +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       #color = 'Model',
       title = "Speciation of Non-Smoke PM2.5 vs Smoke PM2.5") +
  coord_flip() +
  theme_light() 
pm_types_reg_plot

#######################################################################################
# Sensitivity analysis using different FE:
#######################################################################################
# 1. Monitor + year + month of year ----------------------------------------------------
moy_FE = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
                             AL,AS,BR,CA,EC,OC,CHL,CL,CR,
                             CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
                             P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + year + moy, reg_df)

# Save coefficients to table (this is necessary anyways when we make a plot)
moy_coeffs <- coefficients(moy_FE) %>% 
  rename(coef = 'smokePM_pred')

# get standard errors
moy_se <- fixest::se(moy_FE) %>% 
  rename(se = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

# get pvalues
moy_pval <- fixest::pvalue(moy_FE) %>% 
  rename(pval = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

moyFE_res_df <- cbind(moy_coeffs, moy_se, moy_pval) %>% 
  mutate(model= 'Monitor + year + month-of-year FE')

# 2. Monitor + month-of-sample ----------------------------------------------------
month_of_sample_FE = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
                 AL,AS,BR,CA,EC,OC,CHL,CL,CR,
                 CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
                 P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + month, reg_df)

# Save coefficients to table (this is necessary anyways when we make a plot)
month_of_sample_coeffs <- coefficients(month_of_sample_FE) %>% 
  rename(coef = 'smokePM_pred')

# get standard errors
month_of_sample_se <- fixest::se(month_of_sample_FE) %>% 
  rename(se = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

# get pvalues
month_of_sample_pval <- fixest::pvalue(month_of_sample_FE) %>% 
  rename(pval = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

month_of_sample_res_df <- cbind(month_of_sample_coeffs, month_of_sample_se, month_of_sample_pval) %>% 
  mutate(model= 'Monitor + month-of-sample FE')

# 3. Monitor * month of year + year----------------------------------------------------
# 4. Drop negatives -----------------------------------------------------------------
no_negs <- reg_df %>% 
  mutate_if(is.numeric, ~ifelse(. < 0, NA, .)) 

no_negs_reg = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
                             AL,AS,BR,CA,EC,OC,CHL,CL,CR,
                             CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
                             P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + year + month, reg_df)

# Save coefficients to table (this is necessary anyways when we make a plot)
no_negs_coeffs <- coefficients(no_negs_reg) %>% 
  rename(coef = 'smokePM_pred')

# get standard errors
no_negs_se <- fixest::se(no_negs_reg) %>% 
  rename(se = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

# get pvalues
no_negs_pval <- fixest::pvalue(no_negs_reg) %>% 
  rename(pval = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

no_negs_res_df <- cbind(no_negs_coeffs, no_negs_se, no_negs_pval) %>% 
  mutate(model= 'Main, no negatives')

# 5. Change negatives to detection limit values ----------------------------------------------------

##################################################################
# bind all results together and plot
all_results <- bind_rows(main_reg_res_df, moyFE_res_df, 
                         month_of_sample_res_df, no_negs_res_df) # %>% 
  # drop the measures of carbon
  # filter(!lhs %in% c('RC_PM2', 'MF', 'OC'))

all_results$model <- fct_rev(factor(all_results$model, 
                                     levels = c("Main w/ Monitor, Year, Month FE",
                                                "Monitor + year + month-of-year FE",
                                                "Monitor + month-of-sample FE",
                                                "Main, no negatives")))

col <- c('dodgerblue', 'forestgreen', 'salmon', 'plum')

# EXAMPLE PLOT
# plot avg monitor
all_reg_plot <- ggplot(all_results, 
                        aes(x = lhs, y = coef, color=model)) +
  geom_point(aes(shape=model),size=3, alpha = 0.6) +
  scale_color_manual(name="Model",values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(name="Model",values=c(15, 17,18, 19)) +
  geom_hline(yintercept = 10^0, linetype = "dashed", color = 'grey') +
  #geom_linerange(aes(ymin = coef-se, ymax = coef+se)) +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme(axis.text.x=element_text(angle=-45)) +
  ggtitle("Regression Results") +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Model',
       title = "Sensitivity Analysis") +
  coord_flip() +
  theme_light() 
all_reg_plot



# pd <- position_dodge(0.1) # move them .05 to the left and right
# plot <- ggplot(all_results, aes(x = lhs, y = log(coef), color=model)) +
#   geom_point(aes(shape=model),size=3, alpha = 0.6, position=position_jitter(width=.5, seed = 50)) +
#   scale_color_manual(name="Model",values=c("coral","steelblue", "forestgreen", 'plum')) +
#   scale_shape_manual(name="Model",values=c(15, 17,18, 19)) +
#   # scale_x_continuous("model", breaks=1:length(model), labels=model) +
#   # scale_y_continuous("Speciation") +
#   geom_linerange(aes(ymin = coef-se, ymax = coef+se),
#                  position = position_jitter(width = 0.5, seed = 50)) +
#   #geom_errorbar(aes(ymin=coef-se,ymax=coef+se),width=0.1,position=pd) +
#   theme(axis.text.x=element_text(angle=-45)) +
#   ggtitle("Regression Results") +
#   labs(y = 'Coefficient',
#        x = 'Chemical Species',
#        color = 'Model',
#        title = "Speciation of Smoke PM2.5") +
#   coord_flip() +
#   theme_bw() 
# plot



# create scatter of estimates w/standard error bars
# coeff_plot <- ggplot(all_results, 
#                        aes(x = lhs, y = coef)) +
#   geom_point(aes(x = lhs, y = coef, color = model), size = 1, position = "jitter", alpha = .4) +
#   scale_colour_manual(values= col) +
#   
#   # add poulation density base
#   geom_hline(yintercept=0, linetype="dashed", lwd=.5, color='grey') +
#   ggtitle("Regression Results") +
#   labs(y = 'Coefficient',
#        x = 'Chemical Species',
#        color = 'Model',
#        title = "Speciation of Smoke PM2.5") +
#   coord_flip() +
#   theme_bw() 
#   # theme(plot.title =element_text(size=14, face='bold')) +
#   # # theme(plot.subtitle = element_text(size=12)) +
#   # theme(axis.text=element_text(size=11),
#   #       axis.title=element_text(size=12,face="bold")) +
#   # + theme(legend.text=element_text(size=14))
# 
# coeff_plot

return(all_results)

} # end function
