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


# set up main regression for each species:
# Monitor-day Speciesconc = smokePMconc + monitorFE + yearFE + monthFE

# TRY RUNNING THE MAIN MODEL -----------------------------------------------------------------
# run a fixest regression
main_model_reg = feols(c(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4,
                    AL,AS,BR,CA,EC,OC,CHL,CL,CR,
                    CU,FE,K,MG, MN,MO,`NA`,NI,NO3,N2,
                    P,PB,RB,S,SE,SI,SOIL,SO4,SR,V,ZN,ZR) ~ smokePM_pred | SiteCode + year + month, reg_df)

# Save coefficients to table (this is necessary anyways when we make a plot)
# main_model_reg <- summary(main_model_reg, se = "hetero")
main_coeffs <- coefficients(main_model_reg) %>% 
  rename(coef = 'smokePM_pred')

# get standard errors
main_se <- fixest::se(main_model_reg) %>% 
  rename(se = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

# get pvalues
main_pval <- pvalue(main_model_reg) %>% 
  rename(pval = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

main_reg_res_df <- cbind(main_coeffs, main_se, main_pval) %>% 
  mutate(model= 'Main w/ Monitor, Year, Month FE')

# Pooled estimates over the entire country, if we use Marissa’s PM

# Sensitivity analysis using different FE:
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
moy_pval <- pvalue(moy_FE) %>% 
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
month_of_sample_pval <- pvalue(month_of_sample_FE) %>% 
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
no_negs_pval <- pvalue(no_negs_reg) %>% 
  rename(pval = 'smokePM_pred') %>% 
  dplyr::select(-id, -lhs)

no_negs_res_df <- cbind(no_negs_coeffs, no_negs_se, no_negs_pval) %>% 
  mutate(model= 'Main, no negatives')


# 5. Change negatives to detection limit values ----------------------------------------------------

##################################################################
# bind all results together and plot
all_results <- bind_rows(main_reg_res_df, moyFE_res_df, 
                         month_of_sample_res_df, no_negs_res_df) %>% 
  # drop the measures of carbon
  filter(!lhs %in% c('RC_PM2', 'MF', 'OC'))

all_results$model <- fct_rev(factor(all_results$model, 
                                     levels = c("Main w/ Monitor, Year, Month FE",
                                                "Monitor + year + month-of-year FE",
                                                "Monitor + month-of-sample FE",
                                                "Main, no negatives")))

col <- c('dodgerblue', 'forestgreen', 'salmon', 'plum')

# EXAMPLE PLOT
# pd <- position_dodge(0.1) # move them .05 to the left and right

plot <- ggplot(all_results, aes(x = lhs, y = coef, color=model)) +
  geom_point(aes(shape=model),size=3, alpha = 0.6, position=position_jitter(width=.5, seed = 50)) +
  scale_color_manual(name="Model",values=c("coral","steelblue", "forestgreen", 'plum')) +
  scale_shape_manual(name="Model",values=c(15, 17,18, 19)) +
  # scale_x_continuous("model", breaks=1:length(model), labels=model) +
  # scale_y_continuous("Speciation") +
  geom_linerange(aes(ymin = coef-se, ymax = coef+se),
                 position = position_jitter(width = 0.5, seed = 50)) +
  #geom_errorbar(aes(ymin=coef-se,ymax=coef+se),width=0.1,position=pd) +
  theme(axis.text.x=element_text(angle=-45)) +
  ggtitle("Regression Results") +
  labs(y = 'Coefficient',
       x = 'Chemical Species',
       color = 'Model',
       title = "Speciation of Smoke PM2.5") +
  coord_flip() +
  theme_bw() 
plot



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
