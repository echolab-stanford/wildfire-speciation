#library
library(ggplot2)
library(dplyr)
library(lubridate)
library(scales)
library(gridExtra)
library(ggthemes)
library(anytime)
library(tidyr)
library(fixest)
library(plyr)
library(reshape2)
library(ggpmisc)

#setup
setwd("C:/Users/ayako/Documents/StanfordPhD/EchoLab/wildfire-speciation/scripts")
data_path <- "G:/Shared drives/echolab data/wildfire_speciation/intermediate/"
fig_path <- "C:/Users/ayako/Documents/StanfordPhD/EchoLab/wildfire-speciation/figures"

#function to calculate rmse
rmse2 <- function(m, o, na.rm = TRUE){
  res <- sqrt(mean((m-o)^2, na.rm = na.rm)) # o is true value
  return(res)
}


#import data
data <- read.csv(paste(data_path, "pm_plume_speciation_at_sites_tempAK.csv", sep = ''))
#file includes row names, so drop the column 
#data <- data[,-1]
data$smoke_day <- as.factor(data$smoke_day)
#add a short region name
data <- data %>% 
  mutate(abb_region = case_when(
    region == "rocky_mountain" ~ "RM",
    region =="midwest" ~ "MW",
    region =="southwest" ~ "SW",
    region =="noncontiguous" ~ "Non",
    region =="pacific" ~ "PAC",
    region =="northeast" ~ "NE",
    region =="southeast" ~ "SE"
  ))

#convert Date to Date class
data$Date <- as.Date(data$Date)
class(data$Date)

#for all monitors, correlation with smoke PM, monitor FE + year FE + month of year FE
df <- data %>% select(abb_region, SiteCode, year, month, smokePM_pred, MF, RC_PM2, NH4, 
                      tot_metals, ammNO3, ammSO4, AL, AS, BR, 
                      CA, EC, OC, CHL, CL, CR, CU, FE, K, MG, MN, MO, NA., NI,
                      NO3, N2, P, PB, RB, S, SE, SI, SOIL, SO4, SR, V, ZN, ZR)

vars <- data %>% select(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4, AL, AS, BR, 
                        CA, EC, OC, CHL, CL, CR, CU, FE, K, MG, MN, MO, NA., NI,
                        NO3, N2, P, PB, RB, S, SE, SI, SOIL, SO4, SR, V, ZN, ZR
                        )%>% names()

#log transform
df$smokePM_pred <- log(df$smokePM_pred)
df[vars] <- log(df[vars])

#pool all monitors, and monitor + year + month of year FE for each region
#get within R2 for each variable

#filter data for each region

df_RM <- df %>% dplyr::filter(abb_region == "RM")
df_MW <- df %>% dplyr::filter(abb_region == "MW")
df_SW <- df %>% dplyr::filter(abb_region == "SW")
df_Non <- df %>% dplyr::filter(abb_region == "Non")
df_PAC <- df %>% dplyr::filter(abb_region == "PAC")
df_NE <- df %>% dplyr::filter(abb_region == "NE")
df_SE <- df %>% dplyr::filter(abb_region == "SE")

#run for each region

####RM####
datalist <- list()
for(i in vars) {
    formula <- as.formula(paste(i, " ~ smokePM_pred | SiteCode + year + month"))
    within_R2 <- r2(feols(formula, data = df_RM), "war2")
    datalist[[i]] <- within_R2
  }
#create a dataframe
r2 <-  ldply (datalist, data.frame)
colnames(r2) <- c('species','withinR2')
r2$withinR2 <- round(r2$withinR2, 2)

#correlation plots for each species
for(i in vars) {
  r2_var <- r2 %>% filter(species == i)
  adj_r2 <- r2_var$withinR2
g <- ggplot(df_RM, aes(smokePM_pred)) + 
  geom_point(aes(y = get(i), color = "All"), alpha = 0.3) + 
  geom_smooth(aes(y = get(i), color = "All"), method = lm) + 
  stat_poly_eq(aes(y = get(i), color = "All"), formula = y ~ x, parse = T) + 
  #geom_abline() + 
  scale_color_manual(values = "darkblue") + 
  labs(x = "Smoke PM", y = i, title = paste("Correlation plot (log-log), within R2 =", adj_r2, sep = "")) + 
  theme_light() 
  print(g)
  ggsave(path = paste(fig_path,"/11_correlation_plots/RM",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

####MW####
datalist <- list()
for(i in vars) {
  formula <- as.formula(paste(i, " ~ smokePM_pred | SiteCode + year + month"))
  within_R2 <- r2(feols(formula, data = df_MW), "war2")
  datalist[[i]] <- within_R2
}
#create a dataframe
r2 <-  ldply (datalist, data.frame)
colnames(r2) <- c('species','withinR2')
r2$withinR2 <- round(r2$withinR2, 2)

#correlation plots for each species
for(i in vars) {
  r2_var <- r2 %>% filter(species == i)
  adj_r2 <- r2_var$withinR2
  g <- ggplot(df_MW, aes(smokePM_pred)) + 
    geom_point(aes(y = get(i), color = "All"), alpha = 0.3) + 
    geom_smooth(aes(y = get(i), color = "All"), method = lm) + 
    stat_poly_eq(aes(y = get(i), color = "All"), formula = y ~ x, parse = T) + 
    #geom_abline() + 
    scale_color_manual(values = "darkblue") + 
    labs(x = "Smoke PM", y = i, title = paste("Correlation plot (log-log), within R2 =", adj_r2, sep = "")) + 
    theme_light() 
  print(g)
  ggsave(path = paste(fig_path,"/11_correlation_plots/MW",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

####NE####
datalist <- list()
for(i in vars) {
  formula <- as.formula(paste(i, " ~ smokePM_pred | SiteCode + year + month"))
  within_R2 <- r2(feols(formula, data = df_NE), "war2")
  datalist[[i]] <- within_R2
}
#create a dataframe
r2 <-  ldply (datalist, data.frame)
colnames(r2) <- c('species','withinR2')
r2$withinR2 <- round(r2$withinR2, 2)

#correlation plots for each species
for(i in vars) {
  r2_var <- r2 %>% filter(species == i)
  adj_r2 <- r2_var$withinR2
  g <- ggplot(df_NE, aes(smokePM_pred)) + 
    geom_point(aes(y = get(i), color = "All"), alpha = 0.3) + 
    geom_smooth(aes(y = get(i), color = "All"), method = lm) + 
    stat_poly_eq(aes(y = get(i), color = "All"), formula = y ~ x, parse = T) + 
    #geom_abline() + 
    scale_color_manual(values = "darkblue") + 
    labs(x = "Smoke PM", y = i, title = paste("Correlation plot (log-log), within R2 =", adj_r2, sep = "")) + 
    theme_light() 
  print(g)
  ggsave(path = paste(fig_path,"/11_correlation_plots/NE",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}



####Non####
datalist <- list()
for(i in vars) {
  formula <- as.formula(paste(i, " ~ smokePM_pred | SiteCode + year + month"))
  within_R2 <- r2(feols(formula, data = df_Non), "war2")
  datalist[[i]] <- within_R2
}
#create a dataframe
r2 <-  ldply (datalist, data.frame)
colnames(r2) <- c('species','withinR2')
r2$withinR2 <- round(r2$withinR2, 2)

#correlation plots for each species
for(i in vars) {
  r2_var <- r2 %>% filter(species == i)
  adj_r2 <- r2_var$withinR2
  g <- ggplot(df_Non, aes(smokePM_pred)) + 
    geom_point(aes(y = get(i), color = "All"), alpha = 0.3) + 
    geom_smooth(aes(y = get(i), color = "All"), method = lm) + 
    stat_poly_eq(aes(y = get(i), color = "All"), formula = y ~ x, parse = T) + 
    #geom_abline() + 
    scale_color_manual(values = "darkblue") + 
    labs(x = "Smoke PM", y = i, title = paste("Correlation plot (log-log), within R2 =", adj_r2, sep = "")) + 
    theme_light() 
  print(g)
  ggsave(path = paste(fig_path,"/11_correlation_plots/Non",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

####PAC####
datalist <- list()
for(i in vars) {
  formula <- as.formula(paste(i, " ~ smokePM_pred | SiteCode + year + month"))
  within_R2 <- r2(feols(formula, data = df_PAC), "war2")
  datalist[[i]] <- within_R2
}
#create a dataframe
r2 <-  ldply (datalist, data.frame)
colnames(r2) <- c('species','withinR2')
r2$withinR2 <- round(r2$withinR2, 2)

#correlation plots for each species
for(i in vars) {
  r2_var <- r2 %>% filter(species == i)
  adj_r2 <- r2_var$withinR2
  g <- ggplot(df_PAC, aes(smokePM_pred)) + 
    geom_point(aes(y = get(i), color = "All"), alpha = 0.3) + 
    geom_smooth(aes(y = get(i), color = "All"), method = lm) + 
    stat_poly_eq(aes(y = get(i), color = "All"), formula = y ~ x, parse = T) + 
    #geom_abline() + 
    scale_color_manual(values = "darkblue") + 
    labs(x = "Smoke PM", y = i, title = paste("Correlation plot (log-log), within R2 =", adj_r2, sep = "")) + 
    theme_light() 
  print(g)
  ggsave(path = paste(fig_path,"/11_correlation_plots/PAC",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}


####SE####
datalist <- list()
for(i in vars) {
  formula <- as.formula(paste(i, " ~ smokePM_pred | SiteCode + year + month"))
  within_R2 <- r2(feols(formula, data = df_SE), "war2")
  datalist[[i]] <- within_R2
}
#create a dataframe
r2 <-  ldply (datalist, data.frame)
colnames(r2) <- c('species','withinR2')
r2$withinR2 <- round(r2$withinR2, 2)

#correlation plots for each species
for(i in vars) {
  r2_var <- r2 %>% filter(species == i)
  adj_r2 <- r2_var$withinR2
  g <- ggplot(df_SE, aes(smokePM_pred)) + 
    geom_point(aes(y = get(i), color = "All"), alpha = 0.3) + 
    geom_smooth(aes(y = get(i), color = "All"), method = lm) + 
    stat_poly_eq(aes(y = get(i), color = "All"), formula = y ~ x, parse = T) + 
    #geom_abline() + 
    scale_color_manual(values = "darkblue") + 
    labs(x = "Smoke PM", y = i, title = paste("Correlation plot (log-log), within R2 =", adj_r2, sep = "")) + 
    theme_light() 
  print(g)
  ggsave(path = paste(fig_path,"/11_correlation_plots/SE",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}


####SW####
datalist <- list()
for(i in vars) {
  formula <- as.formula(paste(i, " ~ smokePM_pred | SiteCode + year + month"))
  within_R2 <- r2(feols(formula, data = df_SW), "war2")
  datalist[[i]] <- within_R2
}
#create a dataframe
r2 <-  ldply (datalist, data.frame)
colnames(r2) <- c('species','withinR2')
r2$withinR2 <- round(r2$withinR2, 2)

#correlation plots for each species
for(i in vars) {
  r2_var <- r2 %>% filter(species == i)
  adj_r2 <- r2_var$withinR2
  g <- ggplot(df_SW, aes(smokePM_pred)) + 
    geom_point(aes(y = get(i), color = "All"), alpha = 0.3) + 
    geom_smooth(aes(y = get(i), color = "All"), method = lm) + 
    stat_poly_eq(aes(y = get(i), color = "All"), formula = y ~ x, parse = T) + 
    #geom_abline() + 
    scale_color_manual(values = "darkblue") + 
    labs(x = "Smoke PM", y = i, title = paste("Correlation plot (log-log), within R2 =", adj_r2, sep = "")) + 
    theme_light() 
  print(g)
  ggsave(path = paste(fig_path,"/11_correlation_plots/SW",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}
