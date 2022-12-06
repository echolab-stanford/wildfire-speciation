#library
library(ggplot2)
library(dplyr)
library(lubridate)
library(scales)
library(gridExtra)
library(ggthemes)
library(anytime)
library(tidyr)

#setup
setwd("C:/Users/ayako/OneDrive/Documents/echolab/wip-ayako/wildfire-speciation/scripts")
data_path <- "G:/Shared drives/echolab data/wildfire_speciation/intermediate/"
fig_path <- "C:/Users/ayako/OneDrive/Documents/echolab/wip-ayako/wildfire-speciation/figures"

#import data
data <- read.csv(paste(data_path, "pm_plume_speciation_at_sites.csv", sep = ''))
#file includes row names, so drop the column 
data <- data[,-1]
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

data$month_yr <- format(data$Date, "%Y-%m")
monthly <- data %>% select(month_yr, year, month, RC_PM2, NH4, tot_metals, ammNO3,ammSO4,
                        AL, AS, BR, CA, CR, CU, FE, K, MG, MN, NI,
                        NO3, N2, P, PB, RB, S, SE, SI, SOIL, SO4,
                        SR, V, ZN, ZR, MF)

vars <- data %>% select(RC_PM2, NH4, tot_metals, ammNO3,ammSO4,
                        AL, AS, BR, CA, CR, CU, FE, K, MG, MN, NI,
                        NO3, N2, P, PB, RB, S, SE, SI, SOIL, SO4,
                        SR, V, ZN, ZR, MF) %>% names()



#create monthly data, calculate monthly mean
monthly_gb <- monthly %>%
  group_by(month_yr) %>%
  summarise(across(everything(), mean, na.rm = T))

#convert from character to date variable
#monthly_gb$month_yr <- as.Date(as.character(monthly_gb$month_yr))
monthly_gb$year <- as.character(monthly_gb$year)
monthly_gb$month <- as.character(monthly_gb$month)
monthly_gb$month <- factor(monthly_gb$month, levels=c("1", "2", "3",
                                                      "4", "5", "6",
                                                      "7", "8", "9",
                                                      "10", "11", "12"))


#all region, time series of monthly means from 2006-2021
for(i in vars) {
  
  g <- ggplot(monthly_gb, aes(x = month_yr, y = get(i))) +
    geom_point(na.rm=TRUE, color="darkblue", size=1) +
    (scale_x_date(breaks=date_breaks("12 months"),
                  labels=date_format("%y"))) +
    ylab(i) +
    xlab("Year") +
    ggtitle("Monthly means, 2006-2021")
  print(g)
  ggsave(path = paste(fig_path,"/09_monthly_means",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

#all region, seasonality plot of monthly means from 2006-2021
ggplot(monthly_gb, aes(x = month, y = RC_PM2, group = year, color = year)) +
  geom_line()

#all region, time series of monthly means from 2006-2021
for(i in vars) {
  
  g <- ggplot(monthly_gb, aes(x = month, y = get(i), group = year, color = year)) +
    geom_line() +
    ylab(i) +
    xlab("Month") +
    ggtitle("Monthly means, 2006-2021")
  print(g)
  ggsave(path = paste(fig_path,"/10_monthly_means_seasonality",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

