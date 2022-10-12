source("script/functions.R")
source("script/packages.R")

fire <- read_rds("data/hms/hms_fires.rds")

# read in monitor locations
improve_loc <- read_csv("data/IMPROVE_CSN_sites.csv")
improve_loc_sf <- sfheaders::sf_point(improve_loc[,c("Longitude","Latitude")]) %>% mutate(Site = improve_loc$Site)


i  <- 1
firedist <- list()

for(i in 1:length(fire)){

  if(nrow(fire[[i]])>0){
  fire[[i]] <- fire[[i]][unlist(lapply(fire[[i]]$geometry, function(x){is.nan(x)[1]}))==F,]
  test <- st_nn(x = improve_loc_sf, y = fire[[i]],k = 1, sparse = T, returnDist = T)[[2]]%>% unlist()


firedist[[i]] <- data.frame(improve_loc, km2fire = round(test*111), date = as.Date(names(fire)[i], format = "%Y%m%d"))
  }else{
    firedist[[i]] <- data.frame(improve_loc, km2fire = 1e6, date = as.Date(names(fire)[i], format = "%Y%m%d"))
    
    
              }
}


firedist <- data.frame(rbindlist(firedist))
write_rds(firedist, file = "data/clean/improve_daily_dist_fire.rds", compress = "gz")



smoke <- read_rds("data/clean/improve_daily_hms_smoke.rds") %>% left_join(firedist)
smoke$year <- year(smoke$date)
smoke$month <- month(smoke$date)
smoke$day <- mday(smoke$date)





write_rds(smoke, file = "data/clean/improve_daily_hms_smoke_fire.rds", compress = "gz")




us <- tigris::states(); us <- us %>% filter(NAME %in% c("Puerto Rico","Guam","Alaska","Hawaii","American Samoa","Commonwealth of the Northern Mariana Islands","United States Virgin Islands") == F)
ave <- smoke %>% group_by(State, Longitude, Latitude) %>% summarise(low_count = mean(low_count), med_count = mean(med_count), high_count= mean(high_count), dist2fire = mean(km2fire, na.rm = T))


pal <- fields::two.colors(start = "yellow",middle = 'orange', end = 'red3')
int1 <- classInt::classIntervals(ave$low_count)
int2 <- classInt::classIntervals(ave$med_count)
int3 <- classInt::classIntervals(ave$high_count)
col1 <- findColours(int1, pal)
col2 <- findColours(int2, pal)
col3 <- findColours(int3, pal)

pal <- fields::two.colors(start = "red3",middle = 'orange', end = 'yellow')

int4 <- classInt::classIntervals(ave$dist2fire )
col4 <- findColours(int4, pal)

png(file = "figures/smoke_fire_explore_plots.png",width = 800, height = 800)
par(mfrow = c(2,2))
par(mar = c(0,0,4,0))
plot(st_geometry(us), lwd = 0.05)
points(ave$Longitude, ave$Latitude, pch = 21, bg = col1,col = 'white',cex=1.5)
mtext(side = 3, text = "Frequency of low density plumes (red is more frequent)")

plot(st_geometry(us), lwd = 0.05)
points(ave$Longitude, ave$Latitude, pch = 21, bg = col2,col = 'white',cex=1.5)
mtext(side = 3, text = "Frequency of med density plumes (red is more frequent)")

plot(st_geometry(us), lwd = 0.05)
points(ave$Longitude, ave$Latitude, pch = 21, bg = col3,col = 'white',cex=1.5)
mtext(side = 3, text = "Frequency of high density plumes (red is more frequent)")

plot(st_geometry(us), lwd = 0.05)
points(ave$Longitude, ave$Latitude, pch = 21, bg = col4,col = 'white',cex=1.5)
mtext(side = 3, text = "Distance to nearest fire (red is closer)")
dev.off()



