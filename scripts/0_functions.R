# helper functions ----
`%not_in%` <- purrr::negate(`%in%`)


#add a function for adding alphas to colors	
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}	




con2aqi <- function (pollutant, con, type = NULL) 
{
  if (!any(pollutant %in% c("o3", "pm25", "pm10", "co", "so2", 
                            "no2"))) {
    stop("Please input the available pollutant!")
  }
  con2aqi_single = function(s_pollutant, s_con, s_type) {
    pollutant = s_pollutant
    con = s_con
    type = s_type
    il = c(0, 51, 101, 151, 201, 301, 401)
    ih = c(50, 100, 150, 200, 300, 400, 500)
    if (pollutant == "o3") {
      if (type == "8h") {
        bpl = c(0, 0.055, 0.071, 0.086, 0.106)
        bph = c(bpl[-1] - 0.001, 0.2)
        if (con > tail(bph, 1)) {
          stop("Please confirm that the unit of OZONE is ppm!")
        }
      }
      if (type == "1h") {
        bpl = c(0, 0, 0.125, 0.165, 0.205, 0.405, 0.505)
        bph = c(bpl[-1] - 0.001, 0.604)
        if (con > tail(bph, 1)) {
          stop("Please confirm that the unit of OZONE is ppm!")
        }
      }
    }
    if (pollutant == "pm25") {
      bpl = c(0, 12.1, 35.5, 55.5, 150.5, 250.5, 350.5)
      bph = c(bpl[-1] - 0.1, 1000)
      if (con > tail(bph, 1)) {
        stop("Please confirm that the unit of PM2.5 is microgram per cubic metres!")
      }
    }
    if (pollutant == "pm10") {
      bpl = c(0, 55, 155, 255, 355, 425, 505)
      bph = c(bpl[-1] - 1, 604)
      if (con > tail(bph, 1)) {
        stop("Please confirm that the unit of PM10 is microgram per cubic metres!")
      }
    }
    if (pollutant == "co") {
      bpl = c(0, 4.5, 9.5, 12.5, 15.5, 30.5, 40.5)
      bph = c(bpl[-1] - 0.1, 50.4)
      if (con > tail(bph, 1)) {
        stop("Please confirm that the unit of CO is ppm!")
      }
    }
    if (pollutant == "so2") {
      bpl = c(0, 36, 76, 186, 305, 605, 805)
      bph = c(bpl[-1] - 1, 1004)
      if (con > tail(bph, 1)) {
        stop("Please confirm that the unit of SO2 is ppb!")
      }
    }
    if (pollutant == "no2") {
      bpl = c(0, 54, 101, 361, 650, 1250, 1650)
      bph = c(bpl[-1] - 1, 2049)
      if (con > tail(bph, 1)) {
        stop("Please confirm that the unit of SO2 is ppb!")
      }
    }
    rank = which(con <= bph)[1]
    aqi = ceiling((ih[rank] - il[rank])/(bph[rank] - bpl[rank]) * 
                    (con - bpl[rank]) + il[rank])
    names(aqi) = "aqi"
    return(aqi)
  }
  if (length(con) > 1) {
    rst = NULL
    for (i in 1:length(con)) {
      rst = c(rst, con2aqi_single(s_pollutant = pollutant, 
                                  s_con = con[i], s_type = type))
    }
    return(rst)
  }
  if (length(con) == 1) {
    rst = con2aqi_single(s_pollutant = pollutant, s_con = con, 
                         s_type = type)
    names(rst) = "aqi"
    return(rst)
  }
}
