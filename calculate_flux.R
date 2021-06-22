source("find_slopes.R")

calculate_flux_umol <- function(diameter, height, misc_vol, mean_t, slope){
  d  <- diameter        #diameter chamber [cm]
  h  <- height          #height chamber [cm]
  vadd <- misc_vol      #additional volume through collar etc [m^3]
  T1 <- mean_t +273.15  #average air temperature in relevant time window [K]
  m  <- slope           #delta C/delta t [ppm/s]
  
  vc  <- (pi*(d/2)^2*h)/ 10^6 #Volume [m^3]: area[cm^2](pi*diameter/4^2)*height[cm]*10^-6
  v <- vc+vadd                #total volume [m^3]
  v2 <- 22.4*10^-3            #molar volume for ideal gas at 273 K [m^3/mol]
  T2 <- 273.15                #standard temperature [K]
  A  <- pi*(d/2)^2*10^-4     #footprint area chamber [m^2]
  
  # calculate flux
  return(m * (v/(v2*(T2/T1))) * 1/A) #soil CO2 flux [µmol/(m^2*s)]
}




### ⇑ Definitions ⇑ | ⇓ Script body ⇓

filename = "kathy_data.csv"


# Csv file should have columns with the following exact names:
# time  co2  temp  humidity	location	diameter	height	misc_volume  notes
# And values should be provided as the following units
# time  ppm  C     %        text      cm        cm      m^3          text
data = read.csv(filename, dec=".", sep=";")

chamber_diameter = data$diameter[[1]]    # [cm]
chamber_height = data$height[[1]]        # [cm]
chamber_misc_vol = data$misc_volume[[1]] # [m^3]

if("location" %in% colnames(data)){
  location = data$location[[1]]
} else {
  location = "unspecified"
}

if("notes" %in% colnames(data)){
  notes = data$notes[[1]]
} else {
  notes = ""
}

# parse time
data$time <- as.POSIXct(strptime(data$time, "%d.%m.%y %H:%M:%S"))

# Add elapsed seconds column for linear regression
start_time = min(data$time)
data$seconds <- as.numeric(data$time - start_time)


source("find_slopes.R")
spans = find_slopes_periodic(data$seconds, data$co2)

fluxes <- spans %>% lapply(function(span){
  # Calculate mean temperature
  mean_t = mean(data[span$start:span$end, ]$temp)
  # Calculate flux
  flux_umol = calculate_flux_umol(
    chamber_diameter,
    chamber_height,
    chamber_misc_vol,
    mean_t,
    span$slope
  )
  span$flux_umol <- flux_umol
  span$mean_t <- mean_t
  return(span)
})

g <- ggplot(data, aes(x=seconds, y=co2)) + geom_point() + coord_cartesian(ylim=c(min(data$co2)-50,max(data$co2)+50))
for(i in 1:length(fluxes)){
  span = fluxes[[i]]
  g <- g + geom_abline(slope=span$slope, intercept=span$intercept, col="red") +
    geom_point(data=span$points, aes(x=x, y=y), col="red") +
    annotate("text", 
             x=last(span$points$x)-100, 
             y=last(span$points$y)+10, 
             label=sprintf("%f µmol/sm^2", span$flux_umol), 
             angle=70)
}

g













