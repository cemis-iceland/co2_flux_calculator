source("find_slopes.R")
library("tidyverse")
library("lubridate")
library("ggstance")

# used in calculate_flux
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
  
  return(m * (v/(v2*(T2/T1))) * 1/A) #soil CO2 flux [µmol/(m^2*s)]
}

# Main function
# Takes a data frame containing the following columns:
#   seconds : the number of seconds from the start of the measurement
#   co2     : the measured concentration of co2 in ppm
#   temp    : the measured temperature in degrees Celsius
# And takes the following parameters:
#   chamber_diameter    : The diameter of the chamber (used to calculate volume)
#   chamber_height      : The height of the chamber (used to calculate volume)
#   chamber_misc_volume : Volume of the collar and any non-cylindrical things.
# Returns a 
calculate_fluxes <- function(
  data, # seconds, co2, temp
  chamber_diameter,
  chamber_height,
  chamber_misc_vol
){
  # Find measurement periods (rising edges in concentration)
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

  return(fluxes)
}

plot_fluxes <- function(data, fluxes){
  g <- ggplot(data, aes(x=seconds, y=co2)) + geom_point() + coord_cartesian(ylim=c(min(data$co2)-50,max(data$co2)+50))
  for(i in 1:length(fluxes)){
    span = fluxes[[i]]
    g <- g + geom_abline(slope=span$slope, intercept=span$intercept, col="red") +
      geom_point(data=span$points, aes(x=x, y=y), col="red") +
      annotate("text", 
              x=last(span$points$x)-100, 
              y=last(span$points$y)+10, 
              label=sprintf("%f µmol/cm^2", span$flux_umol), 
              angle=70)
  }
  return(g)
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

# Correct for humidity (water vapor dilution)
data$co2 <- data$co2/(1-(0.01*(data$humidity)*((exp(34.494-(4924.99/(data$temp+237.1))))/(data$temp+105)^1.57)*10^-3)/101)

# Correct for temperature 
data$co2 <- (273.15+data$temp)/273.15*data$co2

# End result: ppm of dry molar co2 at 25°C, 1atm

fluxes <- calculate_fluxes(
  data = data,
  chamber_diameter = chamber_diameter,
  chamber_height = chamber_height,
  chamber_misc_vol = chamber_misc_vol
  )

plot_fluxes(data, fluxes)

# save_fluxes_to_csv

row_format_fluxes <- map_dfr(fluxes, function(flux){
  return(list(
    Start_time=data[flux$start,]$time,
    End_time=data[flux$end,]$time,
    Mean_temperature=flux$mean_t,
    Slope=flux$slope,
    Flux_umol=flux$flux_umol))
})

# Get filename without extension
filename_wo_extension = strsplit(filename, ".", fixed=TRUE)[[1]][1]

# Save csv
write_csv(row_format_fluxes, paste(filename_wo_extension, "-fluxes.csv", sep=""))

# Save plot of flux over time
fluxvtime <- ggplot(row_format_fluxes, aes(x=Start_time + (End_time - Start_time) / 2, y=Flux_umol)) +
  geom_smooth(se=F, color="grey", method = "loess", method.args = list(family="symmetric")) +
  geom_point() +
  geom_linerangeh(aes(xmin=Start_time, xmax=End_time)) +
  scale_x_datetime(breaks = pretty_dates(row_format_fluxes$Start_time, 3)) +
  labs(title=paste("Flux over time for", filename), x="Time (middle of measurement)", y="Flux (µmol/scm^2)")
fluxvtime

png(filename=paste(filename_wo_extension, "-flux_over_time.png", sep=""))
print(fluxvtime)
dev.off()











