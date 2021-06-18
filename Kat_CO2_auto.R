#clear workspace from any previous entry
rm(list = ls())

#loading packages
library(ggplot2)
library(Cairo)
library(zoo)
library(dynlm)
library(gdata)
library(stringr)


#set working directory


###loop
list.files()
###load the two necessary scripts
source("local_min_max.R") ###from package spatialEco https://github.com/jeffreyevans/spatialEco
source("CO2M_Licor.R")
source("CO2M.R")

#set the path to your input directory containing your csv files
  inputDir <- "data/"

#run the function, flags are:
#### inputDir = path to data
#### first_m  = how many measurment to exclude at the start of the measurment
#### d        = inner diameter chamber [cm]
#### h        = height chamber [cm]
#### vadd     = additional volume through collar etc [m^3]
###span       = the smoothing parameter from fine (small values) to coarse (high value) this depend of the data so it needs to be empirically tested

###Use CO2M_Licor for licor data or CO2M for the sensor to accomodate for difference in Time formatting
#CO2M_Licor(inputDir = inputDir,first_m = 1, agd = 84.96,vadd = 0.00139339,span=0.3)

###Use CO2_M for your data02
CO2M(inputDir = inputDir,first_m = 0, d = 15.2,h = 4.2,vadd = 0.000335, span = 0.1)

     