source("local_min_max.R")
library(pastecs)
library(dplyr)
library(robustbase)

#process input

filename = "05-07-21_DownloadedLog.csv"

data = read.csv(filename, skip=1)
# next 2 lines deal with bad formatting
names(data) <- c("time", "co2", "", "temp", "", "humidity")
data <- subset(data, select = c("time", "co2", "temp", "humidity"))

# parse time
data$time <- as.POSIXct(strptime(data$time, "%d/%m/%y %H:%M:%S"))

ggplot(data, aes(x=time, y=co2)) + geom_point()

# debug
#find_slopes_periodic <- function(xs, ys, min_span=5){
xs = data$time
ys = data$co2

min_span = 5
#/debug

# We split the data into spans representing individual measurement periods

# Smooth the data and find minima and maxima
smoothed = filtfilt(butter(4, 0.5, "low"), ys) # Complete guess
extrema = turnpoints(smoothed)

# Split the list into sections starting with minima and ending with maxima
if (extrema$firstispeak) { # Drop first if it's a maximum
  extrema$tppos <- extrema$tppos[-1]
}

# If list can't be split into pairs cleanly drop the last index
if (length(extrema$tppos %mod% 2 != 0)){ 
  extrema$tppos <- extrema$tppos[-length(extrema$tppos)]
}

# Split list into (minimum, maximum) pairs
pairs = split(extrema$tppos, ceiling(seq_along(extrema$tppos)/2)) 

# Compute length for each span and convert the list to a data frame
spans = bind_rows(pairs %>% lapply(function(span){ 
    return(list(start=span[1], end=span[2], len=span[2] - span[1]))
  }))

# Exclude short spans to eliminate noise
spans <- filter(spans, len > min_span)

# For each span, find the slope
plot(xs, ys)
for (i in 1:nrow(spans)) {
  # Actually sublist the data
  row <- spans[i, ]
  span <- tibble(x=xs[row$start:row$end], y=ys[row$start:row$end])
  
  # The start and end may have artifacts, so discard them
  range = max(span$y) - min(span$y)
  lower = min(span$y) + 0.10 * range # These percentages can be adjusted to 
  upper = max(span$y) - 0.10 * range # change how much data is discarded
  span <- span %>% filter(y > lower) %>% filter(y < upper)
  
  # Then simply perform a robust linear regression to find the slope of best fit
  
  
  points(span, col="red")
}

#}
