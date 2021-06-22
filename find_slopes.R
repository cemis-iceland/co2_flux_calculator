source("local_min_max.R")
library(pastecs)
library(dplyr)
library(robustbase)



find_slopes_periodic <- function(xs, ys, min_span=6){
  # We split the data into spans representing individual measurement periods
  
  # Smooth the data and find minima and maxima
  smoothed = filtfilt(butter(4, 0.3, "low"), ys - mean(ys)) # Complete guess
  extrema = turnpoints(smoothed)
  plot(smoothed)
  
  
  # Split the list into sections starting with minima and ending with maxima
  if (extrema$firstispeak) { # Drop first if it's a maximum
    extrema$tppos <- extrema$tppos[-1]
  }
  
  # If list can't be split into pairs cleanly drop the last index
  if (length(extrema$tppos) %% 2 != 0){ 
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
  slopes <- spans %>% apply(1, function(row){
    span <- tibble(x=xs[row["start"]:row["end"]], y=ys[row["start"]:row["end"]])
    
    # The start and end may have artifacts, so discard them
    range = max(span$y) - min(span$y)
    lower = min(span$y) + 0.10 * range # These percentages can be adjusted to 
    upper = max(span$y) - 0.00 * range # change how much data is discarded
    span <- span %>% filter(y > lower) %>% filter(y < upper)
    
    # Then simply perform a robust linear regression to find the slope of best fit
    fit = lmrob(y ~ x, span, method="KS2014")
    return(list(
      start=row["start"], 
      end=row["end"],
      points=span,
      fit=fit, 
      slope=fit$coefficients[[2]], 
      intercept=fit$coefficients[[1]]))
  })
  
  return(slopes)
}

