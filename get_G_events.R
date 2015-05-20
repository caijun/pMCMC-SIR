get_G_events <- function(MCMC_params, sample_times, sample_sizes, coal_times, obsv_times) {
  # set up data structure for genealogical data
  # major update 100410: data is now structured by dt interval not events
  
  # All times shoud be in months since tEnd
  sample_times <- sort(sample_times)
  coal_times <- coal_times[coal_times < max(obsv_times)]  # remove coal times that occur before t0
  coal_times <- sort(coal_times)
  obsv_times <- sort(obsv_times)
  lineages <- 0
  tEnd <- max(obsv_times)
  dt_times <- seq(0, tEnd, by = MCMC_params$dt)
  
  event_times <- c(coal_times, obsv_times, sample_times, dt_times)
  event_times <- sort(unique(event_times))
  
  G_lineages <- rep(0, length(event_times))
  G_coal_events <- rep(0, length(event_times))
  G_dt_ref <- rep(0, length(event_times))
  n <- 0
  G_indices <- rep(0, length(event_times))
  for (i in 1:length(event_times)) {
    time <- event_times[i]
    
    # Update G_lingeages
    if (time %in% sample_times) {
      loc <- sample_times == time
      lineages <- lineages + sample_sizes[loc]
    }
    
    # Update G_coal_events
    if (time %in% coal_times) {
      G_coal_events[i] <- 1
      coal_locs <- coal_times == time  # modified to handle more than one coal event at one time
      lineages <- lineages - length(coal_locs)
    }
    
    if (time %in% obsv_times) {
      n <- n + 1
      G_indices[n] <- i
    }
    
    G_lineages[i] <- lineages
    
    # Find closest dt time
    diff <- time -dt_times
    abs_diff <- abs(diff)
    dt_loc <- which.min(abs_diff)
    G_dt_ref <- dt_loc
  }
  
  G_lineages_dt <- rep(0, length(dt_times))
  for (j in 1:length(dt_times)) {
    time <- dt_times[j]
    loc <- event_times == time
    n <- n + 1
    G_lineages_dt[j] <- G_lineages[loc]
  }
  
  # Convert times to times since t0
  G_lineages <- rev(G_lineages)
  G_coal_events <- rev(G_coal_events)
  G_lineages_dt <- rev(G_lineages_dt)
  G_indices <- sort(length(event_times) - G_indices + 1)
  G_dt_ref <- sort(length(dt_times) - G_dt_ref + 1)
  event_times <- sort(tEnd - event_times)
  
  list(G_lineages = G_lineages, 
       G_coal_events = G_coal_events, 
       G_lineages_dt = G_lineages_dt, 
       G_indices = G_indices, 
       G_dt_ref = G_dt_ref, 
       event_times = event_times)
}