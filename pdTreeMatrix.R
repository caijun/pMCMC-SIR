pdTreeMatrix <- function(MCMC_params, theta, lineages, coal_events, event_times, 
                         dt_ref, X_traj, S_traj, P_traj) {
  # Calculate the likelihood of the genealogy given the particle trajectories for all subintervals in current time interval
  mu <- theta[1]
  gamma <- theta[2]
  R0_avg <- theta[3]
  beta <- R0_avg * (gamma + mu)
  alpha <- theta[4]
  
  seas_at_events <- event_times / 12 - floor(event_times / 12)
  beta_at_events <- beta * (1 + alpha * cos(2 * pi * seas_at_events))
  
  dt_ref <- dt_ref - (dt_ref[1] - 1)  # reference to which dt integration time point is closest to genealogical event
  # Get subsinterval lengths
  event_plus_one <- c(0, event_times)
  sub_times <- event_times - head(event_plus_one, -1)
  
  # Set up matrices
  coal_mat <- repmat(coal_events, MCMC_params$J_particles, 1)  # ones where coal events, zeros everywhere else
  no_coal_mat <- 1 - coal_mat  # ones where no coal events, zeros where coal events
  sub_times_mat <- repmat(sub_times, MCMC_params$J_particles, 1)
  
  # Calculate lambda for all particles (rate of coalescence)
  k_mat <- repmat(lineages, MCMC_params$J_particles, 1)  # k is the number of lineages
  k_binom_mat <- k_mat * (k_mat - 1) / 2  # binomial coefficient
  I_mat <- X_traj[, dt_ref]  # this is I for infectious diseases
  I_binom_mat <- I_mat ^ 2 / 2
  S_mat <- S_traj[, dt_ref]
  P_mat <- P_traj[, dt_ref]
  beta_mat <- repmat(beta_at_events, MCMC_params$J_particles, 1)
  lambda_mat <- (k_binom_mat / I_binom_mat) * (beta_mat * S_mat * I_mat / P_mat)
  
  # Get_likelihoods
  coal_prob <- lambda_mat * exp(-lambda_mat * sub_times_mat)
  no_coal_prob <- exp(-lambda_mat * sub_times_mat)
  
  m1 <- coal_prob[, -1] * coal_mat[, -ncol(coal_mat)]  # likelihoods for times with coal events
  m2 <- no_coal_prob[, -1] * no_coal_mat[, -ncol(no_coal_mat)]  # likelihoods for times with no coal events
  
  p_mat <- m1 + m2  # get probabilities for all times in  interval
  w_G <- apply(p_mat, 1, prod) # product over all subintervals
  w_G[is.na(w_G)] <- 0
}