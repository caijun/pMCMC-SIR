get_Likelihoods <- function(theta, data, dataG, MCMC_params, X_start, X_traj_indexes) {
  # Set up for SMC
  X_traj_new <- NULL
  obsv_times <- length(data$y_vals)
  X_F <- array(0, dim = c(length(X_start), MCMC_params$J_particles, obsv_times))  # X_F is a 3D matrix for the state of each particle through time
  X_F[, , 1] <- repmat(X_start, 1, MCMC_params$J_particles)
  w_matrix <- matrix(0, nr = MCMC_params$J_particles, nc = obsv_times - 1)  # w_matrix is for particle weights
  k_matrix <- matrix(0, nr = MCMC_params$J_particles, nc = obsv_times - 1)  # k_matrix holds particle parent indexes
  
  for (n in 1:(obsv_times - 1)) {  # for each time interval
    t_n <- data$t_vals[n]  # start of time interval
    t_n_plus_1 <- data$t_vals[n + 1]  # end of time interval
    start_index <- dataG$indices[n]  # start index in dataG for time interval
    end_index <- dataG$indices[n + 1]  # end index in dataG for time interval
    # Run SMC
    list[X_F[, , n + 1], w_matrix[, n], X_traj, k_matrix[, n]] <- run_SMC(MCMC_params, X_F[, , n], theta, t_n, t_n_plus_1, data$y_vals[n + 1], 
                 dataG$lineages[start_index:end_index], dataG$coal_events[start_index:end_index], 
                 dataG$dt_ref[start_index:end_index], dataG$event_time[start_index:end_index])
    X_traj_new <- cbind(X_traj_new, X_traj[, -1])
  }
  
  # Calculate overall likelihood for each particle
  log_likelihood <- 0
  for (i in 1:(length(data$y_vals) - 1)) {
    sum_likelihood_J <- sum(w_matrix[, i] / MCMC_params$J_particles)
    log_likelihood <- log_likelihood + log(sum_likelihood_J)
  }
  p_theta <- log_likelihood  # prob of y_t given theta
  
  # Track one lineage's trajectory to sample from posterior of the latent variable X_I
  X_traj_1 <- repmat(X_start[2], 1, MCMC_params$J_particles)
  X_traj_final <- cbind(t(X_traj_1), X_traj_new)
  X_traj_sample <- get_sample_traj(k_matrix, X_traj_final, MCMC_params, X_traj_indexes)
  
  list(p_theta = p_theta, X_traj_sample = X_traj_sample)
  
#   ##############################################################################
#   # Calculate effective sampling size (ESS) to check for particle degeneracy
#   W <- colSums(w_matrix, 1)
#   for (i in 1:(obsv_times - 1)) {
#     norm_weights <- w_matrix[, i] / W[i]
#     sum_sqr_weights <- sum(norm_weights ^ 2)
#     ESS[i] <- 1 / sum_sqr_weights
#   }
#   
#   # View particle trajectories
#   I_traj <- matrix(X_F[2, , ], nr = MCMC_params$J_particles, obsv_times)
#   I_traj_var <- apply(I_traj, 2, var)
#   S_traj <- matrix(X[1, , ], MCMC_params$J_particles, obsv_times)
#   ##############################################################################
}