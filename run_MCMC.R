run_MCMC <- function(X_start, data, dataG, MCMC_params, theta_now) {
  # X_traj is the trajectory of the latent variable I
  X_traj_indexes <- seq(1, ((length(data$t_vals) - 1) * (1 / MCMC_params$dt)) + 1, 
                        by = 1 / (MCMC_params$dt * 10))
  X_traj_length <- length(X_traj_indexes)  # (length(data$t_vals) - 1) * (1 / MCMC_params$dt) * 12 + 1
  
  # Get initial likelihoods, samples
  list[p_theta, X_traj_sample] <- get_Likelihoods(theta_now, data, dataG, MCMC_params, X_start, X_traj_indexes)
  p_theta_now <- p_theta
  X_traj_now <- X_traj_sample
  X_traj_samples <- X_traj_now
  theta_samples <- theta_now
  p_samples <- rep(0, MCMC_params$iterations)  # p_sample holds the marginal likelihoods
  p_samples[1] <- p_theta_now
  
  # Set up MCMC sampling
  log_m <- seq(MCMC_params$log_steps, MCMC_params$iterations, by = MCMC_params$log_steps)
  if (MCMC_params$iterations > MCMC_params$save_steps) {
    save_m <- seq(MCMC_params$save_steps, MCMC_params$iterations, by = MCMC_params$save_steps)
  } else {
    save_m <- NULL
  }
  proposals <- NULL
  accept <- 0
  
  ### This is the MCMC part###
  for (m in 2:MCMC_params$iterations) {
    tic
    MCMC_params$m <- m
    print(paste0("MCMC_step = ", m))
    
    # parameter proposals
    theta_new <- rmvnorm(n = 1, mean = theta_now, sigma = MCMC_params$theta_cov)
    while (min(theta_new) < 0 | theta_new[6] > 1) {
      # if any parameter is negative OR if rho is greater than one, propose new theta
      theta_new <- rmvnorm(n = 1, mean = theta_now, sigma = MCMC_params$theta_cov)
    }
    curr_proposal <- theta_new
    
    # Run SMC to get the marginal likelihood for theta_new
    list[p_theta, X_traj_sample] <- get_Likelihoods(theta_new, data, dataG, MCMC_params, X_start, X_traj_indexes)
    p_theta_new <- p_theta
    X_traj_new <- X_traj_sample
    p_samples[m] <- p_theta_new  # these are the marginal likelihoods
    
    # Accept or reject theta proposal
    a <- exp(p_theta_new - p_theta_now)
    z <- runif(1)
    if (z < a) {
      X_traj_now <- X_traj_new
      theta_now <- theta_new
      accept <- accept + 1
      p_theta_now <- p_theta_new
    }
    
    if (m %in% log_m) {
      theta_samples <- cbind(theta_samples, t(theta_now))
      proposals <- cbind(proposals, curr_proposal)
      X_traj_samples <- rbind(X_traj_samples, X_traj_now)
      if (m %in% save_m) {  # save samples
        save(X_traj_samples, theta_samples, file = "MCMC_temp_out")
      }
    }
    toc()
  }
  
  MCMC_out <- list()
  MCMC_out$accept_rate <- accept / m  # acceptance rate
  MCMC_out$p_samples <- p_samples
  MCMC_out$proposals <- proposals
  MCMC_out$X_samples <- X_traj_samples
  
  list(theta_samples = theta_samples, MCMC_out = MCMC_out)
}