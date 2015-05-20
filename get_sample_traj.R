get_sample_traj <- function(k, X_traj_final, MCMC_params, X_traj_indexes) {
  random_lineage <- ceil(runif(1) * MCMC_params$J_particles)
  n_times <- dim(k)[2]
  path <- rep(0, n_times)
  path[n_times] <- random_lineage
  for (j in seq(n_times, 2, by = -1)) {
    parent_index <- path[j]
    path[j - 1] <- k[parent_index, j - 1]
  }
  
  sub_n <- 1 / MCMC_params$dt
  X_traj_sample <- X_traj_final[1, 1]
  for (k in 1:length(path)) {
    start_index <- (k - 1) * sub_n + 2
    end_index <- k * sub_n + 1
    X_traj_new <- X_traj_final[path[k], start_index:end_index]
    X_traj_sample <- c(X_traj_sample, X_traj_new)
  }
  
  X_traj_sample <- X_traj_sample[X_traj_indexes]
}