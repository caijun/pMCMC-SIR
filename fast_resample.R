fast_resample <- function(w, MCMC_params) {
  w_norm <- w / sum(w)
  r <- rmultinom(1, size = MCMC_params$J_particles, prob = w_norm)  # resample by multinomial sampling with replacement
  
  k <- rep(0, MCMC_params$J_particles)
  end_index <- 0
  for (j in 1:MCMC_params$J_particles) {
    reps <- r[j]
    if (reps > 0) {
      curr_index <- end_index + 1
      end_index <- curr_index + reps - 1
      k[curr_index:end_index] <- j
    }
  }
  
  k_shuffle <- sample(MCMC_params$J_particles)
  k <- k[k_shuffle]
}