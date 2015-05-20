get_traces <- function(theta_samples, x) {
  # Make trace plots for all parameters in theta
  theta_new <- theta_samples[x, , drop = F]
  
  nr <- nrow(theta_new)
  par(mfrow = c(nr, 1))
  for (i in 1:nr) {
    plot(theta_new[i, ], type = "l")
  }
  par(mfrow = c(1, 1))
}