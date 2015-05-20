function_f_SIR <- function(MCMC_params, X_F, t_n, t_n_plus_1, theta) {
  t <- seq(t_n, t_n_plus_1, by = MCMC_params$dt)
  time_increments <- length(t)
  
  # Copy over parameters
  dt <- MCMC_params$dt
  mu <- theta[1]
  gamma <- theta[2]
  R0_avg <- theta[3]
  alpha <- theta[4]
  F_noise <- theta[5]
  beta <- R0_avg * (gamma + mu)
  
  # Initial conditions for state variables
  S <- I <- P <- R <- C <- matrix(0, nr = time_increments, nc = MCMC_params$J_particles)
  S[1, ] <- X_F[1, ]
  I[1, ] <- X_F[2, ]
  P[1, ] <- X_F[3, ]
  R[1, ] <- P[1, ] - S[1, ] - I[1, ]
  C[1, ] <- 0
  
  eta_rand <- matrix(rnorm((time_increments - 1) * MCMC_params$J_particles), 
                     nr = time_increments - 1, nc = MCMC_params$J_particles)
  
  for (i in 1:(time_increments - 1)) {  # for each dt time increment
    seas_now <- t[i] / 12 - floor(t[i] / 12)
    beta_now <- beta * (1 + alpha * cos(2 * pi * seas_now))
    
    f_term <- F_noise * beta_now * S[i, ] / P[i, ] * I[i, ]
    dN_SI <- beta_now * S[i, ] / P[i, ] * I[i, ] * dt + f_term * eta_rand[i, ] * dt
    dN_SI[dN_SI < 0] <- 0  # make sure this doesn't go negative!
    dN_IR <- gamma * I[i, ] * dt
    dN_SD <- mu * S[i, ] * dt
    dN_ID <- mu * I[i, ] * dt
    dN_RD <- mu * R[i, ] * dt
    dN_BS <- mu * P[i, ] * dt
    
    dS <- dN_BS - dN_SI + dN_SD
    dI <- dN_SI - dN_IR - dN_ID
    dR <- dN_IR - dN_RD
    
    S[i + 1, ] <- S[i, ] + dS
    I[i + 1, ] <- I[i, ] + dI
    R[i + 1, ] <- R[i, ] + dR
    
    P[i + 1, ] <- S[i + 1, ] + I[i + 1, ] + R[i + 1, ]  # THIS NEEDS TO BE CHANGED!
    C[i + 1, ] <- C[i, ] + dN_SI
  }
  
  end <- time_increments
  X_traj <- t(I)
  S_traj <- t(S)
  P_traj <- t(P)
  X_P <- rbind(S[end, ], I[end, ], P[end, ])
  C_value <- C[end, ]
  list(X_P = X_P, C_value = C_value, X_traj = X_traj, S_traj = S_traj, P_traj = P_traj)
}