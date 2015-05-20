### convert matlab mat file to RData file ###
setwd("~/Desktop/Rasmussen.etal-PLoSComputBiol-2011")

library(TTmisc)
# read matlab mat files
library(R.matlab)
list[] <- readMat("pMCMCSourceCode/SIR_endemic_sim3_mockData.mat", fixNames = F)
save("mockdata", "epi_params", "prevalence", "t_data", 
     file = "pMCMCsrc/SIR_endemic_sim3_mockData.RData")
list[] <- readMat("pMCMCSourceCode/SIR_endemic_sim3_coalTimes.mat", fixNames = F)
save("sample_sizes", "coal_times_back", "sample_times_back", 
     file = "pMCMCsrc/SIR_endemic_sim3_coalTimes.RData")
list[] <- readMat("pMCMCSourceCode/SIR_covmat_041111.mat", fixNames = F)
save("cov_mat", file = "pMCMCsrc/SIR_covmat_041111.RData")
### convert matlab mat file to RData file ###

setwd("~/Desktop/Rasmussen.etal-PLoSComputBiol-2011/pMCMCsrc")
source("get_G_events.R")
source("function_f_SIR.R")
source("pdTreeMatrix.R")
source("fast_resample.R")
source("run_SMC.R")
source("get_sample_traj.R")
source("get_Likelihoods.R")
source("run_MCMC.R")
source("get_traces.R")
library(pracma)  # for repmat, tic, toc functions
library(mvtnorm)  # for rmvnorm function

# Input files
mockData_file <- 'SIR_endemic_sim3_mockData.RData'  # File containing mock time series and true parameters
coalTimes_file <- 'SIR_endemic_sim3_coalTimes.RData'  # File containing genealogical data (sampling and coalescence times)
output_file <- 'SIR_endemic_sim3_MCMCout.RData'  # MCMC output goes here

# Load data files
load(mockData_file)  # time series data
load(coalTimes_file)  # genealogical data: coal times and sampling times in months since present
load("SIR_covmat_041111.RData")  # covariance matrix for MCMC proposal density

# Get time series data
data <- list()
data$t_vals <- t_data - min(t_data)  # observation times (months since t0)
data$y_vals <- mockdata  # case count data (monthly incidence)
data$P <- epi_params$N_init  # population size (constant)

# MCMC set up
MCMC_params <- list()
MCMC_params$J_particles <- 200  # number of particles to use in SMC
MCMC_params$iterations <- 100  # MCMC steps
MCMC_params$log_steps <- 10  # log MCMC samples every x steps
MCMC_params$save_steps <- 1000  # save MCMC samples every x steps
MCMC_params$dt <- epi_params$dt  # integration time step
MCMC_params$burn_in <- 100  # not important

# Get genealogical data
coal_times <- coal_times_back  # coal times in months since present (tEnd)!!
sample_times <- sample_times_back  # sample times in months since present (tEnd)!!
obsv_times <- max(t_data) - t_data  # observation times in months since present (tEnd)!!
obsv_times <- round(obsv_times)
# dataG is the data structure that contains the information used to 
# calculate the likelihood of the genealogy over preset time intervals
list[G_lineages, G_coal_event, G_lineages_dt, G_indices, G_dt_ref, event_times] <- get_G_events(MCMC_params, sample_times, sample_sizes, coal_times, obsv_times)
dataG <- list()
dataG$coal_events <- G_coal_event  # vector that indicates if event is a coalescence events
dataG$lineages <- G_lineages  # vector for number of lineages over time
dataG$indices <- G_indices  # starting and ending index for each observation interval
dataG$dt_ref <- G_dt_ref  # reference to nearest dt integration time
dataG$event_time <- event_times  # all events (sampling, coalescence and observation events)

# View lineags over time plot aligned with time series
par(mfrow = c(1, 2))
plot(G_lineages_dt)
plot(prevalence)
par(mfrow = c(1, 1))

# Initial conditions (can be inferred)
X_I <- rbind(epi_params$S_init, epi_params$I_init, epi_params$N_init)

# Initial parameter values
# mu - host birth/death rate
# gamma - rate of recovery
# R0_avg
# alpha - seasonality
# F_noise - environmental noise
# rho - reporting rate
# tau - observation variance
theta_now <- c(epi_params$mu, epi_params$gamma, epi_params$R0, epi_params$alpha, 
               epi_params$F_noise, epi_params$rho, epi_params$tau)

# Covariances for proposal densities
# theta_var <- c(0, 0.03, 0.04, 0.0003, 0.0004, 0.0001, 0.05)
MCMC_params$theta_cov <- cov_mat  # use diag(theta_var) if no starting covariance matrix

# Run MCMC
list[theta_samples, MCMC_out] <- run_MCMC(X_I, data, dataG, MCMC_params, theta_now)

save(theta_samples, MCMC_out, MCMC_params, data, dataG, file = output_file)