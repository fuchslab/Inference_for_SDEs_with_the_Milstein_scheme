# This R script is used to perform HMC+NUTS sampling using the package rstan
# for the ODE assuming multiplicative measurement error and that the initial time
# point t0 of mRNA realease needs to be estimated.
suppressMessages(library(rstan))

# parameters for the estimation procedure --------------------------------------------
# get arguments handed-over in command line
inputArgs <- commandArgs(trailingOnly = TRUE)

if (length(inputArgs)){
  # folder containing the file with obvervations
  obsFolder <- inputArgs[1]
  # number of observations to be used
  M <- as.numeric(inputArgs[2]) #
  # number of iterations
  numIterations <- as.numeric(inputArgs[3])
  # index of the path
  path_index <- as.numeric(inputArgs[4])

}else{
  # folder containing the file with obvervations
  obsFolder <- "CIR_alpha_1_beta_1_sigma_0.25_x0_3" #  GBM_alpha_1_sigma_2_x0_100 CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10
  # number of observations to be used
  M <- 10 # , 10, 20 ,50
  # number of iterations
  numIterations <- 5e5
  # index of the path
  path_index <- 2
}


# # load data ------------------------------------------------------------------------
# load observations
model_type <- substr(obsFolder,1,3)
load(paste('simulation_study/', obsFolder,'/', model_type, '_obs.data', sep = "")) # contains Y_obs, tau, true_theta, estim_seed
# length of theta
len_theta <- length(true_theta)

# select M observations
index <- seq(from = 1, to = nrow(Y_obs), by = (nrow(Y_obs)-1)/M)
Y_obs <- Y_obs[index, path_index]
tau <- tau[index]
time_step <- unique(round(diff(tau), digits = 10))



# HMC + NUTS sampling using rstan ----------------------------------------------------
rstan_options(auto_write = TRUE)
n_chains <- 4
options(mc.cores = n_chains)

path_stanfit_object <- paste0("simulation_study/",
                               obsFolder , "/output/Stan/stanfit_objects/stanfit_object_M_",
                              M, "_path_", path_index, ".rds", sep = "")
path_stan_optimizing_result <- paste0("simulation_study/",
                              obsFolder , "/output/Stan/optimizing_results/stan_optim_res_M_",
                              M, "_path_", path_index, ".rds", sep = "")

stan_model_file <- paste("simulation_study/Stan/", model_type,
                         "_model.stan", sep = "")

if(model_type == "GBM" ){
  sdata <- list(
    M = M,
    y_obs = Y_obs,
    time_step = time_step
  )

  stan_model_object <- stan_model(file = stan_model_file)
}else if(model_type == "CIR" ){
  sdata <- list(
    M = M,
    y_obs = Y_obs,
    time_step = time_step,
    alpha = true_theta[1]
  )

  stan_model_object <-
    stan_model(file = stan_model_file, model_name = "CIR_model",
               allow_undefined = TRUE,
               includes = paste0('\n#include "', getwd(),
                                 '/simulation_study/Stan/interface.hpp"\n'))
}else{
  stop("In valid model type!")
}

path_sample_file <- paste0("simulation_study/",
                              obsFolder , "/output/Stan/stanfit_objects/stan_sample_file_M_",
                              M, "_path_", path_index, ".rds", sep = "")

# sampling from the posterior
stanfit_object <- sampling(object = stan_model_object,
                           data = sdata,
                           seed = estim_seeds[path_index],
                           chains = n_chains,
                           iter = numIterations)#,
                           #sample_file = path_sample_file)
saveRDS(stanfit_object, file = path_stanfit_object)

# optimizing the posterior
stan_optimizing_result <- rstan::optimizing(stan_model_object, data = sdata)
saveRDS(stan_optimizing_result, file = path_stan_optimizing_result)
