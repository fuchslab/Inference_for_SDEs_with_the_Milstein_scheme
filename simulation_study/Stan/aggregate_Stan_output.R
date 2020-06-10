
suppressMessages(library(boa, quietly = TRUE)) # for highest probability density interval
suppressMessages(library(MCMCpack, quietly = TRUE)) # for the inverse gamma distribution
library(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
suppressMessages(library(mcmcse))
library(stringr)
suppressMessages(library(rstan))

# excute file within container with e.g.
# Rscript --vanilla simulation_study/aggregate_output.R <obsFolder> <M> <m>

# parameters for the estimation procedure --------------------------------------------
# get arguments handed-over in command line
inputArgs <- commandArgs(TRUE)

if (length(inputArgs)){
  # folder containing the file with obvervations
  obsFolder <- inputArgs[1]
  M <- as.numeric(inputArgs[2])
}else{
  # folder containing the file with obvervations
  obsFolder <- "CIR_alpha_1_beta_1_sigma_0.25_x0_3" #GBM_alpha_1_sigma_2_x0_100   CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10  
  M <- 50
}

folderName <- paste('simulation_study/', obsFolder, "/aggregated_output",
                    "/", sep = "")
if (!dir.exists(folderName)){
  dir.create(folderName)
}
setwd(dir = paste("simulation_study/", obsFolder, sep = ""))

# number of observation paths
N <- 100


model_type <- substr(obsFolder,1,3)
if(model_type == "GBM"){
  par_names <- c("alpha", "sigma2")
}else if(model_type == "CIR"){
  par_names <- c("beta", "sigma2")
}

cnames_separate <- c("mean", "median", "variance", "hpd_low", "hpd_up", "ESS_coda", "n_eff", "Rhat", "opt_result")
ncol_separate <- length(cnames_separate)
cnames_overall <- c("multivarESS", "maxDuration", "covariance")
ncol_overall <- length(cnames_overall)

results_separate <- array(rep(NA, length = N * ncol_separate * length(par_names)), 
                           dim = c(N, ncol_separate, length(par_names)),
                          dimnames = list(path_index = 1:N, results=cnames_separate , 
                                          parameter=par_names))

resultsOverall <- matrix(rep(NA, length = N * ncol_overall), ncol = ncol_overall)
colnames(resultsOverall) <- cnames_overall

for(k in 1:N){
  file_name <- paste("stanfit_objects/stanfit_object_", "M_", M, "_path_", k,
                     ".rds", sep = "")
  file_name_optim <- paste("optimizing_results/stan_optim_res_", "M_", M, 
                           "_path_", k,  ".rds", sep = "")

  print(file_name)
  
  if(length(file_name)){
    stanfit_object <- try(readRDS(paste("output/Stan/", file_name, sep = "")))
    stan_optim_res <- try(readRDS(paste("output/Stan/", file_name_optim, sep = "")))
    
    if(class(stanfit_object) != "try-error"){
      sample <- as.matrix(stanfit_object)
      stan_summary <- summary(stanfit_object)$summary
      
      for(i in 1:length(par_names)){
        results_separate[k, "mean", i] <- mean(sample[ ,i])
        results_separate[k, "median", i] <- median(sample[ ,i])
        results_separate[k, "variance", i] <- var(sample[ ,i])
        hpd <- boa.hpd(sample[ ,i], alpha = 0.05)
        results_separate[k, "hpd_low", i] <- hpd[1]
        results_separate[k, "hpd_up", i] <- hpd[2]
        results_separate[k, "ESS_coda", i] <- coda::effectiveSize(sample[ ,i])
        results_separate[k, "n_eff", i] <- stan_summary[i, "n_eff"]
        results_separate[k, "Rhat", i] <- stan_summary[i, "Rhat"]
        results_separate[k, "opt_result", i] <- stan_optim_res$par[i]
      }
      
      resultsOverall[k, "multivarESS"] <- mcmcse::multiESS(sample[,1:2])
      resultsOverall[k, "maxDuration"] <- max(rowSums(get_elapsed_time(stanfit_object)))
      resultsOverall[k, "covariance"] <- cov(sample[,1],sample[ ,2])

      remove(stanfit_object)
    }
  }
}

fileName2 <- paste("true_posterior","M", M, sep = "_")

folderName2 <- "aggregated_output/"

save(results_separate, resultsOverall, M,
     file = paste(folderName2, fileName2, '.data', sep=''))


