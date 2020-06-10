
suppressMessages(library(boa, quietly = TRUE)) # for highest probability density interval
suppressMessages(library(MCMCpack, quietly = TRUE)) # for the inverse gamma distribution
library(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
suppressMessages(library(mcmcse))
library(stringr)

# excute file within container with e.g.
# Rscript --vanilla simulation_study/aggregate_output.R <obsFolder> <M> <m>

# parameters for the estimation procedure --------------------------------------------
# get arguments handed-over in command line
inputArgs <- commandArgs(TRUE)

if (length(inputArgs)){
  # folder containing the file with obvervations
  obsFolder <- inputArgs[1]
  M <- as.numeric(inputArgs[2])
  m <- as.numeric(inputArgs[3])
  methodPathUpdate <- inputArgs[4]
  approxTransDens <-  inputArgs[5]
  approxPropDens <- inputArgs[6]
}else{
  # folder containing the file with obvervations
  obsFolder <- "CIR_alpha_1_beta_1_sigma_0.25_x0_3" #"GBM_alpha_1_sigma_2_x0_100"  CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10
  M <- 20
  m <- 1
  methodPathUpdate <- "leftConditioned" #"DBMilstein"
  approxTransDens <-  "Milstein"
  approxPropDens <- "Milstein"
}

folderName <- paste('simulation_study/', obsFolder, "/aggregated_output",
                    "/", sep = "")
if (!dir.exists(folderName)){
  dir.create(folderName)
}
setwd(dir = paste("simulation_study/", obsFolder, sep = ""))

file_names <- list.files(path = "output/")

# number of observation paths
N <- 100
# number of iterations used as burn-in
burnIn <- 5000

model_type <- substr(obsFolder,1,3)
if(model_type == "GBM"){
  par_names <- c("alpha", "sigma2")
  sample_indices <- c(1,2)
}else if(model_type == "CIR"){
  par_names <- c("beta", "sigma2")
  sample_indices <- c(2,3)
}

cnames_separate <- c("mean", "median", "variance", "hpd_low", "hpd_up", "ESS_coda", "mode")
ncol_separate <- length(cnames_separate)
cnames_overall <- c("numIterations", "multivarESS", "Duration", "ARpath",
                    "ARparam", "covariance", "nNegPointProposals", "nMBSwitchToEuler")
ncol_overall <- length(cnames_overall)

results_separate <- array(rep(NA, length = N * ncol_separate * length(par_names)),
                          dim = c(N, ncol_separate, length(par_names)),
                          dimnames = list(path_index = 1:N, results=cnames_separate ,
                                          parameter=par_names))

resultsOverall <- matrix(rep(NA, length = N * ncol_overall), ncol = ncol_overall)
colnames(resultsOverall) <- cnames_overall


for(k in 1:N){
    fileName_fragment <-
      paste("M_", M, "_m_", m, "_path_", methodPathUpdate, "_td_",
            approxTransDens, "_pd_", approxPropDens, "_", k,
            ".data", sep = "")

  file_name <- file_names %>% str_subset(fileName_fragment)
  print(file_name)

  if(length(file_name)){
    a <- try(load(paste("output/", file_name, sep = "")))

    if(class(a)!= "try-error"){
      attach(results, warn.conflicts = FALSE)

      numIter <- dim(theta)[2]
      sample <- t(theta[sample_indices , (burnIn+1):numIter])

      for(i in 1:length(par_names)){
        results_separate[k, "mean", i] <- mean(sample[ ,i])
        results_separate[k, "median", i] <- median(sample[ ,i])
        results_separate[k, "variance", i] <- var(sample[ ,i])
        hpd <- boa.hpd(sample[ ,i], alpha = 0.05)
        results_separate[k, "hpd_low", i] <- hpd[1]
        results_separate[k, "hpd_up", i] <- hpd[2]
        results_separate[k, "ESS_coda", i] <- coda::effectiveSize(sample[ ,i])
        dens <- density(sample[ ,i])
        results_separate[k, "mode", i] <- dens$x[which.max(dens$y)]
      }

      resultsOverall[k, "multivarESS"] <- mcmcse::multiESS(sample[,1:2])
      resultsOverall[k, "Duration"] <- duration[[1]]
      resultsOverall[k, "covariance"] <- cov(sample[,1],sample[ ,2])
      resultsOverall[k, "nNegPointProposals"] <- numNegPointProposals
      resultsOverall[k, "nMBSwitchToEuler"] <- numMBSwitchToEuler
      resultsOverall[k, "numIterations"] <- numIter
      resultsOverall[k, "ARpath"] <- acceptanceRatePath
      resultsOverall[k, "ARparam"] <- acceptanceRateParam

      remove(results)
    }
  }
}

num_missing <- N - sum(!is.na(results_separate[, "mean", 1]))

if(m != 1){
  fileName2 <- paste("appr_posterior","M", M, "m", m, methodPathUpdate,
                     "td", approxTransDens, "pd", approxPropDens, sep = "_")
}else{
  fileName2 <- paste("appr_posterior","M", M, "m", m, "td", approxTransDens,
                     sep = "_")
}

folderName2 <- "aggregated_output/"
true_theta <- inputs$true_theta

save(results_separate, resultsOverall,
     M, m, N, burnIn, true_theta, num_missing,
     file = paste(folderName2, fileName2, '.data', sep=''))


