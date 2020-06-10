
suppressMessages(library(coda, quietly = TRUE))
suppressMessages(library(MASS, quietly = TRUE))
suppressMessages(library(MCMCpack, quietly = TRUE)) # for the inverse gamma distribution
suppressMessages(library(boa, quietly = TRUE)) # for highest probability density interval



# parameters for the estimation procedure --------------------------------------------
# get arguments handed-over in command line
inputArgs <- commandArgs(TRUE)

if (length(inputArgs)){
  # folder containing the file with obvervations
  obsFolder <- inputArgs[1]
  # number of observations to be used
  M <- as.numeric(inputArgs[2]) # 25, 50, 100, 1000
  # number of imputed points in each inter-observation interval
  m <-  as.numeric(inputArgs[3])  # 2, 10
  # number of iterations
  numIterations <- as.numeric(inputArgs[4])

  methodPathUpdate <- inputArgs[5]  # "leftConditioned", "MB", "DBMilstein
  approxTransDens <-  inputArgs[6] # "Euler",  "Milstein"
  approxPropDens <-  inputArgs[7] # "Euler",  "Milstein"

  # # hyperparameters for the priori distributions of the parameters
  # # use conjugate priors
  # # -> alpha ~ N(alpha_0,rho_2) with
  # alpha_0 <- as.numeric(inputArgs[8])
  # rho_2 <- as.numeric(inputArgs[9])
  # # sigma_2 ~ IG(kappa_0,nu_0) inverse Gamma distribution with
  # kappa_0 <- as.numeric(inputArgs[10])
  # nu_0 <- as.numeric(inputArgs[11])
  #
  # # hyperparameters for the random walk parameter update
  # gamma_alpha <- as.numeric(inputArgs[12])
  # gamma_sigma <- as.numeric(inputArgs[13])

  # index of the path
  path_index <- as.numeric(inputArgs[8])

  # max. computational time in seconds
  comp_time <- as.numeric(inputArgs[9])

}else{
  # folder containing the file with obvervations
  obsFolder <- "CIR_alpha_1_beta_1_sigma_0.25_x0_3" #"GBM_alpha_1_sigma_2_x0_100"  CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10
  # number of observations to be used
  M <- 10 # , 10, 20 ,50
  # number of imputed points in each inter-observation interval
  m <- 5 # , 1, 2, 5
  # number of iterations
  numIterations <- 1e5

  methodPathUpdate <- "MB" # "leftConditioned", "MB", "DBMilstein"
  approxTransDens <-  "Milstein" #  "Euler", "Milstein"
  approxPropDens <- "Milstein" # "Euler", "Milstein"

  # # hyperparameters for the priori distributions of the parameters
  # # use conjugate priors
  # # -> alpha ~ N(alpha_0,rho_2) with
  # alpha_0 <- 0
  # rho_2 <- 1
  # # sigma_2 ~ IG(kappa_0,nu_0) inverse Gamma distribution with
  # kappa_0 <- 2
  # nu_0 <- 2
  #
  # # hyperparameters for the random walk parameter update
  # gamma_alpha <- 0.5
  # gamma_sigma <- 0.5

  # index of the path
  path_index <- 2

  # max. computational time in seconds
  comp_time <- 5
}


# ------------------------------------------------------------------------------------
# for Alg. 7.3. for the choice of the update intervals
#lambda <- 5

# burn-in rate
burnIn <- 0.1

# ------------------------------------------------------------------------------------
# load observations
load(paste('simulation_study/', obsFolder,'/', substr(obsFolder,1,3), '_obs.data', sep = "")) # contains Y_obs, tau, true_theta, estim_seed
# length of theta
len_theta <- length(true_theta)

# select M observations
index <- seq(from = 1, to = nrow(Y_obs), by = (nrow(Y_obs)-1)/M)
Y_obs <- Y_obs[index, path_index]
tau <- tau[index]


# load functions
source(paste("main_functions/", substr(obsFolder,1,3), "_problem_specific_parameter_and_functions.R", sep = ""))
source("main_functions/functions_for_parameter_estimation.R")
source("main_functions/parameter_estimation.R")

# initialize the random number generator
set.seed(estim_seeds[path_index])

# ------------------------------------------------------------------------------------
# gather inputs
inputs <- list(Y_obs = Y_obs, tau = tau, true_theta = true_theta,
               thetaNames = thetaNames, indizes_estim_param = indizes_estim_param,
               methodPathUpdate = methodPathUpdate,
               methodParamUpdate = "RandomWalk",
               approxTransDens = approxTransDens,
               approxPropDens = approxPropDens,
               lambda = lambda,
               numIterations = numIterations, M = M, m = m, seed = estim_seeds[path_index],
               hyperparam = hyperparam, burnIn = burnIn)
# ------------------------------------------------------------------------------------
# generate folder and file names and make sure that the folder exits
folderName <- paste('simulation_study/', obsFolder, "/output",
                    "/", sep = "")
if (!dir.exists(folderName)){
  dir.create(folderName)
}
str_hyperparam <- c()
for(i in 1:length(hyperparam)){
  str_hyperparam <- paste(str_hyperparam, names(hyperparam)[i], hyperparam[[i]], sep = "_")
}
fileName <- paste(substr(obsFolder,1,3), substr(str_hyperparam, 2, nchar(str_hyperparam)),
                  "M", M, "m", m, "path", methodPathUpdate, "td",
                  approxTransDens, "pd", approxPropDens, path_index, sep = "_")
# ------------------------------------------------------------------------------------
# do the parameter estimation
tryCatch(duration <-
           system.time((results <- estimate_parameters(methodPathUpdate = methodPathUpdate,
                                                       methodParamUpdate = "RandomWalk",
                                               approxTransDens = approxTransDens,
                                               approxPropDens = approxPropDens,
                                               numIterations = numIterations,
                                               m = m,
                                               Y_obs, tau,
                                               comp_time = comp_time))),
         error = function(e){
           error_message <- conditionMessage(e)
           traceback_output <- traceback()
           save(inputs, error_message, traceback_output,
                file = paste(folderName, "error_", fileName, ".data", sep = ""))
         })


if(exists("duration")){

  # ------------------------------------------------------------------------------------
  #  save data
  # results list already contains: Y_imp, tpoints, theta, acceptanceRatePath, acceptanceRateParam
  results$duration <- duration
  save(results, inputs, file = paste(folderName, fileName, ".data", sep = ""))

  # # ------------------------------------------------------------------------------------
  # # generate results for some plots
  # cnames <- vector(mode="character", length = 0)
  # for(i in indizes_estim_param){
  #   varName <- thetaNames[i]
  #   cnames <- c(cnames, paste("mean_", thetaNames[i], sep = ""),
  #               paste("hpd_", thetaNames[i], "_l", sep = ""),
  #               paste("hpd_", thetaNames[i], "_u", sep = ""))
  # }
  # # calculate the point estimates for the components of theta and their
  # # highest probability density intervals
  # resultTable1 <- vector(mode = 'numeric')
  # for (i in indizes_estim_param){
  #   resultTable1 <- c(resultTable1,
  #                     mean(results$theta[i, ][floor(burnIn * length(results$theta[i, ])):length(results$theta[i, ])]),
  #                     boa.hpd(results$theta[i, ][floor(burnIn * length(results$theta[i, ])):length(results$theta[i, ])], 0.05))
  # }
  # resultTable1 <- matrix(round(resultTable1,2), nrow = 1, ncol = length(cnames))
  # colnames(resultTable1) <- cnames
  # rownames(resultTable1) <- ""
  #
  #
  # resultTable2 <- c(round(results$acceptanceRatePath,3),
  #                   round(results$acceptanceRateParam,3),
  #                   round(duration[3],4),
  #                   results$numNegPointProposals,
  #                   results$numMBSwitchToEuler)
  # cnames2 <- c("acceptRatePath", "acceptRateParam", "duration", "# of neg. point proposals", "# of switches to MBEuler")
  # resultTable2 <- matrix(resultTable2, nrow = 1, ncol = length(cnames2))
  # colnames(resultTable2) <- cnames2
  # rownames(resultTable2) <- ""
  #
  #
  # # plots
  # pdf(paste(folderName, fileName, '.pdf', sep=''), paper = 'a4', width = 0, height = 0)
  # layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 8, 1), heights = c(1.5,4.4,1.8,1.7,1.7,3.5,3.5,3))
  # par(mar = c(1.5,2,2.5,2), oma = c(1, 1, 1, 1))
  # # paramters
  #
  # tit <- c()
  # for(i in 1:length(true_theta)){
  #   tit <- paste(tit, thetaNames[i], " = ", true_theta[i], ", ", sep = "")
  # }
  #
  # gplots::textplot(paste(tit,
  #                        "M = ", M,", m = ", m,  ", \npath = ", path_index, ", seed = ", estim_seeds[path_index],
  #                        sep = ''))
  # # path
  # if(m > 1){
  #   plot(tau, Y_obs, pch = 18, ylim = c(min(results$Y_imp[, ncol(results$Y_imp)]),
  #                                       max(results$Y_imp[, ncol(results$Y_imp)])))
  #   lines(results$tpoints, results$Y_imp[, ncol(results$Y_imp)], col = 2)
  #   par(xpd = TRUE)
  #   legend(x = 'topright', inset = c(0,-.1), lty = c(NA,1), pch = c(18,NA), lwd=c(1,1),
  #          col=c(1,2), bty='n' , horiz = TRUE,
  #          legend = c('observations', 'last imputed path Y_imp'))
  # }else{
  #   plot(tau, Y_obs, pch = 18, type = 'b')
  #   par(xpd = TRUE)
  #   legend(x = 'topright', inset = c(0,-.1), lty = c(1), pch = c(18), lwd=c(1),
  #          col=c(1,2), bty='n' , horiz = TRUE,
  #          legend = c('observations'))
  # }
  #
  #
  # par(xpd = FALSE)
  #
  #
  # str_hyperparam <- c()
  # for(i in 1:length(hyperparam)){
  #   str_hyperparam <- paste(str_hyperparam, names(hyperparam)[i], " = ", hyperparam[[i]], ", ", sep = "")
  # }
  #
  # gplots::textplot(paste("methodPathUpdate = ", methodPathUpdate,
  #                        ", methodParamUpdate = RandomWalk",
  #                        ",\napproxTransDens = ", approxTransDens,
  #                        ", approxPropDens = ", approxPropDens,
  #                        "\n", str_hyperparam,
  #                        sep = ''))
  # gplots::textplot(resultTable1)
  #
  # gplots::textplot(resultTable2)
  #
  # # alpha
  # plot(results$theta[indizes_estim_param[1], ], type = "l")
  # title(paste("MCMC ", thetaNames[indizes_estim_param[1]]))
  # abline(true_theta[indizes_estim_param[1]], 0, col = 2)
  # abline(mean(results$theta[indizes_estim_param[1], floor(burnIn * ncol(results$theta)):ncol(results$theta)]),
  #        0, col = 3)
  # abline(boa.hpd(results$theta[indizes_estim_param[1], ][floor(burnIn * length(results$theta[indizes_estim_param[1], ])):length(results$theta[indizes_estim_param[1], ])], 0.05)[1], 0, col = 3, lty = 2)
  # abline(boa.hpd(results$theta[indizes_estim_param[1], ][floor(burnIn * length(results$theta[indizes_estim_param[1], ])):length(results$theta[indizes_estim_param[1], ])], 0.05)[2], 0, col = 3, lty = 2)
  # abline(v = floor(burnIn * length(results$theta[indizes_estim_param[1], ])), col = 4)
  # # sigma_2
  # plot(results$theta[indizes_estim_param[2], ], type = "l")
  # title(paste("MCMC ", thetaNames[indizes_estim_param[2]]))
  # abline(true_theta[indizes_estim_param[2]], 0, col = 2)
  # abline(mean(results$theta[indizes_estim_param[2], floor(burnIn * ncol(results$theta)):ncol(results$theta)]),
  #        0, col = 3)
  # abline(boa.hpd(results$theta[indizes_estim_param[2], ][floor(burnIn * length(results$theta[indizes_estim_param[2], ])):length(results$theta[indizes_estim_param[2], ])], 0.05)[1], 0, col = 3, lty = 2)
  # abline(boa.hpd(results$theta[indizes_estim_param[2], ][floor(burnIn * length(results$theta[indizes_estim_param[2], ])):length(results$theta[indizes_estim_param[2], ])], 0.05)[2], 0, col = 3, lty = 2)
  # abline(v = floor(burnIn * length(results$theta[indizes_estim_param[2], ])), col = 4)
  # # log-posterior densities
  # plot(results$log_posterior, type = 'l')
  # title("log-posterior density values")
  #
  # r <- dev.off()
}

print(results$duration)
print(paste("Number of finished iterations: ", results$numInterFinished, sep=""))
