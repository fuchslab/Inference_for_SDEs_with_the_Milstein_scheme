
suppressMessages(library(boa, quietly = TRUE)) # for highest probability density interval
suppressMessages(library(MCMCpack, quietly = TRUE)) # for the inverse gamma distribution
library(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(mcmcse)

v_alpha <- 1#c(1, 5)
v_sigma_2 <- 2#c( 0.25, 2)


# number of observations to be used
v_M <- c(50)
# number of imputed points in each inter-observation interval
v_m <- c(1,2,5)
# number of observation paths
N <- 100
# burn-in rate
burnIn <- 0.1
# number of iterations
numIterations <- 1e5
# acceptance rates are only valid for the corresponding max. numIterations (1e5 or 1e6)
include_acc_rate <- TRUE

methodPathUpdate <- c("leftConditioned", "MB")
methodPathUpdate_short <- c("leftCondi", "MB")
approxTransDens <- c("Euler", "Milstein")
approxPropDens <- c("Euler", "Milstein")

# hyperparameters for the priori distributions of the parameters
# use conjugate priors
# -> alpha ~ N(alpha_0,rho_2) with
alpha_0 <- 0
rho_2 <- c(10)
# sigma_2 ~ IG(kappa_0,nu_0) inverse Gamma distribution with
kappa_0 <- c(2)
nu_0 <- c(2)

# hyperparameters for the random walk parameter update
gamma_alpha <- 0.5
gamma_sigma <- 0.5

for (a in 1:length(v_alpha)){
  for (b in 1:length(v_sigma_2)){

    # theta <- (alpha, sigma^2)
    alpha <- v_alpha[a]
    sigma_2 <- v_sigma_2[b]
    true_theta <- c(alpha, sigma_2)

    folderName <- paste("GBM_alpha_", alpha , "_sigma_", sigma_2,
                        "/", sep = "")
    setwd(dir = paste("simulation_study/", folderName, sep = ""))


    for (i in 1:length(v_M) ){
      for (j in 1:length(v_m)){
        v <- length(approxPropDens) * length(approxTransDens) * length(methodPathUpdate)
        meanAlpha <- matrix(rep(NA, length = N*v), ncol = v)
        meanSigma <- matrix(rep(NA, length = N*v), ncol = v)
        modeAlpha <- matrix(rep(NA, length = N*v), ncol = v)
        modeSigma <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_l_Alpha <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_l_Sigma <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_u_Alpha <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_u_Sigma <- matrix(rep(NA, length = N*v), ncol = v)
        effssAlpha <- matrix(rep(NA, length = N*v), ncol = v) # effective sample size
        effssSigma <- matrix(rep(NA, length = N*v), ncol = v)
        multivarESS <- matrix(rep(NA, length = N*v), ncol = v)
        colnames(meanAlpha)<-colnames(meanAlpha, do.NULL = FALSE)
        colnames(meanSigma)<-colnames(meanSigma, do.NULL = FALSE)
        colnames(modeAlpha)<-colnames(modeAlpha, do.NULL = FALSE)
        colnames(modeSigma)<-colnames(modeSigma, do.NULL = FALSE)
        colnames(hpd_l_Alpha)<-colnames(hpd_l_Alpha, do.NULL = FALSE)
        colnames(hpd_l_Sigma)<-colnames(hpd_l_Sigma, do.NULL = FALSE)
        colnames(hpd_u_Alpha)<-colnames(hpd_u_Alpha, do.NULL = FALSE)
        colnames(hpd_u_Sigma)<-colnames(hpd_u_Sigma, do.NULL = FALSE)
        colnames(effssAlpha)<-colnames(effssAlpha, do.NULL = FALSE)
        colnames(effssSigma)<-colnames(effssSigma, do.NULL = FALSE)
        MAP_alpha <- matrix(rep(NA, length = N), ncol = 1)
        MAP_sigma <- matrix(rep(NA, length = N), ncol = 1)
        ML_alpha <- matrix(rep(NA, length = N), ncol = 1)
        ML_sigma <- matrix(rep(NA, length = N), ncol = 1)
        nNegPointProposals <- matrix(rep(NA, length = N*v), ncol = v)
        nMBSwitchToEuler <- matrix(rep(NA, length = N*v), ncol = v)

        if(include_acc_rate){
          ARpath <- matrix(rep(NA, length = N*v), ncol = v)
          ARparam <- matrix(rep(NA, length = N*v), ncol = v)
          Duration <- matrix(rep(NA, length = N*v), ncol = v)
          colnames(ARpath)<-colnames(ARpath, do.NULL = FALSE)
          colnames(ARparam)<-colnames(ARparam, do.NULL = FALSE)
          colnames(Duration)<-colnames(Duration, do.NULL = FALSE)
        }

        for (p in 1:length(methodPathUpdate)){
          for (q in 1:length(approxTransDens)){
            for (r in 1:length(approxPropDens)){

              current_column <- (p-1)*length(methodPathUpdate)*length(approxTransDens) + (q-1)*length(approxTransDens) + r

              for(k in 1:N){
                fileName <- paste("GBM", "alpha0", alpha_0, "rho2", rho_2, "kappa0",
                                  kappa_0, "nu0", nu_0,"gamma_a", gamma_alpha,
                                  "gamma_s", gamma_sigma, "M", v_M[i], "m", v_m[j], "path", methodPathUpdate[p], "td",
                                  approxTransDens[q], "pd", approxPropDens[r], k, sep = "_")

                a <- try(load(paste("output/", fileName, ".data", sep = "")), silent = TRUE)

                if(class(a)!= "try-error"){

                  attach(results, warn.conflicts = FALSE)

                  MC_alpha <- results$theta[1, (floor(burnIn * numIterations) +1) : numIterations]
                  MC_sigma <- results$theta[2, (floor(burnIn * numIterations) +1) : numIterations]
                  meanAlpha[k, current_column] <- mean(MC_alpha)
                  meanSigma[k, current_column] <- mean(MC_sigma)
                  dens_alpha <- density(MC_alpha)
                  dens_sigma <- density(MC_sigma)
                  modeAlpha[k, current_column] <- dens_alpha$x[which.max(dens_alpha$y)]
                  modeSigma[k, current_column] <- dens_sigma$x[which.max(dens_sigma$y)]
                  hpd_alpha <- boa.hpd(MC_alpha, 0.05)
                  hpd_sigma <- boa.hpd(MC_sigma, 0.05)
                  hpd_l_Alpha[k, current_column] <- hpd_alpha[1]
                  hpd_l_Sigma[k, current_column] <- hpd_sigma[1]
                  hpd_u_Alpha[k, current_column] <- hpd_alpha[2]
                  hpd_u_Sigma[k, current_column] <- hpd_sigma[2]
                  effssAlpha[k, current_column] <- coda::effectiveSize(MC_alpha) # effective sample size
                  effssSigma[k, current_column] <- coda::effectiveSize(MC_sigma)
                  multivarESS[k, current_column] <- mcmcse::multiESS(t(results$theta[, (floor(burnIn * numIterations) +1) : numIterations]))
                  nNegPointProposals <- numNegPointProposals
                  nMBSwitchToEuler <- numMBSwitchToEuler

                  if(include_acc_rate){
                    ARpath[k, current_column] <- results$acceptanceRatePath
                    ARparam[k, current_column] <- results$acceptanceRateParam
                    Duration[k, current_column] <- results$duration[[1]]
                  }

                  if(p==1 & q==1 & r==1){
                    logIncrements <- diff(log(inputs$Y_obs))
                    dt <- diff(inputs$tau)

                    # likelihood function of the parameter
                    objectiveFct1 <- function( theta){
                      -sum(log(dnorm(logIncrements, mean = (theta[1] - 1/2 * theta[2]) * dt, sd = sqrt(theta[2] * dt))))
                    }

                    ML_result <- optim(c(2,2), objectiveFct1, method = "BFGS")
                    ML_alpha[k] <- ML_result$par[1]
                    ML_sigma[k] <- ML_result$par[2]

                    # a posteriori density of the parameter
                    objectiveFct2 <- function( theta){
                      -sum(log(c(dnorm(logIncrements, mean = (theta[1] - 1/2 * theta[2]) * dt, sd = sqrt(theta[2] * dt)),
                                dnorm(theta[1], mean = alpha_0, sd = rho_2 ),
                                dinvgamma(theta[2], shape = kappa_0, scale = nu_0))))
                    }

                    MAP_result <- optim(c(2,2), objectiveFct2, method = "BFGS")
                    MAP_alpha[k] <- MAP_result$par[1]
                    MAP_sigma[k] <- MAP_result$par[2]
                  }
                  remove(results)
                }
              }

              if(v_m[j] != 1){
                colnames(meanAlpha)[current_column] <-
                  colnames(meanSigma)[current_column] <-
                  colnames(modeAlpha)[current_column] <-
                  colnames(modeSigma)[current_column] <-
                  colnames(hpd_l_Alpha)[current_column] <-
                  colnames(hpd_l_Sigma)[current_column] <-
                  colnames(hpd_u_Alpha)[current_column] <-
                  colnames(hpd_u_Sigma)[current_column] <-
                  colnames(effssAlpha)[current_column] <-
                  colnames(effssSigma)[current_column] <-
                  paste(methodPathUpdate_short[p], "\ntd_", approxTransDens[q],
                        "\npd_", approxPropDens[r], "\nn = ", sum(!is.na(meanAlpha[,current_column])), sep = "")

                if(include_acc_rate){
                  colnames(ARpath)[current_column] <-
                    colnames(ARparam)[current_column] <-
                    colnames(Duration)[current_column] <-
                    paste(methodPathUpdate_short[p], "\ntd_", approxTransDens[q],
                          "\npd_", approxPropDens[r], "\nn = ", sum(!is.na(meanAlpha[,current_column])), sep = "")
                }
              }else{
                colnames(meanAlpha)[current_column] <-
                  colnames(meanSigma)[current_column] <-
                  colnames(modeAlpha)[current_column] <-
                  colnames(modeSigma)[current_column] <-
                  colnames(hpd_l_Alpha)[current_column] <-
                  colnames(hpd_l_Sigma)[current_column] <-
                  colnames(hpd_u_Alpha)[current_column] <-
                  colnames(hpd_u_Sigma)[current_column] <-
                  colnames(effssAlpha)[current_column] <-
                  colnames(effssSigma)[current_column] <-
                  paste("td_", approxTransDens[q],
                        "\nn = ", sum(!is.na(meanAlpha[,current_column])), sep = "")

                if(include_acc_rate){
                  colnames(ARpath)[current_column] <-
                    colnames(ARparam)[current_column] <-
                    colnames(Duration)[current_column] <-
                    paste("td_", approxTransDens[q],
                          "\nn = ", sum(!is.na(meanAlpha[,current_column])), sep = "")
                }
              }
            }
          }
        }
        colnames(ML_alpha) <- paste("ML_\nestimate\nn = ", sum(!is.na(ML_alpha)))
        colnames(ML_sigma) <- paste("ML_\nestimate\nn = ", sum(!is.na(ML_sigma)))
        colnames(MAP_alpha) <- paste("MAP_\nestimate\nn = ", sum(!is.na(MAP_alpha)))
        colnames(MAP_sigma) <- paste("MAP_\nestimate\nn = ", sum(!is.na(MAP_sigma)))

        fileName2 <- paste("GBM", "alpha", alpha, "sigma^2", sigma_2,  "M", v_M[i],
                           "m", v_m[j], "nIter", numIterations, sep = "_")

        M <- v_M[i]
        m <- v_m[j]

        if(include_acc_rate){
          save(meanAlpha, meanSigma, modeAlpha, modeSigma, hpd_l_Alpha, hpd_l_Sigma,
             hpd_u_Alpha, hpd_u_Sigma, effssAlpha, effssSigma, multivarESS,
             ARpath, ARparam, Duration,
             ML_alpha, ML_sigma, MAP_alpha, MAP_sigma, nNegPointProposals, nMBSwitchToEuler,
             M, m, N, burnIn, numIterations, alpha_0, rho_2, kappa_0,
             nu_0, gamma_alpha, gamma_sigma, true_theta,
             file = paste("aggregated_output/", fileName2, '.data', sep=''))
        }else{
          save(meanAlpha, meanSigma, modeAlpha, modeSigma, hpd_l_Alpha, hpd_l_Sigma,
               hpd_u_Alpha, hpd_u_Sigma, effssAlpha, effssSigma, multivarESS,
               ML_alpha, ML_sigma, MAP_alpha, MAP_sigma, nNegPointProposals, nMBSwitchToEuler,
               M, m, N, burnIn, numIterations, alpha_0, rho_2, kappa_0,
               nu_0, gamma_alpha, gamma_sigma, true_theta,
               file = paste("aggregated_output/", fileName2, '.data', sep=''))
        }
      }
    }
  }
}
