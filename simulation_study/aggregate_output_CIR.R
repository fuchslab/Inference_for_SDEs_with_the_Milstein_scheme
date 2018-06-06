
suppressMessages(library(boa, quietly = TRUE)) # for highest probability density interval
suppressMessages(library(MCMCpack, quietly = TRUE)) # for the inverse gamma distribution
library(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(mcmcse)

alpha <- 1
v_beta <- 1#c(1, 5)
v_sigma_2 <- 2#c( 0.25, 2)
indizes_estim_param <- c(2,3)


# number of observations to be used
v_M <- c(50)
# number of imputed points in each inter-observation interval
v_m <- c(1,2,5)
# number of observation paths
N <- 100
# burn-in rate
burnIn <- 0.1
# number of iterations
numIterations <- 1e3
# acceptance rates are only valid for the corresponding max. numIterations (1e5 or 1e6)
include_acc_rate <- FALSE

methodPathUpdate <- c("leftConditioned", "MB")
methodPathUpdate_short <- c("leftCondi", "MB")
approxTransDens <- c("Euler", "Milstein")
approxPropDens <- c("Euler", "Milstein")


# hyperparameters for the priori distributions of the parameters
# use conjugate priors
# -> beta ~ IG(kappa_b, nu_b) with
kappa_b <- 3
nu_b <- 3
# sigma_2 ~ IG(kappa_s, nu_s) inverse Gamma distribution with
kappa_s <- 3
nu_s <- 4

# hyperparameters for the random walk parameter update
gamma_b <- .5
gamma_s <- .5





for (a in 1:length(v_beta)){
  for (b in 1:length(v_sigma_2)){

    # theta <- (alpha, beta, sigma^2)
    beta <- v_beta[a]
    sigma_2 <- v_sigma_2[b]
    true_theta <- c(alpha, beta, sigma_2)

    folderName <- paste("CIR_alpha_", true_theta[1], "_beta_", true_theta[2],
                        "_sigma_", true_theta[3], sep='')


    for (i in 1:length(v_M) ){
      for (j in 1:length(v_m)){
        v <- length(approxPropDens) * length(approxTransDens) * length(methodPathUpdate)
        meanBeta <- matrix(rep(NA, length = N*v), ncol = v)
        meanSigma <- matrix(rep(NA, length = N*v), ncol = v)
        modeBeta <- matrix(rep(NA, length = N*v), ncol = v)
        modeSigma <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_l_Beta <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_l_Sigma <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_u_Beta <- matrix(rep(NA, length = N*v), ncol = v)
        hpd_u_Sigma <- matrix(rep(NA, length = N*v), ncol = v)
        effssBeta <- matrix(rep(NA, length = N*v), ncol = v) # effective sample size
        effssSigma <- matrix(rep(NA, length = N*v), ncol = v)
        multivarESS <- matrix(rep(NA, length = N*v), ncol = v)
        colnames(meanBeta)<-colnames(meanBeta, do.NULL = FALSE)
        colnames(meanSigma)<-colnames(meanSigma, do.NULL = FALSE)
        colnames(modeBeta)<-colnames(modeBeta, do.NULL = FALSE)
        colnames(modeSigma)<-colnames(modeSigma, do.NULL = FALSE)
        colnames(hpd_l_Beta)<-colnames(hpd_l_Beta, do.NULL = FALSE)
        colnames(hpd_l_Sigma)<-colnames(hpd_l_Sigma, do.NULL = FALSE)
        colnames(hpd_u_Beta)<-colnames(hpd_u_Beta, do.NULL = FALSE)
        colnames(hpd_u_Sigma)<-colnames(hpd_u_Sigma, do.NULL = FALSE)
        colnames(effssBeta)<-colnames(effssBeta, do.NULL = FALSE)
        colnames(effssSigma)<-colnames(effssSigma, do.NULL = FALSE)
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
                fileName <- paste("CIR", "kappa_b", kappa_b, "nu_b", nu_b, "kappa_s",
                                  kappa_s, "nu_s", nu_s,"gamma_b", gamma_b,
                                  "gamma_s", gamma_s, "M", v_M[i], "m", v_m[j], "path", methodPathUpdate[p], "td",
                                  approxTransDens[q], "pd", approxPropDens[r], k, sep = "_")

                a <- try(load(paste("simulation_study/", folderName, "/output/",
                                    fileName, ".data", sep = "")), silent = TRUE)

                if(class(a)!= "try-error"){

                  attach(results, warn.conflicts = FALSE)

                  MC_beta <- results$theta[2, (floor(burnIn * numIterations) +1) : numIterations]
                  MC_sigma <- results$theta[3, (floor(burnIn * numIterations) +1) : numIterations]
                  meanBeta[k, current_column] <- mean(MC_beta)
                  meanSigma[k, current_column] <- mean(MC_sigma)
                  dens_beta <- density(MC_beta)
                  dens_sigma <- density(MC_sigma)
                  modeBeta[k, current_column] <- dens_beta$x[which.max(dens_beta$y)]
                  modeSigma[k, current_column] <- dens_sigma$x[which.max(dens_sigma$y)]
                  hpd_beta <- boa.hpd(MC_beta, 0.05)
                  hpd_sigma <- boa.hpd(MC_sigma, 0.05)
                  hpd_l_Beta[k, current_column] <- hpd_beta[1]
                  hpd_l_Sigma[k, current_column] <- hpd_sigma[1]
                  hpd_u_Beta[k, current_column] <- hpd_beta[2]
                  hpd_u_Sigma[k, current_column] <- hpd_sigma[2]
                  effssBeta[k, current_column] <- coda::effectiveSize(MC_beta) # effective sample size
                  effssSigma[k, current_column] <- coda::effectiveSize(MC_sigma)
                  multivarESS[k, current_column] <- mcmcse::multiESS(t(results$theta[indizes_estim_param, (floor(burnIn * numIterations) +1) : numIterations]))
                  nNegPointProposals <- numNegPointProposals
                  nMBSwitchToEuler <- numMBSwitchToEuler

                  if(include_acc_rate){
                    ARpath[k, current_column] <- results$acceptanceRatePath
                    ARparam[k, current_column] <- results$acceptanceRateParam
                    Duration[k, current_column] <- results$duration[[1]]
                  }

                  remove(results)
                }
              }

              if(v_m[j] != 1){
                colnames(meanBeta)[current_column] <-
                  colnames(meanSigma)[current_column] <-
                  colnames(modeBeta)[current_column] <-
                  colnames(modeSigma)[current_column] <-
                  colnames(hpd_l_Beta)[current_column] <-
                  colnames(hpd_l_Sigma)[current_column] <-
                  colnames(hpd_u_Beta)[current_column] <-
                  colnames(hpd_u_Sigma)[current_column] <-
                  colnames(effssBeta)[current_column] <-
                  colnames(effssSigma)[current_column] <-
                  paste(methodPathUpdate_short[p], "\ntd_", approxTransDens[q],
                        "\npd_", approxPropDens[r], "\nn = ", sum(!is.na(meanBeta[,current_column])), sep = "")

                if(include_acc_rate){
                  colnames(ARpath)[current_column] <-
                    colnames(ARparam)[current_column] <-
                    colnames(Duration)[current_column] <-
                    paste(methodPathUpdate_short[p], "\ntd_", approxTransDens[q],
                          "\npd_", approxPropDens[r], "\nn = ", sum(!is.na(meanBeta[,current_column])), sep = "")
                }
              }else{
                colnames(meanBeta)[current_column] <-
                  colnames(meanSigma)[current_column] <-
                  colnames(modeBeta)[current_column] <-
                  colnames(modeSigma)[current_column] <-
                  colnames(hpd_l_Beta)[current_column] <-
                  colnames(hpd_l_Sigma)[current_column] <-
                  colnames(hpd_u_Beta)[current_column] <-
                  colnames(hpd_u_Sigma)[current_column] <-
                  colnames(effssBeta)[current_column] <-
                  colnames(effssSigma)[current_column] <-
                  paste("td_", approxTransDens[q],
                        "\nn = ", sum(!is.na(meanBeta[,current_column])), sep = "")

                if(include_acc_rate){
                  colnames(ARpath)[current_column] <-
                    colnames(ARparam)[current_column] <-
                    colnames(Duration)[current_column] <-
                    paste("td_", approxTransDens[q],
                          "\nn = ", sum(!is.na(meanBeta[,current_column])), sep = "")
                }
              }
            }
          }
        }
        fileName2 <- paste("CIR_alpha", true_theta[1], "beta", true_theta[2],
                           "sigma", true_theta[3],  "M", v_M[i],
                           "m", v_m[j], "nIter", numIterations, sep = "_")

        M <- v_M[i]
        m <- v_m[j]

        if(include_acc_rate){
          save(meanBeta, meanSigma, modeBeta, modeSigma, hpd_l_Beta, hpd_l_Sigma,
             hpd_u_Beta, hpd_u_Sigma, effssBeta, effssSigma, multivarESS,
             ARpath, ARparam, Duration,
             nNegPointProposals, nMBSwitchToEuler,
             M, m, N, burnIn, numIterations, kappa_b, nu_b, kappa_s,
             nu_s, gamma_b, gamma_s, true_theta,
             file = paste("simulation_study/aggregated_output/", fileName2, '.data', sep=''))
        }else{
          save(meanBeta, meanSigma, modeBeta, modeSigma, hpd_l_Beta, hpd_l_Sigma,
               hpd_u_Beta, hpd_u_Sigma, effssBeta, effssSigma, multivarESS,
               nNegPointProposals, nMBSwitchToEuler,
               M, m, N, burnIn, numIterations, kappa_b, nu_b, kappa_s,
               nu_s, gamma_b, gamma_s, true_theta,
               file = paste("simulation_study/aggregated_output/", fileName2, '.data', sep=''))
        }
      }
    }
  }
}
