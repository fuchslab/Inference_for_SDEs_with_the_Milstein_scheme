# Parameter estimation for the Itô diffusion process
#
#  which fulfills the SDE dX_t = mu(X_t,theta)*dt + sigma(X_t,theta)*dB_t, X_0=x_0
#
# version : not saving all the path updates
#
# path proposal: "leftConditioned", "MB"
# parameter proposal: "RandomWalk"
# approximation method for transition density (approxTransDens): "Euler", "Milstein"
# approximation method for proposal densities (approxPropDens): "Euler", "Milstein"
# numIterations - number of iterations (MCMC steps)
# m - number of subintervals between two observations (m-1 points are imputed between two observations)
# Y_obs - observed points
# tau - time points of observations


library(MCMCpack, quietly = TRUE) # for the inverse gamma distribution
library(coda, quietly = TRUE)
library(MASS, quietly = TRUE)
library(boa, quietly = TRUE) # for highest probability density interval

estimate_parameters <- function(methodPathUpdate = "leftConditioned",
                     methodParamUpdate = "RandomWalk",
                     approxTransDens = "Euler",
                     approxPropDens = "Euler",
                     numIterations = 10^5, m=2,
                     Y_obs, tau){

  envir <- new.env()
  if (approxTransDens == "Euler" ){
    envir$dtransDens <- dEuler
    envir$rtransDens <- rEuler
  }else if (approxTransDens == "Milstein" ){
    envir$dtransDens <- dMilstein
    envir$rtransDens <- rMilstein
  }else{
    envir$error("The input value of approxTransDens is not valid.")
  }

  if (approxPropDens == "Euler" ){
    envir$dpropDens <- dEuler
    envir$rpropDens <- rEuler
    envir$dMBpropDens <- dMBEuler
    envir$rMBpropDens <- rMBEuler
  }else if (approxPropDens == "Milstein" ){
    envir$dpropDens <- dMilstein
    envir$rpropDens <- rMilstein
    if(m == 2){
      envir$dMBpropDens <- dMBMilstein_unnormalized
    }else{
      envir$dMBpropDens <- dMBMilstein
    }
    envir$rMBpropDens <- rMBMilstein
  }else{
    stop("The input value of approxPropDens is not valid.")
  }

  if( length(Y_obs) != length(tau)){
    stop("The number of observations in Y_obs is not equal to the number of time points in tau.")
  }

  numProposedPaths <- 0 # counts the number of proposed paths
  numAcceptedPaths <- 0 # counts the number of accepted path proposals
  numAcceptedParam <- 0 # counts the number of accepted parameter proposals
  envir$countNegPointProposals <- generateCountingFct() # counts the number of negative proposals during the path update
  envir$countMBSwitchToEuler <- generateCountingFct() # counts the number of times MBEuler had to be used instead of MBMilstein during the path update

  Time <- tau[length(tau)]
  M <- length(tau) - 1 # number of inter-observation intervals
  delta_t <- Time / (M * m) # time step for the augmented path
  tpoints <- seq(from = 0, to = Time, by = delta_t)  # time points for the augmented path
  tpoints[(0:(length(tau)-1))*m+1] <- tau # in order to avoid numerical issues, make
  # sure that tpoints contains the exact same values as tau

  # 1) Initialize Y_imp by linear interpolation
  Y_imp_current <- matrix(approx(tau, Y_obs, tpoints, method = "linear")$y, ncol = 1)

  # 2) Draw initial values for theta
  theta <- rprior_theta(1)

  # memory pre-allocation
  # save numSavedPaths of the accepted paths, the parameter estimates and the log-posterior p(theta | Y_obs, Y_imp)
  numSavedPaths <- 10
  Y_imp <-  matrix( c(Y_imp_current,
                      vector(mode = "numeric",
                             length = (numSavedPaths * nrow(Y_imp_current)))),
                    nrow = nrow(Y_imp_current))
  theta <-  matrix( c(theta, vector(mode = "numeric",
                                   length = (numIterations * length(theta)))),
                    nrow = length(theta))
  log_posterior <- vector(mode = "numeric", length = numIterations+1)
  log_posterior[1] <- sum(log(envir$dtransDens(Y_imp[2:length(Y_imp[,1]), 1],
                                            Y_imp[1:(length(Y_imp[,1])-1),1], theta,
                                            delta_t )))

  # 3) Repeat the path update and parameter update
  while( numProposedPaths < numIterations ) {
    # choice of update interval according to algorithm 7.3 (p. 191 in Fuchs(2013))
    update_intervals <- Alg_7_3(0, M * m, lambda)

    # for each update interval the fct update_path proposes a new path and accepts it
    # with the according probability and then parameter is updated after each update interval
    for ( i in 1:(length(update_intervals) - 1) )    {
      a <- update_intervals[i] + 1 # +1 because the interval limits start at 0
      b <- update_intervals[i+1] + 1

      if ( (b - a > 1) && (numProposedPaths < numIterations) ){ # makes sure, that there are points to be imputed
        if(m > 1){
          # Path update ----------------------------------------------------------------
          updates <- update_path(a, b, delta_t, tau, Y_obs, tpoints, Y_imp_current,
                                 methodPathUpdate, theta[, numProposedPaths + 1],
                                 env = envir)

          # if the proposed path was rejected, updates will be empty
          if (length(updates) != 0)        {
            Y_imp_current <- updates$path
            numProposedPaths <- numProposedPaths + updates$numProposed
            numAcceptedPaths <- numAcceptedPaths + updates$numAccepted
          }
          else{
            numProposedPaths <- numProposedPaths + 1
          }

          # save selection of paths
          if ( (numProposedPaths %% (numIterations / numSavedPaths)) == 0 ){
            Y_imp[, numProposedPaths * numSavedPaths / numIterations + 1] <- Y_imp_current
          }
        }else{# m==1, estimation without data imputation
          numProposedPaths <- numProposedPaths + 1
        }

        # Parameter Update -----------------------------------------------------------
        param_update <- update_parameter(theta[, max(c(numProposedPaths,1))],
                                         gamma_beta, gamma_sigma, beta_0, rho_2,
                                         kappa_0, nu_0, delta_t, Y_imp_current,
                                         methodParamUpdate, env = envir)

        theta[, numProposedPaths + 1] <- param_update$theta
        numAcceptedParam <- numAcceptedParam + param_update$numAccepted
        log_posterior[numProposedPaths + 1] <-  param_update$log_posterior
      }
    }
  }

  acceptanceRatePath <- numAcceptedPaths / numProposedPaths
  acceptanceRateParam <- numAcceptedParam / numIterations

  return( list(theta = theta, Y_imp = Y_imp, tpoints = tpoints,
               acceptanceRatePath = acceptanceRatePath,
               acceptanceRateParam = acceptanceRateParam,
               numNegPointProposals = envir$countNegPointProposals(getCount = TRUE),
               numMBSwitchToEuler = envir$countMBSwitchToEuler(getCount = TRUE),
               log_posterior = log_posterior))
}
