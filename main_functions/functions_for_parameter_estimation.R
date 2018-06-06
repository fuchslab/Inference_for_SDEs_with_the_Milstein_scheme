# Algorithm 7.3. from CF13 p. 191 to choose the update interval
# alteration: arbitrary starting time t_0

Alg_7_3 <- function(t_0, t_end, lambda){
  c <- t_0
  j <- 1
  while (c[j] < t_end){
    Z <- rpois(1, lambda)
    c <- c(c, min(c[length(c)] + Z, t_end))
    j <- j + 1
  }
  return(c)
}
#-------------------------------------------------------------------------------------
# closures to count certain events
generateCountingFct <- function(){
  globalvar <- 0
  countingFct <- function(getCount = FALSE) {
    if(getCount){
      globalvar
    }else{
      globalvar <<- globalvar + 1
    }
  }
  return(countingFct)
}

#-------------------------------------------------------------------------------------
evaluate_quotient <- function(numerator, denominator){
  if (length(numerator) != length(denominator)){
    stop("numerator and denominator are not of the same length")
  }
  # check whether any of the entries are infinite
  if (any(is.infinite(c(numerator, denominator)))){
    #browser()
    # determine the number of Inf and -Inf entries of the numerator and denominator
    numNumerNegInf <- sum(is.infinite(numerator) * (numerator < 0))
    numNumerPosInf <- sum(is.infinite(numerator) * (numerator > 0))
    numDenomNegInf <- sum(is.infinite(denominator) * (denominator < 0))
    numDenomPosInf <- sum(is.infinite(denominator) * (denominator > 0))

    # compare the "net" number of infinite entries
    if ( (numNumerPosInf - numNumerNegInf) > (numDenomPosInf - numDenomNegInf)){
      return(Inf)
    }else if ( (numNumerPosInf - numNumerNegInf) < (numDenomPosInf - numDenomNegInf)){
      return(-Inf)
    }else{
      # net number of infinite entries are the same for numerator and denominator
      # -> exclude all points with infinite entries in either denomitator or numerator
      quotient <- numerator - denominator
      return( sum(quotient[!(is.infinite(quotient) + is.nan(quotient))]) )
    }
  }else{
    return( sum(numerator - denominator) )
  }

}

#-------------------------------------------------------------------------------------
rEuler <- function(n, x, theta, delta_t){
  rnorm(n = n, mean = x + drift_fct(x, theta) * delta_t,
             sd = sqrt(diffusion_fct(x, theta) ^ 2 * delta_t) )
}

dEuler <- function(y, x, theta, delta_t){
  dnorm(y, mean = x + drift_fct(x, theta) * delta_t,
             sd = sqrt(diffusion_fct(x, theta) ^ 2 * delta_t) )
}

rMilstein <- function(n, x, theta, delta_t){
  dB <- rnorm(n = n, mean = 0, sd = sqrt(delta_t)) # increment of a Brownian motion
  diffusion <- diffusion_fct(x, theta)
  return( x + drift_fct(x, theta) * delta_t + diffusion * dB +
            1/2 * diffusion * diffusion_fct_derivative(x, theta) * (dB ^ 2 - delta_t))
}

dMilstein <- function(y, x, theta, delta_t){

  if(length(x) == 1){
    drift <- drift_fct(x, theta)
    diffusion <- diffusion_fct(x, theta)
    diffusionDeriv <- diffusion_fct_derivative(x, theta)

    if(diffusion*diffusionDeriv > 0){
      min_y <- x - 1/2 * diffusion / diffusionDeriv +
        (drift - 1/2 * diffusion * diffusionDeriv) * delta_t
      condition <- (y > min_y)
    }else{
      max_y <- x - 1/2 * diffusion / diffusionDeriv +
        (drift - 1/2 * diffusion * diffusionDeriv) * delta_t
      condition <- (y < max_y)
    }

    result <- vector(mode = 'numeric', length(y)) # contains only zeros for now
    # only the components that fulfill the condition will be replaced
    A <- sqrt( diffusion ^ 2 + 2 * diffusion * diffusionDeriv *
                 (y[condition] - x - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
    B <-  - (diffusion + diffusionDeriv *
               (y[condition] - x - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
    C <- diffusion * diffusionDeriv ^ 2 * delta_t
    dens <- (exp((B - A) / C) + exp((B + A) / C)) / (sqrt(2 * pi * delta_t) * A)

    result[condition] <- dens
    return(result)
  }else{
    if(length(y) != length(x)){
      error("input vectors for the left point and the next point are not of same length!")
    }else{
      result <- vector(mode = "numeric", length = length(x))
      for (i in 1:length(x)){
        result[i] <- dMilstein(y[i], x[i], theta, delta_t)
      }
      return(result)
    }
  }
}

dMBEuler <- function(y, x_k, x_b, theta, delta_t, left_index, right_index){
  return(dnorm(y, mean = x_k + (x_b - x_k) / (right_index - left_index),
              sd = sqrt( (right_index - left_index - 1) / (right_index - left_index) *
                            diffusion_fct(x_k, theta) ^ 2 * delta_t)))
}

rMBEuler <- function(n, x_k, x_b, theta, delta_t, left_index, right_index, precision){
  return(rnorm(n, mean =  x_k + (x_b - x_k) / (right_index - left_index),
               sd = sqrt( (right_index - left_index - 1) / (right_index - left_index) *
                            diffusion_fct(x_k, theta) ^ 2 * delta_t)))
}



#-----------------------------------------------------------------------------------------------------------------
# uses log-transformed densities to calculate the acceptance probabilities
update_path <- function(a, b, delta_t, tau, Y_obs, tpoints, Y_imp, method,
                        theta, env = parent.frame()) {
  current_index <- a
  proposal <- Y_imp
  pointWasProposed <- 0 # binary variable to check whether there was actually a new point proposed
  accept_prob <- 1 # if no point is proposed 'accept_prob' stays equal to 1 and the
  # proposed remains the old path and is always accepted in the end

  if ( b - a > 1){ # makes sure, that there are points to be imputed
    if ( method == "leftConditioned"){
      # initialize variable for the sum of the log transition densities
      # which is needed for the acceptance probability
      numerator <- vector(mode = "numeric", length = 0)
      denominator <- vector(mode = "numeric", length = 0)

      # while there are points to impute between t_a and t_b;
      while ( current_index + 1 < b){
        # if the next point is NOT one of the observed points
        if ( length(which(tpoints[current_index + 1] == tau)) == 0){
          # proposal conditioned on the left point for the next imputed point
          proposal[current_index + 1] <- env$rpropDens(1, proposal[current_index],
                                                       theta, delta_t)

          while(proposal[current_index + 1] < 0){
            env$countNegPointProposals()
            proposal[current_index + 1] <- env$rpropDens(1, proposal[current_index],
                                                         theta, delta_t)
          }

          pointWasProposed <- 1
        }
        else if (pointWasProposed == 1){
          # the next point is an observed point and a new point was
          # previously proposed -> calculate a factor for the acceptance probability
          numerator <- c(numerator, log(env$dtransDens(Y_imp[current_index + 1],
                                              proposal[current_index], theta, delta_t)))
          denominator <- c(denominator, log(env$dtransDens(Y_imp[current_index + 1],
                                                   Y_imp[current_index], theta, delta_t)))
        }
        current_index <- current_index + 1

      }
      # now: current_index+1 = b
      # the last point (current_index) is NOT an observed point and a new point was
      # previously proposed -> calculate a factor for the acceptance probability
      if (pointWasProposed == 1 & ( length(which(tpoints[current_index] == tau)) == 0)){
        numerator <- c(numerator, log(env$dtransDens(Y_imp[current_index + 1],
                                             proposal[current_index], theta, delta_t)))
        denominator <- c(denominator, log(env$dtransDens(Y_imp[current_index + 1],
                                                 Y_imp[current_index], theta, delta_t)))
      }

      # evaluate the quotient taking into account whether any of the entries are infinite
      sum1 <- evaluate_quotient(numerator, denominator)

      # Calculate acceptance probability ("leftConditioned")
      # if no point was proposed sum1 = 0 -> accept_prob = 1
      accept_prob <- min(c(1, exp(sum1)))

    } # end if "leftConditioned"
    else if ( method == "MB"){ # Modified Bridge
      # initialize variable for the terms of the log proposal densities
      # which is needed for the acceptance probability
      denominator_prop_dens <- vector(mode = "numeric", length = 0)
      numerator_prop_dens <- vector(mode = "numeric", length = 0)

      # while there are points to be imputed between t_a and t_b;
      while (current_index + 1 < b){
        # if the next point is NOT one of the observed points, make a new proposal for it
        if (length(which(tpoints[current_index + 1] == tau)) == 0){
          # find the index of the next observation to the right of the point to be proposed,
          # i.e. the right point on which the modified bridge proposal will be conditioned
          b_tau <- min( which(tpoints[current_index + 1] < tau) )
          b_right <- which(tpoints == tau[b_tau])

          # if the next observation is not within the update interval, condition
          # on the right end point b of the update interval
          if ( b_right > b ){
            b_right <- b
          }

          # Modified Bridge proposal for the next imputed point
          proposal[current_index + 1]  <-
            env$rMBpropDens(n = 1, x_k = proposal[current_index], x_b = Y_imp[b_right],
                     theta = theta, delta_t = delta_t, left_index = current_index,
                     right_index = b_right)
          if(!is.na(proposal[current_index + 1])){
            # if the proposed point is negative, make a new proposal
            while(proposal[current_index + 1] < 0){
              env$countNegPointProposals()
              proposal[current_index + 1] <-
                env$rMBpropDens(n = 1, x_k = proposal[current_index], x_b = Y_imp[b_right],
                                theta = theta, delta_t = delta_t, left_index = current_index,
                                right_index = b_right)
            }
            # terms of the quotient of the log proposal densities
            # (to calculate the acceptance probability later)
            numerator_prop_dens <-
              c(numerator_prop_dens,
                log( env$dMBpropDens(y = Y_imp[current_index + 1], x_k = Y_imp[current_index],
                                     x_b = Y_imp[b_right], theta = theta, delta_t = delta_t,
                                     left_index = current_index, right_index = b_right)))
            denominator_prop_dens <-
              c(denominator_prop_dens,
                log( env$dMBpropDens(y = proposal[current_index + 1], x_k = Y_imp[current_index],
                                     x_b = Y_imp[b_right], theta = theta, delta_t = delta_t,
                                     left_index = current_index, right_index = b_right)))
          }else{# if rMBMilstein returns NA, the set of feasible points is empty -> switch to MBEuler
            proposal[current_index + 1] <- rMBEuler(n = 1, x_k = proposal[current_index], x_b = Y_imp[b_right],
                                                    theta = theta, delta_t = delta_t, left_index = current_index,
                                                    right_index = b_right)
            env$countMBSwitchToEuler()
            # if the proposed point is negative, make a new proposal
            while(proposal[current_index + 1] < 0){
              env$countNegPointProposals()
              proposal[current_index + 1] <-
                rMBEuler(n = 1, x_k = proposal[current_index], x_b = Y_imp[b_right],
                         theta = theta, delta_t = delta_t, left_index = current_index,
                         right_index = b_right)
            }
            # terms of the quotient of the log proposal densities
            numerator_prop_dens <-
              c(numerator_prop_dens,
                log( dMBEuler(y = Y_imp[current_index + 1], x_k = Y_imp[current_index],
                              x_b = Y_imp[b_right], theta = theta, delta_t = delta_t,
                              left_index = current_index, right_index = b_right)))
            denominator_prop_dens <-
              c(denominator_prop_dens,
                log( dMBEuler(y = proposal[current_index + 1], x_k = Y_imp[current_index],
                              x_b = Y_imp[b_right], theta = theta, delta_t = delta_t,
                              left_index = current_index, right_index = b_right)))
          }
          pointWasProposed <- 1
        }
        current_index <- current_index + 1
      }

      # Calculate acceptance probability (Modified Bridge)
      # terms of the quotient of the (Euler/Milstein) approximated log transition densities
      numerator_posterior <- log(env$dtransDens(proposal[(a+1):b], proposal[a:(b-1)],
                                        theta, delta_t))
      denominator_posterior <- log(env$dtransDens(Y_imp[(a+1):b], Y_imp[a:(b-1)], theta,
                                          delta_t))

      # evaluate the quotient (/sum in log scale) taking into account whether any of the entries are infinite
      sum1 <- evaluate_quotient(c(numerator_prop_dens, numerator_posterior),
                                  c(denominator_prop_dens, denominator_posterior))

      # acceptance probability (Modified Bridge)
      accept_prob <- min(c(1, exp(sum1) ))

    } # end if Modified Bridge

    if(is.na(accept_prob)|is.null(accept_prob)) browser()

    # accept or reject?
    if ( runif(1) < accept_prob ){
      # path is accepted (accept_prob=1 if no point was proposed)
      # pointWasProposed is equal to 1 if there was actually new points proposed and 0 otherwise
      return(list( path = proposal, numProposed = pointWasProposed, numAccepted = pointWasProposed))
    }else{
      # path is rejected
      return()
    }
  }else{
    # no points to be imputed
    return()
  }
}# end update_path


#-----------------------------------------------------------------------------------------------------------------
update_parameter <- function(theta, gamma_beta, gamma_sigma,
                             beta_0, rho_2, kappa_0, nu_0, delta_t, Y_imp, method,
                             env = parent.frame()){
  # random walk proposal
  proposal_theta <- propose_theta(theta)

  # terms of the quotient of the log posterior densities of the path (using Euler/Milstein approximation)
  numerator <- log(env$dtransDens(Y_imp[2:length(Y_imp)], Y_imp[1:(length(Y_imp)-1)],
                         proposal_theta, delta_t ) )
  denominator <- log(env$dtransDens(Y_imp[2:length(Y_imp)], Y_imp[1:(length(Y_imp)-1)], theta,
                         delta_t ) )
  # log prior density of alpha and sigma_2
  log_prior_theta <- log_dprior_theta(theta)
  log_prior_proposal_theta <- log_dprior_theta(proposal_theta)

  # evaluate the quotient taking into account whether any of the entries are infinite
  sum1 <- evaluate_quotient(c(numerator, log_prior_proposal_theta, log(theta_positive(proposal_theta))),
                            c(denominator, log_prior_theta, log(theta_positive(theta))))

  # calculate the acceptance probability
  accept_prob <- min(c(1, exp(sum1)))

  # accept or reject?
  if ( runif(1) < accept_prob){
    # return the new state, the count that the proposal was accepted, and the current
    # value of the log-posterior
    return(list(theta = proposal_theta, numAccepted = 1,
                log_posterior = sum(c(numerator, log_prior_proposal_theta))))
  }else{
    return(list(theta = theta, numAccepted = 0,
                log_posterior = sum(c(denominator, log_prior_theta))))
  }

}
