# input of the problem specific parameters and functions
#
# Geometric Brownian motion
#
#  which fulfills the SDE dX_t = alpha*X_t*dt + sigma*X_t*dB_t, X_0=x_0

#-------------------------------------------------------------------------------------
# number of parameters (corresponds to the length of theta)
len_theta <- 2
thetaNames <- c("alpha", "sigma^2")
# indizes of the parameters that are to be estimated
indizes_estim_param <- c(1,2)
# value of the parameter that is not estimated
# none
#-------------------------------------------------------------------------------------
drift_fct <- function(x, theta){
  theta[1] * x
}

diffusion_fct <- function(x, theta){
  sqrt(theta[2]) * x
}

diffusion_fct_derivative <- function(x, theta){
  rep(sqrt(theta[2]), length(x))
}


#-------------------------------------------------------------------------------------
theta_positive <- function(theta){
  # returns the components of theta with strictly positve range
  # here: sigma^2
  return(theta[2])
}

# ------------------------------------------------------------------------------------
# for Alg. 7.3. for the choice of the update intervals
lambda <- 5

#-------------------------------------------------------------------------------------
# hyperparameters for the priori distributions of the parameters
# use conjugate priors
# -> alpha ~ N(alpha_0,rho_2) with
alpha0 <- 0
rho2 <- 10
# sigma_2 ~ IG(kappa_0,nu_0) inverse Gamma distribution with
kappa0 <- 2
nu0 <- 2

# hyperparameters for the random walk parameter update
gamma_alpha <- .5
gamma_sigma <- .5

hyperparam <- list(alpha0 = alpha0, rho2 = rho2, kappa0 = kappa0, nu0 = nu0,
                   gamma_alpha = gamma_alpha, gamma_sigma = gamma_sigma)

#-------------------------------------------------------------------------------------
rprior_theta <- function(n){
  # draws a sample of size n from the assumed prior of the parameter theta
  theta <- matrix(vector(mode = "numeric", length = len_theta * n), ncol = n)
  # alpha
  theta[1, ] <- rnorm(n = n, mean = alpha0, sd = sqrt(rho2))
  # sigma^2
  theta[2, ] <- rinvgamma(n = n, shape = kappa0, scale = nu0)

  return(theta)
}

propose_theta <- function(theta){
  proposal_theta <- vector(mode = "numeric", length = length(theta))

  # random walk proposal
  # propose new parameters according to normal and log normal distribution
  #alpha
  proposal_theta[1] <- rnorm(1, mean = theta[1], sd = gamma_alpha)
  # sigma^2
  proposal_theta[2] <- rlnorm(1, mean = log(theta[2]), sd = gamma_sigma)

  return(proposal_theta)
}

log_dprior_theta <- function(theta){
  prior_theta <- vector(mode = "numeric", length = length(theta))
  # prior density of alpha and sigma_2
  # alpha
  prior_theta[1] <- dnorm(theta[1], mean = alpha0, sd = sqrt(rho2) )
  # sigma^2
  prior_theta[2] <- dinvgamma(theta[2], shape = kappa0, scale = nu0)
  return(sum(log(prior_theta)))
}

#-------------------------------------------------------------------------------------
# functions needed for the parameter estimation procedure that combines the modified
# bridge proposal and the Milstein scheme cannot be easily generalized for all processes
# and are therefore implemented individually here
dMBMilstein <- function(y, x_k, x_b, theta, delta_t, left_index, right_index,
                        precision = 10^4){
  if( any(c(y, x_k, x_b) <= 0)){ # GBM-specific
    stop("The function dMBMilstein() is only valid for positive input states y, x_k, x_b.\n",
         " y = ", y,
         ", x_k = ", x_k,
         ", x_b = ", x_b
    )
  }
  stepsToRight <- right_index - left_index - 1
  drift <- drift_fct(x_k, theta)
  diffusion <- diffusion_fct(x_k, theta)
  diffusionDeriv <- diffusion_fct_derivative(x_k, theta)

  ## determine the limits of the support of the density
  # simplified GBM-specific !!!!  limits for y
  # assumes positive state space => diffusion > 0
  min_y  <- max(c( 0, x_k - 1/2 * diffusion / diffusionDeriv +
                     (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
  divisor <- 1/2 + (theta[1] - 1/2 * theta[2]) * delta_t * stepsToRight
  if(divisor >= 0){
    max_y <- x_b / divisor
  }else{ # divisor < 0
    min_y <- max(c(x_b / divisor, min_y))
    max_y <- Inf
  }

  if( min_y > max_y ){
    warning(paste("Set of feasible points y is empty for the given parameters:\n min_y = ",
                  min_y, ", max_y = ", max_y, sep=""))
    return(NA)
  }else if(isTRUE(all.equal(min_y, max_y))){
    warning(paste("The limits of the set of feasible points y are equal to each other:\n min_y = ",
                  min_y, ", max_y = ", max_y, sep=""))
    result <- vector(mode = "numeric", length = length(y))
    for(i in 1:length(y)){ if(isTRUE(all.equal(max_y, y[i]))){
      result[i] <- Inf
    }
    }
    return(result)
  }else{
    dMBMilstein_unnormalized <- function(y, x_k, x_b, theta, delta_t, stepsToRight, scale = 1){

      A_squared <- diffusion^2 + 2 * diffusion * diffusionDeriv *
        (y - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t)
      # due to numercial issues, A_squared may become negative (but small) instead of zero
      cond1 <- A_squared < 0 & abs(A_squared) < 1e-11
      if(any(cond1)){
        A <- vector(mode = "numeric", length = length(y))
        A[cond1] <- 0
        if(any(!cond1)) A[!cond1] <- sqrt(A_squared[!cond1])
      }else{
        A <- sqrt(A_squared)
      }
      B <- - (diffusion + diffusionDeriv *
                (y - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
      C <- diffusion * diffusionDeriv^2 * delta_t

      drift_y <- drift_fct(y, theta)
      diffusion_y <- diffusion_fct(y, theta)
      diffusionDeriv_y <- diffusion_fct_derivative(y, theta)

      # general implementation:
      # A_y_squared <- diffusion_y ^ 2 + 2 * diffusion_y * diffusionDeriv_y *
      #   (x_b - y - (drift_y - 1/2 * diffusion_y * diffusionDeriv_y) *
      #      delta_t * stepsToRight)
      # GBM specific (!) implementation to avoid numerical issues:
      A_y_squared <- y * 2 * theta[2] * x_b -
        y^2 * theta[2] *
        (1 + 2 * (theta[1] - 1/2 * theta[2]) * delta_t * stepsToRight)

      # due to numercial issues, A_y_squared may become negative (but small) instead of zero
      cond2 <- A_y_squared < 0 & abs(A_y_squared) < 1e-11
      if(any(cond2)){
        A_y <- vector(mode = "numeric", length = length(y))
        A_y[cond2] <- 0
        if(any(!cond2)) A_y[!cond2] <- sqrt(A_y_squared[!cond2])
      }else{
        A_y <- sqrt(A_y_squared)
      }
      B_y <- - (diffusion_y + diffusionDeriv_y *
                  (x_b - y - (drift_y - 1/2 * diffusion_y * diffusionDeriv_y) *
                     delta_t * stepsToRight))
      C_y <- diffusion_y * diffusionDeriv_y^2 * delta_t * stepsToRight
      # rescale the result by the factor 'scale' to avoid numerical issues when the function is integrated
      return(scale * (exp((B - A)/C) + exp((B + A)/C)) / (sqrt(2 * pi * delta_t) * A) *
               (exp((B_y - A_y)/C_y) + exp((B_y + A_y)/C_y)) /
               (sqrt(2 * pi * delta_t * stepsToRight) * A_y))
    }

    ## determine the maximum of dMBMilstein_unnormalized() (for rescaling)
    # auxiliary variables to substitute Inf by very large number
    max_y_aux <- min(c(max_y, 1e16))
    if(is.infinite(max_y)){
      # use double-logarithmic scaling
      y_aux <- min_y + exp(exp(seq(from = 0, to =  log(log(max_y_aux - min_y + 1) +1), by = log(log((max_y_aux - min_y))) / precision)) - 1) -1
    }else{
      y_aux <- seq(from = min_y, to = max_y_aux, by = (max_y_aux - min_y) / precision)
    }
    # leave out the interval limits, as the density may become infinite/NA there
    dens_values <- dMBMilstein_unnormalized(y_aux[2:(length(y_aux)-1)], x_k, x_b, theta, delta_t,
                                            stepsToRight)
    maxdens <- max(dens_values, na.rm = TRUE)
    if (is.na(maxdens) | maxdens == 0){
      warning(paste("The maximum of the calculated density is equal to ",
                    maxdens, ",\nwherefore the density is not valid.", sep=""))
      return(NA)
    }else if ((maxdens <= 1e-320) && (var(dens_values[dens_values>0], na.rm = TRUE) == 0)){
      warning(paste("The calculated (unnormalized) density seems to be constant with a value of ",
                    maxdens, "\nwhich indicates that it is not correctly calculated due to numerical issues.", sep=""))
      return(NA)
    }else{
      ## calcualte the normalization constant const
      # if the scale of the unnormalized density is too small, integration will not be correct
      # -> rescale before integration
      try(const <- integrate(f = dMBMilstein_unnormalized, lower = min_y, upper = max_y,  x_k, x_b,
                             theta, delta_t, stepsToRight, scale = 10/maxdens)$value * maxdens/10, silent = TRUE )
      if(!exists("const", mode = 'numeric') || const < 1e-300 || const < maxdens*1e-5){
        index_max <- max(which(dens_values == maxdens))
        # re-adjust the limits of the support (cut off region with very small density value)
        max_y_aux <- y_aux[min(c(length(y_aux), index_max + max(which(dens_values[index_max:length(dens_values)] > maxdens*1e-30)) + 2 ))]
        min_y_aux <- y_aux[max(c(1, min(which(dens_values[1:index_max] > maxdens*1e-30))))]
        const<-integrate(f = dMBMilstein_unnormalized, lower = min_y_aux, upper = max_y_aux,  x_k, x_b,
                         theta, delta_t, stepsToRight, scale = 1)$value
      }

      if (const == 0){
        warning("The function dMBMilstein_unnormalized() integrates to zero.")
        return(NA)
      }else{
        condition <- (y >= min_y) & (y <= max_y) # function values of the points outside of this interval will be zero
        result <- vector(mode = "numeric", length = length(y))
        if(any(condition)){
          result[condition] <- dMBMilstein_unnormalized(y[condition], x_k, x_b, theta,
                                                        delta_t, stepsToRight) / const
        }
        return(result)

      }
    }
  }
}


dMBMilstein_unnormalized <- function(y, x_k, x_b, theta, delta_t, left_index, right_index){
  if( any(c(y, x_k, x_b) <= 0)){
    stop("The function dMBMilstein() is only valid for positive input states y, x_k, x_b.\n",
         " y = ", y,
         ", x_k = ", x_k,
         ", x_b = ", x_b
    )
  }
  stepsToRight <- right_index - left_index - 1
  drift <- drift_fct(x_k, theta)
  diffusion <- diffusion_fct(x_k, theta)
  diffusionDeriv <- diffusion_fct_derivative(x_k, theta)

  ## determine the limits of the support of the density
  # simplified GBM-specific !!!!  limits for y
  # assumes positive state space => diffusion > 0
  min_y  <- max(c(0, x_k - 1/2 * diffusion / diffusionDeriv +
                    (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
  divisor <- 1/2 + (theta[1] - 1/2 * theta[2]) * delta_t * stepsToRight
  if(divisor >= 0){
    max_y <- x_b / divisor
  }else{ # divisor < 0
    min_y <- max(c(x_b / divisor, min_y))
    max_y <- Inf
  }

  if( min_y > max_y ){
    warning(paste("Set of feasible points y is empty for the given parameters:\n min_y = ",
                  min_y, ", max_y = ", max_y, sep=""))
    return(NA)
  }else if(isTRUE(all.equal(min_y, max_y))){
    warning(paste("The limits of the set of feasible points y are equal to each other:\n min_y = ",
                  min_y, ", max_y = ", max_y, sep=""))
    result <- vector(mode = "numeric", length = length(y))
    for(i in 1:length(y)){ if(isTRUE(all.equal(max_y, y[i]))){
      result[i] <- Inf
    }
    }
    return(result)
  }else{
    dMBMilstein_unnormalized_internal <- function(y, x_k, x_b, theta, delta_t, stepsToRight){

      A_squared <- diffusion^2 + 2 * diffusion * diffusionDeriv *
        (y - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t)
      # due to numercial issues, A_squared may become negative (but small) instead of zero
      cond1 <- A_squared < 0 & abs(A_squared) < 1e-11
      if(any(cond1)){
        A <- vector(mode = "numeric", length = length(y))
        A[cond1] <- 0
        if(any(!cond1)) A[!cond1] <- sqrt(A_squared[!cond1])
      }else{
        A <- sqrt(A_squared)
      }
      B <- - (diffusion + diffusionDeriv *
                (y - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
      C <- diffusion * diffusionDeriv^2 * delta_t
      drift_y <- drift_fct(y, theta)
      diffusion_y <- diffusion_fct(y, theta)
      diffusionDeriv_y <- diffusion_fct_derivative(y, theta)
      A_y_squared <- diffusion_y^2 + 2 * diffusion_y * diffusionDeriv_y *
        (x_b - y - (drift_y - 1/2 * diffusion_y * diffusionDeriv_y) *
           delta_t * stepsToRight)
      # due to numercial issues, A_y_squared may become negative (but small) instead of zero
      cond2 <- A_y_squared < 0 & abs(A_y_squared) < 1e-11
      if(any(cond2)){
        A_y <- vector(mode = "numeric", length = length(y))
        A_y[cond2] <- 0
        if(any(!cond2)) A_y[!cond2] <- sqrt(A_y_squared[!cond2])
      }else{
        A_y <- sqrt(A_y_squared)
      }
      B_y <- - (diffusion_y + diffusionDeriv_y *
                  (x_b - y - (drift_y - 1/2 * diffusion_y * diffusionDeriv_y) *
                     delta_t * stepsToRight))
      C_y <- diffusion_y * diffusionDeriv_y^2 * delta_t * stepsToRight
      return((exp((B - A)/C) + exp((B + A)/C)) / (sqrt(2 * pi * delta_t) * A) *
               (exp((B_y - A_y)/C_y) + exp((B_y + A_y)/C_y)) /
               (sqrt(2 * pi * delta_t * stepsToRight) * A_y))
    }

    # # calculate denominator (which is not the normalisation constant!)
    # A_squared_x_b <- diffusion^2 + 2 * diffusion * diffusionDeriv *
    #   (x_b - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t * (right_index - left_index))
    # # due to numercial issues, A_squared may become negative (but small) instead of zero
    # cond3 <- A_squared_x_b < 0 & abs(A_squared_x_b) < 1e-11
    # if(any(cond3)){
    #   A_x_b <- vector(mode = "numeric", length = length(y))
    #   A_x_b[cond3] <- 0
    #   if(any(!cond3)) A_x_b[!cond3] <- sqrt(A_squared_x_b[!cond3])
    # }else{
    #   A_x_b <- sqrt(A_squared_x_b)
    # }
    # B_x_b <- - (diffusion + diffusionDeriv *
    #           (x_b - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t * (right_index - left_index)))
    # C_x_b <- diffusion * diffusionDeriv^2 * delta_t * (right_index - left_index)

    condition <- (y >= min_y) & (y <= max_y) # function values of the points outside of this interval will be zero
    result <- vector(mode = "numeric", length = length(y))
    if(any(condition)){
      result[condition] <- dMBMilstein_unnormalized_internal(y[condition], x_k, x_b, theta,
                                                             delta_t, stepsToRight) #/
      # (exp((B_x_b - A_x_b)/C_x_b) + exp((B_x_b + A_x_b)/C_x_b)) *
      # (sqrt(2 * pi * delta_t *  (right_index - left_index)) * A_x_b)
    }
    return(result)
  }
}


rMBMilstein <- function(n, x_k, x_b, theta, delta_t, left_index, right_index,
                        precision = 10^4){
  if( any(c(x_k, x_b) <= 0)){ # GBM-specific
    stop("The function rMBMilstein() is only valid for positive input states x_k, x_b.\n",
         "x_k = ", x_k,
         ", x_b = ", x_b
    )
  }
  stepsToRight <- right_index - left_index - 1
  drift <- drift_fct(x_k, theta)
  diffusion <- diffusion_fct(x_k, theta)
  diffusionDeriv <- diffusion_fct_derivative(x_k, theta)

  ## determine the limits of the support of the density
  # simplified GBM-specific !!!!  limits for y
  # assumes positive state space => diffusion > 0
  min_y  <- max(c( 0, x_k - 1/2 * diffusion / diffusionDeriv +
                     (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
  divisor <- 1/2 + (theta[1] - 1/2 * theta[2]) * delta_t * stepsToRight
  if(divisor >= 0){
    max_y <- x_b / divisor
  }else{ # divisor < 0
    min_y <- max(c(x_b / divisor, min_y))
    max_y <- Inf
  }

  if( min_y > max_y ){
    warning(paste("Set of feasible points y is empty for the given parameters:\n min_y = ",
                  min_y, ", max_y = ", max_y, sep=""))
    return(NA)
  }else if(isTRUE(all.equal(min_y, max_y))){
    warning(paste("The limits of the set of feasible points y are equal to each other:\n min_y = ",
                  min_y, ", max_y = ", max_y, sep=""))
    return(rep(min_y, times = n))
  }else{
    dMBMilstein_unnormalized1 <- function(y, x_k, x_b, theta, delta_t, stepsToRight){
      A_squared <- diffusion^2 + 2 * diffusion * diffusionDeriv *
        (y - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t)
      # due to numercial issues, A_squared may become negative (but small) instead of zero
      cond1 <- A_squared < 0 & abs(A_squared) < 1e-11
      if(any(cond1)){
        A <- vector(mode = "numeric", length = length(y))
        A[cond1] <- 0
        if(any(!cond1)) A[!cond1] <- sqrt(A_squared[!cond1])
      }else{
        A <- sqrt(A_squared)
      }
      B <- - (diffusion + diffusionDeriv *
                (y - x_k - (drift - 1/2 * diffusion * diffusionDeriv) * delta_t))
      C <- diffusion * diffusionDeriv^2 * delta_t

      drift_y <- drift_fct(y, theta)
      diffusion_y <- diffusion_fct(y, theta)
      diffusionDeriv_y <- diffusion_fct_derivative(y, theta)
      # general implementation:
      # A_y_squared <- diffusion_y ^ 2 + 2 * diffusion_y * diffusionDeriv_y *
      #   (x_b - y - (drift_y - 1/2 * diffusion_y * diffusionDeriv_y) *
      #      delta_t * stepsToRight)
      # GBM specific (!) implementation to avoid numerical issues:
      A_y_squared <- y * 2 * theta[2] * x_b -
        y^2 * theta[2] *
        (1 + 2 * (theta[1] - 1/2 * theta[2]) * delta_t * stepsToRight)
      # due to numercial issues, A_y_squared may become negative (but small) instead of zero
      cond2 <- A_y_squared < 0 & abs(A_y_squared) < 1e-11
      if(any(cond2)){
        A_y <- vector(mode = "numeric", length = length(y))
        A_y[cond2] <- 0
        if(any(!cond2)) A_y[!cond2] <- sqrt(A_y_squared[!cond2])
      }else{
        A_y <- sqrt(A_y_squared)
      }
      B_y <- - (diffusion_y + diffusionDeriv_y *
                  (x_b - y - (drift_y - 1/2 * diffusion_y * diffusionDeriv_y) *
                     delta_t * stepsToRight))
      C_y <- diffusion_y * diffusionDeriv_y^2 * delta_t * stepsToRight
      return((exp((B - A)/C) + exp((B + A)/C)) / (sqrt(2 * pi * delta_t) * A) *
               (exp((B_y - A_y)/C_y) + exp((B_y + A_y)/C_y)) /
               (sqrt(2 * pi * delta_t * stepsToRight) * A_y))
    }
    # auxiliary variables to substitute Inf by very large number
    max_y_aux <- min(c(max_y, 1e16))
    if(is.infinite(max_y)){
      # use logarithmic scaling
      y <- min_y + exp(seq(from = 0, to =  log(max_y_aux - min_y + 1), by = log((max_y_aux - min_y)) / precision*10)) - 1
    }else{
      y <- seq(from = min_y, to = max_y_aux, by = (max_y_aux - min_y) / precision*10)
    }

    # leave out the interval limits, as the density may become infinite/NA there
    dens_values <- dMBMilstein_unnormalized1(y[2:(length(y)-1)], x_k, x_b, theta, delta_t,
                                             stepsToRight)
    maxdens <- max(dens_values, na.rm = TRUE)
    if (is.na(maxdens) | maxdens == 0){
      warning(paste("The maximum of the calculated density is equal to ",
                    maxdens, ",\nwherefore the density is not valid.", sep=""))
      return(NA)
    }else if ((maxdens <= 1e-320) & (var(dens_values, na.rm = TRUE) == 0)){
      warning(paste("The calculated density seems to be constant with a value of ",
                    maxdens, "\nwhich indicates that it is not correctly calculated due to numerical issues.", sep=""))
      return(NA)
    }else{
      index_max <- max(which(dens_values == maxdens))
      # adjust the  limits of the support (cut off region with very small density value)
      max_y_aux <- y[min(c(length(y), index_max + max(which(dens_values[index_max:length(dens_values)] > maxdens*1e-20)) + 1 ))]
      min_y_aux <- y[max(c(1, min(which(dens_values[1:index_max] > maxdens*1e-20)) - 1 ))]
      y <- seq(from = min_y_aux, to = max_y_aux, by = (max_y_aux - min_y_aux) / precision)
      # recalculate the maximum on the adjusted support to increase precision
      dens_values <- dMBMilstein_unnormalized1(y[2:(length(y)-1)], x_k, x_b, theta, delta_t,
                                               stepsToRight)
      maxdens <- max(dens_values)

      # use rejection sampling to draw from dMBMilstein
      naccepted <- 0
      rv <- vector(mode = "numeric", length = n)
      while (naccepted < n){
        y <- runif(1, min_y_aux, max_y_aux)
        dens_y <- runif(1, 0, maxdens)
        if (dens_y <
            (dMBMilstein_unnormalized1(y, x_k, x_b, theta, delta_t, stepsToRight))) {
          naccepted <- naccepted + 1
          rv[naccepted] <- y
        }
      }
      return(rv)
    }
  }
}
