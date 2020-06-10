//
// This Stan program defines a model of the Geometric Brownian Motion
// which fulfills the SDE dX_t = alpha*X_t*dt + sigma*X_t*dB_t, X_0=x_0
// based on the solution process of this SDE.
// The multiplicative increments of the solution process are log-normally distributed.


data {
  int<lower=0> M; // number of points observed
  vector<lower=0>[M + 1] y_obs; // observed states
  real time_step; // time step between observations
}

transformed data{
  vector<lower=0>[M] mult_increments; // multiplicative increments of the data
  for (t in 1:M){
     mult_increments[t] = y_obs[t + 1] / y_obs[t];
  }
}

parameters {
  real alpha; // drift parameter 
  real<lower=0> sigma2; // diffusion parameter
}

model {
  // precalculate the parameters for the log-normal distribution
  real mean_log;
  real sd_log;
  mean_log = (alpha - sigma2 / 2) * time_step;
  sd_log = sqrt(sigma2 * time_step);
  
  // likelihood of the observed path            
  mult_increments ~ lognormal(mean_log, sd_log);

  // priors
  alpha ~ normal(0, 10);
  sigma2 ~ inv_gamma(2, 2);
}
