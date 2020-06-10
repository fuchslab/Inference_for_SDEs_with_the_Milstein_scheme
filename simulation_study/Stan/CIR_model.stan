//
// This Stan program defines a model of the Cox–Ingersoll–Ross process
// described by the SDE dX_t = alpha*(beta - X_t)*dt + sigma*sqrt(X_t)*dB_t, X_0=x_0
// based on the solution process of this SDE.

functions {
  real log_modified_bessel_first_kind(real y, real z); 
  
  real log_transition_density_CIR(real alpha, real beta, real sigma2,
                                  real time_step, real x, real y ){
    real log_dens;
    real c;
    real u;
    real v;
    real order;
    real z;
    c = 2 * alpha / (sigma2 *(1 - exp(-alpha * time_step)));
    u = c * x * exp(-alpha * time_step);
    v = c * y;
    order = 2 * alpha * beta / sigma2 - 1;
    z = 2 * sqrt(u * v);
    log_dens = log(c) + order / 2 * log(v / u) - u - v + 
                log_modified_bessel_first_kind(order, z);
    
    return log_dens;
  }
}

data {
  int<lower=0> M; // number of points observed
  vector<lower=0>[M + 1] y_obs; //  observed states
  real time_step; // time step between observations
  real<lower=0> alpha;
}

parameters {
  real<lower=0> beta; 
  real<lower=0> sigma2;
}

model {
  // likelihood of the observed path
  for (i in 1:M){
    target += log_transition_density_CIR(alpha, beta, sigma2, time_step,
                                          y_obs[i], y_obs[i+1]);
  }

  // priors
  beta ~ inv_gamma(3, 3);
  sigma2 ~ inv_gamma(3, 4);
}
