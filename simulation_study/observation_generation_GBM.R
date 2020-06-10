# observation generation of a
#
# Geometric Brownian motion
#
#  which fulfills the SDE dX_t = alpha*X_t*dt + sigma*X_t*dB_t, X_0=x_0

theta <- matrix(c(c(1, 2), c(1, 4), c(1, 9)), ncol = 3)

set.seed(242)


for (k in 1:ncol(theta)){
  obs_seed <- sample(1:10^4,1)
  set.seed(obs_seed)
  N <- 100 # number of paths

  Time <- 1
  max_M <- 1000
  max_m <- 1
  x_0 <- 100

  delta_t <- Time / (max_M * max_m)
  tau <- seq(from = 0, to = Time, by = delta_t)
  true_theta <- theta[,k]

  Y_obs <- matrix(vector(mode = "numeric", length = length(tau)*N), nrow = length(tau))
  # simulate GBM using the fact that the multiplicative incremetns are log-normally distributed
  Y_obs[1,] <- x_0
  for ( i in 2:length(tau)){
    Y_obs[i,] <- rlnorm(length(Y_obs[i,]), mean = log( Y_obs[i-1,]) +
                         (true_theta[1] - 1/2 * true_theta[2]) * delta_t,
                       sd = sqrt(true_theta[2] * delta_t))
  }


  # seeds for the parameter estiamtion estimation
  estim_seeds <- sample(1:10^4,N)

  # save results
  output_folder <- paste('simulation_study/GBM_alpha_', true_theta[1], 
                         "_sigma_", true_theta[2],  "_x0_", x_0, sep='')
  # data
  if(!dir.exists(output_folder))
  {
    dir.create(output_folder)
  }
  save(tau, Y_obs, obs_seed, true_theta, estim_seeds,
       file = paste(output_folder, '/',  "GBM_obs", '.data', sep=''))

  pdf(paste(output_folder, '/',  "GBM_obs", '.pdf', sep='') )
  plot_color <- rainbow(N)
  plot(tau, Y_obs[,1], type = 'l', col = plot_color[1],  ylab = "Y_obs",
       ylim = c(min(Y_obs), max(Y_obs)))
  title(paste("GBM with alpha = ", true_theta[1], " and sigma^2 = ", true_theta[2],
              "\n x_0 = ", x_0, ", seed = ", obs_seed, sep = ''))
  for(l in 2:N){
    lines(tau, Y_obs[,l], col = plot_color[l])
  }

  plot_color <- rainbow(N)
  plot(tau, Y_obs[,1], type = 'l', col = plot_color[1], ylab = "Y_obs",
       ylim = c(min(Y_obs), max(Y_obs)), log = c("y"))
  title(paste("GBM with alpha = ", true_theta[1], " and sigma^2 = ", true_theta[2],
              "\n x_0 = ", x_0, ", seed = ", obs_seed, "\n log scale", sep = ''))
  for(l in 2:N){
    lines(tau, Y_obs[,l], col = plot_color[l])
  }
  r = dev.off()


  # save plots of individual_paths
  # data
  if(!dir.exists(paste(output_folder, "/individual_paths", sep='')))
  {
    dir.create(paste(output_folder, "/individual_paths", sep=''))
  }
  for( r in 1:N){
    pdf(paste(output_folder, "/individual_paths", '/',  "GBM_obs", r, '.pdf', sep='') )
    plot(tau, Y_obs[,r], type = 'l', col = plot_color[r],  ylab = "Y_obs",
         ylim = c(min(Y_obs), max(Y_obs)))
    title(paste("GBM with alpha = ", true_theta[1], " and sigma^2 = ", true_theta[2],
                "\n x_0 = ", x_0, ", seed = ", obs_seed, sep = ''))
    r = dev.off()
  }
}

