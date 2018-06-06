# observation generation of a
#
# Cox-Ingersoll-Ross process
#
#  which fulfills the SDE dX_t = alpha*(beta - X_t)*dt + sigma*sqrt(X_t)*dB_t, X_0=x_0

theta <- matrix(c(c(1, 1, 2), c(1, 1, 0.25), c(1, 5, 2)), ncol = 3)

set.seed(242)


for (k in 1:ncol(theta)){
  obs_seed <- sample(1:10^4,1)
  set.seed(obs_seed)
  N <- 100 # number of paths

  Time <- 1
  max_M <- 1e6
  max_m <- 1
  x_0 <- 3

  delta_t <- Time / (max_M * max_m)
  tau <- seq(from = 0, to = Time, by = delta_t)
  true_theta <- theta[,k]

  Y_obs <- matrix(vector(mode = "numeric", length = length(tau)*N), nrow = length(tau))
  # simulate CIR process using Euler scheme with very small time steps
  Y_obs[1,] <- x_0
  for ( i in 2:length(tau)){
    Y_obs[i,] <- rnorm(N, mean =  Y_obs[i-1,] +
                         true_theta[1] * (true_theta[2] - Y_obs[i-1,]) * delta_t,
                       sd = sqrt(true_theta[3] * Y_obs[i-1,] * delta_t))
  }



  M <- 1000
  index <- seq(from = 1, to = nrow(Y_obs), by = (nrow(Y_obs)-1)/M)
  Y_obs <- Y_obs[index,]
  tau <- tau[index]

  # seeds for the parameter estiamtion estimation
  estim_seeds <- sample(1:10^4,N)

  # save results
  output_folder <- paste('simulation_study/CIR_alpha_', true_theta[1], "_beta_", true_theta[2], "_sigma_", true_theta[3], sep='')
  # data
  if(!dir.exists(output_folder))
  {
    dir.create(output_folder)
  }
  save(tau, Y_obs, obs_seed, true_theta, estim_seeds,
       file = paste(output_folder, '/',  "CIR_obs", '.data', sep=''))

  pdf(paste( output_folder, '/',  "CIR_obs", '.pdf', sep='') )
  plot_color <- rainbow(N)
  plot(tau, Y_obs[,1], type = 'l', col = plot_color[1], ylim = c(0, max(Y_obs)),
       ylab = "Y_obs")
  title(paste("CIR with alpha = ", true_theta[1], ", beta = ", true_theta[2], ", sigma = ", true_theta[3],
              "\n x_0 = ", x_0, ", seed = ", obs_seed, sep = ''))
  for(l in 2:N){
    lines(tau, Y_obs[,l], col = plot_color[l])
  }
  r = dev.off()


  # data
  if(!dir.exists(paste(output_folder, "/individual paths", sep='')))
  {
    dir.create(paste(output_folder, "/individual paths", sep=''))
  }
  for( r in 1:N){
    pdf(paste(output_folder, "/individual paths", '/',  "CIR_obs", r, '.pdf', sep='') )
    plot(tau, Y_obs[,r], type = 'l', col = plot_color[r], ylim = c(0, max(Y_obs)),
         ylab = "Y_obs")
    title(paste("CIR with alpha = ", true_theta[1], ", beta = ", true_theta[2], ", sigma = ", true_theta[3],
                "\n x_0 = ", x_0, ", seed = ", obs_seed, sep = ''))
    r = dev.off()
  }
}

