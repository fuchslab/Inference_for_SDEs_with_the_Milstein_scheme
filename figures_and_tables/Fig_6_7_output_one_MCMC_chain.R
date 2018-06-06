# create output plots from one MCMC chain


library(coda)
suppressMessages(library(boa, quietly = TRUE)) # for highest probability density interval
suppressMessages(library(MCMCpack, quietly = TRUE)) # for the inverse gamma distribution

#load("simulation_study/GBM_alpha_1_sigma_0.25/output/GBM_alpha0_0_rho2_10_kappa0_2_nu0_2_gamma_a_0.5_gamma_s_0.5_M_100_m_2_path_MB_td_Milstein_pd_Milstein_2.data")
load("simulation_study/GBM_alpha_1_sigma_2/output/GBM_alpha0_0_rho2_10_kappa0_2_nu0_2_gamma_a_0.5_gamma_s_0.5_M_50_m_2_path_MB_td_Milstein_pd_Milstein_1.data")

filetype <- "eps" # "eps", "pdf"

#jpeg("figures_and_tables/Fig_6_MCMC_traceplots.jpeg", width = 7, height = 5.1, units = 'in', res = 1200)

if(filetype == "pdf"){
  pdf("figures_and_tables/Fig_6_MCMC_traceplots.pdf", width = 7, height = 5.1)
}else{
  setEPS()
  postscript("figures_and_tables/Fig_6_MCMC_traceplots.eps", width = 7, height = 5.1)
}

  par(mar=c(3, 4.5, .2, .2))
  layout(matrix(c(1, 2, 3), 3, 1), heights = c(1.7, 1.7, 1.7))
  plot(results$theta[1, ], type = "l", ylab = '', xlab = "", las = 1,
       xaxt = 'n', lwd = .8)
  #title("MCMC alpha")
  abline(inputs$true_theta[1], 0, col = 2)
  abline(mean(results$theta[1, (floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations]),
         0, col = "#3300FF", lwd = 1.5)
  abline(boa.hpd(results$theta[1, ][(floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations], 0.05)[1], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(boa.hpd(results$theta[1, ][(floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations], 0.05)[2], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(v = floor(inputs$burnIn * inputs$numIterations), col = 3)
  axis(1, at=0:5*20000, labels=c(0, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), las = 1)
  title(ylab = expression(alpha), cex.lab = 1.3, line = 3.2)


  # sigma_2
  plot(results$theta[2, ], type = "l", ylab = "", xlab = "", las = 1,
       xaxt = 'n', lwd = .8)
  #title("MCMC sigma^2")
  abline(inputs$true_theta[2], 0, col = 2)
  abline(mean(results$theta[2, (floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations]),
         0, col = "#3300FF", lwd = 1.5)
  abline(boa.hpd(results$theta[2, ][(floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations], 0.05)[1], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(boa.hpd(results$theta[2, ][(floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations], 0.05)[2], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(v = floor(inputs$burnIn * inputs$numIterations), col = 3)
  axis(1, at=0:5*20000, labels=c(0, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), las = 1)
  title(ylab = expression(sigma^2), cex.lab = 1.3, line = 3.2)


  # log-posterior densities
  plot(results$log_posterior,
       type = 'l',
       ylab = '',
       xlab = "",
       las = 1,
       xaxt = 'n', lwd = .8)
  axis(1, at=0:5*20000, labels=c(0, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), las = 1)
  title(xlab = "number of iterations", cex.lab = 1.1, line = 2)
  title(ylab = "log-posterior density", cex.lab = 1.3, line = 3.2)
  abline(v = floor(inputs$burnIn * inputs$numIterations), col = 3)
  #title("log-posterior density values")
dev.off()

if(filetype == "pdf"){
  pdf("figures_and_tables/Fig_7_posterior_densities.pdf", width = 10, height = 3)
}else{setEPS()
  postscript("figures_and_tables/Fig_7_posterior_densities.eps", width = 10, height = 3)
}

  par(mar=c(3.3, 4.5, .2, .2))
  layout(matrix(c(1, 2), 1, 2), heights = c(4))

  logIncrements <- diff(log(inputs$Y_obs))
  dt <- diff(inputs$tau)

  objectiveFct <- function( theta){
    -sum(log(dnorm(logIncrements, mean = (theta[1] - 1/2 * theta[2]) * dt, sd = sqrt(theta[2] * dt))))
  }

  ML_theta <- optim(c(0,1), objectiveFct)$par

  # a posteriori density of the parameter
  objectiveFct2 <- function( theta){
    -sum(log(c(dnorm(logIncrements, mean = (theta[1] - 1/2 * theta[2]) * dt, sd = sqrt(theta[2] * dt)),
               dnorm(theta[1], mean = inputs$alpha_0, sd = inputs$rho_2 ),
               dinvgamma(theta[2], shape = inputs$kappa_0, scale = inputs$nu_0))))
  }

  MAP_theta <- optim(c(2,2), objectiveFct2, method = "BFGS")$par

  plot(density(results$theta[1, ][(floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations]),
       main = "",
       ylab = '',
       xlab = '',
       las = 1)
  abline(v = inputs$true_theta[1], col = 2)
  abline(v = ML_theta[1], col = "#006600")
  abline(v = MAP_theta[1], col = 3)
  abline(v = mean(results$theta[1, (floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations]), col = "#3300FF")
  #abline(v = inputs$alpha_0, col = 8)
  title(ylab = bquote(paste("Density of ", alpha,  sep = "")), cex.lab = 1, line = 2.7)
  title(xlab = expression(alpha), cex.lab = 1, line = 2.2)

  plot(density(results$theta[2, ][(floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations]),
       main = "",
       ylab = '',
       xlab = '',
       las = 1)
  #legend(x = 'topright', legend = c("true value", "MCMC mean estimate", "ML estimate", "MAP estimate", "prior mean"), col = c(2,3,5,6,8), lty = c(1,1,1,1,1), cex = .7)
  legend(x = 'topright', legend = c("true value", "MCMC mean estimate", "ML estimate", "MAP estimate"), col = c(2,"#3300FF","#006600",3), lty = c(1,1,1,1), cex = .7)
  abline(v = inputs$true_theta[2], col = 2)
  abline(v = mean(results$theta[2, (floor(inputs$burnIn * inputs$numIterations) +1) : inputs$numIterations]), col = "#3300FF")
  abline(v = ML_theta[2], col = "#006600")
  abline(v = MAP_theta[2], col = 3)
  #abline(v = inputs$nu_0 / (inputs$kappa_0 - 1), col = 8)
  title(ylab =  bquote(paste("Density of ", sigma^2,  sep = "")), cex.lab = 1, line = 2.2)
  title(xlab = expression(sigma^2), cex.lab = 1, line = 2.2)
dev.off()
