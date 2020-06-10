# create output plots from one MCMC chain


library(coda)
suppressMessages(library(boa, quietly = TRUE)) # for highest probability density interval
suppressMessages(library(MCMCpack, quietly = TRUE)) # for the inverse gamma distribution

load("simulation_study/GBM_alpha_1_sigma_2_x0_100/output/GBM_alpha0_0_rho2_10_kappa0_2_nu0_2_gamma_alpha_0.5_gamma_sigma_0.5_M_20_m_2_path_MB_td_Milstein_pd_Milstein_1.data")
stanfit_object <- readRDS("simulation_study/GBM_alpha_1_sigma_2_x0_100/output/Stan/stanfit_objects/stanfit_object_M_20_path_1.rds")

stan_sample <- as.matrix(stanfit_object)

filetype <- "pdf" # "eps", "pdf"
save_plots <- TRUE


# number of iterations used as burn-in 
burnIn <- 5000
numIterations <- dim(results$theta)[2]


if(save_plots){
  if(filetype == "pdf"){
    pdf("figures_and_tables/Fig_6_MCMC_traceplots.pdf", width = 7, height = 5.1)
  }else{
    setEPS()
    postscript("figures_and_tables/Fig_6_MCMC_traceplots.eps", width = 7, height = 5.1)
  }
}


  par(mar=c(3, 4.6, .2, .2))
  layout(matrix(c(1, 2, 3), 3, 1), heights = c(1.7, 1.7, 1.7))
  plot(results$theta[1, ], type = "l", ylab = '', xlab = "", las = 1,
       xaxt = 'n', lwd = .8)
  #title("MCMC alpha")
  abline(inputs$true_theta[1], 0, col = 2)
  abline(mean(results$theta[1, (burnIn + 1) : numIterations]),
         0, col = "#3300FF", lwd = 1.5)
  abline(boa.hpd(results$theta[1, ][(burnIn + 1) : numIterations], 0.05)[1], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(boa.hpd(results$theta[1, ][(burnIn + 1) : numIterations], 0.05)[2], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(v = burnIn, col = 3)
  #axis(1, at=0:5*20000, labels=c(0, expression('20,000'), expression('40,000'), expression('60,000'), expression('80,000'), expression('100,000')), las = 1)
  axis(1, at=0:6*50000, labels=c(0, expression('50,000'), expression('100,000'), expression('150,000'), expression('200,000'), expression('250,000'), expression('300,000')), las = 1)
  # axis(1, at=0:5*20000, labels=c(0, expression(2%*%10^4), expression(4%*%10^4), expression(6%*%10^4), expression(8%*%10^4), expression(10^5)), las = 1)
  title(ylab = expression(alpha), cex.lab = 1.3, line = 3.2)


  # sigma_2
  plot(results$theta[2, ], type = "l", ylab = "", xlab = "", las = 1,
       xaxt = 'n', lwd = .8)
  #title("MCMC sigma^2")
  abline(inputs$true_theta[2], 0, col = 2)
  abline(mean(results$theta[2, (burnIn + 1) : numIterations]),
         0, col = "#3300FF", lwd = 1.5)
  abline(boa.hpd(results$theta[2, ][(burnIn + 1) : numIterations], 0.05)[1], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(boa.hpd(results$theta[2, ][(burnIn + 1) : numIterations], 0.05)[2], 0, col = "#3300FF", lty = 2, lwd = 1.5)
  abline(v = burnIn, col = 3)
  #axis(1, at=0:5*20000, labels=c(0, expression('20,000'), expression('40,000'), expression('60,000'), expression('80,000'), expression('100,000')), las = 1)
  axis(1, at=0:6*50000, labels=c(0, expression('50,000'), expression('100,000'), expression('150,000'), expression('200,000'), expression('250,000'), expression('300,000')), las = 1)
  title(ylab = expression(sigma^2), cex.lab = 1.3, line = 3.2)


  # log-posterior densities
  plot(results$log_posterior,
       type = 'l',
       ylab = '',
       xlab = "",
       las = 1,
       xaxt = 'n', lwd = .8)
  #axis(1, at=0:5*20000, labels=c(0, expression('20,000'), expression('40,000'), expression('60,000'), expression('80,000'), expression('100,000')), las = 1)
  axis(1, at=0:6*50000, labels=c(0, expression('50,000'), expression('100,000'), expression('150,000'), expression('200,000'), expression('250,000'), expression('300,000')), las = 1)
  title(xlab = "number of iterations", cex.lab = 1.1, line = 2)
  title(ylab = "log-posterior density", cex.lab = 1.3, line = 3.2)
  abline(v = burnIn, col = 3)
  #title("log-posterior density values")
if(save_plots) dev.off()

  
if(save_plots){
  if(filetype == "pdf"){
    pdf("figures_and_tables/Fig_7_posterior_densities.pdf", width = 10, height = 3)
  }else{setEPS()
    postscript("figures_and_tables/Fig_7_posterior_densities.eps", width = 10, height = 3)
  } 
}

  par(mar=c(3.3, 4.5, .2, .2))
  layout(matrix(c(1, 2), 1, 2), heights = c(4))

  col_median <- 2 # "#006600"
  plot(density(results$theta[1, ][(burnIn + 1) : numIterations]),
       main = "",
       ylab = '',
       xlab = '',
       las = 1,
       xlim = c(-3,7.5))
  abline(v = inputs$true_theta[1], col = 1, lty = 3)
  abline(v = mean(results$theta[1, (burnIn + 1) : numIterations]), col = "#3300FF")
  abline(v = median(results$theta[1, (burnIn + 1) : numIterations]), col = col_median)
  abline(v = mean(stan_sample[,1]), col = "#3300FF", lty = 5)
  abline(v = median(stan_sample[,1]), col = col_median, lty = 5)
  title(ylab = bquote(paste("Density of ", alpha,  sep = "")), cex.lab = 1, line = 2.7)
  title(xlab = expression(alpha), cex.lab = 1, line = 2.2)

  plot(density(results$theta[2, ][(burnIn + 1) : numIterations]),
       main = "",
       ylab = '',
       xlab = '',
       las = 1,
       xlim = c(0.85,7))
  legend(x = 'topright', legend = c("true parameter value", "mean estimate from approx. posterior", "median estimate from approx. posterior", "mean estimate from true posterior", "median estimate from true posterior"), col = c(1,"#3300FF",col_median,"#3300FF",col_median), lty = c(3,1,1,5,5), cex = .7)
  abline(v = inputs$true_theta[2], col = 1, lty = 3)
  abline(v = mean(results$theta[2, (burnIn + 1) : numIterations]), col = "#3300FF")
  abline(v = median(results$theta[2, (burnIn + 1) : numIterations]), col = col_median)
  abline(v = mean(stan_sample[,2]), col = "#3300FF", lty = 5)
  abline(v = median(stan_sample[,2]), col = col_median, lty = 5)
  title(ylab =  bquote(paste("Density of ", sigma^2,  sep = "")), cex.lab = 1, line = 2.2)
  title(xlab = expression(sigma^2), cex.lab = 1, line = 2.2)
if(save_plots) dev.off()
