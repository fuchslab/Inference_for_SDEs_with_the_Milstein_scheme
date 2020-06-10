
filetype <- "pdf" # "eps", "pdf"
#
# Geometric Brownian motion
#
#  which fulfills the SDE dX_t = alpha*X_t*dt + sigma*X_t*dB_t, x_k=x_k

#-------------------------------------------------------------------------------------
# number of parameters (corresponds to the length of theta)
len_theta <- 2
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
# source the functions for dEuler(), dMilstein(), etc.
source("main_functions/functions_for_parameter_estimation.R")

Time <- 1
M <- 10
x_k <- 100
delta_t <- Time / M



#-------------------------------------------------------------------------------------
# parameter combination one
theta <- c(1, .25)

y <- seq(from = 50, to = 200, by = .1) #c(seq(from = 50, to = 60, by = .00001), seq(from = 61, to = 200, by = .1))
drift <- drift_fct(x_k, theta)
diffusion <- diffusion_fct(x_k, theta)
diffusionDeriv <- diffusion_fct_derivative(x_k, theta)
min_y <- x_k - 1/2*diffusion/diffusionDeriv +
  (drift - 1/2*diffusion*diffusionDeriv)*delta_t


dens_Euler <- dEuler(y, x_k, theta, delta_t)
dens_Milstein <- dMilstein(y, x_k, theta, delta_t)
dens_true <- dlnorm(y, mean = log( x_k) +
                      (theta[1] - 1/2 * theta[2]) * delta_t,
                    sd = sqrt(theta[2] * delta_t))

if(filetype == "pdf"){
  pdf("figures_and_tables/Fig_2a_trans_dens_alpha_1_sig2_025.pdf", width = 5.5, height = 2.5)
}else{
  setEPS()
  postscript("figures_and_tables/Fig_2a_trans_dens_alpha_1_sig2_025.eps", width = 5.5, height = 2.5)
}
par(mar=c(3, 4, .2, .2))
#ylim_max <- max(c(dens_Euler,dens_Milstein[!is.infinite(dens_Milstein)]))
#ylim_max <- 1.1*max(c(dens_Euler, dens_Milstein, dens_true))
plot(y, dens_Euler, type = 'l', col=1, lwd = 2,
     las = 1,
     yaxt = 'n',
     ylim = c(min(c(dens_Euler,dens_Milstein)), .026),
     ylab = '',
     xlab = '')
axis(2, at=(0:3)/100, labels=c(0,0.01,0.02,0.03), las = 1)

title(xlab = expression("Y"[k+1]), cex.lab = 1.2, line = 2)
title(ylab = bquote(paste(pi," ( Y"[k+1]," | ", "Y"[k], " = ", .(x_k), " )",  sep = "")), cex.lab = 1, line = 2.7)

lines(y, dens_Milstein, type = 'l', col=4, lwd = 2)
lines(y, dens_true, type = 'l', col=2, lwd = 1)
# legend("topright", legend = c("GBM solution", "Euler", "Milstein"), bty = "n", lty = c(1,1,1), lwd = c(1, 2 ,2),
#        col = c(2,1,4))
dev.off()

#-------------------------------------------------------------------------------------
# # parameter combination 2
# theta <- c(.1, 2)
#
# y <- seq(from = -5, to = 30, by = .01)
# drift <- drift_fct(x_k, theta)
# diffusion <- diffusion_fct(x_k, theta)
# diffusionDeriv <- diffusion_fct_derivative(x_k, theta)
# min_y <- x_k - 1/2*diffusion/diffusionDeriv +
#   (drift - 1/2*diffusion*diffusionDeriv)*delta_t
#
#
# dens_Euler <- dEuler(y, x_k, theta, delta_t)
# dens_Milstein <- dMilstein(y, x_k, theta, delta_t)
# dens_true <- dlnorm(y, mean = log( x_k) +
#                       (theta[1] - 1/2 * theta[2]) * delta_t,
#                     sd = sqrt(theta[2] * delta_t))
#
#
# pdf("trans_dens2.pdf", width = 5.5, height = 2.5)
# par(mar=c(3, 4, .2, .2))
# #ylim_max <- max(c(dens_Euler,dens_Milstein[!is.infinite(dens_Milstein)]))
# ylim_max <- .33
# plot(y, dens_Euler, type = 'l', col=1, lwd = 2,
#      las = 1,
#      yaxt = 'n',
#      ylim = c(min(c(dens_Euler,dens_Milstein)), ylim_max),
#      ylab = '',
#      xlab = '')
# axis(2,at=(0:3)/10, labels=c(0,0.1,0.2,0.3), las = 1)
#
# title(xlab = expression("Y"[k+1]), cex.lab = 1, line = 2)
# title(ylab = bquote(paste(pi," ( Y"[k+1]," | ", "Y"[k], " = ", .(x_k), " )",  sep = "")), cex.lab = 1, line = 2.3)
#
# lines(y[y<=min_y], dens_Milstein[y<=min_y], type = 'l', col=4, lwd = 2) # plot separately to avoid connection line from 0 to Inf
# lines(y[y>min_y], dens_Milstein[y>min_y], type = 'l', col=4, lwd = 2)
# lines(y, dens_true, type = 'l', col=2, lwd = 1)
# legend("topright", legend = c("GBM solution", "Euler", "Milstein"), bty = "n", lty = c(1,1,1), lwd = c(1, 2 ,2),
#        col = c(2,1,4))
# dev.off()

#-------------------------------------------------------------------------------------
# parameter combination 3
theta <- c(1, 2)

y <- seq(from = 0, to = 250, by = .1) #c(seq(from = 50, to = 60, by = .00001), seq(from = 61, to = 200, by = .1))
drift <- drift_fct(x_k, theta)
diffusion <- diffusion_fct(x_k, theta)
diffusionDeriv <- diffusion_fct_derivative(x_k, theta)
min_y <- x_k - 1/2*diffusion/diffusionDeriv +
  (drift - 1/2*diffusion*diffusionDeriv)*delta_t


dens_Euler <- dEuler(y, x_k, theta, delta_t)
dens_Milstein <- dMilstein(y, x_k, theta, delta_t)
dens_true <- dlnorm(y, mean = log( x_k) +
                      (theta[1] - 1/2 * theta[2]) * delta_t,
                    sd = sqrt(theta[2] * delta_t))

if(filetype == "pdf"){
  pdf("figures_and_tables/Fig_2b_trans_dens_alpha_1_sig2_2.pdf", width = 4.5, height = 2.5)
}else{
  setEPS()
  postscript("figures_and_tables/Fig_2b_trans_dens_alpha_1_sig2_2.eps", width = 4.5, height = 2.5)
}
par(mar=c(3, 1, .2, .2))
#ylim_max <- max(c(dens_Euler,dens_Milstein[!is.infinite(dens_Milstein)]))
ylim_max <- 1.1*max(c(dens_Euler, dens_Milstein, dens_true))
plot(y, dens_Euler, type = 'l', col=1, lwd = 2,
     las = 1,
     yaxt = 'n',
     ylim = c(min(c(dens_Euler,dens_Milstein)), .026),
     ylab = '',
     xlab = '')
#axis(2, at=(0:2)/100, labels=c(0,0.01,0.02), las = 1)
axis(2, at=(0:2)/100, labels=F, las = 1)

title(xlab = expression("Y"[k+1]), cex.lab = 1.2, line = 2)
#title(ylab = bquote(paste(pi," ( Y"[k+1]," | ", "Y"[k], " = ", .(x_k), " )",  sep = "")), cex.lab = 1, line = 2.7)

lines(y[y<=min_y], dens_Milstein[y<=min_y], type = 'l', col=4, lwd = 2) # plot separately to avoid connection line from 0 to Inf
lines(y[y>min_y], dens_Milstein[y>min_y], type = 'l', col=4, lwd = 2)
lines(y, dens_true, type = 'l', col=2, lwd = 1)
legend("topright", legend = c("GBM solution", "Euler", "Milstein"), bty = "n", lty = c(1,1,1), lwd = c(1, 2 ,2),
       col = c(2,1,4))
dev.off()
