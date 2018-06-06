
# input of the problem specific parameters and functions
#
# Geometric Brownian motion
#
#  which fulfills the SDE dX_t = alpha*X_t*dt + sigma*X_t*dB_t, X_0=x_0

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



N <- 10
Time <- 1
M <- 10
x_0 <- 100
seed <- 4242424
true_theta <- c(1, .25)

delta_t <- Time / M
tau <- seq(from = 0, to = Time, by = delta_t)

set.seed(seed)
Y_exact <- matrix(vector(mode = "numeric", length = length(tau)*N), nrow = length(tau))
# simulate GBM using the fact that the multiplicative incremetns are log-normally distributed
Y_exact[1,] <- x_0
for ( i in 2:length(tau)){
  Y_exact[i,] <- rlnorm(length(Y_exact[i,]), mean = log( Y_exact[i-1,]) +
                          (true_theta[1] - 1/2 * true_theta[2]) * delta_t,
                        sd = sqrt(true_theta[2] * delta_t))
}

set.seed(seed)
Y_Euler <- matrix(vector(mode = "numeric", length = length(tau)*N), nrow = length(tau))
# simulate GBM using the Euler-Maruyama scheme
Y_Euler[1,] <- x_0
for ( i in 2:length(tau)){
  Y_Euler[i,] <- rEuler(N, Y_Euler[i-1,], true_theta, delta_t)
}

set.seed(seed)
Y_Milstein <- matrix(vector(mode = "numeric", length = length(tau)*N), nrow = length(tau))
# simulate GBM using the Milstein scheme
Y_Milstein[1,] <- x_0
for ( i in 2:length(tau)){
  Y_Milstein[i,] <- rMilstein(N, Y_Milstein[i-1,], true_theta, delta_t)
}

# tit <- paste("GMB with alpha = ", true_theta[1], ", sigma^2 = ", true_theta[2],
#              ", left conditioned, \n x_0 = ", x_0, ", M = ", M, ", T = ", Time, sep = "")
#pdf("figures_and_tables/Fig_1_diffusion_approx.pdf", width = 12, height = 3)
setEPS()
postscript("figures_and_tables/Fig_1_diffusion_approx.eps", width = 12, height = 3)
par(mar=c(4, 4, 1, 1))
tit <- ""
plot(tau, Y_exact[,1], type = 'b', col = 2, lwd = 1, lty = 2, cex = 1, pch = 19,
     ylim = c(min(c(Y_exact, Y_Euler, Y_Milstein)), 500),
     main = tit,
     xlab = '',
     ylab = '',
     las = 1,
     xaxt = 'n')
mtext(expression("X"[t]), side=2, las=1, line=3, cex = 1.3)
mtext('time', side=1, line=2.5, cex = 1.3)

axis(1, at=(0:M)/M, labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

lines(tau, Y_Euler[,1], type = 'l', col = 1, lwd = 2)
lines(tau, Y_Milstein[,1], type = 'l', col = 4, lwd = 2)

for(l in 3:4){
  lines(tau, Y_exact[,l], type = 'b', col = 2, lwd = 1, lty = 2, cex = 1, pch = 19)
  lines(tau, Y_Euler[,l], type = 'l', col = 1, lwd = 2)
  lines(tau, Y_Milstein[,l], type = 'l', col = 4, lwd = 2)
}
legend("topleft", legend = c("points of trajectories of the solution",
                             "Euler approximation", "Milstein approximation"),
       bty = "n", lty = c(2,1,1), pch = c(19, NA ,NA), x.intersp=.5,
       col = c(2,1,4), cex = 1.2)
dev.off()
