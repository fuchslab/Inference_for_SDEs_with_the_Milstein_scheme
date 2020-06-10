# illustrate the used example paths of the
#
# Geometric Brownian motion
#
#  which fulfills the SDE dX_t = alpha*X_t*dt + sigma*X_t*dB_t, X_0=x_0

load("simulation_study/GBM_alpha_1_sigma_2_x0_100/GBM_obs.data")
# contains Y_obs, true_theta, tau

# selection of paths:
set.seed(1919495024)
path_indices <- sample.int(100, 10)

M <- 20
max_m <- 1
x_0 <- 100



index <- seq(from = 1, to = nrow(Y_obs), by = (nrow(Y_obs)-1)/M)
Y_obs <- Y_obs[index,]
tau <- tau[index]

Y_mean <- Y_obs[1,1] * exp(true_theta[1] * tau)
Y_min <- apply(Y_obs, 1, min)
Y_max <- apply(Y_obs, 1, max)





pdf("figures_and_tables/Fig_5_sample_paths.pdf", width = 10, height = 5)
#setEPS()
#postscript("figures_and_tables/Fig_5_sample_paths.eps", width = 10, height = 5)

# par(mar=c(4, 5, 1, 1))
# plot(tau, Y_mean, type = 'l',
#      ylim = c(min(Y_obs), max(Y_obs)),
#      main = "",
#      xlab = '',
#      ylab = '',
#      lwd = 4,
#      las = 1)
# mtext(expression("X"[t]), side=2, las=1, line=3.5, cex = 1.7)
# mtext("time", side=1, las=1, line=2, cex = 1.7)
# polygon(c(tau, rev(tau)), c(Y_min,rev(Y_max)), col = "gray", density = 50)
# lines(tau, Y_mean, lwd = 2)
# plot_color <- rainbow(length(path_indices)-1)
# for(l in 1:length(path_indices)){
#   lines(tau, Y_obs[,path_indices[l]], col = plot_color[l], lwd = 1.5)
# }

# in log scale
par(mar=c(4, 5, 1, 1))
plot(tau, log(Y_mean), type = 'l',
     ylim = c(log(min(Y_obs)), log(max(Y_obs))),
     main = "",
     xlab = '',
     ylab = '',
     yaxt = 'n',
     lwd = 4,
     las = 1)
mtext(expression("X"[t]), side=2, las=1, line=3.5, cex = 1.7)
mtext("time", side=1, las=1, line=2, cex = 1.7)
polygon(c(tau, rev(tau)), c(log(Y_min),rev(log(Y_max))), col = "gray", density = 50)
lines(tau, log(Y_mean), lwd = 2)
plot_color <- rainbow(length(path_indices)-1)
for(l in 1:length(path_indices)){
  lines(tau, log(Y_obs[,path_indices[l]]), col = plot_color[l], lwd = 1.8)
}

axis_lab <- c(5, 10, 50, 100, 500, 1000, 5000)
axis(2, at=log(axis_lab), labels = axis_lab, las = 1)

r = dev.off()



