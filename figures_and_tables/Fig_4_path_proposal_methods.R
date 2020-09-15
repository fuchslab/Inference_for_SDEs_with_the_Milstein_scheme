# CF13 p. 178-179
# simulation of the Geometric Brownian motion
#
#  which fulfills the SDE dX_t = alpha*X_t*dt + sigma*X_t*dB_t, x_0=x_0


source("main_functions/functions_for_parameter_estimation.R")
source("main_functions/GBM_problem_specific_parameter_and_functions.R")

filetype <- "pdf" # "eps", "pdf"

save_plots <- TRUE

seed <- 42
theta <- c(0.1, .1)
t_min <- 0
t_max <- 1
x_0 <- 10
x_m <- 25
m <- 10 # number of subintervals
n <- 15 # number of paths

yaxis <- c(5,30)
marg <- .03
xaxis <- c(t_min - marg, t_max + marg)

delta_t <- (t_max - t_min) / m
t <- seq(from = t_min, to = t_max, by = delta_t)

par(mfrow=c(1, 1))
wid <- 2.2
hei <- 2

lab <- c(expression(tau[i]),expression(  tau[i+1])) #axis tiks

###################################################################################################
# Euler
# left-conditioned proposal
set.seed(seed)
Y_euler <- matrix(NA, nrow=n, ncol = m+1) #each row is a path
Y_euler[,1] <- x_0
Y_euler[,m+1] <- x_m

for (k in 1:(m-1)) # time steps
{ Y_euler[,k+1] <- rEuler(n, Y_euler[,k], theta, delta_t)
}

if(save_plots){
  if(filetype == "pdf"){
    pdf("figures_and_tables/Fig_4_c_lc_Euler_proposal.pdf", width = wid, height = hei)
  }else{
    setEPS()
    postscript("figures_and_tables/Fig_4_c_lc_Euler_proposal.eps", width = wid, height = hei)
  }
}

par(mar=c(2, .2, .2, .2))
plot(t, Y_euler[1,], type='l', ylim = yaxis, xlim = xaxis,
     xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
for (j in 2:n) # paths
{ lines(t,Y_euler[j,])
}
lines(c(t_min,t_max), c(x_0,x_m), type = 'p', pch = 20, col = 2)
axis(1, at = c(0,1), labels = lab)
axis(1, at = (1:(m-1)) / m, labels = NA)

if(save_plots)  dev.off()

# Euler - Modified Bridge proposal -----------------------------------------------------------
set.seed(seed)
Y_MB <- matrix(NA, nrow = n, ncol = m+1) #each row is a path
Y_MB[,1] <- x_0
Y_MB[,m+1] <- x_m

for (k in 1:(m-1)) # time steps
{ Y_MB[,k+1] <- rMBEuler(n, Y_MB[,k], Y_MB[,m+1], theta, delta_t, k, m+1)
}

if(save_plots){
  if(filetype == "pdf"){
    pdf("figures_and_tables/Fig_4_e_MB_Euler_proposal.pdf", width = wid, height = hei)
  }else{
    setEPS()
    postscript("figures_and_tables/Fig_4_e_MB_Euler_proposal.eps", width = wid, height = hei)
  }
}
par(mar=c(2, .2, .2, .2))
plot(t, Y_MB[1,], type = 'l', ylim = yaxis, xlim = xaxis, xlab = '', ylab = '',
     yaxt = 'n', xaxt = 'n')
for (j in 2:n) # paths
{ lines(t,Y_MB[j,])
}

lines(c(t_min,t_max),c(x_0,x_m),type='p',pch=20, col=2)
axis(1, at=c(0,1), labels=lab)
axis(1, at=(1:(m-1))/m, labels=NA)

if(save_plots)  dev.off()


###################################################################################################
# Milstein
# left-conditioned proposal
set.seed(seed)
Y_Miltein <- matrix(NA, nrow=n, ncol = m+1) #each row is a path
Y_Miltein[,1] <- x_0
Y_Miltein[,m+1] <- x_m

for (j in 1:n) # paths
{ for (k in 1:(m-1)) # time steps
  { Y_Miltein[j,k+1] <- rMilstein(1, Y_Miltein[j,k], theta, delta_t)
  }
}

if(save_plots){
  if(filetype == "pdf"){
    pdf("figures_and_tables/Fig_4_d_lc_Milstein_proposal.pdf", width = wid, height = hei)
  }else{
    setEPS()
    postscript("figures_and_tables/Fig_4_d_lc_Milstein_proposal.eps", width = wid, height = hei)
  }
}

par(mar=c(2, .2, .2, .2))
plot(t, Y_Miltein[1,], type = 'l', ylim = yaxis, xlim = xaxis, xlab = '',
     ylab = '', yaxt = 'n', xaxt = 'n')
for (j in 2:n) # paths
{ lines(t,Y_Miltein[j,])
}
lines(c(t_min,t_max),c(x_0,x_m),type='p',pch=20, col=2)
axis(1, at=c(0,1), labels=lab)
axis(1, at=(1:(m-1))/m, labels=NA)

if(save_plots) dev.off()

# Milstein - Modified Bridge proposal -----------------------------------------------------------
set.seed(seed)
Y_MBMilstein <- matrix(NA, nrow=n, ncol = m+1) #each row is a path
Y_MBMilstein[,1] <- x_0
Y_MBMilstein[,m+1] <- x_m

for (j in 1:n) # paths
{ for (k in 1:(m-1)) # time steps
  { Y_MBMilstein[j,k+1] <- rMBMilstein(1, Y_MBMilstein[j,k], Y_MBMilstein[j,m+1], theta, delta_t, k, m+1, precision = 10^4)
    #rnorm(1, mean = Y_MBMilstein[j,k] + (Y_MBMilstein[j,m+1] - Y_MBMilstein[j,k])/(m+1-k), sd = sqrt((m-k)/(m+1-k)*sigma*delta_t))
  }
}

if(save_plots){
  if(filetype == "pdf"){
  pdf("figures_and_tables/Fig_4_f_MB_Milstein_proposal.pdf", width = wid, height = hei)
  }else{
  setEPS()
  postscript("figures_and_tables/Fig_4_f_MB_Milstein_proposal.eps", width = wid, height = hei)
  }
}

par(mar=c(2, .2, .2, .2))
plot(t, Y_MBMilstein[1,], type = 'l', ylim = yaxis, xlim = xaxis, xlab = '',
     ylab = '', yaxt = 'n', xaxt = 'n')
for (j in 2:n) # paths
{ lines(t,Y_MBMilstein[j,])
}

lines(c(t_min,t_max),c(x_0,x_m),type='p',pch=20, col=2)
axis(1, at=c(0,1), labels=lab)
axis(1, at=(1:(m-1))/m, labels=NA)

if(save_plots) dev.off()

# Diffusion Bridge Milstein proposal -----------------------------------------------------------
set.seed(seed)
Y_DBMilstein <- matrix(NA, nrow=n, ncol = m+1) #each row is a path
Y_DBMilstein[,1] <- x_0
Y_DBMilstein[,m+1] <- x_m

for (j in 1:n) # paths
{ for (k in 1:(m-1)) # time steps
{ Y_DBMilstein[j,k+1] <- rDBMilstein(1, Y_DBMilstein[j,k], Y_DBMilstein[j,m+1], theta, delta_t, k, m+1, precision = 10^4)
#rnorm(1, mean = Y_MBMilstein[j,k] + (Y_MBMilstein[j,m+1] - Y_MBMilstein[j,k])/(m+1-k), sd = sqrt((m-k)/(m+1-k)*sigma*delta_t))
}
}

if(save_plots){
  if(filetype == "pdf"){
    pdf("figures_and_tables/Fig_4_g_DB_Milstein_proposal.pdf", width = wid, height = hei)
  }else{
    setEPS()
    postscript("figures_and_tables/Fig_4_g_DB_Milstein_proposal.eps", width = wid, height = hei)
  }
}

par(mar=c(2, .2, .2, .2))
plot(t, Y_DBMilstein[1,], type = 'l', ylim = yaxis, xlim = xaxis, xlab = '',
     ylab = '', yaxt = 'n', xaxt = 'n')
for (j in 2:n) # paths
{ lines(t,Y_DBMilstein[j,])
}

lines(c(t_min,t_max),c(x_0,x_m),type='p',pch=20, col=2)
axis(1, at=c(0,1), labels=lab)
axis(1, at=(1:(m-1))/m, labels=NA)

if(save_plots) dev.off()
