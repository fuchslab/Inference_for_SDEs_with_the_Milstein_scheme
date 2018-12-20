# to assess the values of the Kullback-Leibler divergence caluculated between the
# densities estimated from the 100 estimates for each procedure,
# we conduct a small simulation study where we draw 1000 samples of size 100 from a
# normal distribution and calculate the KL divergence between each of the samples

set.seed(26950)
n <- 1000

simulated_estimates_norm <- matrix(rnorm(n * 100, mean = .9, sd = 1), ncol = n)
KLdiv_norm  <- rep(NA, length = n * (n-1))
densities1 <- apply(simulated_estimates_norm, 2, function(x){hist(x,breaks = seq(from=-3.7,to=5.4,by=.001), freq = TRUE)$counts})
densities1 <- densities1[rowSums(densities1)!=0,]
vec <- c()
for(l in 1:(dim(densities1)[1]-1)){
  if(any(densities1[l,]==0)){
    densities1[l+1,] <- densities1[l+1,] + densities1[l,]
    vec <- c(vec, l)
  }
}
densities1 <- densities1[-vec,] / 100
if(any(densities1[dim(densities1)[1],]==0)){
  densities1[dim(densities1)[1]-1,] <-
    densities1[dim(densities1)[1]-1,] +
    densities1[dim(densities1)[1],]
  densities1 <- densities1[-dim(densities1)[1],]
}

simulated_estimates_lnorm <- matrix(rlnorm(n * 100, meanlog = .65, sdlog = .2 ), ncol = n)
KLdiv_lnorm  <- rep(NA, length = n * (n-1))
densities2 <- apply(simulated_estimates_lnorm, 2, function(x){hist(x,breaks = seq(from=.5,to=5,by=.001), freq = TRUE)$counts})
densities2 <- densities2[rowSums(densities2)!=0,]
vec <- c()
for(l in 1:(dim(densities2)[1]-1)){
  if(any(densities2[l,]==0)){
    densities2[l+1,] <- densities2[l+1,] + densities2[l,]
    vec <- c(vec, l)
  }
}
densities2 <- densities2[-vec,] / 100
if(any(densities2[dim(densities2)[1],]==0)){
  densities2[dim(densities2)[1]-1,] <-
    densities2[dim(densities2)[1]-1,] +
    densities2[dim(densities2)[1],]
  densities2 <- densities2[-dim(densities2)[1],]
}

for(j in 2:n){
  for(k in 1:(j-1)){
    KLdiv_norm[(j-2)*(j-1) + (k-1)*2 + 1] <- entropy::KL.empirical(densities1[,j], densities1[,k])
    KLdiv_norm[(j-2)*(j-1) + (k-1)*2 + 2] <- entropy::KL.empirical(densities1[,k], densities1[,j])

    KLdiv_lnorm[(j-2)*(j-1) + (k-1)*2 + 1] <- entropy::KL.empirical(densities2[,j], densities2[,k])
    KLdiv_lnorm[(j-2)*(j-1) + (k-1)*2 + 2] <- entropy::KL.empirical(densities2[,k], densities2[,j])
  }
}

print(paste("mean of KL_div_norm: ", mean(KLdiv_norm), sep = ""))
print(paste("s.d. of KLdiv_norm: ", sd(KLdiv_norm), sep = ""))
print(paste("min of KL_div_norm: ", min(KLdiv_norm), sep = ""))
print(paste("max of KLdiv_norm: ", max(KLdiv_norm), sep = ""))
print(paste("mean of KL_div_lnorm: ", mean(KLdiv_lnorm), sep = ""))
print(paste("s.d. of KLdiv_lnorm: ", sd(KLdiv_lnorm), sep = ""))
print(paste("min of KL_div_lnorm: ", min(KLdiv_lnorm), sep = ""))
print(paste("max of KLdiv_lnorm: ", max(KLdiv_lnorm), sep = ""))

# results:
# [1] "mean of KL_div_norm: 0.132311540116348"
# [1] "s.d. of KLdiv_norm: 0.0573942615654195"
# [1] "min of KL_div_norm: 0.00951635669598365"
# [1] "max of KLdiv_norm: 0.847860532810437"
# [1] "mean of KL_div_lnorm: 0.142202285913364"
# [1] "s.d. of KLdiv_lnorm: 0.0581891205553384"
# [1] "min of KL_div_lnorm: 0.0101386400732847"
# [1] "max of KLdiv_lnorm: 0.562162966354978"




# load(paste("simulation_study/aggregated_output/GBM_alpha_1_sigma^2_2_M_50_m_5_nIter_1e+05.data", sep = ""))
# x <- seq(from = -3, to = 4, by = .1)
# plot(x, dnorm(x, mean=.9, sd=1), type='l', col =2)
# for(i in 1:8){
#   lines(density(meanAlpha[,i]), col=4)
# }
# x <- seq(from = 1, to = 4, by = .01)
# plot(x, dlnorm(x, meanlog = .65, sdlog = .2 ), type='l', col =2)
# for(i in 1:8){
#   lines(density(meanSigma[,i]), col=4)
# }
