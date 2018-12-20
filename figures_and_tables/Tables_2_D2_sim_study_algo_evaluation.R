library(xtable)


v_m <- c(5,2) # reversed order


output <- matrix(rep(NA, length = length(v_m)*8*8), ncol = 8)
rownames(output) <- vector(mode ="character", length = length(v_m)*8)

for(i in 1:length(v_m)){

  #load(paste("simulation_study/aggregated_output/GBM_alpha_1_sigma^2_2_M_50_m_", v_m[i], "_nIter_1e+05.data", sep = ""))

  load(paste("simulation_study/aggregated_output/CIR_alpha_1_beta_1_sigma_0.25_M_50_m_", v_m[i], "_nIter_1e+05.data", sep = ""))

  output[(0:7)*length(v_m)+i,1] <- apply(multivarESS, 2, mean, na.rm = T)
  output[(0:7)*length(v_m)+i,2] <- apply(multivarESS, 2, sd, na.rm = T)
  output[(0:7)*length(v_m)+i,3] <- apply(ARparam, 2, mean, na.rm = T)
  output[(0:7)*length(v_m)+i,4] <- apply(ARparam, 2, sd, na.rm = T)
  output[(0:7)*length(v_m)+i,5] <- apply(ARpath, 2, mean, na.rm = T)
  output[(0:7)*length(v_m)+i,6] <- apply(ARpath, 2, sd, na.rm = T)
  output[(0:7)*length(v_m)+i,7] <- apply(Duration, 2, mean, na.rm = T)
  output[(0:7)*length(v_m)+i,8] <- apply(Duration, 2, sd, na.rm = T)

  rownames(output)[(0:7)*length(v_m)+i] <- apply(matrix(colnames(ARparam),nrow=1), 1,  function(x) paste(x, " m=", v_m[i], sep = ""))


}


colnames(output) <- c( "mean of mESS", "s.d. of mESS",
                       "mean of acc. rate of path", "s.d. of acc. rate of path",
                       "mean of acc. rate of param", "s.d. of acc. rate of param",
                       "mean of comp. time in s", "s.d. of comp. time in s")

output <- output[c(15,16,11,12,13,14,9,10,7,8,3,4,5,6,1,2),]
output <- apply(output, 2, rev)

output <- round(output, 3)

output <- xtable(output, auto = TRUE)
digits(output) <- c(0,0,0,3,3,3,3,1,1)

print(output, booktabs = TRUE)
