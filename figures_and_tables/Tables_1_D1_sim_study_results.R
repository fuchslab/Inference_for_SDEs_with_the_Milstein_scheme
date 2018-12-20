library(xtable)
library(SimDesign)

model <- "GBM" # GBM CIR

v_m <- c(5,2) # reversed order
output <- matrix(rep(NA, length = length(v_m)*8*8), ncol = 8)
rownames(output) <- vector(mode ="character", length = length(v_m)*8)

for(i in 1:length(v_m)){

  if(model == "GBM"){
    load(paste("simulation_study/aggregated_output/GBM_alpha_1_sigma^2_2_M_50_m_", v_m[i], "_nIter_1e+05.data", sep = ""))
    output[(0:7)*length(v_m)+i,1] <- apply(meanAlpha, 2, mean)
    output[(0:7)*length(v_m)+i,2] <- apply(meanAlpha, 2, sd)
    output[(0:7)*length(v_m)+i,3] <- apply(meanAlpha, 2, bias, parameter = true_theta[1])
    output[(0:7)*length(v_m)+i,4] <- apply(meanAlpha, 2, RMSE, parameter = true_theta[1])
    output[(0:7)*length(v_m)+i,5] <- apply(meanSigma, 2, mean)
    output[(0:7)*length(v_m)+i,6] <- apply(meanSigma, 2, sd)
    output[(0:7)*length(v_m)+i,7] <- apply(meanSigma, 2, bias, parameter = true_theta[2])
    output[(0:7)*length(v_m)+i,8] <- apply(meanSigma, 2, RMSE, parameter = true_theta[2])

    colnames(output) <- c("mean of alpha estimates", "s.d. of alpha estimates", "bias of alpha", "RMSE of alpha",
                          "mean of sigma estimates", "s.d. of sigma estimates", "bias of sigma", "RMSE of sigma")
  }else if(model == "CIR"){
    load(paste("simulation_study/aggregated_output/CIR_alpha_1_beta_1_sigma_0.25_M_50_m_", v_m[i], "_nIter_1e+05.data", sep = ""))
    output[(0:7)*length(v_m)+i,1] <- apply(meanBeta, 2, mean)
    output[(0:7)*length(v_m)+i,2] <- apply(meanBeta, 2, sd)
    output[(0:7)*length(v_m)+i,3] <- apply(meanBeta, 2, bias, parameter = true_theta[2])
    output[(0:7)*length(v_m)+i,4] <- apply(meanBeta, 2, RMSE, parameter = true_theta[2])
    output[(0:7)*length(v_m)+i,5] <- apply(meanSigma, 2, mean)
    output[(0:7)*length(v_m)+i,6] <- apply(meanSigma, 2, sd)
    output[(0:7)*length(v_m)+i,7] <- apply(meanSigma, 2, bias, parameter = true_theta[3])
    output[(0:7)*length(v_m)+i,8] <- apply(meanSigma, 2, RMSE, parameter = true_theta[3])

    colnames(output) <- c("mean of beta estimates", "s.d. of beta estimates", "bias of beta", "RMSE of beta",
                          "mean of sigma estimates", "s.d. of sigma estimates", "bias of sigma", "RMSE of sigma")
  }

  rownames(output)[(0:7)*length(v_m)+i] <- apply(matrix(colnames(ARparam),nrow=1), 1,
                                                 function(x){
                                                   x <- gsub("Euler", "E", x)
                                                   x <- gsub("Milstein", "M", x)
                                                   x <- gsub("n = 100", "", x)
                                                   x <- gsub("leftCondi", "lCond", x)
                                                   paste(x, " m=", v_m[i], sep = "")})
}

output <- output[c(2,1,6,5,4,3,8,7,10,9,14,13,12,11,16,15),]
output <- xtable(output, auto = TRUE)
digits(output) <- 3

print(output, booktabs = TRUE)
