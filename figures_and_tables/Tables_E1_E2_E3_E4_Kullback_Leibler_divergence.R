# to compare the estimation results for the different methods, we calculate the
# Kullback-Leibler divergence

library(xtable)
library(entropy)

model <- 'GBM' # GBM CIR
v_m <- c(2,5)

rename_rows <- function(mat){
  apply(matrix(colnames(mat),nrow=1), 1,
        function(x){
          x <- gsub("td_Euler", "tE", x)
          x <- gsub("pd_Euler", "pE", x)
          x <- gsub("td_Milstein", "tM", x)
          x <- gsub("pd_Milstein", "pM", x)
          x <- gsub("n = 100", "", x)
          x <- gsub("leftCondi", "LC", x)
          x <- gsub("\n", " ", x)
          x <- gsub(" ", "", x)
          paste(x, "m", v_m[i], sep = "")})
}
rename_columns <- function(mat){
  apply(matrix(colnames(mat),nrow=1), 1,
        function(x){
          x <- gsub("td_Euler", "tE", x)
          x <- gsub("pd_Euler", "pE", x)
          x <- gsub("td_Milstein", "tM", x)
          x <- gsub("pd_Milstein", "pM", x)
          x <- gsub("n = 100", "", x)
          x <- gsub("\n", " ", x)
          x <- gsub(" ", "", x)
          x <- gsub("leftCondi", "LC", x)})
}


#-- GBM ------------------------------------------------------------------------------
if(model=='GBM'){

  KLDmatrix_alpha <- matrix(rep(NA, length = length(v_m)*8*8), ncol = 8)
  rownames(KLDmatrix_alpha) <- rep(NA, length = length(v_m)*8)
  KLDmatrix_sigma <- matrix(rep(NA, length = length(v_m)*8*8), ncol = 8)
  rownames(KLDmatrix_sigma) <- rep(NA, length = length(v_m)*8)

  for(i in 1:length(v_m)){
    load(paste("simulation_study/aggregated_output/GBM_alpha_1_sigma^2_2_M_50_m_",
               v_m[i], "_nIter_1e+05.data", sep = ""))
    meanAlpha <- meanAlpha[,c(1,3,2,4,5,7,6,8)]
    meanSigma <- meanSigma[,c(1,3,2,4,5,7,6,8)]
    density_matrix_alpha <- apply(meanAlpha, 2,
                                  function(x){
                                    hist(x, breaks = seq(from=-2,to=3.2,by=.01),
                                         freq = TRUE)$counts})
    density_matrix_alpha <- density_matrix_alpha[rowSums(density_matrix_alpha)!=0,]
    vec <- c()
    for(l in 1:(dim(density_matrix_alpha)[1]-1)){
      if(any(density_matrix_alpha[l,]==0)){
        density_matrix_alpha[l+1,] <- density_matrix_alpha[l+1,] + density_matrix_alpha[l,]
        vec <- c(vec, l)
      }
    }
    density_matrix_alpha <- density_matrix_alpha[-vec,] / 100
    if(any(density_matrix_alpha[dim(density_matrix_alpha)[1],]==0)){
      density_matrix_alpha[dim(density_matrix_alpha)[1]-1,] <-
        density_matrix_alpha[dim(density_matrix_alpha)[1]-1,] +
        density_matrix_alpha[dim(density_matrix_alpha)[1],]
      density_matrix_alpha <- density_matrix_alpha[-dim(density_matrix_alpha)[1],]
    }

    density_matrix_sigma <- apply(meanSigma, 2,
                                  function(x){
                                    hist(x, breaks = seq(from=1,to=4,by=.01),
                                         freq = TRUE)$counts})
    density_matrix_sigma <- density_matrix_sigma[rowSums(density_matrix_sigma)!=0,]
    vec <- c()
    for(l in 1:(dim(density_matrix_sigma)[1]-1)){
      if(any(density_matrix_sigma[l,]==0)){
        density_matrix_sigma[l+1,] <- density_matrix_sigma[l+1,] + density_matrix_sigma[l,]
        vec <- c(vec, l)
      }
    }
    density_matrix_sigma <- density_matrix_sigma[-vec,] / 100
    if(any(density_matrix_sigma[dim(density_matrix_sigma)[1],]==0)){
      density_matrix_sigma[dim(density_matrix_sigma)[1]-1,] <-
        density_matrix_sigma[dim(density_matrix_sigma)[1]-1,] +
        density_matrix_sigma[dim(density_matrix_sigma)[1],]
      density_matrix_sigma <- density_matrix_sigma[-dim(density_matrix_sigma)[1],]
    }


    for(j in 1:8){
      for(k in (j):8){
        KLDmatrix_alpha[(i-1)*8 + j, k] <- entropy::KL.plugin(density_matrix_alpha[,j], density_matrix_alpha[,k])
        KLDmatrix_alpha[(i-1)*8 + k, j] <- entropy::KL.plugin(density_matrix_alpha[,k], density_matrix_alpha[,j])
        KLDmatrix_sigma[(i-1)*8 + j, k] <- entropy::KL.plugin(density_matrix_sigma[,j], density_matrix_sigma[,k])
        KLDmatrix_sigma[(i-1)*8 + k, j] <- entropy::KL.plugin(density_matrix_sigma[,k], density_matrix_sigma[,j])
      }
    }

    # rownames and columns names updated for each matrix seperatelxy to make sure, that
    # everything is in the right order
    rownames(KLDmatrix_alpha)[(1:8)+(i-1)*8] <- rename_rows(meanAlpha)
    rownames(KLDmatrix_sigma)[(1:8)+(i-1)*8] <- rename_rows(meanSigma)
  }
  colnames(KLDmatrix_alpha) <- rename_columns(meanAlpha)
  colnames(KLDmatrix_sigma) <- rename_columns(meanSigma)

  KLDmatrix_alpha[KLDmatrix_alpha==0] <- NA
  print(xtable(KLDmatrix_alpha ), digits=3, NA.string = "-")
  KLDmatrix_sigma[KLDmatrix_sigma==0] <- NA
  print(xtable(KLDmatrix_sigma ), digits=3, NA.string = "-")
}


#-- CIR ------------------------------------------------------------------------------
if(model=='CIR'){

  KLDmatrix_beta <- matrix(rep(NA, length = length(v_m)*8*8), ncol = 8)
  rownames(KLDmatrix_beta) <- rep(NA, length = length(v_m)*8)
  KLDmatrix_sigma <- matrix(rep(NA, length = length(v_m)*8*8), ncol = 8)
  rownames(KLDmatrix_sigma) <- rep(NA, length = length(v_m)*8)

  for(i in 1:length(v_m)){
    load(paste("simulation_study/aggregated_output/CIR_alpha_1_beta_1_sigma_0.25_M_50_m_",
               v_m[i], "_nIter_1e+05.data", sep = ""))
    meanBeta <- meanBeta[,c(1,3,2,4,5,7,6,8)]
    meanSigma <- meanSigma[,c(1,3,2,4,5,7,6,8)]
    density_matrix_beta <- apply(meanBeta, 2,
                                 function(x){
                                   hist(x,breaks = seq(from=.5,to=2.5,by=.01),
                                        freq = TRUE)$counts})
    density_matrix_beta <- density_matrix_beta[rowSums(density_matrix_beta)!=0,]
    vec <- c()
    for(l in 1:(dim(density_matrix_beta)[1]-1)){
      if(any(density_matrix_beta[l,]==0)){
        density_matrix_beta[l+1,] <- density_matrix_beta[l+1,] + density_matrix_beta[l,]
        vec <- c(vec, l)
      }
    }
    density_matrix_beta <- density_matrix_beta[-vec,] / 100
    if(any(density_matrix_beta[dim(density_matrix_beta)[1],]==0)){
      density_matrix_beta[dim(density_matrix_beta)[1]-1,] <-
        density_matrix_beta[dim(density_matrix_beta)[1]-1,] +
        density_matrix_beta[dim(density_matrix_beta)[1],]
      density_matrix_beta <- density_matrix_beta[-dim(density_matrix_beta)[1],]
    }

    density_matrix_sigma <- apply(meanSigma, 2,
                                  function(x){
                                    hist(x,breaks = seq(from=.1,to=1,by=.01),
                                         freq = TRUE)$counts})
    density_matrix_sigma <- density_matrix_sigma[rowSums(density_matrix_sigma)!=0,]
    vec <- c()
    for(l in 1:(dim(density_matrix_sigma)[1]-1)){
      if(any(density_matrix_sigma[l,]==0)){
        density_matrix_sigma[l+1,] <- density_matrix_sigma[l+1,] + density_matrix_sigma[l,]
        vec <- c(vec, l)
      }
    }
    density_matrix_sigma <- density_matrix_sigma[-vec,] / 100
    if(any(density_matrix_sigma[dim(density_matrix_sigma)[1],]==0)){
      density_matrix_sigma[dim(density_matrix_sigma)[1]-1,] <-
        density_matrix_sigma[dim(density_matrix_sigma)[1]-1,] +
        density_matrix_sigma[dim(density_matrix_sigma)[1],]
      density_matrix_sigma <- density_matrix_sigma[-dim(density_matrix_sigma)[1],]
    }


    for(j in 1:8){
      for(k in (j):8){
        KLDmatrix_beta[(i-1)*8 + j, k] <- entropy::KL.plugin(density_matrix_beta[,j], density_matrix_beta[,k])
        KLDmatrix_beta[(i-1)*8 + k, j] <- entropy::KL.plugin(density_matrix_beta[,k], density_matrix_beta[,j])
        KLDmatrix_sigma[(i-1)*8 + j, k] <- entropy::KL.plugin(density_matrix_sigma[,j], density_matrix_sigma[,k])
        KLDmatrix_sigma[(i-1)*8 + k, j] <- entropy::KL.plugin(density_matrix_sigma[,k], density_matrix_sigma[,j])
      }
    }
    # rownames and columns names updated for each matrix seperatelxy to make sure, that
    # everything is in the right order
    rownames(KLDmatrix_beta)[(1:8)+(i-1)*8] <- rename_rows(meanBeta)
    rownames(KLDmatrix_sigma)[(1:8)+(i-1)*8] <- rename_rows(meanSigma)
  }
  colnames(KLDmatrix_beta) <- rename_columns(meanBeta)
  colnames(KLDmatrix_sigma) <- rename_columns(meanSigma)

  KLDmatrix_beta[KLDmatrix_beta==0] <- NA
  print(xtable(KLDmatrix_beta ), digits=3, NA.string = "-")
  KLDmatrix_sigma[KLDmatrix_sigma==0] <- NA
  print(xtable(KLDmatrix_sigma ), digits=3, NA.string = "-")
}
