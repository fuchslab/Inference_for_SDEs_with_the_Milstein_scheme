library(xtable)
library(stringr)

# folder containing the file with obvervations
obsFolder <- "GBM_alpha_1_sigma_2_x0_100" #"GBM_alpha_1_sigma_2_x0_100" #   CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10  
M <- 20

v_m <- c(1, 2,5) 

folder_name <- paste0("simulation_study/", obsFolder, "/aggregated_output/")
file_names <- list.files(path = folder_name)

N <- 100
model_type <- substr(obsFolder,1,3)
if(model_type == "GBM"){
  par_names <- c("alpha", "sigma2")
}else if(model_type == "CIR"){
  par_names <- c("beta", "sigma2")
}

aggregated_RMSE_matrix <- function(obsFolder, M, m){
  # loads the aggregated output of the different MCMC methods and stores them in multidim. arrays
  folder_name <- paste("simulation_study/", obsFolder, "/aggregated_output/", sep = "")
  # determine the relevant aggregated_output files
  all_file_names <- list.files(path = folder_name)
  fileName_fragment <-  paste("M_", M, "_m_", m, sep = "")
  file_names <- all_file_names %>% str_subset(fileName_fragment)
  
  num_methods <- length(file_names)
  method_names_long <- sub("^(.[^_]*_+){6}", "", word(file_names, 1, -2, sep = fixed(".data")))
  method_names <- gsub("Milstein", "M", gsub("Euler", "E", method_names_long))
  
  cnames_separate <- c("mean", "median", "variance", "hpd_low", "hpd_up", "ESS_coda", "mode")
  ncol_separate <- length(cnames_separate)
  cnames_overall <- c("numIterations", "multivarESS", "Duration", "ARpath", 
                      "ARparam", "covariance", "nNegPointProposals", "nMBSwitchToEuler")
  ncol_overall <- length(cnames_overall)
  
  results_separate_all <-
    array(rep(NA, length = N * ncol_separate * length(par_names) * num_methods),
          dim = c(N, ncol_separate, length(par_names), num_methods),
          dimnames = list(path_index = 1:N,
                          results = cnames_separate,
                          parameter = par_names,
                          method = method_names))
  
  resultsOverall_all <- array(rep(NA, length = N * ncol_overall * num_methods),
                              dim = c(N, ncol_overall, num_methods),
                              dimnames = list(path_index = 1:N, 
                                              results = cnames_overall, 
                                              method = method_names))
  for (i in 1:num_methods){
    try(load(paste(folder_name, file_names[i] ,sep = "")))
    results_separate_all[ , , , i] <- results_separate
    resultsOverall_all[ , , i] <- resultsOverall
  }
  
  # load Stan aggregated output
  try(load(paste(folder_name, "true_posterior_M_", M, ".data" ,sep = "")))
  # contains: results_separate and resultsOverall
  
  results <- c("mean", "median", "variance")
  num_res <- length(results)
  num_par <- length(par_names)
  mat_RMSE <- matrix(rep(NA, length = num_par * num_methods * num_res), 
                     ncol = num_par * num_res)
  colnames(mat_RMSE) <- apply(as.matrix(expand.grid(results,par_names)), 
                              1, paste0, collapse = "_")
  rownames(mat_RMSE) <- dimnames(results_separate_all)[[4]]
  
  for(k in 1:num_res){
    for(i in 1:num_par){
      # calculate discrepancy between statistic of sample from approximate posterior and
      # and statistic of sample from true posterior (Stan sample)
      res_mat <- results_separate_all[ , results[k], i, ]
      res_Stan <- as.vector(results_separate[, results[k], i])
      discrepancy <- sweep(res_mat, MARGIN = 1, res_Stan)
      mat_RMSE[ ,(i-1) * num_res + k] <- apply(discrepancy, 2, function(x) sqrt(mean(x^2, na.rm = TRUE)))
    }
  }
  
  if(m > 1){
    method_order <- c("MB_td_E_pd_E", "MB_td_M_pd_E",  "MB_td_M_pd_M", "DBM_td_M_pd_M")
    new_names <- c("MBE-E", "MBE-M", "MBM-M", "DBM-M")
  }else{
    method_order <- c("td_E",  "td_M")
    new_names <- c("Euler", "Milstein")
  }
  mat_RMSE <- mat_RMSE[method_order, ]
  rownames(mat_RMSE) <- paste0(new_names, "_m_", m)
  return(mat_RMSE)
}

mat_RMSE <- NULL

for(m in v_m){
  mat <- aggregated_RMSE_matrix(obsFolder = obsFolder, M = M, m = m)
  mat_RMSE <- rbind(mat_RMSE, mat)
}


mat_RMSE <- xtable(mat_RMSE, auto = TRUE)
digits(mat_RMSE) <- 4# 3

print(mat_RMSE, booktabs = TRUE)
