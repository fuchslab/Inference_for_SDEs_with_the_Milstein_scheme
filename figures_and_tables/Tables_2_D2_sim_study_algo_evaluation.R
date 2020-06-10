library(xtable)
library(stringr)

# folder containing the file with obvervations
obsFolder <- "GBM_alpha_1_sigma_2_x0_100" #"GBM_alpha_1_sigma_2_x0_100" #   CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10  
M <- 20

v_m <- c(1,2,5) 

folder_name <- paste0("simulation_study/", obsFolder, "/aggregated_output/")
file_names <- list.files(path = folder_name)

N <- 100

aggregated_performance_matrix <- function(obsFolder, M, m){
  # loads the aggregated output of the different MCMC methods and stores them in multidim. arrays
  folder_name <- paste("simulation_study/", obsFolder, "/aggregated_output/", sep = "")
  # determine the relevant aggregated_output files
  all_file_names <- list.files(path = folder_name)
  fileName_fragment <-  paste("M_", M, "_m_", m, sep = "")
  file_names <- all_file_names %>% str_subset(fileName_fragment)
  
  num_methods <- length(file_names)
  method_names_long <- sub("^(.[^_]*_+){6}", "", word(file_names, 1, -2, sep = fixed(".data")))
  method_names <- gsub("Milstein", "M", gsub("Euler", "E", method_names_long))
  
  cnames_overall <- c("numIterations", "multivarESS", "Duration", "ARpath", 
                      "ARparam", "covariance", "nNegPointProposals", "nMBSwitchToEuler")
  ncol_overall <- length(cnames_overall)
  
  resultsOverall_all <- array(rep(NA, length = N * ncol_overall * num_methods),
                              dim = c(N, ncol_overall, num_methods),
                              dimnames = list(path_index = 1:N, 
                                              results = cnames_overall, 
                                              method = method_names))
  for (i in 1:num_methods){
    try(load(paste(folder_name, file_names[i] ,sep = "")))
    resultsOverall_all[ , , i] <- resultsOverall
  }
  
  measures <-  c("numIterations", "multivarESS", "ARparam", "ARpath")
  measure_array <- resultsOverall_all[ , measures, ]
  
  
  mat_mean_measures <- apply(measure_array, 2:3, mean, na.rm = TRUE)
  rownames(mat_mean_measures) <- paste0(rownames(mat_mean_measures), "_mean")
  mat_sd_measures <- apply(measure_array, 2:3, sd, na.rm = TRUE)
  mat_cv_measures <- apply(measure_array, 2:3, sd, na.rm = TRUE) / mat_mean_measures
  #rownames(mat_sd_measures) <- paste0(rownames(mat_sd_measures), "_sd")
  rownames(mat_cv_measures) <- paste0(rownames(mat_cv_measures), "_cv")
  
  mat_measures <- rbind(mat_mean_measures, mat_cv_measures)
  num_measures <- length(measures)
  new_order <- rep(c(0,num_measures), times = num_measures) + rep(1:num_measures, each = 2)
  if(m > 1){
    method_order <- c("MB_td_E_pd_E", "MB_td_M_pd_E",  "MB_td_M_pd_M", "DBM_td_M_pd_M")
    new_names <- c("MBE-E", "MBE-M", "MBM-M", "DBM-M")
  }else{
    method_order <- c("td_E",  "td_M")
    new_names <- c("Euler", "Milstein")
    mat_measures[c("ARpath_mean", "ARpath_cv"), ] <- NA
  }
  mat_measures <- t(mat_measures[new_order, method_order])
  rownames(mat_measures) <- paste0(new_names, "_m_", m)
  return(mat_measures)
}


performance_mat <- NULL

for(m in v_m){
  mat <- aggregated_performance_matrix(obsFolder = obsFolder, M = M, m = m)
  performance_mat <- rbind(performance_mat, mat)
}

performance_mat <- xtable(performance_mat, auto = TRUE, NA.string = "-")
digits(performance_mat) <- c(0, 0,2,0,2,3,2,3,2)

print.xtable(performance_mat, NA.string = "\\centercell{$-$}", booktabs = TRUE)
