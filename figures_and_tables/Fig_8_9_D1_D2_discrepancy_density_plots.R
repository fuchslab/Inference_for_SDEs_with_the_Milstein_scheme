library(xtable)
library(stringr)


filetype <- "pdf" # "eps", "pdf"
save_plots <- TRUE

l_width <- 1.5

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
  par_labels <- c("alpha", expression(sigma^2))
}else if(model_type == "CIR"){
  par_names <- c("beta", "sigma2")
  par_labels <- c("beta", "sigma^2")
}

if(obsFolder == "GBM_alpha_1_sigma_2_x0_100"){
  true_values <- c(1,2)
}else if(obsFolder == "CIR_alpha_1_beta_1_sigma_0.25_x0_3" ){
  true_values <- c(1,0.25)
}else if(obsFolder == "CIR_alpha_1_beta_1_sigma_2_x0_10"){
  true_values <- c(1,2)
}

results <- c("mean", "median", "variance")
num_res <- length(results)

aggregated_discr_plots <- function(obsFolder, M, m, param_index, xlim){
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
  
  plot_discrep <- function(discrepancy, main, result, xlim){
    # estimate density for each method with the same bandwidth
    bandwidths <- apply(discrepancy, 2, function(x) density(x, na.rm = TRUE)$bw)
    bw <- mean(bandwidths)#mean(bandwidths) # or RMS: sqrt(sum(bandwidths ^ 2) / length(bandwidths))
    densities <- apply(discrepancy, 2, function(x) density(x, na.rm = TRUE, bw = bw))
    x_range <- range(unlist(lapply(densities, function(x) x[["x"]])))
    y_range <- range(unlist(lapply(densities, function(x) x[["y"]])))
    
    if(m > 1){
      v_col <- c(2,1,'#1b9e77',4)
      v_lty <- c(4,3,2, 1)
      v_l_width <- c(l_width, l_width, l_width, 1)
    }else{
      v_col <- c(1,4)
      v_lty <- c(3,1)
      v_l_width <- c(l_width, 1)
    }
    
    # plot the densities of the diffferent methods in the same window
    plot(NULL, xlim = xlim[[result]], ylim = y_range*1.1, main = '', xlab = '', las = 1, ann = FALSE)
    if(m == 1){
      title(main = main, cex.main = 1.4, line = 1.8, font = 2)
    }
    title(xlab = "Deviation", cex.lab = 1, line = 2.2)
    if(result == "mean"){
      title(ylab = paste0("m = ", m), cex.lab = 1.4, line = 3, font = 2, las = 1)
    }
    num_methods <- length(densities)

    for(l in 1:num_methods){
      lines(densities[[l]], col = v_col[l], lty = v_lty[l], xlim = xlim[[result]],
            lwd = v_l_width[l])
    }
    return(densities)
  }
  
  for(k in 1:num_res){
    # calculate discrepancy between statistic of sample from approximate posterior and
    # and statistic of sample from true posterior (Stan sample)
    res_mat <- results_separate_all[ , results[k], param_index, ]
    res_Stan <- as.vector(results_separate[, results[k], param_index])
    discrepancy <- sweep(res_mat, MARGIN = 1, res_Stan)
    if(param_index == 1){
      if(model_type == "GBM"){
        main <- bquote(Posterior ~ .(results[k]) ~ alpha)
      }else if(model_type == "CIR"){
        main <- bquote(Posterior ~ .(results[k]) ~ beta)
      }
    }else if (param_index == 2){
      main <- bquote(Posterior ~ .(results[k]) ~ sigma^2)
    }
    
    densities <- plot_discrep(discrepancy, 
                              main, 
                              results[k], xlim)
  }
  par(xpd=TRUE) # for legend outside plot region
  plot.new()
  if(m > 1){
    order <- c(2, 3, 4, 1)
    #leg_names <- names(densities)[order]
    leg_names <- c("MBE-E", "MBE-M", "MBM-M", "DBM-M")
    v_col <- c(2,1,'#1b9e77',4)[order]
    v_lty <- c(4,3,2, 1)[order]
    v_l_width <- c(l_width, l_width, l_width, l_width)[order]
  }else{
    leg_names <- c("Euler", "Milstein")
    v_col <- c(1,4)
    v_lty <- c(3,1)
    v_l_width <- c(l_width, l_width)
  }
  legend("topleft", legend = leg_names, 
         col = v_col, lty = v_lty,
         bty = "n", cex = 1, seg.len = 2, inset=c(-5,0), lwd = v_l_width)
  par(xpd=FALSE)
  #mtext(paste0("Based on up to ", sum(!(is.na(discrepancy[,1]))), " results."), side=1, cex = .7)
  

}


# save plots -------------------------------------------------------------------------
window_width <- 8.1
window_height <- 5.4

#------plots alpha-------------------------------------------------------------------
if(save_plots){
  if(model_type == "GBM"){
    file_name <- paste0("figures_and_tables/Fig_8_discrepancy_dens_plots_GBM_alpha_M_", M,
                        "_alpha_", true_values[1], "_sigma2_", true_values[2])
  }else if(model_type == "CIR"){
    file_name <- paste0("figures_and_tables/Fig_D1_discrepancy_dens_plots_CIR_alpha_M_", M,
                        "_beta_", true_values[1], "_sigma2_", true_values[2]) 
  }
  if(filetype == "pdf"){
    pdf(paste0(file_name, ".pdf"), width = window_width, height = window_height)
  }else{
    setEPS()
    postscript(paste0(file_name, ".eps"), width = window_width, height = window_height)
  }
}

if(model_type == "GBM"){
  xlim <- list(mean = c(-2,1), median = c(-2,1), variance = c(-3,1))
}else if(model_type == "CIR"){
  xlim <- list(mean = c(-.08,.08), median = c(-.06,.06), variance = c(-.2,.2))
}

par(mar=c(3.1, 4.5, 2.6, 1))
layout(matrix(1:((num_res + 1)*length(v_m)), ncol = num_res + 1, byrow = TRUE),
       widths = c(rep.int(1, num_res), .35))

for(i in 1:length(v_m)){
  aggregated_discr_plots(obsFolder, M, m = v_m[i], param_index = 1, xlim)
}

if(save_plots) dev.off()


#------plots sigma2-------------------------------------------------------------------
if(save_plots){
  if(model_type == "GBM"){
    file_name <- paste0("figures_and_tables/Fig_9_discrepancy_dens_plots_GBM_sigma2_M_", M,
                        "_alpha_", true_values[1], "_sigma2_", true_values[2])
  }else if(model_type == "CIR"){
    file_name <- paste0("figures_and_tables/Fig_D2_discrepancy_dens_plots_CIR_sigma2_M_", M,
                        "_beta_", true_values[1], "_sigma2_", true_values[2]) 
  }
  if(filetype == "pdf"){
    pdf(paste0(file_name, ".pdf"), width = window_width, height = window_height)
  }else{
    setEPS()
    postscript(paste0(file_name, ".eps"), width = window_width, height = window_height)
  }
}

if(model_type == "GBM"){
  xlim <- list(mean = c(-.7,.7), median = c(-.7,.7), variance = c(-0.4,.4))
}else if(model_type == "CIR"){
  xlim <- list(mean = c(-.3,.1), median = c(-.3,.1), variance = c(-.3,.15))
}


par(mar=c(3.1, 4.5, 2.6, 1))
layout(matrix(1:((num_res + 1)*length(v_m)), ncol = num_res + 1, byrow = TRUE),
       widths = c(rep.int(1, num_res), .35))

for(i in 1:length(v_m)){
  aggregated_discr_plots(obsFolder, M, m = v_m[i], param_index = 2, xlim)
}

if(save_plots) dev.off()

