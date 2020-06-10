filetype <- "eps" # "eps", "pdf"
save_plots <- TRUE


obsFolder <- "GBM_alpha_1_sigma_2_x0_100" #GBM_alpha_1_sigma_2_x0_100 CIR_alpha_1_beta_1_sigma_0.25_x0_3 CIR_alpha_1_beta_1_sigma_2_x0_10
M <- 20
N <- 100

num_plots <- 5
indices_2d_dens_plots <- (1:num_plots) * floor(N / num_plots)

model_type <- substr(obsFolder,1,3)
if(model_type == "GBM"){
  par_names <- c("alpha", "sigma2")
  par_labels <- c(expression(alpha), expression(sigma^2))
}else if(model_type == "CIR"){
  par_names <- c("beta", "sigma2")
  par_labels <- c(expression(beta), expression(sigma^2))
}

if(obsFolder == "GBM_alpha_1_sigma_2_x0_100"){
  true_values <- c(1,2)
}else if(obsFolder == "CIR_alpha_1_beta_1_sigma_0.25_x0_3" ){
  true_values <- c(1,0.25)
}else if(obsFolder == "CIR_alpha_1_beta_1_sigma_2_x0_10"){
  true_values <- c(1,2)
}

folderName <- paste('simulation_study/', obsFolder, "/output/Stan/", sep = "")

# 2d-density plots of the parameter sample from the true posterior -------------------

# Multiple plot function
#
# copied from : http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


plots <- list()
for(k in 1:length(indices_2d_dens_plots)){
  file_name <- paste("stanfit_objects/stanfit_object_", "M_", M, "_path_", indices_2d_dens_plots[k],
                     ".rds", sep = "")
  stanfit_object <- try(readRDS(paste(folderName, file_name, sep = "")))
  if(class(stanfit_object)!= "try-error"){
    res <- as.data.frame(stanfit_object)
    # Show the contour only
    p <- ggplot(res, aes_string(x=par_names[1], y=par_names[2])) +
      geom_density_2d()+
      theme(plot.margin = unit(c(0.2,0.1,0.1,0.3), "lines")) +
      ggtitle(paste("path no.", indices_2d_dens_plots[k])) + 
      theme(plot.title = element_text(size = 10, hjust = 0.5, vjust = .5),
            axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      labs(x = par_labels[1], y = par_labels[2])

    plots[[k]] <- p
  }
}

window_width <- 9
window_height <- 1.8

if(save_plots){
  if(model_type == "GBM"){
    file_name <- paste0("figures_and_tables/Fig_E1_2d_dens_plots_GBM_M_", M,
                        "_alpha_", true_values[1], "_sigma2_", true_values[2])
  }else if(model_type == "CIR"){
    file_name <- paste0("figures_and_tables/Fig_E3_2d_dens_plots_CIR_M_", M,
                        "_beta_", true_values[1], "_sigma2_", true_values[2]) 
  }
  if(filetype == "pdf"){
    pdf(paste0(file_name, ".pdf"), width = window_width, height = window_height)
  }else{
    setEPS()
    postscript(paste0(file_name, ".eps"), width = window_width, height = window_height)
  }
}

multiplot(plotlist = plots, cols = 5)

if(save_plots) dev.off()

# histograms of parameter correlations -----------------------------------------------
aggregated_correlations <- function(obsFolder, M, m){
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
  # load aggregated output for the different methods
  for (i in 1:num_methods){
    try(load(paste(folder_name, file_names[i] ,sep = "")))
    results_separate_all[ , , , i] <- results_separate
    resultsOverall_all[ , , i] <- resultsOverall
  }
  
  # load Stan aggregated output
  try(load(paste(folder_name, "true_posterior_M_", M, ".data" ,sep = "")))
  # contains: results_separate and resultsOverall
  
  
  correlation_mat <- matrix(rep(NA, N * (num_methods + 1)), ncol = num_methods + 1 )
  if(m > 1){
    method_order <- c("MB_td_E_pd_E", "MB_td_M_pd_E",  "MB_td_M_pd_M", "DBM_td_M_pd_M")
    new_names <- c("MBE-E", "MBE-M", "MBM-M", "DBM-M")
  }else{
    method_order <- c("td_E",  "td_M")
    new_names <- c("Euler", "Milstein")
  }
  colnames(correlation_mat) <- c("Stan", method_order)
  
  correlation_mat[,1] <- resultsOverall[, "covariance"] / 
    sqrt(results_separate[, "variance", 1] * results_separate[, "variance", 1] )
  
  for(k in method_names){
    correlation_mat[,k] <- resultsOverall_all[, "covariance", k] / 
      sqrt(results_separate_all[, "variance", 1, k] * results_separate_all[, "variance", 1, k] )
  }
  colnames(correlation_mat) <- c("True posterior", new_names)
  correlation_mat
}

m <- 5
correlation_mat <- aggregated_correlations(obsFolder, M, m)


window_width <- 7
window_height <- 2.8

if(save_plots){
  if(model_type == "GBM"){
    file_name <- paste0("figures_and_tables/Fig_E2_hist_corr_GBM_M_", M,
                        "_m_", m, "_alpha_", true_values[1], "_sigma2_", true_values[2])
  }else if(model_type == "CIR"){
    file_name <- paste0("figures_and_tables/Fig_E4_hist_corr_CIR_M_", M,
                        "_m_", m, "_beta_", true_values[1], "_sigma2_", true_values[2]) 
  }
  if(filetype == "pdf"){
    pdf(paste0(file_name, ".pdf"), width = window_width, height = window_height)
  }else{
    setEPS()
    postscript(paste0(file_name, ".eps"), width = window_width, height = window_height)
  }
}

par(mar=c(4, 4.5, 2, 1))
layout.matrix <- matrix(c(1, 2, 4,
                          0, 3, 5), nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = layout.matrix,
       heights = c( 2, 2), 
       widths = c(2, 2, 2)) 
if(model_type == "GBM"){
  xlim <- c(-0.04,0.4)
  breaks <- (-2:20)/50
  ylim <- c(0, 35)
}else if(model_type == "CIR"){
  xlim <- c(-0.05,0.6)
  breaks <- (-2.5:30)/50
  ylim <- c(0, 25)
}

for(k in 1:dim(correlation_mat)[2]){
  main <- dimnames(correlation_mat)[[2]][k]
  hist(correlation_mat[,k], 
       xlim = xlim, ylim = ylim, breaks = breaks,
       main = main, xlab = "",
       las = 1)
  title(xlab = "Parameter correlation", cex.lab = 1, line = 2)
}

if(save_plots) dev.off()