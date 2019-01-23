

library(MCMCpack, quietly = TRUE) # for the inverse gamma distribution
library(reshape2)
library(ggplot2)
library(plyr)
library(grid)
library(gridExtra)
library(gtable)

filetype <- "eps" # "eps", "pdf"

m1 <- m2 <- m5 <- list()

# M <- 100
# nIter <- 1e+05
# theta <- c(1, .25)
# ylim_lower_alpha <- -.5
# ylim_upper_alpha <- 2.25
# ylim_lower_sigma <- 1.75
# ylim_upper_sigma <- 3.75

M <- 50
nIter <- 1e+05
theta <- c(1, 2)
ylim_lower_alpha <- -2
ylim_upper_alpha <- 3.5
ylim_lower_sigma <- 1
ylim_upper_sigma <- 4.2

# M <- 25
# nIter <- 1e+05
# theta <- c(1, 2)
# ylim_lower_alpha <- -7.25
# ylim_upper_alpha <- 3.75
# ylim_lower_sigma <- 0
# ylim_upper_sigma <- 10

# m = 1
load(paste("simulation_study/aggregated_output/GBM_alpha_", theta[1],"_sigma^2_", theta[2],"_M_", M, "_m_1_nIter_", nIter, ".data", sep = ""))
m1$meanAlpha <- meanAlpha
m1$modeAlpha <- modeAlpha
m1$meanSigma <- meanSigma
m1$modeSigma <- modeSigma

# m = 2
load(paste("simulation_study/aggregated_output/GBM_alpha_", theta[1],"_sigma^2_", theta[2],"_M_", M, "_m_2_nIter_", nIter, ".data", sep = ""))
m2$meanAlpha <- meanAlpha
m2$modeAlpha <- modeAlpha
m2$meanSigma <- meanSigma
m2$modeSigma <- modeSigma

# m = 5
load(paste("simulation_study/aggregated_output/GBM_alpha_", theta[1],"_sigma^2_", theta[2],"_M_", M, "_m_5_nIter_", nIter, ".data", sep = ""))
m5$meanAlpha <- meanAlpha
m5$modeAlpha <- modeAlpha
m5$meanSigma <- meanSigma
m5$modeSigma <- modeSigma


alpha <- true_theta[1]
sigma_2 <- true_theta[2]


line_width <- .8
num_left_lines <- 6

window_height <- 16
window_width <- 12
window_ratio <- window_height / window_width

fig_height <- 4


gg_violin <- function(df, method = "all", ylim_lower, ylim_upper){
  if(method=="all"){
    names(df) <- c("", "method", "alpha_estimate")
    df$method <- factor(df$method,
                            levels = c("leftCondi\ntd_Euler\npd_Euler\nn = 100", "leftCondi\ntd_Milstein\npd_Euler\nn = 100",
                                       "leftCondi\ntd_Euler\npd_Milstein\nn = 100", "leftCondi\ntd_Milstein\npd_Milstein\nn = 100",
                                       "MB\ntd_Euler\npd_Euler\nn = 100", "MB\ntd_Milstein\npd_Euler\nn = 100",
                                       "MB\ntd_Euler\npd_Milstein\nn = 100", "MB\ntd_Milstein\npd_Milstein\nn = 100",
                                       "ML_\nestimate\nn =  100", "MAP_\nestimate\nn =  100"), ordered = TRUE)
    # Basic violin plot
    g <- ggplot(df, aes(x=method, y=alpha_estimate)) +
      geom_violin(trim=TRUE, scale = "width") +
      geom_violin(data = df[df$method == "leftCondi\ntd_Euler\npd_Euler\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#990000", lwd = line_width) +
      geom_violin(data = df[df$method == "leftCondi\ntd_Milstein\npd_Euler\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#660099", lwd = line_width) +
      geom_violin(data = df[df$method == "leftCondi\ntd_Euler\npd_Milstein\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#660099", lwd = line_width) +
      geom_violin(data = df[df$method == "leftCondi\ntd_Milstein\npd_Milstein\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#3300FF", lwd = line_width) +
      geom_violin(data = df[df$method == "MB\ntd_Euler\npd_Euler\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#990000", lwd = line_width) +
      geom_violin(data = df[df$method == "MB\ntd_Euler\npd_Milstein\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#660099", lwd = line_width) +
      geom_violin(data = df[df$method == "MB\ntd_Milstein\npd_Euler\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#660099", lwd = line_width) +
      geom_violin(data = df[df$method == "MB\ntd_Milstein\npd_Milstein\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#3300FF", lwd = line_width) +
      geom_violin(data = df[df$method == "ML_\nestimate\nn =  100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#006600", lwd = line_width) +
      geom_violin(data = df[df$method == "MAP_\nestimate\nn =  100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#006600", lwd = line_width) +
      geom_boxplot(width=0.2) +
      coord_cartesian(ylim = c(ylim_lower, ylim_upper))
  }else if(method == "m1"){

    names(df) <- c("", "method", "alpha_estimate")
    df$method <- factor(df$method,
                            levels = c("td_Euler\nn = 100", "td_Euler\nn = 0" ,
                                       "td_Milstein\nn = 0", "td_Milstein\nn = 100" ,
                                       "td_Euler\nn2 = 100", "td_Euler\nn2 = 0" ,
                                       "td_Milstein\nn2 = 100" ,"td_Milstein\nn2 = 0",
                                       "ML_\nestimate\nn =  100", "MAP_\nestimate\nn =  100"),
                            ordered = TRUE)
    # Basic violin plot
    g <- ggplot(df, aes(x=method, y=alpha_estimate)) +
      geom_violin(trim=TRUE, scale = "width") +
      geom_violin(data = df[df$method == "td_Euler\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#990000", lwd = line_width) +
      geom_violin(data = df[df$method == "td_Milstein\nn = 100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#3300FF", lwd = line_width) +
      geom_violin(data = df[df$method == "ML_\nestimate\nn =  100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#006600", lwd = line_width) +
      geom_violin(data = df[df$method == "MAP_\nestimate\nn =  100",], aes(x=method, y=alpha_estimate), trim=TRUE, colour = "#006600", lwd = line_width) +
      geom_boxplot(width=0.2) +
      coord_cartesian(ylim = c(ylim_lower, ylim_upper))
  }
}


Brack <- function(x1,y1,x2,y2,h,orientation = "horizontal", ratio = window_ratio){
  if(orientation == "horizontal"){
    x2 <- x2-x1; y2 <- y2-y1
    v1 <- viewport(x=x1,y=y1,width=sqrt(x2^2+y2^2),
                   height=h,angle=180*atan2(y2,x2)/pi,
                   just=c("left","bottom"),gp=gpar(col="black"))
    #browser()
    pushViewport(v1)
    grid.curve(x2=0,y2=0,x1=.125,y1=.5,curvature=.5)
    grid.move.to(.125,.5)
    grid.line.to(.375,.5)
    grid.curve(x1=.375,y1=.5,x2=.5,y2=1,curvature=.5)
    grid.curve(x2=1,y2=0,x1=.875,y1=.5,curvature=-.5)
    grid.move.to(.875,.5)
    grid.line.to(.625,.5)
    grid.curve(x2=.625,y2=.5,x1=.5,y1=1,curvature=.5)
    popViewport()
  }else{
    x_2 <- y2-y1
    y2 <- x2-x1
    x2 <- x_2 * ratio
    v1 <- viewport(x=x1,y=y1,width=sqrt(x2^2+y2^2),
                   height=h,angle=180*atan2(y2,x2)/pi-90,
                   just=c("right","bottom"),gp=gpar(col="black"))

    pushViewport(v1)
    grid.curve(x2=0,y2=0,x1=.125,y1=.5,curvature=.5)
    grid.move.to(.125,.5)
    grid.line.to(.375,.5)
    grid.curve(x1=.375,y1=.5,x2=.5,y2=1,curvature=.5)
    grid.curve(x2=1,y2=0,x1=.875,y1=.5,curvature=-.5)
    grid.move.to(.875,.5)
    grid.line.to(.625,.5)
    grid.curve(x2=.625,y2=.5,x1=.5,y1=1,curvature=.5)
    popViewport()
  }
}



##------plots alpha-------------------------------------------------------------------
if(filetype == "pdf"){
  pdf(paste("figures_and_tables/Fig_8_violin_plots_alpha_M_", M, "_alpha_", theta[1], "_sigma_", theta[2], "_nIter_", nIter, ".pdf",
            sep = ""), width = window_width, height = window_height)
}else{
  setEPS()
  postscript(paste("figures_and_tables/Fig_8_violin_plots_alpha_M_", M, "_alpha_", theta[1], "_sigma_", theta[2], "_nIter_", nIter, ".eps",
                   sep = ""), width = window_width, height = window_height)

}

# alpha estimate from the mean of the MCMC posteriori-----------------------
# m = 1
colnames(m1$meanAlpha) <- c("td_Euler\nn = 100", "td_Euler\nn = 0" ,
                            "td_Milstein\nn = 100" ,"td_Milstein\nn = 0",
                            "td_Euler\nn2 = 100", "td_Euler\nn2 = 0" ,
                            "td_Milstein\nn2 = 100" ,"td_Milstein\nn2 = 0")
rAlpha <- melt(cbind(m1$meanAlpha, ML_alpha, MAP_alpha))
a <- gg_violin(rAlpha, method = "m1", ylim_lower_alpha, ylim_upper_alpha)  +
  labs(title= "",
       x="", y = bquote(paste( alpha, " estimated by mean", sep = ""))) +
  geom_hline(yintercept = alpha, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 2
rAlpha <- na.omit(melt(cbind(m2$meanAlpha, ML_alpha, MAP_alpha)))
b <- gg_violin(rAlpha, method = "all", ylim_lower_alpha, ylim_upper_alpha) +
  labs(title= "",
       x="", y = bquote(paste( alpha, " estimated by mean", sep = ""))) +
  geom_hline(yintercept = alpha, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 5
rAlpha <- na.omit(melt(cbind(m5$meanAlpha, ML_alpha, MAP_alpha)))
c <- gg_violin(rAlpha, method = "all", ylim_lower_alpha, ylim_upper_alpha) +
  labs(title= "",
       x="", y = bquote(paste( alpha, " estimated by mean", sep = ""))) +
  geom_hline(yintercept = alpha, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# alpha estimate from the mode of the MCMC posteriori-------------------
# m = 1
colnames(m1$modeAlpha) <- c("td_Euler\nn = 100", "td_Euler\nn = 0" ,
                            "td_Milstein\nn = 100" ,"td_Milstein\nn = 0",
                            "td_Euler\nn2 = 100", "td_Euler\nn2 = 0" ,
                            "td_Milstein\nn2 = 100" ,"td_Milstein\nn2 = 0")
rAlpha <- melt(cbind(m1$modeAlpha, ML_alpha, MAP_alpha))
d <- gg_violin(rAlpha, method = "m1", ylim_lower_alpha, ylim_upper_alpha) +
  labs(title= "",
       x="", y = bquote(paste( alpha, " estimated by mode", sep = ""))) +
  geom_hline(yintercept = alpha, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 2
rAlpha <- na.omit(melt(cbind(m2$modeAlpha, ML_alpha, MAP_alpha)))
e <- gg_violin(rAlpha, method = "all", ylim_lower_alpha, ylim_upper_alpha)  +
  labs(title= "",
       x="", y = bquote(paste( alpha, " estimated by mode", sep = ""))) +
  geom_hline(yintercept = alpha, color = "red")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 5
rAlpha <- na.omit(melt(cbind(m5$modeAlpha, ML_alpha, MAP_alpha)))
# Basic violin plot
f <- gg_violin(rAlpha, method = "all", ylim_lower_alpha, ylim_upper_alpha) +
  labs(title= "",
       x="", y = bquote(paste( alpha, " estimated by mode", sep = ""))) +
  geom_hline(yintercept = alpha, color = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

g <- ggplot() + theme_minimal()

grid.arrange(a,b,c,d,e,f,g, ncol = 1, heights = c(unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(2, c("cm"))))


# x axis labels
line_update_method <- .015
line_bracket_update_method <- .035
line_bracket_estimate <- .05
line_trans_dens <- 0.065
line_prop_dens <- .04
sc <- .9
sc_v <- 1
offs <- .1
offs_m <- .07
offs_v_brack <- .05
grid.text("left-conditioned proposal", x = unit(.25*sc + offs, "npc"), y = unit(line_update_method, "npc"))
grid.text("modified bridge proposal", x = unit(.61*sc + offs, "npc"), y = unit(line_update_method, "npc"))
grid.text("Euler", x = unit(.15*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Milstein", x = unit(.33*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Euler", x = unit(.515*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Milstein", x = unit(.7*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Euler", x = unit(.1*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.195*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Euler", x = unit(.29*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.38*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Euler", x = unit(.47*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.56*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Euler", x = unit(.65*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.74*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("ML \nestimate", x = unit(.835*sc + offs, "npc"), y = unit(.0525, "npc"))
grid.text("MAP \nestimate", x = unit(.93*sc + offs, "npc"), y = unit(.0525, "npc"))

Brack(.42*sc + offs,line_bracket_update_method,.07*sc + offs,line_bracket_update_method,.012)
Brack(.78*sc + offs,line_bracket_update_method,.43*sc + offs,line_bracket_update_method,.012)

grid.text("proposal density:", x = unit(.064, "npc"), y = unit(line_prop_dens, "npc"), gp = gpar(fontface="italic"))
grid.text("likelihood function:", x = unit(.06, "npc"), y = unit(line_trans_dens, "npc"), gp = gpar(fontface="italic"))
grid.text("proposal method:", x = unit(.0625 , "npc"), y = unit(line_update_method, "npc"), gp = gpar(fontface="italic"))

# y axis labels
Brack(line_bracket_estimate, 0.45*sc_v + offs_v_brack ,line_bracket_estimate, .05*sc_v + offs_v_brack ,.012, orientation = "vertical")
Brack(line_bracket_estimate, .915*sc_v + offs_v_brack ,line_bracket_estimate, .5*sc_v + offs_v_brack ,.012, orientation = "vertical")
grid.text("m = 1", x = unit(offs_m , "npc"), y = unit(.915, "npc"))
grid.text("m = 2", x = unit(offs_m, "npc"), y = unit(.76, "npc"))
grid.text("m = 5", x = unit(offs_m, "npc"), y = unit(.605, "npc"))
grid.text("m = 1", x = unit(offs_m, "npc"), y = unit(.45, "npc"))
grid.text("m = 2", x = unit(offs_m, "npc"), y = unit(.295, "npc"))
grid.text("m = 5", x = unit(offs_m, "npc"), y = unit(.145, "npc"))

grid.text("mode", x = unit(.02, "npc"), y = unit(.3, "npc"), rot = 90)
grid.text("mean", x = unit(.02, "npc"), y = unit(.76, "npc"), rot = 90)

r <- dev.off()


#
#
##------plots sigma^2----------------------------------------------------------------
if(filetype == "pdf"){
  pdf(paste("figures_and_tables/Fig_9_violin_plots_sigma_M_", M, "_alpha_", theta[1], "_sigma_", theta[2], "_nIter_", nIter, ".pdf",
            sep = ""), width = 12, height = 16)
}else{
  setEPS()
  postscript(paste("figures_and_tables/Fig_9_violin_plots_sigma_M_", M, "_alpha_", theta[1], "_sigma_", theta[2], "_nIter_", nIter, ".eps",
                   sep = ""), width = 12, height = 16)
}


# sigma^2 estimate from the mean of the MCMC posteriori---------------------
# m = 1
colnames(m1$meanSigma) <- c("td_Euler\nn = 100", "td_Euler\nn = 0" ,
                            "td_Milstein\nn = 100" ,"td_Milstein\nn = 0",
                            "td_Euler\nn2 = 100", "td_Euler\nn2 = 0" ,
                            "td_Milstein\nn2 = 100" ,"td_Milstein\nn2 = 0")
rSigma <- melt(cbind(m1$meanSigma, ML_sigma, MAP_sigma))
a <- gg_violin(rSigma, method = "m1", ylim_lower_sigma, ylim_upper_sigma)+
  labs(title= "",
       x="", y = bquote(paste( sigma^2, " estimated by mean", sep = ""))) +
  geom_hline(yintercept = sigma_2, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 2
rSigma <- na.omit(melt(cbind(m2$meanSigma, ML_sigma, MAP_sigma)))
b <- gg_violin(rSigma, method = "all", ylim_lower_sigma, ylim_upper_sigma) +
  labs(title= "",
       x="", y = bquote(paste( sigma^2, " estimated by mean", sep = ""))) +
  geom_hline(yintercept = sigma_2, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 5
rSigma <- na.omit(melt(cbind(m5$meanSigma, ML_sigma, MAP_sigma)))
c <- gg_violin(rSigma, method = "all", ylim_lower_sigma, ylim_upper_sigma) +
  labs(title= "",
       x="", y = bquote(paste( sigma^2, " estimated by mean", sep = ""))) +
  geom_hline(yintercept = sigma_2, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# sigma^2 estimate from the mode of the MCMC posteriori-------------------
# m = 1
colnames(m1$modeSigma) <- c("td_Euler\nn = 100", "td_Euler\nn = 0" ,
                            "td_Milstein\nn = 100" ,"td_Milstein\nn = 0",
                            "td_Euler\nn2 = 100", "td_Euler\nn2 = 0" ,
                            "td_Milstein\nn2 = 100" ,"td_Milstein\nn2 = 0")
rSigma <- melt(cbind(m1$modeSigma, ML_sigma, MAP_sigma))
d <- gg_violin(rSigma, method = "m1", ylim_lower_sigma, ylim_upper_sigma) +
  labs(title= "",
       x="", y = bquote(paste( sigma^2, " estimated by mode", sep = ""))) +
  geom_hline(yintercept = sigma_2, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 2
rSigma <- na.omit(melt(cbind(m2$modeSigma, ML_sigma, MAP_sigma)))
e <- gg_violin(rSigma, method = "all", ylim_lower_sigma, ylim_upper_sigma) +
  labs(title= "",
       x="", y = bquote(paste( sigma^2, " estimated by mode", sep = ""))) +
  geom_hline(yintercept = sigma_2, color = "red")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

# m = 5
rSigma <- na.omit(melt(cbind(modeSigma, ML_sigma, MAP_sigma)))
f <- gg_violin(rSigma, method = "all", ylim_lower_sigma, ylim_upper_sigma)  +
  labs(title= "",
       x="", y = bquote(paste( sigma^2, " estimated by mode", sep = ""))) +
  geom_hline(yintercept = sigma_2, color = "red") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,0,num_left_lines), "lines"))

g <- ggplot() + theme_minimal()

grid.arrange(a,b,c,d,e,f,g, ncol = 1, heights = c(unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(fig_height, c("cm")), unit(2, c("cm"))))



# x axis labels
line_update_method <- .015
line_bracket_update_method <- .035
line_bracket_estimate <- .05
line_trans_dens <- 0.065
line_prop_dens <- .04
sc <- .9
sc_v <- 1
offs <- .1
offs_m <- .07
offs_v_brack <- .05
grid.text("left-conditioned proposal", x = unit(.25*sc + offs, "npc"), y = unit(line_update_method, "npc"))
grid.text("modified bridge proposal", x = unit(.61*sc + offs, "npc"), y = unit(line_update_method, "npc"))
grid.text("Euler", x = unit(.15*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Milstein", x = unit(.33*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Euler", x = unit(.515*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Milstein", x = unit(.7*sc + offs, "npc"), y = unit(line_prop_dens, "npc"))
grid.text("Euler", x = unit(.1*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.195*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Euler", x = unit(.29*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.38*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Euler", x = unit(.47*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.56*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Euler", x = unit(.65*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("Milstein", x = unit(.74*sc + offs, "npc"), y = unit(line_trans_dens, "npc"))
grid.text("ML \nestimate", x = unit(.835*sc + offs, "npc"), y = unit(.0525, "npc"))
grid.text("MAP \nestimate", x = unit(.93*sc + offs, "npc"), y = unit(.0525, "npc"))

Brack(.42*sc + offs,line_bracket_update_method,.07*sc + offs,line_bracket_update_method,.012)
Brack(.78*sc + offs,line_bracket_update_method,.43*sc + offs,line_bracket_update_method,.012)

grid.text("proposal density:", x = unit(.064, "npc"), y = unit(line_prop_dens, "npc"), gp = gpar(fontface="italic"))
grid.text("likelihood function:", x = unit(.06, "npc"), y = unit(line_trans_dens, "npc"), gp = gpar(fontface="italic"))
grid.text("proposal method:", x = unit(.0625 , "npc"), y = unit(line_update_method, "npc"), gp = gpar(fontface="italic"))

# y axis labels
Brack(line_bracket_estimate, 0.45*sc_v + offs_v_brack ,line_bracket_estimate, .05*sc_v + offs_v_brack ,.012, orientation = "vertical")
Brack(line_bracket_estimate, .915*sc_v + offs_v_brack ,line_bracket_estimate, .5*sc_v + offs_v_brack ,.012, orientation = "vertical")
grid.text("m = 1", x = unit(offs_m , "npc"), y = unit(.915, "npc"))
grid.text("m = 2", x = unit(offs_m, "npc"), y = unit(.76, "npc"))
grid.text("m = 5", x = unit(offs_m, "npc"), y = unit(.605, "npc"))
grid.text("m = 1", x = unit(offs_m, "npc"), y = unit(.45, "npc"))
grid.text("m = 2", x = unit(offs_m, "npc"), y = unit(.295, "npc"))
grid.text("m = 5", x = unit(offs_m, "npc"), y = unit(.145, "npc"))

grid.text("mode", x = unit(.02, "npc"), y = unit(.3, "npc"), rot = 90)
grid.text("mean", x = unit(.02, "npc"), y = unit(.76, "npc"), rot = 90)


r <- dev.off()
