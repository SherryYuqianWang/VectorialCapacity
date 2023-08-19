rm(list = ls())
#setwd("~/MariaDecayAnalysis")

source("D:/Basel/master/2 R/VC/function/relchange_entomology.R")
parameter_value = read.csv('D:/Basel/master/2 R/VC/parameter/parameter_value.csv')
    
                    
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(caret)

##############################
set.seed(123)
N_samples = 1000
N_b = 1000 # human population size


value<-parameter_value[parameter_value$species == "sinensis",2:13]

HBI_mean = value[[1]]
HBI_min = value[[2]]
HBI_max = value[[3]]
parous_mean = value[[4]]
parous_min = value[[5]]
parous_max = value[[6]]
sac_mean = value[[7]]
sac_min = value[[8]]
sac_max = value[[9]]
tau_mean = value[[10]]  #the same as minimus(no data available)
tau_min = value[[11]]
tau_max = value[[12]]

plot.data <- data.frame(matrix(ncol = 8, nrow = N_samples))
names(plot.data) = c('x_HBI', "x_parous", 'x_sac', 'x_tau', 'y_HBI','y_parous','y_sac','y_tau')

#HBI
HBI_samples = runif(N_samples, HBI_min, HBI_max)
parous_samples = rep(parous_mean, each = N_samples)
sac_samples = rep(sac_mean, each = N_samples)
tau_samples = rep(tau_mean, each = N_samples)

#plot
plot.data$x_HBI = HBI_samples
plot.data$y_HBI = VCrelchange_entomology(HBI_mean, HBI_samples, parous_mean, parous_samples,
                                                                 sac_mean, sac_samples, tau_mean, tau_samples,
                                                                 N_b,N_samples)

#parous
HBI_samples = rep(HBI_mean, each = N_samples)
parous_samples = runif(N_samples, parous_min, parous_max)
sac_samples = rep(sac_mean, each = N_samples)
tau_samples = rep(tau_mean, each = N_samples)


plot.data$x_parous = parous_samples
plot.data$y_parous = VCrelchange_entomology(HBI_mean, HBI_samples, parous_mean, parous_samples,
                                             sac_mean, sac_samples, tau_mean, tau_samples,
                                             N_b,N_samples)

#sac
HBI_samples = rep(HBI_mean, each = N_samples)
parous_samples = rep(parous_mean, each = N_samples)
sac_samples = runif(N_samples, sac_min, sac_max)
tau_samples = rep(tau_mean, each = N_samples)

plot.data$x_sac = sac_samples
plot.data$y_sac = VCrelchange_entomology(HBI_mean, HBI_samples, parous_mean, parous_samples,
                                                 sac_mean, sac_samples, tau_mean, tau_samples,
                                                 N_b,N_samples)
#tau
HBI_samples = rep(HBI_mean, each = N_samples)
parous_samples = rep(parous_mean, each = N_samples)
sac_samples = rep(sac_mean, each = N_samples)
tau_samples = floor(runif(N_samples, tau_min, tau_max+1))

plot.data$x_tau = tau_samples
plot.data$y_tau = VCrelchange_entomology(HBI_mean, HBI_samples, parous_mean, parous_samples,
                                                 sac_mean, sac_samples, tau_mean, tau_samples,
                                                 N_b,N_samples)


process <- preProcess(as.data.frame(plot.data[,1:4]), method=c("range"))
norm_scale <- predict(process, as.data.frame(plot.data[,1:4]))
plot.data[,1:4]<-norm_scale


plot_tidy<-plot.data %>% 
  pivot_longer(
    cols = everything(), 
    cols_vary = "slowest",
    names_to = c(".value", "parameter"), 
    names_sep = "_"
  )


p = ggplot(data=plot_tidy) +
  geom_point(aes(y = y, x = x,color=parameter)) +
  #facet_wrap(~parameter,nrow=2,scales = "free_x")+
  scale_y_continuous(trans = pseudo_log_trans(), breaks=c(-1,0,10, 100,1000))+
  #  geom_point(aes(y = Output, x = parous)) +
  #  geom_point(aes(y = Output, x = sac)) +
  #  geom_point(aes(y = Output, x = tau)) +
  #ylim(-1,10)+
   xlab('Normalized parameter value') +
   ylab('Relative change in VC') +
 #  ylab('Relative change in VC (pseudo log transform)') +
  theme_bw()
p


ggsave("dirus_normalized.png", width = 4, height = 4)
ggsave("minimus_normalized.png", width = 4, height = 4)
ggsave("sinensis_normalized.png", width = 4, height = 4)
ggsave("sundaicus_normalized.png", width = 4, height = 4)
ggsave("maculatus_normalized.png", width = 4, height = 4)



