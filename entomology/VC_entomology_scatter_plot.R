rm(list = ls())
#setwd("~/MariaDecayAnalysis")

source("D:/Basel/master/2 R/VC/function/relchange_entomology.R")
parameter_value = read.csv('D:/Basel/master/2 R/VC/parameter/parameter_value.csv')
    
                    
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
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

plot.data <- data.frame(matrix(ncol = 5, nrow = N_samples))
names(plot.data) = c('HBI', "parous", 'sac', 'tau', 'Output')
#mean
#HBI_samples = rep(HBI_mean, each = N_samples)
#parous_samples = rep(parous_mean, each = N_samples)
#sac_samples = rep(sac_mean, each = N_samples)
#tau_samples = rep(tau_mean, each = N_samples)
#range
set.seed(123)
HBI_samples = runif(N_samples, HBI_min, HBI_max)
parous_samples = runif(N_samples, parous_min, parous_max)
sac_samples = runif(N_samples, sac_min, sac_max)
tau_samples = floor(runif(N_samples, tau_min, tau_max+1))

#plot
plot.data$HBI = HBI_samples
plot.data$parous = parous_samples
plot.data$sac = sac_samples
plot.data$tau = tau_samples
plot.data$Output = VCrelchange_entomology(HBI_mean, HBI_samples, parous_mean, parous_samples,
                                          sac_mean, sac_samples, tau_mean, tau_samples,
                                          N_b,N_samples)

plot_tidy <- plot.data %>%
  pivot_longer(
    cols = c('HBI','parous','sac','tau'),
    names_to = 'parameter',
    values_to = 'value')

summary(plot.data$Output)

p = ggplot(data=plot_tidy) +
  geom_point(aes(y = Output, x = value),size=0.3,alpha=0.3) +
  facet_wrap(~parameter,nrow=2,scales = "free_x")+
  #scale_y_continuous(trans = pseudo_log_trans(), breaks=c(0,10, 100,1000))+
  #  geom_point(aes(y = Output, x = parous)) +
  #  geom_point(aes(y = Output, x = sac)) +
  #  geom_point(aes(y = Output, x = tau)) +
   xlab('parameter value') +
   ylab('Relative reduction in vectorial capacity') +
  theme_bw()
p


ggsave("dirus_ento.png", width = 4, height = 4)
ggsave("minimus_ento.png", width = 4, height = 4)
ggsave("sinensis_ento.png", width = 4, height = 4)
ggsave("sundaicus_ento.png", width = 4, height = 4)
ggsave("maculatus_ento.png", width = 4, height = 4)