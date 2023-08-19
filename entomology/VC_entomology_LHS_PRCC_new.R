rm(list = ls())
#setwd("~/MariaDecayAnalysis")


source("~/malaria/VC/function/relchange_entomology_GSA_new.R")
parameter_value = read.csv('~/malaria/VC/parameter/parameter_value.csv')

library(ggplot2)
library(lhs)
library(sensitivity)
library(tidyverse)

library(foreach)
library(doParallel)
library(doSNOW)

##############################
set.seed(123)
n_params = 7
repeats = 50
species = 4

#default parameters
N_samples = 500
N_b = 1000 # human population size

value<-parameter_value[parameter_value$species == "all",2:13]

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
theta_d_min = 5/24
theta_d_max = 12/24
theta_s_min = 9
theta_s_max = 15
N_b_min = 500
N_b_max = 2000
#prcc<-data.frame(matrix(ncol = repeats, nrow = n_params))


#cl<-makeCluster(4)
#registerDoSNOW(cl)

prcc<-foreach(r = 1:repeats, .combine='cbind')%dopar%{
  A<-randomLHS(n_samples,n_params)
  X<-data.frame(
    HBI_samples=qunif(A[,1],min = HBI_mean, max = HBI_max),
    parous_samples=qunif(A[,2],min = parous_min, max = parous_max),
    sac_samples=qunif(A[,3],min = sac_min, max = sac_max),
    tau_samples=1/2*floor(qunif(A[,4],min = tau_min*2, max = tau_max*2+1)),
    theta_s_samples=qunif(A[,5],min = theta_s_min, max = theta_s_max),
    theta_d_samples=1/2*floor(qunif(A[,6],min = theta_d_min*2, max = theta_d_max*2+1))
  )
  
  # sensitivity analysis
  y=rep(0,n_samples)
  
  #Calculate vectorial capacity
  for (s in 1:n_samples){
    y[s] <- with(X[s,], VCrelchange(HBI_mean, HBI_samples, parous_mean, parous_samples,
                                    sac_mean, sac_samples, tau_mean, tau_samples,theta_d_samples,
                                    theta_s_samples,N_b))
    
    
    #Calculate PRCCs for each parameter set as PRCC
    x <- pcc(X, y, rank=TRUE)
    x$PRCC
  }
}



prcc<-data.frame(matrix(ncol = repeats ,nrow = n_params))

for (r in 1:repeats){
#prcc<-foreach(r = 1:2,.combine=cbind)%dopar%{
  A<-randomLHS(N_samples, n_params)
  X<-data.frame(
     HBI_samples=qunif(A[,1],min = HBI_min, max = HBI_max),
     parous_samples=qunif(A[,2],min = parous_min, max = parous_max),
     sac_samples=qunif(A[,3],min = sac_min, max = sac_max),
     tau_samples=floor(qunif(A[,4],min = tau_min, max = tau_max+1)),
     theta_d_samples=qunif(A[,5],min = theta_d_min, max = theta_d_max),
     theta_s_samples=floor(qunif(A[,6],min = theta_s_min, max = theta_s_max+1)),
     N_b_samples=floor(qunif(A[,7],min = N_b_min, max = N_b_max+1))
  )

  
  # sensitivity analysis
  y=rep(0,N_samples)
  
  #Calculate PRCC
  for (s in 1:N_samples){
  y[s] <- with(X[s,],VCrelchange_entomology_GSA(HBI_mean, HBI_samples, parous_mean, parous_samples,
                              sac_mean, sac_samples, tau_mean, tau_samples, theta_d_samples, theta_s_samples, N_b_samples))
  
  #Calculate PRCCs for each parameter set as PRCC
    x <- pcc(X, y, rank=TRUE)
    prcc[,r]<-x$PRCC
  }
}


output<-cbind(X,y)
plot(output$sac_samples,output$y)

prcc_new = as.data.frame(t(prcc))
colnames(prcc_new) = c('HBI','parous','sac', 'tau','theta_d','theta_s','N')

prcc_plot<-data.frame(matrix(ncol = 4, nrow = n_params))
prcc_plot[,1]<-sapply(prcc_new,quantile, probs=0.975)
prcc_plot[,2]<-sapply(prcc_new,mean)
prcc_plot[,3]<-sapply(prcc_new,quantile, probs=0.025)
prcc_plot[,4]<-c("HBI","parous","sac","tau","theta_d","theta_s","N")
names(prcc_plot) = c('high.ci','mean','low.ci', 'parameters')


saveRDS(prcc_plot,file="~/malaria/VC/output/ent_prcc_new")

#calculate non significant t-values

t_bound<-qt((0.05/n_params)/2,N_samples-2-(n_params-1))

lower_bound<--sqrt(t_bound^2/((N_samples-2-(n_params-1))+(t_bound^2)))
upper_bound<-sqrt(t_bound^2/((N_samples-2-(n_params-1))+(t_bound^2)))


p<-ggplot(prcc_plot,aes(x=parameters, y=mean))+
  geom_bar(stat="identity",
           fill=c("#00467D"), # Use black outlines,
           linewidth=0.05,
           width = 0.5) +      # Thinner lines
  #geom_rect(aes(xmin=0, xmax=Inf, ymin=lower_bound, ymax=upper_bound), fill="grey",alpha=0.1) +  
  geom_hline(yintercept=upper_bound, linetype="dashed", color = "red")+
  geom_hline(yintercept=lower_bound, linetype="dashed", color = "red")+
  geom_hline(yintercept=0, color = "black",size=0.4)+
#  annotate("text", x=1, y=upper_bound, label="non significant", size=4, color="blue") +
#  annotate("text", x=1, y=lower_bound, label="non significant", size=4, color="blue") +
  geom_errorbar(aes(ymin=low.ci, ymax=high.ci), size= 0.4, width=0.15)+
  xlab("Parameter") +
  ylab("PRCC") +
  #scale_fill_hue(name="Indices type", # Legend label, use darker colors
  #               breaks=c("first", "total"),
  #               labels=c("first-order", "total-order")) +
  #scale_fill_manual(values=c("#08306b")) +
  ylim(-1,1)+
  theme_bw()

p


ggsave("dirus_PRCC.png", width = 3, height = 4)
ggsave("minimus_PRCC.png", width = 3, height = 4)
ggsave("sinensis_PRCC.png", width = 3, height = 4)
ggsave("sundaicus_PRCC.png", width = 3, height = 4)
ggsave("maculatus_PRCC.png", width = 3, height = 4)
