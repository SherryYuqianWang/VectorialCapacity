rm(list = ls())
#setwd("~/MariaDecayAnalysis")

source("~/malaria/VC/function/relchange_intervention2_GSA.R")
parameter_value = read.csv('~/malaria/VC/parameter/parameter_value.csv')
int2 = read.csv('~/malaria/VC/parameter/int2_value.csv')

library(ggplot2)
library(lhs)
library(sensitivity)
library(tidyr)
library(dplyr)
library(abind)
library(profvis)
library(epiR)
library(Rmisc)
library(ggpubr)
library(data.table)

##############################
#set.seed(123)
n_samples = 300 #number of samples generated from LHS
n_params = 4 #number of parameters
repeats = 50 #number of repeats for averaging

#default parameters
N_b=1000 

value<-parameter_value[parameter_value$species == "all",2:13]

HBI_mean = value[[1]]
parous_mean = value[[4]]
sac_mean = value[[7]]
tau_mean = value[[10]]  


betar_min = int2$betar_min
betar_max = int2$betar_max
betam_min = int2$betam_min
betam_max = int2$betam_max
betad_min = int2$betad_min
betad_max = int2$betad_max
xi_min = int2$xi_min
xi_max = int2$xi_max
#rho_min = int2$rho_min
#rho_max = int2$rho_max


sac_range=seq(0.16,0.88,by=0.01)
coverage_sample=c(0.1,0.3,0.5)

###calculate PRCC
#prcc<-list()

###foreach
for(v in 1:length(HBI_range)){
  HBI_mean=HBI_range[v]
  prcc[[v]]<-data.frame(matrix(ncol = repeats, nrow = n_params))
  
  #generate sample  
  prcc[[v]]<- foreach (r = 1:repeats,.combine = cbind)%dopar%{
    A<-randomLHS(n_samples,n_params)
    X<-data.frame(
      pi_samples=qunif(A[,1],min = pi_min, max = pi_max),
      k_samples=qunif(A[,2],min = 0.2, max = 2.2),
      xi_samples=qunif(A[,3],min = 0, max = 0.4),
      omega_samples=qunif(A[,4],min = 0.1, max = 0.9))
    
    
    

    #Calculate relative change of vectorial capacity
    y=rep(0,n_samples)
    
    for (s in 1:n_samples){
      y[s] <- with(X[s,], VCrelchange(HBI_mean, pi_samples, parous_mean, k_samples,
                                      sac_mean, xi_samples, tau_mean, omega_samples,
                                      N_b,c))
    }
    
    #Calculate PRCCs for each parameter set as PRCC
    x <- pcc(X, y, rank=TRUE)
    x$PRCC
  }
}







t_bound<-qt((0.05/n_params)/2,n_samples-2-(n_params-1))
lower_bound<--sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))
upper_bound<-sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))

DT_plot <- data.table("mean" = numeric(),    # Create an empty data.table
                      "low.ci" = numeric(),
                      "high.ci" = numeric(),
                      "parameter" = character(),
                      "sac" = numeric(),
                      "coverage" = numeric())
#only use for loop
for (i in 1:3){
  
  coverage=coverage_sample[i]
  
for(v in 1:length(sac_range)){
  sac_mean=sac_range[v]
  prcc<-data.frame(matrix(ncol = repeats, nrow = n_params))

  #generate sample  
  for (r in 1:repeats){
    A<-randomLHS(n_samples,n_params)
    X<-data.table(
      betar=qunif(A[,1],min = betar_min, max = betar_max),
      betam=qunif(A[,2],min = betam_min, max = betam_max),
      betad=qunif(A[,3],min = betad_min, max = betad_max),
      xi=qunif(A[,4],min = xi_min, max = xi_max))
      #rho=qunif(A[,5],min = rho_min, max = rho_max))
  
    #Calculate relative change of vectorial capacity
    y=rep(0,n_samples)
    
    for (s in 1:n_samples){
      y[s] <- with(X[s,], VCrelchange_intervention2_GSA(HBI_mean, betar, parous_mean, betam, betad,
                                                        sac_mean, xi, tau_mean,
                                                        N_b,coverage))
    }
    
    #Calculate PRCCs for each parameter set as PRCC
    x <- pcc(X, y, rank=TRUE)
    prcc[,r]<-x$PRCC
    
  }
  DT_prcc<-data.table("mean" = apply(prcc,1,mean),
                       "low.ci" = apply(prcc,1, quantile,probs=0.025),
                       "high.ci" = apply(prcc, 1, quantile, probs=0.975),
                       "parameter" = c("betar","betam","betad","xi"),
                       "sac" = sac_mean,
                       "coverage" = coverage)
  DT_plot<-rbind(DT_plot,DT_prcc)
 }
}

saveRDS(DT_plot,file="~/malaria/VC/output/int_sac_PRCC")

#plot
p<-ggplot(DT_plot,aes(x=sac, y=mean,color=parameter))+
  #geom_rect(aes(xmin=min(parous_range), xmax=max(parous_range), ymin=lower_bound, ymax=upper_bound), fill="grey", color=NA,alpha=0.01) +
  geom_line(size=0.7)+
  #geom_hline(yintercept=upper_bound, linetype="dashed", color = "grey32",)+
  #geom_hline(yintercept=lower_bound, linetype="dashed", color = "grey32")+
  #annotate(geom="text", x=0.3, y=upper_bound-0.1, label="non statistically significant", size=3.5, color="grey32")+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  facet_grid(coverage~.)+
  xlab('Sac rate') +
  ylab("PRCC") +
  theme_bw()
p


ggsave("all_prcc_int2_sac.png", width = 8, height = 7)
