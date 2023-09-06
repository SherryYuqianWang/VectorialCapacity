rm(list = ls())

source("~/1_function/relchange_entomology_GSA.R") #load the function
parameter_value = read.csv('~/2_parameter/bionomics_value.csv')

library(lhs)
library(sensitivity)

##############################
set.seed(123)


#default parameters

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
tau_mean = value[[10]]
tau_min = value[[11]]
tau_max = value[[12]]

#LHS-PRCC
n_params = 4  #number of parameters
N_samples = 500  # number of sample set
repeats = 50   # number of repeat

prcc<-data.frame(matrix(ncol = repeats ,nrow = n_params))

for (r in 1:repeats){
  A<-randomLHS(N_samples,n_params)
  X<-data.frame(
     HBI_samples=qunif(A[,1],min = HBI_min, max = HBI_max),
     parous_samples=qunif(A[,2],min = parous_min, max = parous_max),
     sac_samples=qunif(A[,3],min = sac_min, max = sac_max),
     tau_samples=floor(qunif(A[,4],min = tau_min, max = tau_max+1))
  )


  y=rep(0,N_samples)
  
  #Calculate PRCC
  for (s in 1:N_samples){
  y[s] <- with(X[s,],VCrelchange_entomology_GSA(HBI_mean, HBI_samples, parous_mean, parous_samples,
                              sac_mean, sac_samples, tau_mean, tau_samples,
                              N_b))
  
  #Calculate PRCCs for each parameter set as PRCC
    x <- pcc(X, y, rank=TRUE)
    prcc[,r]<-x$PRCC
  }
}


prcc_new = as.data.frame(t(prcc))
colnames(prcc_new) = c('HBI','parous','sac', 'tau')

prcc_plot<-data.frame(matrix(ncol = 4, nrow = n_params))
prcc_plot[,1]<-sapply(prcc_new,quantile, probs=0.975)
prcc_plot[,2]<-sapply(prcc_new,mean)
prcc_plot[,3]<-sapply(prcc_new,quantile, probs=0.025)
prcc_plot[,4]<-c("HBI","parous","sac","tau")
names(prcc_plot) = c('high.ci','mean','low.ci', 'parameters')

saveRDS(prcc_plot,file="~/5_output/ent_prcc")
