rm(list = ls())

source("~/1_function/relchange_intervention_GSA.R")
parameter_value = read.csv('~/2_parameter/bionomics_value.csv')
int = read.csv('~/2_parameter/intervention_value.csv')

library(lhs)
library(sensitivity)
library(data.table)

##############################
#set.seed(123)


#default parameters
N_b=1000 #total number of hosts

value<-parameter_value[parameter_value$species == "all",2:13]

HBI_mean = value[[1]]
parous_mean = value[[4]]
sac_mean = value[[7]]
tau_mean = value[[10]]  


betar_min = int$betar_min
betar_max = int$betar_max
betam_min = int$betam_min
betam_max = int$betam_max
betad_min = int$betad_min
betad_max = int$betad_max
xi_min = int$xi_min
xi_max = int$xi_max

sac_range=seq(0.16,0.88,by=0.01)
coverage_sample=c(0.1,0.3,0.5)

#LHS-PRCC
n_samples = 500 #number of samples generated from LHS
n_params = 4 #number of parameters
repeats = 50 #number of repeats for averaging

DT_plot <- data.table("mean" = numeric(),    # Create an empty data.table
                      "low.ci" = numeric(),
                      "high.ci" = numeric(),
                      "parameter" = character(),
                      "sac" = numeric(),
                      "coverage" = numeric())

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
    
    #Calculate relative change of vectorial capacity
    y=rep(0,n_samples)
    
    for (s in 1:n_samples){
      y[s] <- with(X[s,], VCrelchange_intervention_GSA(HBI_mean, betar, parous_mean, betam, betad,
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

saveRDS(DT_plot,file="~/5_output/int_sac_PRCC")
