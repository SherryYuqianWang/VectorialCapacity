rm(list = ls())

source("~/1_function/relchange_intervention2_GSA.R")
parameter_value = read.csv('~/2_parameter/parameter_value.csv')
int = read.csv('~/2_parameter/intervention_value.csv')

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
library(qrng)


##############################
#set.seed(123)
N_samples = 10000 #number of samples generated from 
n_params = 4 #number of parameters
repeats = 10 #number of repeats for averaging

#default parameters
N_b=1000 
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

sobol<-list()
DT_plot <- data.table("mean" = numeric(),    # Create an empty data.table
                   "low.ci" = numeric(),
                   "high.ci" = numeric(),
                   "order" = character(),
                   "parameter" = character(),
                   "sac" = numeric(),
                   "coverage" = numeric()
                   )


N <- 1000
params <- c("betar","betam","betad","xi")
order <- "first"
R<-50


X<-sobol_matrices(
  matrices = c("A", "B", "AB"),
  N=N,
  params=params,
  order = order,
  type = "QRN",
)


X[, 1] <- qunif(X[,1], min = betar_min, max = betar_max)
X[, 2] <- qunif(X[,2],min = betam_min, max = betam_max)
X[, 3] <- qunif(X[,3],min = betad_min, max = betad_max)
X[, 4] <- qunif(X[,4],min = xi_min, max = xi_max)

#datatable for output
X1<-data.table(X)

for (i in 1:3){
  coverage=coverage_sample[i]
 
for(v in 1:length(sac_range)){
  sac_mean=sac_range[v]
  y=rep(0,nrow(X))
  for (s in 1:nrow(X)){
    y[s] <- with(X1[s,], VCrelchange_intervention_GSA(HBI_mean, betar, parous_mean, betam, betad,
                                                   sac_mean, xi, tau_mean,
                                                   N_b,coverage))
  }
  
    #Calculate sobols for each parameter set as sobol
  ind <- sobol_indices(Y = y, N = N, params = params, boot = TRUE, R = R)
  sobol<-data.table(ind$results)
  DT_sobol<-data.table("mean" = sobol$original,
                   "low.ci" = sobol$low.ci,
                   "high.ci" = sobol$high.ci,
                   "order" = sobol$sensitivity,
                   "parameter" = sobol$parameters,
                   "sac" = sac_mean,
                   "coverage" = coverage
                   )
  DT_plot<-rbind(DT_plot,DT_sobol)
}

}

saveRDS(DT_plot,file="~/5_output/int_sac_sobol")