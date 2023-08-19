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
#coverage=0.3 #coverage

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


parous_range=seq(0.39,0.78,by=0.01)
coverage_sample=c(0.1,0.3,0.5)

###calculate sobol
#sobol<-list()


for(v in 1:length(HBI_range)){
  HBI_mean=HBI_range[v]
  sobol[[v]]<-data.frame(matrix(ncol = repeats, nrow = n_params))
  
  #generate sample  
#  sobol[[v]]<- foreach (r = 1:repeats,.combine = cbind)%dopar%{
    A<-randomLHS(n_samples,n_params)
    X<-data.frame(
      pi_samples=qunif(A[,1],min = 0.3, max = 1.0),
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
    
    #Calculate sobols for each parameter set as sobol
    x <- pcc(X, y, rank=TRUE)
    x$sobol
  }
}


sobol<-list()
DT_plot <- data.table("mean" = numeric(),    # Create an empty data.table
                   "low.ci" = numeric(),
                   "high.ci" = numeric(),
                   "order" = character(),
                   "parameter" = character(),
                   "parous" = numeric(),
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

#matrix for scatter plot
X[, 1] <- qunif(X[,1], min = betar_min, max = betar_max)
X[, 2] <- qunif(X[,2],min = betam_min, max = betam_max)
X[, 3] <- qunif(X[,3],min = betad_min, max = betad_max)
X[, 4] <- qunif(X[,4],min = xi_min, max = xi_max)
#X[, 5] <- qunif(X[,5],min = rho_min, max = rho_max)

#datatable for output
X1<-data.table(X)

#only use for loop
for (i in 1:3){
#  sobol<-list()
  coverage=coverage_sample[i]
 
for(v in 1:length(parous_range)){
  parous_mean=parous_range[v]
  #generate sample
  y=rep(0,nrow(X))
  for (s in 1:nrow(X)){
    y[s] <- with(X1[s,], VCrelchange_intervention2_GSA(HBI_mean, betar, parous_mean, betam, betad,
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
                   "parous" = parous_mean,
                   "coverage" = coverage
                   )
  #l = list(DT_plot,DT_sobol)
  #DT_plot<-rbindlist(l)
  DT_plot<-rbind(DT_plot,DT_sobol)
}

}

saveRDS(DT_plot,file="~/malaria/VC/output/int_parous_sobol")

p<-ggplot(DT_plot[order=="Si"],aes(x=parous, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=0, ymax=0.05), fill="grey", color=NA,alpha=0.05) +
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1)+
  facet_grid(coverage~.)+
  #xlab('parous') +
  ylab("Sobol's first order indices") +
  theme_bw()
p

p<-ggplot(DT_plot[order=="Ti"],aes(x=parous, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=0, ymax=0.05), fill="grey", color=NA,alpha=0.05) +
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1)+
  facet_grid(coverage~.)+
  #xlab('parous') +
  ylab("Sobol's total order indices") +
  theme_bw()
p


ggsave("all_first_int2_parous.png", width = 8, height = 7)
ggsave("all_total_int2_parous.png", width = 8, height = 7)



###LHS

sobol<-data.frame(matrix(ncol = 6, nrow = n_params*length(parous_range)*length(coverage)))
#sobol_plot<-list()
#p<-list()
#only use for loop
for (c in 1:3){
  #  sobol<-list()
  coverage_sample=coverage[c]
  
  for(v in 1:length(parous_range)){
    parous_mean=parous_range[v]
    #generate sample
    A1<-sobol(n_samples, d = 4,randomize = "digital.shift")
    X1<-data.frame(
      pi_samples=qunif(A1[,1],min = 0.3, max = 1.0),
      theta_samples=qunif(A1[,2],min = 0.1, max = 2.5), 
      xi_samples=qunif(A1[,3],min = 0, max = 0.4),
      omega_samples=qunif(A1[,4],min = 0.1, max = 0.9))
    
    X1$theta_samples=X1$pi_samples*X1$theta_samples
    names(X1)[names(X1) == 'theta_samples'] <- 'k_samples'
    
    A2<-sobol(n_samples, d = 4,randomize = "digital.shift")
    X2<-data.frame(
      pi_samples=qunif(A2[,1],min = 0.3, max = 1.0),
      theta_samples=qunif(A2[,2],min = 0.1, max = 2.5), 
      xi_samples=qunif(A2[,3],min = 0, max = 0.4),
      omega_samples=qunif(A2[,4],min = 0.1, max = 0.9))
    
    X2$theta_samples=X2$pi_samples*X2$theta_samples
    names(X2)[names(X2) == 'theta_samples'] <- 'k_samples'
    
    #Calculate relative change of vectorial capacity
    x <- sobolSalt(model = NULL, X1, X2, nboot = 1000)
    X <- rbind(X1, X2)
    for (i in 1:n_params) {
      Xb <- X1
      Xb[,i] <- X2[,i]
      X <- rbind(X, Xb) 
    } 
    
    y=rep(0,nrow(X))
    for (s in 1:nrow(X)){
      y[s] <- with(X[s,], VCrelchange(HBI_mean, pi_samples, parous_mean, k_samples,
                                      sac_mean, xi_samples, tau_mean, omega_samples,
                                      N_b,coverage_sample))
    }
    
    tell(x,y)
    #Calculate sobols for each parameter set as sobol
    
    sobol[(n_params*length(HBI_range)*(c-1))+(n_params*v-3):(n_params*v),1]<-x$S[,1] #mean
    sobol[(n_params*length(HBI_range)*(c-1))+(n_params*v-3):(n_params*v),2]<-x$S[,4] #low
    sobol[(n_params*length(HBI_range)*(c-1))+(n_params*v-3):(n_params*v),3]<-x$S[,5] #high
    sobol[(n_params*length(HBI_range)*(c-1))+(n_params*v-3):(n_params*v),4]<-HBI_range[v]
    sobol[(n_params*length(HBI_range)*(c-1))+(n_params*v-3):(n_params*v),5]<-coverage_sample
    sobol[(n_params*length(HBI_range)*(c-1))+(n_params*v-3):(n_params*v),6]<-c("pi","k","xi","omega")
    
  }
}


names(sobol) = c('mean','lower','upper', 'HBI','coverage','parameter')
#sobol$coverage<-as.character(sobol$coverage)

p<-ggplot(sobol,aes(x=HBI, y=mean,col=parameter))+
  geom_line()+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  #xlab('HBI') +
  ylab('sobol') +
  facet_grid(coverage~.)+
  theme_minimal()
p



#plot
sobol_new<-list()

sobol_plot[[c]]<-data.frame(matrix(ncol = 5, nrow = length(HBI_range)*n_params))
#calculate mean and t-values
#sobol_means=rep(0,n_params)
#t_values=rep(0,n_params)
for(v in 1:length(HBI_range)){
  sobol_new[[v]] = data.frame(t(sobol[[v]]))
  
  sobol_plot[[c]][(n_params*v-3):(n_params*v),1]<-sapply(sobol_new[[v]],quantile, probs=0.975)
  sobol_plot[[c]][(n_params*v-3):(n_params*v),2]<-sapply(sobol_new[[v]],mean)
  sobol_plot[[c]][(n_params*v-3):(n_params*v),3]<-sapply(sobol_new[[v]],quantile, probs=0.025)
  sobol_plot[[c]][(n_params*v-3):(n_params*v),4]<-HBI_range[v]
  sobol_plot[[c]][(n_params*v-3):(n_params*v),5]<-c("pi","k","xi","omega")
    #t_values[i] <- sobol_means[i]*sqrt((n_samples-2)/(1-sobol_means[i]))

}

names(sobol_plot[[c]]) = c('upper','mean','lower', 'HBI','parameter')
#sobol_tidy <- sobol_means %>%
#  pivot_longer(
#    cols = c('pi','k','xi','omega'),
#    names_to = 'parameter',
#    values_to = 'sobol')

p[[c]]<-ggplot(sobol_plot[[c]],aes(x=HBI, y=mean,col=parameter))+
  geom_line()+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  #xlab('HBI') +
  ylab('sobol') +
  theme_minimal()

}

ggarrange(p[[1]],p[[2]],p[[3]] + rremove("x.text"), 
          labels = c("A", "B", "C","D"),common.legend = TRUE,legend = "right",
          ncol = 1, nrow = 3)

#



library(randtoolbox)

sobol(100, d = 4,randomize = "digital.shift")


# Test case : the non-monotonic Sobol g-function
# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]
library(sensitivity)
library(boot)
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))
# sensitivity analysis
x <- soboljansen(model = sobol.fun, X1, X2, nboot = 100)
x$X[3001,]
print(x)
plot(x)
library(ggplot2)
ggplot(x)
# Only for demonstration purposes: a model function returning a matrix
sobol.fun_matrix <- function(X){
  res_vector <- sobol.fun(X)
  cbind(res_vector, 2 * res_vector)
}
x_matrix <- soboljansen(model = sobol.fun_matrix, X1, X2)
plot(x_matrix, y_col = 2)
title(main = "y_col = 2")
# Also only for demonstration purposes: a model function returning a
# three-dimensional array
sobol.fun_array <- function(X){
  res_vector <- sobol.fun(X)
  res_matrix <- cbind(res_vector, 2 * res_vector)
  array(data = c(res_matrix, 5 * res_matrix),
        dim = c(length(res_vector), 2, 2))
}
x_array <- soboljansen(model = sobol.fun_array, X1, X2)
plot(x_array, y_col = 2, y_dim3 = 2)
title(main = "y_col = 2, y_dim3 = 2")




 #scatter plot
N_samples = 1000
N_b = 1000 # human population size
#intervention
pi_mean=0.65
pi_min=0.3
pi_max=1.0
k_mean=1.4
k_min=0.2
k_max=2.6
xi_mean=0.28
xi_min=0
xi_max=0.56
omega_mean=0.46
omega_min=0
omega_max=0.92


X<-data.frame(
  pi_samples = rep(pi_mean, each = N_samples),
  k_samples = rep(k_mean, each = N_samples),
  xi_samples = rep(xi_mean, each = N_samples),
  omega_samples = rep(omega_mean, each = N_samples)
)
#range
X[,1] = runif(N_samples, pi_min, pi_max)
X[,2] = runif(N_samples, k_min, k_max)
X[,3] = runif(N_samples, xi_min, xi_max)
X[,4] = runif(N_samples, omega_min, omega_max)



#calculate VC
y=rep(0,N_samples)
for (s in 1:N_samples){
  y[s] <- with(X[s,], VCrelchange(HBI_mean, pi_samples, parous_mean, k_samples,
                                  sac_mean, xi_samples, tau_mean, omega_samples,
                                  N_b))
  
}

#plot
plot.data <- data.frame(matrix(ncol = 5, nrow = N_samples))
names(plot.data) = c('pi', "k", 'xi', 'omega', 'Output')

plot.data$pi = X[,1]
plot.data$k = X[,2]
plot.data$xi = X[,3]
plot.data$omega = X[,4]
plot.data$Output = y



p = ggplot(data=plot.data) +
  #geom_point(aes(y = Output, x = pi)) +
  #geom_point(aes(y = Output, x = k)) +
  geom_point(aes(y = Output, x = xi)) +
  #geom_point(aes(y = Output, x = omega)) +
  # xlab('Coverage including adherence (%)') +
  # ylab('Relative reduction in vectorial capacity') +
  theme_minimal()
p


pi=runif(1000,0.3,1.0)
r=runif(1000,0.1,2.2)
k=pi*r
k=runif(100,0.2,2.2)
y=k/pi
hist(k)
summary()

df<-data.frame(x1=c(0.74,0.53,0.87,0.73,0.97,0.74,0.96,0.82,0.83,0.40,0.6,0.42,0.84,0.57,0.98,0.99,0.41,0.29),
x2=c(0.48,0.24,1.65,0.44,0.85,0.43,0.84,0.57,0.74,0.79,0.55,0.51,2.16,1.99,0.42,0.69,0.91,0.34),
x3=c(0.02,0.01,0.11,0.02,0.08,0.06,0.07,0.01,0.02,0.01,0.01,0.12,0.13,0.41,0.11,0.26,0.07,0.07),
x4=c(0.65,0.91,0.05,0.63,0.16,0.56,0.52,0.73,0.51,0.17,0.43,0.39,0.30,0.32,0.78,0.49,0.5,0.5),
x5=rep(1,18),
x6=rep(1,18))
names(df) = c('pi', "k", 'post_kill', 'omega', 'disarm','pre_kill')

df$x5=df$x2*df$x4
df$x6=df$x2*(1-df$x4)

plot(df$post_kill,df$pre_kill)
plot(df$post_kill,df$disarm)
plot(df$pi,df$pre_kill)
plot(df$pi,df$k)
Corr <- cor(sapply(df, as.numeric),
            use = "pairwise.complete.obs", method = "spearman")
corrplot::corrplot(Corr, method = "square", type = "upper",
                   tl.col = "black")
summary(df)
