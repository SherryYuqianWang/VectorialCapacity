rm(list = ls())
#setwd("~/MariaDecayAnalysis")

source("~/malaria/VC/function/relchange_entomology_GSA.R")
parameter_value = read.csv('~/malaria/VC/parameter/parameter_value.csv')


library(ggplot2)
library(lhs)
library(sensitivity)
library(tidyr)
library(dplyr)
#library(doSNOW)
library(sensobol)
library(data.table)

##############################
set.seed(123)
n_params = 5
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
tau_mean = value[[10]] 
tau_min = value[[11]]
tau_max = value[[12]]

#sobol<-data.frame(matrix(ncol = repeats*2, nrow = 4))

N <- 10000
params <- c("HBI","parous","sac","tau")
order <- "first"
R<-50
type <- "percent"
conf <- 0.95


X<-sobol_matrices(
  matrices = c("A", "B", "AB"),
  N=N,
  params=params,
  order = order,
  type = "QRN",
)


#matrix for scatter plot
X[, 1] <- qunif(X[, 1], min = HBI_min, max = HBI_max)
X[, 2] <- qunif(X[,2],min = parous_min, max = parous_max)
X[, 3] <- qunif(X[,3],min = sac_min, max = sac_max)
X[, 4] <- floor(qunif(X[,4],min = tau_min, max = tau_max+1))

#datatable for output
X1<-data.table(X)

y=rep(0,nrow(X1))
for (s in 1:nrow(X1)){
  y[s] <- with(X1[s,], VCrelchange_entomology_GSA(HBI_mean, HBI, parous_mean, parous,
                                                 sac_mean, sac, tau_mean, tau,
                                                 N_b))
}

#plot_scatter(data = X, N = N, Y = y, params = params)
##plot_multiscatter(data = X, N = N, Y = y, params = params, smpl = 2^14)


ind <- sobol_indices(Y = y, N = N, params = params, order=order, boot = TRUE, R = R)

plot<-data.frame(ind$results)

ind.dummy <- sobol_dummy(Y = y, N = N, params = params, boot = TRUE,
                         R = R)

saveRDS(plot,file="~/malaria/VC/output/ent_sobol")
saveRDS(ind.dummy,file="~/malaria/VC/output/ent_dummy")

p<-ggplot(plot,aes(x=parameters, y=original, fill=sensitivity))+
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines
           width = 0.9,
           linewidth=0.05) +      # Thinner lines
  geom_errorbar(aes(ymin=low.ci, ymax=high.ci), size= 0.4, width=0.2, position=position_dodge(0.9))+
  geom_hline(yintercept=ind.dummy$original[2], linetype="dashed", color = "blue")+
  geom_hline(yintercept=ind.dummy$original[1], linetype="dashed", color = "red")+
  xlab("Parameter") +
  ylab("Sobol' indices") +
  #scale_fill_hue(name="Indices type", # Legend label, use darker colors
  #               breaks=c("first", "total"),
  #               labels=c("first-order", "total-order")) +
  scale_fill_manual(values=c("#c6dbef","#08306b"),
                   labels=c('First-order', 'Total-order'),
                    name="Sobol's index") +
  #scale_fill_manual(values=c("#c6dbef","#9ecae1","#08306b")) +
  ylim(-0.01,1)+
  theme_bw()
  #theme(legend.position = "none")
  
p

ggsave("ento_sobol.png", width = 8, height = 7)

ggsave("dirus_sobol.png", width = 3, height = 4)
ggsave("minimus_sobol.png", width = 3, height = 4)
ggsave("sinensis_sobol.png", width = 3, height = 4)
ggsave("sundaicus_sobol.png", width = 3, height = 4)
ggsave("maculatus_sobol.png", width = 4, height = 4)
