rm(list = ls())

#manuscript plot
prcc_plot<-readRDS("~/malaria/VC/output/ent_prcc")
plot<-readRDS("~/malaria/VC/output/ent_sobol")
ind.dummy<-readRDS("~/malaria/VC/output/ent_dummy")

#extend discussion plot
prcc_plot<-readRDS("~/malaria/VC/output/ent_prcc_new")
plot<-readRDS("~/malaria/VC/output/ent_sobol_new")
ind.dummy<-readRDS("~/malaria/VC/output/ent_dummy_new")



library(ggpubr)
library(cowplot)


n_params=7
n_samples=500

t_bound<-qt((0.05/n_params)/2,n_samples-2-(n_params-1))
lower_bound<--sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))
upper_bound<-sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))


p_ent_prcc<-ggplot(prcc_plot,aes(x=parameters, y=mean))+
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



p_ent_sobol<-ggplot(plot,aes(x=parameters, y=original, fill=sensitivity))+
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
  ylim(-0.01,1.26)+
  theme_bw()
#theme(legend.position = "none")

#manuscript plot
ggarrange(p_ent_prcc,p_ent_sobol,
          labels = c("A", "B"),
          align = "h",
          widths = c(0.45, 0.58),
          ncol = 2, nrow = 1)
ggsave("ent.png", width = 8, height = 3.5)


#extend discussion plot
ggarrange(p_ent_prcc,p_ent_sobol,
          labels = c("A", "B"),
          align = "h",
          widths = c(0.45, 0.58),
          ncol = 2, nrow = 1)
ggsave("ent_new.png", width = 8, height = 3)
