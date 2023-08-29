rm(list = ls())

#load the data for ploting
HBI_prcc<-readRDS("~/malaria/VC/output/int_HBI_PRCC")
parous_prcc<-readRDS("~/malaria/VC/output/int_parous_PRCC")
sac_prcc<-readRDS("~/malaria/VC/output/int_sac_PRCC")
HBI_sobol<-readRDS("~/malaria/VC/output/int_HBI_sobol")
parous_sobol<-readRDS("~/malaria/VC/output/int_parous_sobol")
sac_sobol<-readRDS("~/malaria/VC/output/int_sac_sobol")


library(ggpubr)
library(cowplot)

coverage.labs <- c("10% coverage", "30% coverage", "50% coverage")
names(coverage.labs) <- c("0.1", "0.3", "0.5")

n_params=4
n_samples=500

t_bound<-qt((0.05/n_params)/2,n_samples-2-(n_params-1))
lower_bound<--sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))
upper_bound<-sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))

p_HBI_PRCC<-ggplot(HBI_prcc,aes(x=HBI, y=mean,col=parameter))+
  geom_line(size=0.7)+
  geom_hline(yintercept=upper_bound, linetype="dashed", color = "grey32",)+
  geom_hline(yintercept=lower_bound, linetype="dashed", color = "grey32")+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  facet_grid(coverage~., labeller=labeller(coverage = coverage.labs))+
  xlab('Human blood index') +
  ylab("PRCC") +
  theme_bw()+
  theme(legend.position = "none")

p_parous_PRCC<-ggplot(parous_prcc,aes(x=parous, y=mean,col=parameter))+
  geom_line(size=0.7)+
  geom_hline(yintercept=upper_bound, linetype="dashed", color = "grey32",)+
  geom_hline(yintercept=lower_bound, linetype="dashed", color = "grey32")+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  facet_grid(coverage~., labeller=labeller(coverage = coverage.labs))+
  xlab('Parity proportion') +
  ylab("PRCC") +
  theme_bw()+
  theme(legend.position = "none")

p_sac_PRCC<-ggplot(sac_prcc,aes(x=sac, y=mean,col=parameter))+
  geom_line(size=0.7)+
  geom_hline(yintercept=upper_bound, linetype="dashed", color = "grey32",)+
  geom_hline(yintercept=lower_bound, linetype="dashed", color = "grey32")+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  facet_grid(coverage~., labeller=labeller(coverage = coverage.labs))+
  xlab('Sac proportion') +
  ylab("PRCC") +
  theme_bw()+
  theme(legend.position = "none")

p_HBI_sobol<-ggplot(HBI_sobol[order=="Ti"],aes(x=HBI, y=mean,col=parameter))+
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1.06)+
  facet_grid(coverage~., labeller=labeller(coverage = coverage.labs))+
  ylab("Sobol's total indices") +
  theme_bw()+
  scale_fill_discrete(labels=c('Disarm', 'Pre-kill','Repel','Post-kill'),name="Effect")+
  scale_color_hue(labels=c('Disarm', 'Pre-kill','Repel','Post-kill'),name="Effect")

p_parous_sobol<-ggplot(parous_sobol[order=="Ti"],aes(x=parous, y=mean,col=parameter))+
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1.06)+
  facet_grid(coverage~., labeller=labeller(coverage = coverage.labs))+
  xlab('Parity proportion') +
  ylab("Sobol's total indices") +
  theme_bw()+
  scale_fill_discrete(labels=c('Disarm', 'Pre-kill','Repel','Post-kill'),name="Effect")+
  scale_color_hue(labels=c('Disarm', 'Pre-kill','Repel','Post-kill'),name="Effect")


p_sac_sobol<-ggplot(sac_sobol[order=="Ti"],aes(x=sac, y=mean,col=parameter))+
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1.06)+
  facet_grid(coverage~., labeller=labeller(coverage = coverage.labs))+
  xlab('Sac proportion') +
  ylab("Sobol's total indices") +
  theme_bw()+  
  scale_fill_discrete(labels=c('Disarm', 'Pre-kill','Repel','Post-kill'),name="Effect")+
  scale_color_hue(labels=c('Disarm', 'Pre-kill','Repel','Post-kill'),name="Effect")


ggarrange(p_HBI_PRCC,p_HBI_sobol,
          labels = c("A", "B"),
          align = "h",
          widths = c(0.45, 0.58),
          ncol = 2, nrow = 1)
ggsave("int_HBI.png", width = 8, height = 4)

ggarrange(p_parous_PRCC,p_parous_sobol,
          labels = c("A", "B"),
          align = "h",
          widths = c(0.45, 0.58),
          ncol = 2, nrow = 1)
ggsave("int_parous.png", width = 8, height = 4)

ggarrange(p_sac_PRCC,p_sac_sobol,
          labels = c("A", "B"),
          align = "h",
          widths = c(0.45, 0.58),
          ncol = 2, nrow = 1)
ggsave("int_sac.png", width = 8, height = 4)
