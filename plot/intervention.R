rm(list = ls())
HBI_prcc<-readRDS("~/malaria/VC/output/int_HBI_PRCC")
parous_prcc<-readRDS("~/malaria/VC/output/int_parous_PRCC")
sac_prcc<-readRDS("~/malaria/VC/output/int_sac_PRCC")
HBI_sobol<-readRDS("~/malaria/VC/output/int_HBI_sobol")
parous_sobol<-readRDS("~/malaria/VC/output/int_parous_sobol")
sac_sobol<-readRDS("~/malaria/VC/output/int_sac_sobol")




library(ggpubr)
library(cowplot)


n_params=4
n_samples=500

t_bound<-qt((0.05/n_params)/2,n_samples-2-(n_params-1))
lower_bound<--sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))
upper_bound<-sqrt(t_bound^2/((n_samples-2-(n_params-1))+(t_bound^2)))

p_HBI_PRCC<-ggplot(HBI_prcc,aes(x=HBI, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=lower_bound, ymax=upper_bound), fill="grey", color=NA,alpha=0.05) +
  geom_line(size=0.7)+
  #geom_smooth(method = "loess",size=0.7)+
  geom_hline(yintercept=upper_bound, linetype="dashed", color = "grey32",)+
  geom_hline(yintercept=lower_bound, linetype="dashed", color = "grey32")+
  #annotate(geom="text", x=0.2, y=upper_bound-0.1, label="non statistically significant", size=3.5, color="grey32")+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  facet_grid(coverage~.)+
  xlab('HBI') +
  ylab("PRCC") +
  theme_bw()+
  theme(legend.position = "none")

p_parous_PRCC<-ggplot(parous_prcc,aes(x=parous, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=lower_bound, ymax=upper_bound), fill="grey", color=NA,alpha=0.05) +
  geom_line(size=0.7)+
  #geom_smooth(method = "loess",size=0.7)+
  geom_hline(yintercept=upper_bound, linetype="dashed", color = "grey32",)+
  geom_hline(yintercept=lower_bound, linetype="dashed", color = "grey32")+
  #annotate(geom="text", x=0.2, y=upper_bound-0.1, label="non statistically significant", size=3.5, color="grey32")+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  facet_grid(coverage~.)+
  xlab('Parity proportion') +
  ylab("PRCC") +
  theme_bw()+
  theme(legend.position = "none")

p_sac_PRCC<-ggplot(sac_prcc,aes(x=sac, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=lower_bound, ymax=upper_bound), fill="grey", color=NA,alpha=0.05) +
  geom_line(size=0.7)+
  #geom_smooth(method = "loess",size=0.7)+
  geom_hline(yintercept=upper_bound, linetype="dashed", color = "grey32",)+
  geom_hline(yintercept=lower_bound, linetype="dashed", color = "grey32")+
  #annotate(geom="text", x=0.2, y=upper_bound-0.1, label="non statistically significant", size=3.5, color="grey32")+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(-1,1)+
  facet_grid(coverage~.)+
  xlab('Sac proportion') +
  ylab("PRCC") +
  theme_bw()+
  theme(legend.position = "none")

p_HBI_sobol<-ggplot(HBI_sobol[order=="Ti"],aes(x=HBI, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=0, ymax=0.05), fill="grey", color=NA,alpha=0.05) +
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1.06)+
  facet_grid(coverage~.)+
  #xlab('parous') +
  ylab("Sobol's total indices") +
  theme_bw()+
  #theme(legend.position = "none")
  scale_fill_discrete(labels=c('disarm', 'pre-kill','repell','post-kill'),name="effect")+
  scale_color_hue(labels=c('disarm', 'pre-kill','repell','post-kill'),name="effect")

p_parous_sobol<-ggplot(parous_sobol[order=="Ti"],aes(x=parous, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=0, ymax=0.05), fill="grey", color=NA,alpha=0.05) +
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1.06)+
  facet_grid(coverage~.)+
  xlab('Parity proportion') +
  ylab("Sobol's total indices") +
  theme_bw()+
  #theme(legend.position = "none")
  scale_fill_discrete(labels=c('disarm', 'pre-kill','repell','post-kill'),name="effect")+
  scale_color_hue(labels=c('disarm', 'pre-kill','repell','post-kill'),name="effect")


p_sac_sobol<-ggplot(sac_sobol[order=="Ti"],aes(x=sac, y=mean,col=parameter))+
  #geom_rect(aes(xmin=min(HBI_range), xmax=max(HBI_range), ymin=0, ymax=0.05), fill="grey", color=NA,alpha=0.05) +
  geom_line()+
  geom_ribbon(aes(ymin=low.ci,ymax=high.ci,fill=parameter),alpha=0.2,colour = NA)+
  ylim(0,1.06)+
  facet_grid(coverage~.)+
  xlab('Sac proportion') +
  ylab("Sobol's total indices") +
  theme_bw()+  
  #theme(legend.position = "none")
  scale_fill_discrete(labels=c('disarm', 'pre-kill','repell','post-kill'),name="effect")+
  scale_color_hue(labels=c('disarm', 'pre-kill','repell','post-kill'),name="effect")


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

ggdraw() +
  draw_plot(p_sac_PRCC, x = 0, y = 0, width = .4, height = .99) +
  draw_plot(p_sac_sobol, x = 0.5, y = 0, width = .5, height = .99) +
  #draw_plot(p_parous_PRCC, x = 0, y = .33, width = .2, height = .33) +
  #draw_plot(p_parous_sobol, x = .5, y = .33, width = .2, height = .33) +
  #draw_plot(p_sac_PRCC, x = 0, y = 0, width = .2, height = .33) +
  #draw_plot(p_sac_sobol, x = .5, y = 0, width = .2, height = .33) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))

ggsave("int_sac.png", width = 8, height = 4)
