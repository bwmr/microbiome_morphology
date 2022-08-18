rm(list = ls())

library('tidyverse')
library('cowplot')

# Import Distance vs. LCL dataset
distance_df <- read.csv('results/distance_vis_vs_genome_vs_LCL.csv', row.names = 1)

# Plot by LCL
lcl <- ggplot(distance_df, aes(fct_inorder(LCL), cryoEM_Distance, fill = fct_relevel(LCL,"Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain")))+
  geom_boxplot(outlier.size=1, outlier.alpha = 0.5, outlier.color = NULL)+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme_half_open(font_size = 10)+
  background_grid()+
  scale_x_discrete(limits = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain"), labels = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain","Between Domains"))+
  scale_fill_brewer(type = "seq", direction = 1)+
  xlab(element_blank())+
  ylab('Euclidean Distance')+  
  theme(legend.position="none", axis.text = element_text(size = 10, family = "sans", face = "plain"), axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8), axis.text.y = element_text(vjust = 0.5))

ggsave('Fig5-Diversity-PanelA.pdf', plot = lcl, device = 'pdf', width = 70, height = 80, units = c("mm"), dpi = 300, limitsize = FALSE)
