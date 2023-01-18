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

ggsave('Fig5-Diversity-PanelA.png', plot = lcl, device = 'png', width = 70, height = 70, units = c("mm"), dpi = 300, limitsize = FALSE)

# Calculate values for Panel B/C
keep <- c("Bacteroides_caccae_DSM_19024 Bacteroides_cellulosyliticus_CRE21_DSM_14838","Bacteroides_caccae_DSM_19024 Bacteroides_thetaiotaomicron_DSM_2079","Bacteroides_cellulosyliticus_CRE21_DSM_14838 Bacteroides_thetaiotaomicron_DSM_2079")
bacteroides <- df[(df$comp %in% keep),]
bact_agg <- aggregate.data.frame(bacteroides$cryoEM_Distance, by = list(bacteroides$comp), FUN = mean)

keep <- c("Butyrivibrio_fibrisolvens_CF3 Lachnospira_multipara_G6_ATCC_19207","Butyrivibrio_fibrisolvens_CF3 Roseburia_intestinalis_L1-82_DSM_14610",
          "Lachnospira_multipara_G6_ATCC_19207 Roseburia_intestinalis_L1-82_DSM_14610")
lachno <- df[(df$comp %in% keep),]
lachno_agg <- aggregate.data.frame(lachno$cryoEM_Distance, by = list(lachno$comp), FUN = mean)
