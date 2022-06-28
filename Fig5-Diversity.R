rm(list = ls())

library('tidyverse')
library('cowplot')
library('ggtext')

# Import Distance by Family Long
pairwise <- read.csv('results/distance_within_vs_between_families_upper.csv', row.names = 1)

# Import Genome Distance Long
distance_df <- read.csv('results/distance_vis_vs_genome_vs_LCL.csv', row.names = 1)

# Import Bf/Ar measurements
bf_ar <- read.csv('results/distance_between_within_Ar_Bf.csv', row.names = 1)

# Plot by LCL
lcl <- ggplot(distance_df, aes(fct_inorder(LCL), cryoEM_Distance, fill = fct_relevel(LCL,"Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain")))+
  geom_boxplot(outlier.size=1, outlier.alpha = 0.5, outlier.color = NULL)+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme_half_open(font_size = 10)+
  background_grid()+
  scale_x_discrete(limits = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain"), labels = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain","> Domain"))+
  scale_fill_brewer(type = "seq", direction = 1)+
  xlab(element_blank())+
  ylab('Visual Distance')+  
  theme(legend.position="none", axis.text = element_text(size = 10, family = "sans", face = "plain"), axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8), axis.text.y = element_text(vjust = 0.5))+
  theme(plot.margin = unit(c(1,0.5,0,1), "cm"))

# Plot heatmap full
heatmap <- ggplot(pairwise, aes(x = fct_relevel(X), y = fct_relevel(Y)))+
  geom_tile(aes(fill = Mean))+
  theme_minimal_grid(10)+
  scale_fill_fermenter(name='**Mean Dist.**', palette = "BuPu", direction = 1, na.value = "#b3b3b3")+
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), legend.title = element_markdown(vjust = 1), legend.text = element_text(angle = 60, hjust = 0.6))+
  theme(plot.margin = unit(c(0,0,0,2), "cm"))

legend_heatmap <- get_legend(heatmap)
rm('heatmap')

heatmap_nolegend <- ggplot(pairwise, aes(x = fct_relevel(X), y = fct_relevel(Y)))+
  geom_tile(aes(fill = Mean))+
  theme_minimal_grid(10)+
  scale_x_discrete(position = "top")+
  scale_fill_fermenter(palette = "BuPu", direction = 1, na.value = "#b3b3b3")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank())+
  xlab(element_blank())+
  ylab(element_blank())+
  theme(legend.position = "none", text = element_text(size = 10, family = "sans", face = "plain", color = "black"), axis.text.x = element_text(angle = 60, hjust = 0), axis.text.y = element_text(hjust = 1))+
  theme(plot.margin = unit(c(0,0.5,0,1), "cm"))

# Plot heatmap Bf / Ar

heatmap_bf_ar <- ggplot(bf_ar, aes(x = fct_relevel(X), y = fct_relevel(Y)))+
  geom_tile(aes(fill = Mean))+
  theme_minimal_grid(10)+
  scale_fill_fermenter(name='Mean Dist.', palette = "OrRd", direction = 1, na.value = "#b3b3b3")+
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), legend.title = element_blank(), legend.text = element_text(angle = 60, hjust = 0.6))+
  theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))

legend_bf_ar <- get_legend(heatmap_bf_ar)
rm(heatmap_bf_ar)

bf_ar_nolegend <- ggplot(bf_ar, aes(x = fct_relevel(X), y = fct_relevel(Y)))+
  geom_tile(aes(fill = Mean))+
  theme_minimal_grid(10)+
  scale_x_discrete(position = "top", breaks = c("Agathobacter_ruminis_DSM_29029","Butyrivibrio_fibrisolvens_CF3"), 
                   labels = element_blank())+
  scale_y_discrete(breaks = c("Butyrivibrio_fibrisolvens_CF3","Agathobacter_ruminis_DSM_29029"), 
                   labels = element_blank())+
  scale_fill_fermenter(palette = "OrRd", direction = 1, na.value = "#b3b3b3")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank())+
  xlab(element_blank())+
  ylab(element_blank())+
  theme(legend.position = "none", text = element_text(size = 10, family = "sans", face = "plain", color = "black"), axis.text.x = element_markdown(angle = 60, hjust = 0), axis.text.y = element_markdown(hjust = 1))+
  theme(plot.margin = unit(c(4,1,0,0), "cm"))

# fuse together

panel_b <- plot_grid(heatmap_nolegend, bf_ar_nolegend, legend_heatmap, legend_bf_ar, rel_widths = c(3,1), rel_heights = c(4,1), nrow = 2, ncol = 2)
  
final <- plot_grid(lcl, NULL, panel_b, NULL, nrow = 2, ncol = 2, labels = c("a","c","b","d"), label_size = 18, label_fontfamily = "sans")

ggsave('Fig5-Diversity.pdf', plot = final, device = 'pdf', width = 210, height = 148, units = c("mm"), dpi = 300, limitsize = FALSE)
