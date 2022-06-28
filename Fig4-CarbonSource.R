rm(list = ls())

library('tidyverse')
library('cowplot')
library('ggtext')
library('ggpubr')

# Re-Read scaled data, separate measurements from metadata, remove measurements for PM - OM:
mini_repl <- read.csv("results/cryoem_scaled_mini_replicates.csv", row.names = 1, header = T)
values <- read.csv("data/cryoem_measurements_replicates.csv")

# Calculate PCA
mini_pca <- prcomp(mini_repl[2:6], center = FALSE, scale. = FALSE)
mini_pca$results <- as.data.frame(mini_pca$x)
mini_pca$results$Sample_ID <- mini_repl$Sample_ID 

# Plot PCA with different carbon sources highlighted
btheta_df <- mini_pca$results[mini_pca$results$Sample_ID == "200714-1" | mini_pca$results$Sample_ID == "201019-1" ,]
rbromii_df <- mini_pca$results[mini_pca$results$Sample_ID == "211130-1" | mini_pca$results$Sample_ID == "211130-3" | mini_pca$results$Sample_ID == "211130-2",]
ct_df <- mini_pca$results[mini_pca$results$Sample_ID == "200120-1" | mini_pca$results$Sample_ID == "180628-1",]

bt_plot <- ggplot(btheta_df, aes(PC1, PC2, color = Sample_ID))+
  geom_point(size = 3, na.rm=TRUE)+
  scale_color_manual(name='*B. thetaiotaomicron*  \n DSM 2079',
                     labels = c('On Glucose','On Starch/Maltose'),
                     values = c('200714-1'="#f768a1",'201019-1'="#ae017e"))+
  theme(legend.title = element_markdown(size = 10), legend.key = element_rect(fill = 'white'))
bt_legend <- cowplot::get_legend(bt_plot)
rm('bt_plot')

rb_plot <- ggplot(rbromii_df, aes(PC1, PC2, color = Sample_ID))+
  geom_point(size = 3, na.rm=TRUE)+
  scale_color_manual(name='*R. bromii* L2.63',
                     labels = c('On Fructose','On Pullulan','On Starch'),
                     values = c('211130-1'='#b2e2e2','211130-3' = '#66c2a4','211130-2'='#238b45'))+
  theme(legend.title = element_markdown(size = 10), legend.key = element_rect(fill = 'white'))
rb_legend <- cowplot::get_legend(rb_plot)
rm('rb_plot')

ct_plot <- ggplot(ct_df, aes(PC1, PC2, color = Sample_ID))+
  geom_point(size = 3, na.rm=TRUE)+
  scale_color_manual(name='*C. thermocellum*  \n DSM 1313',
                     labels = c('On Cellobiose','On MCC'),
                     values = c('200120-1'="#fe9929",'180628-1'="#cc4c02"))+
  theme(legend.title = element_markdown(size = 10), legend.key = element_rect(fill = 'white'))
ct_legend <- cowplot::get_legend(ct_plot)
rm('ct_plot')

# Plot Full PCA
full_plot <- ggplot(mini_pca$results, aes(PC1, PC2))+
  theme_minimal_grid(12)+
  geom_point(color = "#b3b3b3", size=3, alpha = 0.8, na.rm=TRUE, show.legend = FALSE)+
  xlab("PC1 (37.4%)")+
  ylab("PC2 (30.2%)")+
  geom_point(data=btheta_df, aes(PC1,PC2, color = Sample_ID), alpha = 0.9, size=3, na.rm=TRUE, show.legend = FALSE)+
  geom_point(data=rbromii_df, aes(PC1,PC2, color = Sample_ID), alpha = 0.9, size=3, na.rm=TRUE, show.legend = FALSE)+
  geom_point(data=ct_df, aes(PC1,PC2, color = Sample_ID), alpha = 0.9, size=3, na.rm=TRUE, show.legend = FALSE)+
  scale_color_manual(values = c('200714-1'="#f768a1",'201019-1'="#ae017e",
                                '211130-1'='#b2e2e2','211130-2'='#238b45', '211130-3' = '#66c2a4',
                                '200120-1'="#fe9929",'180628-1'="#cc4c02"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5),"cm"))
rm('btheta_df',"ct_df","rbromii_df")

# Plot B. theta comparison
btheta_df <- values[values$Sample_ID == "200714-1" | values$Sample_ID == "201019-1" ,]

bt_cw <- ggplot(btheta_df, aes(y = CW, x = Sample_ID))+
  theme_minimal_hgrid(12)+
  geom_jitter(aes(color = Sample_ID), show.legend = FALSE, size = 3, alpha = 0.9, width = 0.1, height = 0)+
  scale_x_discrete(labels = c('Glucose','Starch/Maltose'))+
  scale_color_manual(labels = c('On Glucose','On Starch/Maltose'),
                     values = c('200714-1'="#f768a1",'201019-1'="#ae017e"))+
  xlab(element_blank())+
  ylab("CW Thickness [nm]")+
  scale_y_continuous(expand = expansion(mult = c(0.2)))+
  stat_compare_means(method = "t.test", comparisons = list(c("200714-1","201019-1")), 
                     paired = FALSE, tip.length=0.03, vjust = -0.2, label = "p.signif", inherit.aes = FALSE,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  theme(plot.margin = unit(c(0.5,1,0,1), "cm"))

# Plot R. bromii comparison diameter
rbromii_df <- values[values$Sample_ID == "211130-1" | values$Sample_ID == "211130-3" | values$Sample_ID == "211130-2",]

rb_dia <- ggplot(rbromii_df, aes(y = Diameter, x = fct_relevel(Sample_ID,"211130-1","211130-3","211130-2")))+
  theme_minimal_hgrid(12)+
  geom_jitter(aes(color = Sample_ID), show.legend = FALSE, size = 3, alpha = 0.9, width = 0.1, height = 0)+
  scale_x_discrete(labels = c('Fructose','Pullulan','Starch'))+
  scale_color_manual(labels = c('On Fructose','On Pullulan','On Starch'),
                     values = c('211130-1'='#b2e2e2','211130-3'='#66c2a4','211130-2'='#238b45'))+
  xlab(element_blank())+
  ylab("Diameter [um]")+
  scale_y_continuous(expand = expansion(mult = c(0.4)))+
  stat_compare_means(method = "t.test", comparisons = list(c("211130-1","211130-3"),c("211130-1","211130-2")), 
                     paired = FALSE, tip.length=c(0.05,0.1), label = "p.signif", interhit.aes = FALSE, vjust = -0.1, step.increase = c(0.3),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  theme(plot.margin = unit(c(0,1,0.5,1), "cm"))

# Bring it all together
legends <- plot_grid(ct_legend, bt_legend, rb_legend, nrow = 1, ncol = 3, labels = NULL, align = "hv", greedy = FALSE, axis = "tl")
pca <- plot_grid(full_plot, legends, nrow = 2, ncol = 1, labels = c("a",NULL), label_size = 18, rel_heights = c(2,1), align = "v", axis = "l", greedy = FALSE)

comparisons <- plot_grid(bt_cw, rb_dia, NULL, nrow = 3, labels = c("b","c"), label_size = 18, align = "hv")

final <- plot_grid(pca, comparisons, nrow = 1, ncol = 2, labels = NULL, rel_widths = c(3,2))

ggsave('Fig4-CarbonSource.pdf', plot = final, device = 'pdf', width = 210, height = 148, units = c("mm"), dpi = 300)

