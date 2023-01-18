rm(list = ls())

library('tidyverse')
library('vegan')
library('gdata')
library('cowplot')

# Import values for cryoEM, select only intersect
scaled.mini_repl <- read.csv('results/cryoem_scaled_mini_replicates.csv', sep = ",", header = T, row.names = 1)

include = c("Marvinbryantia_formatexigens","Bacteroides_cellulosyliticus_CRE21","Bacteroides_caccae","Eubacterium_siraeum")
vis_cleaned <- scaled.mini_repl[(scaled.mini_repl$Strain %in% include), ] 

rm('scaled.mini_repl','include')

# Import WGS distance matrix, select only intersect
wgs_dist <- read.csv('data/wgs_distance_jaccard.csv',sep = ",", header = TRUE, row.names = 1, check.names = FALSE)

include = c("Marvinbryantia_formatexigens_DSM_14469","Bacteroides_cellulosyliticus_CRE21_DSM_14838","Bacteroides_caccae_DSM_19024","Eubacterium_siraeum_DSM_15702")
wgs_cleaned <- wgs_dist[(colnames(wgs_dist) %in% include),(rownames(wgs_dist) %in% include)]

rm('wgs_dist','include')

# Import metabolomic fold-change distance matrix as published in Han and Van Treuren et al 2021, select only intersect
tax_to_dn <- read.csv('data/from_Han_vanTreuren_2021/taxonomy_to_display_name.txt', sep = "\t", header = T)

include = c("Marvinbryantia formatexigens DSMZ 14469","Bacteroides cellulosilyticus DSMZ 14838","Bacteroides caccae ATCC 43185","Eubacterium siraeum DSMZ 15702")
legend <- tax_to_dn[(tax_to_dn$taxonomy %in% include), ] %>% mutate(taxonomy = str_replace_all(taxonomy," ","_"))

metab_dist <- read.csv('data/from_Han_vanTreuren_2021/metab_dm.txt', sep = '\t', header = T, row.names = 1)
metab_cleaned <- metab_dist[(rownames(metab_dist) %in% legend$disp_name),(colnames(metab_dist) %in% legend$disp_name)]

rm('include','tax_to_dn','metab_dist')

# Replace display names with taxonomies in metabolomics distance matrix
colnames(metab_cleaned) <- legend$taxonomy
rownames(metab_cleaned) <- legend$taxonomy

rm('legend')

# Summarize cryo-EM measurements, calculate distance matrix
vis_sum <- aggregate.data.frame(vis_cleaned[2:6], by = list(vis_cleaned$Strain), FUN = mean)
vis_sum <- vis_sum %>% column_to_rownames(var="Group.1")

cryoem_dist <- as.data.frame(as.matrix(vegdist(vis_sum, method = "euclidean", upper = TRUE)))

rm('vis_cleaned','vis_sum')

# Calculate Mantel statistic to check for correlation
vis_metab_test <- mantel(cryoem_dist,metab_cleaned,method = "pearson")
vis_metab_test

wgs_metab_test <- mantel(wgs_cleaned,metab_cleaned,method = "pearson")
wgs_metab_test

vis_wgs_test <- mantel(cryoem_dist,wgs_cleaned,method = "pearson")
vis_wgs_test

# Make long tables for plotting
lowerTriangle(cryoem_dist) <- NA
cryoem_dist$X <- rownames(cryoem_dist)

cryoem_long <- cryoem_dist %>% 
  as.data.frame() %>%
  pivot_longer(-c(X), names_to = "Y", values_to = "cryoEM_Distance")

lowerTriangle(metab_cleaned) <- NA
metab_cleaned$X <- rownames(metab_cleaned)

metab_long <- metab_cleaned %>% 
  as.data.frame() %>%
  pivot_longer(-c(X), names_to = "Y", values_to = "Metab_Distance")

rm('cryoem_dist','metab_cleaned')

# Merge Tables, add (manually determined) LCL
cryoem_long$Metab_Distance <- metab_long$Metab_Distance
cryoem_long$LCL <- c("Species","Genus","Domain","Domain","Genus","Species","Domain","Domain","Domain","Domain","Species","Order","Domain","Domain","Order","Species")

cryoem_long <- na.omit(cryoem_long)

# Plot
ggplot(cryoem_long, aes(x = cryoEM_Distance, y = Metab_Distance, color = fct_relevel(LCL,"Species","Genus","Order","Domain")))+
  theme_cowplot(12)+
  geom_point(size = 3, na.rm=TRUE, alpha = 0.6)+
  scale_color_brewer(name = 'Shared Taxonomy', palette = "Dark2", direction = 1)+
  geom_smooth(inherit.aes = F, aes(x = cryoEM_Distance, y = Metab_Distance, fill = 'Linear Fit / 95% CI'),method='lm')+
  scale_fill_manual(name = "Linear Fit", values = c("grey80"))+
  xlab("Structural Similarity")+
  ylab("Metabolic Similarity")

cor.test(cryoem_long$cryoEM_Distance, cryoem_long$Metab_Distance, method = "pearson")