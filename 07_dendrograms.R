rm(list = ls())

library('tidyverse')
library('ape')

# Import Scaled Measurements for clustering
vis <- read.csv('results/cryoem_scaled_mini_replicates.csv', header = T, row.names = 1)

# Import WGS Distance Matrix
wgs_dist <- read.csv('results/wgs_dist_cleaned.csv',sep = ",", header = TRUE, row.names = 1, check.names = FALSE)

# Export WGS of replicates
wgs_dist <- as.dist(wgs_dist)
wgs_tree <- as.phylo(hclust(wgs_dist, method = "complete"))

write.tree(wgs_tree,"results/wgs_dist_tree_complete.nwe")

# Exclude complex carbon sources, calculate average values by Sample ID, calculate visual tree
exclude = c("211130-2", "211130-3", "201019-1","201211-1","180628-1")
vis_cleaned <- vis[!(vis$Sample_ID %in% exclude), ]
vis_cleaned <- vis_cleaned %>% mutate(concat = case_when(!is.na(Identifier) ~ paste0(Strain, "_", Identifier),T ~ Strain))

vis_sum <- aggregate.data.frame(vis_cleaned[2:6], by = list(vis_cleaned$concat), FUN = mean)
vis_sum <- vis_sum %>% column_to_rownames(var="Group.1")

vis_dist <- dist(vis_sum, method = "euclidean")
vis_tree <- as.phylo(hclust(vis_dist, method = "complete"))

write.tree(vis_tree,"results/vis_dist_tree_complete.nwe")
