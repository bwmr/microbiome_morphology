rm(list = ls())

library('tidyverse')
library('ape')

# Import Scaled Measurements
vis <- read.csv('results/cryoem_scaled_mini_replicates.csv', header = T, row.names = 1)

# Import WGS Distance Matrix
wgs_dist <- read.csv('results/wgs_dist_cleaned.csv',sep = ",", header = TRUE, row.names = 1, check.names = FALSE)

# Exclude complex carbon sources, calculate average values by Sample ID, calculate visual tree
exclude = c("211130-2", "211130-3", "201019-1","201211-1","180628-1")
vis_cleaned <- vis[!(vis$Sample_ID %in% exclude), ]
vis_cleaned <- vis[vis$Family=='Lachnospiraceae',]
vis_cleaned <- vis_cleaned %>% mutate(concat = case_when(!is.na(Identifier) ~ paste0(Strain, "_", Identifier),T ~ Strain))

vis_sum <- aggregate.data.frame(vis_cleaned[2:6], by = list(vis_cleaned$concat), FUN = mean)
vis_sum <- vis_sum %>% column_to_rownames(var="Group.1")

# Filter & Export WGS of replicates
wgs_cleaned <- wgs_dist[(rownames(wgs_dist) %in% rownames(vis_sum)),(colnames(wgs_dist) %in% rownames(vis_sum))]
wgs_cleaned <- as.dist(wgs_cleaned)
wgs_tree <- as.phylo(hclust(wgs_cleaned, method = "ward.D2"))

write.tree(wgs_tree,"results/wgs_dist_tree_lachno_ward.nwe")
