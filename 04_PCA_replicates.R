rm(list = ls())

library('tidyverse')
library('vegan')

# Re-Read scaled data, separate measurements from metadata, remove measurements for PM - OM:
scaled.mini_repl <- read.csv('results/cryoem_scaled_mini_replicates.csv', sep = ",", header = T, row.names = 1)

# Calculate PCA
mini_pca <- prcomp(scaled.mini_repl[2:6], center = FALSE, scale. = FALSE)
mini_pca$results <- as.data.frame(mini_pca$x)
mini_pca$results$Sample_ID <- scaled.mini_repl$Sample_ID 

# Scree Plot / Loadings
summary(mini_pca)

var_explained_df <- data.frame(PC= paste0("PC",1:5),
                               var_explained=(mini_pca$sdev)^2/sum((mini_pca$sdev)^2))

var_explained_df %>%
  ggplot(aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4)+
  ylab('Variation Explained')+
  geom_line()+
  labs(title="Scree plot: PCA")

rm(var_explained_df)

biplot(mini_pca, scale = 1, pc.biplot = TRUE)
