rm(list = ls())

library('tidyverse')
library('vegan')
library('parallel')

# Re-Read scaled data, separate measurements from metadata, remove measurements for PM - OM:
scaled.mini_repl <- read.csv('results/cryoem_scaled_mini_replicates.csv', sep = ",", header = T, row.names = 1)

# Calculate overall ANOSIM by Sample
overall <- anosim(scaled.mini_repl[2:6], scaled.mini_repl$Sample, permutations = 9999, distance = "euclidean", strata = NULL, parallel = 2)
summary(overall)

# Calculate Pairwise ANOSIM for all Genera which have multiple examples in dataset
df <- scaled.mini_repl %>% group_by(Family) %>% filter(n() > 5)

genera <- unique(df$Family)

pairwise <- as.data.frame(matrix(data = NA, ncol = 4, nrow = 0))
colnames(pairwise) <- c("X","Y","ANOSIM_R","signif")

counter = 1
i = 1

for (i in 1:(length(genera)-1)) {
  # Print appropriate Genus name
  print(genera[i])
  
  # Generate counter for all outstanding objects
  counter = i + 1
  for (j in counter:length(genera)) {
    # Subset dataframe including these two genera
    temp <- subset(scaled.mini_repl, (Family == genera[i] | Family == genera[j]))
    
    # Calculate ANOSIM on this pair
    temp_anosim <- anosim(temp[2:6], temp$Family, permutations = 9999, distance = "euclidean", strata = NULL, parallel = 2)
    
    # Save to table
    pairwise[nrow(pairwise)+1,] <- c("","","","")
    
    pairwise$X[nrow(pairwise)] <- as.character(genera[i])
    pairwise$Y[nrow(pairwise)] <- as.character(genera[j])
    
    pairwise$ANOSIM_R[nrow(pairwise)] <- as.numeric(temp_anosim$statistic)
    pairwise$signif[nrow(pairwise)] <- as.numeric(temp_anosim$signif)
    
    # Save reverse direction to table
    pairwise[nrow(pairwise)+1,] <- c("","","","")
    pairwise$X[nrow(pairwise)] <- as.character(genera[j])
    pairwise$Y[nrow(pairwise)] <- as.character(genera[i])
    
    pairwise$ANOSIM_R[nrow(pairwise)] <- as.numeric(temp_anosim$statistic)
    pairwise$signif[nrow(pairwise)] <- as.numeric(temp_anosim$signif)

  }
  
}

pairwise$signif <- as.numeric(pairwise$signif)
pairwise$ANOSIM_R <- as.numeric(pairwise$ANOSIM_R)

rm('df','temp','temp_anosim','counter','genera','i','j')

# Remove all nonsignificant values
pairwise$ANOSIM_R <- ifelse(pairwise$signif > 0.05, NA, pairwise$ANOSIM_R)

write.csv(pairwise,'results/anosim_pairwise_family_long_replicates.csv')

# plot as heatmap  
ggplot(pairwise, aes(X, Y))+
  geom_tile(aes(fill = ANOSIM_R))+
  scale_fill_fermenter(palette = "BuPu", direction = 1, na.value = "#b3b3b3")+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  xlab(element_blank())+
  ylab(element_blank())+
  theme(legend.title = element_blank())
