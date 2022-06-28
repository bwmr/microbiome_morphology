rm(list = ls())

library('tidyverse')
library('vegan')
library('gdata')

# Import Data, exclude multiple carbon sources, retain only those datasets where a minimum number of samples per Family etc. is present
scaled.mini_repl <- read.csv('results/cryoem_scaled_mini_replicates.csv', sep = ",", header = T, row.names = 1) %>%
  mutate(concat = case_when(!is.na(Identifier) ~ paste0(Strain, "_", Identifier),T ~ Strain))

exclude = c("211130-2", "211130-3", "201019-1","201211-1","180628-1")
vis_cleaned <- scaled.mini_repl[!(scaled.mini_repl$Sample_ID %in% exclude), ] 
vis_cleaned <- vis_cleaned %>% group_by(Family) %>% filter(n() > 5)

# Calculate Distances
cryoem_dist <- as.data.frame(as.matrix(vegdist(vis_cleaned[2:6], method = "euclidean", upper = TRUE)))

colnames(cryoem_dist) <- vis_cleaned$concat

cryoem_dist_sorted <- cryoem_dist[ ,order(colnames(cryoem_dist))]
colnames(cryoem_dist_sorted) <- str_split_fixed(colnames(cryoem_dist_sorted), '\\.', 2)[,1]

cryoem_dist_sorted$X <- vis_cleaned$concat
cryoem_dist_sorted <- cryoem_dist_sorted[order(cryoem_dist_sorted$X),]

rm('cryoem_dist','scaled.mini_repl','exclude')

# Turn into long table
cryoem_long <- cryoem_dist_sorted %>% 
  as.data.frame() %>%
  pivot_longer(-c(X), names_to = "Y", values_to = "cryoEM_Distance")

rm('cryoem_dist_sorted')

# Add Family_X / Family_Y columns
info <- vis_cleaned %>% select(concat, Family) %>% distinct(.keep_all = TRUE)

cryoem_extended <- merge(cryoem_long, info, by.x = 'X', by.y = 'concat', all.x = TRUE, all.y = FALSE) %>% rename(Family_X = Family)
cryoem_extended <- merge(cryoem_extended, info, by.x = 'Y', by.y = 'concat', all.x = TRUE, all.y = FALSE) %>% rename(Family_Y = Family)

rm('cryoem_long','info','vis_cleaned')

# Calculate average distances within each Family, average distance between families

todo <- unique(cryoem_extended$Family_X)

distance <- as.data.frame(matrix(data = NA, ncol = 3, nrow = 0))
colnames(distance) <- c("X","Y","Mean")

for (i in 1:(length(todo))) {

  # Subset table for Group_X = group name
  temp <- subset(cryoem_extended, (Family_X == todo[i]))
  
  # Summarize by Group_Y, make data frame, also add inverted
  temp_sum <- aggregate.data.frame(temp$cryoEM_Distance, by = list(temp$Family_Y), FUN = mean) %>% rename(Y = Group.1, Mean = x)
  temp_sum$X <- c(todo[i],todo[i],todo[i],todo[i],todo[i],todo[i],todo[i],todo[i])
  distance <- rbind(distance,temp_sum)
  
}

distance$Mean <- as.numeric(distance$Mean)

rm('temp','temp_sum','todo','i')

# Pivot to wide, remove lower triangle
distance_wide <- as.data.frame(pivot_wider(distance, id_cols = Y, names_from = X,values_from = Mean)) 

distance_wide <- column_to_rownames(distance_wide, var = "Y")

distance_wide <- distance_wide[order(rownames(distance_wide)), ]
distance_wide <- distance_wide[ ,order(colnames(distance_wide))]


lowerTriangle(distance_wide, diag = FALSE) <- NA

# Pivot to long
distance_long <- distance_wide %>% 
  as.data.frame() %>%
  rownames_to_column("X") %>%
  pivot_longer(-c(X), names_to = "Y", values_to = "Mean") %>%
  na.omit

rm('distance','distance_wide')

# Export
write.csv(distance_long,'results/distance_within_vs_between_families_upper.csv')

# Extract Butyrivibrio fibrosolvens and Agathobacter ruminis data

todo <- c("Butyrivibrio_fibrisolvens_CF3","Agathobacter_ruminis_DSM_29029")

distance_bf_ar <- cryoem_extended[(cryoem_extended$X %in% todo), ]
distance_bf_ar <- distance_bf_ar[(distance_bf_ar$Y %in% todo), ]

mean_dist_bf_ar <- distance_bf_ar %>%
  group_by(X,Y) %>%
  dplyr::summarize(Mean = mean(cryoEM_Distance)) %>%
  as.data.frame()

rm('cryoem_extended','distance_bf_ar','todo')

# Remove Bf / Ar duplicate
mean_dist_bf_ar <- mean_dist_bf_ar[-c(3),]

write.csv(mean_dist_bf_ar,'results/distance_between_within_Ar_Bf.csv')
