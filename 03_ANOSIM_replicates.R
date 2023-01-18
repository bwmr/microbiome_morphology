rm(list = ls())

library('tidyverse')
library('vegan')
library('parallel')

# Re-Read scaled data, separate measurements from metadata, remove measurements for PM - OM:
scaled.mini_repl <- read.csv('results/cryoem_scaled_mini_replicates.csv', sep = ",", header = T, row.names = 1)

# Calculate ANOSIM at all levels
levels = c("Sample_ID","Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain")

anosim_result <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(anosim_result) <- c("Grouping","ANOSIM R","p")

for (level in levels) {
  cmd <- sprintf("anosim(scaled.mini_repl[2:6], scaled.mini_repl$%s, permutations = 9999, distance = 'euclidean', strata = NULL, parallel = 4)",level)
  result <- eval(parse(text=cmd))
  anosim_result[nrow(anosim_result) + 1,] = c(level,result$statistic, result$signif)
}

write.csv(anosim_result,"results/anosim_by_grouping_new.csv")
