rm(list = ls())

library('tidyverse')

# Read replicate-level data, merge with metadata:
mini_repl <- read.csv('data/cryoem_measurements_replicates.csv', sep = ",", header = T)
meta <- read.csv('results/cryoem_metadata_enriched.csv', sep = ",", header = T, row.names = 1)

mini_meta <- merge(mini_repl, meta, by.x = 'Sample_ID', by.y = 'row.names')

rm('mini_repl','meta')

# Scale and center data
scaled.mini_repl <- mini_meta
scaled.mini_repl <- mini_meta %>% mutate_if(is.numeric,scale)

write.csv(scaled.mini_repl,"results/cryoem_scaled_mini_replicates.csv")
