rm(list = ls())

library('tidyverse')
library('vegan')

# Re-Read scaled data, separate measurements from metadata, remove measurements for PM - OM:
mini_repl <- read.csv('data/cryoem_measurements_replicates.csv', sep = ",", header = T)
meta <- read.csv('results/cryoem_metadata_enriched.csv', sep = ",", header = T, row.names = 1)

mini_meta <- merge(mini_repl, meta, by.x = 'Sample_ID', by.y = 'row.names')

rm('mini_repl','meta')

# Scale data, make scaled dataset
scaled.mini_repl <- mini_meta
scaled.mini_repl <- mini_meta %>% mutate_if(is.numeric,scale)

scaled.mini_repl$OM_Present <- as.character(mini_meta$OM_Present)
scaled.mini_repl[scaled.mini_repl == 'Monoderm'] <- '0'
scaled.mini_repl[scaled.mini_repl == 'Archaeum'] <- '0'
scaled.mini_repl[scaled.mini_repl == 'Diderm'] <- '1'
scaled.mini_repl$OM_Present <- as.numeric(scaled.mini_repl$OM_Present)

write.csv(scaled.mini_repl,"results/cryoem_scaled_mini_replicates.csv")