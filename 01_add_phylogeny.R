rm(list = ls())

library('tidyverse')
library('taxize')
library('stringr')

# Import Metadata
metadata_cryoem <- read.csv('data/cryoem_metadata.csv', sep = ",", header = T, row.names = 1)

# Extract Genus from name
metadata_cryoem$Genus <- as.character(map(strsplit(metadata_cryoem$Strain, split = "_"),1))

# Extract Species from name
metadata_cryoem$Species <- str_c(str_split_fixed(metadata_cryoem$Strain, "_", 6)[,1], "_", str_split_fixed(metadata_cryoem$Strain, "_", 6)[,2])
  
# Pull background info from NCBI
class_object = classification(metadata_cryoem$Genus,db="ncbi")

# Inserting the Information
metadata_cryoem$Family = class_object %>% map(1) %>% map(7)
metadata_cryoem$Order = class_object %>% map(1) %>% map(6)
metadata_cryoem$Class = class_object %>% map(1) %>% map(5)
metadata_cryoem$Phylum = class_object %>% map(1) %>% map(4)
metadata_cryoem$Kingdom = class_object %>% map(1) %>% map(3)
metadata_cryoem$Domain = class_object %>% map(1) %>% map(2)

# Turn contents to proper text:
metadata_cryoem_new <- apply(metadata_cryoem,2,as.character)

# Bring back row names:
row.names(metadata_cryoem_new) <- row.names(metadata_cryoem)

# Save:
write.csv(metadata_cryoem_new,"results/cryoem_metadata_enriched.csv")
