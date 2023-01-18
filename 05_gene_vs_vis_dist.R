rm(list = ls())

library('tidyverse')
library('vegan')
library('gdata')

# Re-Read scaled data, and WGS Distance.
scaled.mini_repl <- read.csv('results/cryoem_scaled_mini_replicates.csv', sep = ",", header = T, row.names = 1) %>%
  mutate(concat = case_when(!is.na(Identifier) ~ paste0(Strain, "_", Identifier),T ~ Strain))
wgs_dist <- read.csv('data/wgs_distance_jaccard.csv',sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
meta <- read.csv('results/cryoem_metadata_enriched.csv',sep = ",", header = TRUE, row.names = 1) %>%
  mutate(concat = case_when(!is.na(Identifier) ~ paste0(Strain, "_", Identifier),T ~ Strain))

# Check that all samples, for which cryo-EM measurements exist, are part of WGS table
setdiff(scaled.mini_repl$concat,rownames(wgs_dist))
setdiff(scaled.mini_repl$concat,colnames(wgs_dist))

# Remove distance measurements from WGS for which no 3D cryo-EM data exists
row.names.remove <- setdiff(rownames(wgs_dist),scaled.mini_repl$concat)
wgs_dist_cleaned <- wgs_dist[!(row.names(wgs_dist) %in% row.names.remove), ]
wgs_dist_cleaned <- wgs_dist_cleaned[,!names(wgs_dist_cleaned) %in% row.names.remove]

write.csv(wgs_dist_cleaned,"results/wgs_dist_cleaned.csv")

rm('wgs_dist','row.names.remove')

# Sort alphabetically, use only upper triangle of WGS matrix
wgs_dist_cleaned <- wgs_dist_cleaned[order(rownames(wgs_dist_cleaned)), ]
wgs_dist_cleaned <- wgs_dist_cleaned[ ,order(colnames(wgs_dist_cleaned))]

lowerTriangle(wgs_dist_cleaned) <- NA

# Calculate replicate-level cryo-EM distance, sort, retain only upper triangle
cryoem_dist <- as.data.frame(as.matrix(vegdist(scaled.mini_repl[2:6], method = "euclidean", upper = TRUE)))

colnames(cryoem_dist) <- scaled.mini_repl$concat

cryoem_dist_sorted <- cryoem_dist[ ,order(colnames(cryoem_dist))]
colnames(cryoem_dist_sorted) <- str_split_fixed(colnames(cryoem_dist_sorted), '\\.', 2)[,1]

cryoem_dist_sorted$X <- scaled.mini_repl$concat
cryoem_dist_sorted <- cryoem_dist_sorted[order(cryoem_dist_sorted$X),]

write.csv(cryoem_dist,'results/cryoem_dist_replicates.csv')

lowerTriangle(cryoem_dist_sorted[1:225]) <- NA

rm('cryoem_dist')

# Turn into long tables
wgs_long <- wgs_dist_cleaned %>% 
  as.data.frame() %>%
  rownames_to_column("X") %>%
  pivot_longer(-c(X), names_to = "Y", values_to = "WGS_Distance")
wgs_long$comp <- paste(wgs_long$X,wgs_long$Y)
wgs_long <- na.omit(wgs_long)

cryoem_long <- cryoem_dist_sorted %>% 
  as.data.frame() %>%
  pivot_longer(-c(X), names_to = "Y", values_to = "cryoEM_Distance")
cryoem_long$comp <- paste(cryoem_long$X,cryoem_long$Y)
cryoem_long <- na.omit(cryoem_long)

rm('wgs_dist_cleaned','cryoem_dist_sorted','scaled.mini_repl')

# Introduce information on lowest common level (LCL)
find_lcl <- function(lhs, rhs, df) {
  # Define Phylogeny Vectors as [Genus, Family, Order, Class, Phylum, Kindom, Domain]
  
  ret <- c('Domain','Kingdom','Phylum','Class','Order','Family','Genus','Species')
  
  df_lhs <- df[df$concat == lhs,]
  df_rhs <- df[df$concat == rhs,]
  
  v_lhs <- c(df_lhs$Species,df_lhs$Genus,df_lhs$Family,df_lhs$Order,df_lhs$Class,df_lhs$Phylum,df_lhs$Kingdom,df_lhs$Domain)
  v_rhs <- c(df_rhs$Species,df_rhs$Genus,df_rhs$Family,df_rhs$Order,df_rhs$Class,df_rhs$Phylum,df_rhs$Kingdom,df_rhs$Domain)
  
  if(length(intersect(v_lhs,v_rhs)) > 0) {
    lcl <- c(ret[length(intersect(v_lhs,v_rhs))])
  } else {
    lcl <- c('>Domain')
  }
  
  
  return(lcl)
} 

wgs_long$LCL <- NA

for (i in 1:length(wgs_long$X)) {
  wgs_long$LCL[i] <- find_lcl(wgs_long$X[i],wgs_long$Y[i],meta)
}

rm('meta','i','find_lcl')

# Merge tables
df <- merge(cryoem_long, wgs_long, by = 'comp', all.x = TRUE)
df <- subset(df, select = -c(2,3,5,6))
df <- na.omit(df)

write.csv(df,'results/distance_vis_vs_genome_vs_LCL.csv')

# Plot Vis vs. WGS
ggplot(df, aes(WGS_Distance, cryoEM_Distance, color = fct_relevel(LCL,"Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain")))+
  geom_point(size = 1, alpha = 0.5)+
  scale_fill_manual(values = c('#000000'))+
  theme(legend.title = element_blank())+
  xlab('Genome Distance')+
  ylab('Visual Distance')

# Plot Vis vs. LCL
ggplot(df, aes(fct_inorder(LCL), cryoEM_Distance, fill = fct_relevel(LCL,"Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain")))+
  geom_boxplot(outlier.size=1, outlier.alpha = 0.5, outlier.fill = NULL, outlier.color = NULL)+
  scale_x_discrete(limits = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain"))+
  scale_fill_brewer(type = "seq", direction = 1)+
  xlab('Last Branching')+
  ylab('Visual Distance')+  
  theme(legend.position="none")

# Plot WGS vs LCL
ggplot(df, aes(fct_inorder(LCL), WGS_Distance, fill = fct_relevel(LCL,"Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain")))+
  geom_boxplot(outlier.size=1, outlier.alpha = 0.5, outlier.fill = NULL, outlier.color = NULL)+
  scale_x_discrete(limits = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain",">Domain"))+
  scale_fill_brewer(type = "seq", direction = 1)+
  xlab('Last Branching')+
  ylab('Genomic Distance')+  
  theme(legend.position="none")

# Calculate ANOVA Visual Dist vs. LCL
sign <- aov(cryoEM_Distance ~ LCL, data = df)
summary(sign)

# Calculate correlation WGS / cryoEM for all LCL Species / Genus / Family
keep <- c("Species","Genus","Family")
species_genus_family <- df[(df$LCL %in% keep),]
cor.test(species_genus_family$cryoEM_Distance, species_genus_family$WGS_Distance, method = "pearson")

other <- df[!(df$LCL %in% keep),]
cor.test(other$cryoEM_Distance, other$WGS_Distance, method = "pearson")
