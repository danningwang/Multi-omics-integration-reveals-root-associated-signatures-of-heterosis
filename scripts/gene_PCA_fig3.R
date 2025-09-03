################################################################################
# RNA-seq: PCA, Dispersion, PERMANOVA, MANOVA, and PCA Plot
# Author: Danning Wang
# Project: Multi-omics data analysis for heterosis 
################################################################################

library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(DESeq2)
library(vegan)
library(phyloseq)
library(dplyr)
library(reshape2)


setwd('H:/Research Projects/Heterosis/')

mytheme = theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), 
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 8)) 


# --- Load metadata and gene expression data ---
md <- read.table('./intermediate_data/rna/rna_sample_metadata.txt', header = TRUE)

# Subset only CK treatment samples, and repeatly do this for other soil
md.ck <- md[md$Treatment == 'CK', ]

# Load variance-stabilized normalized gene expression table
gene.tab <- readRDS('./intermediate_data/rna/filtered_gene_table_vst_normlized_CK_corrected.RDS')


# =========================================================
# 1. DISPERSION ANALYSIS
# =========================================================

# Calculate Euclidean distance matrix between samples
gene.matrix <- dist(t(gene.tab), method = "euclidean")

# Test dispersion (betadisper) between hybrids and inbreds
dispersion <- betadisper(gene.matrix, md.ck$Hybrid_inbred, type = 'centroid')
anova(dispersion)  # tests if group dispersions differ significantly


# =========================================================
# 2. PERMANOVA TEST
# =========================================================

# Tests whether centroids of hybrid vs inbred groups differ
p.test <- adonis2(gene.matrix ~ Hybrid_inbred, data = md.ck, permutations = 999)
p.test


# =========================================================
# 3. PCA ANALYSIS
# =========================================================

# Run PCA on transposed gene expression matrix
pcaDat <- prcomp(t(gene.tab))

# Screeplot of first 10 PCs
screeplot(pcaDat, type = "l", npcs = 10, main = "Screeplot of the first 10 PCs")

# Calculate explained variance (%) for each PC
percentage <- round(pcaDat$sdev^2 / sum(pcaDat$sdev^2) * 100, 2)  
percentage <- paste(colnames(pcaDat$x), "(", percentage, "%)", sep = "")

# Convert PCA result into dataframe
pcs.df <- as.data.frame(pcaDat$x)
pcs.df$Genotype <- md.ck$Genotype
pcs.df$Hybrid   <- md.ck$Hybrid_inbred


# =========================================================
# 4. MANOVA TEST
# =========================================================

# Multivariate test of differences across all PCs for hybrids vs inbred lines
manova_res <- manova(as.matrix(pcs.df[, 1:83]) ~ Hybrid, data = pcs.df)
summary(manova_res)  # Pillai’s trace is a robust test


# =========================================================
# 5. CALCULATE CENTROIDS (group averages in PCA space)
# =========================================================

centroids <- pcs.df %>%
  group_by(Genotype, Hybrid) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), .groups = "drop")


# =========================================================
# 6. COLOR SETTINGS (inbreds vs hybrids)
# =========================================================

# 7 distinct colors for inbreds
inbred_colors <- brewer.pal(7, "Dark2")

# Hybrids mapped to lighter shades of parental colors
hybrid_colors <- c(
  alpha(inbred_colors[1], 0.8), alpha(inbred_colors[1], 0.5), alpha(inbred_colors[1], 0.3),
  alpha(inbred_colors[2], 0.8), alpha(inbred_colors[2], 0.5), alpha(inbred_colors[2], 0.3),
  alpha(inbred_colors[3], 0.8), alpha(inbred_colors[3], 0.5), alpha(inbred_colors[3], 0.3),
  alpha(inbred_colors[4], 0.8), alpha(inbred_colors[4], 0.7), alpha(inbred_colors[4], 0.6),
  alpha(inbred_colors[4], 0.5), alpha(inbred_colors[4], 0.4), alpha(inbred_colors[4], 0.3)
)


# =========================================================
# 7. PCA PLOT WITH CONVEX HULLS
# =========================================================

# Function to compute convex hull per group
find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

# Apply hull function by genotype
hulls <- pcs.df %>% group_by(Genotype) %>% do(find_hull(.))

# PCA plot
pca.pl <- ggplot(pcs.df, aes(PC1, PC2, color = Genotype, shape = Hybrid)) +
  geom_point(size = 3) +
  geom_polygon(data = hulls, aes(fill = Genotype), alpha = 0.4) +  # convex hulls
  scale_fill_manual(values = c(hybrid_colors, inbred_colors)) +
  scale_color_manual(values = c(hybrid_colors, inbred_colors)) +
  geom_text(data = centroids, aes(label = Genotype), color = 'black', size = 2.5) +
  xlab(percentage[1]) + ylab(percentage[2]) +
  coord_fixed(ratio = 1) +
  mytheme

# Save PCA plot
ggsave(pca.pl, file = "./results/plots/PCA_and_PCOA/gene-PCA-plot-CK.pdf", 
       width = 20, height = 15, units = "cm")



