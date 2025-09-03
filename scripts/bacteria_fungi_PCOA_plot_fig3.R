################################################################################
# Analysis of bacterial/fungal community variation in hybrids vs inbreds
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


mytheme = theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), 
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 8)) 



# --- Setup ---
setwd('H:/Research Projects/Heterosis/heterosis/')

# Load data
bac.vst.ck.corr <- readRDS("./intermediate_data/bacteria/bac.LN.vst.hidden.batch.corrected.RDS")
meta <- as.data.frame(sample_data(bac.vst.ck.corr))

# ============================================================
# 1. Beta-dispersion (distance to centroid) test
# ============================================================

# Calculate Bray-Curtis distance
bac.dist <- phyloseq::distance(bac.vst.ck.corr, method = "bray")

# Test differences in dispersion (homogeneity of variance)
dispersion <- betadisper(bac.dist, meta$hybrid, type = "centroid")
anova(dispersion)  


# ============================================================
# 2. PERMANOVA (adonis2)
# ============================================================

p.test <- adonis2(as.matrix(bac.dist) ~ hybrid, data = meta, permutations = 999)
print(p.test)

# ============================================================
# 3. PCoA ordination
# ============================================================

PCoA.bray <- ordinate(bac.vst.ck.corr, method = "PCoA", distance = "bray")

# Extract eigenvectors & metadata
pcoa.df <- cbind.data.frame(PCoA.bray$vectors, meta)

# Percent variance explained
pct <- PCoA.bray$values$Relative_eig[1:2]  # first 2 axes
percentage <- paste0(c("PCoA1", "PCoA2"), " (", round(pct * 100, 2), "%)")

# Factor order
pcoa.df$hybrid <- factor(pcoa.df$hybrid, levels = c("hybrid", "inbred"))

# ============================================================
# 4. MANOVA test on all axes
# ============================================================

manova_res <- manova(as.matrix(pcoa.df[, grep("^Axis", names(pcoa.df))]) ~ hybrid, data = pcoa.df)
summary(manova_res)

# =========================================================
# 5. COLOR SETTINGS (inbreds vs hybrids)
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


# ============================================================
# 6. Plot PCoA
# ============================================================

# Centroids for each genotype
centroids <- pcoa.df %>%
  group_by(genotype, hybrid) %>%
  summarise(PC1 = mean(Axis.1), PC2 = mean(Axis.2), .groups = "drop")

# Convex hulls
find_hull <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
hulls <- pcoa.df %>% group_by(genotype) %>% do(find_hull(.))

# Final plot
pcoa.pl <- ggplot(pcoa.df, aes(x = Axis.1, y = Axis.2, color = genotype, shape = hybrid)) +
  geom_point(size = 3) +
  geom_polygon(data = hulls, aes(fill = genotype), alpha = 0.4) +
  scale_fill_manual(values = c(hybrid_colors, inbred_colors)) +
  scale_color_manual(values = c(hybrid_colors, inbred_colors)) +
  xlab(percentage[1]) + ylab(percentage[2]) +
  coord_fixed(ratio = 1) +
  mytheme

ggsave(pcoa.pl, file = "./results/plots/PCA and PCOA/bac-PCoA-plot-LN.pdf", 
       width = 20, height = 15, units = "cm")

## repeat above steps for other soil conditions and fungi



