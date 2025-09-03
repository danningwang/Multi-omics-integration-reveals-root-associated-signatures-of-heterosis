################################################################################
# WGCNA Analysis of fungal community 
# Author: Danning Wang
# Project: Multi-omics data analysis for heterosis 
################################################################################

library(WGCNA)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(writexl)
library(readxl)




# --- Setup --------------------------------------------------------------------

setwd("H:/Research Projects/Heterosis/heterosis/")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Input data
f.ck.corrected <- readRDS("./intermediate_data/fungi/fungi.CK.vst.hidden.batch.corrected.RDS")

#===============================================================================
# 1. Prepare expression data & sample clustering
#===============================================================================

# Convert OTU table to expression matrix (samples x ASVs)
datExpr.fungi <- t(as.data.frame(otu_table(f.ck.corrected)))
dim(datExpr.fungi)

# MPH (mid-parent heterosis) table
MPH.fungi <- read.csv("./results/MPH_correlation/CK/fungi_ASVs_MPH_percent_CK.csv", header = TRUE)

# Reshape to wide format (Hybrids x Traits)
MPH.fungi.wide <- MPH.fungi %>%
  select(Hybrid, Trait, MPH.pct) %>%
  pivot_wider(names_from = Trait, values_from = MPH.pct)

rownames(MPH.fungi.wide) <- MPH.fungi.wide$Hybrid
datExpr.fungi <- MPH.fungi.wide[, -1]   # remove Hybrid column

# Cluster samples
sampleTree <- hclust(dist(datExpr.fungi), method = "average")

pdf("./results/fungiteria/WGCNA/1-sampleClustering.pdf", width = 30, height = 9)
par(cex = 0.8, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", xlab = "", sub = "")
dev.off()

#===============================================================================
# 2. Choose soft-threshold power
#===============================================================================

powers <- c(1:10, seq(12, 30, 2))

sft <- pickSoftThreshold(
  datExpr.fungi, 
  powerVector = powers, 
  verbose = 5, 
  networkType = "unsigned"
)

pdf("./results/fungiteria/WGCNA/2-softThreshold.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))

# Scale-free topology fit index
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale independence", type = "n")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.8, col = "red")

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = "red")
dev.off()

#===============================================================================
# 3. Construct Topological Overlap Matrix (TOM)
#===============================================================================

power <- sft$powerEstimate
TOM <- TOMsimilarityFromExpr(datExpr.fungi, power = power, TOMType = "unsigned")
saveRDS(TOM, "./results/fungi/WGCNA/TOM_power.RDS")

dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

pdf("./results/fungi/WGCNA/3-geneClustering.pdf", width = 12, height = 9)
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

#===============================================================================
# 4. Module identification
#===============================================================================

dynamicMods <- cutreeDynamic(
  dendro = geneTree, 
  distM = dissTOM,
  deepSplit = 4,
  pamRespectsDendro = FALSE, 
  minClusterSize = 5
)

dynamicColors <- labels2colors(dynamicMods)

pdf("./results/fungi/WGCNA/4-moduleTree.pdf", width = 8, height = 6)
plotDendroAndColors(
  geneTree, dynamicColors, "Dynamic Tree Cut",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

#===============================================================================
# 5. Merge close modules
#===============================================================================

merge <- mergeCloseModules(
  datExpr.fungi, dynamicColors, cutHeight = 0.25, verbose = 3
)

mergedColors <- merge$colors
mergedMEs <- merge$newMEs

pdf("./results/fungi/WGCNA/5-mergedModules.pdf", width = 20, height = 9)
plotDendroAndColors(
  geneTree, cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

#===============================================================================
# 6. Export networks for Cytoscape
#===============================================================================

outdir <- "./results/fungi/WGCNA/cytoscape/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (mod in unique(mergedColors)) {
  modOTUs <- colnames(datExpr.fungi)[mergedColors == mod]
  modTOM <- TOM[modOTUs, modOTUs]
  dimnames(modTOM) <- list(modOTUs, modOTUs)
  
  exportNetworkToCytoscape(
    modTOM,
    edgeFile = file.path(outdir, paste0("edges-", mod, ".txt")),
    nodeFile = file.path(outdir, paste0("nodes-", mod, ".txt")),
    weighted = TRUE, threshold = -1,
    nodeNames = modOTUs, nodeAttr = mergedColors[mergedColors == mod]
  )
}

#===============================================================================
# 7. Module eigengenes heatmap
#===============================================================================

rownames(mergedMEs) <- rownames(MPH.fungi.wide)

pdf("./results/fungi/WGCNA/6-moduleEigengenes.pdf", width = 15, height = 40)
pheatmap(mergedMEs,
         cluster_cols = TRUE, cluster_rows = FALSE,
         show_rownames = FALSE, fontsize = 20,
         breaks = seq(-0.5, 0.5, 0.01))
dev.off()

#===============================================================================
# 8. Module-trait correlations
#===============================================================================

pheno.MPH <- readxl::read_xlsx("./results/phenotype/CK/MPH/all_traits_MPH_value_new.xlsx")

moduleTraitCor <- cor(
  mergedMEs[match(pheno.MPH$Hybrid, rownames(mergedMEs)), ],
  pheno.MPH[, -1],
  use = "p", method = "spearman"
)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(mergedMEs))

# Save results
write.table(moduleTraitCor, "./results/fungi/WGCNA/module_trait_correlation.txt", sep = "\t")
write.table(moduleTraitPvalue, "./results/fungi/WGCNA/module_trait_pvalue.txt", sep = "\t")

# Heatmap
pdf("./results/fungi/WGCNA/7-moduleTraitHeatmap.pdf", width = 12, height = 15)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = colnames(mergedMEs),
               ySymbols = colnames(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = "fungi module - trait relationships")
dev.off()



