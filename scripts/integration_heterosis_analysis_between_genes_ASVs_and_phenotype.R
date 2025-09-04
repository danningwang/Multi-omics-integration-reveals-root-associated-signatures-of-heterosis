################################################################################
## Multi-omics MPH correlation and network construction
## Author: Danning Wang
## Project: Multi-omics data analysis of heterosis 
################################################################################

library(psych)
library(phyloseq)
library(tidyverse)
library(readxl)
library(igraph)
library(writexl)

setwd('H:/Research Projects/Heterosis/heterosis/')

#-----------------------------------------------#
# 1. Function: Filter MPH data
#-----------------------------------------------#
filter_MPH <- function(MPH_data, het_info, min_hybrids = 6){
  # wide-format
  MPH_wide <- MPH_data %>%
    select(Hybrid, Trait, MPH.pct) %>%
    pivot_wider(names_from = Trait, values_from = MPH.pct)
  
  # Filter traits with sufficient heterosis
  het_traits <- het_info %>% filter(nhb >= min_hybrids) %>% pull(Trait)
  common_cols <- intersect(colnames(MPH_wide), het_traits)
  MPH_wide <- MPH_wide[, c("Hybrid", common_cols)]
  
  return(MPH_wide)
}

#-----------------------------------------------#
# 2. Function: Spearman correlation with FDR
#-----------------------------------------------#
calc_corr <- function(X, Y, method="spearman", adjust="fdr", threshold_rho=0.5, threshold_p=0.05){
  # Ensure same row order
  if(is.data.frame(X) & is.data.frame(Y)){
    stopifnot(nrow(X) == nrow(Y))
  }
  
  cor_res <- corr.test(X, Y, method=method, adjust=adjust, ci=F, normal=F)
  
  cor_mat <- as.data.frame(as.table(cor_res$r))
  cor_p <- as.data.frame(as.table(cor_res$p))
  cor_padj <- as.data.frame(as.table(cor_res$p.adj))
  
  corr_df <- cbind.data.frame(cor_mat, padj = cor_padj$Freq, p = cor_p$Freq)
  colnames(corr_df)[1:3] <- c("Source", "Target", "rho")
  
  # Filter by thresholds
  corr_df_sig <- corr_df[abs(corr_df$rho) >= threshold_rho & corr_df$padj < threshold_p, ]
  return(corr_df_sig)
}

#----------------------------------------------------#
# 3. Function: Create network from correlation
#----------------------------------------------------#
create_network <- function(edges_df, taxa_table=NULL){
  # Ensure edges have Source, Target, rho
  stopifnot(all(c("Source","Target","rho") %in% colnames(edges_df)))
  
  # Create igraph object
  g <- graph_from_data_frame(edges_df, directed=FALSE)
  
  # Add edge weight
  E(g)$weight <- abs(edges_df$rho)
  
  # Extract node info
  nodes <- data.frame(
    ID = V(g)$name,
    hubScore = hub_score(g)$vector,
    Degree = degree(g),
    closeness = closeness(g))
  
  nodes$type <- substr(nodes$ID, 1, 4)
  
  # Add taxonomy info if provided
  if(!is.null(taxa_table)){
    
    nodes <- left_join(nodes, taxa_table, by = c("ID"="row.names"))
  }
  
  return(list(graph = g, nodes = nodes))
}

#--------------------------------------------------------------#
# 4. Correlations between MPH of phenotype and bacteria
#--------------------------------------------------------------#

# Load data
MPH_bac <- read.csv('./results/MPH_correlation/CK/bac_ASVs_MPH_percent_CK.csv')
het_bac <- read_xlsx('./heterosis/results/bacteria/heterosis_test/bac_vst_CK_num_hybrids.xlsx')

pheno_MPH <- read.csv('./heterosis/results/phenotype/CK/MPH/all_traits_MPH_value.csv')

# Filter MPH data
MPH_bac_wide <- filter_MPH(MPH_bac, het_bac, min_hybrids = 6)

# Match Hybrid order
MPH_bac_wide <- MPH_bac_wide[match(pheno_MPH$Hybrid, MPH_bac_wide$Hybrid), ]

# Correlation phenotype vs ASVs
corr_bac_pheno <- calc_corr(MPH_bac_wide[,-1], pheno_MPH[, -1])

# Save results
write_xlsx(corr_bac_pheno, './results/MPH_correlation/CK/corr_bac_pheno.xlsx')

#--------------------------------------------------------------#
# 5. Correlations between MPH of phenotype and fungi
#--------------------------------------------------------------#

# Load data
MPH_fungi <- read.csv('./heterosis/results/MPH_correlation/CK/fungi_ASVs_MPH_percent_CK.csv')
het_fungi <- read_xlsx('./heterosis/results/fungi/heterosis_test/fungi_vst_CK_num_hybrids.xlsx')

#pheno_MPH <- read.csv('./heterosis/results/phenotype/CK/MPH/all_traits_MPH_value.csv')

# Filter MPH data
MPH_fungi_wide <- filter_MPH(MPH_fungi, het_fungi, min_hybrids = 6)

# Match Hybrid order
MPH_fungi_wide <- MPH_fungi_wide[match(pheno_MPH$Hybrid, MPH_fungi_wide$Hybrid), ]

# Correlation phenotype vs ASVs
corr_fungi_pheno <- calc_corr(MPH_fungi_wide[,-1], pheno_MPH[, -1])

# Save results
write_xlsx(corr_fungi_pheno, './results/MPH_correlation/CK/corr_fungi_pheno.xlsx')

#----------------------------------------------------------#
# 6. Correlations between root genes and phenotypes
#----------------------------------------------------------#

MPH_gene <- read.csv('./results/MPH_correlation/CK/Genes_MPH_percent_CK.csv')
het_gene <- read_xlsx('./results/rna_seq/heterotic_genes/CK/heterotic_genes_CK.xlsx')

MPH_gene_wide <- filter_MPH(MPH_gene, het_gene, min_hybrids = 6)

# Match Hybrid order
MPH_gene_wide <- MPH_gene_wide[match(pheno_MPH$Hybrid, MPH_gene_wide$Hybrid), ]

corr_gene_pheno <- calc_corr(MPH_gene_wide[, -1], pheno_MPH[, -1])
write_xlsx(corr_gene_pheno, './results/MPH_correlation/CK/corr_gene_pheno.xlsx')

#-----------------------------------------------#
# 7. Construct network
#-----------------------------------------------#
edges_all <- bind_rows(corr_bac_pheno, corr_fungi_pheno, corr_gene_pheno)
net <- create_network(edges_all)
all.net <- net$graph
nodes_df <- net$nodes

write_xlsx(nodes_df, './results/MPH_correlation/CK/network_nodes.xlsx')
write_xlsx(edges_all, './results/MPH_correlation/CK/network_edges.xlsx')













