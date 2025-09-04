################################################################################
# RDA Analysis 
# Author: Danning Wang
# Project: Multi-omics data analysis for heterosis 
################################################################################

# Load packages
library(vegan)
library(ggplot2)
library(readxl)
library(dplyr)
library(corrplot)
library(psych)


mytheme = theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), 
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 8)) 


setwd('H:/Research Projects/Heterosis')

#==============================================
# 1. Load and preprocess data
#==============================================

# --- Transcriptome data ---
ck.gene.tab <- readRDS("./intermediate_data/rna/filtered_gene_table_vst_normlized_CK_corrected.RDS")
ck.gene.tab <- t(as.data.frame(ck.gene.tab))

md <- read.table("./intermediate_data/rna/rna_sample_metadata_correct.txt", header = TRUE)
md.ck <- md[md$Treatment == "CK", ]

stopifnot(all(rownames(ck.gene.tab) == md.ck$SampleID))
ck.gene.md <- cbind.data.frame(ck.gene.tab, md.ck)

# Aggregate gene expression per genotype
ck.gene.mean <- aggregate(. ~ Genotype, data = ck.gene.md[, c(1:32047, 32049)], mean)


# --- Bacteria data ---
ck.bac <- readRDS("./intermediate_data/bacteria/bac.ck.vst.hidden.batch.corrected.RDS")
ck.bac.tab <- t(as.data.frame(otu_table(ck.bac)))
ck.bac.tab <- cbind.data.frame(ck.bac.tab, sample_data(ck.bac)[, "genotype"])
ck.bac.mean <- aggregate(. ~ genotype, data = ck.bac.tab, mean)


# --- Fungi data ---
ck.fungi <- readRDS("./intermediate_data/fungi/fungi.CK.vst.hidden.batch.corrected.RDS")
ck.fungi.tab <- t(as.data.frame(otu_table(ck.fungi)))
ck.fungi.tab <- cbind.data.frame(ck.fungi.tab, sample_data(ck.fungi)[, "genotype"])
ck.fungi.mean <- aggregate(. ~ genotype, data = ck.fungi.tab, mean)


# --- Phenotype data ---
pheno <- read.table("./results/phenotype/all_phenotpe_traits_mean.txt", header = TRUE)
# Standardize phenotype data (exclude IDs)
pheno.scaled <- apply(pheno[, -c(1,2)], 2, scale)
pheno.scaled <- cbind.data.frame(pheno[, 1:2], pheno.scaled)


#==========================
# 2. Fit RDA model
#==========================

# Check alignment between genotype and phenotype
stopifnot(all(ck.gene.mean$genotype == pheno.scaled$Genotype))

# Full model with all phenotype variables, repeat this for each omic data
gene.rda <- rda(ck.gene.mean[, -1] ~ ., data = pheno.scaled[, c(2:5, 7, 15:16, 21:22)])

# Forward selection
fwd.sel <- ordiR2step(
  rda(ck.gene.mean[, -1] ~ 1, data = pheno.scaled[, c(2:5, 7, 15:16, 21:22)]),
  scope = formula(gene.rda),
  direction = "forward",
  R2scope = TRUE,
  pstep = 1000,
  trace = FALSE
)
print(fwd.sel$call)


#================================
# 3. Evaluate model
#================================

gene.rda.signif <- rda(ck.gene.mean[, -1] ~ ., data = pheno.scaled[, c(2:5, 7, 15:16, 21:22)])

# Adjusted R² and ANOVA
print(RsquareAdj(gene.rda.signif)$adj.r.squared)
anova.cca(gene.rda.signif, step = 1000)               # global test
anova.cca(gene.rda.signif, step = 1000, by = "term")  # per-variable test


#============================================
# 4. ggplot RDA plot with arrows
#============================================

# Extract sample scores
rda.df <- scores(gene.rda.signif, display = "sites", scaling = 1) %>% as.data.frame()
rda.df$name <- ck.gene.mean$Genotype
rda.df$type <- factor(ifelse(grepl("x", rda.df$name), "Hybrid", "Inbred"), levels = c("Inbred", "Hybrid"))

# Extract biplot (arrows)
pl.arrows <- data.frame(RDA1 = gene.rda.signif$CCA$biplot[, 1],
                        RDA2 = gene.rda.signif$CCA$biplot[, 2],
                        label = rownames(gene.rda.signif$CCA$biplot))

# ggplot
f4 <- ggplot(rda.df, aes(x = RDA1, y = RDA2, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#8498AB", "#EDC66A")) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_segment(data = pl.arrows, aes(x = 0, y = 0, xend = 2*RDA1, yend = 2*RDA2),
               arrow = arrow(type = "open", length = unit(0.2, "cm")),
               color = "royalblue", size = 0.5) +
  geom_text(data = pl.arrows, aes(x = 2*RDA1, y = 2*RDA2, label = label),
            color = "royalblue", size = 3, vjust = -0.5) +
  ggtitle("Transcriptome RDA (CK)") +
  coord_fixed() +
  mytheme

ggsave(f4, file = "./results/RDA/CK_gene_phenotype_RDA_plot.pdf",
       width = 12, height = 8, units = "cm")


################################################################################
##  Plot of MPH of phenotype and correlation heatmap
################################################################################

#=======================================
# 1. Load and preprocess data
#=======================================

# Read in anatomy dataset (sheet 2 of the Excel file)
anatomy <- read_xlsx("./data/heterosis root triats clean.xlsx", sheet = 2) %>%
  as.data.frame()

# Root anatomy traits
vars.names.root <- colnames(anatomy)[4:6]

# Identify hybrids
hybrids <- unique(anatomy$Genotype[grep("x", anatomy$Genotype)])

#=================================================
# 2. Calculate heterosis (MPH) for each zone
#=================================================

het.res.anatomy <- do.call(rbind, lapply(unique(anatomy$Zone), function(z) {
  het.res <- MPH.pct(hybrids, vars.names.root, anatomy[anatomy$Zone == z, ])
  cbind(het.res, Zone = z)
}))

#======================================================
# 3. Adjust p-values and classify heterosis type
#======================================================

heterosis.results.anatomy <- het.res.anatomy %>%
  group_by(Zone, Hybrid) %>%
  mutate(MPH.padj = p.adjust(MPH.p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    HetType = case_when(
      MPH.padj < 0.05 ~ "MPH",
      TRUE ~ "none"),
    HetType = factor(HetType, levels = c("none", "MPH")),
    comb = paste(Hybrid, Trait, sep = "_")
  )

# Save results
write_xlsx(heterosis.results.anatomy,
           "./results/phenotype/CK/MPH/anatomy_traits_MPH_pct_results.xlsx")

#=================================================
# 4. Plot significant MPH results
#=================================================

fig4a <- ggplot(
  data = heterosis.results.anatomy %>% filter(HetType == "MPH"),
  aes(x = Zone, y = MPH.pct, color = Hybrid, shape = Trait)) +
  geom_point(size = 3, alpha = 0.8, stroke = 0) +
  geom_line(aes(group = comb)) +
  geom_text(aes(label = Hybrid), color = "black", size = 1, hjust = -1) +
  mytheme

# Save plot
ggsave(fig4a,
       file = "./results/plots/main fig/fig4/anatomical_traits_MPH_plot.pdf",
       width = 12, height = 12, units = "cm")



################################################################################
##  Correlation between rhizosheath and other phenotypes
################################################################################

#===============================
# 1.  Load phenotype data
#===============================

pheno <- read.table("./results/phenotype/all_phenotpe_traits_mean.txt",
  header = TRUE)

RZsheath <- read_xlsx("./heterosis/data/heterosis root triats clean.xlsx", sheet = 5) %>%
  pivot_longer(
    cols = colnames(.)[2:23],
    names_to = "Genotype",
    values_to = "Rhizosheath")

# Calculate mean rhizosheath per genotype
rhizo.mean <- aggregate(Rhizosheath ~ Genotype, data = RZsheath, mean)

# Match order with pheno$Genotype
rhizo.mean <- rhizo.mean[match(pheno$Genotype, rhizo.mean$Genotype), ]

# Add rhizosheath to phenotype table
pheno <- cbind.data.frame(pheno, Rhizosheath = rhizo.mean$Rhizosheath)

#============================================================
# 2.  Correlation analysis (Spearman, FDR-adjusted)
#============================================================

traits.cor <- corr.test(
  pheno[, c(3:8, 10:12, 14:16, 18:21)], 
  method = "spearman",
  adjust = "fdr",
  ci = FALSE,
  normal = FALSE)

# Extract correlation and adjusted p-values
cor.mat   <- as.data.frame(traits.cor$r)
cor.padj  <- as.data.frame(traits.cor$p.adj)


#===================================================================
# 3.  Helper function: Flatten p-values into symmetric matrix
#===================================================================

flatten_to_matrix <- function(vec, n) {
  mat <- matrix(NA, n, n)
  mat[lower.tri(mat)] <- vec
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  diag(mat) <- 0
  return(mat)
}

# Build p-value matrix with correct dimensions
p_mat <- flatten_to_matrix(cor.padj$`traits.cor$p.adj`, ncol(cor.mat))
colnames(p_mat) <- colnames(cor.mat)
rownames(p_mat) <- rownames(cor.mat)

#==================================
# 4.  Visualization
#==================================

# Corrplot with significance masking
pdf(file = "./heterosis/results/plots/main fig/fig4/rhizosheath_corr_other_traits.pdf", 
  width = 6, height = 6)

corrplot(as.matrix(cor.mat), 
  type = "upper", 
  tl.col = "black",
  p.mat = as.matrix(p_mat),  # significance matrix
  sig.level = 0.05,          # cutoff
  insig = "blank", 
  diag = FALSE)
dev.off()


################################################################################
# Function for calculate mid/parent heterosis percent 
################################################################################

MPH.pct <- function(hb.lev, vars.names, het.input ){
  heterosis <- c()
  for (h in hb.lev){ 
    # Identify parents
    parts <- strsplit(h, "x")[[1]]
    mat <- parts[1]
    pat <- parts[2]
    for (v in vars.names) { 
      # Subset sample data to only include the 3 genotypes from this cross
      print(v)
      mat.data <- het.input[het.input$Genotype==mat, v] 
      pat.data <- het.input[het.input$Genotype==pat, v]
      hyb.data <- het.input[het.input$Genotype==h, v]
      
      # Calculate mean for hybrid, maternal, and paternal
      hyb.Mean <- mean(hyb.data, na.rm=TRUE) 
      mat.Mean <- mean(mat.data, na.rm=TRUE) 
      pat.Mean <- mean(pat.data, na.rm=TRUE)
      
      # Calculate mid-parent value
      midparent <- sum(mat.Mean, pat.Mean)/2 

      if(hyb.Mean > midparent){
        tt.MPH <- t.test(hyb.data, mu = midparent, alternative = 'greater')
      }else{
        tt.MPH <- t.test(hyb.data, mu = midparent, alternative = 'less')
      }
       
      
      # Store results
      tt.results <- data.frame('Maternal' = mat,'Paternal' = pat,'Hybrid' = h,
                              'Trait' = v, 'MPH.p' = tt.MPH$p.value,
                              'MPH.pct' = (hyb.Mean-midparent)/midparent)
      heterosis <- rbind(heterosis, tt.results)
      
      
      
    }
    
  }
  return(heterosis)
}



