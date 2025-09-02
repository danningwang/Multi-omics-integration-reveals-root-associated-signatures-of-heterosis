################################################################################
# Phenotype stats & plots + MPH calculation for heterosis
# Author: Danning Wang
# Project: Multi-omics data analysis realated to heterosis 
################################################################################


library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(rcompanion)   # cldList
library(FSA)          # dunnTest
library(RColorBrewer) # palettes
library(scales)       # alpha()


mytheme <- theme_bw() +
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),
    plot.title  = element_text(size = 10),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line   = element_line(colour = "black"),
    axis.text   = element_text(colour = "black", size = 8)
  )

setwd('H:/Research Projects/Heterosis/')

## ==================================================================
## 1. Prep data & comparison between hybrids and inbred lines
## ==================================================================

pheno22 <- read_xlsx('./data/phenotype/exploratory_exp_pheno_data.xlsx', sheet = 4) |>
  as.data.frame()

# Tidy factors used downstream
pheno22$Soil <- factor(pheno22$Soil, levels = c('Unsterilized', 'Sterilized'))
pheno22$Type <- factor(pheno22$Type, levels = c('Inbred', 'Hybrid'))
pheno22$Genotype <- gsub(' ', '', pheno22$Genotype)        
pheno22$comb <- paste(pheno22$Soil, pheno22$Type, sep = '_')

# Subset for PRL boxplot across Soil x Type, do this also for other phenotype
pheno.sub <- pheno22[, c('PRL', 'comb', 'Soil', 'Type')]

# Nonparametric group test + pairwise 
kruskal.test(PRL ~ comb, data = pheno.sub)

dunn.res <- dunnTest(PRL ~ comb, data = pheno.sub, method = "bh")
dunn.tab <- dunn.res$res

# Compact letter display per Soil x Type
label.dunn <- cldList(P.adj ~ Comparison, data = dunn.tab, threshold = 0.05)[, 1:2]
names(label.dunn) <- c("Group", "Letter")
label.dunn$Soil  <- sub('_.*', '', label.dunn$Group)
label.dunn$Type  <- sub('^[^_]*_', '', label.dunn$Group)
label.dunn$Soil  <- factor(label.dunn$Soil, levels = c('Unsterilized', 'Sterilized'))
label.dunn$Type  <- factor(label.dunn$Type, levels = c('Inbred', 'Hybrid'))

## ----- PRL boxplot with Dunn letters ----- ##
fig2a <- ggplot(data = pheno.sub, aes(x = Soil, y = PRL, fill = Type)) +
  geom_boxplot(outlier.colour = NA, width = 0.75) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
    alpha = 0.6, size = 2, shape = 1) +
  scale_fill_manual(values = c("#8498AB", "#EDC66A")) +
  geom_text(
    data = label.dunn,
    aes(x = Soil, y = rep(22, nrow(label.dunn)), label = Letter),
    position = position_dodge(width = 0.75),
    vjust = 0, size = 3) +
  labs(x = "Soil", y = "Primary root length") +
  mytheme

ggsave(fig2a, file = "./results/phenotype/Primary_root_length_boxplot.pdf",
       width = 12, height = 8, units = "cm")

## ===========================================================
## 2. MPH calculation function
## ===========================================================

# mid-parent heterosis percentage for a set of hybrids and traits
MPH.pct <- function(hb.geno, vars.names, het.input) {
  heterosis <- c()
  # do it for each hybrid-parents pairs
  for (h in hb.geno) {
    parts <- strsplit(h, 'x')[[1]]
    mat <- parts[1]  # maternal inbred
    pat <- parts[2]  # paternal inbred
    # do this for each phenotype
    for (v in vars.names) {
      mat.data <- het.input[het.input$Genotype == mat, v, drop = TRUE]
      pat.data <- het.input[het.input$Genotype == pat, v, drop = TRUE]
      hyb.data <- het.input[het.input$Genotype == h,   v, drop = TRUE]
      
      hyb.Mean <- mean(hyb.data, na.rm = TRUE)
      mat.Mean <- mean(mat.data, na.rm = TRUE)
      pat.Mean <- mean(pat.data, na.rm = TRUE)
      
      midparent <- (mat.Mean + pat.Mean) / 2
      
      t.res <- t.test(hyb.data, mu = midparent, alternative = 'two.sided')
      
      res <- data.frame('Maternal' = mat,'Paternal' = pat,'Hybrid' = h,
                       'Trait' = v, 'test.p' = t.res$p.value,
                       'MPH.pct'=(hyb.Mean-midparent)/midparent)
      heterosis = rbind(heterosis, res)

    }
  }
  return(heterosis)
}

## ==================================================
## 3. Compute MPH under each soil condition
## ==================================================

# Traits to evaluate 
vars.names <- colnames(pheno22)[c(2, 5, 7)]  # PRL, PRW, Standard_rhizosheath

# Define hybrids 
hybrids <- unique(pheno22$Genotype[grepl('x', pheno22$Genotype, fixed = FALSE)])

het.res.pct <- c()

for (s in unique(pheno22$Soil)) {
  pheno22_sub <- pheno22[pheno22$Soil == s, ]
  mph <- MPH.pct(hybrids, vars.names, pheno22_sub)
  het.res.pct = rbind(het.res.pct, cbind(mph, Soil = s))
}

het.res.pct$Soil <- factor(het.res.pct$Soil, levels = c('Unsterilized', 'Sterilized'))

## ----- tests on MPH by Soil for each trait -----
w.res.PRL  <- t.test(MPH.pct ~ Soil, data = subset(het.res.pct, Trait == 'PRL'))
w.res.PRW  <- t.test(MPH.pct ~ Soil, data = subset(het.res.pct, Trait == 'PRW'))
w.res.RHZ  <- t.test(MPH.pct ~ Soil, data = subset(het.res.pct, Trait == 'Standard_rhizosheath'))

## ----- MPH boxplot -----
fig2d <- ggplot(data = het.res.pct,
  aes(x = Trait, y = MPH.pct, fill = Soil)) +
  geom_boxplot(outlier.colour = NA, width = 0.7) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
    alpha = 0.6, size = 2, shape = 1) +
  scale_fill_manual(values = c("#D9DEE7", "#C9A1CA")) +
  ylim(NA, 1) +
  labs(x = "Trait", y = "MPH %") +
  mytheme

ggsave(fig2d, file = "./results/phenotype/MPH_boxplot_soil_22geno.pdf",
       width = 12, height = 8, units = "cm")

## ==========================================================
## 4. Figure 2e: prep phenotype data and boxplot of MPH
## ==========================================================

pheno.3tr <- read_xlsx('./data/phenotype/exploratory_exp_pheno_data_3_treatments.xlsx', sheet = 2)

# Long -> remove trailing 1..4 replicate numbers in column names, then widen back
pheno.3tr.long <- pheno.3tr |>
  pivot_longer(cols = 3:14, names_to = 'pheno', values_to = 'value') |>
  mutate(pheno = gsub('[1234]', '', pheno))

pheno.3tr.df <- pheno.3tr.long |>
  group_by(Treatment, Genotypes, pheno) |>
  mutate(replicate = row_number()) |>
  ungroup() |>
  pivot_wider(names_from = pheno, values_from = value) |>
  as.data.frame()

# Tidy labels
pheno.3tr.df$Genotypes <- gsub(' ', '', pheno.3tr.df$Genotypes)
pheno.3tr.df$Type <- ifelse(grepl('x', pheno.3tr.df$Genotypes), 'Hybrids', 'Inbredlines')
pheno.3tr.df$Type <- factor(pheno.3tr.df$Type, levels = c('Inbredlines', 'Hybrids'))
pheno.3tr.df$comb <- paste(pheno.3tr.df$Treatment, pheno.3tr.df$Type, sep = '_')
colnames(pheno.3tr.df)[2] <- 'Genotype'

# Traits to evaluate 
vars.names.3tr <- colnames(pheno.3tr.df)[4:6]
hybrids.3tr    <- unique(pheno.3tr.df$Genotype[grepl('x', pheno.3tr.df$Genotype)])

# MPH per treatment
het.res.pct.3tr <- c()

for (treat in unique(pheno.3tr.df$Treatment)) {
  pheno.sub <- subset(pheno.3tr.df, Treatment == treat)
  mph   <- MPH.pct(hybrids.3tr, vars.names.3tr, pheno.sub)
  het.res.pct.3tr = rbind(het.res.pct.3tr, cbind(mph, Treatment = treat))
}


# do this for each trait
dunn.res <- dunnTest(MPH.pct~Treatment, data = het.res.pct[het.res.pct$Trait == 'Normsheath', ], method = "bh")
print(dunn.res$res)

## ----- Boxplot across treatments ----- ##
fig2e <- ggplot(het.res.pct.3tr, aes(x = Trait, y = MPH.pct, fill = Treatment)) +
  geom_boxplot(outlier.colour = NA, width = 0.75) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
    alpha = 0.6, size = 2, shape = 1
  ) +
  scale_fill_manual(values = c("grey", "#9FDAF7", "#8FB4DC")) +
  labs(x = "Trait", y = "MPH") +
  mytheme

ggsave(fig2e, file = "./results/phenotype/pheno_MPH_boxplot_3_treatments.pdf",
       width = 18, height = 10, units = "cm")




