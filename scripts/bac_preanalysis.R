################################################################################
## Bacteria 16S rRNA Analysis 
## Author: Danning Wang
## Project: Multi-omics data analysis of heterosis 
################################################################################


library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(DESeq2)
library(vegan)
library(phyloseq)
library(readr)
library(dplyr)
library(reshape2)
library(qiime2R)
library(biomformat)
library(lmerTest)
library(rcompanion)
library(FSA)
library(viridis)


mytheme <- theme_bw() +
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),
    plot.title = element_text(size = 10),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 8)
  )

#===========================================================
# 1. Import data
#===========================================================
setwd("H:/Research Projects/Heterosis/")

## ASV table from output of QIIME
biom.table <- read_biom("./data/feature-table.biom")
raw.bac.table <- as.matrix(biom_data(biom.table))
dim(raw.bac.table)    # check dimensions

## FeatureID -> ASV mapping
featureID2asv <- data.frame(
  FeatureID = rownames(raw.bac.table),
  ASV = paste0("bASV", seq_len(nrow(raw.bac.table))))
write.table(featureID2asv,
            "./intermediate_data/feature_ID_to_ASV_file.txt",
            row.names = FALSE, sep = "\t", col.names = TRUE)

asv.table <- raw.bac.table
rownames(asv.table) <- featureID2asv$ASV

## Taxonomy from output of QIIME
bac.taxa <- read.table("./data/taxonomy.tsv", sep = "\t", header = TRUE)
stopifnot(all(rownames(raw.bac.table) == bac.taxa$Feature.ID))

bac.taxa.tab <- bac.taxa %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(
    Kingdom = gsub("d__", "", Kingdom),
    Phylum  = gsub(" p__", "", Phylum),
    Class   = gsub(" c__", "", Class),
    Order   = gsub(" o__", "", Order),
    Family  = gsub(" f__", "", Family),
    Genus   = gsub(" g__", "", Genus),
    Species = gsub(" s__", "", Species)
  )

rownames(bac.taxa.tab) <- featureID2asv$ASV
bac.taxa.mat <- as.matrix(bac.taxa.tab[, -1])
bac.taxa.mat[is.na(bac.taxa.mat)] <- "Unidentified"

## Metadata
bac.metadata <- read.table(
  "./intermediate_data/metadata.txt",
  sep = "\t", header = TRUE)
stopifnot(all(bac.metadata$sample.id == colnames(raw.bac.table)))
rownames(bac.metadata) <- bac.metadata$sample.id

#============================================================
# 2. Build phyloseq object
#============================================================

bac.OTU <- otu_table(asv.table, taxa_are_rows = TRUE)
bac.TAX <- tax_table(bac.taxa.mat)
bac.samples <- sample_data(bac.metadata)

bac.ps <- phyloseq(bac.OTU, bac.TAX, bac.samples)
saveRDS(bac.ps, './intermediate_data/bac_ps_raw.RDS')

#==============================================================
# 3. Filtering
#==============================================================
ps.prefilter <- subset_taxa(bac.ps,
                            Kingdom == "Bacteria" &
                              Phylum != "Unidentified" &
                              Family != "Mitochondria" &
                              Order != "Chloroplast")

## Relative abundance filtering: keep ASVs >0.05% in ≥5% samples
ps.RA <- transform_sample_counts(ps.prefilter, function(x) x / sum(x))
ps.keep <- filter_taxa(ps.RA, function(x) sum(x > 0.0005) >= 0.05 * nsamples(ps.RA), TRUE)

ps.abund <- prune_taxa(taxa_names(ps.keep), ps.prefilter)

## remove samples with <2000 reads
ps.abund <- subset_samples(ps.abund, sample_sums(ps.abund) > 2000)
saveRDS(ps.abund, './intermediate_data/ps_abund_filtered.RDS')

#=============================================================
# 4. Alpha diversity
#=============================================================

## Rarefaction
ps.rare <- rarefy_even_depth(ps.abund, 15000, rngseed = 123555, replace = FALSE)

## Diversity measures
adiv <- estimate_richness(ps.rare, measures = c("Observed", "Shannon"))
adiv <- cbind(adiv, sample_data(ps.rare))
write.table(adiv,
            "./intermediate_data/alpha_diversity_table_15000.txt",
            sep = "\t", quote = FALSE)

## Kruskal-wallis test and boxplot
adiv$CompTre <- paste(adiv$treatment, adiv$hybrid, sep = "_")
k.res <- kruskal.test(Shannon ~ CompTre, data = adiv)
dunn.res <- dunnTest(Shannon ~ CompTre, data = adiv, method = "bh")

label.dunn <- cldList(P.adj ~ Comparison, data = dunn.res$res, threshold = 0.05)[, 1:2]
label.dunn$treatment <- gsub("_.*$", "", label.dunn$Group)
label.dunn$hybrid <- sub("^[^_]*_", "", label.dunn$Group)

alpha.plot <- ggplot(adiv, aes(x = treatment, y = Shannon, fill = hybrid)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(color = "black",
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
    alpha = 0.4, size = 2) +
  scale_fill_brewer(palette = "Dark2") +
  ylim(6, 6.3) +
  mytheme +
  geom_text(data = label.dunn,
            aes(y = 6.3, label = Letter),
            position = position_dodge(width = 1))

ggsave(alpha.plot,
       file = "./results/bacteria/pre_process/alpha_div/Shannon_div_rhizo.pdf",
       width = 12, height = 8, units = "cm")

#=============================================================
# 5. Beta diversity
#=============================================================

## Rhizosphere under soil condition CK subset
ps.rz <- subset_samples(ps.abund, compartment == "R" & treatment == "CK")

## DESeq2 VST normalization
bac.deseq <- phyloseq_to_deseq2(ps.rz, ~genotype)
bac.deseq <- estimateSizeFactors(bac.deseq, type = "poscount")
bac.vst <- varianceStabilizingTransformation(bac.deseq, blind = FALSE)

bac.vst.ck = ps.rz
otu_table(bac.vst.ck) = otu_table(assay(bac.vst), taxa_are_rows = T)
bac.vst.ck = transformSampleCounts(bac.vst.ck,function(x) ifelse(x<0,0,x)) 

# PCoA (Bray-Curtis)
bray.dist <- phyloseq::distance(bac.vst.ck, method = "bray")
ord <- ordinate(bac.vst.ck, method = "PCoA", distance = bray.dist)

pcoa.df <- cbind(ord$vectors[, 1:2], sample_data(bac.vst.ck))
pct <- ord$values$Relative_eig[1:2]
percentage <- paste0(c("PCoA1", "PCoA2"), " (", round(pct * 100, 2), "%)")

ggplot(pcoa.df, aes(Axis.1, Axis.2, color = genotype, shape = hybrid)) +
  geom_point(size = 3) +
  geom_text(aes(label = genotype), nudge_y = -0.00025, size = 2) +
  mytheme +
  xlab(percentage[1]) + ylab(percentage[2]) +
  coord_fixed()

# PERMANOVA test
adonis2(bray.dist ~ treatment + genotype,
        data = as(sample_data(bac.vst.ck), "data.frame"),
        permutations = 1999, by = "margin")



