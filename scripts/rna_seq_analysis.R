################################################################################
## Root RNA-Seq Data Analysis 
## Author: Danning Wang
## Project: Multi-omics data analysis of heterosis 
################################################################################



library(DESeq2)
library(Rsubread)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(readr)
library(glue)


#----------------------------#
#   Utility Functions
#----------------------------#

## 1. Read featureCounts output
read_counts <- function(file) {
  counts <- read.delim(file, header = TRUE, row.names = 1)
  counts <- counts[, 6:ncol(counts)]   # keep count columns only
  counts
}

## 2. Filter low-expression genes
filter_counts <- function(counts, min_cpm = 10, min_samples = 3) {
  keep <- rowSums(counts >= min_cpm) >= min_samples
  counts[keep, ]
}

## 3. Build DESeq2 object for one treatment
make_dds <- function(counts, metadata, treatment, design = ~Hybrid_inbred) {
  md <- metadata %>% filter(Treatment == treatment)
  counts <- counts[, md$SampleID]
  stopifnot(all(colnames(counts) == md$SampleID))
  DESeqDataSetFromMatrix(countData = counts,
                         colData = md,
                         design = design)
}

## 4. Run DESeq2 and extract results
run_deseq <- function(dds, contrast, alpha = 0.05, lfc = 1) {
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast, alpha = alpha, lfcThreshold = lfc)
  res[order(res$padj), ]
}

## 5. Volcano plot
plot_volcano <- function(res, title = "Volcano Plot") {
  resdf <- as.data.frame(res) %>% 
    mutate(GeneID = rownames(res),
           status = case_when(
             log2FoldChange >= 1 & padj < 0.05 ~ "Up",
             log2FoldChange <= -1 & padj < 0.05 ~ "Down",
             TRUE ~ "NS"
           ))
  
  ggplot(resdf, aes(x = log2FoldChange, y = -log10(padj), fill = status)) +
    geom_point(shape = 21, colour = "black", size = 2) +
    scale_fill_manual(values = c("Up" = "#ffad73", "Down" = "#26b3ff", "NS" = "grey")) +
    theme_classic() +
    labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
    theme(legend.position = "top")
}

## 6. PCA plot
plot_pca <- function(dds, ntop = 500, intgroup = "Hybrid_inbred", title = "PCA") {
  vsd <- vst(dds, blind = FALSE)
  pca <- prcomp(t(assay(vsd)[order(rowVars(assay(vsd)), decreasing = TRUE)[1:ntop], ]))
  percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pcaData <- data.frame(pca$x, colData(dds))
  
  ggplot(pcaData, aes(PC1, PC2, color = intgroup, label = SampleID)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3) +
    labs(x = glue("PC1 ({percentVar[1]}%)"),
         y = glue("PC2 ({percentVar[2]}%)"),
         title = title) +
    theme_classic()
}

## 7. Heatmap of top DEGs
plot_heatmap <- function(dds, res, topN = 50, ann_col = NULL) {
  vsd <- vst(dds, blind = FALSE)
  topGenes <- head(order(res$padj), topN)
  mat <- assay(vsd)[topGenes, ]
  pheatmap(mat, annotation_col = ann_col)
}

################################################################################
##  Process bam files and basic downstream analysis
################################################################################

#-----------------------------#
# 1. Define file paths
#-----------------------------#
bam_dir   <- "/home/data/Precision-5820_4TB/heterosis/bam"
gtf_file  <- "~/emmynoether/data/zm_ref_v5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.53.chr.gtf"
out_dir   <- "~/Heterosis/rna_seqs"

# Ensure output folders exist
dir.create(file.path(out_dir, "intermediate_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

#-----------------------------#
# 2. Collect BAM files
#-----------------------------#
all_bam_files <- list.files(path = bam_dir, pattern = "\\.bam$", full.names = TRUE)
stopifnot(length(all_bam_files) > 0)  # sanity check


#-----------------------------#
# 3. Run featureCounts
#-----------------------------#
fc <- featureCounts(
  files = all_bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  countMultiMappingReads = FALSE,
  useMetaFeatures = TRUE,
  countChimericFragments = FALSE,
  nthreads = 18
)

#-----------------------------#
# 4. Save raw counts
#-----------------------------#
raw_counts <- fc$counts
write.csv(raw_counts,
          file = file.path(out_dir, "intermediate_data/raw_gene_counts.csv"),
          row.names = TRUE, quote = FALSE)

#--------------------------------------------------#
# 5. PCA and identify DEGs for each soil condition
#--------------------------------------------------#

## Read metadata
metadata_file <- "intermediate_data/rna_sample_metadata_correct.txt"
metadata <- read.table(metadata_file, header = TRUE)

# Filter genes
counts_filt <- filter_counts(raw_counts)

# Run for each treatment
for (trt in c("CK", "LP", "LN")) {
  message(glue("=== Running DESeq2 for {trt} ==="))
  
  # Make DESeq object
  dds <- make_dds(counts_filt, metadata, treatment = trt)
  
  # Run DESeq2
  res <- run_deseq(dds, contrast = c("Hybrid_inbred", "Hybrid", "Inbred"))
  
  # Save results
  write.csv(as.data.frame(res),
            glue("results/deseq2/{trt}_DESeq2_results.csv"),
            quote = FALSE)
  
  # Volcano plot
  ggsave(glue("results/plots/{trt}_volcano.pdf"),
         plot_volcano(res, title = glue("{trt}: Hybrid vs Inbred")),
         width = 6, height = 5)
  
  # PCA
  ggsave(glue("results/plots/{trt}_PCA.pdf"),
         plot_pca(dds, title = glue("{trt}: PCA")),
         width = 6, height = 5)
  
}










