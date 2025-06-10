# Differential Gene Expression analysis using 'pasilla' dataset

# Source: Drosophila melanogaster (fruit fly)

# Design: Comparison of cells treated with pasilla RNAi (to knock down a splicing regulator) vs untreated

# Goal: Identify genes whose expression changes due to RNAi knockdown


# installing and loading packages

if (!requireNamespace("BiocManager",quietly = T))
  install.packages("BiocManager")

BiocManager::install(c("pasilla","SummarizedExperiment", force= T))
BiocManager::install(c("AnnotationDbi","org.Dm.eg.db"))
library(DESeq2)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(tidyverse)
library(AnnotationDbi)
library(org.Dm.eg.db)

# cheching data found in pasilla
library(pasilla)
data(package =  "pasilla")

# LOADING DATA
# The dataset is included inside the package’s extdata directory

counts_file   <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
anno_file     <- system.file("extdata/pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)

counts_data <- read_tsv(counts_file, col_types = cols()) %>% 
  column_to_rownames("gene_id") 
sample_info <- read_csv(anno_file, col_types = cols())

# Clean and Prepare Sample Metadata to match conts_data

sample_info <- sample_info %>% 
  filter(type == 'single-read') %>% 
  mutate(condition = factor(condition), 
         file = gsub("fb$","",file))

counts_data <- counts_data[, sample_info$file]

# checking consistent order
sample_info <- as.data.frame(sample_info)
rownames(sample_info) <- sample_info$file
all(colnames(counts_data) %in% rownames(sample_info))

# creating DESeq data set and running DESeq2

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = sample_info,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
view(res)

# gene annotation
library(org.Dm.eg.db)
res$symbol <- mapIds(org.Dm.eg.db,
                         keys = rownames(res),
                         column = "SYMBOL",
                         keytype = "FLYBASE",
                         multiVals = "first")
view(res)  
write.table(res, file = "pasilla_DESeqDataSet.csv")
# visualization
# 1. MA plot
plotMA(res, ylim = c(-5,5))

summary(res)

# volcano plot

 EnhancedVolcano(res,
                 lab = res$symbol,
                 x = "log2FoldChange" ,
                 y = "pvalue",
                legendLabels= c("p>0.1 & fc<1(not sig.)",
                   "p>0.1","p<0.1 & fc<1","p<0.1 & fc>1 (sig.)")

# Genes like Kal1, Sox100B, CG34330 are highly upregulated with strong statistical support
# Genes like Ant2, Ama, Prip, Hsp23 appear on the far left, showing reduced expression post-siRNA.
# The bulk of genes are gray or green and lie near the center → indicating no meaningful change, which is typical in RNA-seq.








