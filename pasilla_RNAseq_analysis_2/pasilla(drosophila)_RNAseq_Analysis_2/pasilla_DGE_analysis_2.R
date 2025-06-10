
# Pasilla gene (the Drosophila homologue of the mammalian splicing regulators Nova-1 and Nova-2 proteins) 
# the Pasilla (PS) gene in Drosophila melanogaster is depleted by RNA interference (RNAi)
# Total RNA was then isolated and used to prepare both single-end and paired-end reads
# 4 untreated samples: GSM461176, GSM461177, GSM461178, GSM461182
# 3 treated samples (Pasilla gene depleted by RNAi): GSM461179, GSM461180, GSM461181 
# two treated and two untreated samples are paired-end reads, remaining samples single-end reads.
# Comparing RNA-Seq data for the treated and the untreated samples to find effects of PS gene depletion on GE.

# loading libraries
library(DESeq2)
BiocManager::install("pheatmap")
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)

# loading the datasets

count_data <- read.csv("count_matrix.csv",header = T, row.names = 1)

sample_info <- read.csv("design.csv",header = T,row.names = 1)

colnames(sample_info)

# setting factor levels (treatment, sequencing type)

sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$Sequencing <- factor(sample_info$Sequencing)

# creating a DESeq object and importing  count data and samp info

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ Sequencing+Treatment)

# setting reference for treatment factor
# first one in c() becomes reference
# "untreated" – becomes the reference level (also called the baseline)
# "treated" – this is the level that will be compared against the reference

dds$Treatment <- factor(dds$Treatment, 
                        levels = c("untreated","treated"))
# filtering genes with read counts > 10
# Keep genes that are expressed (count > 10) in at least 3 samples (because the smallest group has 3 samples).

keep <- rowSums(counts(dds) > 10) >= min(table(sample_info$Treatment))

dds <- dds[keep,]

# performing statistical tests to find DEG
dds <- DESeq(dds)
dseq_res <- results(dds)
view(dseq_res)

# changing dseq results into R dataframe
dseq_res <- as.data.frame(dseq_res)

# ordering results by increasing p value
res_ordered <- dseq_res[order(dseq_res$pvalue),]

# making some queries
dseq_res["FBgn0037363",] # down regulated & diff. expressed

# filtering genes by p <0.05  logFC (<1 & >1)
filtered <- dseq_res %>% 
  filter(dseq_res$padj < 0.05)
filtered <- filtered %>% 
  filter(abs(filtered$log2FoldChange)>1)
 dim(dseq_res)
 dim(filtered)
 
# saving results dataset
write.csv(dseq_res,"dds_results2.csv")
write.csv(filtered,"filtered_dds_results2.csv")

# saving the normalized counts
# raw read counts not adjusted for differences in sequencing depth or sample composition.
# normalization accounts for Library size differences (samples with more reads overall) and Sample-to-sample variability
norm_counts <- counts(dds, normalized = T)
write.csv(norm_counts,"normalizd_counts2.csv")

# Visualization
# dispersion plot

plotDispEsts(dds)

#How well does DESeq2 model the variability in my gene expression data?
#Are the black dots mostly near the red line? ✅ That means ,data fits the model well.
#Are the blue dots not wildly far from the others? ✅ That means DESeq2 did a good job estimating variability.

# PCA plot: dimentionality reduction technique
# used to explain variance in gene expression 

vsd <- vst(dds, blind = F)
plotPCA(vsd,intgroup = c("Sequencing","Treatment"))

#Each point = a sample (here 7 samples)
#Samples with similar gene activity land near each other.
#It helps you see which factor (treatment, batch, etc.) is driving the biggest differences in your data.
#PC1 is the first principal component — treated vs untreated
#it's the direction that captures the most variation across all your samples.
# PC2 : Single-end vs paied-end

#Heatmap
# heatmap of sample-to-sample distance matrix
# distance matrix
sample_dist <- dist(t(assay(vsd)))
sample_dis_matrix <- as.matrix(sample_dist)
colnames(sample_dis_matrix)
# generating heatmap
pheatmap(sample_dis_matrix, 
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist)
# dark blue: less distance between samples- similar expression
# red: large distance between samples- distinct expression (different groups)

# MA plot
plotMA(dds,ylim = c(-2,2))

# blue points have adj_p < 0.1
# to remove noice
# resLFC <- lfcShrink(dds, coef="Treatment_treated_vs_untreated",type="apeglm)
# pllotMA(resLFC,ylim=c(-2,2))

# annotating the genes 
library(org.Dm.eg.db)
dseq_res$symbol <- mapIds(org.Dm.eg.db,
                          keys = rownames(dseq_res),
                          column = "SYMBOL",
                          keytype = 'FLYBASE',
                          multivals = 'first')
view(dseq_res)
summary(dseq_res)
# volcano plot
library(EnhancedVolcano)
EnhancedVolcano(dseq_res,
                lab = dseq_res$symbol,
                x = "log2FoldChange",
                y = "pvalue",
                legendLabels = c("p>0.1 & FC<1","FC>1","p<0.1","p<0.1 & FC>1 Significant"))

# Ant2 and Ama genes are statistically sig. and are upregulated
# Kal1 & SPARC are down regulated
































