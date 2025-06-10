

#1. loading libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)

#2. reading data
sample_info <- read.csv("design.tsv",sep = "\t",row.names = 1)
count_table <- read.csv("raw_counts.tsv",sep = '\t',row.names=1)
view(count_table)

dim(count_table)
dim(sample_info)

#3. setting factor levels control vs. tgf-beta

sample_info$Group <- factor(sample_info$Group, levels = c("control","tgf-beta"))

sample_info$Group

# creating DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = sample_info,
                              design = ~Group)
# filtering low count genes
# keep genes with atleast N countrs >=10, where N = size of smallest group

keep <- rowSums(counts(dds) >=10 ) >= min(table(sample_info$Group))

dds <- dds[keep,]

# performing statistical tests
dds <- DESeq(dds,test = "Wald",sfType = "poscount")
dsq_res <- results(dds)
view(dsq_res)
dsq_res <- as.data.frame(dsq_res)

# here in results gene names are set as index, changing into a column
dsq_res$gene_name <- row.names(dsq_res)
view(dsq_res)

# subsetting order of cols

dsq_res <- dsq_res %>% 
  select("gene_name","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean")
view(dsq_res)

# saving results
write.table(dsq_res,file = "dsq_res_TGF.tsv",row.names = F,sep = "\t")

# filter gens with padj < 0.05 and LFC <=-1 or >=1
# abs() consider both -1 and 1
deg <- subset(dsq_res,padj<0.05 & abs(log2FoldChange)>=1)
dim(deg)
dim(dsq_res)

# ordering deg by padj
deg <- deg[order(deg$padj),]


# visualization
# dispersion plot with title
plotDispEsts(dds,main = "GSE203159 Dispersion Estimates")

# histogram
hist(dsq_res$padj, breaks = seq(0,1,length=21))

# PCA plot : dimentions reduction plot
vsd <- vst(dds,blind = F)
colnames(vsd)
plotPCA(vsd, intgroup = "Group")
 # 99% of variance is due to presence of TGF-beta gene in test than control

# heat map of sample to sample distance/similarity

samp_dist <- dist(t(assay(vsd)))
samp_dist_matrx <- as.matrix(samp_dist)
pheatmap(samp_dist_matrx)

# there is a clear clustering between control and TGF-beta samples

# MAplot
plotMA(dds, ylim = c(-3,3))

# blue points are with p<0.05

# volcano plot
library(EnhancedVolcano)
EnhancedVolcano(dsq_res,
                lab = row.names(dsq_res),
                x = "log2FoldChange",
                y = "pvalue",
                legendLabels = c("NO sig.","FC>1","p<0.05","p<0.05 & FC>1 Sig."))





