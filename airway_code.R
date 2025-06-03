# Bioconductor 'airway' data set
# Gene Expression analysis with DESeq2

BiocManager::install("DESeq2")
BiocManager::install("airway")

library(airway)
library(DESeq2)
library(tidyverse)

# loading dta

data(airway)

# isolating necessary data colData and assay

sample_info<-as.data.frame(colData(airway))

sample_info %>% 
  select(cell,dex) %>% 
  setNames(c("cellLine","dexamethasone")) %>% 
  mutate(dexamethasone = fct_recode(dexamethasone, 
                                    "treated" = "trt", "untreated" = "untrt")) %>% 
  write.table(file = "sample_info.csv", sep="," , col.names = T, row.names = T,quote = F)

# preparation of counts dataset from airway (assay)

countsData<-assay(airway) %>% 
  write.table(file ="count_data.csv",sep =",",col.names =T,row.names =T,quote = F)

# calling in sample_info and counts data

counts_data<- read.csv("count_data.csv")
colData<-read.csv("sample_info.csv")

# checking row names and order in sampleinfo matching with col names in counts data

all(colnames(counts_data)%in% rownames(colData)) # present 
all(colnames(counts_data) == rownames(colData)) # order

# creating a DESeq object and running DESeq2

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~dexamethasone)
str(dds)

# filter out genes with less read counts into an object

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# setting levels treated and untreated. reference would be untreated
# checking the DE in treated w.r.t untreated

dds$dexamethasone<-relevel(dds$dexamethasone, ref = "untreated")

dds <- DESeq(dds)
# we get a table which can be seen with results function
# table contains, log2fold change, p_vals, padjusted vals

res <- results(dds)

summary(res)
# in summary the p value = o.1
# changing the p value to 0.01
# p=0.01 means that the DE of genes by random chance is less than 1% 

res2 <- results(dds,alpha = 0.01)
summary(res2)

# visualization
# we need a data base for annotation to ensemble ids

BiocManager::install('org.Hs.eg.db')

library(AnnotationDbi)
library(org.Hs.eg.db)

# saving the results as csv dataset for future use

write.table(res2,file = "airway.csv",quote = T, sep = ",")

# loading the results airway file backin

rnaseq <- read.csv('airway.csv')
View(rnaseq)

library(dplyr)
rnaseq <- rnaseq %>% 
  mutate(symbol = mapIds(org.Hs.eg.db, key = rownames(rnaseq),
                         column = ('SYMBOL'), keytype = "ENSEMBL")) %>% 
  data.frame() %>% 
  View()
# gene names are annotated w.r.t emsembl ids (rows)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(rnaseq,
                lab = rnaseq$symbol,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 0.01,
                legendLabels = c("p>0.01 & fc<1(not sig.)",
                                 "p>0.01","p<0.01 & fc<1","p<0.01 & fc>1 (sig.)"))

# The volcano plot displays the differential expression of 22,369 genes. 
# Genes labeled at the top-left and top-right represent the most significantly downregulated and upregulated genes, respectively, under treatment. 
# For instance, FKBP5 and DUSP1 are highly upregulated, while COL1A1 and ZBTB16 are strongly downregulated, 
# making them potential candidates for downstream biological validation.

















