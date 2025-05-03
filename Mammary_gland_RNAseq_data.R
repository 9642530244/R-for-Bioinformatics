library(tidyverse)
sampl <- read.csv("GSE60450_filtered_metadata-1.csv")
cnt <- read.csv("GSE60450_GeneLevel_NormalizedCPM.and_.TMM_data.csv")

#1) reading files into vars
sampl = read.csv("filename.cas") # file is in cwd 
cnt = read.csv(filename.csv)

#2) observe the tables, assign col names if not theme_replace
colnames(sampl)[1] <- "sampl_id"
colnames(cnt)[1] <- "gene_id"

#3) formatting wide table into long table # ideal for plotting with ggplot
seqdata = pivot_longer(cnt, cols = starts_with("GSM"), 
                       names_to = "sample",values_to = "count")

#4) joining the two tables having same column with full join 
allinfo = full_join(seqdata,sampl, by = c("sample" = "sampl_id"))


view(head(allinfo))
# labelling the 6 charecters to simple names for plotting
allinfo <- mutate(allinfo,group = case_when(
  str_detect(characteristics, "basal.*virgin") ~ "bvirg",
  str_detect(characteristics, "basal.*preg") ~ "bpreg",
  str_detect(characteristics, "basal.*lact") ~ "blact",
  str_detect(characteristics,"luminal.*virg") ~ "lvirg",
  str_detect(characteristics,"luminal.*preg") ~ "lpreg",
  str_detect(characteristics,"luminal.*lact") ~ "llact"
))
view(head(allinfo,n=3))
mygenes2 <- c("Csn1s2a", "Csn1s1", "Csn2", "Glycam1", "COX1", "Trf", "Wap", "Eef1a1")
mygen2_cnt <- filter(allinfo, gene_symbol %in% mygenes2)
mygen2_cnt
ggplot(data =mygen2_cnt,mapping = aes(x=group,y=log2(count+1),colour=group))+
  geom_boxplot()+facet_wrap(~gene_symbol)
#point plot
ggplot(data=mygen2_cnt,mapping=aes(x=group,y=log2(count+1),colour=group))+
  geom_point()+facet_wrap(~gene_symbol)

#The luminal compartment is transcriptionally active during pregnancy and lactation, aligning with its role in milk synthesis.

#The basal compartment shows little expression of milk protein genes, supporting its structural, non-secretory role.

#Gene expression data reflects developmental programming of the mammary gland in response to hormonal cues.

#This pattern can be used to identify differentially expressed genes (DEGs) and infer functional specialization in mammary epithelial subtypes.

#Would you like a short report version of this summary for your project or
























# another method i tried
mygenes <- allinfo %>% 
  group_by(gene_symbol) %>% 
  summarise(totg_count = sum(count)) %>% 
  arrange(desc(totg_count)) 
head(mygenes,n=8)
ggplot(data = head(mygenes,n=8),
       mapping = aes(x=gene_symbol,y=log2(totg_count+1),fill =gene_symbol ))+
  geom_boxplot()+ facet_wrap(~ gene_symbol)
# end






  
#
 