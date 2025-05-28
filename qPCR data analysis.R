library(tidyverse)

#--------------------TUTORIAL PROJECT----------------------------
# U2OS cells (a human osteosarcoma cell line).
# Six samples (5 treatment + 1 control_siscr) treated with 5 different siRNAs
# House Keeping gene: GAPDH
# Target gene OGG1 gene. DNA repair gene
# 3 different primers used in PCR aplification OGG1_1, OGG1_2, OGG1_3
# RT-qPCR to measure OGG1 gene expression in U2OS cells after treatment with different siRNAs.


res = read.delim("qPCR.txt","\t",
                  skip = 10)
view(res)
# Step.1: selecting desired columns into a variable
# data cleaning
impres<-res %>% 
  select(sample='Sample Name',primer='Detector Name',ct='Ct') %>% 
  drop_na() %>% 
  filter(primer!='Detector Name',ct!='Ct') %>%
  filter(ct!='Undetermined') %>%
  mutate_at('ct',as.numeric)
  # changing dtype of ct to umeric
  # mutate_at used to do operation on a column
view(impres)

# Step.2: Calculating housekeeping gene avd_ct per sample

HKgene<-impres %>%
  filter(primer=='GAPDH') %>% 
  group_by(sample) %>% 
  summarize(HK_ct=mean(ct)) %>% 
  view()


# Step.3 Merging and calculating ΔCt

results<-impres %>% 
  filter(primer!='GAPDH') %>% 
  left_join(HKgene,by='sample') %>%
  mutate(dct= ct - HK_ct) %>% 
  view()

# Step.4: geting control ΔCt  (for ΔΔCt  baseline)

control_dct<-results %>% 
  filter(sample=='U2OS_siscr') %>% 
  group_by(primer) %>% 
  summarise(ctrl_dct = mean(dct)) %>% 
  view()

# Step 5: Add ΔΔCt and fold change to the same tibble

results<-results %>% 
  left_join(control_dct,by='primer') %>% 
  mutate(ddct = dct - ctrl_dct,
         fold_change = 2^(-ddct))
view(results)  

# Step 6: Summary table of mean and SD per group

summary_results<-results %>% 
  group_by(sample,primer) %>% 
  summarise(avg_dct = mean(dct),
         sd_dct = sd(dct),
         avg_ddct = mean(ddct),
         sd_ddct = sd(ddct),
         sd_fc = sd(fold_change),
         avg_fc = mean(fold_change),.groups = 'drop')
view(summary_results)

# Visualization
# ΔCt Plot

ggplot(summary_results,aes(x=sample,y=avg_dct,fill = primer))+
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin =avg_dct-sdev ,
                    ymax =avg_dct+sdev),
                width=.2,position = position_dodge(.9))

# ΔΔCt Plot for all samples hilighting 3 primers used

ggplot(summary_results,aes(x=sample,y=avg_ddct,fill = primer))+
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin = avg_ddct-sd_ddct,ymax = avg_ddct+sd_ddct),
                width=.2,
                position = position_dodge(.9))+
  ggtitle('ΔΔCt Plot')
  
# Fold-Change plot

ggplot(summary_results,aes(x=sample,y=avg_fc,fill = primer))+
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin = avg_fc-sd_fc, ymax = avg_fc+sd_fc),
                width=.2,position = position_dodge(.9))+
  scale_y_log10()+ggtitle('GE Fold Change')
  


