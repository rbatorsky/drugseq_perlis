LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
select <- dplyr::select
library(openxlsx)
library(DESeq2)
library(data.table)
library(ComplexHeatmap)
library(clusterProfiler)
library(data.table)

# read data
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# read in data ----
mat = readRDS("rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.rds")
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")  %>%
  mutate(drug = paste0(Drug_Treatment, "_", Concentration)) 
nrow(meta)

# analyze distribution of random with best DMSO ----
DE_count_500_sig = fread("rd1_rd2_analysis/de/dmso_500_batch_step2_padj_0.1_abslfc_1_run3.csv")
#DE_count_500_sig = fread("rd1_rd2_analysis/de/dmso_500_batch_step2_padj_0.05_abslfc_1_run3.csv")

n_de_run = DE_count_500_sig %>%
  group_by(run) %>%
  summarise(count=n()) %>%
  ungroup() 

DE_count_500_sig = DE_count_500_sig %>%
  left_join(n_de_run, by="run")

h = hist(n_de_run$count,breaks=100, 
         col="gray", labels = FALSE, main = 'Frequency of number DE genes per rand_run')

### 4.2.2 DE counts statistical summary
quantiles_table<-data.frame(t(round(quantile(n_de_run$count,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),stringsAsFactors = F,check.names = F) 
quantiles_table%>% data.table(caption = "Second 500 random run quantile of numDEGs summary") 

batches<-DE_count_500_sig %>% 
  dplyr::select(batch) %>% 
  distinct() %>% 
  arrange(batch) %>% 
  .$batch

quantiles_stats<-NULL
for(h in 1:length(batches)) {
  quantiles_stats<-rbind(quantiles_stats,data.frame(batch_id=batches[h],t(round(quantile(subset(DE_count_500_sig,DE_count_500_sig$batch==batches[h],select='count')$count,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),stringsAsFactors = F,check.names = F))
}
data.table(quantiles_stats,caption = "Second 500 run quantile of numDEGs summary by batch")

# analyze the active ----
mat = readRDS("rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.rds")
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")  %>%
  mutate(drug = paste0(Drug_Treatment, "_", Concentration)) 
#thresh = 171 # p < 0.1
thresh = 94 # p<0.05

nde = fread("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run3_ndeg.csv")
colnames(nde)
res_filter = fread("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run3.csv")

hist(nde$sig_padj_0.1_abslfc_1,breaks=100, 
         col="gray", labels = FALSE, main = 'Frequency of number DE genes per active') 

abline(v = thresh, col = 'red', lwd = 2, lty = 'dashed')

active_res = meta %>%
  left_join(nde, by="drug") %>%
  mutate(active = ifelse(sig_padj_0.1_abslfc_1 > thresh,1,0))

t = as.data.frame.matrix(table(active_res$drug, active_res$active))
t$treatment = rownames(t)
write.xlsx(t, "rd1_rd2_analysis/de/active_compound_list_9jul24.xlsx")

