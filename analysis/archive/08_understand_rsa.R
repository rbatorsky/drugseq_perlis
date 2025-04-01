LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(data.table)


# read data
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')


# save and plot Step1 results -----
DE_count_500 = read.csv("rd1_rd2_analysis/de/dmso_sig_100_batch.csv")

nde = DE_count_500 %>%
  group_by(run) %>%
  summarise(count=n()) %>%
  ungroup() 


nde

DE_count_500 = DE_count_500 %>%
  left_join(nde, by="run") %>%
  dplyr::rename(run_no = run)


head(DE_count_500)

h = hist(DE_count_500$count,breaks=100, 
         col="gray", labels = FALSE, main = 'Frequency of number DE genes per rand_run+batch (N=1500)_Step_1')


# DE counts statistical summary ----
RC_to_AC_num_DE_seprow<-separate_rows(data = DE_count_500,
                                      'group_a',sep = '_')

view(RC_to_AC_num_DE_seprow)
RC_to_AC_num_DE_seprow$batch_plate_well<-paste0(RC_to_AC_num_DE_seprow$batch,'_',RC_to_AC_num_DE_seprow$group_a)

DE_summary1 <- data.frame(t(round(quantile(DE_count_500$count,
                                           probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),
                          stringsAsFactors = F,
                          check.names = F) 
DE_summary1 %>% 
  write.csv("rd1_rd2_analysis/de/summary_first_500run2.csv",row.names = FALSE)

DE_summary1 %>% 
  data.table(caption = "First 500 random run quantile of numDEGs summary")


batches<-DE_count_500 %>% 
  select(batch) %>% 
  distinct() %>% 
  arrange(batch) %>% 
  .$batch
batches

quantiles_stats<-NULL
for(h in 1:length(batches)) {
  quantiles_stats<-rbind(quantiles_stats,
                         data.frame(batch_id=batches[h],
                                    t(round(quantile(subset(DE_count_500,
                                                            DE_count_500$batch==batches[h],select='count')$count,
                                                     probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),
                                    stringsAsFactors = F,
                                    check.names = F))
}

quantiles_stats %>%  data.table(caption = "First 500 random run quantile of numDEGs summary by batch")

# select bad and best DMSO-----

rsa_res<-NULL
source('scripts/RNAi-package.r')

for(h in 1:length(batches)) {
  print(h)
  batchid<-batches[h]
  #RC_to_AC_num_DE_seprow$batch<-as.integer(RC_to_AC_num_DE_seprow$batch)
  RC_to_AC_num_DE_seprow_batch <- RC_to_AC_num_DE_seprow %>% 
    filter(batch==batchid) %>% 
    mutate(run_plate_well=paste0(batch_plate_well,"_",run_no)) %>% 
    mutate(log2FC=count/median(count,na.rm = T)) %>%  # log2FC is the 
    mutate(rz_score=(log2FC-median(log2FC,na.rm = T))/mad(log2FC,na.rm = T)) %>% # normalize to median
    arrange(count)
  
  #RSA on ranked values, rz_score transformed
  rsa_res_batch<-score.screen.ZY.RQz(data = RC_to_AC_num_DE_seprow_batch,
                                     gene.column = 'batch_plate_well',
                                     value.column = 'rz_score',
                                     log.column = 'log2FC')                       
  rsa_res<-rbind(rsa_res,data.frame(batch_id=batches[h],rsa_res_batch,stringsAsFactors = F))
}