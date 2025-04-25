LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(data.table)
# read data
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# read in random runs ----
DE_count_500 = read.csv("rd1_rd2_analysis/de/dmso_500_batch_notrim_padj_0.1_abslfc_1.csv") 

# which genes are often deg in dmso -----
gene_count = DE_count_500 %>%
  group_by(gene) %>%
  summarise(count=n()) %>%
  ungroup() 

# which samples have most de per random run -----
n_de_run = DE_count_500 %>%
  group_by(run) %>%
  summarise(count=n()) %>%
  ungroup() 

h = hist(n_de_run$count,breaks=100, 
         col="gray", labels = FALSE, main = 'Frequency of number DE genes per rand_run')


quantiles_table<-data.frame(t(round(quantile(n_de_run$count,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),stringsAsFactors = F,check.names = F) 
quantiles_table%>% data.table(caption = "Second 500 random run quantile of numDEGs summary") 


DE_count_500 = DE_count_500 %>%
  left_join(n_de_run, by="run") %>%
  dplyr::rename(run_no = run)

h = hist(DE_count_500$count,breaks=100, 
     col="gray", labels = FALSE, main = 'Frequency of number DE genes per rand_run')

# DE counts statistical summary ----
RC_to_AC_num_DE_seprow<-separate_rows(data = DE_count_500,
                                      'group_a',sep = '_')
colnames(RC_to_AC_num_DE_seprow)

data = RC_to_AC_num_DE_seprow %>%
  select(run_no, group_a, count) %>%
  distinct()

n_rand_run_per_well = data %>%
  group_by(group_a) %>%
  summarise(n_rand_run=n())

nde = data %>%
  group_by(group_a) %>%
  summarise(n_de_total=sum(count)) %>%
  left_join(n_rand_run_per_well, by="group_a") %>%
  mutate(nde_per_run = n_de_total/n_rand_run) %>% 
  arrange(-nde_per_run)

RC_to_AC_num_DE_seprow$batch_plate_well<-paste0(RC_to_AC_num_DE_seprow$batch,'_',RC_to_AC_num_DE_seprow$group_a)

head(RC_to_AC_num_DE_seprow)

quantile(nde$nde_per_run,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1))
         
DE_summary1 <- data.frame(t(round(quantile(nde$n_rand_run,
                                           probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),
                          stringsAsFactors = F,
                          check.names = F) 
DE_summary1

DE_summary1 %>% 
  write.csv("rd1_rd2_analysis/de/summary_first_500run2_23aug24.csv",row.names = FALSE)

# old = read.csv("rd1_rd2_analysis/de/summary_first_500run2.csv")
# old

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
    mutate(log2FC=count/median(count,na.rm = T)) %>% 
    mutate(rz_score=(log2FC-median(log2FC,na.rm = T))/mad(log2FC,na.rm = T)) %>% # normalize to median
    arrange(count)
  
  #RSA on ranked values, rz_score transformed
  rsa_res_batch<-score.screen.ZY.RQz(data = RC_to_AC_num_DE_seprow_batch,
                                     gene.column = 'batch_plate_well',
                                     value.column = 'rz_score',
                                     log.column = 'log2FC')                       
  rsa_res<-rbind(rsa_res,data.frame(batch_id=batches[h],rsa_res_batch,stringsAsFactors = F))
}

rsa_res[!is.finite(rsa_res$logP_RSA_Down),'logP_RSA_Down']<-0
rsa_res[!is.finite(rsa_res$logP_RSA_Up),'logP_RSA_Up']<-0

#prepare the plate_barcode and well
rsa_res<-separate(rsa_res,col = 'batch_plate_well',into =c('batch','sample_name'),sep = '_',remove = F)

#rsa_res<-separate(rsa_res,col = 'rest',into =c('plate_barcode','plate_well'),sep = ' ',remove = T)
#rsa_res$batch<-NULL
#rsa_res$plate_barcode<-paste0('VH02001',rsa_res$plate_barcode)
#rsa_res$sector<-1

# for(i in 5:15){
#  rsa_res[,i]<-as.numeric(rsa_res[,i])
# }

dim(subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3))) 
rsa_res_rmv<-subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3))

rsa_res_rmv_summr<-rsa_res_rmv %>% 
  group_by(batch) %>% 
  dplyr::summarize(num_DMSO_wells_rmv_per_batch=n())

rsa_res_rmv_summr

# best DMSO
# for(i in 5:15){
#   rsa_res[,i]<-as.numeric(rsa_res[,i])
# }
dim(rsa_res) 


rsa_res<-rsa_res[order(rsa_res$logP_RSA_Down,
                       decreasing = F),]

#rsa_res<-rsa_res[order(rsa_res$plate_barcode),]
rsa_res_best_DMSO_wells<-NULL

# for each BATCH, select the 3 wells with best
for(h in 1:length(batches)) {
  print(h)
  batchid<-batches[h]
  rsa_res_1plate<-subset(rsa_res,
                         rsa_res$batch==batchid)
  
  rsa_res_1plate$order <- seq(1,nrow(rsa_res_1plate),by=1)
  rsa_res_1plate<-rsa_res_1plate[1:3,] # 3 per plate
  rsa_res_best_DMSO_wells<-rbind(rsa_res_best_DMSO_wells,rsa_res_1plate)
}

write.csv(rsa_res_rmv, 
          file = 'rd1_rd2_analysis/de/bad_DMSO_500remove.csv', 
          row.names = FALSE)
write.csv(rsa_res_best_DMSO_wells,
          file ='rd1_rd2_analysis/de/best_DMSO_500keep.csv',
          row.names = FALSE) 

rsa_res_rmv
rsa_res_best_DMSO_wells


# # process the results -----
source('scripts/RNAi-package.r')

DE_count_500 = read.csv("rd1_rd2_analysis/de/dmso_sig_100_batch.csv") %>%
  mutate(rz_score=(log2FoldChange-median(log2FoldChange,na.rm = T))/mad(log2FoldChange,na.rm = T))

head(DE_count_500)
rsa_res=NULL
for (b in unique(DE_count_500$batch)){
  
  DE_count_500_b = DE_count_500 %>%
    dplyr::filter(batch == b)
  
  rsa_res_batch<-score.screen.ZY.RQz(data = DE_count_500_b,
                                     gene.column = 'X',
                                     value.column = 'rz_score',
                                     log.column = 'log2FoldChange')
  
  rsa_res_batch[!is.finite(rsa_res_batch$logP_RSA_Down),'logP_RSA_Down']<-0
  rsa_res_batch[!is.finite(rsa_res_batch$logP_RSA_Up),'logP_RSA_Up']<-0
  
  rsa_res<-rbind(rsa_res,data.frame(batch=b,rsa_res_batch,stringsAsFactors = F))
}

view(rsa_res)


#prepare the plate_barcode and well
rsa_res<-separate(rsa_res,col = 'batch_plate_well',into =c('batch','rest'),sep = '_',remove = F)

rsa_res<-separate(rsa_res,col = 'rest',into =c('plate_barcode','plate_well'),sep = ' ',remove = T)
rsa_res$batch<-NULL
rsa_res$plate_barcode<-paste0('VH02001',rsa_res$plate_barcode)
rsa_res$sector<-1


for(i in 5:15){
  rsa_res[,i]<-as.numeric(rsa_res[,i])
}
dim(rsa_res)

dim(subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3)))
rsa_res_rmv<-subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3))


rsa_res_rmv_summr<-rsa_res_rmv %>% group_by(batch_id,sector)%>% dplyr::summarize(num_DMSO_wells_rmv_per_batch=n())
rsa_res_rmv_summr

# best DMSO
for(i in 5:15){
  rsa_res[,i]<-as.numeric(rsa_res[,i])
}
dim(rsa_res)


rsa_res<-rsa_res[order(rsa_res$logP_RSA_Down,decreasing = F),]
rsa_res<-rsa_res[order(rsa_res$plate_barcode),]
rsa_res_best_DMSO_wells<-NULL
# for each plate, select the 3 wells with best

uniq_plate_barcodes<-unique(rsa_res$plate_barcode)
for(h in 1: length(uniq_plate_barcodes)){
  rsa_res_1plate<-subset(rsa_res,rsa_res$plate_barcode==uniq_plate_barcodes[h])
  rsa_res_1plate$order <- seq(1,nrow(rsa_res_1plate),by=1)
  rsa_res_1plate<-rsa_res_1plate[1:3,] # 3 per plate
  rsa_res_best_DMSO_wells<-rbind(rsa_res_best_DMSO_wells,rsa_res_1plate)
}

write.csv(rsa_res_rmv,file = file.path(out_data_dir,'bad_DMSO_500remove.csv'),row.names = FALSE)
write.csv(rsa_res_best_DMSO_wells,file =file.path(out_data_dir,'best_DMSO_500keep.csv'),row.names = FALSE)



## 2.3 select bad and best DMSO 
rsa_res<-NULL
source(file.path(analysis_dir,'RNAi-package.r'))
for(h in 1:length(batches)) {
  #print(h)
  batchid<-batches[h]
  RC_to_AC_num_DE_seprow$batch<-as.integer(RC_to_AC_num_DE_seprow$batch)
  RC_to_AC_num_DE_seprow_batch <- RC_to_AC_num_DE_seprow %>% 
    filter(batch==batchid) %>% 
    mutate(run_plate_well=paste0(batch_plate_well,"_",run_no)) %>% 
    mutate(log2FC=count/median(count,na.rm = T)) %>% 
    mutate(rz_score=(log2FC-median(log2FC,na.rm = T))/mad(log2FC,na.rm = T)) %>% 
    arrange(count)
  
  #RSA on ranked values, rz_score transformed
  rsa_res_batch<-score.screen.ZY.RQz(data = RC_to_AC_num_DE_seprow_batch,
                                     gene.column = 'batch_plate_well',
                                     value.column = 'rz_score',
                                     log.column = 'log2FC')                       
  rsa_res<-rbind(rsa_res,data.frame(batch_id=batches[h],rsa_res_batch,stringsAsFactors = F))
  
  
  
}

rsa_res[!is.finite(rsa_res$logP_RSA_Down),'logP_RSA_Down']<-0
rsa_res[!is.finite(rsa_res$logP_RSA_Up),'logP_RSA_Up']<-0

#prepare the plate_barcode and well
rsa_res<-separate(rsa_res,col = 'batch_plate_well',into =c('batch','rest'),sep = '_',remove = F)
rsa_res<-separate(rsa_res,col = 'rest',into =c('plate_barcode','plate_well'),sep = ' ',remove = T)
rsa_res$batch<-NULL
rsa_res$plate_barcode<-paste0('VH02001',rsa_res$plate_barcode)
rsa_res$sector<-1


for(i in 5:15){
  rsa_res[,i]<-as.numeric(rsa_res[,i])
}
dim(rsa_res) 

dim(subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3))) 
rsa_res_rmv<-subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3))


rsa_res_rmv_summr<-rsa_res_rmv %>% group_by(batch_id,sector)%>% dplyr::summarize(num_DMSO_wells_rmv_per_batch=n())
rsa_res_rmv_summr

# best DMSO
for(i in 5:15){
  rsa_res[,i]<-as.numeric(rsa_res[,i])
}
dim(rsa_res) 


rsa_res<-rsa_res[order(rsa_res$logP_RSA_Down,decreasing = F),]
rsa_res<-rsa_res[order(rsa_res$plate_barcode),]
rsa_res_best_DMSO_wells<-NULL
# for each plate, select the 3 wells with best

uniq_plate_barcodes<-unique(rsa_res$plate_barcode)
for(h in 1: length(uniq_plate_barcodes)){
  rsa_res_1plate<-subset(rsa_res,rsa_res$plate_barcode==uniq_plate_barcodes[h])
  rsa_res_1plate$order <- seq(1,nrow(rsa_res_1plate),by=1)
  rsa_res_1plate<-rsa_res_1plate[1:3,] # 3 per plate
  rsa_res_best_DMSO_wells<-rbind(rsa_res_best_DMSO_wells,rsa_res_1plate)
}

write.csv(rsa_res_rmv,file = file.path(out_data_dir,'bad_DMSO_500remove.csv'),row.names = FALSE)
write.csv(rsa_res_best_DMSO_wells,file =file.path(out_data_dir,'best_DMSO_500keep.csv'),row.names = FALSE) 
