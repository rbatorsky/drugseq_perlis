LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
library(openxlsx)
library(DESeq2)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# read data
mat = readRDS("rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.rds")
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name

dmso = meta %>%
  dplyr::filter(Drug_Treatment=="DMSO")

# step 1 analyze dmso ---- 
n_run=500
set.seed(1)
full_res = NULL
for (b in unique(meta$batch)){
  for (i in 1:n_run){
    
    dmso_batch <- dmso %>% 
      filter(batch == b) 
    
    group_a = sample(dmso_batch$sample_name, 4)
    print(group_a)
    
    group_b = dmso_batch %>%
      dplyr::filter(!(sample_name %in% group_a))
    
    group_b = sample(group_b$sample_name, 4)
    
    mat_i = mat[,c(group_a,group_b)]
    meta_i = meta %>%
      dplyr::filter(sample_name %in% c(group_a,group_b)) %>%
      mutate(group = ifelse(sample_name %in% group_a,'a','b'))
    
    mat_i = mat_i[,rownames(meta_i)]
    
    dds <- DESeqDataSetFromMatrix(countData = mat_i, colData = meta_i, design = ~ group)
    dds$group <- relevel(dds$group, "a")
    dds <- DESeq(dds)
    res <- results(dds)
    
    res = data.frame(res) 
    
    if(nrow(res)>0){
      res$batch = b
      res$run = i
      res$group_a = paste(group_a, collapse="_")
      res = res %>%
        rownames_to_column("gene")
      res = res %>%
       dplyr::select(gene, everything())
      if(is.null(full_res)){
        full_res = res
      }else{
        full_res = rbind(full_res, res)
      }
    }
  }
}

write.csv(full_res, "rd1_rd2_analysis/de/dmso_500_batch.csv")

library(data.table)
library(tidyverse)
full_res = fread("rd1_rd2_analysis/de/dmso_500_batch.csv")

full_res_sig = full_res %>%
  dplyr::filter(padj<0.1, abs(log2FoldChange) > 1)

write.csv(full_res_sig, "rd1_rd2_analysis/de/dmso_500_batch_padj_0.1_abslfc_1.csv", row.names = F)

