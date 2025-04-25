LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
BiocManager::install("apeglm")
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(apeglm)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# read data
mat= readRDS("rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.rds")
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name

# what are the best and worst well ----
# 1st time, based on 100 reps
# best_wells = c("C08|DMSO|rd2","F05|DMSO|rd2","F11|DMSO|rd1","E09|DMSO|rd1")
# worst_wells=c("H01|DMSO|rd1","E05|DMSO|rd2")

# 2nd time, based on 500 reps
best_wells = c("F05|DMSO|rd2","C02|DMSO|rd2","H07|DMSO|rd1","C01|DMSO|rd1")
worst_wells=c("C07|DMSO|rd2","G01|DMSO|rd1")

dmso <- meta %>%
  dplyr::filter(!Drug_Treatment == "Empty well") %>%
  dplyr::filter( Drug_Treatment == "DMSO" ) %>%
  dplyr::filter(!(sample_name %in% c(worst_wells, best_wells)))


dmso_best <- meta %>%
  dplyr::filter(!Drug_Treatment == "Empty well") %>%
  dplyr::filter( Drug_Treatment == "DMSO" ) %>%
  dplyr::filter(sample_name %in% best_wells) %>%
  .$sample_name

# step 1 analyze dmso ---- 
n_run=500
set.seed(1)
full_res = NULL
for (b in unique(dmso$batch)){
  for (i in 1:n_run){
    
    dmso_batch <- dmso %>% 
      filter(batch == b) 
    
    group_a = sample(dmso_batch$sample_name, 4)
    print(group_a)
    
    group_b = best_wells
    
    mat_i = mat[,c(group_a,group_b)]
    
    meta_i = meta %>%
      dplyr::filter(sample_name %in% c(group_a,group_b)) %>%
      mutate(group = ifelse(sample_name %in% group_a,'a','b'))
    
    mat_i = mat_i[,rownames(meta_i)]
    
    dds <- DESeqDataSetFromMatrix(countData = mat_i, colData = meta_i, design = ~ group)
    dds$group <- relevel(dds$group, "a")
    
    dds <- DESeq(dds)
    #res <- results(dds)
    #resultsNames(dds)
    res <- lfcShrink(dds, coef="group_b_vs_a", type="apeglm")
    
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

write.csv(full_res, "rd1_rd2_analysis/de/dmso_500_batch_step2_run3.csv")

# relaxed threshold
full_res_sig = full_res %>%
  dplyr::filter(padj<0.1, abs(log2FoldChange) > 1)

write.csv(full_res_sig, "rd1_rd2_analysis/de/dmso_500_batch_step2_padj_0.1_abslfc_1_run3.csv", row.names = F)

# more stringent threshold
full_res_sig = read.csv("rd1_rd2_analysis/de/dmso_500_batch_step2_padj_0.1_abslfc_1_run3.csv")

full_res_sig = full_res_sig %>%
  dplyr::filter(padj<0.05)

write.csv(full_res_sig, "rd1_rd2_analysis/de/dmso_500_batch_step2_padj_0.05_abslfc_1_run3.csv", row.names = F)
