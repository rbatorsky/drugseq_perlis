library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(openxlsx)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# read in -----
mat1 = readRDS("20240312_Joshua_Bowen/analysis/processed/umi.notrim.nondedup.counts.rds")
mat2 = readRDS("20240530_Joshua_Bowen/analysis/processed/rd2.umi.nondedup.counts.noempty.rds")

rename_df = data.frame(n1 = colnames(mat1), n2 = 'rd1')
rename_df$n3 = paste0(rename_df$n1,'|',rename_df$n2)
head(rename_df)
colnames(mat1) = rename_df$n3

rename_df = data.frame(n1 = colnames(mat2), n2 = 'rd2')
rename_df$n3 = paste0(rename_df$n1,'|',rename_df$n2)
colnames(mat2) = rename_df$n3

mat = cbind(mat1,mat2)
saveRDS(mat, "rd1_rd2_analysis/rd1_rd2_umi.notrim.nondedup.counts.rds")

meta1 = read.xlsx("20240312_Joshua_Bowen/Drug-Seq_Rd1_Treatment+Barcode_format.xlsx")
meta2 = read.xlsx("20240530_Joshua_Bowen/Drug-Seq_Rd2_Treatment+Barcode_format.xlsx")
meta1$batch = "rd1"
meta2$batch = "rd2"

meta = rbind(meta1,meta2)
meta$sample_name = paste0(meta$sample_name,"|",meta$batch)
head(meta)
write.xlsx(meta, "rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
