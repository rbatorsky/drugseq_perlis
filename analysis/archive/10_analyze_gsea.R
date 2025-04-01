LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
library(openxlsx)
library(DESeq2)
library('org.Hs.eg.db')
library(clusterProfiler)
library(data.table)
library(msigdbr)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# meta ----
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")  %>%
  mutate(drug = paste0(Drug_Treatment, "_", Concentration)) 

# gsea ----
gsea = readRDS("rd1_rd2_analysis/de/gsea/all_readable.rds")

view(gsea)


dotplot(ck_readable)
