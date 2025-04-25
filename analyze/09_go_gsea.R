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

# read data
res = fread("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2.csv")

res   = res  %>%
  separate(gene, into = c("ens","gene"), sep="\\|", remove=F) %>%
  mutate(dir = ifelse(log2FoldChange < 0,"dn","up"))

unique(res$drug)
min(abs(res$log2FoldChange))


ck<- compareCluster(geneCluster = gene~drug + dir,
                    data = res,
                    OrgDb = org.Hs.eg.db,
                    keyType="SYMBOL",
                    fun = "enrichGO",
                    ont="BP")

saveRDS(ck, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp_dir_29jul24.rds")


# active only  ---- 
ck = readRDS("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp.rds")
write.xlsx(ck, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp.xlsx")


active_res   = active_res  %>%
  separate(gene, into = c("ens","gene"), sep="\\|", remove=F)

ck<- compareCluster(geneCluster = gene~drug,
                    data = all_res,
                    OrgDb = org.Hs.eg.db,
                    keyType="SYMBOL",
                    fun = "enrichGO",
                    ont="BP")

saveRDS(ck, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp.rds")
ck = readRDS("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp.rds")
write.xlsx(ck, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp.xlsx")


# reactome ----
convert = clusterProfiler::bitr(active_res$ens, fromType="ENSEMBL", toType="ENTREZID",OrgDb = org.Hs.eg.db)
head(convert)

active_res = active_res %>%
  left_join(convert, by=c("ens" = "ENSEMBL")) %>%
  dplyr::filter(!is.na(ENTREZID))

head(active_res)

ck<- compareCluster(geneCluster = ENTREZID~dir,
                    data = pell_deg_filter ,
                    organism = "human",
                    fun = "enrichPathway",
                    readable = T)



# gsea -----
all_gene_sets = msigdbr(species = "Homo sapiens")%>%
  dplyr::select(gs_name, human_ensembl_gene)

hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, human_ensembl_gene)


# data(geneList, package="DOSE")
# head(geneList)
# 
# y <- GSEA(geneList, 
#           TERM2GENE = C3_t2g)


ff = list.files(path = "rd1_rd2_analysis/de/individual_drug/", pattern = "*_run2.csv")
# https://stephenturner.github.io/deseq-to-fgsea/

dd = gsub("_run2.csv","",ff)

all_list <- vector(mode="list", length=length(dd))
names(all_list) <- dd
ff
for (f in ff){
  d = gsub("_run2.csv","",f)
  res = fread(paste0("rd1_rd2_analysis/de/individual_drug/",f)) %>%
    separate(V1, into = c("ens","gene"), sep="\\|", remove=F) %>%
    dplyr::select(ens,log2FoldChange, pvalue) %>%
    na.omit() %>%
    mutate(score = log2FoldChange*-10*log(pvalue)) %>%
    arrange(-score) %>%
    dplyr::select('ens','score') %>%
    distinct()
  ranks <- deframe(res)
  all_list[d] = list(ranks)
  
  # y <- GSEA(ranks, 
  #           TERM2GENE =  all_gene_sets)
  # 
  # y_readable = setReadable(y, OrgDb=org.Hs.eg.db, keyType = "ENSEMBL")
  # saveRDS(y_readable, paste0("rd1_rd2_analysis/de/gsea/",d, "_all.rds"))
}


ck<- compareCluster(geneClusters=all_list,
                    fun = "GSEA",
                    TERM2GENE =  all_gene_sets)

saveRDS(ck, "rd1_rd2_analysis/de/gsea/all.rds")
