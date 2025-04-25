library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(openxlsx)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/')

# read in -----
# ma_nodedup = read.table("analysis/star_rd2/Solo.out/Gene/processed/umi.nondedup.counts.txt")
# head(ma_nodedup)

mat = read.table("analysis/star_rd3/Solo.out/Gene/processed/umi.1mm_all_dedup.counts.txt")
meta = read.xlsx("Drug-Seq_Rd3_Treatment+Barcode.xlsx")

# format meta ----
table(meta$Drug_Treatment, meta$Concentration)

meta = meta %>%
  mutate(Drug_Treatment = gsub('\\)','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\(','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\-','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\+','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('_$','',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('^_','',Drug_Treatment)) %>%
  mutate(sample_name = paste0(Plate_Well,"|",Drug_Treatment)) 


rownames(meta) = meta$sample_name

head(meta)
table(meta$Drug_Treatment)

# format mat -----
length(colnames(mat))
colnames(mat)
head(meta)
head(mat) 

# change the column names from barcode to sample_id|treatment
convert_names = data.frame(Barcode = colnames(mat))
head(convert_names)

convert_names = convert_names %>%
  left_join(meta, by="Barcode")


head(convert_names)

colnames(mat) = convert_names$sample_name
head(mat)

# remove empty well for better visualization
rm = which(colnames(mat) %in% c("H12|Empty well", "G12|Empty well"))
rm

meta = meta %>%
  dplyr::filter(!Drug_Treatment == "Empty well")
mat =mat[, -rm]

read_per_sample = colSums(mat)

hist(read_per_sample)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(mat)
gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(gene_list, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

saveRDS(mat_w_genes, "analysis/star_rd3/Solo.out/Gene/processed/umi.1mm_all_dedup.counts.noempty.rds")



# start deseq2 -----
all(colnames(mat_w_genes) %in% rownames(meta))
mat_w_genes = mat_w_genes[,rownames(meta)]
all(colnames(mat_w_genes) == rownames(meta))
ncol(mat_w_genes)

dds <- DESeqDataSetFromMatrix(countData = mat_w_genes, colData = meta, design = ~ Drug_Treatment)
dds$Drug_Treatment <- relevel(dds$Drug_Treatment, "DMSO")

dds <- DESeq(dds)
saveRDS(dds, "analysis/star_rd2/Solo.out/Gene/processed/umi.1mm_all_dedup.counts.noempty.rds")
vst <- vst(dds, blind=TRUE)

# dimension reduction ----
dds = readRDS("analysis/processed/umi.1mm_all_dedup.counts.noempty.rds")

# with vst
vst_mat <- assay(vst)

## Compute pairwise correlation values
vst_cor <- cor(vst_mat)   

library(ComplexHeatmap)

meta_heatmap = meta %>%
  dplyr::select(Drug_Treatment)

top_ha = HeatmapAnnotation(df = meta_heatmap)

Heatmap(vst_cor,
        top_annotation = top_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))

pcaData <- plotPCA(vst, intgroup=c("Drug_Treatment","Plate_Well"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData = pcaData %>%
  separate(name, into=c("Plate_Well","Drug_Treatment"), remove=F, sep="\\|")

head(pcaData)
ggplot(pcaData, aes(PC1, PC2, color=Drug_Treatment, label=name)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + 
  geom_text(hjust=0, vjust=0, size=2)

# dpi=300
# pdf("analysis/plots/hc_annotate_legend.pdf",
#     width = 12,
#     height = 10 )
# dev.off()

# with normalized counts
# normalized_counts <- data.frame(counts(dds, normalized=TRUE))
# normalized_counts_cor <- cor(normalized_counts)   
# 
# library(ComplexHeatmap)
# 
# meta_heatmap = meta %>%
#   dplyr::select(Drug_Treatment)
# 
# top_ha = HeatmapAnnotation(df = meta_heatmap)
# 
# Heatmap(normalized_counts_cor,
#         top_annotation = top_ha,
#         column_names_gp = grid::gpar(fontsize = 6),
#         row_names_gp = grid::gpar(fontsize = 6))
# 
# 


# highly variable genes -----

var_genes <- apply(vst_mat, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:75]
head(select_var)
vst_hv <- vst_mat[select_var,]
dim(vst_hv)
head(vst_hv)
rownames(vst_hv) = gsub(".*\\|","",rownames(vst_hv))

meta_heatmap = meta %>%
  dplyr::select(Drug_Treatment)

bottom_ha = HeatmapAnnotation(df = meta_heatmap)

head(vst_hv)
vst_hv_scaled = t(scale(t(vst_hv)))

Heatmap(vst_hv_scaled,
        bottom_annotation = bottom_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))




# back to analysis -----
normalized_counts <- data.frame(counts(dds, normalized=TRUE))

# look at two example conditions that group well -----
c1 = "GSK805"
n1 = paste0("Drug_Treatment_",c1, '_vs_DMSO')
r1 <- results(dds, name = n1)
r1a <- lfcShrink(dds, coef = n1, type="apeglm")

c2 = "Vorinostat"
n2 = paste0("Drug_Treatment_",c2, '_vs_DMSO')
r2 <- results(dds, name = n2)
r2a <- lfcShrink(dds, coef = n2, type="apeglm")

r1adf = data.frame(r1a) %>%
  rownames_to_column("gene") 

write.xlsx(r1adf, "analysis/processed/r1adf.xlsx")

r2adf = data.frame(r2a) %>%
  rownames_to_column("gene") 

write.xlsx(r2adf, "analysis/processed/r2adf.xlsx")


r1adf_sig = r1adf %>%
  dplyr::filter(padj< 0.05 & abs(log2FoldChange) > 1.5)%>%
  mutate(treatment = c1)

r2adf_sig = r2adf %>%
  dplyr::filter(padj< 0.05 & abs(log2FoldChange) > 1.5) %>%
  mutate(treatment = c2)

results = rbind(r1adf_sig, r2adf_sig) %>%
  separate(gene, into = c("ens","gene"), sep="\\|", remove=F)

table(results$treatment)
head(results)

library(clusterProfiler)
library(org.Hs.eg.db)

ck<- compareCluster(geneCluster = gene~treatment,
                    data = results ,
                    OrgDb = org.Hs.eg.db,
                    keyType="SYMBOL",
                    fun = "enrichGO",
                    ont="BP")

saveRDS(ck, "analysis/processed/go_bp_test.rds")

dotplot(ck, show = 10,font.size = 10)



# all results ----
results_names = resultsNames(dds)
results_names = results_names[-1]

all_results = NULL
for (r in results_names){
  res <- results(dds, name = r)
  res = data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(name = r)
  
  if(is.null(all_results)){
    all_results = res
  }else{
    all_results = rbind(all_results, res)
  }
}

write.xlsx(all_results, file = "analysis/processed/results_vs_dmso.nodedup.xlsx")


all_results$name = gsub("Drug_Treatment_","",all_results$name)

r1 = all_results %>%
  dplyr::filter(name == paste0(c1, '_vs_DMSO'))

r1_sig = r1 %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1.5)

nrow(r1_sig)


# # results, first round  ----
# 
# # compare genotype only in three conditions
# ko_un_b6_un <- results(dds, contrast=c("cond_geno","unstim_T169A","unstim_B6"), alpha = 0.05)
# ko_lps_b6_lps  <- results(dds, contrast=c("cond_geno","LPS_T169A","LPS_B6"), alpha = 0.05)
# ko_lps5z7_b6_lps5z7  <- results(dds, contrast=c("cond_geno","LPS5z7_T169A","LPS5z7_B6"), alpha = 0.05)
# 
# # compare condition only in two genotypes
# ko_lps5z7_ko_lps  <- results(dds, contrast=c("cond_geno","LPS5z7_T169A","LPS_T169A"), alpha = 0.05)
# b6_lps5z7_b6_lps  <- results(dds, contrast=c("cond_geno","LPS5z7_B6","LPS_B6"), alpha = 0.05)
# 
# # compare condition only in two genotypes with unstim 
# ko_lps5z7_ko_unstim  <- results(dds, contrast=c("cond_geno","LPS5z7_T169A","unstim_T169A"), alpha = 0.05)
# ko_lps_ko_unstim  <- results(dds, contrast=c("cond_geno","LPS_T169A","unstim_T169A"), alpha = 0.05)
# b6_lps5z7_b6_unstim  <- results(dds, contrast=c("cond_geno","LPS5z7_B6","unstim_B6"), alpha = 0.05)
# b6_lps_b6_unstim  <- results(dds, contrast=c("cond_geno","LPS_B6","unstim_B6"), alpha = 0.05)
# 
# 
# #shrink
# library('ashr')
# ko_un_b6_un_s <- lfcShrink(dds, contrast=c("cond_geno","unstim_T169A","unstim_B6"), res=ko_un_b6_un ,type='ashr')
# ko_lps_b6_lps_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS_T169A","LPS_B6"), res=ko_lps_b6_lps,type='ashr')
# ko_lps5z7_b6_lps5z7_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS5z7_T169A","LPS5z7_B6"), res=ko_lps5z7_b6_lps5z7,type='ashr')
# ko_lps5z7_ko_lps_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS5z7_T169A","LPS_T169A"), res=ko_lps5z7_ko_lps,type='ashr')
# b6_lps5z7_b6_lps_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS5z7_B6","LPS_B6"), res=b6_lps5z7_b6_lps,type='ashr')
# 
# ko_lps5z7_ko_unstim_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS5z7_T169A","unstim_T169A"), res=ko_lps5z7_ko_unstim ,type='ashr')
# ko_lps_ko_unstim_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS_T169A","unstim_T169A"), res=ko_lps_ko_unstim ,type='ashr')
# b6_lps5z7_b6_unstim_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS5z7_B6","unstim_B6"), res=b6_lps5z7_b6_unstim ,type='ashr')
# b6_lps_b6_unstim_s  <- lfcShrink(dds, contrast=c("cond_geno","LPS_B6","unstim_B6"), res=b6_lps_b6_unstim ,type='ashr')
# 
# #format
# ko_un_b6_un_s$type = "ko_un_b6_un"
# ko_lps_b6_lps_s$type = "ko_lps_b6_lps"
# ko_lps5z7_b6_lps5z7_s$type = "ko_lps5z7_b6_lps5z7"
# ko_lps5z7_ko_lps_s$type = "ko_lps5z7_ko_lps"
# b6_lps5z7_b6_lps_s$type = "b6_lps5z7_b6_lps"
# 
# ko_lps5z7_ko_unstim_s$type  <- "ko_lps5z7_ko_un"
# ko_lps_ko_unstim_s$type  <- "ko_lps_ko_un"
# b6_lps5z7_b6_unstim_s$type  <- "b6_lps5z7_b6_un"
# b6_lps_b6_unstim_s$type  <- "b6_lps_b6_un"
# 
# ko_un_b6_un_s= data.frame(ko_un_b6_un_s) %>% rownames_to_column("gene_id")
# ko_lps_b6_lps_s= data.frame(ko_lps_b6_lps_s) %>% rownames_to_column("gene_id")
# ko_lps5z7_b6_lps5z7_s= data.frame(ko_lps5z7_b6_lps5z7_s) %>% rownames_to_column("gene_id")
# ko_lps5z7_ko_lps_s= data.frame(ko_lps5z7_ko_lps_s) %>% rownames_to_column("gene_id")
# b6_lps5z7_b6_lps_s= data.frame(b6_lps5z7_b6_lps_s) %>% rownames_to_column("gene_id")
# 
# ko_lps5z7_ko_unstim_s= data.frame(ko_lps5z7_ko_unstim_s) %>% rownames_to_column("gene_id")
# ko_lps_ko_unstim_s= data.frame(ko_lps_ko_unstim_s) %>% rownames_to_column("gene_id")
# b6_lps5z7_b6_unstim_s= data.frame(b6_lps5z7_b6_unstim_s) %>% rownames_to_column("gene_id")
# b6_lps_b6_unstim_s= data.frame(b6_lps_b6_unstim_s) %>% rownames_to_column("gene_id")
# 
# 
# # analysis with interaction term -----
# 
# # interaction one at a time -----
# rm(meta_int)
# rm(counts_int)
# rm(dds)
# rm(int_s)
# # remove unstim
# meta_int = meta %>%
#   dplyr::filter(cond != "unstim")
# 
# counts_int = counts[,rownames(meta_int)]
# 
# ## check sample order
# all(colnames(counts_int) %in% rownames(meta_int))
# counts_int = counts_int[,rownames(meta_int)]
# all(colnames(counts_int) == rownames(meta_int))
# 
# dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = meta_int, design = ~ geno + cond + cond:geno)
# dds$geno = factor(dds$geno, levels=c("B6","T169A"))
# dds$cond = factor(dds$cond, levels=c("LPS","LPS5z7"))
# 
# dds <- DESeq(dds)
# rld <- rlog(dds, blind=TRUE)
# 
# resultsNames(dds)
# int <- results(dds, name="genoT169A.condLPS5z7")
# 
# #shrink
# library('ashr')
# int_s <- lfcShrink(dds, 
#                    coef="genoT169A.condLPS5z7", 
#                    res=int)
# 
# int_s_1 = data.frame(int_s) %>%
#   mutate(type = "genoT169A.condLPS5z7_LPS") %>%
#   rownames_to_column("gene_id") 
# 
# 
# int_s_1 %>% dplyr::filter(padj < 0.05)
# 
# # remove LPS5z7
# rm(meta_int)
# rm(counts_int)
# rm(dds)
# rm(int_s)
# 
# meta_int = meta %>%
#   dplyr::filter(cond != "LPS5z7")
# 
# counts_int = counts[,rownames(meta_int)]
# 
# ## check sample order
# all(colnames(counts_int) %in% rownames(meta_int))
# counts_int = counts_int[,rownames(meta_int)]
# all(colnames(counts_int) == rownames(meta_int))
# 
# dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = meta_int, design = ~ geno + cond + cond:geno)
# dds$geno = factor(dds$geno, levels=c("B6","T169A"))
# dds$cond = factor(dds$cond, levels=c("unstim","LPS"))
# dds <- DESeq(dds)
# rld <- rlog(dds, blind=TRUE)
# 
# # results ----
# 
# resultsNames(dds)
# int <- results(dds, name="genoT169A.condLPS")
# 
# #shrink
# library('ashr')
# int_s <- lfcShrink(dds, coef="genoT169A.condLPS", res=int)
# 
# 
# int_s_2 = data.frame(int_s) %>%
#   mutate(type = "genoT169A.condLPS_un") %>%
#   rownames_to_column("gene_id") 
# 
# int_s_2 %>% dplyr::filter(padj<0.2)
# 
# # remove LPS
# rm(meta_int)
# rm(counts_int)
# rm(dds)
# rm(int_s)
# 
# meta_int = meta %>%
#   dplyr::filter(cond != "LPS")
# 
# counts_int = counts[,rownames(meta_int)]
# 
# ## check sample order
# all(colnames(counts_int) %in% rownames(meta_int))
# counts_int = counts_int[,rownames(meta_int)]
# all(colnames(counts_int) == rownames(meta_int))
# 
# dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = meta_int, design = ~ geno + cond + cond:geno)
# dds$geno = factor(dds$geno, levels=c("B6","T169A"))
# dds$cond = factor(dds$cond, levels=c("unstim","LPS5z7"))
# dds <- DESeq(dds)
# rld <- rlog(dds, blind=TRUE)
# 
# # results ----
# 
# resultsNames(dds)
# int <- results(dds, name="genoT169A.condLPS5z7")
# 
# #shrink
# library('ashr')
# int_s <- lfcShrink(dds, coef="genoT169A.condLPS5z7", res=int)
# 
# int_s_3 = data.frame(int_s) %>%
#   mutate(type = "genoT169A.condLPS5z7_un") %>%
#   rownames_to_column("gene_id") 
# 
# 
# int_s_3 %>% dplyr::filter(padj<0.2)
# 
# 
# # combine and write out ----
# combined=rbind(ko_un_b6_un_s,
#                ko_lps_b6_lps_s,
#                ko_lps5z7_b6_lps5z7_s,
#                ko_lps5z7_ko_lps_s,
#                b6_lps5z7_b6_lps_s,
#                ko_lps5z7_ko_unstim_s,
#                ko_lps_ko_unstim_s,
#                b6_lps5z7_b6_unstim_s,
#                b6_lps_b6_unstim_s,
#                int_s_1,
#                int_s_2,
#                int_s_3
# )
# 
# combined %>% dplyr::filter(gene_name == "Ift57")
# 
# 
# combined = combined %>% 
#   left_join(gtf, by="gene_id")
# 
# write.xlsx(combined,"/cluster/tufts/poltoraklab/rbator01/ripk1_may2023/analysis/deseq2/deg_all_pairwise_and_int_26may23.xlsx")
# 
# combined_sig=data.frame(combined) %>% 
#   dplyr::filter(padj<0.05) 
# 
# write.xlsx(combined_sig,"/cluster/tufts/poltoraklab/rbator01/ripk1_may2023/analysis/deseq2/deg_all_pairwise_and_int_padj_0.05_26may23.xlsx")
# 
# combined_less_sig=data.frame(combined) %>% 
#   dplyr::filter(padj<0.2) 
# 
# write.xlsx(combined_less_sig,"/cluster/tufts/poltoraklab/rbator01/ripk1_may2023/analysis/deseq2/deg_all_pairwise_and_int_padj_0.2_26may23.xlsx")
# 
# table(combined_less_sig$type)
# 
# #combined_old_sig = read.xlsx("/cluster/tufts/poltoraklab/rbator01/ripk1_may2023/analysis/deseq2/deg_all_pairwise_and_int_padj_0.05.xlsx")
# #table(combined_old_sig$type)
# 
# 
# # make a heatmap of the sig DEG ----
# # top20 per condition
# 
# head(combined_sig)
# sex_genes=c("Ddx3y","Eif2s3y","Kdm5d","Uty","Xist")
# combined_sig_nona = combined_sig %>%
#   dplyr::filter(!(gene_name %in% sex_genes)) 
# 
# head(combined_sig)
# #combined_sig_nona = combined_sig 
# 
# unstim = combined_sig_nona %>%
#   dplyr::filter(type == "ko_un_b6_un") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# lps = combined_sig_nona %>%
#   dplyr::filter(type == "ko_lps_b6_lps") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# lps5z7 = combined_sig_nona %>%
#   dplyr::filter(type == "ko_lps5z7_b6_lps5z7") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# ko_1 = combined_sig_nona %>%
#   dplyr::filter(type == "ko_lps5z7_ko_lps") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# b6_1 = combined_sig_nona %>%
#   dplyr::filter(type == "b6_lps5z7_b6_lps") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# 
# ko_2 = combined_sig_nona %>%
#   dplyr::filter(type == "ko_lps_ko_un") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# b6_2 = combined_sig_nona %>%
#   dplyr::filter(type == "b6_lps_b6_un") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# ko_3 = combined_sig_nona %>%
#   dplyr::filter(type == "ko_lps5z7_ko_un") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# b6_3 = combined_sig_nona %>%
#   dplyr::filter(type == "b6_lps5z7_b6_un") %>%
#   mutate(abs_lfc = abs(log2FoldChange))%>%
#   slice_max(order_by = abs_lfc, n = 20)
# 
# 
# normalized_counts_conv = normalized_counts %>%
#   rownames_to_column("gene_id") %>%
#   left_join(gtf,by="gene_id")
# 
# norm_sig <- normalized_counts_conv %>%
#   dplyr::filter(gene_name %in% b6_3$gene_name) %>%
#   column_to_rownames("gene_name")
# 
# norm_sig$gene_type = NULL
# norm_sig$gene_id = NULL
# 
# head(norm_sig)
# 
# annotation <- meta %>%
#   dplyr::select(sample_name, geno,cond) %>%
#   data.frame(row.names = "sample_name")
# 
# heat_colors <- brewer.pal(6, "YlOrRd")
# 
# dpi=300
# #png("/cluster/tufts/poltoraklab/rbator01/trif_june2021//plots/heatmap.png",width = dpi*15, height = dpi*10, units = "px",res = dpi,type='cairo')
# pheatmap(norm_sig ,
#          color = heat_colors,
#          cluster_rows = T,
#          show_rownames = T,
#          annotation = annotation,
#          border_color = NA,
#          fontsize = 10,
#          scale = "row",
#          fontsize_row = 10,
#          height = 20, 
#          main="B6:LPS5z7vsUnstim")
# 
# #dev.off()
# 
# ## clusterprofiler
# combined_sig$type = factor(combined_sig$type, levels=c("ko_un_b6_un", "ko_lps_b6_lps","ko_lps5z7_b6_lps5z7","ko_lps5z7_ko_lps","b6_lps5z7_b6_lps", 
#                                                        "ko_lps_ko_un","b6_lps_b6_un","ko_lps5z7_ko_un","b6_lps5z7_b6_un",
#                                                        "genoT169A.condLPS5z7_LPS", "genoT169A.condLPS_un","genoT169A.condLPS5z7_un"))
# table(combined_sig$type)
# 
# combined_sig$ensembl = gsub('\\..+$', '', combined_sig$gene_id)
# 
# ck <- compareCluster(geneCluster = ensembl~type,
#                      data = combined_sig,
#                      OrgDb = org.Mm.eg.db,
#                      keyType="ENSEMBL",
#                      fun = "enrichGO",
#                      ont="BP",
#                      readable=T)
# #keytypes(org.Mm.eg.db)
# 
# saveRDS(ck,"/cluster/tufts/poltoraklab/rbator01/ripk1_may2023/analysis/deseq2/ck_pairwise_and_int_padj_0.05_26may23.rds")
# write.xlsx(ck,"/cluster/tufts/poltoraklab/rbator01/ripk1_may2023/analysis/deseq2/ck_pairwise_and_int_padj_0.05_26may23.xlsx")
# 
# 
# p = clusterProfiler::dotplot(ck, x =~Cluster,showCategory =6, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# 
# # include all treatment levels
# # this seems to get the same genes as one at a time, and that is clearer to interpret?
# # 
# # dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ geno + cond + cond:geno)
# # dds$geno = factor(dds$geno, levels=c("B6","T169A"))
# # dds$cond = factor(dds$cond, levels=c("unstim","LPS","LPS5z7"))
# # 
# # dds <- DESeq(dds)
# # rld <- rlog(dds, blind=TRUE)
# # 
# # resultsNames(dds)
# # int_lps_unstim <- results(dds, name="genoT169A.condLPS")
# # 
# # int_lps_unstim_s <- lfcShrink(dds, 
# #                               coef="genoT169A.condLPS", 
# #                               res=int_lps_unstim)
# # 
# # int_lps_unstim_s_res = data.frame(int_lps_unstim_s) %>%
# #   mutate(type = "genoT169A.condLPS_unstim") %>%
# #   rownames_to_column("gene_id") 
# # 
# # 
# # 
# # 
# # int_lps5z7_unstim <- results(dds, name="genoT169A.condLPS5z7")
# # 
# # int_lps5z7_unstim_s <- lfcShrink(dds, 
# #                                  coef="genoT169A.condLPS5z7", 
# #                                  res=int_lps5z7_unstim)
# # 
# # int_lps5z7_unstim_s_res = data.frame(int_lps5z7_unstim_s) %>%
# #   mutate(type = "genoT169A.condLPS5z7_unstim") %>%
# #   rownames_to_column("gene_id") 
# # 
# # combined = combined %>% 
# #   left_join(gtf, by="gene_id")
# # 
# # combined_sig=data.frame(combined) %>% 
# #   dplyr::filter(padj<0.05) 
# # 
# # int_lps5z7_lps <- results(dds, name="genoT169A.condLPS5z7_LPS")
# # 
# 
# 
# 
# # lots of plots ----
# 
# plot_gene <- function(name, gtf, dds){
#   
#   get_ens = gtf %>% dplyr::filter(gene_name == name)
#   ens=get_ens$gene_id
#   
#   data = plotCounts(dds, ens, intgroup = "cond_geno", returnData = T)
#   
#   p = ggplot(data) +
#     geom_point(aes(x=cond_geno, y=count)) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#     ggtitle(name) +
#     xlab("") +
#     scale_y_continuous(trans='log10')
#   
#   ggsave(p, file = paste0("analysis/plots/",name, "_counts.pdf"),
#          width = 3,
#          height = 3,
#          units = "in")
# }
# 
# # for (gene in sasha_gene){
# #   plot_gene(gene, gtf, dds)
# # }
# 
# sex_genes=c("Ddx3y","Eif2s3y","Kdm5d","Uty","Xist")
# 
# for (gene in sex_genes){
#   plot_gene(gene, gtf, dds)
# }
# 
# plot_gene("Ido2", gtf, dds)
# plot_gene("Nrg1", gtf, dds)
# 
# name="Atp6v0c-ps2"
# get_ens = gtf %>% dplyr::filter(gene_name == name)
# get_ens
# ens="ENSMUSG00000080242.6"
# 
# data = plotCounts(dds, ens, intgroup = "cond_geno", returnData = T)
# 
# 
