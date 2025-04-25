library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(openxlsx)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/')

# analyze rd1 only, this was done before the second batch was generated

# read in -----
mat_w_genes = readRDS("analysis/star_rd2/Solo.out/Gene/processed/umi.1mm_all_dedup.counts.noempty.rds")
meta = read.xlsx("Drug-Seq_Rd2_Treatment+Barcode.xlsx")

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

