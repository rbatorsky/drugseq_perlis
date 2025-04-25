LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(openxlsx)

# rd3 analysis
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

mat = readRDS("analysis/processed/rd3.umi.nondedup.counts.rds")

meta = read.xlsx("Drug-Seq_Rd3_Treatment+Barcode_format.xlsx")

rownames(meta) = meta$sample_name

meta = meta %>%
  mutate(drug = paste0(Drug_Treatment, "_", `LPS.y/n`, "_", Cell.Type))


table(meta$Drug_Treatment, meta$`LPS.y/n`, meta$Treatment.Time, meta$Cell.Type)


# start deseq2 -----
all(colnames(mat) %in% rownames(meta))
mat = mat[,rownames(meta)]
all(colnames(mat) == rownames(meta))
ncol(mat)

unique(meta$drug)

# make a meta column of all the variables
dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = ~ drug)
dds$drug <- relevel(dds$drug, "VEH_DMSO_n_PBMC-iMG")

dds <- DESeq(dds)
saveRDS(dds, "rd3_analysis/de/umi.nondedup.counts.dds.rds")

# dimension reduction ----
dds = readRDS("rd3_analysis/de/umi.nondedup.counts.dds.rds")
vst <- vst(dds, blind=FALSE)
vst_mat <- assay(vst)

## Compute pairwise correlation values
vst_cor <- cor(vst_mat)   

library(ComplexHeatmap)

colnames(meta)
meta_heatmap = meta %>%
  mutate(drug = paste0(meta$Drug_Treatment, "_", meta$Concentration)) %>%
  dplyr::select(drug, "LPS.y/n", "Treatment.Time","Cell.Type")

top_ha = HeatmapAnnotation(df = meta_heatmap)

Heatmap(vst_cor,
        top_annotation = top_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))


pcaData <- plotPCA(vst, intgroup=c("LPS.y/n", "Treatment.Time","Cell.Type","Drug_Treatment","Plate_Well"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData = pcaData %>%
  mutate(celltype_time_lps = paste0(Cell.Type, "_", Treatment.Time,"_",`LPS.y.n`)) %>%
  separate(name, into=c("Plate_Well","Drug_Treatment"), remove=F, sep="\\|")

head(pcaData)
ggplot(pcaData, aes(PC1, PC2, color=celltype_time_lps, label=name)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + 
  geom_text(hjust=0, vjust=0, size=2)

# PBMC-iMG only  -----

dds_pbmc <- dds[, dds$Cell.Type %in% c("PBMC-iMG")]

vst <- vst(dds_pbmc, blind=TRUE)
vst_mat <- assay(vst)

## Compute pairwise correlation values
vst_cor <- cor(vst_mat)   

library(ComplexHeatmap)

colnames(meta)
meta_heatmap = meta %>%
  mutate(drug = paste0(meta$Drug_Treatment, "_", meta$Concentration)) %>%
  dplyr::filter(Cell.Type == "PBMC-iMG") %>%
  dplyr::select(drug, "LPS.y/n", "Treatment.Time","Cell.Type")

top_ha = HeatmapAnnotation(df = meta_heatmap)

Heatmap(vst_cor,
        top_annotation = top_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))


pcaData <- plotPCA(vst, intgroup=c("LPS.y/n", "Treatment.Time","Cell.Type","Drug_Treatment","Plate_Well"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData = pcaData %>%
  mutate(celltype_time_lps = paste0(Cell.Type, "_", Treatment.Time,"_",`LPS.y.n`)) %>%
  separate(name, into=c("Plate_Well","Drug_Treatment"), remove=F, sep="\\|")

head(pcaData)
ggplot(pcaData, aes(PC1, PC2, color=Drug_Treatment, shape=celltype_time_lps, label=name)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + 
  geom_text(hjust=0, vjust=0, size=2)

# iPSC-iMG only, no LPS  -----

dds_ipsc <- dds[, (dds$Cell.Type %in% c("iPSC-iMG") & dds$`LPS.y/n` == 'n')]

vst <- vst(dds_ipsc, blind=TRUE)
vst_mat <- assay(vst)

## Compute pairwise correlation values
vst_cor <- cor(vst_mat)   

library(ComplexHeatmap)

colnames(meta)
meta_heatmap = meta %>%
  mutate(drug = paste0(meta$Drug_Treatment, "_", meta$Concentration)) %>%
  dplyr::filter(Cell.Type == "iPSC-iMG" & `LPS.y/n` == 'n') %>%
  dplyr::select(drug, "LPS.y/n", "Treatment.Time","Cell.Type")

head(meta_heatmap)

top_ha = HeatmapAnnotation(df = meta_heatmap)

Heatmap(vst_cor,
        top_annotation = top_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))


dds_ipsc_lps <- dds[, (dds$Cell.Type %in% c("iPSC-iMG") & dds$`LPS.y/n` == 'y')]

vst <- vst(dds_ipsc_lps, blind=TRUE)
vst_mat <- assay(vst)

## Compute pairwise correlation values
vst_cor <- cor(vst_mat)   

library(ComplexHeatmap)

colnames(meta)
meta_heatmap = meta %>%
  mutate(drug = paste0(meta$Drug_Treatment, "_", meta$Concentration)) %>%
  dplyr::filter(Cell.Type == "iPSC-iMG" & `LPS.y/n` == 'y') %>%
  dplyr::select(drug, "LPS.y/n", "Treatment.Time","Cell.Type")

head(meta_heatmap)

top_ha = HeatmapAnnotation(df = meta_heatmap)

Heatmap(vst_cor,
        top_annotation = top_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))


# umap, looks the same as pca ----
# # https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
# library(umap)
# dds = readRDS("rd3_analysis/de/umi.nondedup.counts.dds.rds")
# vst <- vst(dds, blind=FALSE)
# vst_mat <- assay(vst)
# 
# ## Compute pairwise correlation values
# vst_cor <- cor(vst_mat)   
# 
# umap_results <- umap::umap(t(vst_mat))
# head(umap_results)
# 
# colnames(data.frame(umap_results$layout))
# umap_plot_df <- data.frame(umap_results$layout) %>%
#   rownames_to_column("sample_name") %>%
#   inner_join(meta, by = "sample_name") %>%
#   dplyr::rename(UMAP1 = X1) %>%
#   dplyr::rename(UMAP2 = X2)
# 
# 
# 
# # Plot using `ggplot()` function
# ggplot(
#   umap_plot_df,
#   aes(
#     x = UMAP1,
#     y = UMAP2,
#     label= drug,
#     color=`LPS.y/n`,
#     shape=Cell.Type
#   )
# ) +
#   geom_point(size = 3, position = position_dodge(width = 1)) +
#   geom_text(position = position_dodge(width = 1), 
#             vjust = 0.4, 
#             hjust = -0.1, 
#             size = 2.5, 
#             show.legend = FALSE, 
#             color='black')



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
  dplyr::select(Cell.Type, `LPS.y/n`, Drug_Treatment)

bottom_ha = HeatmapAnnotation(df = meta_heatmap)

head(vst_hv)
vst_hv_scaled = t(scale(t(vst_hv)))

Heatmap(vst_hv_scaled,
        bottom_annotation = bottom_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))

# analyze dmso var genes ------
dmso_meta = meta %>% 
  filter(Drug_Treatment=="DMSO") 

vst_batch_dmso = vst_batch[,dmso_meta$sample_name]

var_genes <- apply(vst_batch_dmso, 1, var)

select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]

vst_hv <- vst_batch_dmso[select_var,]

rownames(vst_hv) = gsub(".*\\|","",rownames(vst_hv))

meta_heatmap = dmso_meta %>%
  dplyr::select(batch, Drug_Treatment)

bottom_ha = HeatmapAnnotation(df = meta_heatmap)

vst_hv_scaled = t(scale(t(vst_hv)))

Heatmap(vst_hv_scaled,
        bottom_annotation = bottom_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))

# dds object
dmso_meta = meta %>% 
  filter(Drug_Treatment=="DMSO") 
mat_dmso = mat[,dmso_meta$sample_name]

all(colnames(mat_dmso) %in% rownames(dmso_meta))

dds_dmso <- DESeqDataSetFromMatrix(countData = mat_dmso, colData = dmso_meta, design = ~ batch)
dds_dmso <- DESeq(dds_dmso)
saveRDS(dds_dmso, "rd1_rd2_analysis/de/rd1_rd2_dmso.nondedup.counts.dds.rds")
vst_dmso <- vst(dds_dmso, blind=FALSE)
vst_dmso_batch <- limma::removeBatchEffect(assay(vst_dmso), vst_dmso$batch)

pcaData <- plotPCA(vst_dmso , intgroup=c("batch","Plate_Well"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

head(pcaData)
ggplot(pcaData, aes(PC1, PC2, color=batch, label=Plate_Well)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + 
  geom_text(hjust=0, vjust=0, size=2)


library(umap)
umap_results <- umap::umap(t(vst_dmso_batch))
head(umap_results)

colnames(data.frame(umap_results$layout))
umap_plot_df <- data.frame(umap_results$layout) %>%
  rownames_to_column("sample_name") %>%
  inner_join(meta, by = "sample_name") %>%
  dplyr::rename(UMAP1 = X1) %>%
  dplyr::rename(UMAP2 = X2)



# Plot using `ggplot()` function
ggplot(
  umap_plot_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    label= Plate_Well,
    color=batch
  )
) +
  geom_point(size = 3, position = position_dodge(width = 1)) +
  geom_text(position = position_dodge(width = 1), 
            vjust = 0.4, 
            hjust = -0.1, 
            size = 2.5, 
            show.legend = FALSE, 
            color='black')



# SWITCH TO PIPELINE 
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

