LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(openxlsx)

# add rd2, second batch, to the rd1 analysis

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')
#setwd('~/Box/0my_projects/draft_pres/perlis_prep/drugseq_mar24/')

mat = readRDS("rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.rds")

meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name
meta = meta %>%
  dplyr::filter(!(Drug_Treatment == "Empty well")) %>%
  mutate(drug = paste0(Drug_Treatment, "_", Concentration))

head(meta)

# start deseq2 -----
all(colnames(mat) %in% rownames(meta))
mat = mat[,rownames(meta)]
all(colnames(mat) == rownames(meta))
ncol(mat)

dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = ~ batch + Drug_Treatment)
dds$Drug_Treatment <- relevel(dds$Drug_Treatment, "DMSO")

dds <- DESeq(dds)
saveRDS(dds, "rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.dds.rds")

# dimension reduction ----
dds = readRDS("rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.dds.rds")
vst <- vst(dds, blind=FALSE)
vst_mat <- assay(vst)
vst_batch <- limma::removeBatchEffect(assay(vst), vst$batch)
head(vst_batch)

## Compute pairwise correlation values
vst_cor <- cor(vst_mat)   

library(ComplexHeatmap)

meta_heatmap = meta %>%
  mutate(drug = paste0(meta$Drug_Treatment, "_", meta$Concentration)) %>%
  dplyr::select(batch, drug)

top_ha = HeatmapAnnotation(df = meta_heatmap)

Heatmap(vst_cor,
        top_annotation = top_ha,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))


pcaData <- plotPCA(vst, intgroup=c("batch","Drug_Treatment","Plate_Well"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData = pcaData %>%
  separate(name, into=c("Plate_Well","Drug_Treatment","batch"), remove=F, sep="\\|")

head(pcaData)
ggplot(pcaData, aes(PC1, PC2, color=Drug_Treatment, label=name)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + 
  geom_text(hjust=0, vjust=0, size=2)


# batch correct ----
plotPCA(vst, "batch")
plotPCA(vst_batch, "batch")

pcaData <- plotPCA(vst , intgroup=c("batch","drug","Plate_Well"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData = pcaData %>%
  separate(name, into=c("Plate_Well","drug","batch"), remove=F, sep="\\|")

head(pcaData)
ggplot(pcaData, aes(PC1, PC2, color=drug, label=name)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + 
  geom_text(hjust=0, vjust=0, size=2)

# umap
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
library(umap)
umap_results <- umap::umap(t(vst_batch))
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
    label= drug,
    color=drug,
    shape=batch
  )
) +
  geom_point(size = 3, position = position_dodge(width = 1)) +
  geom_text(position = position_dodge(width = 1), 
            vjust = 0.4, 
            hjust = -0.1, 
            size = 2.5, 
            show.legend = FALSE, 
            color='black')



# highly variable genes -----

# test
var_genes=c('ACTB, ACTG1')

rownames(vst_batch)[grepl('ACTG1', rownames(vst_batch))]


#select_var= c("ENSG00000075624|ACTB", "ENSG00000184009|ACTG1")
# test
#var_genes <- apply(vst_batch, 1, var)
#head(var_genes)
#select_var <- names(sort(var_genes, decreasing=TRUE))[1:75]
#head(select_var)
vst_hv <- vst_mat[select_var,]
dim(vst_hv)
head(vst_hv)
rownames(vst_hv) = gsub(".*\\|","",rownames(vst_hv))

meta_heatmap = meta %>%
  dplyr::select(batch, Drug_Treatment) %>%
  filter(Drug_Treatment %in% c('AEE788', 'Vorinostat','DMSO'))

bottom_ha = HeatmapAnnotation(df = meta_heatmap)

head(vst_hv)
vst_hv_scaled = t(scale(t(vst_hv)))

head(vst_hv_scaled)
vst_hv_scaled = vst_hv_scaled[,rownames(meta_heatmap)]
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

