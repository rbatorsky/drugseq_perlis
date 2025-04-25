LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(openxlsx)
library(ggforce)
library('org.Hs.eg.db')
select = dplyr::select
rename = dplyr::rename

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# export counts for liam 1/15 ----
mat = readRDS("rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.rds")
write.xlsx(mat,"rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.xlsx" )

dds = readRDS("rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.dds.rds")
counts = counts(dds, normalized=FALSE)
write.xlsx(counts,"rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.xlsx" )
norm_counts = counts(dds, normalized=TRUE)
write.xlsx(norm_counts,"rd1_rd2_analysis/rd1_rd2_umi.nondedup.norm_counts.xlsx" )

# read in meta data ----
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name
meta = meta %>%
  dplyr::filter(!(Drug_Treatment == "Empty well")) %>%
  mutate(drug = paste0(Drug_Treatment, "_", Concentration))

# add the rna activity, filter for <50% and RNA active ---
activity = read.xlsx("rd1_rd2_analysis/de/all_compounds_rnaactivity_functionalscreen_23aug24.xlsx") %>%
  select(Compound_Name, ndeg) %>%
  mutate(ndeg = ifelse(is.na(ndeg),0,ndeg)) %>%
  mutate(rna_active = factor(ifelse(ndeg > 193,1,0), levels=c(0,1))) %>%
  filter(rna_active == 1)

nrow(activity)

scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format_add_rna_activity.xlsx")

scr2_50 = scr2 %>%
 filter(PI_below_50percent_calc == 1 & active == 1)

meta = meta %>%
  filter(Concentration != '5uM') %>%
  inner_join(scr2_50, by=c("Drug_Treatment" = "Compound_Name"))

# read in count matrix -----
mat = readRDS("rd1_rd2_analysis/rd1_rd2_umi.notrim.nondedup.counts.rds")
mat = mat[,rownames(meta)]
rownames(meta)
all(colnames(mat) %in% rownames(meta))
mat = mat[,rownames(meta)]
all(colnames(mat) == rownames(meta))
ncol(mat)

dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = ~ batch + Drug_Treatment)
dds <- DESeq(dds)
saveRDS(dds, "rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.notrim.dds.2scr_50_active.rds")
#
# dimension reduction ----
dds = readRDS("rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.dds.2scr_50_active.rds")

# reads per sample
counts = dds@assays@data$counts

cs= colSums(counts)
# hist(cs)
# median(cs)
# mean(cs)
# sd(cs)

vst <- vst(dds, blind=FALSE)
vst_mat <- assay(vst)
vst_batch <- limma::removeBatchEffect(assay(vst), vst$batch)

# dds_batch_correct = dds
# vst_batch_correct = vst
# assay(vst_batch_correct) = vst_batch
#
#
# ## Compute pairwise correlation values ----

colnames = colnames(vst_batch)

colnames =  str_replace(colnames, "^(.*?)\\|(.*?)\\|rd(\\d+)$", "\\2|\\1|\\3") |> 
  str_replace_all("_hydrochloride","_HCL") |> 
  str_replace_all("_Hydrochloride","_HCL")|> 
  str_replace_all("_acid","_Acid")|> 
  str_replace_all("_ditosylate","_Ditosylate")|> 
  str_replace_all("_tosylate","_Tosylate")|> 
  str_replace_all("\\|","_")

colnames(vst_batch) = colnames

head(vst_batch)
vst_cor <- cor(vst_batch)

library(ComplexHeatmap)
colnames(meta)
meta_heatmap = meta %>%
  select(sample_name, Drug_Treatment) %>%
  mutate(sample_name = str_replace(sample_name, "^(.*?)\\|(.*?)\\|rd(\\d+)$", "\\2|\\1|\\3")) %>%
  mutate(sample_name = gsub("_hydrochloride","_HCL",sample_name)) %>%
  mutate(sample_name = gsub("_Hydrochloride","_HCL",sample_name)) %>%
  mutate(sample_name = gsub("_acid","_Acid",sample_name)) %>%
  mutate(sample_name = gsub("_ditosylate","_Ditosylate",sample_name)) %>%
  mutate(sample_name = gsub("_tosylate","_Tosylate",sample_name)) %>%
  mutate(sample_name = gsub("_","",sample_name)) %>%
  column_to_rownames("sample_name")

# some are duplicate, that's why the batch name is neede

head(meta_heatmap)

dim(vst_cor)
dim(meta_heatmap)

top_ha = HeatmapAnnotation(df = meta_heatmap)

colnames(vst_cor)
head(meta_heatmap)

pdf("rd1_rd2_analysis/paper_plots/corr_heatmap.pdf",height=8, width=10)
p = Heatmap(vst_cor,
        top_annotation = top_ha,
        column_names_side = "top",                # Put column names at the top
        column_dend_side = "top",  
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = 6))
print(p)
dev.off()


# PCA ----
# pcaData <- plotPCA(vst_batch_correct, intgroup=c("batch","Drug_Treatment","Plate_Well"), returnData=TRUE)
# 
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# 
# pcaData = pcaData %>%
#   separate(name, into=c("Plate_Well","Drug_Treatment","batch"), remove=F, sep="\\|")
# 
# head(pcaData)
# ggplot(pcaData, aes(PC1, PC2, color=Drug_Treatment, label=name)) +
#   geom_point(size=1) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed() +
#   geom_text(hjust=0, vjust=0, size=2)

# UMAP -----
# # https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
library(umap)

custom.config <- umap.defaults
custom.config$random_state <- 123
umap_results <- umap::umap(t(vst_batch))

umap_plot_df <- data.frame(umap_results$layout) %>%
  rownames_to_column("sample_name") %>%
  inner_join(meta, by = "sample_name") %>%
  dplyr::rename(UMAP1 = X1) %>%
  dplyr::rename(UMAP2 = X2) %>%
  mutate(drug = gsub("_10uM","", drug)) %>%
  mutate(drug = gsub("_2uM","", drug)) %>%
  mutate(drug = gsub("_hydrochloride","_HCL",drug)) %>%
  mutate(drug = gsub("_Hydrochloride","_HCL",drug)) %>%
  mutate(drug = gsub("_acid","_Acid",drug)) %>%
  mutate(drug = gsub("_ditosylate","_Ditosylate",drug)) %>%
  mutate(drug = gsub("_tosylate","_Tosylate",drug)) %>%
  mutate(drug = gsub("_","",drug))


# Plot using `ggplot()` function
p = ggplot(
  umap_plot_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    color=drug
  )
) +
  geom_point(size = 3) +
  theme_bw() +
  theme(legend.position="none") +
  ggforce::geom_mark_ellipse(aes(color = drug, label=drug))
ggsave(p, filename="rd1_rd2_analysis/paper_plots/4_b_umap.pdf", height=10, width=10)

# ggplot(
#   umap_plot_df,
#   aes(
#     x = UMAP1,
#     y = UMAP2,
#     label= drug,
#     color=drug
#   )
# ) +
#   geom_point(size = 3) +
#   geom_text_repel(col="black", size=3)  +
#   theme_bw(base_size = 12) +
#   guides(col= guide_legend(title= "Treatment"))+
#   theme(legend.position="none")+
#   stat_ellipse(level = 0.9)

# which are the active ones

dont_group = c("Lazertinib","Lapatinib","LapatinibDitosylate","SRX246","AlectinibHCL","SRT2104","BF227")
dont_group_df = umap_plot_df %>%
  filter(drug %in% dont_group) %>%
  select(drug, Drug_Treatment)

activity = read.xlsx("rd1_rd2_analysis/de/all_compounds_rnaactivity_functionalscreen_23aug24.xlsx") %>%
  select(Compound_Name, ndeg) %>%
  filter(Compound_Name %in% dont_group_df$Drug_Treatment)

activity

select = umap_plot_df %>%
  filter(!(drug %in% dont_group))

unique(select$drug)

ggplot(
  select,
  aes(
    x = UMAP1,
    y = UMAP2,
    color=drug
  )
) +
  geom_point(size = 3) +
  theme_bw() +
  theme(legend.position="none") +
  ggforce::geom_mark_ellipse(aes(color = drug, label=drug),
                             label.fontsize=10,
                             con.type="straight",
                             con.size = 0.5)


# # highly variable genes -----
# 
# var_genes <- apply(vst_batch, 1, var)
# select_var <- names(sort(var_genes, decreasing=TRUE))[1:75]
# vst_hv <- vst_mat[select_var,]
# vst_hv
# 
# rownames(vst_hv) = gsub(".*\\|","",rownames(vst_hv))
# 
# meta_heatmap = meta %>%
#   dplyr::select(batch, Drug_Treatment) 
# 
# bottom_ha = HeatmapAnnotation(df = meta_heatmap)
# 
# vst_hv_scaled = t(scale(t(vst_hv)))
# 
# vst_hv_scaled = vst_hv_scaled[,rownames(meta_heatmap)]
# Heatmap(vst_hv_scaled,
#         bottom_annotation = bottom_ha,
#         column_names_gp = grid::gpar(fontsize = 6),
#         row_names_gp = grid::gpar(fontsize = 6))
# 
# head(deg)
# 
# # Median over the replicates ----
# meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
# rownames(meta) = meta$sample_name
# meta = meta %>%
#   dplyr::filter(!(Drug_Treatment == "Empty well")) %>%
#   mutate(drug = paste0(Drug_Treatment, "_", Concentration))
# 
# scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format_add_rna_activity.xlsx") 
# 
# scr2_50 = scr2 %>%
#   filter(PI_below_50percent_calc == 1) 
# 
# meta = meta %>%
#   filter(Concentration != '5uM') %>%
#   inner_join(scr2_50, by=c("Drug_Treatment" = "Compound_Name"))
# 
# phenos=c('mean_normalized_PI','well_mean_eccentricity','well_mean_solidity',"well_mean_normalized_IBA1_intensity")
# 
# vst_batch_long = vst_batch %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(names_to="sample_name", cols = colnames(vst_batch)) %>%
#   left_join(meta %>% as.data.frame() %>% select(sample_name, Drug_Treatment), by="sample_name") %>%
#   group_by(gene, Drug_Treatment) %>%
#   #summarise(med = median(value, na.rm = T))%>%
#   summarise(mean_AIF1_expression = mean(value, na.rm = T))
# 
# write.xlsx(vst_batch_long, "rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.dds.2scr_50_active.rds")
# head(vst_batch_long)
# # AIF1 corraltion ------
# 
# aif1 = vst_batch_long %>%
#   filter(gene == "ENSG00000204472|AIF1") %>%
#   left_join(scr2 %>%
#               select(Compound_Name, well_mean_normalized_IBA1_intensity), by=c("Drug_Treatment" = "Compound_Name"))
# 
# head(aif1)
# 
# ggscatter(aif1, x = "mean_AIF1_expression", y = "well_mean_normalized_IBA1_intensity",
#           add = "reg.line",                               
#           conf.int = TRUE,
#           label = "Drug_Treatment", repel = TRUE,
#           # Add confidence interval
#           add.params = list(color = "blue",
#                             fill = "lightgray")
# )+
#   stat_cor(method = "pearson", label.x = 10, label.y = .05)  # Add correlation coefficient
# 
# 
# # Not using anything below ----
# vst_batch_med = vst_batch_long %>%
#   pivot_wider(names_from = "Drug_Treatment", values_from="med", id_cols="gene") %>%
#   column_to_rownames("gene") %>%
#   as.matrix()
# 
# vst_batch_cor <- cor(vst_batch_med)   
# 
# library(ComplexHeatmap)
# 
# head(vst_batch_cor)
# meta_heatmap = meta %>%
#   select(Drug_Treatment, phenos) %>%
#   distinct() %>%
#   column_to_rownames("Drug_Treatment") 
# 
# meta_heatmap = meta_heatmap[rownames(vst_batch_cor),]
# colnames(meta_heatmap )
# library(circlize)
# 
# # bottom_ha = HeatmapAnnotation(df = meta_heatmap, col = list(MOA = c(
# #   "Cell_Cycle_Proliferation_Inhibitors" = "red",
# #   "Kinase_Inhibitors" = "blue",
# #   "DNA_Damage_Apoptosis_Regulators" = "black",
# #   "Immunomodulation_Inflammation" = "cyan",
# #   "Signal_Transduction_Modulators"  = "pink",
# #   "Neuroprotection_Neurogenesis" ="green",
# #   "Epigenetic_Modulators" = "purple"
# # )))
# 
# col_fun = colorRamp2(c(0.90, 0.95, 1), c("blue", "yellow", "red"))
# 
# Heatmap(vst_batch_cor,
#         column_names_gp = grid::gpar(fontsize = 10),
#         row_names_gp = grid::gpar(fontsize = 10),
#         name = "sample-sample correlation",
#         col = col_fun)
# 
# # no using
# # library(ggfortify)
# # 
# # vst_batch_med_t = t(vst_batch_med) %>%
# #   as.data.frame() %>%
# #   rownames_to_column("Drug_Treatment") %>%
# #   left_join(meta_heatmap %>%
# #               rownames_to_column("Drug_Treatment"), by="Drug_Treatment")
# # 
# # pca_res <- prcomp(t(vst_batch_med))
# # 
# # library(ggrepel)
# # 
# # dev.off()
# # autoplot(pca_res, 
# #          data = vst_batch_med_t, 
# #          color = "MOA", 
# #          label=T,
# #          label.size = 5) + xlim(-1,0.4)
# # 
# # 
# # autoplot(pca_res,    
# #          x = 3,    
# #          y = 4, 
# #          data = vst_batch_med_t, 
# #          color = "MOA", 
# #          label=T,
# #          label.size = 5) + xlim(-0.75,0.75)
# # 
# # autoplot(pca_res,    
# #          x = 5,    
# #          y = 6, 
# #          data = vst_batch_med_t, 
# #          color = "MOA", 
# #          label=T,
# #          label.size = 5) 
# # 
# # pca_res.weights <- data.frame(pca_res$rotation)
# # 
# # pc1 <- pca_res.weights[order(pca_res.weights[, 1], decreasing = TRUE),][, 1, drop = FALSE]
# # pc2 <- pca_res.weights[order(pca_res.weights[, 2], decreasing = TRUE),][, 2, drop = FALSE]
# # 
# # pc1_top = data.frame(gene = rownames(head(pc1,100)))
# # pc2_top = data.frame(gene = rownames(head(pc2,100)))
# # 
# # pc1_top$type = "pc1"
# # pc2_top$type = "pc2"
# # 
# # pc1_top
# # pc2_top
# # 
# # pc = rbind(pc1_top, pc2_top) %>%
# #   separate(gene, into=c("ens", "sym"), sep="\\|", remove=T)
# # 
# # ck<- compareCluster(geneCluster = ens~type,
# #                     data = pc,
# #                     OrgDb = org.Hs.eg.db,
# #                     keyType="ENSEMBL",
# #                     fun = "enrichGO",
# #                     ont="BP", 
# #                     readable=T)
# # 
# # dotplot(ck,  show=2, by="Count")+ 
# #   theme(text = element_text(size = 10),
# #         axis.text.x = element_text(angle = 90, size=10, hjust = 1),
# #         axis.text.y = element_text(size=10)) +
# #   scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) 
# # 
# # 
# # view(ck)
# # saveRDS(ck, "rd1_rd2_analysis/dim_red/pc1_pc2_gobp.rds")
# # 
# # 
# # # umap
# # # https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
# # library(umap)
# # umap_results <- umap::umap(t(vst_batch_med))
# # 
# # colnames(data.frame(umap_results$layout))
# # umap_plot_df <- data.frame(umap_results$layout) %>%
# #   rownames_to_column("Drug_Treatment") %>%
# #   inner_join(meta_heatmap %>%
# #                rownames_to_column("Drug_Treatment"), by = "Drug_Treatment") %>%
# #   dplyr::rename(UMAP1 = X1) %>%
# #   dplyr::rename(UMAP2 = X2)
# # 
# # # Plot using `ggplot()` function
# # ggplot(
# #   umap_plot_df,
# #   aes(
# #     x = UMAP1,
# #     y = UMAP2,
# #     label= Drug_Treatment,
# #     color=MOA
# #   )
# # ) +
# #   geom_point(size = 3, position = position_dodge(width = 1)) +
# #   geom_text(position = position_dodge(width = 1), 
# #             vjust = 0.4, 
# #             hjust = -0.1, 
# #             size = 2.5, 
# #             show.legend = FALSE, 
# #             color='black')
# 
# 
# # highly variable genes -----
# var_genes <- apply(vst_batch_med, 1, var)
# select_var <- names(sort(var_genes, decreasing=TRUE))[1:50]
# vst_hv <- vst_batch_med[select_var,]
# vst_hv
# 
# rownames(vst_hv) = gsub(".*\\|","",rownames(vst_hv))
# 
# head(meta)
# meta_heatmap = meta %>%
#   dplyr::select(Drug_Treatment, phenos) %>%
#   distinct() %>%
#   column_to_rownames("Drug_Treatment")
# 
# # bottom_ha = HeatmapAnnotation(df = meta_heatmap, col = list(MOA = c(
# #   "Cell_Cycle_Proliferation_Inhibitors" = "red",
# #   "Kinase_Inhibitors" = "blue",
# #   "DNA_Damage_Apoptosis_Regulators" = "black",
# #   "Immunomodulation_Inflammation" = "cyan",
# #   "Signal_Transduction_Modulators"  = "pink",
# #   "Neuroprotection_Neurogenesis" ="green",
# #   "Epigenetic_Modulators" = "purple"
# # )))
# 
# vst_hv_scaled = t(scale(t(vst_hv)))
# 
# Heatmap(vst_hv_scaled,
#         column_names_gp = grid::gpar(fontsize = 10),
#         row_names_gp = grid::gpar(fontsize = 8))
# 
# 
# # analyze dmso var genes ------
# dmso_meta = meta %>% 
#   filter(Drug_Treatment=="DMSO") 
# 
# vst_batch_dmso = vst_batch[,dmso_meta$sample_name]
# 
# var_genes <- apply(vst_batch_dmso, 1, var)
# 
# select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]
# 
# vst_hv <- vst_batch_dmso[select_var,]
# 
# rownames(vst_hv) = gsub(".*\\|","",rownames(vst_hv))
# 
# meta_heatmap = dmso_meta %>%
#   dplyr::select(batch, Drug_Treatment)
# 
# bottom_ha = HeatmapAnnotation(df = meta_heatmap)
# 
# vst_hv_scaled = t(scale(t(vst_hv)))
# 
# Heatmap(vst_hv_scaled,
#         bottom_annotation = bottom_ha,
#         column_names_gp = grid::gpar(fontsize = 6),
#         row_names_gp = grid::gpar(fontsize = 6))
# 
# # dds object
# dmso_meta = meta %>% 
#   filter(Drug_Treatment=="DMSO") 
# mat_dmso = mat[,dmso_meta$sample_name]
# 
# all(colnames(mat_dmso) %in% rownames(dmso_meta))
# 
# dds_dmso <- DESeqDataSetFromMatrix(countData = mat_dmso, colData = dmso_meta, design = ~ batch)
# dds_dmso <- DESeq(dds_dmso)
# saveRDS(dds_dmso, "rd1_rd2_analysis/de/rd1_rd2_dmso.nondedup.counts.dds.rds")
# vst_dmso <- vst(dds_dmso, blind=FALSE)
# vst_dmso_batch <- limma::removeBatchEffect(assay(vst_dmso), vst_dmso$batch)
# 
# pcaData <- plotPCA(vst_dmso , intgroup=c("batch","Plate_Well"), returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# 
# head(pcaData)
# ggplot(pcaData, aes(PC1, PC2, color=batch, label=Plate_Well)) +
#   geom_point(size=1) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed() + 
#   geom_text(hjust=0, vjust=0, size=2)
# 
# 
# library(umap)
# umap_results <- umap::umap(t(vst_dmso_batch))
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
#     label= Plate_Well,
#     color=batch
#   )
# ) +
#   geom_point(size = 3, position = position_dodge(width = 1)) +
#   geom_text(position = position_dodge(width = 1), 
#             vjust = 0.4, 
#             hjust = -0.1, 
#             size = 2.5, 
#             show.legend = FALSE, 
#             color='black')
# 
# 
# 
# # SWITCH TO PIPELINE 
# # back to analysis -----
# normalized_counts <- data.frame(counts(dds, normalized=TRUE))
# 
# # look at two example conditions that group well -----
# c1 = "GSK805"
# n1 = paste0("Drug_Treatment_",c1, '_vs_DMSO')
# r1 <- results(dds, name = n1)
# r1a <- lfcShrink(dds, coef = n1, type="apeglm")
# 
# c2 = "Vorinostat"
# n2 = paste0("Drug_Treatment_",c2, '_vs_DMSO')
# r2 <- results(dds, name = n2)
# r2a <- lfcShrink(dds, coef = n2, type="apeglm")
# 
# r1adf = data.frame(r1a) %>%
#   rownames_to_column("gene") 
# 
# write.xlsx(r1adf, "analysis/processed/r1adf.xlsx")
# 
# r2adf = data.frame(r2a) %>%
#   rownames_to_column("gene") 
# 
# write.xlsx(r2adf, "analysis/processed/r2adf.xlsx")
# 
# 
# r1adf_sig = r1adf %>%
#   dplyr::filter(padj< 0.05 & abs(log2FoldChange) > 1.5)%>%
#   mutate(treatment = c1)
# 
# r2adf_sig = r2adf %>%
#   dplyr::filter(padj< 0.05 & abs(log2FoldChange) > 1.5) %>%
#   mutate(treatment = c2)
# 
# results = rbind(r1adf_sig, r2adf_sig) %>%
#   separate(gene, into = c("ens","gene"), sep="\\|", remove=F)
# 
# table(results$treatment)
# head(results)
# 
# library(clusterProfiler)
# library(org.Hs.eg.db)
# 
# ck<- compareCluster(geneCluster = gene~treatment,
#                     data = results ,
#                     OrgDb = org.Hs.eg.db,
#                     keyType="SYMBOL",
#                     fun = "enrichGO",
#                     ont="BP")
# 
# saveRDS(ck, "analysis/processed/go_bp_test.rds")
# 
# dotplot(ck, show = 10,font.size = 10)
# 
# 
# 
# # all results ----
# results_names = resultsNames(dds)
# results_names = results_names[-1]
# 
# all_results = NULL
# for (r in results_names){
#   res <- results(dds, name = r)
#   res = data.frame(res) %>%
#     rownames_to_column("gene") %>%
#     mutate(name = r)
#   
#   if(is.null(all_results)){
#     all_results = res
#   }else{
#     all_results = rbind(all_results, res)
#   }
# }
# 
# write.xlsx(all_results, file = "analysis/processed/results_vs_dmso.nodedup.xlsx")
# 
# 
# all_results$name = gsub("Drug_Treatment_","",all_results$name)
# 
# r1 = all_results %>%
#   dplyr::filter(name == paste0(c1, '_vs_DMSO'))
# 
# r1_sig = r1 %>%
#   dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1.5)
# 
# nrow(r1_sig)
# 
