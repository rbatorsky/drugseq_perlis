LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
select <- dplyr::select
library(openxlsx)
library(DESeq2)
library('org.Hs.eg.db')
library(clusterProfiler)
library(data.table)
library(msigdbr)
library(ggfortify)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# scr2  ----
scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format_add_rna_activity.xlsx") 

scr2_50 = scr2 %>%
  filter(PI_below_50percent_calc == 1)

nrow(scr2_50) # 28 are PI <50% DMSO

scr2_50_active = scr2_50 %>%
  filter(active == 1)

nrow(scr2_50_active) # 16 are active

# make a bar plot of the ndeg per compound -----

scr2_50_deg = scr2_50 %>%
  select(c(Compound_Name, ndeg, active)) %>%
  mutate(ndeg = ifelse(is.na(ndeg),0,ndeg)) %>%
  mutate(rna_active = factor(ifelse(active == 1,1,0), levels=c(0,1))) %>%
  mutate(Compound_Name = gsub("_10uM","", Compound_Name)) %>%
  mutate(Compound_Name = gsub("_2uM","", Compound_Name)) %>%
  mutate(Compound_Name = gsub("_hydrochloride","_HCL",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_Hydrochloride","_HCL",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_acid","_Acid",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_ditosylate","_Ditosylate",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_tosylate","_Tosylate",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_","",Compound_Name))%>%
  rename(`RNA Active` = rna_active)


scr2_50_deg$Compound_Name

p = ggplot(data=scr2_50_deg, aes(x=reorder(Compound_Name, ndeg), y=ndeg, fill=`RNA Active`)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Number of DEG vs. DMSO") +
  xlab("")+ coord_flip() + theme_classic()+ 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))+
  scale_fill_manual(values=c("grey", "black"))

ggsave(p, filename = "rd1_rd2_analysis/paper_plots/5_a_deg_barplot.pdf",
       device = cairo_pdf,
       width = 5, height = 7,
       units = "in")

show(p)

# What if the threshold was p<0.05

deg = fread("rd1_rd2_analysis/de/all_res_run3.csv") %>%
  mutate(Compound_Name = gsub("_10uM", "", drug)) %>%
  mutate(Compound_Name = gsub("_2uM", "", Compound_Name)) %>%
  mutate(Compound_Name = gsub("_run3.csv", "", Compound_Name)) 

deg_stringent = deg %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

table(deg_stringent$drug)


# analyze AIF1 expression ----
deg = fread("rd1_rd2_analysis/de/all_res_run3.csv") %>%
  mutate(Compound_Name = gsub("_10uM", "", drug)) %>%
  mutate(Compound_Name = gsub("_2uM", "", Compound_Name)) %>%
  mutate(Compound_Name = gsub("_run3.csv", "", Compound_Name)) 


# aif1 correlation
aif1_deg = deg %>%
  filter(gene == "ENSG00000204472|AIF1") %>%
  select(log2FoldChange, Compound_Name) %>%
  inner_join(scr2 %>%
              filter(active == 1) %>%
               select(Compound_Name, well_mean_normalized_IBA1_intensity), by="Compound_Name")

library(ggpubr)
ggscatter(aif1_deg, x = "log2FoldChange", y = "well_mean_normalized_IBA1_intensity",
          add = "reg.line",                               
          conf.int = TRUE,
          label = "Compound_Name", repel = TRUE,
          # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = -2, label.y = .06)  + # Add correlation coefficient
  xlab("AIF1 Log2FoldChange vs. DMSO")


# select these genes from all gene expression matrix
deg_select = deg %>%
  filter(Compound_Name %in% unique(scr2_50$Compound_Name)) %>%
  filter(gene  == "ENSG00000204472|AIF1") %>%
  select("gene","log2FoldChange","padj", "Compound_Name")

lfc_select = deg_select %>%
  dplyr::select("gene","log2FoldChange", "Compound_Name") %>%
  spread("Compound_Name","log2FoldChange", fill=0)


rownames(lfc_select) = lfc_select$gene
lfc_select$gene = NULL

p_select = deg_select %>%
  dplyr::select("gene","padj", "Compound_Name") %>%
  spread("Compound_Name","padj", fill=1)

rownames(p_select) = p_select$gene
p_select$gene = NULL

heatmap_max = max(lfc_select)
heatmap_min = min(lfc_select)

library(circlize)
library(ComplexHeatmap)

annot = scr2_50 %>%
  select(Compound_Name, mean_normalized_PI, well_mean_eccentricity, well_mean_solidity,well_mean_normalized_IBA1_intensity) %>%
  distinct() %>%
  column_to_rownames("Compound_Name") %>%
  as.data.frame() 

lfc_select = lfc_select[,rownames(annot)]
p_select = p_select[,rownames(annot)]

p_select_mat = as.matrix(p_select)
lfc_select_mat = as.matrix(lfc_select)

ph = Heatmap(lfc_select_mat,
             col = circlize::colorRamp2(c(-10, 0, 10), c("Darkblue", "white", "red")),
             cluster_columns=T,
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(p_select_mat[i, j] < 0.1 ) {
                 grid.text("*", x, y)
               }
             },
             heatmap_legend_param = list(title = "Scaled log\nfold change"),
             column_names_gp = grid::gpar(fontsize = 8),
             row_names_gp = grid::gpar(fontsize = 8) ,
             column_title = "AIF1",
             top_annotation = HeatmapAnnotation(df = annot,
                                                annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
                                                annotation_name_gp= gpar(fontsize = 8)))


# make a pca of the fold changes -----
deg_filter = deg %>%
  filter(Compound_Name %in% unique(scr2_50_active$Compound_Name)) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 1) %>%
  group_by(Compound_Name) %>%
  slice_max(n=200, order_by=abs(log2FoldChange))

head(deg_filter)
table(deg_filter$Compound_Name)

deg_filter_wide = deg %>%
  filter(Compound_Name %in% unique(scr2_50_active$Compound_Name)) %>%
  filter(gene %in% unique(deg_filter$gene)) %>%
  select(gene, Compound_Name, log2FoldChange) %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange),0, log2FoldChange)) %>%
  pivot_wider(id_cols="Compound_Name", names_from="gene", values_from = "log2FoldChange", values_fill = 0) 


deg_filter_wide_meta = deg_filter_wide %>%
  left_join(scr2_50_active, by="Compound_Name") %>%
  column_to_rownames("Compound_Name")

deg_filter_wide  = deg_filter_wide  %>%
  column_to_rownames("Compound_Name")

# PCA, really dominated by strong compounds -----
pca_res <- prcomp(deg_filter_wide)

autoplot(pca_res,
         data = deg_filter_wide_meta, label=T,label.size = 3)

pca_res.weights <- data.frame(pca_res$rotation)

pc1 <- pca_res.weights[order(pca_res.weights[, 1], decreasing = TRUE),][, 1, drop = FALSE]
pc2 <- pca_res.weights[order(pca_res.weights[, 2], decreasing = TRUE),][, 2, drop = FALSE]

pc1$pc = "pc1"
pc1$weight = pc1$PC1
pc1$PC1=NULL
pc2$pc = "pc2"
pc2$weight = pc2$PC2
pc2$PC2=NULL
pc = rbind(pc1, pc2)

head(pc)
pc_genes = pc %>%
  rownames_to_column("feature") %>%
  separate(feature, into=c("ens","sym"),sep="\\|") %>%
  group_by(pc) %>%
  slice_max(weight, n=200)

# GO Enrichment
ck<- compareCluster(geneCluster = ens ~ pc,
                    data = pc_genes,
                    OrgDb = org.Hs.eg.db,
                    keyType="ENSEMBL",
                    fun = "enrichGO",
                    ont="BP",
                    readable=T)

dotplot(ck)
view(ck)
saveRDS(ck, "rd1_rd2_analysis/de/pc1_pc2_gobp_top200.rds")


# umap
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
library(umap)
umap_results <- umap::umap(deg_filter_wide)

head(scr2_50)
colnames(data.frame(umap_results$layout))

umap_plot_df <- data.frame(umap_results$layout) %>%
  rownames_to_column("Compound_Name") %>%
  inner_join(scr2_50 , by = "Compound_Name") %>%
  dplyr::rename(UMAP1 = X1) %>%
  dplyr::rename(UMAP2 = X2)


plot="well_mean_area"                     
plot="well_mean_eccentricity"
plot="well_mean_solidity"                 
plot="well_mean_normalized_IBA1_intensity"
plot="mean_normalized_PI"    

ggplot(
  umap_plot_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    label= Compound_Name,
    color = mean_normalized_PI
  )
) +
  geom_point(size = 3, position = position_dodge(width = 1)) +
  geom_text(position = position_dodge(width = 1), 
            vjust = 0.4, 
            hjust = -0.2, 
            size = 2.5, 
            show.legend = FALSE, 
            color='black')


# filter -----
deg_filter = deg %>%
  filter(!gene == "") %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 1) %>%
  filter(Compound_Name %in% unique(scr2_50$treatment))

# pick which genes to plot
deg_to_plot = deg %>%
  filter(Compound_Name %in% unique(scr2_50_active$Compound_Name)) %>%
  filter(!gene == "") %>%
  filter(padj < 0.001) %>%
  slice_max(order_by = abs(log2FoldChange), n = 3, by = Compound_Name)

# select these genes from all gene expression matrix
deg_select = deg %>%
  filter(Compound_Name %in% unique(scr2_50_active$Compound_Name)) %>%
  filter(gene %in% unique(deg_to_plot$gene)) %>%
  arrange(factor(Compound_Name, levels=unique(deg$Compound_Name))) %>%
  select("gene","log2FoldChange","padj", "Compound_Name")

lfc_select = deg_select %>%
  dplyr::select("gene","log2FoldChange", "Compound_Name") %>%
  spread("Compound_Name","log2FoldChange", fill=0)

rownames(lfc_select) = lfc_select$gene
lfc_select$gene = NULL
lfc_select_mat = as.matrix(lfc_select)

p_select = deg_select %>%
  dplyr::select("gene","padj", "Compound_Name") %>%
  spread("Compound_Name","padj", fill=1)

rownames(p_select) = p_select$gene
p_select$gene = NULL
p_select_mat = as.matrix(p_select)

heatmap_max = max(lfc_select)
heatmap_min = min(lfc_select)
heatmap_max
library(circlize)
library(ComplexHeatmap)

annot = scr2_50_active %>%
  select(Compound_Name, mean_normalized_PI, well_mean_eccentricity, well_mean_solidity) %>%
  distinct() %>%
  column_to_rownames("Compound_Name") %>%
  as.data.frame() 

dim(lfc_select_mat)
lfc_select_mat = lfc_select_mat[,rownames(annot)]
p_select_mat = p_select_mat[,rownames(annot)]

ph = Heatmap(lfc_select_mat,
             col = circlize::colorRamp2(c(-10, 0, 10), c("Darkblue", "white", "red")),
             cluster_columns=T,
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(p_select_mat[i, j] < 0.1 ) {
                 grid.text("*", x, y)
               }
             },
             heatmap_legend_param = list(title = "Scaled log\nfold change"),
             column_names_gp = grid::gpar(fontsize = 8),
             row_names_gp = grid::gpar(fontsize = 8) ,
             column_title = "top 3 DEG per treatment by abs(log2FoldChange)",
             top_annotation = HeatmapAnnotation(df = annot,
                                                annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
                                                annotation_name_gp= gpar(fontsize = 8)))


pdf("rd1_rd2_analysis/paper_plots/top3_deg_heatmap.pdf", height=16, width=7)
print(ph)
dev.off()

