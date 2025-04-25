LIB="/cluster/tufts/patralab/rbator01/R_libs/4.4.0/"
.libPaths(LIB)
library(openxlsx)
library(tidyverse)
select <- dplyr::select
library(clusterProfiler)
library(org.Hs.eg.db)
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# wgcna genes heatmap -----

module_genes = read.xlsx("rd1_rd2_analysis/wgcna/module_genes_2sd.xlsx") %>%
  separate(gene, into=c("ens","gene_name"), sep="\\|", remove=F)

head(module_genes)
table(module_genes$module_color)

ck = clusterProfiler::compareCluster(ens ~ module_color,
                                     data = module_genes,
                                     fun="enrichGO",
                                     OrgDb='org.Hs.eg.db',
                                     keyType = "ENSEMBL",
                                     ont="BP",
                                     readable = T)

saveRDS(ck, "rd1_rd2_analysis/wgcna/go_bp_module_genes_2sd.rds")
write.xlsx(ck,"rd1_rd2_analysis/wgcna/go_bp_module_genes_2sd.xlsx")

ck = readRDS("rd1_rd2_analysis/wgcna/go_bp_module_genes_2sd.rds")

dotplot(ck, show=5, size="Count") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 1000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


colnames(ck@compareClusterResult)
ck_corr = ck %>%
  filter(module_color %in% c("black","green", "red"))

dotplot(ck_corr, show=10, size="Count") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 1000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# plot the module genes
c="green"
module_gene_select = module_genes %>%
  dplyr::filter(module_color == c)

module_gene_select
library(pheatmap)
expr_normalized_1 = read.csv("rd1_rd2_analysis/wgcna/expr_used_for_wgcna.csv") %>%
  rename(gene = X) %>%
  filter(gene %in% module_gene_select$gene) %>%
  separate(gene, into=c("ens","gene_name"), sep="\\|", remove=T) %>%
  select(-ens) %>%
  mutate(gene_name = make.unique(gene_name)) %>%
  column_to_rownames("gene_name")

# make heatmap
pheatmap(expr_normalized_1,
         show_rownames = F,
         fontsize_col = 8,
         scale = "row", 
         cluster_rows = T,
         main = c)

# # overlap with the DEG
# 
# exo_deg = read.xlsx('../differential_abundance/mixed_effect_exo_combinedexopell_rm_mut5_6may24.xlsx') %>%
#   mutate(type="exosome") %>%
#   dplyr::select(c(accession, gene_name,Estimate , p.value, p.adjust, type, dir))%>%
#   mutate(gene_name = gsub(";.*","", gene_name))
# 
# pell_deg = read.xlsx('../differential_abundance/mixed_effect_pell_combinedexopell_rm_mut6_7may24.xlsx') %>%
#   mutate(type = "pellet") %>%
#   dplyr::select(c(accession, gene_name,Estimate , p.value, p.adjust, type, dir)) %>%
#   mutate(gene_name = gsub(";.*","", gene_name))
# 
# # subtract pellet DEG and contaminants -----
# padj_cut=0.1
# lfc_cut=0
# 
# deg_filter = pell_deg  %>%
#   dplyr::filter(p.adjust <padj_cut & abs(Estimate) > lfc_cut)
# 
# module_genes_deg = module_genes %>%
#   dplyr::filter(gene_name %in% deg_filter$gene_name)
# 
# table(module_genes_deg$colors)
# 
# ck = clusterProfiler::compareCluster(accession ~ colors,
#                                      data = module_genes_deg,
#                                      fun="enrichGO",
#                                      OrgDb='org.Hs.eg.db',
#                                      keyType = "UNIPROT",
#                                      ont="BP",
#                                      readable = T)
# 
# dotplot(ck, show=5, size="Count") +
#   scale_y_discrete(labels = function(x) str_wrap(x, width = 1000)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# saveRDS(ck, paste0(save_string,"/go_bp_module_genes_deg_28may24.rds"))
# write.xlsx(ck, paste0(save_string,"/go_bp_module_genes_deg_28may24.xlsx"))
# 
# ck = readRDS(paste0(save_string,"/go_bp_module_genes_deg_28may24.rds"))
# 
# view(ck)
# 
# # plot the module genes
# 
# unique(module_genes$colors)
# c="greenyellow"
# module_gene_select = module_genes_deg  %>%
#   dplyr::filter(colors == c)
# 
# library(pheatmap)
# 
# exo_data_ro_norm = read.xlsx("exosome/data_ro_norm_exo_combinedexopell_rm_mut5_6may24.xlsx") %>%
#   mutate(sample = paste0("exo_", sample))%>%
#   mutate(type = "exosome") %>%
#   dplyr::select(c(accession, gene_name, sample, differentiation, genotype, type, type_genotype, normalised_intensity_log2))
# 
# 
# data_ro_norm_filtered_wide = exo_data_ro_norm %>%
#   filter(gene_name %in% module_gene_select$gene_name) %>%
#   dplyr::filter(type  == "exosome") %>%
#   mutate(normalised_intensity_log2 = ifelse(is.na(normalised_intensity_log2),0,normalised_intensity_log2)) %>%
#   pivot_wider(id_cols = gene_name,
#               names_from=sample,
#               values_from=normalised_intensity_log2,
#               values_fill = 0) %>%
#   column_to_rownames("gene_name")
# 
# # annotation
# head(module_gene_select)
# 
# plot_annotation = exo_data_ro_norm %>%
#   dplyr::filter(type  == "exosome") %>%
#   dplyr::select(c(sample, genotype, differentiation)) %>%
#   unique() %>%
#   as.data.frame()
# 
# rownames(plot_annotation) = plot_annotation$sample
# plot_annotation$sample = NULL
# 
# ann_colors = list(
#   differentiation = c("1" = "red","2"="blue","3"="yellow", "4"="green"),
#   genotype = c("MUT" = "black", "WT" = "pink")
# )
# 
# # make heatmap
# pheatmap(data_ro_norm_filtered_wide,
#          annotation_col = plot_annotation,
#          annotation_colors = ann_colors,
#          show_rownames = T,
#          fontsize_row = 10,
#          scale = "row")
# 
# 
# 
# 
# 
# 
# # overlap with the DEG -----
# module_genes_deg = module_genes %>%
#   dplyr::filter(gene_name %in% exo_deg_filter_subtract$gene_name)
# 
# table(module_genes_deg$color)/table(module_genes$color)
# 
# 
