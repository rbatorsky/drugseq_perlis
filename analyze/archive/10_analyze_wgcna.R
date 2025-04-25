# run wgcna on the cluster, that's where the impute library is installed

library(tidyverse)
library(clusterProfiler)
library(openxlsx)
library(org.Hs.eg.db)

# wgcna genes heatmap -----
setwd("/cluster/tufts/patralab/rbator01/degterev/necroptosis_resistant_ht29_3x_5x_jul24//")
filename="analysis/gene_modules.txt"
module_genes = read.table(filename, sep="\t", header=T) 
ids=module_genes$gene_id
sym=mapIds(org.Hs.eg.db, keys = ids, keytype="ENSEMBL", column = "SYMBOL")
module_genes$symbol = sym


blue_gene = module_genes %>% 
  dplyr::filter(colors=='blue')



ck = clusterProfiler::compareCluster(gene_id ~ colors,
                                     data = module_genes,
                                     fun="enrichGO",
                                     OrgDb='org.Hs.eg.db',
                                     keyType = "ENSEMBL",
                                     ont="BP",
                                     readable = T)

dotplot(ck, show=5, size="Count") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 1000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

saveRDS(ck, "analysis/wgcna/go_bp_module_genes_14jul24.rds")
write.xlsx(ck, "analysis/wgcna/go_bp_module_genes_14jul24.xlsx")

ck = readRDS("analysis/wgcna/go_bp_module_genes_14jul24.rds")

head(ck)

methyl = ck %>% 
  dplyr::filter(grepl('meth',Description, ignore.case=T))

view(methyl)

ck_blue = ck %>%
  dplyr::filter(colors == "blue") 

dotplot(ck_blue, show = 10)

methyl = ck_blue %>% 
  dplyr::filter(grepl('meth',Description, ignore.case=T))


head(methyl)

ck_blue_ripk3 = ck_blue %>%
  dplyr::filter(grepl('RIPK3', geneID))


test = ck %>%
  dplyr::filter(grepl('DNMT1', geneID))

view(test)


dotplot(ck_blue_ripk3, show = 10)

view(ck_blue_ripk3)
# plot the module genes
c="purple"
module_gene_select = module_genes %>%
  dplyr::filter(colors == c)

library(pheatmap)

exo_data_ro_norm = read.xlsx("exosome/data_ro_norm_exo_combinedexopell_rm_mut5_6may24.xlsx") %>%
  mutate(sample = paste0("exo_", sample))%>%
  mutate(type = "exosome") %>%
  dplyr::select(c(accession, gene_name, sample, differentiation, genotype, type, type_genotype, normalised_intensity_log2))


data_ro_norm_filtered_wide = exo_data_ro_norm %>%
  filter(gene_name %in% module_gene_select$gene_name) %>%
  dplyr::filter(type  == "exosome") %>%
  mutate(normalised_intensity_log2 = ifelse(is.na(normalised_intensity_log2),0,normalised_intensity_log2)) %>%
  pivot_wider(id_cols = gene_name,
              names_from=sample,
              values_from=normalised_intensity_log2,
              values_fill = 0) %>%
  column_to_rownames("gene_name")

# annotation
head(module_gene_select)

plot_annotation = exo_data_ro_norm %>%
  dplyr::filter(type  == "exosome") %>%
  dplyr::select(c(sample, genotype, differentiation)) %>%
  unique() %>%
  as.data.frame()

rownames(plot_annotation) = plot_annotation$sample
plot_annotation$sample = NULL

ann_colors = list(
  differentiation = c("1" = "red","2"="blue","3"="yellow", "4"="green"),
  genotype = c("MUT" = "black", "WT" = "pink")
)

# make heatmap
pheatmap(data_ro_norm_filtered_wide, 
         annotation_col = plot_annotation, 
         annotation_colors = ann_colors,
         show_rownames = T,
         fontsize_row = 10, 
         scale = "row")










# overlap with the DEG

exo_deg = read.xlsx('../differential_abundance/mixed_effect_exo_combinedexopell_rm_mut5_6may24.xlsx') %>%
  mutate(type="exosome") %>%
  dplyr::select(c(accession, gene_name,Estimate , p.value, p.adjust, type, dir))%>%
  mutate(gene_name = gsub(";.*","", gene_name))

pell_deg = read.xlsx('../differential_abundance/mixed_effect_pell_combinedexopell_rm_mut6_7may24.xlsx') %>%
  mutate(type = "pellet") %>%
  dplyr::select(c(accession, gene_name,Estimate , p.value, p.adjust, type, dir)) %>%
  mutate(gene_name = gsub(";.*","", gene_name))

# subtract pellet DEG and contaminants -----
padj_cut=0.1
lfc_cut=0

deg_filter = pell_deg  %>%
  dplyr::filter(p.adjust <padj_cut & abs(Estimate) > lfc_cut) 

module_genes_deg = module_genes %>%
  dplyr::filter(gene_name %in% deg_filter$gene_name)

table(module_genes_deg$colors)

ck = clusterProfiler::compareCluster(accession ~ colors,
                                     data = module_genes_deg,
                                     fun="enrichGO",
                                     OrgDb='org.Hs.eg.db',
                                     keyType = "UNIPROT",
                                     ont="BP",
                                     readable = T)

dotplot(ck, show=5, size="Count") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 1000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

saveRDS(ck, paste0(save_string,"/go_bp_module_genes_deg_28may24.rds"))
write.xlsx(ck, paste0(save_string,"/go_bp_module_genes_deg_28may24.xlsx"))

ck = readRDS(paste0(save_string,"/go_bp_module_genes_deg_28may24.rds"))

view(ck)

# plot the module genes

unique(module_genes$colors)
c="greenyellow"
module_gene_select = module_genes_deg  %>%
  dplyr::filter(colors == c)

library(pheatmap)

exo_data_ro_norm = read.xlsx("exosome/data_ro_norm_exo_combinedexopell_rm_mut5_6may24.xlsx") %>%
  mutate(sample = paste0("exo_", sample))%>%
  mutate(type = "exosome") %>%
  dplyr::select(c(accession, gene_name, sample, differentiation, genotype, type, type_genotype, normalised_intensity_log2))


data_ro_norm_filtered_wide = exo_data_ro_norm %>%
  filter(gene_name %in% module_gene_select$gene_name) %>%
  dplyr::filter(type  == "exosome") %>%
  mutate(normalised_intensity_log2 = ifelse(is.na(normalised_intensity_log2),0,normalised_intensity_log2)) %>%
  pivot_wider(id_cols = gene_name,
              names_from=sample,
              values_from=normalised_intensity_log2,
              values_fill = 0) %>%
  column_to_rownames("gene_name")

# annotation
head(module_gene_select)

plot_annotation = exo_data_ro_norm %>%
  dplyr::filter(type  == "exosome") %>%
  dplyr::select(c(sample, genotype, differentiation)) %>%
  unique() %>%
  as.data.frame()

rownames(plot_annotation) = plot_annotation$sample
plot_annotation$sample = NULL

ann_colors = list(
  differentiation = c("1" = "red","2"="blue","3"="yellow", "4"="green"),
  genotype = c("MUT" = "black", "WT" = "pink")
)

# make heatmap
pheatmap(data_ro_norm_filtered_wide, 
         annotation_col = plot_annotation, 
         annotation_colors = ann_colors,
         show_rownames = T,
         fontsize_row = 10, 
         scale = "row")






# overlap with the DEG -----
module_genes_deg = module_genes %>%
  dplyr::filter(gene_name %in% exo_deg_filter_subtract$gene_name)

table(module_genes_deg$color)/table(module_genes$color)


