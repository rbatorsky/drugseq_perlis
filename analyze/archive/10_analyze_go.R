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

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# meta ----
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")  %>%
  mutate(drug = paste0(Drug_Treatment, "_", Concentration)) 

# moa ----
moa = read.xlsx("steve_moa_format.xlsx")
meta_moa = meta %>% 
  full_join(moa, by = "Drug_Treatment")

meta_moa = meta_moa %>% 
  mutate(match_found = ifelse(is.na(moa) | is.na(sample_name),'n','y')) %>%
  dplyr::select(-name_found)

#write.xlsx(meta_moa, "rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format_moa_16jul24.xlsx")

# active res ----
active = read.xlsx("rd1_rd2_analysis/de/active_compound_list_addscreen2_29jul24.xlsx") %>%
  filter(active_or_confirmed_2sd == 1)
head(active)

meta_moa_active = meta_moa %>% 
  mutate(active = ifelse(drug %in% active$treatment,1,0)) %>%
  dplyr::filter(active == 1)

head(meta_moa_active)

# with dir
ck = readRDS("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp_dir_29jul24.rds")

dotplot(ck, x="drug", show=1, by="Count")+ facet_grid(~dir) + 
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) 

# separate dotplot for the two groups, what are they doing differently?
g1 = groups_df %>%
  dplyr::filter(groups == 1)

g2 = groups_df %>%
  dplyr::filter(groups == 2)

ck_g1 = ck %>%
  dplyr::filter(drug %in% g1$treatment)

ck_g2 = ck %>%
  dplyr::filter(drug %in% g2$treatment)


dotplot(ck_g1, x="drug", show=1, by="Count")+ facet_grid(~dir) + 
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) 


dotplot(ck_g2, x="drug", show=1, by="Count")+ facet_grid(~dir) + 
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) 

# which specific processes are different in the two groups
ck_groups = ck@compareClusterResult %>%
  left_join(groups_df, by=c("drug" = "treatment"))

write.xlsx(ck_groups, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp_dir_29jul24_addgroups.xlsx")





# go ----
ck = readRDS("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp_29jul24.rds") %>%
  dplyr::filter(Cluster %in% active$treatment)

length(unique(ck@compareClusterResult$Cluster))

ck_active = ck %>%
  dplyr::filter(Cluster %in% active$treatment)

write.xlsx(ck_active, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp_active_or_confirmed_29jul24.xlsx")

dotplot(ck_active, show=3, by="Count")+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

# cluster treatments by enrichment results using GeneRatio ------
# this is really noisy since so many are missing
colnames(ck_active@compareClusterResult)
ck_active_select = ck_active@compareClusterResult %>%
  select(drug, Description, GeneRatio) %>%
  separate(GeneRatio, into=c("num","den"), sep="/") %>%
  mutate(GeneRatio = as.numeric(num)/as.numeric(den)) %>%
  pivot_wider(id_cols=Description, values_from=GeneRatio, names_from=drug, values_fill = 0) %>%
  column_to_rownames("Description") 

dim(ck_active_select)

# correlation
cormat <- cor(ck_active_select)
pheatmap(cormat)
 
# pheatmap(as.matrix(ck_active_select_g2), cluster_rows = F, cluster_cols=F)

# which enrichment results correlate the most with the phagocytosis score ------
ck = readRDS("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2_gobp_29jul24.rds") 

all_res = read.xlsx("rd1_rd2_analysis/de/active_compound_list_addscreens_29jul24.xlsx") 

ck_active_select = ck_active@compareClusterResult %>%
  select(drug, ID, GeneRatio) %>%
  separate(GeneRatio, into=c("num","den"), sep="/") %>%
  mutate(GeneRatio = as.numeric(num)/as.numeric(den)) %>%
  pivot_wider(id_cols=drug, values_from=GeneRatio, names_from=ID, values_fill = 0)  %>%
  left_join(all_res %>%
              select(treatment, Compound_mean, mean_phago, ndeg), by=c("drug" = "treatment")) %>%
  filter(!is.na(Compound_mean)) %>%
  column_to_rownames("drug")

rownames(ck_active_select)
ck_active_select$Compound_mean

cs = colSums(ck_active_select > 0)
ck_active_select_g = ck_active_select[ , which(cs>15)]

tt = lm(Compound_mean ~ ., data=ck_active_select_g) 
summary(tt)


# We want to find a process that is predictive of phagocytosis scores
# do the # of DEG in the phagocytosis pathway or GeneRatio in phagocytosis pathway correlate with the score?
# If so, try on other pathways, if not then it's not a good metric
# Go back to the sample grouping over genes, cut the tree to reveal samples over variable genes
# Are there any genes that correlate with phagocytosis scores? We tried this with genes in the phagocytosis pathway and didn't find too much


# make the phagocytosis plots ------
deg = fread("rd1_rd2_analysis/de/all_res_run2.csv") %>%
  separate(gene, into = c("ens","gene"), sep="\\|", remove=F)

head(deg)
go_phago_run1 = ck %>%
  dplyr::filter(grepl('phago',Description, ignore.case=T)) %>%
  dplyr::filter(!(grepl('run2', Cluster)))

# go_phago_run2 = ck %>%
#   dplyr::filter(grepl('phago',Description, ignore.case=T)) %>%
#   dplyr::filter(grepl('run2', Cluster))


dotplot(go_phago_run2, by="Count")+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

# phagocytosis heatmap
term="phagocytosis"
ck_select = ck %>% 
  dplyr::filter(Description == term)
ck_select
ck_select_genes = data.frame(genes = ck_select@compareClusterResult$geneID)
ck_select_genes$split = gsub("/",",",ck_select_genes$genes)
ck_select_genes_list = unique(unlist(strsplit(ck_select_genes$split,",")))

ck_select_genes_list

head(deg)
deg_select = deg %>%
  dplyr::filter(gene %in% ck_select_genes_list) %>%
  arrange(factor(drug, levels=unique(deg$drug))) %>%
  dplyr::select("gene","log2FoldChange","padj", "drug")

head(deg_select)

#deg_select$cell_type= factor(deg_select$cell_type, levels=select_cluster)

lfc_select = deg_select %>%
  dplyr::select("gene","log2FoldChange", "drug") %>%
  spread("drug","log2FoldChange")

rownames(lfc_select) = lfc_select$gene
lfc_select$gene = NULL
lfc_select_mat = as.matrix(lfc_select)
lfc_select_mat[is.na(lfc_select_mat)] <- 0

p_select = deg_select %>%
  dplyr::select("gene","padj", "drug") %>%
  spread("drug","padj")

rownames(p_select) = p_select$gene
p_select$gene = NULL
p_select_mat = as.matrix(p_select)
p_select_mat[is.na(p_select_mat)] <- 1

heatmap_max = max(deg_select$log2FoldChange)
heatmap_min = min(deg_select$log2FoldChange)
heatmap_max
library(circlize)
library(ComplexHeatmap)

lfc_select_mat = t(lfc_select_mat)
p_select_mat = t(p_select_mat)
ph = Heatmap(lfc_select_mat,
             col = circlize::colorRamp2(c(-10, 0, 10), c("Darkblue", "white", "red")),
             cluster_columns=T,
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(p_select_mat[i, j] < 0.05 & abs(lfc_select_mat[i,j])>0.2) {
                 grid.text("*", x, y)
               }
             },
             heatmap_legend_param = list(title = "Scaled log\nfold change"),
             column_names_gp = grid::gpar(fontsize = 8),
             row_names_gp = grid::gpar(fontsize = 8) ,
             column_title = term)

pdf("rd1_rd2_analysis/de/phago_heatmap.pdf", height=16, width=7)
print(ph)
dev.off()

# phago_index 
# move this over to the 11_compare file
pi = read.csv("phenotypic_screen_data/primary_screen_functional_means.csv") %>%
  rename(Drug_Treatment = Compound_Name) %>%
  mutate(Drug_Treatment = gsub('\\)','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\(','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\-','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\+','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('_$','',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('^_','',Drug_Treatment)) 

head(pi)

cor(lfc_select_mat)


g="YES1"
lfc = lfc_select_mat %>%
  as.data.frame() %>%
  rownames_to_column("drug")  %>%
  left_join(meta %>% dplyr::select("Drug_Treatment","drug"), by="drug") %>%
  inner_join(pi %>% dplyr::select(Drug_Treatment, Compound_mean), by="Drug_Treatment") %>%
  dplyr::select(-Drug_Treatment) %>%
  dplyr::select(drug, Compound_mean, everything()) %>%
  distinct() %>%
  column_to_rownames("drug")

rn = rownames(lfc)
lfc <- sapply( lfc, as.numeric )
rownames(lfc) = rn
head(lfc)

lfc = data.frame(lfc)
cm = lfc[,1]
head(cm)
dat = lfc[,-1]
head(dat)
tt = lm(Compound_mean ~ ., data=lfc) 

summary(tt) 

# why are there correlations
tmp <- cor(lfc)

sort(tmp[,"Compound_mean"])

#LYST         SPON2          AZU1 CLEC7A           BTK         RAB34


tt = lm(Compound_mean ~ LYST+SPON2+AZU1+CLEC7A+BTK+RAB34, data=lfc) 
summary(tt)


library(ggpubr)

ggscatter(lfc, x = "SPON2", y = "Compound_mean",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n")
)


#%>%
#  dplyr::select(-Drug_Treatment) %>%
#  column_to_rownames("drug")

head(lfc)  

