# ARCHIVE below and move to IPA ----
# d <- dist(t(lfc_select_mat))
# hc <- hclust(d)
# plot(hc)
# rect.hclust(hc, k = 2)
# 
# groups = cutree(hc, k = 2)
# groups_df = data.frame(groups) %>%
#   rownames_to_column("treatment")
# groups_df
# # Is the phago stuff different between groups?
# 
# active = read.xlsx("rd1_rd2_analysis/de/active_compound_list_addscreens_29jul24.xlsx") %>%
#   filter(active_or_confirmed_2sd == 1)
# 
# active_groups = active %>%
#   inner_join(groups_df, by="treatment") 
# 
# head(active_groups)
# ggboxplot(active_groups, x = "groups", y = "ndeg",
#           color = "groups", palette = "jco",
#           add = "jitter") + stat_compare_means() + 
#   ggtitle("screen2")
# 
# 
# 
# # ck_active = ck %>%
# #   dplyr::filter(Cluster %in% active$treatment)
# # 
# # go_phago = ck_active %>%
# #   dplyr::filter(grepl('phago',Description, ignore.case=T)) 
# # 
# # dotplot(go_phago, by="Count")+
# #   theme(text = element_text(size = 10),
# #         axis.text.x = element_text(angle = 90, size=10, hjust = 1),
# #         axis.text.y = element_text(size=10)) +
# #   scale_y_discrete(labels = function(x) str_wrap(x, width = 100))
# # 
# # unique(go_phago@compareClusterResult$Description)
# 
# # phagocytosis heatmap
# term="phagocytosis"
# #term="positive regulation of cytokine production"
# 
# ck_select = ck_active %>% 
#   dplyr::filter(Description == term)
# ck_select
# ck_select_genes = data.frame(genes = ck_select@compareClusterResult$geneID)
# ck_select_genes
# ck_select_genes$split = gsub("/",",",ck_select_genes$genes)
# ck_select_genes_list = unique(unlist(strsplit(ck_select_genes$split,",")))
# 
# ck_select_genes_list
# 
# deg_select = deg %>%
#   dplyr::filter(gene %in% ck_select_genes_list) %>%
#   arrange(factor(Compound_Name, levels=unique(deg$Compound_Name))) %>%
#   dplyr::select("gene","log2FoldChange","padj", "Compound_Name")
# 
# lfc_select = deg_select %>%
#   dplyr::select("gene","log2FoldChange", "Compound_Name") %>%
#   spread("Compound_Name","log2FoldChange")
# 
# rownames(lfc_select) = lfc_select$gene
# lfc_select$gene = NULL
# lfc_select_mat = as.matrix(lfc_select)
# lfc_select_mat[is.na(lfc_select_mat)] <- 0
# 
# p_select = deg_select %>%
#   dplyr::select("gene","padj", "Compound_Name") %>%
#   spread("Compound_Name","padj")
# 
# rownames(p_select) = p_select$gene
# p_select$gene = NULL
# p_select_mat = as.matrix(p_select)
# p_select_mat[is.na(p_select_mat)] <- 1
# 
# heatmap_max = max(deg_select$log2FoldChange)
# heatmap_min = min(deg_select$log2FoldChange)
# heatmap_max
# library(circlize)
# library(ComplexHeatmap)
# 
# lfc_select_mat = t(lfc_select_mat)
# p_select_mat = t(p_select_mat)
# ph = Heatmap(lfc_select_mat,
#              col = circlize::colorRamp2(c(-10, 0, 10), c("Darkblue", "white", "red")),
#              cluster_columns=T,
#              cell_fun = function(j, i, x, y, w, h, fill) {
#                if(p_select_mat[i, j] < 0.05 & abs(lfc_select_mat[i,j])>0.2) {
#                  grid.text("*", x, y)
#                }
#              },
#              heatmap_legend_param = list(title = "Scaled log\nfold change"),
#              column_names_gp = grid::gpar(fontsize = 8),
#              row_names_gp = grid::gpar(fontsize = 8) ,
#              column_title = term)
# 
# #pdf("rd1_rd2_analysis/de/phago_heatmap.pdf", height=16, width=7)
# print(ph)
# #dev.off()
# 
# head(lfc_select_mat)
# d <- dist(lfc_select_mat)
# hc <- hclust(d)
# plot(hc)
# rect.hclust(hc, k = 2)
# 
# groups = cutree(hc, k = 2)
# groups_df = data.frame(groups) %>%
#   rownames_to_column("treatment")
# groups_df
# # Is the phago stuff different between groups?
# 
# active = read.xlsx("rd1_rd2_analysis/de/active_compound_list_addscreens_29jul24.xlsx") %>%
#   filter(active_or_confirmed_2sd == 1)
# 
# active_groups = active %>%
#   inner_join(groups_df, by="treatment") 
# 
# head(active_groups)
# ggboxplot(active_groups, x = "groups", y = "ndeg",
#           color = "groups", palette = "jco",
#           add = "jitter") + stat_compare_means() + 
#   ggtitle("screen2")
# 
# 
# 
# 
# # phago_index -----
# # move this over to the 11_compare file
# pi = read.csv("phenotypic_screen_data/primary_screen_functional_means.csv") %>%
#   rename(Drug_Treatment = Compound_Name) %>%
#   mutate(Drug_Treatment = gsub('\\)','_',Drug_Treatment))%>%
#   mutate(Drug_Treatment = gsub('\\(','_',Drug_Treatment))%>%
#   mutate(Drug_Treatment = gsub('\\-','_',Drug_Treatment))%>%
#   mutate(Drug_Treatment = gsub('\\+','_',Drug_Treatment))%>%
#   mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
#   mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
#   mutate(Drug_Treatment = gsub('_$','',Drug_Treatment))%>%
#   mutate(Drug_Treatment = gsub('^_','',Drug_Treatment)) 
# 
# head(pi)
# 
# cor(lfc_select_mat)
# 
# 
# g="YES1"
# lfc = lfc_select_mat %>%
#   as.data.frame() %>%
#   rownames_to_column("drug")  %>%
#   left_join(meta %>% dplyr::select("Drug_Treatment","drug"), by="drug") %>%
#   inner_join(pi %>% dplyr::select(Drug_Treatment, Compound_mean), by="Drug_Treatment") %>%
#   dplyr::select(-Drug_Treatment) %>%
#   dplyr::select(drug, Compound_mean, everything()) %>%
#   distinct() %>%
#   column_to_rownames("drug")
# 
# rn = rownames(lfc)
# lfc <- sapply( lfc, as.numeric )
# rownames(lfc) = rn
# head(lfc)
# 
# lfc = data.frame(lfc)
# cm = lfc[,1]
# head(cm)
# dat = lfc[,-1]
# head(dat)
# tt = lm(Compound_mean ~ ., data=lfc) 
# 
# summary(tt) 
# 
# # why are there correlations
# tmp <- cor(lfc)
# 
# sort(tmp[,"Compound_mean"])
# 
# #LYST         SPON2          AZU1 CLEC7A           BTK         RAB34
# 
# 
# tt = lm(Compound_mean ~ LYST+SPON2+AZU1+CLEC7A+BTK+RAB34, data=lfc) 
# summary(tt)
# 
# 
# library(ggpubr)
# 
# ggscatter(lfc, x = "SPON2", y = "Compound_mean",
#           color = "black", shape = 21, size = 3, # Points color, shape and size
#           add = "reg.line",  # Add regressin line
#           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#           conf.int = TRUE, # Add confidence interval
#           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
#           cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n")
# )
# 
# 
# #%>%
# #  dplyr::select(-Drug_Treatment) %>%
# #  column_to_rownames("drug")
# 
# head(lfc)  
# 
