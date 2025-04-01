LIB="/cluster/tufts/patralab/rbator01/R_libs/4.4.0/"
.libPaths(LIB)
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')
#setwd('~/Box/0my_projects/draft_pres/perlis_prep/drugseq_mar24/')

library(openxlsx)
library(tidyverse)
library(data.table)
library(circlize)
library(ComplexHeatmap)
select = dplyr::select
rename = dplyr::rename
library(ggpubr)

# combine ipa results -----
files <- Sys.glob("rd1_rd2_analysis/ipa/*cp.txt")
files
df =NULL
files
for (f in files) {
  print(f)
  new=fread(f)
  new$rna_name = gsub("_cp.txt","",basename(f))
  if(nrow(new)>0)
    if(is.null(df))
      df <- new
  else
    df = rbind(df, new)
}

df = df %>%
  rename(ipa_canonical = `Ingenuity Canonical Pathways`) %>%
  rename(z.score = `z-score`)

# df  = df %>%
#   mutate(ipa_canonical = ifelse(ipa_canonical == 'Fc Receptor-mediated Phagocytosis in Macrophages and Monocytes',
#                                 'Fc_gamma Receptor-mediated Phagocytosis in Macrophages and Monocytes',ipa_canonical))

df$count = sapply(strsplit(as.character(df$Molecules), ","), length)
write.xlsx(df, "rd1_rd2_analysis/ipa/combined_cp_scr2_50.xlsx")

##Fig 5C select steve and rebecca pathways ----
all = read.xlsx("rd1_rd2_analysis/ipa/combined_cp_scr2_50_active_select_v5_sds.xlsx") %>%
  mutate(edited_name = gsub("TopotecanHydrochloride", "TopotecanHCL", edited_name))

# make it so that all selected have 1
steve_paths=all %>%
  filter(steve_select == 1) %>%
  .$ipa_canonical

rebecca_paths=all %>%
  filter(rebecca_select == 1) %>%
  .$ipa_canonical
rebecca_paths
select_paths = unique(c(steve_paths, rebecca_paths))

df_top = all %>%
  filter(ipa_canonical %in% select_paths)

df_top_heat = df_top %>%
  select(ipa_canonical, edited_name, z.score) %>%
  pivot_wider(names_from=edited_name, values_from=z.score, id_cols=ipa_canonical) %>%
  column_to_rownames("ipa_canonical")

df_top_heat[is.na(df_top_heat)] = 0
head(df_top_heat)
phenos=c('mean_normalized_PI','well_mean_eccentricity')

annot = df_top %>%
  select(edited_name, all_of(phenos)) %>%
  distinct() %>%
  column_to_rownames("edited_name") %>%
  as.data.frame()

head(df_top)
row_annot = df_top %>%
  select(ipa_canonical, color) %>%
  distinct() %>%
  column_to_rownames("ipa_canonical") %>%
  as.data.frame()

head(row_annot)

min_heat = min(df_top$z.score, na.rm=T)
max_heat = max(df_top$z.score, na.rm=T)
min_heat
max_heat
col_fun = colorRamp2(c(min_heat, 0, max_heat), c("blue", "white", "red"))

head(annot)
bot_annot = HeatmapAnnotation(`Phagocytoc Index` = anno_barplot(annot$mean_normalized_PI),
                              annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
                              annotation_name_gp= gpar(fontsize = 8))

#df_top_heat_abs_rs = rowSums(abs(df_top_heat))

#rt_annot = rowAnnotation(category = row_annot)

row_ha = rowAnnotation(category = row_annot$color,
                       col = list(category = c("actin" = "black",
                                          "activation" = "pink",
                                          "cell cycle" = "cyan",
                                          "cell death" = "darkkhaki",
                                          "immune response" = "green",
                                          "lysozome" = "darkorange",
                                          "metabolism" = "purple",
                                          "phagocytosis" = "yellow")),
                       annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
                       annotation_name_gp= gpar(fontsize = 0))


#ha = rowAnnotation(foo = anno_empty(border = FALSE,
#                                    width = max_text_width(unlist(text_list)) + unit(4, "mm")))

pdf("rd1_rd2_analysis/paper_plots/5_c_steve_select.pdf", height=6, width=9)
p = Heatmap(df_top_heat,
            show_heatmap_legend = T,
            col = col_fun,
            name="IPA z-score",
            column_names_gp = grid::gpar(fontsize = 10),
            row_names_gp = grid::gpar(fontsize = 10),
            heatmap_legend_param = list(title_gp = gpar( fontsize = 8)),
            row_names_max_width = unit(12, "cm"),
            bottom_annotation = bot_annot,
            right_annotation = row_ha,
            row_dend_side = "left")
print(p)
dev.off()

# # gene heatmap of phago cat -----

#select_cat ="Fcgamma receptor (FCGR) dependent phagocytosis"
select_cat ="Phagosome Formation"

df_top = all %>%
  filter(ipa_canonical == select_cat)

df_top_genes = data.frame(genes = df_top$Molecules)
df_top_genes_list = unique(unlist(strsplit(df_top_genes$genes,",")))
length(df_top_genes_list)
df_top_genes_list[grep("ACT",df_top_genes_list)]

deg = fread("rd1_rd2_analysis/de/all_res_run3.csv") %>%
  separate(gene, into=c("ens","gene"), sep="\\|") 

test =  fread("rd1_rd2_analysis/de/all_res_run3.csv") %>%
  dplyr::filter(gene %in% df_top_genes_list)


  
  # mutate(Compound_Name = gsub("_10uM", "", drug)) %>%
  # mutate(Compound_Name = gsub("_2uM", "", Compound_Name)) %>%
  # mutate(Compound_Name = gsub("_run3.csv", "", Compound_Name)) %>%
  # separate(gene, into=c("ens","gene"), sep="\\|") %>%
  #mutate(drug = gsub("Topotecan_Hydrochloride", "TopotecanHCL", drug))

top_deg = deg %>%
  filter(padj < 0.000000001 & abs(log2FoldChange)>1.5)

genes_to_plot = intersect(df_top_genes_list, top_deg$gene)
length(unique(genes_to_plot))
deg_select_all = deg %>%
  dplyr::filter(gene %in% genes_to_plot) %>%
  mutate(drug = gsub("_10uM","", drug)) %>%
  mutate(drug = gsub("_2uM","", drug)) %>%
  mutate(drug = gsub("_acid","_Acid",drug)) %>%
  mutate(drug = gsub("_ditosylate","_Ditosylate",drug)) %>%
  mutate(drug = gsub("_tosylate","_Tosylate",drug)) %>%
  mutate(drug = gsub("Topotecan_Hydrochloride", "TopotecanHCL", drug)) %>%
  mutate(drug = gsub("_","",drug)) %>%
  filter(drug %in% colnames(df_top_heat)) %>%
  dplyr::select("gene","log2FoldChange","padj", "drug")

unique(deg_select_all$drug)

setdiff(colnames(df_top_heat), unique(deg_select_all$drug))
setdiff( unique(deg_select_all$drug),colnames(df_top_heat))

# now choose the limits
heatmap_max = max(deg_select_all$log2FoldChange, na.rm = T)
heatmap_min = min(deg_select_all$log2FoldChange, na.rm = T)

lfc_select = deg_select_all %>%
  dplyr::select("gene","log2FoldChange", "drug") %>%
  spread("drug","log2FoldChange")

rownames(lfc_select) = lfc_select$gene
lfc_select$gene = NULL
lfc_select_mat = as.matrix(lfc_select)
lfc_select_mat[is.na(lfc_select_mat)] <- 0
p_select = deg_select_all %>%
  dplyr::select("gene","padj", "drug") %>%
  spread("drug","padj")
rownames(p_select) = p_select$gene
p_select$gene = NULL
p_select_mat = as.matrix(p_select)
p_select_mat[is.na(p_select_mat)] <- 1

dim(lfc_select_mat)
library(circlize)
library(ComplexHeatmap)
pdf("rd1_rd2_analysis/paper_plots/phago_form_genes.pdf", height=9.5, width=5)
ph = Heatmap(lfc_select_mat,
             col = circlize::colorRamp2(c(heatmap_min, 0, heatmap_max), c("Darkblue", "white", "red")),
             cluster_columns=T,
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(p_select_mat[i, j] < 0.05 & abs(lfc_select_mat[i,j])>0.2) {
                 grid.text("*", x, y)
               }
             },
             heatmap_legend_param = list(title = "Log2FoldChange\nCompound vs. DMSO"),
             column_names_gp = grid::gpar(fontsize = 10),
             row_names_gp = grid::gpar(fontsize = 10))


print(ph)
dev.off()


# Phagosome formation vs. Neuroinflammation -----
df_phago =  all %>%
  filter(ipa_canonical %in% c("Fcgamma receptor (FCGR) dependent phagocytosis","Neuroinflammation Signaling Pathway")) %>%
  select(edited_name, ipa_canonical, z.score) %>%
  pivot_wider(id_cols = "edited_name", values_from="z.score", names_from = "ipa_canonical")

colnames(df_phago) = make.names(colnames(df_phago))

library(ggpubr)
p = ggscatter(df_phago, x = "Fcgamma.receptor..FCGR..dependent.phagocytosis", y = "Neuroinflammation.Signaling.Pathway",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          label = "edited_name",
          font.label = 8,
          repel=T
)+ stat_cor(method = "pearson", label.x = -2, label.y = -1) +
  ylab('Neuroinflammation Signaling Pathway') +
  xlab('Fc-Gamma receptor dependent phagocytosis')
show(p)
ggsave(p,filename="rd1_rd2_analysis/paper_plots/5d.pdf" ,width=5,height=5)

df_phago =  all %>%
  filter(ipa_canonical %in% c("Phagosome Formation","Neuroinflammation Signaling Pathway")) %>%
  select(edited_name, ipa_canonical, z.score) %>%
  pivot_wider(id_cols = "edited_name", values_from="z.score", names_from = "ipa_canonical")

colnames(df_phago) = make.names(colnames(df_phago))

library(ggpubr)
p = ggscatter(df_phago, x = "Phagosome.Formation", y = "Neuroinflammation.Signaling.Pathway",
              add = "reg.line",
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              label = "edited_name",
              font.label = 8,
              repel=T
)+ stat_cor(method = "pearson", label.x =-1, label.y = -1) +
  ylab('Neuroinflammation Signaling Pathway') +
  xlab('Phagosome Formation')
show(p)
ggsave(p,filename="rd1_rd2_analysis/paper_plots/5d_phago_form.pdf" ,width=5,height=5)

# linear modeling -----
# problem with this is it's trying to find pathways that explain the loss of phago in all
# we suspect there are multiple

df = read.xlsx("rd1_rd2_analysis/ipa/combined_cp_scr2_50.xlsx") %>%
  mutate(edited_name = gsub("_hydrochloride","_HCL",rna_name)) %>%
  mutate(edited_name = gsub("_acid","_Acid",edited_name)) %>%
  mutate(edited_name = gsub("_ditosylate","_Ditosylate",edited_name)) %>%
  mutate(edited_name = gsub("_tosylate","_Tosylate",edited_name)) %>%
  mutate(edited_name = gsub("_","",edited_name))

activity = read.xlsx("rd1_rd2_analysis/de/all_compounds_rnaactivity_functionalscreen_23aug24.xlsx") %>%
  select(Compound_Name, ndeg)

scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format_add_rna_activity.xlsx") %>%
  inner_join(activity, by="Compound_Name") %>%
  rename(treatment = Compound_Name)


# scr2 = read.xlsx("phenotypic_screen_data/2024.8.29_SecondaryScreen_CombinedData_format.xlsx") %>%
#   inner_join(activity, by="Compound_Name") %>%
#   rename(treatment = Compound_Name) 

scr2_50 = scr2 %>%
  filter(PI_below_50percent == "yes" & ndeg > 200)

#%>%
#  left_join(moa %>% select(treatment, MOA), by="treatment")

# do we have all 2sd confirmed in the ipa results
intersect(unique(scr2_50$treatment), unique(df$rna_name))

setdiff(unique(df$rna_name), unique(scr2_50$treatment))
setdiff(unique(scr2_50$treatment), unique(df$rna_name))
# "Lapatinib"         

df = df %>%
  inner_join(scr2_50,
             by=c("rna_name" = "treatment"))

head(df)
phenos=c('mean_normalized_PI','well_mean_eccentricity','well_mean_normalized_IBA1_intensity')


df_significant = df %>%
  filter(`-log(p-value)` > 1.3) %>%
  filter(abs(z.score) > 2) %>%
  group_by(ipa_canonical) %>%
  summarize(count = n()) %>%
  filter(count > 2)

paths = unique(df_significant$ipa_canonical)

res_df = NULL
for (pheno in phenos){
  
  if(pheno =='well_mean_normalized_IBA1_intensity'){
    ymax = 0.1
  }else{ymax = 1}
  
  for (p in paths){
    
    path_i = df %>%
      filter(ipa_canonical == p)
    
    nvals = length(which(!is.na(path_i$z.score)))
    
    if(nvals > 5){
      
      lm = lm(paste0(pheno, '~', 'z.score'),data=path_i) 
      
      s = summary(lm)
      
      result = data.frame(Estimate = s$coefficients['z.score','Estimate'], pval = s$coefficients['z.score','Pr(>|t|)'], ipa_canonical = p)
      result$phenotype = pheno
      
      if(is.null(res_df)){
        res_df = result
      }else{
        res_df = rbind(res_df, result)
      }
      
      if(result$pval < 0.05){
        
        fig = ggscatter(path_i, x = "z.score", y = pheno,
                        add = "reg.line",  
                        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                        conf.int = TRUE, # Add confidence interval
                        label = "rna_name",
                        font.label = 8,
                        repel=T
        )+ stat_cor(method = "pearson", label.x = .2, label.y = ymax) + 
          ggtitle(p)
        
        p_name = gsub(' ','_',p)
        p_name = gsub('\\/','_',p_name)
        p_name = gsub('-','_',p_name)
        p_name = gsub('\\.','_',p_name)
        
        ggsave(fig, filename=paste0('rd1_rd2_analysis/ipa/',p_name,'_',pheno,'correlation.pdf'), height=5, width=5)
      }
    }
  }
}

res_df$p.adjust <- p.adjust(res_df$pval, method = 'BH')

write.xlsx(res_df,"rd1_rd2_analysis/ipa/combined_cp_correlation_with_phenos_50dmso_active_11nov24.xlsx")


# # read in ipa and add functional data -----
# 
# # df = read.xlsx("rd1_rd2_analysis/ipa/combined_cp_scr2_50.xlsx") %>%
# #   mutate(edited_name = gsub("_hydrochloride","_HCL",rna_name)) %>%
# #   mutate(edited_name = gsub("_acid","_Acid",edited_name)) %>%
# #   mutate(edited_name = gsub("_ditosylate","_Ditosylate",edited_name)) %>%
# #   mutate(edited_name = gsub("_tosylate","_Tosylate",edited_name)) %>%
# #   mutate(edited_name = gsub("_","",edited_name))
# 
# df = read.xlsx("rd1_rd2_analysis/ipa/combined_cp_scr2_50_active_select_v4_sds.xlsx") %>%
#   mutate(edited_name = gsub("TopotecanHydrochloride", "TopotecanHCL", edited_name)) 
# 
# 
# scr2_50_active = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format_add_rna_activity.xlsx")  %>%
#   filter(PI_below_50percent_calc == 1)%>%
#   filter(active == 1) %>%
#   mutate(Compound_Name = gsub("_10uM","", Compound_Name)) %>%
#   mutate(Compound_Name = gsub("_2uM","", Compound_Name)) %>%
#   mutate(Compound_Name = gsub("_hydrochloride","_HCL",Compound_Name)) %>%
#   mutate(Compound_Name = gsub("_Hydrochloride","_HCL",Compound_Name)) %>%
#   mutate(Compound_Name = gsub("_acid","_Acid",Compound_Name)) %>%
#   mutate(Compound_Name = gsub("_ditosylate","_Ditosylate",Compound_Name)) %>%
#   mutate(Compound_Name = gsub("_tosylate","_Tosylate",Compound_Name)) %>%
#   mutate(Compound_Name = gsub("_","",Compound_Name)) %>%
#   rename(treatment = Compound_Name)
# 
# 
# # do we have all 2sd confirmed in the ipa results
# intersect(unique(scr2_50_active$treatment), unique(df$edited_name))
# 
# setdiff(unique(df$edited_name), unique(scr2_50_active$treatment))
# setdiff(unique(scr2_50_active$treatment), unique(df$edited_name))
# 
# df = df %>%
#   inner_join(scr2_50_active,
#             by=c("edited_name" = "treatment"))
# 
# 
# # compare top categories -----
# top_cat = c()
# for(c in unique(df$rna_name)){
#   
#   df_c = df %>%
#     filter(rna_name == c) %>%
#     filter(!(ipa_canonical %in% top_cat)) %>%
#     filter(`-log(p-value)` > 1.3) %>%
#     arrange(-abs(z.score)) %>%
#     head(3)
#   
#   top_cat = c(top_cat, df_c$ipa_canonical)
# }
# 
# # heatmap of top categories -----
# 
# library(ComplexHeatmap)
# df_top = df %>%
#   filter(ipa_canonical %in% top_cat) 
# 
# #write.xlsx(df_top, "rd1_rd2_analysis/ipa/top3_cp_percompound.xlsx")
# 
# head(df_top)
# df_top_heat = df_top %>%
#   select(ipa_canonical, edited_name, z.score) %>%
#   pivot_wider(names_from=edited_name, values_from=z.score, id_cols=ipa_canonical) %>%
#   column_to_rownames("ipa_canonical")
# 
# df_top_heat[is.na(df_top_heat)] = 0
# head(df_top)
# phenos=c('mean_normalized_PI','well_mean_eccentricity','well_mean_solidity')
# annot = df_top %>%
#   select(edited_name, all_of(phenos)) %>%
#   distinct() %>%
#   column_to_rownames("edited_name") %>%
#   as.data.frame()
# 
# min_heat = min(df_top$z.score, na.rm=T)
# max_heat = max(df_top$z.score, na.rm=T)
# min_heat
# max_heat
# col_fun = colorRamp2(c(min_heat, 0, max_heat), c("blue", "white", "red"))
# 
# col_fun
# 
# head(annot)
# head(df_top_heat)
# 
# df_top_heat = as.matrix(df_top_heat)
# Heatmap(df_top_heat,
#         col = col_fun,
#         name="IPA z-score",
#         column_names_gp = grid::gpar(fontsize = 8),
#         row_names_max_width = unit(12, "cm"),
#         row_names_gp = grid::gpar(fontsize = 8),
#         heatmap_legend_param = list(title_gp = gpar( fontsize = 8)),
#         bottom_annotation = HeatmapAnnotation(df = annot,
#                                            annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
#                                            annotation_name_gp= gpar(fontsize = 8)))
# 
# 
# 
#  # study vorinostat vs. sr3306 ----
# 
# select_compound=c("Vorinostat","SR_3306")
# two_path = df %>%
#   filter(rna_name %in% select_compound) 
# 
# to_plot = two_path %>% 
#   filter(!is.na(z.score)) %>%
#   select(all_of(c("ipa_canonical", "z.score","rna_name"))) %>%
#   pivot_wider( values_from="z.score", names_from = "rna_name") %>%
#   filter(!is.na(Vorinostat) & !is.na(SR_3306)) %>%
#   filter(sign(Vorinostat) != sign(SR_3306)) %>%
#   filter(abs(Vorinostat) > 4 | abs(SR_3306) > 4) %>%
#   column_to_rownames("ipa_canonical")
# 
# to_plot[is.na(to_plot)] = 0
# 
# head(to_plot)
# phenos=c('mean_normalized_PI')
# head(scr2)
# annot = scr2 %>%
#   filter(treatment %in% select_compound) %>%
#   select(treatment, all_of(phenos)) %>%
#   distinct() %>%
#   column_to_rownames("treatment") %>%
#   as.data.frame()
# 
# min_heat = min(c(to_plot$SR_3306, to_plot$Vorinostat), na.rm=T)
# max_heat = max(c(to_plot$SR_3306, to_plot$Vorinostat), na.rm=T)
# min_heat
# max_heat
# col_fun = colorRamp2(c(min_heat, 0, max_heat), c("blue", "white", "red"))
# 
# col_fun
# 
# 
# to_plot = as.matrix(to_plot)
# 
# bot_annot = HeatmapAnnotation(`Phagocytic Index` = anno_barplot(annot$mean_normalized_PI),
#                   annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
#                   annotation_name_gp= gpar(fontsize = 8))
# 
# 
# Heatmap(to_plot,
#         col = col_fun,
#         name="IPA z-score",
#         column_names_gp = grid::gpar(fontsize = 8),
#         row_names_gp = grid::gpar(fontsize = 8),
#         heatmap_legend_param = list(title_gp = gpar( fontsize = 8)),
#         row_names_max_width = unit(12, "cm"),
#         bottom_annotation = bot_annot)
# 
# 
# # select pathways ----
# 
# phago_path = df %>%
#   filter(grepl('phago', ipa_canonical, ignore.case=TRUE)) %>%
#   select(ipa_canonical) %>%
#   distinct()
# 
# phago_path$ipa_canonical
# df_top = df %>%
#   filter(ipa_canonical %in% unique(phago_path$ipa_canonical)) 
# 
# df_top_heat = df_top %>%
#   select(ipa_canonical, rna_name, z.score) %>%
#   pivot_wider(names_from=rna_name, values_from=z.score, id_cols=ipa_canonical) %>%
#   column_to_rownames("ipa_canonical")
# 
# df_top_heat[is.na(df_top_heat)] = 0
# head(df_top)
# phenos=c('mean_normalized_PI','well_mean_area','well_mean_eccentricity','well_mean_solidity')
# annot = df_top %>%
#   select(rna_name, all_of(phenos)) %>%
#   distinct() %>%
#   column_to_rownames("rna_name") %>%
#   as.data.frame()
# 
# col_fun = colorRamp2(c(-6, 0, 2), c("blue", "white", "red"))
# 
# Heatmap(df_top_heat,
#         col = col_fun,
#         column_names_gp = grid::gpar(fontsize = 8),
#         row_names_gp = grid::gpar(fontsize = 8),
#         heatmap_legend_param = list(title_gp = gpar( fontsize = 8)),
#         bottom_annotation = HeatmapAnnotation(df = annot,
#                                               annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
#                                               annotation_name_gp= gpar(fontsize = 8)))
# 
# 
# 
# # types of categories -----
# 
# phago_cp = phago_path$`Canonical Pathways`
# phago_cp
# # which pathways have at least one compound > 1.5
# df_top_select = df %>%
#   filter(ipa_canonical %in% phago_cp) %>%
#   filter(abs(z.score) > 1.5)
# 
# df_top = df %>%
#   filter(ipa_canonical %in% unique(df_top_select$ipa_canonical)) 
# 
# df_top_heat = df_top %>%
#   select(ipa_canonical, rna_name, z.score) %>%
#   pivot_wider(names_from=rna_name, values_from=z.score, id_cols=ipa_canonical) %>%
#   column_to_rownames("ipa_canonical")
# 
# df_top_heat[is.na(df_top_heat)] = 0
# 
# annot = df_top %>%
#   select(rna_name, mean_normalized_PI, well_mean_eccentricity, well_mean_solidity) %>%
#   distinct() %>%
#   column_to_rownames("rna_name") %>%
#   as.data.frame()
# 
# col_fun = colorRamp2(c(-6, 0, 2), c("blue", "white", "red"))
# 
# Heatmap(df_top_heat,
#         column_names_gp = grid::gpar(fontsize = 8),
#         row_names_gp = grid::gpar(fontsize = 6),
#         heatmap_legend_param = list(title_gp = gpar( fontsize = 8)),
#         top_annotation = HeatmapAnnotation(df = annot,
#                                            annotation_legend_param = list(title_gp = gpar( fontsize = 8)),
#                                            annotation_name_gp= gpar(fontsize = 8)),
#         col = col_fun) 
# 
# 
# # Phagosome formation vs. Neuroinflammation -----
# df_phago = df %>% 
#   filter(ipa_canonical %in% c("Fcgamma receptor (FCGR) dependent phagocytosis","Neuroinflammation Signaling Pathway")) %>%
#   select(rna_name, ipa_canonical, z.score) %>%
#   pivot_wider(id_cols = "rna_name", values_from="z.score", names_from = "ipa_canonical") 
# 
# colnames(df_phago) = make.names(colnames(df_phago))
# 
# colnames(df_phago)
# library(ggpubr)
# ggscatter(df_phago, x = "Fcgamma.receptor..FCGR..dependent.phagocytosis", y = "Neuroinflammation.Signaling.Pathway",
#           add = "reg.line",  
#           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#           conf.int = TRUE, # Add confidence interval
#           label = "rna_name",
#           font.label = 8,
#           repel=T
# )+ stat_cor(method = "pearson", label.x = -2, label.y = -1) + 
#   ylab('Neuroinflammation Signaling Pathway') + 
#   xlab('Fc-Gamma receptor dependent phagocytosis')
# 
# 
# # compare morph, function, ipa features -----
# 
# phenos=c('mean_normalized_PI','well_mean_area','well_mean_eccentricity','well_mean_solidity')
# path_to_select=c("Fcgamma receptor (FCGR) dependent phagocytosis","Neuroinflammation Signaling Pathway",'SRP-dependent cotranslational protein targeting to membrane',
#  "Macrophage Classical Activation Signaling Pathway", "Cell Cycle Checkpoints")
# 
# colnames(df_top)
# df_top %>% filter(grepl('Classical', ipa_canonical))
#      
# multi_feature = df_top %>% 
#   filter(ipa_canonical %in% path_to_select) %>%
#   select(rna_name, ipa_canonical, z.score, phenos)  %>%
#   pivot_wider(id_cols = c("rna_name", phenos), values_from="z.score", names_from = "ipa_canonical") 
# 
# multi_feature[is.na(multi_feature)] = 0
# 
# multi_feature = multi_feature %>%
#   column_to_rownames("rna_name") 
# 
# multi_feature_scale = scale(multi_feature)
# 
# pheatmap(multi_feature_scale, scale="none")


