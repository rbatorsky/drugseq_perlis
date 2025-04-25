LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
select <- dplyr::select
rename <- dplyr::rename
library(openxlsx)
library(DESeq2)
library(data.table)
library('org.Hs.eg.db')
library(clusterProfiler)
library(ggrepel)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# drug seq meta ----
remove_from_meta = c('DMSO', 'ketamine', 'LSD', 'Empty well')
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx") %>%
  filter(!Drug_Treatment %in% remove_from_meta)

# format newest 2nd screen data ----
scr2 = read.csv("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table.csv")%>%
  mutate(Compound_Name = gsub('\\)','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\(','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\-','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\+','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('_$','',Compound_Name))%>%
  mutate(Compound_Name = gsub('^_','',Compound_Name))

#write.xlsx(scr2, "phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format.xlsx")

# calculate activity, ndeg calculated in 09_de_compounds_vs_best_dmso.R

scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format.xlsx")

determine_activity = read.csv("rd1_rd2_analysis/de/all_res_padj_0.05_abslfc_1_run3_ndeg.csv")

scr2_add_rna_activity = scr2 %>%
  left_join(determine_activity, by=c("Compound_Name" = "drug"))

old_table_transfer_info = read.xlsx("phenotypic_screen_data/2024.8.29_SecondaryScreen_CombinedData_format_add_rna_activity.xlsx") %>%
  select(Compound_Name, label,PI_below_50percent)

scr2_add_rna_activity = scr2_add_rna_activity %>% 
  full_join(old_table_transfer_info, by="Compound_Name") %>%
  mutate(PI_round = round(mean_normalized_PI,2)) %>%
  mutate(PI_below_50percent_calc = ifelse(PI_round <= 0.5,1,0))

#write.xlsx(scr2_add_rna_activity,   "phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format_add_rna_activity.xlsx")

# read back in ----
scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format_add_rna_activity.xlsx") %>%
  mutate(edited_name = gsub("_hydrochloride","_HCL",Compound_Name)) %>%
  mutate(edited_name = gsub("_acid","_Acid",edited_name)) %>%
  mutate(edited_name = gsub("_ditosylate","_Ditosylate",edited_name)) %>%
  mutate(edited_name = gsub("_tosylate","_Tosylate",edited_name)) %>%
  mutate(edited_name = gsub("JQ_1","(+)-JQ-1",edited_name)) %>%
  mutate(edited_name = gsub("_","",edited_name)) %>%
  mutate(name = edited_name) %>%
  column_to_rownames("edited_name") %>%
  rename(Eccentricity = well_mean_eccentricity) %>%
  rename(Solidity = well_mean_solidity)%>%
  mutate(`Phagocytic Index norm to DMSO` = mean_normalized_PI) %>%
  rename(`Mean IBA1 Intensity` = well_mean_normalized_IBA1_intensity) %>%
  mutate(ll = ifelse(PI_below_50percent_calc==1,name,"")) %>%
  mutate(label = ifelse(!is.na(label), name, label))%>%
  mutate(cc = ifelse(label == "DMSO", "darkgreen", 
                     ifelse(!is.na(label), 'blue', 'black')))%>%
  mutate(ss = ifelse(!is.na(label), 3, 2)) %>%
  mutate(sh = factor(ifelse(Compound_Name=="DMSO", 4,1)))%>% 
  filter(Compound_Name=="DMSO" | PI_below_50percent_calc==1) %>% 
  filter(Compound_Name!= "Lapatinib_ditosylate")

scr2

# summarize screen
table(scr2$PI_below_50percent_calc)

# do we have all 2sd confirmed in the ipa results
intersect(unique(scr2$Compound_Name), unique(meta$Drug_Treatment))
setdiff(unique(meta$Drug_Treatment), unique(scr2$Compound_Name))
setdiff(unique(scr2$Compound_Name), unique(meta$Drug_Treatment)) 
  
scr2_select = scr2 %>%
  select(c(Eccentricity, Solidity, `Phagocytic Index norm to DMSO`, `Mean IBA1 Intensity`))

scr2_select  <- sapply( scr2_select , as.numeric )
scr2_select
library(ggfortify)

pca = prcomp(scr2_select, scale=T)

library(ggpubr)
library(cowplot)

p1 = ggscatter(scr2, x = "Phagocytic Index norm to DMSO", y = "Eccentricity",
               shape="sh")+ 
  theme_bw() + 
  geom_point(color=scr2$cc, size=scr2$ss, shape=scr2$sh) + 
  geom_text_repel(label=scr2$label, color=scr2$cc) +
  theme(legend.position="none")

print(p1)
ggsave(p1, 
       filename = "rd1_rd2_analysis/paper_plots/4_b_pi50.pdf", 
       height=5,
       width=5.5)



p2 = ggscatter(scr2, x = "Mean IBA1 Intensity", y = "Eccentricity",
               shape="sh")+ 
  theme_bw() + 
  geom_point(color=scr2$cc, size=scr2$ss, shape=scr2$sh) + 
  geom_text_repel(label=scr2$label, color=scr2$cc) +
  theme(legend.position="none")

show(p2)
ggsave(p2, 
       filename = "rd1_rd2_analysis/paper_plots/4_c_pi50.pdf", 
       height=5,
       width=5.5)

p3 = autoplot(pca,
         data = scr2,
         shape="sh",
         loadings = TRUE, 
         loadings.colour = 'red', 
         loadings.label = T,  
         loading.label.color = 'red',
         loadings.label.repel=T) + 
  theme_bw() + 
  geom_text_repel(label=scr2$label, color=scr2$cc,nudge_x = .01) + 
  theme(legend.position="none")

show(p3)
ggsave(p3, 
       filename = "rd1_rd2_analysis/paper_plots/4_d_pi50.pdf", 
       height=5,
       width=5.5)


p4 = ggscatter(scr2, x = "Phagocytic Index norm to DMSO", y = "Solidity",
               shape="sh")+ 
  theme_bw() + 
  geom_point(color=scr2$cc, size=scr2$ss, shape=scr2$sh) + 
  geom_text_repel(label=scr2$label, color=scr2$cc) +
  theme(legend.position="none")

print(p4)
ggsave(p4, 
       filename = "rd1_rd2_analysis/paper_plots/SI_4a_pi50.pdf", 
       height=5,
       width=5.5)



p5 = ggscatter(scr2, x = "Mean IBA1 Intensity", y = "Solidity",
               shape="sh")+ 
  theme_bw() + 
  geom_point(color=scr2$cc, size=scr2$ss, shape=scr2$sh) + 
  geom_text_repel(label=scr2$label, color=scr2$cc)+
  theme(legend.position="none")

show(p5)
ggsave(p5, 
       filename = "rd1_rd2_analysis/paper_plots/SI_4b_pi50.pdf", 
       height=5,
       width=5.5)


#final = plot_grid(p1,p2,p3, ncol=3,
#          labels=c("B","C","D"))
#show(final)
#ggsave(final, 
#    filename = "rd1_rd2_analysis/paper_plots/4_bcd.pdf", 
#    height=5,
#    width=17)

# summarize the screen
table(scr2$PI_below_50percent)

view(meta)
# do we have all 2sd confirmed in the ipa results
intersect(unique(scr2$Compound_Name), unique(meta$Drug_Treatment))
setdiff(unique(meta$Drug_Treatment), unique(scr2$Compound_Name))
setdiff(unique(scr2$Compound_Name), unique(meta$Drug_Treatment))


# cluster the phenos

colnames(scr2)
scr2_select = scr2 %>%
  select("Compound_Name"                  ,"well_mean_area"             ,"well_mean_eccentricity"    , "well_mean_solidity"      ,   "mean_normalized_PI") %>%
  column_to_rownames("Compound_Name")

scr2_select_scale = t(scale(t(scr2_select)))
head(scr2_select_scale)
pca_res <- prcomp(scr2_select, scale=T)

autoplot(pca_res, 
         label=T,
         label.size = 3,
         label.repel=T,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         loadings.label.repel=T) 



pca_res.weights <- data.frame(pca_res$rotation)

pc1 <- pca_res.weights[order(pca_res.weights[, 1], decreasing = TRUE),][, 1, drop = FALSE]
pc2 <- pca_res.weights[order(pca_res.weights[, 2], decreasing = TRUE),][, 2, drop = FALSE]

# compare screen 1 and screen 2 -------

# all overlapping compounds
scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format.xlsx") 

scr1 = read.csv("phenotypic_screen_data/primary_screen_functional_means_corrected20240816.csv") %>%
  mutate(Compound_Name = gsub('\\)','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\(','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\-','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\+','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('_$','',Compound_Name))%>%
  mutate(Compound_Name = gsub('^_','',Compound_Name)) %>%
  select(Compound_Name, Compound_mean)

extra_data = data.frame(Compound_Name = c("ML218","GSK356278"), 
                        Compound_mean = c(0.190991556,0.266467147))

scr1 = rbind(scr1, extra_data)

setdiff(scr2$Compound_Name, scr1$Compound_Name)

# do the scores correlate scr1 to scr2?
colnames(scr1)
colnames(scr2)

scr1_scr2 = scr2 %>%
  select(Compound_Name, mean_normalized_PI) %>%
  rename(`Secondary screen mean normalized PI` = mean_normalized_PI) %>%
  inner_join(scr1 %>%
               select(Compound_Name, Compound_mean) %>%
               rename(`Primary screen mean normalized PI` = Compound_mean) , by="Compound_Name") %>%
  mutate(Compound_Name = gsub("_hydrochloride","_HCL",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_Hydrochloride","_HCL",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_acid","_Acid",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_ditosylate","_Ditosylate",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_tosylate","_Tosylate",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_","",Compound_Name))

library(ggpubr)
p = ggscatter(scr1_scr2, x = "Primary screen mean normalized PI", y = "Secondary screen mean normalized PI",
              add = "reg.line",  
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              label = "Compound_Name",
              font.label = 10,
              repel=T
)+ stat_cor(method = "pearson", label.x = 0.01, label.y = 1) + ylim(0,1) + xlim(0,0.5)

show(p)

ggsave(p, filename = "rd1_rd2_analysis/paper_plots/s1_s2_corr_all.pdf", height = 6, width = 6)

# confirmed only
scr2 = read.xlsx("phenotypic_screen_data/2024.9.24_secondary_results_supplement_table_format.xlsx") %>%
  filter(mean_normalized_PI <=0.55 )

scr1 = read.csv("phenotypic_screen_data/primary_screen_functional_means_corrected20240816.csv") %>%
  mutate(Compound_Name = gsub('\\)','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\(','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\-','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\+','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('_$','',Compound_Name))%>%
  mutate(Compound_Name = gsub('^_','',Compound_Name)) 

setdiff(scr2$Compound_Name, scr1$Compound_Name)
setdiff(scr2$Compound_Name, scr1$Compound_Name)

scr1$Compound_Name

# do the scores correlate scr1 to scr2?
colnames(scr1)
colnames(scr2)

scr1_scr2 = scr2 %>%
  select(Compound_Name, mean_normalized_PI) %>%
  rename(`Secondary screen mean normalized PI` = mean_normalized_PI) %>%
  inner_join(scr1 %>%
  select(Compound_Name, Compound_mean) %>%
  rename(`Primary screen mean normalized PI` = Compound_mean) , by="Compound_Name") %>%
  mutate(Compound_Name = gsub("_hydrochloride","_HCL",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_Hydrochloride","_HCL",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_acid","_Acid",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_ditosylate","_Ditosylate",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_tosylate","_Tosylate",Compound_Name)) %>%
  mutate(Compound_Name = gsub("_","",Compound_Name))


nrow(scr1_scr2)

library(ggpubr)
p = ggscatter(scr1_scr2, x = "Primary screen mean normalized PI", y = "Secondary screen mean normalized PI",
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          label = "Compound_Name",
          font.label = 10,
          repel=T
)+ stat_cor(method = "pearson", label.x = 0.05, label.y = 0.6) + ylim(0,0.7)

show(p)

ggsave(p, filename = "rd1_rd2_analysis/paper_plots/s1_s2_corr.pdf", height = 6, width = 6)
# plot only secondary 2sd confirmed -----



# archive below -----
rna = read.xlsx("rd1_rd2_analysis/de/active_compound_list_9jul24.xlsx")
view(rna)
head(scr2)
head(rna)

all_compounds = active %>%
  mutate(ndeg = ifelse(`0` == 0, `1`,`0`)) %>%
  mutate(active = ifelse(`1` == 0, 0,1)) %>%
  filter(!grepl('5uM', treatment)) %>%
  mutate(Compound_Name = gsub('_10uM','',treatment)) %>%
  mutate(Compound_Name = gsub('_2uM','',Compound_Name)) %>%
  rename(rna_name = treatment) %>%
  full_join(scr2, 
            by=c("Compound_Name")) %>%
  mutate(active_or_confirmed_2sd = ifelse(active == 1 | raw_PI_Below_2SD_threshold == "yes", 1, 0)) %>%
  mutate(active_and_confirmed_2sd = ifelse(active == 1 & raw_PI_Below_2SD_threshold == "yes", 1, 0))

view(all_compounds)

# which are we missing in the functional
all_compounds %>% filter(is.na(mean_normalized_PI))

#Clomipramine_hydrochloride_10uM
#MK_28_10uM
#ketamine (expected)
#lsd (expected)

# which are we missing in the rna
all_compounds %>% filter(is.na(ndeg))

#Lapatinib
#Zorifertinib

write.xlsx(all_compounds, "rd1_rd2_analysis/de/all_compounds_rnaactivity_functionalscreen_23aug24.xlsx")

confirmed_2sd = all_compounds %>%
  filter(raw_PI_Below_2SD_threshold == "yes")

view(confirmed_2sd)

active = all_compounds %>%
  filter(active == 1) 

# Archive below

# Compare Venn overlap archive ----
# overlap <- list(
#   scr2_1SD = confirmed_scr2_1sd$Compound_Name,
#   scr2_2SD = confirmed_scr2_2sd$Compound_Name ,
#   active = unique(active$Drug_Treatment)
# )
# 
# setdiff(unique(active_compounds$Drug_Treatment), confirmed_scr2_1sd$Compound_Name)
# 
# dpi=300
# #png(paste0("analysis/deg/venn_vic_padj_0.05_lfc_0.5deg.png"),width = dpi*5, height = dpi*5, units = "px",res = dpi, type="cairo")
# 
# library(ggvenn)
# p=ggvenn(
#   overlap,
#   stroke_size = 0.25, set_name_size = 4
# ) + ggtitle("Treatment Comparison")
# 
# print(p)
# #dev.off()
# 



# meta ----
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
# screen data ----
scr1 = read.csv("phenotypic_screen_data/primary_screen_functional_means.csv") %>%
  mutate(Compound_Name = gsub('\\)','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\(','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\-','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\+','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('_$','',Compound_Name))%>%
  mutate(Compound_Name = gsub('^_','',Compound_Name))  

scr1

scr2 = read.csv("phenotypic_screen_data/2024.8.13_Secondary_screen_morphology.csv") %>%
  mutate(Compound_Name = gsub('\\)','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\(','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\-','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('\\+','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('__','_',Compound_Name))%>%
  mutate(Compound_Name = gsub('_$','',Compound_Name))%>%
  mutate(Compound_Name = gsub('^_','',Compound_Name))  

view(scr2)

write.xlsx(scr2, "phenotypic_screen_data/2024.8.13_Secondary_screen_morphology_format.xlsx")

confirmed_scr2_2sd = scr2 %>%
  filter(raw_PI_Below_2SD_threshold == "yes") 

# do the scores correlate scr1 to scr2?
scr1_scr2 = scr1 %>%
  select(Compound_Name, Compound_mean) %>%
  rename(Compound_mean.s1 = Compound_mean) %>%
  full_join(scr2 %>%
              select(Compound_Name, mean_normalized_PI) %>%
              rename(mean_normalized_PI.s2 = mean_normalized_PI), by="Compound_Name")


library(ggpubr)
ggscatter(scr1_scr2, x = "Compound_mean.s1", y = "mean_normalized_PI.s2",
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          label = "Compound_Name",
          font.label = 8
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 2) + 
  ggtitle('correlation of phagocytosis scores')



# Does the ndeg correlate with the phago score? -----
tt = lm(mean_phago ~ ndeg, data=all_compounds) 
summary(tt)
# NO

# Analyze active and confirmed samples -----
active_or_confirmed = all_compounds %>%
  filter(active_or_confirmed_2sd == 1) %>%
  select(treatment, mean_phago)

head(active_or_confirmed)
# make the phagocytosis plots ------
deg = fread("rd1_rd2_analysis/de/all_res_run2.csv") 

# get top 200 DEG by lfc -----
deg_filter = deg %>% 
  dplyr::filter(!is.na(padj) & padj < 0.00001) %>%
  arrange(-abs(log2FoldChange)) %>%
  head(200)

head(deg_filter)
deg_select = deg %>%
  dplyr::filter(gene %in% deg_filter$gene) %>%
  arrange(factor(drug, levels=unique(deg$drug))) %>%
  dplyr::select("gene","log2FoldChange","padj", "drug")

lfc_select = deg_select %>%
  select("gene","log2FoldChange", "drug") %>%
  pivot_wider(id_cols="gene",names_from="drug", values_from="log2FoldChange", values_fill=0)

head(lfc_select)
rownames(lfc_select) = lfc_select$gene
lfc_select$gene = NULL
lfc_select_mat = as.matrix(lfc_select)
lfc_select_mat[is.na(lfc_select_mat)] <- 0

p_select = deg_select %>%
  select("gene","padj", "drug") %>%
  pivot_wider(id_cols="gene",names_from="drug", values_from="padj", values_fill=1)

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

annot = df_top %>%
  mutate(mean_phago = ifelse(treatment == "DMSO_10uM", NA, mean_phago)) %>%
  mutate(below_2SD_threshold = ifelse(treatment == "DMSO_10uM", NA, below_2SD_threshold)) %>%
  select(treatment, mean_phago, below_2SD_threshold) %>%
  distinct() %>%
  column_to_rownames("treatment") %>%
  as.data.frame()


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
             column_title = "top 200 DEG")

pdf("rd1_rd2_analysis/de/top_100_heatmap.pdf", height=16, width=7)
print(ph)
dev.off()

# HC and dendrogram

d <- dist(lfc_select_mat)
hc <- hclust(d)
plot(hc)
rect.hclust(hc, k = 2)

groups = cutree(hc, k = 2)
groups_df = data.frame(groups) %>%
  rownames_to_column("treatment")

write.xlsx(groups_df, "rd1_rd2_analysis/de/two_hc_groups.xlsx")

# Is the phago stuff different between groups?

active_or_confirmed_groups = active_or_confirmed %>%
  inner_join(groups_df, by="treatment") 

# ndeg
ggboxplot(active_or_confirmed_groups, x = "groups", y = "mean_phago",
          color = "groups", palette = "jco",
          add = "jitter") + stat_compare_means() + 
  ggtitle("screen2")


# What are these genes doing?

deg_select_group = deg_select %>% 
  left_join(groups_df, by=c("drug" = "treatment")) %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(gene, groups) %>%
  summarise_at(vars(log2FoldChange),
               list(Mean_log2FoldChange = mean)) %>%
  mutate(dir = ifelse(Mean_log2FoldChange > 0,"up","dn")) %>%
  separate(gene, into=c("ens","sym"), sep="\\|") 

head(deg_select_group)
ck<- compareCluster(geneCluster = ens~groups+dir,
                    data = deg_select_group,
                    OrgDb = org.Hs.eg.db,
                    keyType="ENSEMBL",
                    fun = "enrichGO",
                    ont="BP",
                    readable = T)


saveRDS(ck, "rd1_rd2_analysis/de/deg_hcgroups_dir_gobp.rds")
write.xlsx(ck, "rd1_rd2_analysis/de/deg_hcgroups_dir_gobp.xlsx")

ck = readRDS("rd1_rd2_analysis/de/deg_hcgroups_dir_gobp.rds")
dotplot(ck, show=30, by="Count")+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

ck_noion = ck %>%
  filter(!grepl(" ion", Description))


view(ck_noion)

dotplot(ck_noion, show=30, by="Count")+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

# Active AND confirmed, let's understand these first -------

active_2sd = all_compounds %>%
  filter(active == 1 & below_2SD_threshold == "yes" ) %>%
  select(treatment, mean_phago)

dmso_row = data.frame(treatment = "DMSO_10uM", mean_phago = NA)
active_2sd = rbind(active_2sd, dmso_row)

annot = active_2sd %>%
  column_to_rownames("treatment") %>%
  as.data.frame()

# do once only
#deg = fread("rd1_rd2_analysis/de/all_res_run2.csv") 
#deg = deg %>%
#  separate(gene, into=c("ens","sym"), sep="\\|")

deg_active_2sd = deg %>%
  filter(drug %in% active_2sd$treatment)

head(deg_active_2sd)

dmso = fread("rd1_rd2_analysis/de/dmso_500_batch_step2_run2_unique_genes.csv") %>%
  left_join(deg %>% 
              select(ens, sym), by="ens") %>%
  mutate(drug = "DMSO_10uM")

deg_active_2sd_dmso = rbind(deg_active_2sd[, colnames(dmso)], dmso) %>%
  distinct() %>%
  mutate(gene=paste0(sym,"|", ens)) 


# get top 10 DEG per sample -----
deg_filter = deg_active_2sd_dmso %>% 
  filter(!sym == "") %>%
  dplyr::filter(!is.na(padj) & padj < 0.01) %>%
  slice_max(order_by=abs(log2FoldChange), n=5, by=drug)

deg_select = deg_active_2sd_dmso %>%
  dplyr::filter(gene %in% deg_filter$gene) %>%
  arrange(factor(drug, levels=unique(deg$drug))) %>%
  dplyr::select("gene","sym","ens","log2FoldChange","padj", "drug")

lfc_select = deg_select %>%
  select("sym","log2FoldChange", "drug") %>%
  pivot_wider(id_cols="sym",names_from="drug", values_from="log2FoldChange", values_fill=0) 

lfc_select = lfc_select %>%
  as.data.frame()
rownames(lfc_select) = lfc_select$sym
lfc_select$sym = NULL
#lfc_select_mat = as.matrix(lfc_select)
lfc_select[is.na(lfc_select)] <- 0

dim(lfc_select)
p_select = deg_select %>%
  select("sym","padj", "drug") %>%
  pivot_wider(id_cols="sym",names_from="drug", values_from="padj", values_fill=1) 

rownames(p_select) = p_select$sym
p_select$sym = NULL
p_select_mat = as.matrix(p_select)
p_select_mat[is.na(p_select_mat)] <- 1

heatmap_max = max(deg_select$log2FoldChange)
heatmap_min = min(deg_select$log2FoldChange)
library(circlize)
library(ComplexHeatmap)


#lfc_select_mat = t(lfc_select_mat)
#p_select_mat = t(p_select_mat)
ph = Heatmap(lfc_select,
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
             column_title = "top 10 DEG per active + 2SD treatment",
             top_annotation = HeatmapAnnotation(df = annot))

#pdf("rd1_rd2_analysis/de/top_100_heatmap.pdf", height=16, width=7)
print(ph)
#dev.off()

# HC and dendrogram

d <- dist(t(lfc_select))
hc <- hclust(d)
plot(hc)
rect.hclust(hc, k = 2)

groups = cutree(hc, k = 2)
groups_df = data.frame(groups) %>%
  rownames_to_column("treatment")

write.xlsx(groups_df, "rd1_rd2_analysis/de/two_hc_groups_active_2sd.xlsx")

# Is the phago stuff different between groups?

active_2sd_groups = active_2sd %>%
  inner_join(groups_df, by="treatment") 

# ndeg
ggboxplot(active_2sd_groups, x = "groups", y = "mean_phago",
          color = "groups", palette = "jco",
          add = "jitter") + stat_compare_means() + 
  ggtitle("screen2")


# What are these genes doing?

deg_select_group = deg_select %>% 
  left_join(groups_df, by=c("drug" = "treatment")) %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(ens, groups) %>%
  summarise_at(vars(log2FoldChange),
               list(Mean_log2FoldChange = mean)) %>%
  mutate(dir = ifelse(Mean_log2FoldChange > 0,"up","dn")) 
head(deg_select)

deg_filter = deg_active_2sd_dmso %>% 
  dplyr::filter(!is.na(padj) & padj < 0.01) %>%
  slice_max(order_by=abs(log2FoldChange), n=50, by=drug)

head(deg_filter)

ck<- compareCluster(geneCluster = ens~drug,
                    data = deg_filter,
                    OrgDb = org.Hs.eg.db,
                    keyType="ENSEMBL",
                    fun = "enrichGO",
                    ont="BP",
                    readable = T)


saveRDS(ck, "rd1_rd2_analysis/de/deg_active_2sd_hcgroups_dir_gobp.rds")
write.xlsx(ck, "rd1_rd2_analysis/de/deg_active_2sd_hcgroups_dir_gobp.xlsx")

ck = readRDS("rd1_rd2_analysis/de/deg_active_2sd_hcgroups_dir_gobp.rds")
dotplot(ck, show=30, by="Count")+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))

ck_noion = ck %>%
  filter(!grepl(" ion", Description))


view(ck_noion)

dotplot(ck_noion, show=30, by="Count")+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1),
        axis.text.y = element_text(size=10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))


