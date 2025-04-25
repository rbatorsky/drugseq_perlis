LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))
library(tidyverse)
select <- dplyr::select
rename <- dplyr::rename
library(openxlsx)
library(DESeq2)
library(data.table)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# read data
mat = readRDS("rd1_rd2_analysis/rd1_rd2_umi.nondedup.counts.rds")
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")

table(meta$Drug_Treatment, meta$Concentration)

meta$drug = paste0(meta$Drug_Treatment, "_", meta$Concentration)
table(meta$drug)

rownames(meta) = meta$sample_name

# 1st time, based on 100 reps only
#best_wells = c("C08|DMSO|rd2","F05|DMSO|rd2","F11|DMSO|rd1","E09|DMSO|rd1")
#worst_wells = c("H01|DMSO|rd1","E05|DMSO|rd2")

#2nd time - use this
best_wells = c("F05|DMSO|rd2","C02|DMSO|rd2","H07|DMSO|rd1","C01|DMSO|rd1")
worst_wells=c("C07|DMSO|rd2","G01|DMSO|rd1")

dmso <- meta %>%
  dplyr::filter( Drug_Treatment == "DMSO" ) %>%
  dplyr::filter(sample_name %in% best_wells)

meta = meta %>%
  dplyr::filter(!(Drug_Treatment == "Empty well") & !(Drug_Treatment == 'DMSO')) 


# step 2, analyze compounds -----
all_res = NULL
for (d in unique(meta$drug)){
  
  active = meta %>%
    dplyr::filter(drug == d)
  
  # check that both batches are present
  if(length(unique(active$batch))>1){
    
    mat_i = mat[,c(active$sample_name, dmso$sample_name)]
    
    meta_i = meta %>%
      dplyr::filter(sample_name %in% c(active$sample_name))
    
    meta_i = rbind(meta_i, dmso)
    
    mat_i = mat_i[,rownames(meta_i)]
    
    
    dds <- DESeqDataSetFromMatrix(countData = mat_i, colData = meta_i, design = ~ batch + drug)
    
    #select_var= c("ENSG00000075624|ACTB", "ENSG00000184009|ACTG1")
    #plotCounts(dds, gene="ENSG00000075624|ACTB", intgroup="drug")
  }else{
    
    b = unique(active$batch)
    dmso_b = dmso %>%
      dplyr::filter(batch == b)
    
    mat_i = mat[,c(active$sample_name, dmso_b$sample_name)]
    
    meta_i = meta %>%
      dplyr::filter(sample_name %in% c(active$sample_name))
    
    meta_i = rbind(meta_i, dmso_b)
    
    mat_i = mat_i[,rownames(meta_i)]
    
    dds <- DESeqDataSetFromMatrix(countData = mat_i, colData = meta_i, design = ~ drug)
  }
  
  dds$drug <- relevel(dds$drug, "DMSO_10uM")
  
  dds <- DESeq(dds)
  name=paste0("drug_", d,"_vs_DMSO_10uM")
  #res <- results(dds)
  res <- lfcShrink(dds, coef=name, type="apeglm")
  
  res = data.frame(res)
  res$gene = rownames(res)
  rownames(res) = NULL
  
  write.csv(res, paste0("rd1_rd2_analysis/de/individual_drug/",d,"_run3.csv"))
  
  res_sig = res %>%
    dplyr::filter(padj<0.1, abs(log2FoldChange) > 1)
  
  if(nrow(res_sig) > 0){
    res_sig$drug = d
    
    if(is.null(all_res)){
      all_res = res_sig
    }else{
      all_res = rbind(all_res, res_sig)
    }
  }
}

# combine files -----
ff = list.files(path = "rd1_rd2_analysis/de/individual_drug/", pattern = "*_run3.csv")
ff
all_res=NULL
for (f in ff){
  d = gsub("_run4.notrim.csv","",f)
  res = fread(paste0("rd1_rd2_analysis/de/individual_drug/",f))%>%
    select(-V1) %>%
    filter(padj<0.1, abs(log2FoldChange) > 1) %>%
    mutate(drug = d)
  all_res = rbind(all_res, res)
}

write.csv(all_res, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run3.csv")

# no cutoff
ff = list.files(path = "rd1_rd2_analysis/de/individual_drug/", pattern = "*_run3.csv")
ff
all_res=NULL
for (f in ff){
  d = gsub("_run3.csv","",f)
  res = fread(paste0("rd1_rd2_analysis/de/individual_drug/",f)) %>%
    select(-V1) %>%
    mutate(drug = d)
  all_res = rbind(all_res, res)
}

write.csv(all_res, "rd1_rd2_analysis/de/all_res_run4.notrim.csv")

# look at the res ---
all_res = read.csv("rd1_rd2_analysis/de/all_res_run3.csv") %>%
  mutate(drug = gsub("_run4.notrim.csv","", drug)) %>%
  filter(!(grepl('LSD',drug))) %>%
  filter(!(grepl('ketamine',drug)))%>%
  filter(!(grepl('5uM',drug)))%>%
  mutate(drug = gsub("_10uM","",drug)) %>%
  mutate(drug = gsub("_2uM","",drug))

all_res_filter = all_res %>%
  mutate(sig = ifelse(padj < 0.1 & abs(log2FoldChange) > 1,1,0))

determine_activity = as.data.frame.matrix(table(all_res_filter$drug, all_res_filter$sig))  %>%
  rename(ndeg = `1`) %>%
  select(-`0`) %>%
  rownames_to_column("drug") %>%
  mutate(active = ifelse(ndeg > 179,1,0))

write.csv(determine_activity, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run3_ndeg.csv")

# investigate more stringent threshold ---- 

all_res = read.csv("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run3.csv") 

head(all_res)

all_res_filter = all_res %>%
  mutate(sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1,1,0))

determine_activity = as.data.frame.matrix(table(all_res_filter$drug, all_res_filter$sig))  %>%
  rename(ndeg = `1`) %>%
  select(-`0`) %>%
  rownames_to_column("drug") %>%
  mutate(active = ifelse(ndeg > 94,1,0))

table(determine_activity$active)

write.csv(determine_activity, "rd1_rd2_analysis/de/all_res_padj_0.05_abslfc_1_run3_ndeg.csv")


# determine rna activity - need to compare to functional (so this is a bit out of order) 
# do this here, remove from other places 10/8/24 -----

all_res = read.csv("rd1_rd2_analysis/de/all_res_run3.csv") %>%
  mutate(drug = gsub("_run4.notrim.csv","", drug)) %>%
  filter(!(grepl('LSD',drug))) %>%
  filter(!(grepl('ketamine',drug)))%>%
  filter(!(grepl('5uM',drug)))%>%
  mutate(drug = gsub("_10uM","",drug)) %>%
  mutate(drug = gsub("_2uM","",drug))

all_res_filter = all_res %>%
  mutate(sig = ifelse(padj < 0.1 & abs(log2FoldChange) > 1,1,0))

determine_activity = as.data.frame.matrix(table(all_res_filter$drug, all_res_filter$sig))  %>%
  rename(ndeg = `1`) %>%
  select(-`0`) %>%
  rownames_to_column("drug") %>%
  mutate(active = ifelse(ndeg > 179,1,0))

write.csv(determine_activity, "rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run3_ndeg_26oct24.csv")

head(determine_activity)
active_only = determine_activity %>%
  filter(active == 1) %>%
  .$drug

#compare=read.csv("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run3_ndeg.csv")
#compare = cbind(t,compare)
#view(compare)

#  IPA format 1 - take active or 2scr 2sd confirm only using screen results from july ----

all_compounds = read.xlsx("rd1_rd2_analysis/de/active_compound_list_addscreens_29jul24.xlsx")

active_2sd = all_compounds %>%
  filter(active == 1 & below_2SD_threshold == "yes")

active_2sd_filter = all_res %>%
  left_join(meta %>%
              mutate(drug = paste0(Drug_Treatment,"_", Concentration)) %>%
              select(drug, Drug_Treatment) %>%
              distinct(), by="drug") %>%
  filter(Drug_Treatment %in% active_2sd$Drug_Treatment) %>%
  dplyr::filter(padj<0.1)

unique(active_2sd_filter$drug)
max(min(abs(active_2sd_filter$log2FoldChange)))

length(unique(active_2sd$Drug_Treatment))

active_2sd_filter_wide = active_2sd_filter %>%
  separate(gene, into=c("ens","sym"), sep="\\|") %>%
  select(ens, log2FoldChange, padj, drug) %>%
  distinct() %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange),0,log2FoldChange)) %>%
  mutate(padj = ifelse(is.na(padj),1,padj)) %>%
  pivot_wider(id_cols = ens,
              names_from=drug,
              values_from = c("log2FoldChange", "padj"),
              values_fill = list(log2FoldChange = 0, padj = 1),
              names_glue = "{drug}_{.value}")


write.csv(active_2sd_filter_wide, "rd1_rd2_analysis/de/deg_active_and_confirmed_2sd_run3_padj_0.1_ipa.csv", row.names=F)

# IPA format 2 -- add more from the updated 2sd sheet + non active but 2sd -----

# old screen results
old = read.xlsx("phenotypic_screen_data/2024.8.13_Secondary_screen_morphology_format.xlsx")

old_2sd = old %>%
  filter(raw_PI_Below_2SD_threshold == "yes")

# new screen results
new = read.xlsx("phenotypic_screen_data/2024.8.16_Secondary_screen_morphology_format.xlsx")

new_2sd = new %>%
  filter(raw_PI_Below_2SD_threshold == "yes")

setdiff(new$Compound_Name,old$Compound_Name)
# [1] "Palbociclib_hydrochloride" "Lazertinib"                "SRT_2104"                  "JQ_1"                     
# [5] "Zorifertinib" 

setdiff(new_2sd$Compound_Name,old_2sd$Compound_Name)
# no new 2sd

# old screen + rna_activity results
all_compounds = read.xlsx("rd1_rd2_analysis/de/active_compound_list_addscreens_29jul24.xlsx")

# meta to link names from rna to screen
meta = read.xlsx("rd1_rd2_analysis/Drug-Seq_rd1_rd2_Treatment+Barcode_format.xlsx")
meta$drug = paste0(meta$Drug_Treatment, "_", meta$Concentration)

# all deg results
all_res = read.csv("rd1_rd2_analysis/de/all_res_run3.csv") 

all_res = all_res %>%
  mutate(sig = ifelse(padj < 0.1,1,0))

table(all_res$drug, all_res$sig)

# IPA FORMAT 3- these are the ones we put into ipa before (active and 2sd)
already_in_ipa = all_compounds %>%
  filter(active == 1 & below_2SD_threshold == "yes")

add_to_ipa = setdiff(new_2sd$Compound_Name, already_in_ipa$Drug_Treatment)
add_to_ipa

add_to_ipa_res = all_res %>%
  left_join(meta %>%
              mutate(drug = paste0(Drug_Treatment,"_", Concentration)) %>%
              select(drug, Drug_Treatment) %>%
              distinct(), by="drug") %>%
  filter(Drug_Treatment %in% add_to_ipa) %>%
  dplyr::filter(padj<0.1)

unique(add_to_ipa_res$drug)
max(min(abs(add_to_ipa_res$log2FoldChange)))

length(unique(add_to_ipa_res$Drug_Treatment))

add_to_ipa_res_wide = add_to_ipa_res %>%
  separate(gene, into=c("ens","sym"), sep="\\|") %>%
  select(ens, log2FoldChange, padj, drug) %>%
  distinct() %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange),0,log2FoldChange)) %>%
  mutate(padj = ifelse(is.na(padj),1,padj)) %>%
  pivot_wider(id_cols = ens,
              names_from=drug,
              values_from = c("log2FoldChange", "padj"),
              values_fill = list(log2FoldChange = 0, padj = 1),
              names_glue = "{drug}_{.value}")


write.csv(add_to_ipa_res_wide, "rd1_rd2_analysis/de/deg_notactive_and_confirmed_2sd_run3_padj_0.1_ipa.csv", row.names=F)

# IPA format 3- add the samples with new threshold -----

scr2 = read.xlsx("phenotypic_screen_data/2024.8.29_SecondaryScreen_CombinedData_format.xlsx") %>%
  filter(PI_below_50percent == "yes")

scr2$Compound_Name

alread_in_ipa_1 = read.csv("rd1_rd2_analysis/de/deg_notactive_and_confirmed_2sd_run3_padj_0.1_ipa.csv")
alread_in_ipa_2 = read.csv("rd1_rd2_analysis/de/deg_active_and_confirmed_2sd_run3_padj_0.1_ipa.csv")

already_in_ipa = c(colnames(alread_in_ipa_1), colnames(alread_in_ipa_2))
already_in_ipa = gsub("_padj","",already_in_ipa)
already_in_ipa = gsub("_log2FoldChange","",already_in_ipa)
already_in_ipa = gsub("_10uM","",already_in_ipa)
already_in_ipa = gsub("_5uM","",already_in_ipa)
already_in_ipa = unique(already_in_ipa)
already_in_ipa

setdiff(already_in_ipa, scr2$Compound_Name)

add_to_ipa = setdiff(scr2$Compound_Name,already_in_ipa)
add_to_ipa
all_res = read.csv("rd1_rd2_analysis/de/all_res_run3.csv")

all_res = all_res %>%
  mutate(sig_padj_0.1_abslfc_1 = ifelse(padj < 0.1 & abs(log2FoldChange)>1, 1,0))

add_to_ipa_res = all_res %>%
  left_join(meta %>%
              mutate(drug = paste0(Drug_Treatment,"_", Concentration)) %>%
              select(drug, Drug_Treatment) %>%
              distinct(), by="drug") %>%
  filter(Drug_Treatment %in% add_to_ipa) %>%
  dplyr::filter(padj<0.1)

unique(add_to_ipa_res$drug)
max(min(abs(add_to_ipa_res$log2FoldChange)))

length(unique(add_to_ipa_res$Drug_Treatment))

add_to_ipa_res_wide = add_to_ipa_res %>%
  separate(gene, into=c("ens","sym"), sep="\\|") %>%
  select(ens, log2FoldChange, padj, drug) %>%
  distinct() %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange),0,log2FoldChange)) %>%
  mutate(padj = ifelse(is.na(padj),1,padj)) %>%
  pivot_wider(id_cols = ens,
              names_from=drug,
              values_from = c("log2FoldChange", "padj"),
              values_fill = list(log2FoldChange = 0, padj = 1),
              names_glue = "{drug}_{.value}")


write.csv(add_to_ipa_res_wide, "rd1_rd2_analysis/de/deg_add_50percentdmsocutoff_run3_padj_0.1_ipa.csv", row.names=F)

add_to_ipa_res_wide = read.csv("rd1_rd2_analysis/de/deg_add_50percentdmsocutoff_run3_padj_0.1_ipa.csv")
colnames(add_to_ipa_res_wide)

# IPA format 4 - add the notrim reads -----
# I decided not to go back to the notrim reads
scr2_50_active = read.xlsx("phenotypic_screen_data/2024.8.29_SecondaryScreen_CombinedData_format_add_rna_activity.xlsx") %>%
  filter(PI_below_50percent == "yes") %>%
  filter(active == 1) %>%
  rename(drug = Compound_Name)

head(scr2_50_active)

all_res = read.csv("rd1_rd2_analysis/de/all_res_run4.notrim.csv") %>%
  mutate(drug = gsub("_run4.notrim.csv","", drug)) %>%
  filter(!(grepl('LSD',drug))) %>%
  filter(!(grepl('ketamine',drug)))%>%
  filter(!(grepl('5uM',drug)))%>%
  mutate(drug = gsub("_10uM","",drug)) %>%
  mutate(drug = gsub("_2uM","",drug))

# these are the res filtered down to 50 and active drug
res_50_active = all_res %>%
  inner_join(scr2_50_active, by="drug")

# these are the genes to consider, filter padj for simplicity
res_50_active_filter = res_50_active %>%
  filter(padj < 0.1 )

# then plot these genes in these drugs
res_50_active_filter_wide = res_50_active %>%
  filter(gene %in% res_50_active_filter$gene) %>%
  separate(gene, into=c("ens","sym"), sep="\\|") %>%
  select(ens, log2FoldChange, padj, drug) %>%
  distinct() %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange),0,log2FoldChange)) %>%
  mutate(padj = ifelse(is.na(padj),1,padj)) %>%
  pivot_wider(id_cols = ens,
              names_from=drug,
              values_from = c("log2FoldChange", "padj"),
              values_fill = list(log2FoldChange = 0, padj = 1),
              names_glue = "{drug}_{.value}")


write.csv(res_50_active_filter_wide, "rd1_rd2_analysis/de/deg_active_and_confirmed_2sd_run4_notrim_padj_0.1_ipa.csv", row.names=F)
