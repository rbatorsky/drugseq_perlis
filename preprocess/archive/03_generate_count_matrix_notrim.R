library(data.table)
library(Matrix)
library(openxlsx)
library(tidyverse)
# this contains a code chunk for generating the count matrix for rd1, rd2, rd3

# rd1 -------
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/')

# generate count matrix -----
matrix_dir <- "analysis/star/Solo.out/Gene/raw/"

# the non-dedup -----
f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "analysis/processed/umi.notrim.nondedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# the dedup -----
f <- file(paste0(matrix_dir, "umiDedup-1MM_All.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "analysis/processed/umi.notrim.1mm_all_dedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# combine with meta data -----
mat = read.table("analysis/processed/umi.notrim.1mm_all_dedup.counts.txt")
meta = read.xlsx("Drug-Seq_Rd1_Treatment+Barcode.xlsx")

# format meta ----
table(meta$Drug_Treatment, meta$Concentration)

meta = meta %>%
  mutate(Drug_Treatment = gsub('\\)','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\(','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\-','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\+','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('_$','',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('^_','',Drug_Treatment)) %>%
  mutate(sample_name = paste0(Plate_Well,"|",Drug_Treatment)) 

write.xlsx(meta, "Drug-Seq_Rd1_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name

head(meta)
table(meta$Drug_Treatment)

# format mat -----
length(colnames(mat))
colnames(mat)
head(meta)
head(mat) 

# change the column names from barcode to sample_id|treatment
convert_names = data.frame(Barcode = colnames(mat))
head(convert_names)

convert_names = convert_names %>%
  left_join(meta, by="Barcode")


head(convert_names)

colnames(mat) = convert_names$sample_name
head(mat)

read_per_sample = colSums(mat)

hist(read_per_sample)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(mat)
gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(gene_list, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

saveRDS(mat_w_genes, "analysis/processed/umi.notrim.1mm_all_dedup.counts.rds")


# combine with meta data -----
mat = read.table("analysis/processed/umi.notrim.nondedup.counts.txt")
meta = read.xlsx("Drug-Seq_Rd1_Treatment+Barcode.xlsx")

# format meta ----
table(meta$Drug_Treatment, meta$Concentration)

meta = meta %>%
  mutate(Drug_Treatment = gsub('\\)','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\(','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\-','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\+','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('_$','',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('^_','',Drug_Treatment)) %>%
  mutate(sample_name = paste0(Plate_Well,"|",Drug_Treatment)) 

write.xlsx(meta, "Drug-Seq_Rd1_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name

head(meta)
table(meta$Drug_Treatment)

# format mat -----
length(colnames(mat))
colnames(mat)
head(meta)
head(mat) 

# change the column names from barcode to sample_id|treatment
convert_names = data.frame(Barcode = colnames(mat))
head(convert_names)

convert_names = convert_names %>%
  left_join(meta, by="Barcode")


head(convert_names)

colnames(mat) = convert_names$sample_name
head(mat)

read_per_sample = colSums(mat)

hist(read_per_sample)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(mat)
gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(gene_list, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

saveRDS(mat_w_genes, "analysis/processed/umi.notrim.nondedup.counts.rds")



#  rd2 ------
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/')

# generate count matrix -----
matrix_dir <- "analysis/star_rd2/Solo.out/Gene/raw/"

# the non-dedup -----
f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "analysis/processed/rd2.umi.nondedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# the dedup -----
f <- file(paste0(matrix_dir, "umiDedup-1MM_All.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "analysis/processed/rd2.umi.1mm_all_dedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)



# combine with meta data -----
mat = read.table("analysis/processed/rd2.umi.1mm_all_dedup.counts.txt")
meta = read.xlsx("Drug-Seq_Rd2_Treatment+Barcode.xlsx")

# format meta ----
table(meta$Drug_Treatment, meta$Concentration)

meta = meta %>%
  mutate(Drug_Treatment = gsub('\\)','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\(','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\-','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\+','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('_$','',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('^_','',Drug_Treatment)) %>%
  mutate(sample_name = paste0(Plate_Well,"|",Drug_Treatment)) 

write.xlsx(meta, "Drug-Seq_Rd2_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name

head(meta)
table(meta$Drug_Treatment)

# format mat -----
length(colnames(mat))
colnames(mat)
head(meta)
head(mat) 

# change the column names from barcode to sample_id|treatment
convert_names = data.frame(Barcode = colnames(mat))
head(convert_names)

convert_names = convert_names %>%
  left_join(meta, by="Barcode")


head(convert_names)

colnames(mat) = convert_names$sample_name
head(mat)

# remove empty well for better visualization
rm = which(colnames(mat) %in% c("H12|Empty well", "G12|Empty well"))
rm

meta_rm = meta %>%
  dplyr::filter(!Drug_Treatment == "Empty well")
mat =mat[, -rm]

read_per_sample = colSums(mat)

hist(read_per_sample)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(mat)
gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(gene_list, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

saveRDS(mat_w_genes, "analysis/processed/rd2.umi.1mm_all_dedup.counts.noempty.rds")

# nodedup ----

mat = read.table("analysis/processed/rd2.umi.nondedup.counts.txt")
meta = read.xlsx("Drug-Seq_Rd2_Treatment+Barcode_format.xlsx")

# change the column names from barcode to sample_id|treatment
convert_names = data.frame(Barcode = colnames(mat))
convert_names = convert_names %>%
  left_join(meta, by="Barcode")
convert_names
colnames(mat) = convert_names$sample_name

colnames(mat)
# remove empty well for better visualization
rm = which(colnames(mat) %in% c("H12|Empty well", "G12|Empty well"))
rm

meta_rm = meta %>%
  dplyr::filter(!Drug_Treatment == "Empty well")
mat =mat[, -rm]

read_per_sample = colSums(mat)

hist(read_per_sample)

mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(gene_list, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

saveRDS(mat_w_genes, "analysis/processed/rd2.umi.nondedup.counts.noempty.rds")

#  rd3 ------
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/')

# generate count matrix -----
matrix_dir <- "analysis/star_rd3/Solo.out/Gene/raw/"

# the non-dedup -----
f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "analysis/processed/rd3.umi.nondedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# the dedup -----
f <- file(paste0(matrix_dir, "umiDedup-1MM_All.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "analysis/processed/rd3.umi.1mm_all_dedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)


# combine with meta data -----
mat = read.table("analysis/processed/rd3.umi.1mm_all_dedup.counts.txt")
meta = read.xlsx("Drug-Seq_Rd3_Treatment+Barcode.xlsx")

# format meta ----
table(meta$Drug_Treatment, meta$Concentration)

unique(meta$Drug_Treatment)
meta = meta %>%
  mutate(Drug_Treatment = gsub('\\)','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\(','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\-','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('\\+','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('__','_',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('_$','',Drug_Treatment))%>%
  mutate(Drug_Treatment = gsub('^_','',Drug_Treatment)) %>%
  mutate(sample_name = paste0(Plate_Well,"|",Drug_Treatment)) 

head(meta)
write.xlsx(meta, "Drug-Seq_Rd3_Treatment+Barcode_format.xlsx")
rownames(meta) = meta$sample_name

head(meta)
table(meta$Drug_Treatment)

# format mat -----
length(colnames(mat))
colnames(mat)
head(meta)
head(mat) 

# change the column names from barcode to sample_id|treatment
convert_names = data.frame(Barcode = colnames(mat))
head(convert_names)

convert_names = convert_names %>%
  left_join(meta, by="Barcode")


head(convert_names)

colnames(mat) = convert_names$sample_name
head(mat)

read_per_sample = colSums(mat)

hist(read_per_sample)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(mat)
gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(gene_list, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

saveRDS(mat_w_genes, "analysis/processed/rd3.umi.1mm_all_dedup.counts.rds")

# nodedup ----

mat = read.table("analysis/processed/rd3.umi.nondedup.counts.txt")
meta = read.xlsx("Drug-Seq_Rd3_Treatment+Barcode_format.xlsx")

# change the column names from barcode to sample_id|treatment
convert_names = data.frame(Barcode = colnames(mat))
convert_names = convert_names %>%
  left_join(meta, by="Barcode")
convert_names
colnames(mat) = convert_names$sample_name

colnames(mat)

read_per_sample = colSums(mat)

hist(read_per_sample)

mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(gene_list, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

saveRDS(mat_w_genes, "analysis/processed/rd3.umi.nondedup.counts.rds")
