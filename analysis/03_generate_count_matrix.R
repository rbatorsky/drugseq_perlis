library(data.table)
library(Matrix)
library(openxlsx)
library(tidyverse)
library(biomaRt)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/')

# Function to process count matrices
process_count_matrix <- function(matrix_dir, file_prefix, output_file) {
  # Read count matrix (dedup and non-dedup)
  f <- file(paste0(matrix_dir, file_prefix, ".mtx"), "r")
  mat <- as.data.frame(as.matrix(readMM(f)))
  close(f)
  
  # Read feature and barcode names
  feature.names <- fread(paste0(matrix_dir, "features.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
  barcode.names <- fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
  
  # Assign row and column names
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V1
  
  # Save matrix to file
  fwrite(mat, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

# Function to format metadata
format_metadata <- function(meta) {
  # Clean up Drug_Treatment names and create a sample name
  meta <- meta %>%
    mutate(Drug_Treatment = gsub('[^[:alnum:]_]', '_', Drug_Treatment)) %>%
    mutate(Drug_Treatment = gsub('_+', '_', Drug_Treatment)) %>%
    mutate(Drug_Treatment = gsub('^_|_$', '', Drug_Treatment)) %>%
    mutate(sample_name = paste0(Plate_Well, "|", Drug_Treatment))
  
  rownames(meta) <- meta$sample_name
  return(meta)
}

# Function to update count matrix column names
update_matrix_colnames <- function(mat, meta) {
  convert_names <- data.frame(Barcode = colnames(mat))
  convert_names <- left_join(convert_names, meta, by = "Barcode")
  colnames(mat) <- convert_names$sample_name
  return(mat)
}

# Function to remove empty wells
remove_empty_wells <- function(mat, meta, empty_wells = c("H12|Empty well", "G12|Empty well")) {
  rm <- which(colnames(mat) %in% empty_wells)
  meta <- meta %>% filter(!Drug_Treatment %in% "Empty well")
  mat <- mat[, -rm]
  return(list(mat = mat, meta = meta))
}

# Function to add gene annotations
add_gene_annotations <- function(mat, genes) {
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = genes, mart = mart)
  
  mat_w_genes <- mat
  mat_w_genes$ensembl_gene_id <- rownames(mat_w_genes)
  
  mat_w_genes <- mat_w_genes %>%
    left_join(gene_list, by = "ensembl_gene_id") %>%
    mutate(ens_gene = paste0(ensembl_gene_id, "|", hgnc_symbol)) %>%
    column_to_rownames("ens_gene") %>%
    select(-c(ensembl_gene_id, hgnc_symbol))
  
  return(mat_w_genes)
}

# Path for data
matrix_dir <- "analysis/star/Solo.out/Gene/raw/"
meta_file <- "Drug-Seq_Rd1_Treatment+Barcode.xlsx"

# Process non-dedup count matrix (Rd1)
process_count_matrix(matrix_dir, "umiDedup-NoDedup", "analysis/processed/umi.notrim.nondedup.counts.txt")

# Process dedup count matrix (Rd1)
process_count_matrix(matrix_dir, "umiDedup-1MM_All", "analysis/processed/umi.notrim.1mm_all_dedup.counts.txt")

# Read meta and format it
meta <- read.xlsx(meta_file)
meta <- format_metadata(meta)

# Save formatted metadata
write.xlsx(meta, "Drug-Seq_Rd1_Treatment+Barcode_format.xlsx")

# Load count matrix and update column names
mat <- read.table("analysis/processed/umi.notrim.1mm_all_dedup.counts.txt")
mat <- update_matrix_colnames(mat, meta)

# Remove empty wells for better visualization
result <- remove_empty_wells(mat, meta)
mat <- result$mat
meta <- result$meta

# Add gene annotations
mat_w_genes <- add_gene_annotations(mat, rownames(mat))

# Save processed matrix with gene annotations
saveRDS(mat_w_genes, "analysis/processed/umi.notrim.1mm_all_dedup.counts.rds")

# Process similar steps for Rd2 and Rd3 using the same functions...