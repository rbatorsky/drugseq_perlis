LIB <- '/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("", LIB))

library(jsonlite)
library(TidyMultiqc)
library(tidyverse)
library(openxlsx)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/analysis/qualimap/multiqc_data/')

# Function to read and format metadata
load_metadata <- function(meta_file) {
  meta <- read.xlsx(meta_file)
  return(meta)
}

# Function to load and process multiqc data
load_multiqc_data <- function(multiqc_data_path, meta) {
  df <- TidyMultiqc::load_multiqc(multiqc_data_path)
  
  # Rename and join metadata
  df <- df %>%
    rename(Plate_Well = metadata.sample_id) %>%
    left_join(meta, by = c("Plate_Well"))
  
  return(df)
}

# Function to create bar plot
create_bar_plot <- function(df) {
  ggbarplot(df, x = "sample_name", y = "general.reads_aligned", fill = "sample_name") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
    guides(fill = "none")
}

# Function to extract intergenic fraction from multiqc data
extract_intergenic_fraction <- function(multiqc_data_path) {
  df <- TidyMultiqc::load_multiqc(multiqc_data_path, 
                                  sections = 'plot', 
                                  plots = "qualimap_genomic_origin")
  
  df_int <- df %>%
    as.data.frame %>%
    tidyr::unnest(cols = plot.qualimap_genomic_origin) %>%
    column_to_rownames("metadata.sample_id")
  
  # Sum of all values for each sample
  df_sum <- rowSums(df_int)
  
  # Calculate intergenic fraction
  intergenic_fraction <- data.frame(type = "drugseq", 
                                    intergenic_fraction = df_int$intergenic / df_sum)
  
  return(intergenic_fraction)
}

# Main execution
main <- function() {
  # Set file paths
  meta_file <- "/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/Drug-Seq_Rd1_Treatment+Barcode_format.xlsx"
  multiqc_data_path <- "multiqc_data.json"
  
  # Load metadata
  meta <- load_metadata(meta_file)
  
  # List available plots
  TidyMultiqc::list_plots(multiqc_data_path)
  
  # Load and merge multiqc data with metadata
  df <- load_multiqc_data(multiqc_data_path, meta)
  
  # Create and display the bar plot
  create_bar_plot(df)
  
  # Extract intergenic fraction
  intergenic_drugseq <- extract_intergenic_fraction(multiqc_data_path)
  
  # Optionally save the plot or data
  # ggsave("intergenic_reads_fraction.png")
  # write.xlsx(intergenic_drugseq, "intergenic_fraction.xlsx")
}

# Run the main function
main()