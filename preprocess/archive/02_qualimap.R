LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))

library("jsonlite")
library(TidyMultiqc)
library(tidyverse)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/analysis/qualimap/multiqc_data/')

meta = read.xlsx("/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/Drug-Seq_Rd1_Treatment+Barcode_format.xlsx")
head(meta)

multiqc_data_path = "multiqc_data.json"
TidyMultiqc::list_plots(multiqc_data_path)

df = TidyMultiqc::load_multiqc(multiqc_data_path)

df = df %>%
  rename(Plate_Well = metadata.sample_id) %>%
  left_join(meta, by=c("Plate_Well")) 

head(df)

ggbarplot(df, x = "sample_name", y = "general.reads_aligned", fill="sample_name") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  guides(fill="none")



read_json(multiqc_data_path)

df = TidyMultiqc::load_multiqc(multiqc_data_path, 
  sections = 'plot',
  plots = "qualimap_genomic_origin")

df_int = df %>%
  as.data.frame %>%
  tidyr::unnest(cols = plot.qualimap_genomic_origin) %>%
  column_to_rownames("metadata.sample_id")

df_sum = rowSums(df_int)

intergenic_drugseq = data.frame(type = "drugseq", 
                                intergenic_fraction = df_int$intergenic/df_sum)


# multiqc_data_path = "~/poltoraklab/rbator01/ripk1_may2023/analysis/qualimap/multiqc_data/multiqc_data.json"
# TidyMultiqc::list_plots(multiqc_data_path)
# 
# df = TidyMultiqc::load_multiqc(
#   multiqc_data_path, 
#   sections = 'plot',
#   plots = "qualimap_genomic_origin"
# )
# 
# df_rna = df %>%
#   as.data.frame %>%
#   tidyr::unnest(cols = plot.qualimap_genomic_origin) %>%
#   column_to_rownames("metadata.sample_id")
# 
# head(df_rna)
# df_sum = rowSums(df_rna)
# 
# intergenic_rna = data.frame(type = "bulk_rnaseq", intergenic_fraction = df_rna$intergenic/df_sum)
# 
# 
# to_plot = rbind(intergenic_drugseq, intergenic_rna)
# 
# 
# ggplot(to_plot , aes(x=type,y=intergenic_fraction, fill=type)) +
#   geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) + 
#   ggtitle('intergenic reads fraction') + 
#   
# 
