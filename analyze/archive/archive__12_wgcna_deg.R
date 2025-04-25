LIB="/cluster/tufts/patralab/rbator01/R_libs/4.4.0/"
.libPaths(LIB)
library(openxlsx)
library(tidyverse)
library(flashClust)
library(matrixStats)
library(WGCNA)
allowWGCNAThreads()  
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')
save_dir = "rd1_rd2_analysis/wgcna"

# NOT WORKING

# active res ----
active = read.xlsx("rd1_rd2_analysis/de/active_compound_list_addscreen2_29jul24.xlsx") %>%
  filter(active_or_confirmed_2sd == 1)
head(active)

# run wgcna -----
dds = readRDS("rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.dds.rds")
colData(dds)

meta = colData(dds) 

meta_active = meta %>%
  as.data.frame() %>%
  mutate(drug = paste0(Drug_Treatment,"_", Concentration)) %>%
  filter(drug %in% active$treatment)

deg = fread("rd1_rd2_analysis/de/all_res_padj_0.1_abslfc_1_run2.csv") %>%
  dplyr::filter(padj < 0.01)
deg = unique(deg$gene)
length(deg)
dds = dds[deg,rownames(meta_active)]

input_mat = counts(dds, normalized=TRUE)
dim(input_mat)
head(input_mat)

# scale free topology
powers = seq(from = 1, to = 40, by = 2)

sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
# saveRDS(sft, "rd1_rd2_analysis/wgcna/sft.rds")
# sft = readRDS("rd1_rd2_analysis/wgcna/sft.rds")

dev.off()
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 10
# no block ----

adj= adjacency(input_mat,
               type = "unsigned", 
               power = picked_power);

adj[1:10,1:10]

TOM=TOMsimilarityFromExpr(input_mat,
                          networkType = "unsigned", 
                          TOMType = "unsigned");

TOM[1:10,1:10]

dissTOM=1-TOM
geneTree = flashClust(as.dist(dissTOM),method="average");

dev.off()
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize = 100;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  
                            distM = dissTOM,
                            method="tree",
                            minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, 
#                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

module_colors= setdiff(unique(dynamicColors), "grey")
module_genes = NULL
for (color in module_colors){
  module=data.frame(gene = gene.names[which(dynamicColors==color)], module_color = color)
  if(is.null(module_genes)){
    module_genes = module
  }else{
    module_genes = rbind(module_genes, module)
  }
}

head(module_genes)
write.xlsx(module_genes, "rd1_rd2_analysis/wgcna/module_genes.xlsx")

# module-module correlations
MEList = moduleEigengenes(input_mat, colors = dynamicColors)
MEs = MEList$eigengenes
MEs

# correlation with genotype ------

#Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, colors = dynamicColors)$eigengenes
MEs <- orderMEs(MEs0)
head(MEs)

# Reorder modules so similar modules are next to each other
module_order = names(MEs) %>% gsub("ME","", .)

# Add treatment names

drug = data.frame(sample = row.names(MEs)) %>%
  left_join(meta_active %>%
              as.data.frame() %>%
              rownames_to_column("sample") %>%
              select(sample, drug), by="sample") %>%
  left_join(active %>%
              select(treatment, mean_phago), by=c("drug" = "treatment")) %>%
  column_to_rownames("sample")

nsamples = length(drug$mean_phago)

moduleTraitCor = cor(MEs, drug$mean_phago, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nsamples)

names(MEs)
moduleTraitPvalue
data.frame(moduleTraitCor[,1])
data.frame(moduleTraitPvalue[,1])

write.xlsx(data.frame(moduleTraitCor[,1]),paste0("rd1_rd2_analysis/wgcna/module_cor.xlsx"))
write.xlsx(data.frame(moduleTraitPvalue[,1]),paste0("rd1_rd2_analysis/wgcna/module_p.xlsx"))


textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")",  sep = "")
textMatrix
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor, xLabels = "mean_phago", yLabels = module_order, 
               ySymbols = names(module_order), colorLabels = FALSE, colors = greenWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1), 
               main = paste("Module-trait relationships"))



