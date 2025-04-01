LIB="/cluster/tufts/patralab/rbator01/R_libs/4.4.0/"
.libPaths(LIB)
library(openxlsx)
library(tidyverse)
library(DESeq2)
select <- dplyr::select
renames <- dplyr::rename

library(flashClust)
library(matrixStats)

library(WGCNA)
allowWGCNAThreads()  

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# meta data -----
scr2 = read.xlsx("phenotypic_screen_data/2024.8.16_Secondary_screen_morphology_format.xlsx") 

scr2_2sd = scr2 %>%
  filter(raw_PI_Below_2SD_threshold == "yes")

# run wgcna -----
dds = readRDS("rd1_rd2_analysis/de/rd1_rd2_umi.nondedup.counts.dds.rds")
colData(dds)

# do we have all 2sd confirmed
meta = colData(dds) 

meta_2sd = meta %>%
  as.data.frame() %>%
  filter(!Concentration == "5uM") %>%
  mutate(drug = paste0(Drug_Treatment,"_", Concentration)) %>%
  inner_join(scr2_2sd, by=c("Drug_Treatment" = "Compound_Name")) 

meta_2sd[1:10,1:10]

rownames(meta_2sd) = meta_2sd$sample_name

numericMeta = meta_2sd

dds = dds[,rownames(meta_2sd)]

input_mat = counts(dds, normalized=TRUE)

topVarGenes <- head(order(rowVars(input_mat), decreasing = TRUE), 5000)
cleanDat  <- input_mat[ topVarGenes, ]
write.csv(cleanDat, "rd1_rd2_analysis/wgcna/expr_used_for_wgcna.csv")
dim(cleanDat)
gene.names=rownames(cleanDat)


# scale free topology ------
library("doParallel")
parallelThreads=24
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)

powers <- seq(4,28,by=1)  #initial power check -- try to get SFT.R.sq to go > 0.80
sft <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")

#plot initial SFT.R.sq vs. power curve
tableSFT<-sft[[2]]
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R?")


#choose power at elbow of SFT R? curve approaching asymptote near or ideally above 0.80
power=16 #SFT reached
minModSize=50
enforceMMS=TRUE

## Run an automated network analysis (ds=4 and mergeCutHeight=0.07, more liberal)
# choose parameters deepSplit and mergeCutHeight to get respectively more modules and more stringency sending more low connectivity genes to grey (not in modules).
net <- blockwiseModules(t(cleanDat),power=power,deepSplit=4,minModuleSize=minModSize,
                        mergeCutHeight=0.07,TOMdenom="mean", #detectCutHeight=0.9999,        #TOMdenom="mean" may get more small modules here.
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)       #maxBlockSize always more than the number of rows in cleanDat
#blockwiseModules can take 30 min+ for large numbers of gene products/proteins (10000s of rows); much quicker for smaller proteomic data sets

nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
modules
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

net.ds4.mms10<-net

#we will explore the blockwiseModules() function-built network with parameter deepSplit=4, minimum (initial) module size of 10
net<-net.ds4.mms10

# If necessary, return module members of small modules below size minSize=X to grey
if (enforceMMS) {
  removedModules<-orderedModules[which(modules<minModSize),"Color"]
  for(i in removedModules) { net$colors[net$colors==i] <- "grey" }
  for(i in removedModules) { net$MEs[,paste0("ME",i)] <- NULL }

  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
  as.data.frame(cbind(orderedModules,Size=modules))
}

#calculate kME table up front, in case we need to correct color assignments
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
net$MEs <- MEs
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")
table(net$colors)["grey"]
# 64
paste0(round(table(net$colors)["grey"]/nrow(cleanDat)*100,2),"% grey")
# 1.89% grey


##ITERATIVE until condition met that all module membes are at least 0.28 kMEintramodule.
#Go back and do final algorithm fix of module colors (remove kMEintramodule<0.28 members, reassign grey with kMEintramodule>0.35; max difference from kMEmax<0.10)

retry=TRUE;
kMEmaxDiff=0.1
reassignIfGT=0.30
greyIfLT=0.30
iter=1;
while (retry) {
  cat(paste0("\nkME table Cleanup, processing iteration ",iter,"..."))
  colorVecFixed<-colorVecBackup<-net$colors
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEintramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]) )  #all sig digits (no rounding), so max will be unique.
  colorVecFixed[kMEintramoduleVector<greyIfLT]<-"grey"
  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEmaxColorsVec<-apply( as.data.frame(cbind(kMEmaxVec,kMEdat)),1, function(x) gsub("kME","",colnames(kMEdat)[which(x==x[1])[2]-1]) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffTooBig<-(kMEmaxVec-kMEintramoduleVector) >= kMEmaxDiff
  colorVecFixed[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )] <- kMEmaxColorsVec[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )]
  net$colors<-colorVecFixed

#  table(net$colors)["grey"]  #decreased to x


# Are colors still in rank order? -- put them in order by recoloring modules that changed rank
  sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

  oldcolors <- names(sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"])
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==oldcolors[i]]<-paste0("proxy",labels2colors(i))
  }
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==paste0("proxy",labels2colors(i))]<-labels2colors(i)
  }

# one can check that colors are in order by size now
  #sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

# recalculate kME table, since we have corrected color assignments
  MEs<-tmpMEs<-data.frame()
  MEList = moduleEigengenes(t(cleanDat), colors = net$colors, verbose=0)
  MEs = orderMEs(MEList$eigengenes)
  net$MEs <- MEs
  colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
  rownames(MEs)<-rownames(numericMeta)

  tmpMEs <- MEs #net$MEs
  colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
  MEs[,"grey"] <- NULL
  tmpMEs[,"MEgrey"] <- NULL

  kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")

# recheck min kMEintramodule and max diff from kMEmax
  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEsIntramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { 1 } ) #grey proteins set to dummy value of 1 (ignore)

  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffCalc<- kMEmaxVec-kMEintramoduleVector
  if (min(kMEsIntramoduleVector)>=greyIfLT & max(kMEmaxDiffCalc)<=kMEmaxDiff) { cat(paste0("\nkME table 'clean' in ",iter," iterations.")); retry=FALSE; }
  iter=iter+1
  if (iter>30) break; #**
}
#** breaks after iteration 30 if did not reach criteria.

nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

#outputtabs<-outputfigs<-rootdir
FileBaseName="initialExo18samp_Net_toFindMEtoRegress_enforce_mms50"
########################################
#Write Module Membership/kME table
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDat)),rownames(cleanDat),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
write.table(kMEtable,file=paste0("rebecca_calc/2.ModuleAssignments-",FileBaseName,".txt"),sep="\t",row.names=FALSE)
#(load above file in excel and apply green-yellow-red conditional formatting heatmap to the columns with kME values); then save as excel with Tesco Exo multiple sheet data.

write.table(kMEtable,file=paste0("rebecca_calc/2.ModuleAssignments-",FileBaseName,".txt"),sep="\t",row.names=FALSE)
saveRDS(net, paste0("rebecca_calc/2.net-",FileBaseName,".rds"))