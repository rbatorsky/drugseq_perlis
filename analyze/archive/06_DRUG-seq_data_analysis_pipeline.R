library(stringr)
library(DT)
library(tidyr)
library(dplyr)
library(foreach)
library(limma)
library(edgeR)
library(umap)
library(data.table)
library(ggplot2)
library(GGally)
library(plotly)
library(ggrepel)
library(viridis)
library(RANN)
library(igraph)
library(RColorBrewer)
library(patchwork)

#functions----

get_union_toptable <- function(fit3) {
  contrast_names  <-  colnames(coef(fit3))
  n_contrasts  <-  length(contrast_names)
  
  tt <-  foreach(i=1:n_contrasts, .combine=rbind) %do%
  { topTable(fit3, coef = i, 
             p.value = 1, number = Inf) %>% 
      mutate(GeneId=rownames(.), contrast=colnames(coef(fit3))[i]) }
  
  tt$contrast  <-  factor(tt$contrast, levels=contrast_names)
  return(tt)
}

DE.limma_trend<-function(Exp=Exp,CTRL,meta,best_DMSO=NULL,random_run=FALSE){
  if(!is.null(best_DMSO)){ 
    print('Starting DE calculations using best DMSO.')
    }else{print('Starting DE calculations.')
    }
  
  de<-lapply(1:length(Exp), function(i){
    b<-names(Exp)[i]
    meta_sub<- meta %>% filter(batch==b) %>% 
      mutate(plate_barcode_short = substr(plate_barcode, nchar(plate_barcode)-2, nchar(plate_barcode)))
    if(!is.null(best_DMSO)) {
      meta_sub <- meta_sub %>% 
        filter(treatment!=CTRL|barcode_well %in% best_DMSO$barcode_well) #only use best DMSO as control
    } 

    rownames(meta_sub)<-meta_sub$plate_id
    UMI_batch<-lapply(1:length(Exp[[i]]), function(j){
      #print(paste0('i ',i,' / j ',j))
      x<-Exp[[i]][[j]]
      mat<-x$UMI.counts
      mat<-mat[!rownames(mat) %in% grep("ERCC-",rownames(mat),value = T),]
      mat<-mat[order(row.names(mat)),]
      return(mat)
      })
    # active checking that the objects in the list have exactly the same names
    unames <- unique(unlist(lapply(UMI_batch, rownames)))
    stopifnot( all( sapply(UMI_batch, function(x) identical(rownames(x), unames)) ) )
    UMI_batch<-do.call(cbind,UMI_batch)


    UMI_batch<-UMI_batch[,meta_sub$plate_id]
    clean<-UMI_batch[apply(UMI_batch,1,max)>quantile(apply(UMI_batch,1,max),0.75),] 
    clean<-clean[,colSums(clean)>0] 

    d<-DGEList(counts=clean)
    d <- calcNormFactors(d, method="TMM")
    clean.TMM<-log2(edgeR::cpm(d,  normalized.lib.sizes=T,log=F)+1) 

    meta_sub<-meta_sub[colnames(clean.TMM),]
    meta_sub$Sample<-sapply(meta_sub$Sample, function(k) gsub("-",".",k))



    design<-model.matrix(~0+factor(meta_sub$Sample)+factor(meta_sub$replicate))
    colnames(design)<-c(levels(factor(meta_sub$Sample)),levels(factor(meta_sub$replicate))[-1])
    rownames(design)<-colnames(clean)
    samples<-unique(meta_sub$Sample)
    ctrl<-grep(CTRL,samples,value = T)
    #print(ctrl)
    samples<-samples[!samples %in% ctrl]

    contrasts<-sapply(samples, function(ctrs)paste0(ctrs,"-",ctrl))
    contrast.matrix<-limma::makeContrasts(contrasts=contrasts, levels = design)

    fit = limma::lmFit(clean.TMM, design)
    fit2 = limma::contrasts.fit(fit,contrast.matrix)
    fit3 = limma::eBayes(fit2,  trend = T)

    DE.table <- get_union_toptable(fit3) %>% 
    separate(contrast,c('contrast','ctrl_arm'),sep="-") %>% 
    select(-ctrl_arm) %>% 
    merge(meta_sub %>% 
            filter(!treatment ==CTRL) %>% 
            select(plate_barcode_short,treatment,treatment_dose,timepoint,Sample,plate_well,batch) %>% 
            distinct() %>% 
            group_by(treatment,treatment_dose,timepoint,Sample,batch) %>% 
            summarise(plate_well = paste(paste(plate_barcode_short,plate_well,sep = " "), collapse=",")),
          sort=FALSE, by.x='contrast',by.y='Sample') %>% mutate(Sample=contrast) # combine metadata with toptable
    
  })
  if (random_run==TRUE){
    de<-do.call(rbind,de) %>% 
      filter(Sample=="RC_to_SA") %>% 
      group_by(batch,plate_well) %>%
      filter(abs(logFC)>1,adj.P.Val<0.1) %>%
      summarise(count=n()) %>%
      mutate(run_no=n,contrast='RC_to_SA') %>% 
      ungroup() 
    }else{
    de<-do.call(rbind,de) 
    }
  return(de)
}

vpPlot_DElabel <- function(this_contrast_toptabel,batch_select=1,this_arms){
  
  pcut<-0.1
  #color= black,red, blue
  #        0     1      2   
  color_custum<-c(rgb(0,0,0,0.5),rgb(1,0,0,1),rgb(0,0,1,1))
  
  this_contrast_toptabel <- this_contrast_toptabel %>% 
    filter(treatment %in% this_arms,batch==batch_select) %>% 
    mutate(label=0)
  
  

  high_concentration<-this_contrast_toptabel %>% select(treatment_dose) %>% max() 
  
  this_contrast_toptabel_plot<-NULL

  for (SampleName_lib_temp in this_arms) { 
    upregulate_gene <- this_contrast_toptabel %>% 
      filter(treatment == SampleName_lib_temp) %>% 
      filter(treatment_dose==high_concentration, adj.P.Val<pcut, logFC>1) %>% 
      select(GeneId) %>% 
      .$GeneId
    downregulate_gene <-this_contrast_toptabel %>% 
    filter(treatment == SampleName_lib_temp) %>%
    filter(treatment_dose==high_concentration, adj.P.Val<pcut, logFC< -1) %>% 
    select(GeneId) %>% 
    .$GeneId

    this_contrast_toptabel_plot <- rbind(this_contrast_toptabel_plot,this_contrast_toptabel 
                                      %>%  filter(treatment==SampleName_lib_temp) %>% 
                                        mutate(label=ifelse(GeneId %in% upregulate_gene,1,label)) %>% 
                                        mutate(label=ifelse(GeneId %in% downregulate_gene,2,label)))
  }
  
  this_contrast_toptabel_plot$label<-as.factor(this_contrast_toptabel_plot$label)
  levels(this_contrast_toptabel_plot$label)<-c("not significant ","upregulated genes","downregulated genes")
  this_contrast_toptabel_plot$treatment<-factor(this_contrast_toptabel_plot$treatment)

 
  
  vp<-this_contrast_toptabel_plot %>% arrange(label) %>% 
    ggplot(aes(x=logFC, y=-log10(adj.P.Val), color=label))+
    scale_color_manual(values = color_custum)+
    geom_point(fill="white",shape=21, size=1, stroke=0.5) +
    geom_abline(slope = 0, intercept = -log10(pcut), size=0.3, linetype='dashed') +
    geom_vline(xintercept = -1, size=0.3, linetype='dashed') +
    geom_vline(xintercept = 1, size=0.3, linetype='dashed') +
    xlim(-8,8)+
    theme_bw() +
    theme(axis.text.x = element_text(face="bold",size=14),
          axis.text.y = element_text(face="bold", size=10),
          axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12) ) +
    facet_grid(treatment~treatment_dose)+
    theme(strip.text = element_text(face="bold", size=12))

  return(vp)
  }

## works up to 3 groups comparison now
VD_plot <- function(this_contrast_toptable,batch_select=1){
  this_contrast_toptable_plot <- this_contrast_toptable %>% 
  mutate(sig=ifelse(abs(logFC)>1 & adj.P.Val<0.1,-1,0)) %>% 
    mutate(sig=ifelse(logFC>1 & adj.P.Val<0.1,1,sig)) %>%
    mutate(contrast=paste0(treatment_dose,'uM')) %>% 
    select(contrast,GeneId,sig) %>% 
    reshape2::dcast(GeneId~contrast,value.var="sig")
  
  
  rownames(this_contrast_toptable_plot) <- this_contrast_toptable_plot$GeneId
  return(this_contrast_toptable_plot %>% select(-GeneId) %>% vennDiagram(circle.col=c("turquoise", "salmon","palegreen"),cex=1))

}

pca_umap <- function(counts=NULL,center_pca=T,metric="euclidean", scale_pca=T,n_neighbors=15,n_pca=20){
  
  pca<-prcomp(t(counts), scale = scale_pca, center = center_pca)
  rownames(pca$x)  <-  colnames(counts)
  umap  <-  umap(pca$x[,1:n_pca],method = "umap-learn",metric=metric, n_neighbors=n_neighbors)
  rownames(umap$layout)  <-  colnames(counts)
  return(umap)

  
}

ggpairs_custom<-function(toptable_meta,this_arms,p_cut,corr_size=3,batch_s=1,strip_size=NULL,label_size=NULL,exclusive_gene=NULL){
  
  toptable_meta=toptable_meta %>% mutate(contrast=paste(treatment,treatment_dose,timepoint,sep="_"))
  this_contrast=toptable_meta %>% 
    filter(Sample %in% this_arms,batch==batch_s) %>% 
    arrange(treatment_dose) %>% 
    select(contrast) %>% 
    distinct() %>%.$contrast
  p_cut <- p_cut

  
if (is.null(exclusive_gene)) {
valid_genes <- toptable_meta %>% 
  filter(Sample %in% this_arms) %>%
  filter(batch==batch_s) %>%
  filter(adj.P.Val<p_cut) %>%
  .$GeneId %>%
  unique()
} else {
  
  valid_genes <- toptable_meta %>% 
    filter(Sample %in% this_arms) %>%
    filter(batch==batch_s) %>%
    filter(adj.P.Val<p_cut) %>%
  .$GeneId %>%
  unique() %>% 
    setdiff(exclusive_gene)
  
}

table_to_plot1<-toptable_meta %>% 
  filter(Sample %in% this_arms) %>%
  filter(batch==batch_s) %>%
  filter(GeneId %in% valid_genes) %>%
  dplyr::select(GeneId, contrast, logFC) %>%
  setDT %>%
  data.table::dcast(GeneId ~ contrast, value.var = 'logFC') %>%
  as.data.frame()

table_to_plot1<-table_to_plot1[,c("GeneId",this_contrast)]
rownames(table_to_plot1)<-table_to_plot1$GeneId

table_to_plot2<-toptable_meta %>% 
  filter(Sample %in% this_arms) %>%
  filter(batch==batch_s) %>%
  filter(adj.P.Val<p_cut) %>%
  dplyr::select(GeneId, contrast, logFC) %>%
  setDT %>%
  data.table::dcast(GeneId ~ contrast, value.var = 'logFC') %>%
  as.data.frame()

missing_col<-setdiff(colnames(table_to_plot1),colnames(table_to_plot2))
if (!identical(missing_col, character(0))) {
  missing_mat<-data.frame(matrix(0, nrow = nrow(table_to_plot2), ncol =as.numeric( length(missing_col))))
  missing_mat[missing_mat==0]=NA
  colnames(missing_mat)<-missing_col
  table_to_plot2<-cbind(table_to_plot2,missing_mat)
}

table_to_plot2<-table_to_plot2[,c("GeneId",this_contrast)]

data_plot<-inner_join(table_to_plot1,table_to_plot2,by="GeneId")
  
  

  rownames(data_plot)<-data_plot$GeneId
  data_plot<-data_plot %>% select(-"GeneId")
  
  axis_limt<-max(abs(data_plot),na.rm = TRUE)
  
  diagfun <- function(data,mapping){
    ggplot(data = data, mapping = mapping)+
      geom_density(color="darkblue", fill="lightblue") +
       scale_x_continuous(limits = c(-axis_limt, axis_limt)) 
      # scale_y_continuous(limits = c(-axis_limt, axis_limt))
    }
  
  corr_heatmap <- function(data, mapping, method="p", use="pairwise", ...){

              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

              # calculate correlation
              corr <- cor(x, y, method=method, use=use)
              colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
              fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
              ggally_cor(data = data, mapping = mapping, ...) + 
                theme_classic()+
                theme(panel.background = element_rect(fill=fill,colour = fill))
  }
  
  
  pm<-ggpairs(data_plot,c(1:as.numeric(length(colnames(data_plot))/2)),columnLabels=colnames(table_to_plot1)[-1],diag=list(continuous = wrap(diagfun)),
              upper = list(continuous = corr_heatmap),progress = F)+theme(axis.text = element_text(size = label_size),strip.placement = "outside", 
                                                                          text = element_text(size = strip_size) )
  
  method<-"pearson"
  pm2 <- pm
  color_custum<-c(rgb(0,0,0,0.5),rgb(0,1,0,0.7),rgb(0,0,1,0.7),rgb(1,0,0,1))
  label_list<-c('not significant','x axis DE genes','y axis DE genes','shared DE genes')
  for(i in 2:pm$nrow) {
       for(j in 1:(i-1)) {
         
         x <- GGally::eval_data_col(pm[i,j]$data, pm[i,j]$mapping$x)
         y <- GGally::eval_data_col(pm[i,j]$data, pm[i,j]$mapping$y)
         corr<- cor(x,y,method=method,use="na.or.complete")
         col_name_x_x<-paste0(gsub("\\.[^\\.]*$","",pm[i,j]$labels$x),".x")
         col_name_y_x<-paste0(gsub("\\.[^\\.]*$","",pm[i,j]$labels$y),".x")
         col_name_x_y<-paste0(gsub("\\.[^\\.]*$","",pm[i,j]$labels$x),".y")
         col_name_y_y<-paste0(gsub("\\.[^\\.]*$","",pm[i,j]$labels$y),".y")
         col_name<-c(col_name_x_x,col_name_y_x,col_name_x_y,col_name_y_y)
         data_sub<-pm$data[,col_name]
         
         data_sub<-data_sub %>% 
           mutate(label = ifelse(!is.na(!!sym(col_name_x_y)) & 
                                   abs(!!sym(col_name_x_y))>1,1,0)) %>%
          mutate(label = ifelse(!is.na(!!sym(col_name_y_y)) &
                                  abs(!!sym(col_name_y_y))>1,2,label)) %>%
          mutate(label = ifelse(!is.na(!!sym(col_name_x_y)) &
                        abs(!!sym(col_name_x_y))>1 &
                          !is.na(!!sym(col_name_y_y))&
                           abs(!!sym(col_name_y_y))>1 ,3,label)) %>% 
           arrange(label)
         
          unique_label<-sort(unique(data_sub$label))
          data_sub$label<-as.factor(data_sub$label)
          levels(data_sub$label)<-label_list[unique_label+1]
          
          p<-ggplot(data_sub,aes(x=!!sym(pm[i,j]$labels$x),y=!!sym(pm[i,j]$labels$y),color=label))+
               geom_point(size=2,shape=21,stroke=0.5,fill="white")+theme_bw()+scale_color_manual(values = color_custum[unique_label+1])+ scale_shape(solid = TRUE) +
            scale_x_continuous(limits = c(-axis_limt, axis_limt)) +
            scale_y_continuous(limits = c(-axis_limt, axis_limt)) +
            geom_abline(slope = 0, intercept = -1, size=0.3, linetype='dashed') +
            geom_abline(slope = 0, intercept =  1, size=0.3, linetype='dashed') +
            geom_vline(xintercept = -1, size=0.3, linetype='dashed') +
            geom_vline(xintercept = 1, size=0.3, linetype='dashed')
          pm2[i,j]=p+geom_label(data = data.frame(xlabel = -axis_limt,
                                          ylabel = max(axis_limt, na.rm = TRUE),
                                          lab = round(corr, digits = 3)),
                                mapping = ggplot2::aes(x = xlabel, y = ylabel, label = lab),
                                hjust = 0, vjust = 1,size = corr_size, fontface = "bold",
                          inherit.aes = FALSE # do not inherit anything from the ...
                          ) 
          
          
       }
   }
  return(pm2)
}


#UMI count table preparation -----
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/')

# input_data_dir <- file.path('20240530_Joshua_Bowen/')
# out_data_dir <- file.path('20240530_Joshua_Bowen/analysis/drugseq_pipeline/')
# 
# meta <- read.xlsx("20240530_Joshua_Bowen/Drug-Seq_Rd2_Treatment+Barcode_format.xlsx")
# colnames(meta)
# meta %>% dplyr::filter(Drug_Treatment == "DMSO")
# Exp = readRDS("20240530_Joshua_Bowen/analysis/star_rd2/Solo.out/Gene/processed/umi.1mm_all_dedup.counts.noempty.rds")
# head(Exp)
# Exp<-lapply(unique(meta$batch), function(x){
#   b<-unique(meta[meta$batch==x,]$plate_barcode)
#   res<-lapply(b, function(y){
#     load(file.path(input_data_dir, paste0("flowcell_4000_UMI_decode_",y,".RData")),verbose=T)
#     annot<-meta[meta$plate_barcode==y,]
#     rownames(annot) <- apply( annot[ , c('plate_barcode','well_index') ] , 1 , paste , collapse = "_" )
#     UMI_decode<-UMI_decode[!rownames(UMI_decode) %in% grep("ERCC-",rownames(UMI_decode),value = T),]
#     UMI_annot<-list(UMI_decode,annot)
#     names(UMI_annot)<-c("UMI.counts", "Annotation")
#     return(UMI_annot)
#   })
#   names(res)<-b
#   return(res)
# })
# 
# names(Exp)<-unique(meta$batch)
# 
# save(Exp, file=file.path(out_data_dir,"Exp_Init.RData"),compress = TRUE)
#
# 2 Step1
## 2.1 run first 500 random run
#
# not used
# step1_func <- function(n){
#   
#   CTRL="DMSO"
#   rand_runs=TRUE
#   
#   meta <- read.csv(file.path(input_data_dir,'meta.csv'),stringsAsFactors = F) %>%  
#     mutate(replicate=paste0('rep',plate_replicate),plate_id=paste0(plate_barcode,"_",well_index))
#   load(file.path(out_data_dir, 'Exp_Init.RData'),verbose = T)
#   
#   DMSO_sample <- meta %>% filter(treatment==CTRL) %>% select(batch,plate_barcode,plate_well) %>%
#     group_by(plate_barcode) %>% 
#     do(sample_n(.,1)) %>% 
#     mutate(Sample="RC_to_SA",treatment="RC_to_SA")
# 
# 
#   meta <- meta %>% 
#     left_join(DMSO_sample,by=c( 'batch','plate_barcode','plate_well'),suffix=c("","_randR")) %>%
#     mutate(Sample=ifelse(is.na(Sample_randR),Sample,Sample_randR)) %>% 
#     mutate(treatment=ifelse(is.na(treatment_randR),treatment,treatment_randR))
#   
#   DE<-DE.limma_trend(Exp=Exp, CTRL=CTRL, meta=meta, random_run =rand_runs)
# }
# DE_count_500 <- foreach(n=1:500,.packages=c("dplyr","limma","tidyr","foreach","edgeR"),.combine=rbind) %do% step1_func(n)
#
#
## 2.2 save and plot Step1 results
### 2.2.1 DE gene distribution (give it a new path to save)

pdf(file = file.path(out_data_dir,'first_500run.pdf') ,width =11 ,height =6 )
  H<-hist(DE_count_500$count,breaks=100, col="gray", labels = FALSE,ylim=c(0,900), main = 'Frequency of number DE genes per rand_run+batch (N=1500)_Step_1')
  text(x = H$mids, y = H$counts, labels = H$counts, cex = 0.75, pos=1,srt=90,offset = -1)
dev.off()

write.csv(DE_count_500,file=file.path(out_data_dir,'first_500run_DEcount.csv'))


### 2.2.2 DE counts statistical summary 

RC_to_AC_num_DE_seprow<-separate_rows(data = DE_count_500,'plate_well',sep = ',')
RC_to_AC_num_DE_seprow$batch_plate_well<-paste0(RC_to_AC_num_DE_seprow$batch,'_',RC_to_AC_num_DE_seprow$plate_well)
DE_summary1 <- data.frame(t(round(quantile(DE_count_500$count,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),stringsAsFactors = F,check.names = F) 
DE_summary1 %>% write.csv(file.path(out_data_dir,"summary_first_500run2.csv"),row.names = FALSE)
DE_summary1%>% datatable(caption = "First 500 random run quantile of numDEGs summary")


batches<-DE_count_500 %>% select(batch) %>% distinct() %>% arrange(batch) %>% .$batch
  quantiles_stats<-NULL
  for(h in 1:length(batches)) {
    quantiles_stats<-rbind(quantiles_stats,data.frame(batch_id=batches[h],t(round(quantile(subset(DE_count_500,DE_count_500$batch==batches[h],select='count')$count,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),stringsAsFactors = F,check.names = F))
  }
  
quantiles_stats %>%  datatable(caption = "First 500 random run quantile of numDEGs summary by batch")


## 2.3 select bad and best DMSO 

rsa_res<-NULL
source(file.path(analysis_dir,'RNAi-package.r'))
for(h in 1:length(batches)) {
    #print(h)
  batchid<-batches[h]
  RC_to_AC_num_DE_seprow$batch<-as.integer(RC_to_AC_num_DE_seprow$batch)
  RC_to_AC_num_DE_seprow_batch <- RC_to_AC_num_DE_seprow %>% 
    filter(batch==batchid) %>% 
    mutate(run_plate_well=paste0(batch_plate_well,"_",run_no)) %>% 
    mutate(log2FC=count/median(count,na.rm = T)) %>% 
    mutate(rz_score=(log2FC-median(log2FC,na.rm = T))/mad(log2FC,na.rm = T)) %>% 
    arrange(count)
  
  #RSA on ranked values, rz_score transformed
  rsa_res_batch<-score.screen.ZY.RQz(data = RC_to_AC_num_DE_seprow_batch,
                                     gene.column = 'batch_plate_well',
                                     value.column = 'rz_score',
                                     log.column = 'log2FC')                       
  rsa_res<-rbind(rsa_res,data.frame(batch_id=batches[h],rsa_res_batch,stringsAsFactors = F))

  
  
}

rsa_res[!is.finite(rsa_res$logP_RSA_Down),'logP_RSA_Down']<-0
rsa_res[!is.finite(rsa_res$logP_RSA_Up),'logP_RSA_Up']<-0

#prepare the plate_barcode and well
rsa_res<-separate(rsa_res,col = 'batch_plate_well',into =c('batch','rest'),sep = '_',remove = F)
rsa_res<-separate(rsa_res,col = 'rest',into =c('plate_barcode','plate_well'),sep = ' ',remove = T)
rsa_res$batch<-NULL
rsa_res$plate_barcode<-paste0('VH02001',rsa_res$plate_barcode)
rsa_res$sector<-1


for(i in 5:15){
rsa_res[,i]<-as.numeric(rsa_res[,i])
}
dim(rsa_res) 

dim(subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3))) 
rsa_res_rmv<-subset(rsa_res,!is.finite(rsa_res$logP_RSA_Up) | rsa_res$logP_RSA_Up<(-3))


rsa_res_rmv_summr<-rsa_res_rmv %>% group_by(batch_id,sector)%>% dplyr::summarize(num_DMSO_wells_rmv_per_batch=n())
rsa_res_rmv_summr

# best DMSO
for(i in 5:15){
  rsa_res[,i]<-as.numeric(rsa_res[,i])
}
dim(rsa_res) 


rsa_res<-rsa_res[order(rsa_res$logP_RSA_Down,decreasing = F),]
rsa_res<-rsa_res[order(rsa_res$plate_barcode),]
rsa_res_best_DMSO_wells<-NULL
# for each plate, select the 3 wells with best

uniq_plate_barcodes<-unique(rsa_res$plate_barcode)
for(h in 1: length(uniq_plate_barcodes)){
  rsa_res_1plate<-subset(rsa_res,rsa_res$plate_barcode==uniq_plate_barcodes[h])
  rsa_res_1plate$order <- seq(1,nrow(rsa_res_1plate),by=1)
    rsa_res_1plate<-rsa_res_1plate[1:3,] # 3 per plate
  rsa_res_best_DMSO_wells<-rbind(rsa_res_best_DMSO_wells,rsa_res_1plate)
}

write.csv(rsa_res_rmv,file = file.path(out_data_dir,'bad_DMSO_500remove.csv'),row.names = FALSE)
write.csv(rsa_res_best_DMSO_wells,file =file.path(out_data_dir,'best_DMSO_500keep.csv'),row.names = FALSE) 


#3 DE gene analysis using best DMSO 
meta <- read.csv(file.path(input_data_dir, 'meta.csv'),stringsAsFactors = F) %>%  
  mutate(replicate=paste0('rep',plate_replicate),plate_id=paste0(plate_barcode,"_",well_index),barcode_well=paste0(plate_barcode,"_",plate_well))
load(file.path(out_data_dir, 'Exp_Init.RData'),verbose = T)

CTRL<-'DMSO'
rsa_res_best_DMSO_wells <- read.csv(file.path(out_data_dir,'best_DMSO_500keep.csv'),stringsAsFactors = F) 
rsa_res_best_DMSO_wells1 <- rsa_res_best_DMSO_wells %>% mutate(barcode_well=paste0(plate_barcode,"_",plate_well)) 

DE<-DE.limma_trend(Exp=Exp,CTRL=CTRL,meta=meta,best_DMSO = rsa_res_best_DMSO_wells1,random_run = FALSE)
write.csv(DE,file.path(out_data_dir,"toptable_on_best_DMSO_500.csv"),row.names = FALSE)


## 3.2 data visulization volcano plot
treatment_to_plot <- c("QC-05-UB63","KA-73-NB69")
vp <- vpPlot_DElabel(DE,batch_select=2,this_arms=treatment_to_plot)
vp


## 3.3 Venn Diagram 
# works for up to 3 groups comparison
DE_selected  <- DE %>% filter(treatment=="QC-05-UB63",batch==1,treatment_dose %in% c(10.0000,3.16456,1.00144)) 
VD_plot(DE_selected)

## 3.4 pairplot
this_arms <- DE %>% 
  filter(treatment=="QC-05-UB63",batch==1) %>%
  arrange(treatment_dose) %>% 
  select(contrast) %>% 
  distinct() %>% 
  .$contrast

pm <- ggpairs_custom(DE,this_arms,0.1,3,batch_s=1,strip_size=8)
pm

this_arms <- DE %>% 
  filter(treatment=="KA-73-NB69",batch==1) %>%
  arrange(treatment_dose) %>% 
  select(contrast) %>% 
  distinct() %>% 
  .$contrast

pm <- ggpairs_custom(DE,this_arms,0.05,3,batch_s=1,strip_size=8)
pm



# 4 Step2 random run for UMAP plot 
## 4.1 generate SA from DMSO excluding bad DMSO and best DMSO
step2_func<-function(n){
  
  CTRL<-"DMSO"
  rand_runs<-TRUE
  
  meta<-read.csv(file.path(input_data_dir,'meta.csv'),stringsAsFactors = F) %>%  
    mutate(replicate=paste0('rep',plate_replicate),plate_id=paste0(plate_barcode,"_",well_index),barcode_well=paste0(plate_barcode,"_",plate_well))
  load(file.path(out_data_dir, 'Exp_Init.RData'),verbose = T)
  
  rsa_res_rmv<-read.csv(file.path(out_data_dir,'bad_DMSO_500remove.csv'),stringsAsFactors = F)
  rsa_res_best_DMSO_wells<-read.csv(file.path(out_data_dir,'best_DMSO_500keep.csv'),stringsAsFactors = F)


  rand_runs<-TRUE
  CTRL<-'DMSO'
  
  
  rsa_res_rmv1<-rsa_res_rmv %>% mutate(barcode_well=paste0(plate_barcode,"_",plate_well))
  rsa_res_best_DMSO_wells1<-rsa_res_best_DMSO_wells %>% mutate(barcode_well=paste0(plate_barcode,"_",plate_well)) 

  
  
    
  DMSO_sample<-meta %>% filter(treatment==CTRL,!barcode_well %in% c(rsa_res_rmv1$barcode_well,rsa_res_best_DMSO_wells1$barcode_well)) %>% 
    group_by(plate_barcode) %>% 
    do(sample_n(.,1)) %>% 
    mutate(Sample="RC_to_SA",treatment="RC_to_SA") %>% 
    select(batch,plate_barcode,plate_well,Sample,treatment)
  
  
  meta<-meta %>% left_join(DMSO_sample,by=c( 'batch','plate_barcode','plate_well'),suffix=c("","_randR")) %>%
    mutate(Sample=ifelse(is.na(Sample_randR),Sample,Sample_randR)) %>% 
    mutate(treatment=ifelse(is.na(treatment_randR),treatment,treatment_randR))
  
  DE<-DE.limma_trend(Exp=Exp,CTRL=CTRL,meta=meta, best_DMSO=rsa_res_best_DMSO_wells1,random_run =rand_runs)
}
DE_count_500_2<- foreach(n=1:500,.packages=c("dplyr","limma","tidyr","foreach","edgeR"),.combine=rbind) %do% step2_func(n)




## 4.2 save and plot Step4 results
### 4.2.1 DE gene distribution 
pdf(file<-file.path(out_data_dir,'second_500run.pdf') ,width =11 ,height =6 )
H<-hist(DE_count_500_2$count,breaks=100, col="gray", labels = FALSE,ylim=c(0,500), main = 'Frequency of number DE genes per rand_run+batch (N1500)_Step2')
text(x = H$mids, y = H$counts, labels = H$counts, cex = 0.75, pos=1,srt=90,offset = -1)
dev.off()

write.csv(DE_count_500_2,file=file.path(out_data_dir,"second_500run_DEcounts.csv"))

### 4.2.2 DE counts statistical summary
RC_to_AC_num_DE_seprow<-separate_rows(data = DE_count_500_2,'plate_well',sep = ',')
RC_to_AC_num_DE_seprow$batch_plate_well<-paste0(RC_to_AC_num_DE_seprow$batch,'_',RC_to_AC_num_DE_seprow$plate_well)
quantiles_table<-data.frame(t(round(quantile(DE_count_500_2$count,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),stringsAsFactors = F,check.names = F) 
quantiles_table %>% write.csv(file = file.path(out_data_dir,'summary_second_500run.csv'),row.names = FALSE)
quantiles_table%>% datatable(caption = "Second 500 random run quantile of numDEGs summary") 


batches<-DE_count_500_2 %>% select(batch) %>% distinct() %>% arrange(batch) %>% .$batch
  quantiles_stats<-NULL
  for(h in 1:length(batches)) {
    quantiles_stats<-rbind(quantiles_stats,data.frame(batch_id=batches[h],t(round(quantile(subset(DE_count_500_2,DE_count_500_2$batch==batches[h],select='count')$count,probs = c(0,0.01,0.05,0.1,0.5,0.9,0.95,0.99,1)))),stringsAsFactors = F,check.names = F))
  }
datatable(quantiles_stats,caption = "Second 500 run quantile of numDEGs summary by batch")



# 4.3 UMAP ----
### 4.3.1 active compounds and gene list filter
meta<-read.csv(file.path(input_data_dir, 'meta.csv'),stringsAsFactors = F) %>%
  mutate(barcode_well=paste0(plate_barcode,"_",plate_well),plate_id=paste0(plate_barcode,"_",well_index),batch_Sample=paste(batch,Sample,sep = "_"))
threshold_numDEGs=round(quantile(DE_count_500_2$count,probs = 0.95))
log2FC<-1
adj_pval<-0.1

DE.results<-DE %>% filter(abs(logFC)>log2FC,adj.P.Val<adj_pval) %>% 
  group_by(batch,Sample) %>% 
  summarise(N_DEgenes=n()) %>% 
  separate(Sample,c("treatment","dose","timepoint"),sep="_") %>% 
  mutate(treatment=str_replace_all(treatment,"\\.","-"),Sample=paste(treatment,dose,timepoint,sep="_")) %>% 
  select(-c(treatment,dose,timepoint)) %>% 
  mutate(active=ifelse(N_DEgenes>threshold_numDEGs,1,0))  

keep_sample<-DE.results %>% filter(active==1) %>% mutate(batch_Sample=paste(batch,Sample,sep = "_")) %>% 
  select(batch_Sample) %>% .$batch_Sample
keep_sample_meta<-meta %>% filter(batch_Sample %in% keep_sample)

rsa_res_best_DMSO_wells=read.csv(file.path(out_data_dir,'best_DMSO_500keep.csv'),stringsAsFactors = F) %>% 
  mutate(barcode_well=paste0(plate_barcode,"_",plate_well)) %>% select(barcode_well) %>% .$barcode_well
keep_DMSO_meta<-meta %>% filter(barcode_well %in% rsa_res_best_DMSO_wells)

keep_meta<-rbind(keep_sample_meta,keep_DMSO_meta)

# gene list filter
keep_genes<-DE %>% 
  filter(abs(logFC)>log2FC,adj.P.Val<adj_pval) %>% 
  select(GeneId) %>% 
  distinct() %>% 
  .$GeneId

### 4.3.2 prepare count table for UMAP
ump_colname<-keep_meta %>% select(plate_id) %>% .$plate_id
load(file.path(out_data_dir, 'Exp_Init.RData'),verbose = T)
UMI_batch_umap<-lapply(1:length(Exp), function(i){
  b<-names(Exp)[i]
  UMI_batch<-lapply(1:length(Exp[[i]]), function(j){
      print(paste0('i ',i,' / j ',j))
      x<-Exp[[i]][[j]]
      mat<-x$UMI.counts
      mat<-mat[!rownames(mat) %in% grep("ERCC-",rownames(mat),value = T),]
      mat<-mat[order(row.names(mat)),]
      return(mat)
    })
UMI_batch<-do.call(cbind,UMI_batch)
return(UMI_batch)

})

unames <- unique(unlist(lapply(UMI_batch_umap, rownames)))
stopifnot( all( sapply(UMI_batch_umap, function(x) identical(rownames(x), unames)) ) )
UMI_batch_umap<-do.call(cbind,UMI_batch_umap)
UMI_batch_umap<-UMI_batch_umap[keep_genes,ump_colname]
d<-DGEList(counts=UMI_batch_umap)
d <- calcNormFactors(d, method="TMM")
mat_ump_wbatchE<-log2(edgeR::cpm(d,  normalized.lib.sizes=T,log=F)+1)    
    

### 4.3.3 UMAP without remove batch correction
set.seed(87654)
umap_10<-pca_umap(counts = mat_ump_wbatchE, center_pca=T,metric ="euclidean",scale_pca = T, n_neighbors=20)

rownames(keep_meta)<-keep_meta$plate_id
color<-keep_meta %>% select(Sample,plate_barcode,treatment,batch,treatment_dose)
rownames(color)=keep_meta$plate_id



### 4.3.3.1 UMAP labeled by plate barcode
colnames(umap_10$layout)<-c('umap1','umap2')
umap_wbatchE_layout<-umap_10$layout 
umap_wbatchE_plot<-cbind(umap_wbatchE_layout,color)
umap_wbatchE_plot$Sample<-factor(umap_wbatchE_plot$Sample)
umap_wbatchE_plot$plate_barcode<-factor(umap_wbatchE_plot$plate_barcode)
umap_wbatchE_plot$treatment<-factor(umap_wbatchE_plot$treatment)
umap_wbatchE_plot$batch<-factor(umap_wbatchE_plot$batch)


n_color<-umap_wbatchE_plot$treatment %>% unique() %>% length()
p<-umap_wbatchE_plot%>% plot_ly(x = ~umap1, y = ~umap2,color = ~treatment,type = "scatter",mode = 'markers', size = ~treatment_dose,fill = ~'',text = ~Sample,marker = list(opacity = 0.5,sizemode='area'),colors=colorRampPalette(brewer.pal(8, "Set2"))(n_color))
htmlwidgets::saveWidget(as_widget(p), file.path(out_data_dir,"wtRBE_500plotly.html"))

p<-umap_wbatchE_plot%>% plot_ly(x = ~umap1, y = ~umap2,color = ~batch,type = "scatter",mode = 'markers',marker = list(size = 5),text = ~Sample)
htmlwidgets::saveWidget(as_widget(p), file.path(out_data_dir,"wtRBE_batch_500plotly.html"))


###4.3.4 UMAP with remove batch effect
groups<-factor(keep_meta$Sample)
groups<-relevel(groups,unique(grep(CTRL,groups,value = T)))
seq_batch<-keep_meta$plate_barcode
design<-model.matrix(~groups)
mat_ump_RBE <- removeBatchEffect(mat_ump_wbatchE, batch=seq_batch, design=design)

set.seed(87654)
umap_RBE<-pca_umap(counts = mat_ump_RBE, center_pca=T,metric ="euclidean",scale_pca = T, n_neighbors=20)


colnames(umap_RBE$layout)<-c('umap1','umap2')
color<-keep_meta %>% select(Sample,plate_barcode,treatment,batch,treatment_dose,compound)
rownames(color)<-keep_meta$plate_id
umap_RBE_layout<-umap_RBE$layout #%>% unlist() %>% as.data.frame()
umap_RBE_plot<-cbind(umap_RBE_layout,color)
umap_RBE_plot$Sample<-factor(umap_RBE_plot$Sample)
umap_RBE_plot$plate_barcode<-factor(umap_RBE_plot$plate_barcode)
umap_RBE_plot$batch<-factor(umap_RBE_plot$batch)
umap_RBE_plot$label<- paste0(umap_RBE_plot$compound,'_',umap_RBE_plot$treatment_dose,'um')
umap_RBE_plot$label<-factor(umap_RBE_plot$label)


n_color<-umap_wbatchE_plot$treatment %>% unique() %>% length()
p<-umap_RBE_plot%>% plot_ly(x = ~umap1, y = ~umap2,color = ~compound,type = "scatter", mode = 'markers',fill = ~'', size = ~treatment_dose,text = ~label,marker = list(opacity = 0.5,sizemode='area'),colors=colorRampPalette(brewer.pal(8, "Set2"))(n_color))
htmlwidgets::saveWidget(as_widget(p), file.path(out_data_dir,"RBE_500plotly.html"))


p<-umap_RBE_plot%>% plot_ly(x = ~umap1, y = ~umap2,color = ~batch,type = "scatter",fill = ~'',mode = 'markers', marker = list(size = 5),text = ~label)

htmlwidgets::saveWidget(as_widget(p), file.path(out_data_dir,"RBE_batch_500plotly.html"))



## 4.4 Umap with louvain clustering

n_pc<-20
pca<-prcomp(t(mat_ump_RBE), scale = T, center = T)
knn.info <- RANN::nn2(pca$x[,1:n_pc], k=35)
n_matrix <- dim(pca$x)[1]
knn <- knn.info$nn.idx
adj <- matrix(0, n_matrix, n_matrix)
rownames(adj) <- colnames(adj) <- rownames(pca$x[,1:n_pc])
for(i in seq_len(nrow(pca$x[,1:n_pc]))) {
    adj[i,rownames(pca$x[,1:n_pc])[knn[i,]]] <- 1
}

g <- igraph::graph.adjacency(adj, mode="undirected")
g <- simplify(g) ## remove self loops

clusterlouvain <- cluster_louvain(g)

### 4.4.1labeled by clustering number
cl<-data.frame(membership=clusterlouvain$membership)
n_color<-cl$membership %>% unique() %>% length()
rownames(cl)<-rownames(pca$x)
umap_RBE_plot_cl<-cbind(umap_RBE_layout,color,cluster=cl$membership)
umap_RBE_plot_cl$plate_barcode<-factor(umap_RBE_plot_cl$plate_barcode)
umap_RBE_plot_cl$batch<-factor(umap_RBE_plot_cl$batch)
umap_RBE_plot_cl$cluster<-factor(umap_RBE_plot_cl$cluster)
umap_RBE_plot_cl%>%  ggplot(aes(x = umap1, y = umap2,color=cluster))+
  geom_point(aes(size=treatment_dose),alpha = 1/5)+
  scale_fill_viridis()+
  theme_bw()+
  ggtitle("RBE_louvain_clustering") 
ggsave(file.path(out_data_dir,"louvain_clustering_labeled_by_clusters_500run.pdf"))

### 4.4.2 labeled by treatment
cl<-data.frame(membership=clusterlouvain$membership)
n_color<-cl$membership %>% unique() %>% length()
rownames(cl)<-rownames(pca$x)
umap_RBE_plot_cl<-cbind(umap_RBE_layout,color,cluster=cl$membership)
umap_RBE_plot_cl$compound<-factor(umap_RBE_plot_cl$compound)
umap_RBE_plot_cl$plate_barcode<-factor(umap_RBE_plot_cl$plate_barcode)
umap_RBE_plot_cl$batch<-factor(umap_RBE_plot_cl$batch)
umap_RBE_plot_cl$cluster<-factor(umap_RBE_plot_cl$cluster)
umap_RBE_plot_cl%>%  ggplot(aes(x = umap1, y = umap2,color=compound))+
  geom_point(aes(size=treatment_dose),alpha = 1/5)+
  scale_fill_viridis()+
  theme_bw()+
  ggtitle("RBE_louvain_clustering") 
ggsave(file.path(out_data_dir,"louvain_clustering_labeled_by_compound_500run.pdf"))


#### 4.3.4.2 clustering matrix
cl<-data.frame(membership=clusterlouvain$membership)
rownames(cl)<-rownames(pca$x)
confusion_matrix<-cbind(color,cluster=cl$membership) %>% select(batch,compound,cluster) %>% group_by(batch,compound,cluster) %>% summarise(count=n())
wide_matrix<-reshape2::dcast(confusion_matrix, batch+compound ~ cluster, value.var="count")
wide_matrix[is.na(wide_matrix)] <- 0

batch1<-wide_matrix %>% filter(batch==1)
rownames(batch1)<-batch1$compound
batch1<-batch1 %>% select(-c(compound,batch))
batch1<-batch1/rowSums(batch1)

batch1<-batch1 %>% mutate(compound=rownames(.))
mvar<-colnames(batch1 %>% select(-compound))
batch1<-batch1%>% reshape2::melt(id.vars = "compound", measure.vars = mvar) 
g1=ggplot(data = batch1, aes(x=variable, y=compound, fill=value))+ geom_tile()+labs(x = "batch1 clusters")+scale_x_discrete(position = "top")+scale_fill_gradient2(low = "white", high = "blue", limit = c(0,1), space = "Lab",guide = FALSE)


batch2<-wide_matrix %>% filter(batch==2)
rownames(batch2)<-batch2$compound
batch2<-batch2 %>% select(-c(compound,batch))
batch2<-batch2/rowSums(batch2)

batch2<-batch2 %>% mutate(compound=rownames(.))
mvar<-colnames(batch2 %>% select(-compound))
batch2<-batch2%>% reshape2::melt(id.vars = "compound", measure.vars = mvar)
g2=ggplot(data = batch2, aes(x=variable, y=compound, fill=value))+ geom_tile()+labs(x = "batch2 clusters")+scale_x_discrete(position = "top")+scale_fill_gradient2(low = "white", high = "blue", limit = c(0,1), space = "Lab",guide = FALSE)+theme(axis.title.y=element_blank())

batch3<-wide_matrix %>% filter(batch==3)
rownames(batch3)<-batch3$compound
batch3<-batch3 %>% select(-c(compound,batch))
batch3<-batch3/rowSums(batch3)

batch3<-batch3 %>% mutate(compound=rownames(.))
mvar<-colnames(batch3 %>% select(-compound))
batch3<-batch3%>% reshape2::melt(id.vars = "compound", measure.vars = mvar) 
g3=ggplot(data = batch3, aes(x=variable, y=compound, fill=value))+ geom_tile()+labs(x = "batch3 clusters")+scale_x_discrete(position = "top")+scale_fill_gradient2(low = "white", high = "blue", limit = c(0,1), space = "Lab", 
    name="clustering proportion")+theme(axis.title.y=element_blank())
g1+g2+g3 
ggsave(file.path(out_data_dir,"clustering_500run.pdf"))
