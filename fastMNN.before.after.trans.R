#' ---
#' title: "fastMNN before and after trans"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library
#R4.3
rm(list=ls())
#condaENV <- "/home/chenzh/miniconda3/envs/R4.3"
#LBpath <- paste0(condaENV ,"/lib/R/library")
#.libPaths(LBpath)

#https://www.nature.com/articles/s41586-019-1168-5
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scran)
  library(batchelor)
  library(Seurat)
  library(SeuratWrappers)
  library(cowplot)
  #library(scuttle)
  #library(SeuratDisk)
  #library(uwot)
})


# working directory
DIR <- "~/My_project/sc_pan"
knitr::opts_knit$set(root.dir=DIR)
setwd(DIR)


#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")


suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="Oct_2023"

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter
options(digits = 4)
options(future.globals.maxSize= 3001289600)

load("tmp_data/gene.meta.Rdata",verbose=T)


if (file.exists(paste0("tmp_data/",TD,"/before.after.trans.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.ob.umap.rds"))
  
  HS980.before.fm.mk <- readRDS(paste0("tmp_data/",TD,"/HS980.before.trans.fm.mk.rds"))
  HS980.after.fm.mk <- readRDS(paste0("tmp_data/",TD,"/HS980.after.trans.fm.mk.rds"))
  
}else{
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("HS980_trans","HS980_notrans_CM310"))
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))[,meta.filter$cell]
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  if (file.exists(paste0("tmp_data/",TD,"/before.after.trans.fastMNN.lognormExp.mBN.rds"))) {
    lognormExp.mBN <- readRDS(file=paste0("tmp_data/",TD,"/before.after.trans.fastMNN.lognormExp.mBN.rds"))
  }else{
    expG.set <- list()
    for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) {
      temp.cell <- meta.filter %>% filter(pj==b) %>% pull(cell)
      expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
    }
    sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
    length(sel.expG )
    
    sce.ob <- list()
    for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) {
      print(b)
      temp.M <- meta.filter %>% filter(pj==b)
      temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
      sce.ob[[b]] <- temp.sce
    }
    
    mBN.sce.ob <- multiBatchNorm(sce.ob$HS980_notrans_CM310,sce.ob$HS980_trans)
    lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
    saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/before.after.trans.fastMNN.lognormExp.mBN.rds"))
    
  }
  sel.od <- c("HS980_notrans_CM310","HS980_trans")
  #sel.od <- c("HS980_trans","HS980_notrans_CM310")
  temp.M <- meta.filter
  temp.sel.expG <- rownames(lognormExp.mBN)
  
  data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)#%>% CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes) 
  
  data.merge@assays$RNA$data <- as.matrix(lognormExp.mBN[temp.sel.expG,rownames(data.merge@meta.data)])
  data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
  data.spt <- data.spt[sel.od]
  
  nGene=2000;pc=25;
  mnn.VGs <- SelectIntegrationFeatures(object.list = data.spt, nfeatures = nGene)
  for (b in names(data.spt)) {
    VariableFeatures(data.spt[[b]]) <- mnn.VGs
  }
  
  #' release memory
  rm(data.merge)
  #rm(counts.filter)
  rm(sce.ob)
  rm(mBN.sce.ob)
  #rm(lognormExp.mBN)
  
  data.spt <- data.spt[sel.od]
  set.seed(123)
  data.ob <- RunFastMNN(data.spt,verbose=F,features=length(mnn.VGs)) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc)#nn.method="annoy",annoy.metric="cosine") 
  data.temp <- data.ob %>% FindClusters(reso=1,verbose=F)
  
  data.temp <- data.ob %>%  FindClusters(resolution = 0.6,verbose = FALSE) 
  plot_grid(
    DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
    FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
    FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
    DimPlot(data.temp,label=T,group.by = "pj")+NoAxes()+NoLegend()
  )
  plot_grid(plotlist = FunFP_plot(data.temp,unlist(main.mk)))
  
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% select(cell,SID:mt.perc) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")  # 
  data.ob.umap <- data.ob.umap  %>% mutate(EML=recode(SC,'C0'="late_beta",'C1'="late_beta",'C2'="early_beta",'C6'="early_beta",'C3'="SCEC","C4"="alpha","C5"="polyhormonal",'C8'="psc",'C9'="prolif",'C7'="delta"))
  
  Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"])
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend()
  
  data.temp <- JoinLayers(data.ob)
  data.temp <- AddModuleScore(data.temp, features = list(prolifSig=mature.score.mk), name = "matureSig")
  data.ob.umap$matureSig <- data.temp@meta.data[data.ob.umap$cell,"matureSig1"]
  
  saveRDS(data.ob,paste0("tmp_data/",TD,"/before.after.trans.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/before.after.trans.data.ob.umap.rds"))
  
  #' create the object for separate DEG detection
  # before
  temp.M <- meta.filter %>% filter(devTime=="HS980_notrans_CM310")
  temp.sel.expG <- rownames(lognormExp.mBN)
  
  data.before.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.before.ob@assays$RNA$data <- as.matrix(lognormExp.mBN[temp.sel.expG,rownames(data.before.ob@meta.data)])
  Idents(data.before.ob) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.before.ob@meta.data),"EML"])
  
  #after 
  temp.M <- meta.filter %>% filter(devTime=="HS980_trans")
  temp.sel.expG <- rownames(lognormExp.mBN)
  
  data.after.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.after.ob@assays$RNA$data <- as.matrix(lognormExp.mBN[temp.sel.expG,rownames(data.after.ob@meta.data)])
  Idents(data.after.ob) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.after.ob@meta.data),"EML"])
  
  
  HS980.before.fm.mk <- list()
  
  
  data.deg <- subset(data.before.ob,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.before.fm.mk$separate <- FindAllMarkers(data.deg) %>% tbl_df() %>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% mutate(cluster=recode(cluster,"early-beta"="early_beta","late-beta"="late_beta"))  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  data.deg <- subset(data.before.ob,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.deg <- RenameIdents(data.deg,"early_beta"="beta","late_beta"="beta")
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.before.fm.mk$merge <- FindAllMarkers(data.deg) %>% tbl_df() %>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene) %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  
  data.deg <- data.before.ob
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"EML"])
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.before.fm.mk$wp <- FindAllMarkers(data.deg) %>% tbl_df()%>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  data.deg <- subset(data.before.ob,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.deg <- RenameIdents(data.deg,"early_beta"="beta","late_beta"="beta")
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.before.fm.mk$merge <- FindAllMarkers(data.deg) %>% tbl_df()%>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  
  #' after
  HS980.after.fm.mk <- list()
  
  data.deg <- subset(data.after.ob,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.after.fm.mk$separate <- FindAllMarkers(data.deg) %>% tbl_df() %>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% mutate(cluster=recode(cluster,"early-beta"="early_beta","late-beta"="late_beta"))  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  data.deg <- subset(data.after.ob,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.deg <- RenameIdents(data.deg,"early_beta"="beta","late_beta"="beta")
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.after.fm.mk$merge <- FindAllMarkers(data.deg) %>% tbl_df() %>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene) %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  
  data.deg <- data.after.ob
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"EML"])
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.after.fm.mk$wp <- FindAllMarkers(data.deg) %>% tbl_df()%>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  data.deg <- subset(data.after.ob,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.deg <- RenameIdents(data.deg,"early_beta"="beta","late_beta"="beta")
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  HS980.after.fm.mk$merge <- FindAllMarkers(data.deg) %>% tbl_df()%>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  
  
  saveRDS(HS980.before.fm.mk,paste0("tmp_data/",TD,"/HS980.before.trans.fm.mk.rds"))
  saveRDS(HS980.after.fm.mk,paste0("tmp_data/",TD,"/HS980.after.trans.fm.mk.rds"))
  
  
}




#' resolution
data.temp <- data.ob %>%  FindClusters(resolution = 0.6,verbose = FALSE) 
data.temp@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"]
plot_grid(
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
  DimPlot(data.temp,label=T,group.by = "EML")+NoAxes()+NoLegend(),
  DimPlot(data.temp,label=T,group.by = "EML")+NoAxes()+NoLegend()
)


#' check markers
data.temp <- JoinLayers(data.ob)
temp.mk <- HS980BA.fm.mk$separate %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup()

temp.mk <- temp.mk %>% bind_rows(HS980BA.fm.mk$merge %>% filter(cluster=="beta") %>% filter(p_val_adj < 0.05) %>% filter(!gene %in% temp.mk$gene) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup())

temp.mk <- temp.mk %>% bind_rows(HS980BA.fm.mk$wp %>% filter(cluster=="polyhormonal") %>% filter(p_val_adj < 0.05) %>% filter(!gene %in% temp.mk$gene) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup())
temp.mk <-  (temp.mk %>% split(.,.$cluster))[c(mk.od,"psc")] %>% lapply(function(x){x %>% arrange(p_val_adj) %>% pull(gene)})


temp.M <- data.ob.umap %>% mutate(SID=paste(devTime,EML,sep=":")) %>% split(SID) %>% lapply(function(x){FunMaSF(200)}) %>% do.call("bind_rows",.) %>% select(cell,devTime,EML) %>% mutate(od=factor(EML,c(EML.od,"psc"),ordered = T)) %>% arrange(od) %>% select(-od)

temp.exp <- data.temp@assays$RNA$data[unlist(temp.mk) ,temp.M$cell]
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
zs.limit <- 2.5
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
temp.anno <- temp.M  %>% tibble::column_to_rownames("cell")
pheatmap::pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=4,border_color="NA",gaps_row =unlist(lapply( temp.mk,function(x){return(length(x))})) %>% cumsum(),main="",raster = TRUE) #%>%ggplotify::as.ggplot()



#‘ 
#library("reticulate")
#use_python('/home/chenzh/miniconda3/envs/R4.3/bin/python', require=T)
#' need to be run in terminal
check_terminal <- FALSE
if (check_terminal) {
  data.lend.temp <- data.ob %>%  FindClusters(resolution = 0.6,verbose = FALSE,algorithm = "4")
  saveRDS(data.lend.temp,paste0("tmp_data/",TD,"/before.after.trans.data.Leident.from.Terminal.rds"))
}
data.lend.temp <- readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.Leident.from.Terminal.rds"))
data.lend.temp <- RenameIdents(data.lend.temp,"1"="0","2"="1","3"="2","4"="3","5"="0","6"="4","7"="5","8"="6","9"="7","10"="8","11"="9")
data.temp <- data.ob %>%  FindClusters(resolution = 0.6,verbose = FALSE)




cowplot::plot_grid(
  DimPlot(data.temp,label=T) +NoLegend()+ggtitle("Louvain(Default)")+theme(plot.title = element_text(hjust=0.5)),
  DimPlot(data.lend.temp,label=T) +NoLegend()+ggtitle("Leiden")+theme(plot.title = element_text(hjust=0.5))
)



table(Idents(RenameIdents(data.temp,"0"="C0_1","1"="C0_1"))==Idents(RenameIdents(data.lend.temp,"0"="C0_1","1"="C0_1")))
