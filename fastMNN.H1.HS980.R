#' ---
#' title: "fastMNN H1 and HS980"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library
#R4.0
rm(list=ls())
condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

#https://www.nature.com/articles/s41586-019-1168-5
suppressPackageStartupMessages({
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


if (file.exists(paste0("tmp_data/",TD,"/H1_HS980.trans.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/H1_HS980.trans.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/H1_HS980.trans.data.ob.umap.rds"))
}else{
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("H1","HS980_notrans_CM310"))
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))[,meta.filter$cell]
  
  #' update annotation
  meta.filter <- meta.filter %>% rows_update(readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.umap.rds")) %>% select(cell,EML),by="cell")
  meta.filter <- meta.filter %>% rows_update(readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.ob.umap.rds")) %>% filter(devTime %in% c("HS980_notrans_CM310")) %>% select(cell,EML),by="cell")
  
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  if (file.exists(paste0("tmp_data/",TD,"/H1_HS980.fastMNN.lognormExp.mBN.rds"))) {
    lognormExp.mBN <- readRDS(file=paste0("tmp_data/",TD,"/H1_HS980.fastMNN.lognormExp.mBN.rds"))
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
    
    mBN.sce.ob <- multiBatchNorm(sce.ob$HS980_notrans_CM310,sce.ob$H1)
    lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
    saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/H1_HS980.fastMNN.lognormExp.mBN.rds"))
    
  }
  sel.od <- c("HS980_notrans_CM310","H1")
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
  rm(counts.filter)
  rm(sce.ob)
  rm(mBN.sce.ob)
  rm(lognormExp.mBN)
  
  data.spt <- data.spt[sel.od]
  set.seed(123)
  data.ob <- RunFastMNN(data.spt,verbose=F,features=length(mnn.VGs)) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc)#nn.method="annoy",annoy.metric="cosine") 
  data.temp <- data.ob %>% FindClusters(reso=1,verbose=F)
  
  data.temp <- data.ob %>%  FindClusters(resolution = 0.4,verbose = FALSE) 
  plot_grid(
    DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
    DimPlot(data.temp,label=T,group.by = "EML")+NoAxes()+NoLegend(),
    #FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
    #FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
    DimPlot(data.temp,label=T,group.by = "pj")+NoAxes()+NoLegend()
  )
  
  DimPlot(data.temp,label=T,group.by = "EML",split.by = "pj")+NoAxes()+NoLegend()
  plot_grid(plotlist = FunFP_plot(data.temp,unlist(main.mk)))
  
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% select(cell,SID:mt.perc) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")  # 
  #data.ob.umap <- data.ob.umap  %>% mutate(EML=recode(SC,'C0'="early_beta",'C1'="early_beta",'C2'="late_beta",'C3'="SCEC","C4"="alpha","C5"="polyhormonal",'C6'="psc",'C7'="prolif",'C8'="delta"))
  
  Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"])
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend()
  
  saveRDS(data.ob,paste0("tmp_data/",TD,"/H1_HS980.trans.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/H1_HS980.trans.data.ob.umap.rds"))
}

#' resolution
data.temp <- data.ob %>%  FindClusters(resolution = 0.4,verbose = FALSE) 
data.temp@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"]
plot_grid(
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
  FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
  DimPlot(data.temp,label=T,group.by = "EML")+NoAxes()+NoLegend()
)

