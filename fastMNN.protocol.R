#' ---
#' title: "check IT with different protocols"
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


if (file.exists(paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.umap.rds"))
}else{
  metas.list <- list()
  counts.list <- list()
  
  metas.list$H1 <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj=="H1")
  metas.list$H1 <-  metas.list$H1 %>% select(-EML) %>% left_join(readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.umap.rds")) %>% select(cell,EML,prolifSig),by="cell")
  counts.list$H1 <- readRDS(file=paste0("tmp_data/",TD,"/counts.filter.rds"))[,metas.list$H1$cell]
  
  metas.list$Veres_2019 <- readRDS(paste0("tmp_data/","GSE114412_Veres_2019","/Veres_2019.meta.filter.rds")) %>% filter(devTime=="week4") %>% mutate(subCT=ifelse(EML %in% c("acinar_like","ductal_like","early_exo","late_exo","mesenchymal"),EML,subCT))%>% mutate(EML=ifelse(EML %in% c("acinar_like","ductal_like","early_exo","late_exo","mesenchymal"),"exo",EML)) %>% mutate(EML=recode(EML,"repl"="prolif","sc_alpha"="alpha","sc_beta"="beta","sc_ec"="SCEC","sst"="delta"))
  counts.list$Veres_2019 <- readRDS(file=paste0("tmp_data/","GSE114412_Veres_2019","/Veres_2019.counts.filter.rds"))[,metas.list$Veres_2019$cell]
  
  metas.list$BalBoa_2022 <- readRDS(paste0("tmp_data/","BalBoa_GSE167880","/BalBoa_GSE167880.meta.filter.rds")) %>% filter(grepl("S7.d20",cell)) #%>% filter(!grepl("S7.d20",cell)) %>% filter(!grepl("S7.d25",cell))  %>% filter(!grepl("S7.d18",cell))
  counts.list$BalBoa_2022 <- readRDS(file=paste0("tmp_data/","BalBoa_GSE167880","/BalBoa_GSE167880.counts.filter.rds"))[,metas.list$BalBoa_2022$cell]
  
  metas.list$Augsor_2022 <- readRDS(paste0("tmp_data/","Augsor_GSM4567006","/Augsor_GSM4567006.meta.filter.rds")) 
  counts.list$Augsor_2022 <- readRDS(file=paste0("tmp_data/","Augsor_GSM4567006","/Augsor_GSM4567006.counts.filter.rds"))[,metas.list$Augsor_2022$cell]
  
  meta.filter <- metas.list %>% do.call("bind_rows",.)
  lapply(counts.list,nrow)
  ov.genes <- rownames(counts.list$H1) %>% intersect(rownames(counts.list$Veres_2019)) %>% intersect(rownames(counts.list$BalBoa_2022)) %>% intersect(rownames(counts.list$Augsor_2022))
  counts.filter <- counts.list$H1[ov.genes,] %>% cbind(counts.list$Veres_2019[ov.genes,]) %>% cbind(counts.list$BalBoa_2022[ov.genes,]) %>% cbind(counts.list$Augsor_2022[ov.genes,])
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  if (file.exists(paste0("tmp_data/",TD,"/DP.fastMNN.lognormExp.mBN.rds"))) {
    lognormExp.mBN <- readRDS(file=paste0("tmp_data/",TD,"/DP.fastMNN.lognormExp.mBN.rds"))
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
    
    mBN.sce.ob <- multiBatchNorm(sce.ob$H1,sce.ob$Veres_2019,sce.ob$BalBoa_2022,sce.ob$Augsor_2022)
    lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
    saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/DP.fastMNN.lognormExp.mBN.rds"))
    
  }
  sel.od <- c("BalBoa_2022","Augsor_2022","H1","Veres_2019")
  temp.M <- meta.filter
  temp.sel.expG <- rownames(lognormExp.mBN)
  
  data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)#%>% CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes) 
  
  data.merge@assays$RNA$data <- as.matrix(lognormExp.mBN[temp.sel.expG,rownames(data.merge@meta.data)])
  data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
  data.spt <- data.spt[sel.od]
  
  nGene=2500;pc=25;
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
  data.temp <- data.ob %>% FindClusters(reso=0.6,verbose=F)
  
  DimPlot(data.temp,label=T,split.by = "pj",group.by = "EML")+NoLegend()
  DimPlot(data.temp,label=T)+NoLegend()
  plot_grid(plotlist = FunFP_plot(data.temp,unlist(main.mk)))

  
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "pj",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% select(cell,SID:mt.perc,EML,prolifSig) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")  # 
  data.ob.umap <- data.ob.umap  %>% mutate(cluster_EML=recode(SC,'C0'="alpha",'C1'="early_beta",'C2'="late_beta",'C3'="polyhormonal","C4"="delta","C5"="SCEC",'C6'="alpha",'C7'="SCEC",'C8'="exo",'C9'="SCEC",'C10'="prolif",'C11'="alpha",'C12'="early_beta",'C13'="exo",'C14'="exo"))
  
  Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"cluster_EML"])
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend()
  FeaturePlot(data.temp,"prolifSig")
  
  saveRDS(data.ob,paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.umap.rds"))
  
}

#'
data.temp <- data.ob %>% FindClusters(reso=1,verbose=F)
Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"cluster_EML"])
DimPlot(data.temp,label=T)+NoAxes()+NoLegend()



