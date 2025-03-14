#' ---
#' title: "QC"
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

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scran)
  #library(batchelor)
  library(Seurat)
  #library(SeuratWrappers)
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
#source("src/local.quick.fun.R")


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

qc.nGene.min <- 1000
qc.nGene.max <- 7000
qc.mt.perc <- 0.25
qc.MH.nGene.FC <- 1/2.25
if (file.exists(paste0("tmp_data/",TD,"/meta.filter.rds"))) {
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  meta.all <- meta.all %>% mutate(MH.nGene.FC =mouse.nGene/nGene)
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
}else{
  #' loading data
  
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  load(paste0("tmp_data/","/gene.meta.Rdata"),verbose=T)
  
  meta.all <- meta.all %>% mutate(MH.nGene.FC =mouse.nGene/nGene)
  
  
  meta.filter <- meta.all %>% filter( MH.nGene.FC < qc.MH.nGene.FC & nGene < qc.nGene.max & mt.perc < qc.mt.perc & nGene > qc.nGene.min) %>% mutate_all(as.vector) 
  
  counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene),meta.filter$cell]
  
  
  #' get the expressed genes ( expressed in at least 1 datatsets)
  expG.set <- list()
  for (b in unique(meta.filter$devTime  %>% unique() %>% as.vector())) {
    temp.cell <- meta.filter %>% filter(devTime==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()

  sce.ob <- list()
  for (b in unique(meta.filter$devTime  %>% unique() %>% as.vector())) {
    print(b)
    temp.M <- meta.filter %>% filter(devTime==b)
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    sce.ob[[b]] <- temp.sce
  }

  mBN.sce.ob <- batchelor::multiBatchNorm(sce.ob$H1,sce.ob$HS980_notrans_CM310,sce.ob$HS980_notrans_CM311,sce.ob$HS980_trans)
  lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)

  saveRDS(meta.all,file=paste0("tmp_data/",TD,"/meta.all.rds"))
  saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/counts.filter.rds"))
  saveRDS(meta.filter %>% mutate(pj=EML),file=paste0("tmp_data/",TD,"/meta.filter.rds"))
  saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
}



#' #### check mouse.nGene and human.nGene distribution
#+ fig.width=9,fig.height=9
print(
  meta.all %>% ggplot()+geom_point(mapping=aes(x=nGene,y=mouse.nGene,color=EML),size=0.5)+geom_abline(slope=qc.MH.nGene.FC,intercept = 0,linetype="dashed")+facet_wrap(~EML)#+geom_abline(slope=1,intercept = 0,linetype="dashed")
)

print(
  meta.all %>% ggplot()+geom_histogram(mapping=aes(x=mouse.nGene/nGene,fill=EML),bins=100)+facet_wrap(~EML)+geom_vline(xintercept = qc.MH.nGene.FC,linetype="dashed")
)

#' #### check the general distribution
#+ fig.width=9,fig.height=9
print(
  meta.all%>% filter( MH.nGene.FC < qc.MH.nGene.FC ) %>% ggplot+geom_histogram(mapping=aes(x=nGene,fill=devTime),bins=100)+geom_vline(xintercept=qc.nGene.min)+geom_vline(xintercept=qc.nGene.max)+facet_wrap(.~devTime)+ylab("Number of cells")
)
#' #### check the mt.dis
#+ fig.width=9,fig.height=9
print(
  meta.all %>% filter( MH.nGene.FC < qc.MH.nGene.FC )%>% ggplot+geom_histogram(mapping=aes(x=mt.perc,fill=devTime),bins=100)+geom_vline(xintercept=qc.mt.perc)+facet_wrap(.~devTime)+ylab("Number of cells")
)
meta.all %>% filter( MH.nGene.FC < qc.MH.nGene.FC ) %>% ggplot+geom_point(mapping=aes(x=nGene,y=mt.perc,col=devTime))+facet_wrap(.~devTime)


#' #### number of cells
#+ fig.width=9,fig.height=9
print(
  meta.all %>% group_by(devTime) %>% summarise(BeforeQCnCellWithM=n()) %>% inner_join(meta.all%>% filter( MH.nGene.FC < qc.MH.nGene.FC )  %>% group_by(devTime) %>% summarise(BeforeQCnCellNoM=n()) )%>% inner_join(meta.filter  %>% group_by(devTime) %>% summarise(AfterQCnCell=n()),by="devTime") %>% gather(QC,nCell,-devTime) %>% ggplot+geom_bar(mapping=aes(x=devTime,fill=QC,y=nCell),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90))
)

#' check the cell number during QC
table(meta.all$devTime)
table(meta.all %>% filter(MH.nGene.FC < qc.MH.nGene.FC) %>% pull(devTime))
table(meta.filter$devTime)



