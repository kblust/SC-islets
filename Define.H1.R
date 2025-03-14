#' ---
#' title: "define the H1"
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
  library(cowplot)
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
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))


heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)

if (file.exists(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.rds"))) {
  data.ob <-  readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.rds"))
  data.ob.umap<- readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.umap.rds"))
}else{

  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  # H1 cells
  temp.M <- meta.filter %>% filter(pj =="H1") 
  set.seed(123)
  temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >=5]
  data.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)  
  
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
  temp.norm <- scuttle::normalizeCounts(temp.sce)
  data.ob@assays$RNA$data <- as.sparse(temp.norm[temp.sel.expG,colnames(data.ob)])
  
  data.ob <- data.ob %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:25,verbose=F) %>% FindNeighbors( dims = 1:25,verbose = FALSE) #  #vars.to.regress=c("mt.perc","S.Score", "G2M.Score")#,nn.method="annoy",annoy.metric="cosine"
  
  data.ob <- AddModuleScore(data.ob, features = list(prolifSig=prolif.mk), name = "prolifSig")
  data.temp <- data.ob %>%  FindClusters(resolution = 0.4,verbose = FALSE) 
  plot_grid(
    DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
    FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
    FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend()
    #DimPlot(data.temp,label=T,group.by = "Phase")+NoAxes()+NoLegend()
  )
  plot_grid(plotlist = FunFP_plot(data.temp,unlist(main.mk)))
  
  data.ob.umap <-  data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,EML:Phase,prolifSig1)) %>% rename(prolifSig=prolifSig1) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))%>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell")
  data.ob.umap <- data.ob.umap %>% mutate(EML=recode(seurat_clusters,'C0'="early_beta",'C1'="late_beta",'C2'="alpha",'C3'="early_beta",'C4'="SCEC",'C5'="early_beta",'C6'="polyhormonal",'C7'="prolif","C8"="delta"))
  
  
  Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"EML"])
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend()
  
  #H1.mk <- FunRF_FindAllMarkers_para(data.temp)
  H1.fm.mk <- list()
  
  data.deg <- subset(data.temp,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"EML"])
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  H1.fm.mk$separate <- FindAllMarkers(data.deg) %>% tbl_df() %>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% mutate(cluster=recode(cluster,"early-beta"="early_beta","late-beta"="late_beta"))  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))

  data.deg <- subset(data.temp,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.deg <- RenameIdents(data.deg,"early_beta"="beta","late_beta"="beta")
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  H1.fm.mk$merge <- FindAllMarkers(data.deg) %>% tbl_df() %>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene) %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  
  data.deg <- data.temp
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"EML"])
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  H1.fm.mk$wp <- FindAllMarkers(data.deg) %>% tbl_df()%>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
  
  
  
  
  
  data.deg <- subset(data.temp,cell=(data.ob.umap %>% filter(EML!="polyhormonal") %>% pull(cell)))
  data.deg <- RenameIdents(data.deg,"early_beta"="beta","late_beta"="beta")
  data.ave.exp <- log1p(AverageExpression(data.deg)$RNA) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df()
  H1.fm.mk$merge <- FindAllMarkers(data.deg) %>% tbl_df()%>% filter(p_val_adj <0.05 & avg_log2FC > 0.25)  %>% inner_join(data.ave.exp %>% gather(cluster,ave_exp,-gene)  %>% group_by(gene) %>% top_n(1,ave_exp),by = c("cluster", "gene"))
  
 
  
  saveRDS(data.ob,file=paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.rds"))
  saveRDS(data.ob.umap,file=paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.umap.rds"))
  saveRDS(H1.fm.mk,file=paste0("tmp_data/",TD,"/sc.pan.H1.fm.mk.rds"))

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
FeaturePlot(data.temp,"prolifSig1")
plot_grid(plotlist = c(
          FunFP_plot(data.temp,unlist(main.mk)),
          list(umap=(DimPlot(data.temp,label=T)+NoAxes()+NoLegend())))
)

temp.mk <- H1.fm.mk$separate %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup()
 
temp.mk <- temp.mk %>% bind_rows(H1.fm.mk$merge %>% filter(cluster=="beta") %>% filter(p_val_adj < 0.05) %>% filter(!gene %in% temp.mk$gene) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup())

temp.mk <- temp.mk %>% bind_rows(H1.fm.mk$wp %>% filter(cluster=="polyhormonal") %>% filter(p_val_adj < 0.05) %>% filter(!gene %in% temp.mk$gene) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup())
temp.mk <-  (temp.mk %>% split(.,.$cluster))[mk.od] %>% lapply(function(x){x %>% arrange(p_val_adj) %>% pull(gene)})


temp.M <- data.ob.umap %>% select(cell,EML) %>% mutate(od=factor(EML,EML.od,ordered = T)) %>% arrange(od) %>% select(-od)
temp.exp <- data.ob@assays$RNA$data[unlist(temp.mk) ,temp.M$cell]
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
zs.limit <- 2.5
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
temp.anno <- temp.M  %>% tibble::column_to_rownames("cell")
pheatmap::pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=4,border_color="NA",gaps_row =unlist(lapply( temp.mk,function(x){return(length(x))})) %>% cumsum(),main="batch1_2,marker miRNA") #%>%ggplotify::as.ggplot()
# EPHA5， TMEM190， TLE4