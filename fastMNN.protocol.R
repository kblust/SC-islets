#' ---
#' title: "protocal fastMNN added Rajaei 2025"
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
  
  metas.list$Rajaei_2025 <- readRDS(paste0("tmp_data/","Rajaei_2025","/Rajaei_2025.meta.filter.rds")) 
  counts.list$Rajaei_2025 <- readRDS(file=paste0("tmp_data/","Rajaei_2025","/Rajaei_2025.counts.filter.rds"))[,metas.list$Rajaei_2025$cell]
  
  
  meta.filter <- metas.list %>% do.call("bind_rows",.)
  lapply(counts.list,nrow)
  ov.genes <- rownames(counts.list$H1) %>% intersect(rownames(counts.list$Veres_2019)) %>% intersect(rownames(counts.list$BalBoa_2022)) %>% intersect(rownames(counts.list$Augsor_2022)) %>% intersect(rownames(counts.list$Rajaei_2025))
  counts.filter <- counts.list$H1[ov.genes,] %>% cbind(counts.list$Veres_2019[ov.genes,]) %>% cbind(counts.list$BalBoa_2022[ov.genes,]) %>% cbind(counts.list$Augsor_2022[ov.genes,])%>% cbind(counts.list$Rajaei_2025[ov.genes,])
  
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
    
    mBN.sce.ob <- multiBatchNorm(sce.ob$H1,sce.ob$Veres_2019,sce.ob$BalBoa_2022,sce.ob$Augsor_2022,sce.ob$Rajaei_2025)
    lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
    saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/DP.fastMNN.lognormExp.mBN.rds"))
    
  }
  # sel.od <- c("Rajaei_2025","BalBoa_2022","Augsor_2022","H1","Veres_2019")
  #sel.od <- c("Rajaei_2025","Veres_2019","BalBoa_2022","Augsor_2022","H1")
  #sel.od <- c("BalBoa_2022","Rajaei_2025","Veres_2019","Augsor_2022","H1") #bad 
  #sel.od <- c("Rajaei_2025","Veres_2019","BalBoa_2022","Augsor_2022","H1")
  sel.od <- c("Rajaei_2025","Augsor_2022","H1","BalBoa_2022","Veres_2019")
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
  data.temp <- data.ob %>% FindClusters(reso=0.8,verbose=F)
  
  DimPlot(data.temp,label=T,split.by = "pj",group.by = "EML")+NoLegend()
  DimPlot(data.temp,label=T,split.by = "pj")+NoLegend()
  plot_grid(plotlist = FunFP_plot(data.temp,unlist(main.mk)))

  
  cowplot::plot_grid(
    DimPlot(data.temp,label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "EML",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "pj",label=T)+NoLegend()+NoAxes(),
    DimPlot(data.temp,group.by = "devTime",label=T)+NoLegend()+NoAxes()
  )
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% select(cell,SID:mt.perc,EML,prolifSig) %>% mutate(SC=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")  # 
  data.ob.umap <- data.ob.umap  %>% mutate(cluster_EML=recode(SC,'C0'="alpha","C1"="alpha",'C2'="late_beta","C3"="early_beta", 'C4'="SCEC",'C5'="SCEC",'C6'="alpha",'C7'="exo",'C8'="alpha",'C9'="exo",'C10'="exo",'C11'="polyhormonal",'C12'="SCEC",'C13'="prolif",'C14'="delta",'C15'="SCEC",'C16'="polyhormonal",'C17'="exo",'C18'="exo",'C19'="exo",'C20'="exo",'C21'="NeuroEndo"))
  
  Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"cluster_EML"])
  DimPlot(data.temp,label=T)+NoAxes()+NoLegend()
  FeaturePlot(data.temp,"prolifSig")
  
  saveRDS(data.ob,paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.umap.rds"))
  
}

#'
data.temp <- data.ob %>% FindClusters(reso=0.8,verbose=F)
Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"cluster_EML"])
DimPlot(data.temp,label=T)+NoAxes()+NoLegend()

data.ob.umap %>% group_by(pj,cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% spread(cluster_EML,nCell)
#' data.temp <- JoinLayers(data.temp)
#' temp.mk=FindMarkers(data.temp,ident.1 = "12") %>% tibble::rownames_to_column("gene")
#' DimPlot(data.ob,group.by = "EML",label=T,split.by="pj")+NoAxes()+NoLegend()
#' temp.SG <- c("LIN28A","CDH26","KLF9","GYPC","CDX2","ANGPT2","SPOCK1","DKK2","CDH7","SOX9","LCP1","KRT19","KRT17","MUC1","CFTR","PTF1A","CPA1","CPA2","CTRB2","VIM","SPARC","COL3A1","COL1A1","PDGFRB") 
#' 
#' temp.plot <- list()
#' temp.plot[["UMAP"]] <- DimPlot(data.ob,group.by = "EML",label=T)+NoAxes()+NoLegend()
#' for (g in temp.SG) {
#'   if (g %in% rownames(data.ob@assays$RNA@features)) {
#'     temp.plot[[g]] <-  FeaturePlot(data.ob,g)+NoAxes()+NoLegend()
#'   }
#' }
#' cowplot::plot_grid(plotlist = temp.plot)
#' 
#' #' check cell number 
#' data.ob.umap %>% mutate(EML=new_EML) %>% group_by(devTime,EML) %>% summarise(nCell=n_distinct(cell)) %>% spread(EML,nCell)%>% replace(.,is.na(.),0)
#' print(
#'   data.ob.umap %>% mutate(EML=new_EML)  %>% group_by(devTime,EML) %>% summarise(nCell=n_distinct(cell))  %>% spread(EML,nCell)%>% replace(.,is.na(.),0) %>% gather(EML,nCell,-devTime) %>% group_by(devTime) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot()+geom_bar(mapping=aes(x=EML,y=prop,fill=devTime),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90))+xlab("")+ggtitle("Veres et al")+FunTitle()
#' )
#' 
#' 
#' #' check the beta_dis and beta  (beta_dis ribo gene highly expressed)
#' 
#' temp.DEG <- DEG.mk$detail %>% filter(BG=="beta" & EG=="beta_dis") %>% filter(power > 0.4) %>% filter(pct.2 < 0.3) %>% mutate(UpDown="UpRe") %>% bind_rows(DEG.mk$detail %>% filter(EG=="beta" & BG=="beta_dis") %>% filter(power > 0.4) %>% mutate(UpDown="DownRe")) %>% group_by(UpDown) %>% top_n(15,power)#%>% filter(!gene %in% ribo.gene) 
#' temp.M <- data.ob.umap %>% filter(new_EML %in% c("beta","beta_dis")) %>% arrange(new_EML,seurat_clusters,devTime) %>% select(cell,new_EML,seurat_clusters,devTime)
#' FunPreheatmapNoLog(data.ob@assays$RNA@data,temp.DEG$gene,temp.M$cell) %>% pheatmap::pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,annotation_col = (temp.M %>% tibble::column_to_rownames("cell")))
#' 
#' 
#' 
#' data.temp <- data.ob %>% FindClusters(reso=1,verbose=F)
#' data.temp@meta.data$new_EML <- (data.ob.umap %>% column_to_rownames("cell"))[rownames(data.temp@meta.data),"new_EML"]
#' VlnPlot(data.temp,group.by = "new_EML",c("ABCC8","MIAT","PLXNA2","PNISR"))
#' plot_grid(
#'   DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
#'   DimPlot(data.temp,group.by="EML",label=T)+NoAxes()+NoLegend(),
#'   DimPlot(data.temp,group.by="new_EML",label=T)+NoAxes()+NoLegend(),
#'   DimPlot(data.temp,group.by="devTime",label=T)+NoAxes()+NoLegend(),
#'   FeaturePlot(data.temp,"nGene")+NoAxes()+NoLegend(),
#'   #DimPlot(data.temp,group.by="batch",label=T)+NoAxes()+NoLegend(),
#'   FeaturePlot(data.temp,"mt.perc")+NoAxes()+NoLegend(),
#'   DimPlot(data.temp,label=T,group.by = "Phase")+NoAxes()+NoLegend(),
#'   FeaturePlot(data.temp,"beta_PSDT")+NoAxes()+NoLegend()
#' )
#' APPJ_devTime(data.temp,"Veres_2019")
#' APPJ(data.temp,"Veres_2019")
#' DimPlot(data.temp,group.by="new_EML",label=T)+NoAxes()+NoLegend()
#' 
#' plot_grid(
#'   plotlist = FunFP_plot(data.temp,unlist(human.mk.list))
#' )
#' # plot_grid(plotlist = FunFP_plot(data.temp,c("GAP43","ONECUT3","NEUROG3"))) ### GAP43 markers
#' # plot_grid(plotlist = FunFP_plot(data.temp,c("RPL3","PCP4","QDPR","MAP1B","DCX","RTN1","FXYD6"))) ### chec some sub markers
#' 
#' temp.plot <- list()
#' data.sub.temp <- data.ob %>% subset(cell=(data.ob.umap %>% filter(new_EML %in% c("beta")) %>% pull(cell)))
#' for (g in unlist(beta.stage.mk)) {
#'   temp.plot[[g]] <- VlnPlot(data.sub.temp,g,group.by = "devTime")+NoLegend()
#' }
#' plot_grid(plotlist = temp.plot)
#' 
#' 
#' # sub beta important
#' plot_grid(plotlist = FunFP_plot(data.temp,c("MIAT","ABCC8","WSB1","IGFBP5","INS","GCG","FEV","GAP43","CHGA"))) ##### disfunctional beta cells
#' plot_grid(plotlist = FunFP_plot(data.temp,c("SST","TOP2A","CDK1","NEUROG3","FOXJ1","ARX","TTR","CHGA","PCSK2"))) ##### other interesting to check
#' 
#' #' check major pancreatic cell type
#' plot_grid(
#'   plotlist = FunFP_plot(data.temp,c(unlist(human.mk.list[c(1:3,5,6,7)]),"ARX"))
#' )


temp.M <- data.ob.umap %>% filter(pj=="Veres_2019") %>% pull(subCT)

temp.plot <- list()
for ( n in unique(data.ob.umap %>% filter(pj=="Veres_2019") %>% pull(subCT) %>% unique())) {
  temp.plot[[n]] <-DimPlot(data.temp,cells.highlight = (data.ob.umap %>% filter(pj=="Veres_2019") %>% filter(subCT==n) %>% pull(cell)))+ggtitle(n)+FunTitle()+NoAxes()+NoLegend()
}
print(cowplot::plot_grid(plotlist=temp.plot))

temp.plot <- list()
for ( n in unique(data.ob.umap %>% filter(pj=="Veres_2019" & EML=="other") %>% pull(subCT) %>% unique())) {
  temp.plot[[n]] <-DimPlot(data.temp,cells.highlight = (data.ob.umap %>% filter(pj=="Veres_2019") %>% filter(subCT==n) %>% pull(cell)))+ggtitle(n)+FunTitle()+NoAxes()+NoLegend()
}
print(cowplot::plot_grid(plotlist=temp.plot))
