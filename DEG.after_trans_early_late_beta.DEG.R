#' ---
#' title: "DEG late vs early beta in after and early trans"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library
#R4.3
rm(list=ls())
# condaENV <- "/home/chenzh/miniconda3/envs/R4.3"
# LBpath <- paste0(condaENV ,"/lib/R/library")
# .libPaths(LBpath)

suppressPackageStartupMessages({
  library(clusterProfiler)
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


data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.ob.umap.rds"))

PATH_list <- list()
for (tp in c("KEGG_2021_Human")) {
  temp <- read.delim(paste("~/Genome_new/enrichr_library/modified/",tp,".mod.txt",sep=""),head=F,stringsAsFactors=F)
  colnames(temp) <- c("TERM","GENE")
  #rownames(temp) <- temp$TERM
  PATH_list[[tp]] <- temp
}

gsea.results <- list()
pt.results <- list()
DEG.results <- list()
fgsea.results <- list()
if (file.exists(paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.pt.results.rds"))) {
  DEG.results <- readRDS(paste0("tmp_data/",TD,"/trans.after.DEG.cluster.results.rds"))
  pt.results <- readRDS(paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.pt.results.rds"))
  gsea.results <- readRDS(paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.gsea.results.rds"))
  fgsea.results <- readRDS(paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.fgsea.results.rds"))
  
}else{
  load("tmp_data/gene.meta.Rdata",verbose=T)
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(devTime=="HS980_trans")
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))[,meta.filter$cell]
  
  lognormExp.mBN <- readRDS(file=paste0("tmp_data/",TD,"/before.after.trans.fastMNN.lognormExp.mBN.rds"))
  
  
  #' create psd object
  temp.sel.expG <- rownames(lognormExp.mBN)
  temp.M <- meta.filter %>% filter(devTime=="HS980_trans")
  
  data.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.ob@assays$RNA$data <- as.matrix(lognormExp.mBN[temp.sel.expG,rownames(data.ob@meta.data)])
  Idents(data.ob) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"])
  
  
  
  library(topGO)
  
  genename.bg <- ALL_gene <- rownames(counts.filter)
  
  geneID2GO=inverseList(readMappings(file = "~/Genome_new/Human/GO/Human.GO2geneID_ALL.refseq.map"))
  
  data.deg <- data.ob
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"EML"])
  
  temp.DEG <- FindMarkers(data.deg,ident.1 = "late_beta",ident.2="early_beta") %>% tibble::rownames_to_column("gene") %>% tbl_df()
  temp.compair <- "late_beta_vs_early_beta"
  DEG.results[[temp.compair]] <- list()
  DEG.results[[temp.compair]][["DEG.all.results"]] <- temp.DEG 
  DEG.results[[temp.compair]][["DEG.result.up"]] <- temp.DEG %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25 & pct.1 >=1/4)
  DEG.results[[temp.compair]][["DEG.result.down"]] <- temp.DEG %>% filter(p_val_adj < 0.05 & avg_log2FC <  -0.25  & pct.2 >=1/4 )
  DEG.results[[temp.compair]][["DEG.result"]] <- DEG.results[[temp.compair]][["DEG.result.up"]] %>% bind_rows(DEG.results[[temp.compair]][["DEG.result.down"]])
  
  temp.up.ID <- DEG.results[[temp.compair]][["DEG.result.up"]]$gene
  temp.down.ID <- DEG.results[[temp.compair]][["DEG.result.up"]]$gene
  DEG.results[[temp.compair]] [["DEG.up.GO.result"]] <- topGO_enrichment2(temp.up.ID,ALL_gene,0.05,geneID2GO)
  DEG.results[[temp.compair]] [["DEG.down.GO.result"]] <- topGO_enrichment2(temp.down.ID,ALL_gene,0.05,geneID2GO)
  
  
  genename.bg <-  ALL_gene ### background genes
  
  pt.results <- list()
  for(n in names(DEG.results)) {
    print(n)
    temp.compair <- n
    pt.results[[temp.compair]] <- list()
    
    temp.DEG <- DEG.results[[temp.compair]]$DEG.result
    up.DEG <- DEG.results[[temp.compair]]$DEG.result.up$gene
    down.DEG <- DEG.results[[temp.compair]]$DEG.result.down$gene
    for (tp in c("KEGG_2021_Human")) {
      pt.results[[temp.compair]][[tp]]  <- list()
      
      suppressMessages(up.DEG.pt <- enricher(up.DEG, pvalueCutoff = 0.05, pAdjustMethod = "none",genename.bg ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=PATH_list[[tp]]))
      suppressMessages(down.DEG.pt <- enricher(down.DEG, pvalueCutoff = 0.05, pAdjustMethod = "none",genename.bg ,minGSSize = 5,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=PATH_list[[tp]]))
      
      pt.results[[temp.compair]][[tp]]$up <- up.DEG.pt
      pt.results[[temp.compair]][[tp]]$down <- down.DEG.pt
      
      #' GSEA results
      #' 
      temp.all.DEG <- DEG.results[[temp.compair]]$DEG.all.results
      temp.fc <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(avg_log2FC)
      names(temp.fc) <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(gene)
      temp.gsea <- GSEA(temp.fc,pvalueCutoff = 0.25, TERM2GENE =PATH_list[[tp]] )
      temp.gsea <- GSEA(temp.fc,pvalueCutoff = 0.25, TERM2GENE =PATH_list[[tp]] )
      gsea.results[[temp.compair]][[tp]] <- temp.gsea
      
      
      temp.sel.tp.path <- PATH_list[[tp]] %>% split(.,.$TERM) %>% lapply(function(x){x$GENE})
      temp.fgsea.out <-fgsea::fgsea(pathways = temp.sel.tp.path , stats = temp.fc)
      fgsea.results[[temp.compair]][[tp]] <- temp.fgsea.out 
      
    }
  }
  
  
  
  
  saveRDS(DEG.results,paste0("tmp_data/",TD,"/trans.after.DEG.cluster.results.rds"))
  saveRDS(pt.results,paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.pt.results.rds"))
  saveRDS(gsea.results,paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.gsea.results.rds"))
  saveRDS(fgsea.results,paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.fgsea.results.rds"))
  
}