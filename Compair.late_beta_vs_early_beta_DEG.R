#' ---
#' title: ""
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
  library(mascarade)
  library(ComplexHeatmap)
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


heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)


meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))

#' for H1 related results
data.H1.ob.umap<- readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.umap.rds"))

#' for  HS980 before and after transplantation
data.BA.trans.umap <-  readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.ob.umap.rds"))


DEG.results <- list()
#' loading late vs early beta cells
DEG.results[["H1"]] <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.results.rds"))
DEG.results[["BT"]] <- readRDS(paste0("tmp_data/",TD,"/trans.before.DEG.cluster.results.rds"))
DEG.results[["AT"]] <- readRDS(paste0("tmp_data/",TD,"/trans.after.DEG.cluster.results.rds"))

#' fgsea results 
H1.fgsea.results <-  readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.fgsea.results.rds"))
BT.fgsea.results <- readRDS(paste0("tmp_data/",TD,"/trans.before.DEG.cluster.MSigR.fgsea.results.rds"))
AT.fgsea.results <- readRDS(paste0("tmp_data/",TD,"/trans.after.DEG.cluster.MSigR.fgsea.results.rds"))


#' pathway
PATH_list <- list()
for (tp in c("KEGG_2021_Human")) {
  temp <- read.delim(paste("~/Genome_new/enrichr_library/modified/",tp,".mod.txt",sep=""),head=F,stringsAsFactors=F)
  colnames(temp) <- c("TERM","GENE")
  #rownames(temp) <- temp$TERM
  PATH_list[[tp]] <- temp
}

plot.results <- list()


ggvenn::ggvenn(
  list(
  H1=DEG.results[["H1"]]$late_beta_vs_early_beta$DEG.result %>% mutate(UpDown=ifelse(avg_log2FC >0,"UpRe","DownRe")) %>% mutate(SID=paste(gene,UpDown,sep=":")) %>% pull(SID),
  BT=DEG.results[["BT"]]$late_beta_vs_early_beta$DEG.result %>% mutate(UpDown=ifelse(avg_log2FC >0,"UpRe","DownRe")) %>% mutate(SID=paste(gene,UpDown,sep=":")) %>% pull(SID),
  AT=DEG.results[["AT"]] $late_beta_vs_early_beta$DEG.result %>% mutate(UpDown=ifelse(avg_log2FC >0,"UpRe","DownRe")) %>% mutate(SID=paste(gene,UpDown,sep=":")) %>% pull(SID)
  )
)

#' check
temp.genes <- c(DEG.results[["H1"]]$late_beta_vs_early_beta$DEG.result$gene,DEG.results[["BT"]]$late_beta_vs_early_beta$DEG.result$gene,DEG.results[["AT"]]$late_beta_vs_early_beta$DEG.result$gene)


temp.input <- DEG.results[["H1"]]$late_beta_vs_early_beta$DEG.all.results %>% filter(gene %in% temp.genes) %>% mutate(H1=avg_log2FC) %>% select(gene,H1) %>% inner_join(DEG.results[["BT"]]$late_beta_vs_early_beta$DEG.all.results %>% filter(gene %in% temp.genes) %>% mutate(BT=avg_log2FC) %>% select(gene,BT),by="gene") %>% inner_join(DEG.results[["AT"]]$late_beta_vs_early_beta$DEG.all.results %>% filter(gene %in% temp.genes) %>% mutate(AT=avg_log2FC) %>% select(gene,AT),by="gene") %>% tibble::column_to_rownames("gene")
#temp.input[temp.input >3] <- 3
#temp.input[temp.input < -3] <- -3
temp.input[] %>% pheatmap::pheatmap(scale="none",cluster_rows = T,cluster_cols = F,show_rownames = F,col=heat.col)


# check trans influence
temp.genes <- c(DEG.results[["BT"]]$late_beta_vs_early_beta$DEG.result$gene,DEG.results[["AT"]]$late_beta_vs_early_beta$DEG.result$gene)
temp.input <- DEG.results[["BT"]]$late_beta_vs_early_beta$DEG.all.results %>% filter(gene %in% temp.genes) %>% mutate(BT=avg_log2FC) %>% select(gene,BT) %>% inner_join(DEG.results[["AT"]]$late_beta_vs_early_beta$DEG.all.results %>% filter(gene %in% temp.genes) %>% mutate(AT=avg_log2FC) %>% select(gene,AT),by="gene") %>% mutate(BT_type=ifelse(BT>0,"Up","Down"), AT_type=ifelse(AT>0,"Up","Down")) %>% mutate(SID=paste(BT_type,AT_type,sep=":")) #%>% gather(trans,log2FC,-gene)

temp.sel.genes <- c("SLC30A8","ASCL1","MAFA","SCG2","CHGB","C7","SST","IAPP","VGF","WNT4","CDH1","MAPK1","ROCK1")

plot.results$BT.AT.DEG.log2FC <- temp.input  %>% ggplot()+geom_point(mapping = aes(x=BT,y=AT,col=SID),size=0.5)+ggtitle("log2FC")+theme_classic()+FunTitle()+ylim(-4,5)+xlim(-4,7)+scale_color_manual(values=c("Up:Up"=as.vector(EML.lineage.col.set["late_beta"]),"Down:Down"=as.vector(EML.lineage.col.set["early_beta"]),"Up:Down"=as.vector(EML.lineage.col.set["exo"]),"Down:Up"=as.vector(EML.lineage.col.set["psc"])))+ggrepel::geom_label_repel(temp.input %>% filter(gene %in% temp.sel.genes),mapping=aes(x=BT,y=AT,label=gene,col=SID),label.size=NA,fill=NA,max.overlaps=25,box.padding = 0.5)
plot.results$BT.AT.DEG.log2FC 

temp.input %>% filter(BT*AT < 0) %>% pull(gene)
temp.input %>% filter(BT*AT < 0 & AT >0) %>% pull(gene)
temp.input %>% filter(BT*AT < 0 & AT <0) %>% pull(gene)
temp.input %>% filter(BT*AT > 0) %>% pull(gene)

#' #### Pathway enrichment analysis in BT
temp.compair <- "late_beta_vs_early_beta"
temp.all.DEG <- DEG.results[["BT"]][[temp.compair]]$DEG.all.results
temp.fc <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(avg_log2FC)
names(temp.fc) <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(gene)
temp.sc <- c("Hippo signaling pathway","TGF-beta signaling pathway")
temp.fgsea.out <- BT.fgsea.results[[temp.compair]]$KEGG_2021_Human
temp.sel.tp.path <- PATH_list[[tp]] %>% split(.,.$TERM) %>% lapply(function(x){x$GENE})

temp.plot <- list()
for (temp.sel.id in c(temp.sc)) {
  temp.SID <- temp.fgsea.out %>% as.data.frame() %>% tbl_df() %>% filter(pathway==temp.sel.id)  %>% mutate(SID=paste0("NES:",round(NES,3),",","Pvalue:",round(pval,3))) %>% pull(SID)
  temp.plot[[temp.sel.id]] <- fgsea::plotEnrichment(temp.sel.tp.path[[temp.sel.id]],temp.fc) + ggtitle(paste(temp.sel.id,temp.SID,"",sep="\n"))+theme(plot.title = element_text(hjust =0.5))
}
cowplot::plot_grid(plotlist = temp.plot)
plot.results$DEG.BT.EL.beta.fgsea=temp.plot



#' #### Pathway enrichment analysis in AT
temp.compair <- "late_beta_vs_early_beta"
temp.all.DEG <- DEG.results[["AT"]][[temp.compair]]$DEG.all.results
temp.fc <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(avg_log2FC)
names(temp.fc) <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(gene)
temp.sc <- c("Hippo signaling pathway","TGF-beta signaling pathway")
temp.fgsea.out <- AT.fgsea.results[[temp.compair]]$KEGG_2021_Human
temp.sel.tp.path <- PATH_list[[tp]] %>% split(.,.$TERM) %>% lapply(function(x){x$GENE})

temp.plot <- list()
for (temp.sel.id in c(temp.sc)) {
  temp.SID <- temp.fgsea.out %>% as.data.frame() %>% tbl_df() %>% filter(pathway==temp.sel.id)  %>% mutate(SID=paste0("NES:",round(NES,3),",","Pvalue:",round(pval,3))) %>% pull(SID)
  temp.plot[[temp.sel.id]] <- fgsea::plotEnrichment(temp.sel.tp.path[[temp.sel.id]],temp.fc) + ggtitle(paste(temp.sel.id,temp.SID,"",sep="\n"))+theme(plot.title = element_text(hjust =0.5))
}
cowplot::plot_grid(plotlist = temp.plot)
plot.results$DEG.AT.EL.beta.fgsea=temp.plot




# check conserved fgsea results
BT.fgsea.results[[temp.compair]]$KEGG_2021_Human %>% tbl_df() %>% filter(pval <0.05) %>% select(pathway,pval,NES) %>% inner_join(H1.fgsea.results[[temp.compair]]$KEGG_2021_Human %>% tbl_df() %>% filter(pval <0.05) %>% select(pathway,pval,NES) ,by="pathway")



#‘ some other check 
# check trans influence
temp.genes <- c(DEG.results[["BT"]]$late_beta_vs_early_beta$DEG.result$gene,DEG.results[["H1"]]$late_beta_vs_early_beta$DEG.result$gene)
temp.input <- DEG.results[["BT"]]$late_beta_vs_early_beta$DEG.all.results %>% filter(gene %in% temp.genes) %>% mutate(BT=avg_log2FC) %>% select(gene,BT) %>% inner_join(DEG.results[["H1"]]$late_beta_vs_early_beta$DEG.all.results %>% filter(gene %in% temp.genes) %>% mutate(H1=avg_log2FC) %>% select(gene,H1),by="gene") #%>% gather(trans,log2FC,-gene)
temp.input  %>% ggplot()+geom_point(mapping = aes(x=BT,y=H1))+ggtitle("log2FC")+FunTitle()

