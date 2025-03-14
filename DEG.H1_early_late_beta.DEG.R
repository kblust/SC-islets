#' ---
#' title: "DEG late vs early beta in H1"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library
#R4.0
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
if (file.exists(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.pt.results.rds"))) {
  DEG.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.results.rds"))
  pt.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.pt.results.rds"))
  gsea.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.gsea.results.rds"))
  fgsea.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.fgsea.results.rds"))
  
}else{
  load("tmp_data/gene.meta.Rdata",verbose=T)
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(devTime=="H1")
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))[,meta.filter$cell]
  
  #' loading dataset
  heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
  
  if (file.exists(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.rds"))) {
    data.ob <-  readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.rds"))
    data.ob.umap<- readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.umap.rds"))
  }
  

  library(topGO)
  
  genename.bg <- ALL_gene <- rownames(readRDS(paste0("tmp_data/",TD,"/counts.filter.rds")))

  geneID2GO=inverseList(readMappings(file = "~/Genome_new/Human/GO/Human.GO2geneID_ALL.refseq.map"))
  
  data.deg <- data.ob
  Idents(data.deg) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.deg@meta.data),"EML"])
  
  temp.DEG <- FindMarkers(data.deg,ident.1 = "late_beta",ident.2="early_beta") %>% tibble::rownames_to_column("gene") %>% tbl_df()
  temp.compair <- "late_beta_vs_early_beta"
  DEG.results[[temp.compair]] <- list()
  DEG.results[[temp.compair]][["DEG.all.results"]] <- temp.DEG 
  DEG.results[[temp.compair]][["DEG.result.up"]] <- temp.DEG %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25 & pct.1 >=1/4 )
  DEG.results[[temp.compair]][["DEG.result.down"]] <- temp.DEG %>% filter(p_val_adj < 0.05 & avg_log2FC <  -0.25 & pct.2 >=1/4)
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
  

 
  
  saveRDS(DEG.results,paste0("tmp_data/",TD,"/H1.DEG.cluster.results.rds"))
  saveRDS(pt.results,paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.pt.results.rds"))
  saveRDS(gsea.results,paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.gsea.results.rds"))
  saveRDS(fgsea.results,paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.fgsea.results.rds"))
  
}

sel.pathway <- c("Diabetic cardiomyopathy","Maturity onset diabetes of the young","Insulin secretion","Thermogenesis","Oxidative phosphorylation","Ferroptosis","Maturity onset diabetes of the young","Carbohydrate digestion and absorption","Insulin resistance","Type II diabetes mellitus","AMPK signaling pathway","cAMP signaling pathway","Cholesterol metabolism","PPAR signaling pathway","Protein processing in endoplasmic reticulum")
sel.genes <- c("IAPP","NEFM","SIX2","LDHA","NEUROD1","PAX6","PDX1","INS","ABCC8","CACNA1D","G6PC2","SLC37A4","IGF1R","PKM","PPARA","ADCYAP1","ADCY2","RAP1B","RAP1A","CD36","APOC3","ABCA1","SCD","ERO1B","P4HB","MAP1LC3A","ASPH", "ASCL1","BACE2", "PCDH7") #"HOPX","IGF2","NEFM","LDHA"


#' #### Pathway enrichment analysis
temp.compair <- "late_beta_vs_early_beta"
temp.all.DEG <- DEG.results[[temp.compair]]$DEG.all.results
temp.fc <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(avg_log2FC)
names(temp.fc) <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(gene)

temp.sc <- c("Glucagon signaling pathway","Insulin resistance","MAPK signaling pathway","FoxO signaling pathway","AMPK signaling pathway","JAK-STAT signaling pathway","PPAR signaling pathway","Oxidative phosphorylation","Thermogenesis")

temp.gsea.out <- fgsea.results[[temp.compair]]$KEGG_2021_Human
temp.plot <- list()
for (temp.sel.id in c(temp.sc)) {
    temp.SID <- temp.gsea.out %>% as.data.frame() %>% tbl_df() %>% filter(pathway==temp.sel.id)  %>% mutate(SID=paste0("NES:",round(NES,3),",","Pvalue:",round(pval,3))) %>% pull(SID)
    temp.plot[[temp.sel.id]] <- fgsea::plotEnrichment(temp.sel.tp.path[[temp.sel.id]],temp.fc) + ggtitle(paste(temp.sel.id,temp.SID,"",sep="\n"))+theme(plot.title = element_text(hjust =0.5))
}
cowplot::plot_grid(plotlist = temp.plot)
plot.results$DEG.EL.beta.fgsea=temp.plot




pdf("tmp_data/temp.DEG.fgsea.plot.pdf",8,8)
temp.plot <- plot.results$DEG.EL.beta.fgsea
cowplot::plot_grid(plotlist = temp.plot)
dev.off()

#' output tables
DEG.results$late_beta_vs_early_beta$DEG.result.up %>% mutate(UpDown="Up-Re") %>% bind_rows(DEG.results$late_beta_vs_early_beta$DEG.result.down %>% mutate(UpDown="Down-Re") ) %>% mutate(Comparison="Late_beta_vs_early_beta") %>% write.table("tmp_data/late_beta_vs_early_beta/DEG",quote=F,sep="\t",col.names = T,row.names = F)

pt.results$late_beta_vs_early_beta$KEGG_2021_Human$up %>% as.data.frame() %>% tbl_df() %>% mutate(UpDown="Up-Re") %>% bind_rows(pt.results$late_beta_vs_early_beta$KEGG_2021_Human$down %>% as.data.frame() %>% tbl_df() %>% mutate(UpDown="Down-Re")) %>% select(ID,GeneRatio,pvalue,geneID,Count,UpDown)  %>% mutate(Comparison="Late_beta_vs_early_beta") %>% filter(pvalue < 0.05) %>% rename(Pathway=ID)%>% write.table("tmp_data/late_beta_vs_early_beta/enriched.KEGG",quote=F,sep="\t",col.names = T,row.names = F)

fgsea.results$late_beta_vs_early_beta$KEGG_2021_Human%>% as.data.frame() %>% tbl_df() %>% mutate(Pathway=pathway,pvalue=pval) %>% select(Pathway,pvalue,NES)  %>% filter(pvalue < 0.05)  %>% mutate(Comparison="Late_beta_vs_early_beta")%>% write.table("tmp_data/late_beta_vs_early_beta/GSEA",quote=F,sep="\t",col.names = T,row.names = F)


data.temp <- readRDS("data/seurat_ob_from_kelly/comp_trans_int.rds")
DimPlot(data.temp)



#perl  ~/PC/code/txt2excel-1.0.pl -f DEG,enriched.KEGG,GSEA -X sup.lateVSearlyBeta.xls
###################### don't run 
# if (FALSE) {
#   #GS_db =  readRDS("big_doc/Human.msigdbr.rds")
#   for (m in names(DEG.results)) {
#     temp.compair <- m
#     temp.DEG <- DEG.results[[temp.compair]]$DEG.result  %>% arrange(desc(avg_log2FC)) ### mast
#     temp.fc <-  temp.DEG$avg_log2FC
#     names(temp.fc) <- temp.DEG$gene
#     for (tp in c("CGP","CP","CP:BIOCARTA", "CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS", "TFT:GTRD")) {
#       print(paste(m,tp))
#       temp.path <- GS_db %>% filter(gs_subcat == tp) %>% mutate(TERM=gs_name,GENE=gene_symbol) %>% select(TERM,GENE) %>% unique() %>% as.data.frame()
#       suppressMessages(pt.results[[temp.compair]][[tp]]$up  <- enricher(DEG.results[[temp.compair]]$DEG.result.up$gene, pvalueCutoff = 1, pAdjustMethod = "BH",genename.bg ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=temp.path))
#       suppressMessages(pt.results[[temp.compair]][[tp]]$down  <- enricher(DEG.results[[temp.compair]]$DEG.result.down$gene, pvalueCutoff = 1, pAdjustMethod = "BH",genename.bg ,minGSSize = 2,maxGSSize =1000,qvalueCutoff = 1, TERM2GENE=temp.path))
#       suppressMessages(gsea.results[[temp.compair]][[tp]] <-  GSEA(temp.fc,pvalueCutoff = 0.5, pAdjustMethod = "BH",minGSSize = 2,maxGSSize = 1000, TERM2GENE=temp.path))
#     }
#   }
#   saveRDS(DEG.results,paste0("tmp_data/",TD,"/DEG.cluster.results.rds"))
#   saveRDS(pt.results,paste0("tmp_data/",TD,"/DEG.cluster.MSigR.pt.results.rds"))
#   saveRDS(gsea.results,paste0("tmp_data/",TD,"/DEG.cluster.MSigR.gsea.results.rds"))
#   #+ fig.width=15,fig.height=9
#   pdf("tmp_data/pt.way.temp.pdf",15,9)
#   for (temp.compair in names(DEG.results)) {
#     for (tp in names(pt.results[[names(DEG.results)]] )) {
#       up.DEG.pt <- pt.results[[temp.compair]][[tp]]$up
#       down.DEG.pt <- pt.results[[temp.compair]][[tp]]$down
#       if (! is.null(up.DEG.pt) ) {
#         if (nrow(up.DEG.pt)!=0) {
#           up.DEG.pt.sig=subset(as.data.frame(up.DEG.pt),up.DEG.pt$pvalue <0.05)[,c(1,3,4,5,8,9)]
#           print(paste("up-regulated DEGs in",temp.compair,paste("(",tp,")",sep="")))
#           print(dotplot(up.DEG.pt,color="pvalue",title=paste("up-regulated DEGs in",temp.compair,paste("(",tp,")",sep=""))))
#         }else{
#           print(paste("No significant item in up-regulated genes"))
#         }
#       }else{
#         print(paste("No significant item in up-regulated genes"))
#       }
#       if (! is.null(down.DEG.pt) ) {
#         if (nrow(down.DEG.pt)!=0) {
#           down.DEG.pt.sig=subset(as.data.frame(down.DEG.pt),down.DEG.pt$pvalue <0.05)[,c(1,3,4,5,8,9)]
#           print(paste("down-regulated DEGs in",temp.compair,paste("(",tp,")",sep="")))
#           print(dotplot(down.DEG.pt,color="pvalue",title=paste("down-regulated DEGs in",temp.compair,paste("(",tp,")",sep=""))))
#         }else{
#           print(paste("No significant item in down-regulated genes"))
#         }
#         
#       }else{
#         print(paste("No significant item in down-regulated genes"))
#       }
#     }
#   }
#   dev.off()
# }


#' #### Pathway enrichment analysis
# temp.compair <- "late_beta_vs_early_beta"
# temp.gsea <- gsea.results[[temp.compair]]$KEGG_2021_Human
# temp.all.DEG <- DEG.results[[temp.compair]]$DEG.all.results
# temp.fc <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(avg_log2FC)
# names(temp.fc) <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(gene)
# temp.fc.sel <- temp.fc[names(temp.fc) %>% intersect(DEG.results[[temp.compair]]$DEG.result$gene)]
# temp.gsea.out <- temp.gsea  %>% as.data.frame() %>% tbl_df()
# 
# temp.sc <- c("Glucagon signaling pathway","Oxidative phosphorylation","Insulin resistance","AGE-RAGE signaling pathway in diabetic complications")
#gseaplot(temp.gsea, by = "all", title = "Glucagon signaling pathway", geneSetID = 1)
#dotplot(temp.gsea, showCategory=10, split=".sign")
#emapplot(temp.gsea, showCategory = 10)
#cnetplot(temp.gsea, categorySize="pvalue", foldChange=temp.fc.sel*4, showCategory = temp.sc,cex.gene=0.25)
