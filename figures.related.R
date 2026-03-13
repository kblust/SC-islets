#' ---
#' title: "generate the figures"
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
  #library(mascarade)
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




meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))

#' for H1 related results
data.H1.ob.umap<- readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.umap.rds"))

H1.fm.mk <- readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.fm.mk.rds"))

#' for protocal comparison
data.DP.ob.umap <- readRDS(paste0("tmp_data/",TD,"/DP.fastMNN.data.ob.umap.rds"))

#' for H1 with HS980
data.H1_H9S980.umap <- readRDS(paste0("tmp_data/",TD,"/H1_HS980.trans.data.ob.umap.rds"))

#' for  HS980 before and after transplantation
data.BA.trans.umap <-  readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.ob.umap.rds"))


#' loading late vs early beta cells
H1.DEG.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.results.rds"))
H1.pt.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.pt.results.rds"))
H1.gsea.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.gsea.results.rds"))
H1.fgsea.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.MSigR.fgsea.results.rds"))

H1.DEG.results <- readRDS(paste0("tmp_data/",TD,"/H1.DEG.cluster.results.rds"))



###
plot.results <- list()


### extra heavy loading
data.H1.ob <-  readRDS(paste0("tmp_data/",TD,"/sc.pan.H1.data.ob.rds"))
PATH_list <- list()
for (tp in c("KEGG_2021_Human")) {
  temp <- read.delim(paste("~/Genome_new/enrichr_library/modified/",tp,".mod.txt",sep=""),head=F,stringsAsFactors=F)
  colnames(temp) <- c("TERM","GENE")
  #rownames(temp) <- temp$TERM
  PATH_list[[tp]] <- temp
}



#' for H1 cells
#' plot raw clustre
temp.umap <- data.H1.ob.umap %>% select(cell,EML,seurat_clusters,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2) %>%mutate(EML=seurat_clusters)
temp.text.pos <- temp.umap %>% group_by(EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) #%>% mutate(UMAP_2=ifelse(rename_EML=="PE",UMAP_2+0.5,UMAP_2))
plot.results$UMAP.H1.SC <- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.25,)+geom_text(data=temp.text.pos,mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=4.5 ) + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+NoLegend()+xlim(-6,8.5)+ylim(-7.25,5.5)
plot.results$UMAP.H1.SC



#' plot EML
temp.umap <- data.H1.ob.umap %>% select(cell,EML,seurat_clusters,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2)
temp.text.pos <- temp.umap %>% group_by(EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) #%>% mutate(UMAP_2=ifelse(rename_EML=="PE",UMAP_2+0.5,UMAP_2))
plot.results$UMAP.H1.EML <- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.25,)+geom_text(data=temp.text.pos,mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=4.5 ) + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+NoLegend()+scale_color_manual(values=EML.lineage.col.set)+xlim(-6,8.5)+ylim(-7.25,5.5)
plot.results$UMAP.H1.EML
plot.results$UMAP.H1.EML.legend<- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=1.5)+geom_text(data=temp.text.pos,mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=4 ) + theme_classic()+NoAxes()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=EML.lineage.col.set)
plot.results$UMAP.H1.EML.legend <- plot.results$UMAP.H1.EML.legend %>% ggpubr::get_legend() %>% ggpubr::as_ggplot()
plot.results$UMAP.H1.EML.legend

#temp.hl.genes <- c("INS","HADH","DLK1","ASCL1","GNAS","BACE2","IAPP","PCDH7","CDH8","PRSS23","PDX1","CADM1","CD81","CD99","ARX","ISL1","GCG","GC","CLU","GLS","MKI67","CDK1","TOP2A","PBK","CCNB2","CCNA2")

# shown marker expression 
temp.mk <- H1.fm.mk$separate %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup()

temp.mk <- temp.mk %>% bind_rows(H1.fm.mk$merge %>% filter(cluster=="beta") %>% filter(p_val_adj < 0.05) %>% filter(!gene %in% temp.mk$gene) %>% group_by(cluster) %>% top_n(50,-1*p_val_adj) %>% ungroup())


temp.mk.od <- c("early_beta","late_beta","beta","alpha","delta","SCEC","prolif")
temp.mk.sel <-  ((temp.mk %>% group_by(cluster) %>% top_n(15,-1*p_val_adj) %>% ungroup()) %>% split(.,.$cluster))[temp.mk.od] %>% lapply(function(x){x %>% arrange(p_val_adj) %>% pull(gene)})

temp.M <- data.H1.ob.umap %>% select(cell,EML) %>% mutate(od=factor(EML,EML.od,ordered = T)) %>% arrange(od) %>% select(-od)
temp.exp <- data.H1.ob@assays$RNA$data[unlist(temp.mk.sel) ,temp.M$cell]
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
zs.limit <- 2.5
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
temp.anno <- temp.M  %>% tibble::column_to_rownames("cell")
plot.results$H1.mk.ph <- pheatmap::pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,show_rownames=T,color=heat.col, fontsize_row=4,annotation_colors =list(EML=EML.lineage.col.set[unique(temp.anno$EML)]), border_color="NA",gaps_row =unlist(lapply( temp.mk.sel,function(x){return(length(x))})) %>% cumsum(),main="H1 markers") %>%ggplotify::as.ggplot()


#' Dot plot 
temp.lm.od <- c("early_beta","late_beta","beta","alpha","delta","polyhormonal","SCEC","prolif")
temp.M <-  data.H1.ob.umap  %>% mutate(cluster_EML=EML) %>% filter(cluster_EML %in% temp.lm.od) %>% mutate(od=factor(cluster_EML,temp.lm.od,ordered = T)) %>% arrange(od)
temp.mk.sel <- temp.mk %>% filter(gene %in% c("INS","HADH","ASCL1","BACE2","IAPP","PCDH7","PDX1","NKX6-1","ISL1","GCG","ARX","IRX2","HHEX","SST","SUCNR1","TPH1","COL5A2","FEV","MKI67","CDK1","TOP2A")) %>% mutate(od=factor(cluster,temp.lm.od,ordered = T)) %>% arrange(od) %>% select(-od)
#c("INS","HADH","DLK1","ASCL1","GNAS","BACE2","IAPP","PCDH7","CDH8","PRSS23","PDX1","CADM1","CD81","ARX","ISL1","GCG","GC","CLU","GLS","MKI67","CDK1","TOP2A","PBK","CCNB2","CCNA2")
nrow(temp.mk.sel)

temp.mk.sel$gene %>% duplicated()%>% table()

temp <- data.H1.ob@assays$RNA$data[temp.mk.sel$gene,temp.M$cell] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% gather(cell,logExp,-gene) %>% inner_join(temp.M %>% select(cell,cluster_EML) %>% dplyr::rename(cell.cluster=cluster_EML),by="cell") %>% inner_join(temp.mk.sel %>% select(gene,cluster) %>% dplyr::rename(gene.cluster=cluster) ,by="gene")
temp.input <- temp %>% group_by(cell.cluster,gene.cluster,gene) %>% summarise(meanlogExp=mean(logExp),nCell=n_distinct(cell)) %>% left_join(temp %>% group_by(cell.cluster,gene.cluster,gene) %>% filter(logExp >0)%>% summarise(nExpCell=n_distinct(cell)),by=c("cell.cluster","gene.cluster","gene"))  %>% mutate(nExpCell=ifelse(is.na(nExpCell),0,nExpCell))%>% mutate(nExp.ct=nExpCell/nCell) %>% group_by(gene) %>% mutate(sv=(meanlogExp-mean(meanlogExp))/sd(meanlogExp)) %>% mutate(sv=ifelse(sv> 1.75,1.75,sv)) %>% mutate(sv=ifelse(sv< -1.75,-1.75,sv)) %>% mutate(nExp.ct= as.vector(cut(nExp.ct*100,c(0,5,25,50,75,100),label=c(0,25,50,75,100)))) %>% replace(.,is.na(.),"0") %>% mutate(nExp.ct=as.numeric(nExp.ct))


plot.results$H1.exp.dot <- temp.input %>% mutate(gene=factor(gene,(temp.mk.sel$gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,temp.lm.od,ordered = T))%>% mutate(gene.cluster=factor(gene.cluster,temp.lm.od,ordered = T)) %>%  ggplot(aes(y=cell.cluster, x = gene, color =sv, size = nExp.ct)) + geom_point() + scale_color_gradient2(low="lightgrey",mid="#9595DF",high="blue") + cowplot::theme_cowplot() +theme(axis.text.x = element_text(angle =90, vjust = 1, hjust=1)) +ylab('')  +scale_size(name="",breaks = c(0,25,50,75,100),range = c(0.01, 4))+xlab("")+theme(axis.text.y = element_text(size = 9))+theme(axis.text.x = element_text(size = 7))+ theme(legend.position="top")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())
plot.results$H1.exp.dot 

#' featurePlot 
temp.ft.genes <- c("HADH","BACE2","INS","GCG","SST","TPH1","TOP2A","CHGA","PPY","GPC5-AS1","THSD7A","GHRL","ACSL1","NPY1R","PTF1A","CPA1","DES", "PDGFRB","PDGFRA","PECAM1","CD34","VWF","TH","CHAT","SLC18A3")

#DESMIN, PDGFRb, PDGFRa), endothelial (PECAM-1, CD34, vWF) and neuron (TH, ChAT, VACHT
temp.plot <- list()
data.temp <- data.H1.ob
data.temp@meta.data$SS <- "SS"
for (g in temp.ft.genes ) {
  if (g %in% rownames(data.H1.ob@assays$RNA$counts)) {
    temp.plot[[g]] <- FeaturePlot(data.temp,g,pt.size = 0.1)+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoLegend()+NoAxes()
  }else{
    temp.plot[[g]] <- DimPlot(data.temp,group.by = "SS",cols = "lightgrey",pt.size = 0.1)+ggtitle(paste(g))+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoAxes()+NoLegend()
  }
}
plot.results$H1.ft.plot <- temp.plot
plot.results$H1.ft.plot.legend <- FeaturePlot(data.H1.ob,"TOP2A")+ggtitle("")



#' check the late vs early beta cells
sel.pathway <- c("Diabetic cardiomyopathy","Maturity onset diabetes of the young","Insulin secretion","Thermogenesis","Oxidative phosphorylation","Ferroptosis","Maturity onset diabetes of the young","Carbohydrate digestion and absorption","Insulin resistance","Type II diabetes mellitus","AMPK signaling pathway","cAMP signaling pathway","Cholesterol metabolism","PPAR signaling pathway","Protein processing in endoplasmic reticulum")
sel.genes <- c("IAPP","CPE","MAFA","SIX2","LDHA","NEUROD1","PAX6","PDX1","INS","G6PC2","IGF1R","ADCYAP1","ADCY2","RAP1B","RAP1A","APOC3","P4HB","ASPH", "ASCL1","BACE2", "PCDH7","CDH8","CACNA2D1") #"HOPX","IGF2","NEFM","LDHA"

 #c("HOPX","UCN","IAPP","CPE","SIX3","BACE2","MAFA","FXYD2")

#' plot late beta vs early beta genes scatter 
n <- "late_beta_vs_early_beta"
temp <- H1.DEG.results[[n]]$DEG.all.results %>% mutate(UpDown=ifelse(p_val_adj < 0.05 & avg_log2FC > 0.25 , "UpRe","notDEG"))%>% mutate(UpDown=ifelse(p_val_adj < 0.05 & avg_log2FC < -0.25 , "DownRe",UpDown)) %>% mutate(p_val_adj=ifelse(p_val_adj < 10^-40,10^-40,p_val_adj)) %>% arrange(avg_log2FC)

temp.s1 <- "late_beta"
temp.s2 <- "early_beta"

plot.results$DEG.EL.beta.scatter <- ggplot()+geom_point(temp,mapping=aes(x=avg_log2FC,y=-log10(p_val_adj),col=UpDown),size=0.5)+ggrepel::geom_label_repel(temp %>% filter(gene %in% sel.genes & UpDown %in% c("UpRe","DownRe")),mapping=aes(x=avg_log2FC,y=-log10(p_val_adj),col=UpDown,label=gene),fill="#ffffff33",label.size=NA,max.overlaps=25)+theme_classic()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=c("UpRe"=as.vector(EML.lineage.col.set[temp.s1]),"DownRe"=as.vector(EML.lineage.col.set[temp.s2]), "notDEG"="lightgrey"))+xlim(-8.5,8.5)+geom_vline(xintercept = -0.25 ,linetype="dashed")+geom_vline(xintercept = 0.25 ,linetype="dashed")+geom_hline(yintercept = -log10(0.05) ,linetype="dashed")+NoLegend()
plot.results$DEG.EL.beta.scatter
plot.results$DEG.EL.beta.scatter.legend <-  ggplot()+geom_point(temp,mapping=aes(x=avg_log2FC,y=-log10(p_val_adj),col=UpDown),size=5)+theme_classic()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=c("UpRe"=as.vector(EML.lineage.col.set[temp.s1]),"DownRe"=as.vector(EML.lineage.col.set[temp.s2]), "notDEG"="lightgrey"))+xlim(-8.5,8.5)+geom_vline(xintercept = -0.25 ,linetype="dashed")+geom_vline(xintercept = 0.25 ,linetype="dashed")+geom_hline(yintercept = -log10(0.05) ,linetype="dashed")
plot.results$DEG.EL.beta.scatter.legend <-plot.results$DEG.EL.beta.scatter.legend %>% ggpubr::get_legend() %>% ggpubr::as_ggplot()
plot.results$DEG.EL.beta.scatter.legend

#' #### Pathway enrichment analysis
temp.compair <- "late_beta_vs_early_beta"
temp.all.DEG <- H1.DEG.results[[temp.compair]]$DEG.all.results
temp.fc <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(avg_log2FC)
names(temp.fc) <-  temp.all.DEG %>% arrange(desc(avg_log2FC)) %>% pull(gene)

#temp.sc <- c("Glucagon signaling pathway","Insulin resistance","MAPK signaling pathway","FoxO signaling pathway","AMPK signaling pathway","JAK-STAT signaling pathway","PPAR signaling pathway","Oxidative phosphorylation","Thermogenesis")
#temp.sc <- c("FoxO signaling pathway","TGF-beta signaling pathway")
temp.sc <- c("Hippo signaling pathway","TGF-beta signaling pathway")

temp.fgsea.out <- H1.fgsea.results[[temp.compair]]$KEGG_2021_Human
temp.sel.tp.path <- PATH_list[[tp]] %>% split(.,.$TERM) %>% lapply(function(x){x$GENE})

temp.plot <- list()
for (temp.sel.id in c(temp.sc)) {
  temp.SID <- temp.fgsea.out %>% as.data.frame() %>% tbl_df() %>% filter(pathway==temp.sel.id)  %>% mutate(SID=paste0("NES:",round(NES,3),",","Pvalue:",round(pval,3))) %>% pull(SID)
  temp.plot[[temp.sel.id]] <- fgsea::plotEnrichment(temp.sel.tp.path[[temp.sel.id]],temp.fc) + ggtitle(paste(temp.sel.id,temp.SID,"",sep="\n"))+theme(plot.title = element_text(hjust =0.5))
}
cowplot::plot_grid(plotlist = temp.plot)
plot.results$DEG.H1.EL.beta.fgsea=temp.plot


#' output DEG tables
H1.DEG.results$late_beta_vs_early_beta$DEG.result.up %>% mutate(UpDown="Up-Re") %>% bind_rows(H1.DEG.results$late_beta_vs_early_beta$DEG.result.down %>% mutate(UpDown="Down-Re") ) %>% mutate(Comparison="Late_beta_vs_early_beta") %>% write.table("tmp_data/late_beta_vs_early_beta/DEG",quote=F,sep="\t",col.names = T,row.names = F)

H1.pt.results$late_beta_vs_early_beta$KEGG_2021_Human$up %>% as.data.frame() %>% tbl_df() %>% mutate(UpDown="Up-Re") %>% bind_rows(H1.pt.results$late_beta_vs_early_beta$KEGG_2021_Human$down %>% as.data.frame() %>% tbl_df() %>% mutate(UpDown="Down-Re")) %>% select(ID,GeneRatio,pvalue,geneID,Count,UpDown)  %>% mutate(Comparison="Late_beta_vs_early_beta") %>% filter(pvalue < 0.05) %>% rename(Pathway=ID)%>% write.table("tmp_data/late_beta_vs_early_beta/enriched.KEGG",quote=F,sep="\t",col.names = T,row.names = F)

H1.fgsea.results$late_beta_vs_early_beta$KEGG_2021_Human%>% as.data.frame() %>% tbl_df() %>% mutate(Pathway=pathway,pvalue=pval) %>% select(Pathway,pvalue,NES)  %>% filter(pvalue < 0.05)  %>% mutate(Comparison="Late_beta_vs_early_beta")%>% write.table("tmp_data/late_beta_vs_early_beta/GSEA",quote=F,sep="\t",col.names = T,row.names = F)
#perl  ~/PC/code/txt2excel-1.0.pl -f DEG,enriched.KEGG,GSEA -X sup.lateVSearlyBeta.xls


#' check the UMAP for DP integration
temp.umap <- data.DP.ob.umap %>% select(cell,EML,cluster_EML,SC,prolifSig,pj,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2)  %>% mutate(EML=cluster_EML) 

temp.text.pos <- temp.umap %>% group_by(EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) #%>% mutate(UMAP_2=ifelse(rename_EML=="PE",UMAP_2+0.5,UMAP_2))
temp.maskTable <- generateMask( dims=as.data.frame(temp.umap[,c("UMAP_1","UMAP_2")]), cluster=recode(temp.umap$EML,"late_beta"="beta","early_beta"="beta"), minDensity = 5,smoothSigma = 0.05)

plot.results$UMAP.DP.EML <-ggplot()+geom_point(temp.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.1,alpha=0.75)+geom_path( temp.maskTable,mapping=aes(x=UMAP_1,y=UMAP_2,group=group),linewidth=0.5,linetype = 2)+geom_text(data=temp.text.pos,mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=4.5 )  + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+NoLegend()+scale_color_manual(values=EML.lineage.col.set)+xlim(-6,12)+ylim(-8,9)
plot.results$UMAP.DP.EML
plot.results$UMAP.DP.EML.legend<- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=1.5)+geom_text(data=temp.text.pos,mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=4 ) + theme_classic()+NoAxes()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=EML.lineage.col.set)
plot.results$UMAP.DP.EML.legend <- plot.results$UMAP.DP.EML.legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()
plot.results$UMAP.DP.EML.legend

#' split DP plot
temp.umap <- data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% select(cell,EML,cluster_EML,SC,prolifSig,pj,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2) %>% mutate(EML=cluster_EML)

plot.results$UMAP.DP.EML.split <- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.03,alpha=0.75) + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+NoLegend()+scale_color_manual(values=EML.lineage.col.set)+xlim(-6,12)+ylim(-8,9)+facet_wrap(~factor(pj,pj.od,ordered = T),ncol=7)
plot.results$UMAP.DP.EML.split 

# ' highlight prolif
# temp.umap <- data.DP.ob.umap %>% select(cell,EML,cluster_EML,SC,prolifSig,pj,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2) %>% mutate(cluster_EML=recode(cluster_EML,"late_beta"="beta","early_beta"="beta")) %>% mutate(EML=cluster_EML)
# plot.results$UMAP.DP.HL.prolif <- ggplot()+geom_point(temp.umap %>% filter(prolifSig > 0) ,mapping=aes(x=UMAP_1,y=UMAP_2),size=0.03,color="firebrick3") +geom_point(temp.umap %>% filter(prolifSig <= 0) ,mapping=aes(x=UMAP_1,y=UMAP_2),size=0.05,alpha=0.75,color="lightgrey")+ theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+NoLegend()+xlim(-6,12)+ylim(-8,9)+facet_wrap(~factor(pj,pj.od,ordered = T),ncol=5)
# plot.results$UMAP.DP.HL.prolif

# ' highlight prolif ft
temp.umap <- data.DP.ob.umap%>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% select(cell,EML,cluster_EML,SC,prolifSig,pj,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2) %>% mutate(cluster_EML=recode(cluster_EML,"late_beta"="beta","early_beta"="beta")) %>% mutate(EML=cluster_EML) %>% arrange(prolifSig)
plot.results$UMAP.DP.ft.prolif <- ggplot()+geom_point(temp.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,col=prolifSig),size=0.03) + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+NoLegend()+xlim(-6,12)+ylim(-8,9)+facet_wrap(~factor(pj,pj.od,ordered = T),ncol=7)+scale_color_gradient2(low="lightgrey",mid="lightgrey",high="blue")
plot.results$UMAP.DP.ft.prolif.legend <- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,col=prolifSig),size=0.03) + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+xlim(-6,12)+ylim(-8,9)+facet_wrap(~factor(pj,pj.od,ordered = T),ncol=4)+scale_color_gradient2(low="lightgrey",mid="lightgrey",high="blue")
plot.results$UMAP.DP.ft.prolif.legend <- plot.results$UMAP.DP.ft.prolif.legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()
plot.results$UMAP.DP.ft.prolif.legend

#' check the prolif percentage 
data.DP.ob.umap  %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj))%>% filter(prolifSig > 0) %>% group_by(pj) %>% summarise(pp_Cell=n_distinct(cell)) %>% inner_join(data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% group_by(pj) %>% summarise(total_Cell=n_distinct(cell)),by="pj") %>% mutate(prop=round(pp_Cell/total_Cell*100,2)) 


#data.DP.ob.umap  %>% mutate(pj=ifelse(pj=="Rajaei_2025" ,devTime,pj))%>% filter(prolifSig > 0 & cluster_EML!="prolif") %>% group_by(pj) %>% summarise(pp_Cell=n_distinct(cell)) %>% inner_join(data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% group_by(pj) %>% summarise(total_Cell=n_distinct(cell)),by="pj") %>% mutate(prop=round(pp_Cell/total_Cell*100,2)) 


data.DP.ob.umap  %>% mutate(pj=ifelse(pj=="Rajaei_2025" ,devTime,pj))%>% filter(prolifSig > 0 & cluster_EML=="exo") %>% group_by(pj) %>% summarise(pp_Cell=n_distinct(cell)) %>% inner_join(data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% group_by(pj) %>% summarise(total_Cell=n_distinct(cell)),by="pj") %>% mutate(prop=round(pp_Cell/total_Cell*100,2)) 

data.DP.ob.umap  %>% mutate(pj=ifelse(pj=="Rajaei_2025" ,devTime,pj))%>% filter(prolifSig > 0 & cluster_EML %in% c("alpha","delta","early_beta","late_beta","NeuroEndo","polyhormonal","prolif","SCEC")) %>% group_by(pj) %>% summarise(pp_Cell=n_distinct(cell)) %>% inner_join(data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% group_by(pj) %>% summarise(total_Cell=n_distinct(cell)),by="pj") %>% mutate(prop=round(pp_Cell/total_Cell*100,2)) 


#' show the barplot number
temp.umap <- data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% select(cell,EML,cluster_EML,SC,prolifSig,pj,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2) %>% mutate(EML=cluster_EML) %>% mutate(cluster_EML=recode(cluster_EML,"late_beta"="beta","early_beta"="beta"))
temp.input1 <- temp.umap %>% group_by(pj,cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) %>% select(-nCell) %>% spread(cluster_EML,prop)%>% replace(.,is.na(.),0) %>% gather(cluster_EML,prop,-pj) %>% mutate(pj=factor(pj,pj.od,ordered = T)) 
temp.input2 <- temp.umap %>% group_by(pj,EML,cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) %>% unite(cluster_EML,c(EML,cluster_EML),sep=":") %>% select(-nCell) %>% spread(cluster_EML,prop)%>% replace(.,is.na(.),0) %>% gather(cluster_EML,prop,-pj) %>% separate(cluster_EML,c("EML","cluster_EML"),sep=":") %>% mutate(pj=factor(pj,pj.od,ordered = T)) 
plot.results$H1.bar.left <- ggplot()+geom_bar(data=temp.input1  %>% filter(cluster_EML %in% c("beta","alpha","SCEC"))  %>% mutate(cluster_EML=factor(cluster_EML,c("beta","alpha","SCEC"),ordered = T)) , mapping=aes(x=cluster_EML,y=prop,fill=pj),stat="identity",position="dodge",alpha=0.75,width=0.75,col="black")+geom_bar(data=temp.input2 %>% filter(EML %in% c("early_beta","alpha","SCEC")) %>% mutate(cluster_EML=factor(cluster_EML,c("beta","alpha","SCEC"),ordered = T)),mapping=aes(x=cluster_EML,y=prop,fill=pj),stat="identity",position="dodge",width=0.75,col="grey33")+scale_fill_manual(values=pj.col.set)+xlab("")+theme_classic()+ylim(0,0.75)+ylab("")
plot.results$H1.bar.left

plot.results$H1.bar.right <- ggplot()+geom_bar(data=temp.input1 %>% filter(cluster_EML %in% c("delta","polyhormonal","NeuroEndo","prolif","exo"))  %>% mutate(cluster_EML=factor(cluster_EML,c("delta","polyhormonal","NeuroEndo","prolif","exo"),ordered = T)) , mapping=aes(x=cluster_EML,y=prop,fill=pj),stat="identity",position="dodge",width=0.75,col="black")+scale_fill_manual(values=pj.col.set)+xlab("")+theme_classic()+ylim(0,0.41)+ylab("")+ggbreak::scale_y_break(c(0.18, 0.35))
plot.results$H1.bar.right

data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% group_by(cluster_EML,pj) %>% summarise(nCell=n_distinct(cell))%>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) %>% filter(cluster_EML=="SCEC")

data.DP.ob.umap %>% group_by(cluster_EML,pj) %>% summarise(nCell=n_distinct(cell))%>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) %>% filter(cluster_EML=="SCEC")



data.DP.ob.umap %>% mutate(pj=ifelse(pj=="Rajaei_2025",devTime,pj)) %>% group_by(cluster_EML,pj) %>% summarise(nCell=n_distinct(cell))%>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) %>% filter(cluster_EML=="exo")

data.DP.ob.umap %>% group_by(cluster_EML,pj) %>% summarise(nCell=n_distinct(cell))%>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) %>% filter(cluster_EML=="exo")
temp.umap %>%select(-SC)  %>% rename(annotation=EML,cluster_annotation=cluster_EML)%>% write.table("tmp_data/temp.source.data.fig4G_J.tsv",col.names = T,row.names = F,quote = F,sep="\t")

#' show DP integration based on orignal annotation
temp.umap <- data.DP.ob.umap %>% filter(pj=="H1") %>% bind_rows(data.DP.ob.umap %>% filter(pj=="Veres_2019") %>% mutate(EML=ifelse(EML=="exo",subCT,EML)) %>%mutate(EML=ifelse(subCT=="other__gap43","NeuroEndo",EML)) %>% mutate(EML=recode(EML,"neurog3"="other","acinar_like"="acinar","ductal_like"="ductal","early_exo"="exo","late_exo"="exo") ))  %>% bind_rows(data.DP.ob.umap %>% filter(pj=="Augsor_2022") %>% mutate(EML=recode(EML,"prolif_alpha"="prolif")) ) %>% bind_rows(data.DP.ob.umap %>% filter(pj=="BalBoa_2022") %>% mutate(EML=recode(EML,"prolif_alpha"="prolif","not_endo"="exo")))  %>% bind_rows(data.DP.ob.umap %>% filter(pj=="Rajaei_2025") %>% mutate(pj=devTime) %>% mutate(EML=recode(EML,"early_duct"="ductal","duct"="ductal","GAP43pNE"="NeuroEndo")))  %>% rename(UMAP_1=umap_1,UMAP_2=umap_2) 

temp.plot <- list()
temp.plot.legend <- list()
for ( p in c("H1","Augsor_2022","BalBoa_2022","Veres_2019","Rajaei_not_enriched","Rajaei_enriched")) {
  temp.plot[[p]] <- ggplot()+geom_point(temp.umap %>% filter(pj==p),mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.2,alpha=0.75) + theme_classic()+ggtitle(p)+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=c(EML.lineage.col.set,EML.lineage.col.extra.set))+xlim(-6,12)+ylim(-8,9)+NoLegend()
  temp.plot.legend [[p]] <- ggplot()+geom_point(temp.umap %>% filter(pj==p),mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=2) + theme_classic()+ggtitle(p)+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=c(EML.lineage.col.set,EML.lineage.col.extra.set))+xlim(-6,12)+ylim(-8,9)
  temp.plot.legend [[p]] <- temp.plot.legend [[p]]%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()
}
cowplot::plot_grid(plotlist = temp.plot)

plot.results$UMAP.DP.raw.EML.split <- temp.plot
plot.results$UMAP.DP.raw.EML.split.legend <- temp.plot.legend




#' H1 HS980
temp.umap <- data.H1_H9S980.umap %>% select(cell,EML,pj,umap_1,umap_2,pj) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2)  
temp.text.pos <- temp.umap %>% group_by(EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) #%>% mutate(UMAP_2=ifelse(rename_EML=="PE",UMAP_2+0.5,UMAP_2))
plot.results$UMAP.H1980.EML <-ggplot()+geom_point(temp.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.3,alpha=0.75)+geom_text(data=temp.text.pos,mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=4.5 )  + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=EML.lineage.col.set)
plot.results$UMAP.H1980.EML

plot.results$UMAP.H1980.pj <- ggplot()+geom_point(temp.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=pj),size=0.3,alpha=0.75) + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=pj.col.set)
plot.results$UMAP.H1980.pj


#' show the barplot number
temp.input1 <- data.H1_H9S980.umap %>% filter(pj=="HS980_notrans_CM310")  %>% mutate(cluster_EML=recode(EML,"late_beta"="beta","early_beta"="beta")) %>% group_by(pj,cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) # %>% mutate(EML=factor(EML, EML.od,ordered = T)) 

temp.input2 <- data.H1_H9S980.umap %>% filter(pj=="HS980_notrans_CM310")  %>% mutate(cluster_EML=recode(EML,"late_beta"="beta","early_beta"="beta")) %>% group_by(pj,cluster_EML,EML) %>% summarise(nCell=n_distinct(cell)) %>% group_by(pj) %>% mutate(prop=nCell/sum(nCell)) 


plot.results$HS980.bar <- ggplot()+geom_bar(data=temp.input1 %>% mutate(cluster_EML=factor(cluster_EML,c("beta",EML.od[3:7]),ordered = T)) , mapping=aes(x=cluster_EML,y=prop*100,fill=pj),stat="identity",position="dodge",alpha=0.75,width=0.75,col="black")+geom_bar(data=temp.input2 %>% filter(!EML %in% c("late_beta")) %>% mutate(cluster_EML=factor(cluster_EML,c("beta",EML.od[3:7]),ordered = T)),mapping=aes(x=cluster_EML,y=prop*100,fill=pj),stat="identity",position="dodge",width=0.75,col="grey33")+scale_fill_manual(values=pj.col.set)+xlab("")+theme_classic()+ylab("Proportion(%)")+NoLegend()+ggtitle("HS980")+FunTitle()


#' shown polyhormonal difference
data.temp <- JoinLayers(readRDS(paste0("tmp_data/",TD,"/H1_HS980.trans.data.ob.rds")))
temp.input <- data.temp@assays$RNA$data[c("INS","GCG"),] %>% as.data.frame()%>% tibble::rownames_to_column("gene") %>% tbl_df()%>% gather(cell,logExp,-gene) %>% spread(gene,logExp)%>% inner_join(data.H1_H9S980.umap,by="cell")
rm(data.temp)

plot.results$HS980.co.check <- temp.input %>% filter(pj=="HS980_notrans_CM310") %>% ggplot()+geom_point(mapping=aes(x=INS,y=GCG,color=EML),size=0.75)+theme_classic()+scale_color_manual(values=EML.lineage.col.set)
plot.results$HS980.co.check

plot.results$H1.co.check <- temp.input %>% filter(pj=="H1") %>% ggplot()+geom_point(mapping=aes(x=INS,y=GCG,color=EML),size=0.75)+theme_classic()+scale_color_manual(values=EML.lineage.col.set)
plot.results$H1.co.check


#' shown before and after trans results
temp.umap <- data.BA.trans.umap  %>% select(cell,devTime,EML,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2)
temp.text.pos <- temp.umap %>% group_by(EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) #%>% mutate(UMAP_2=ifelse(rename_EML=="PE",UMAP_2+0.5,UMAP_2))
plot.results$UMAP.BA.EML <- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.25,alpha=0.75)+geom_text(data=temp.text.pos,mapping=aes(x=UMAP_1,y=UMAP_2,label=EML),size=4.5 ) + theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=EML.lineage.col.set)+xlim(-5,13)+ylim(-7.5,6)
plot.results$UMAP.BA.EML

plot.results$UMAP.BA.EML.split <- ggplot()+geom_point(temp.umap  ,mapping=aes(x=UMAP_1,y=UMAP_2,color=EML),size=0.25,alpha=0.75)+ theme_classic()+ggtitle("Annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=EML.lineage.col.set)+NoLegend()+facet_wrap(~devTime)
plot.results$UMAP.BA.EML.split+xlim(-5,13)+ylim(-7.5,6)

# check markers
if (generate_dot_plot) {
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))[,meta.filter$cell]
  BA.lognormExp.mBN <- readRDS(file=paste0("tmp_data/",TD,"/before.after.trans.fastMNN.lognormExp.mBN.rds"))
  

  for (sa in c("HS980_notrans_CM310","HS980_trans")) {
    temp.M <- data.BA.trans.umap %>% filter(devTime %in% sa) %>% select(cell:mt.perc)
    temp.sel.expG <- rownames(BA.lognormExp.mBN)
    
    data.temp <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
    data.temp@assays$RNA$data <- as.matrix(BA.lognormExp.mBN[temp.sel.expG,rownames(data.temp@meta.data)])
    data.temp <- data.temp  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:25,verbose=F)
    data.temp@reductions$umap@cell.embeddings[,1] <- (data.BA.trans.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"umap_1"]
    data.temp@reductions$umap@cell.embeddings[,2] <- (data.BA.trans.umap  %>% as.data.frame()%>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"umap_2"]
    Idents(data.temp) <- factor(data.temp@meta.data$EML)
    data.temp@meta.data$SS <- "SS"
    #temp.ft.genes <- c("HADH","IAPP","INS","GCG","SST","TPH1","TOP2A","CHGA","PPY","GPC5-AS1","THSD7A","GHRL","ACSL1","NPY1R","PTF1A","CPA1","DCN","LGALS1","MAFA","SIX3","UCN3")#,"CALB2","BACE2",
    #temp.ft.genes <- c("HADH","BACE2","INS","GCG","SST","TPH1","TOP2A","CHGA","PPY","GPC5-AS1","THSD7A","GHRL","ACSL1","NPY1R","PTF1A","CPA1","DCN","LGALS1","MAFA","SIX3","UCN3","DES", "PDGFRB","PDGFRA","PECAM1","CD34","VWF","TH","CHAT","SLC18A3")
    temp.ft.genes <- c("HADH","BACE2","INS","GCG","SST","TPH1","TOP2A","CHGA","PPY","GPC5-AS1","THSD7A","GHRL","ACSL1","NPY1R","PTF1A","CPA1","DCN","LGALS1","MAFA","FXYD2","G6PC2","SCGN","PCSK1N","RBP4","DES", "PDGFRB","PDGFRA","PECAM1","CD34","VWF","TH","CHAT","SLC18A3","IAPP","FXYD2","G6PC2","SCGN","PCSK1N","RBP4","NEUROG3")
    temp.ft.genes %>% setdiff(rownames( counts.filter))
    temp.plot <- list()
    for (g in temp.ft.genes ) {
      if (!g %in% rownames(data.temp@assays$RNA$counts) & g %in% rownames(counts.filter)) {
        print(g)
        temp.plot[[g]] <- DimPlot(data.temp,group.by = "SS",cols = "lightgrey",pt.size = 0.1)+ggtitle(paste(g))+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoAxes()+NoLegend()+xlim(-5,13)+ylim(-7.5,6)
          ggplot()+theme_void()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="plain"))
      }else if ( sum(data.temp@assays$RNA$counts[g,])==0) {
        temp.plot[[g]] <- FeaturePlot(data.temp,g,pt.size = 0.001,cols =c("lightgrey", "lightgrey"))+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoLegend()+NoAxes()+xlim(-5,13)+ylim(-7.5,6)
      }else{
        temp.plot[[g]] <- FeaturePlot(data.temp,g,pt.size = 0.001)+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoLegend()+NoAxes()+xlim(-5,13)+ylim(-7.5,6)
        
      }
    }
    cowplot::plot_grid(plotlist = temp.plot)
    plot.results$HS980.BA.ft.plot[[sa]] <- temp.plot
  }
  #' featurePlot of mature markers
  temp.M <- data.BA.trans.umap  %>% select(cell:mt.perc)
  temp.sel.expG <- rownames(BA.lognormExp.mBN)
  
  data.temp <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.temp@assays$RNA$data <- as.matrix(BA.lognormExp.mBN[temp.sel.expG,rownames(data.temp@meta.data)])
  data.temp <- data.temp  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:25,verbose=F)
  data.temp@reductions$umap@cell.embeddings[,1] <- (data.BA.trans.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"umap_1"]
  data.temp@reductions$umap@cell.embeddings[,2] <- (data.BA.trans.umap  %>% as.data.frame()%>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"umap_2"]
  Idents(data.temp) <- factor(data.temp@meta.data$EML)
  data.temp@meta.data$SS <- "SS"
  temp.ft.genes <- c("INS", "IAPP","SIX3","MAFA","FXYD2","G6PC2", "SCGN", "PCSK1N", "HADH", "FXYD2", "RBP4")
  
  temp.plot <- list()
  for (g in temp.ft.genes ) {
    temp.plot[[g]] <- FeaturePlot(data.temp,g,pt.size = 0.001,split.by = "devTime",order = T)
  }
  cowplot::plot_grid(plotlist = temp.plot)
  plot.results$HS980.BA.mature.ft.plot <- temp.plot
  
  
  
  
  
  
  
  temp.lm.od <- c("early_beta","late_beta","beta","alpha","delta","polyhormonal","SCEC","prolif","psc")
  temp.M <-  data.BA.trans.umap %>% mutate(cluster_EML=EML) %>% filter(cluster_EML %in% temp.lm.od) %>% mutate(od=factor(cluster_EML,temp.lm.od,ordered = T)) %>% arrange(od)
  temp.mk.sel <- data.frame(gene=c("HADH","INS","ASCL1","BACE2","IAPP","PCDH7","PDX1","ISL1","NKX6-1","ARX","IRX2","GCG","HHEX","SUCNR1","SST","TPH1","COL5A2","FEV","MKI67","CDK1","TOP2A","DCN","LGALS1","CAV1"),cluster=c("early_beta","early_beta","early_beta","late_beta","late_beta","late_beta","beta","beta","beta","alpha","alpha","alpha","delta","delta","delta","SCEC","SCEC","SCEC","prolif","prolif","prolif","psc","psc","psc")) %>% tbl_df()
  
  nrow(temp.mk.sel)
  
  temp.mk.sel$gene %>% duplicated()%>% table()
  
  temp <- BA.lognormExp.mBN[temp.mk.sel$gene,temp.M$cell] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% gather(cell,logExp,-gene) %>% inner_join(temp.M %>% select(cell,devTime,cluster_EML) %>% dplyr::rename(cell.cluster=cluster_EML),by="cell") %>% inner_join(temp.mk.sel %>% select(gene,cluster) %>% dplyr::rename(gene.cluster=cluster) ,by="gene")
  temp.input <- temp %>% group_by(cell.cluster,gene.cluster,gene,devTime) %>% summarise(meanlogExp=mean(logExp),nCell=n_distinct(cell)) %>% left_join(temp %>% group_by(cell.cluster,gene.cluster,gene,devTime) %>% filter(logExp >0)%>% summarise(nExpCell=n_distinct(cell)),by=c("cell.cluster","gene.cluster","gene","devTime"))  %>% mutate(nExpCell=ifelse(is.na(nExpCell),0,nExpCell))%>% mutate(nExp.ct=nExpCell/nCell) %>% group_by(gene,devTime) %>% mutate(sv=(meanlogExp-mean(meanlogExp))/sd(meanlogExp)) %>% mutate(sv=ifelse(sv> 1.75,1.75,sv)) %>% mutate(sv=ifelse(sv< -1.75,-1.75,sv)) %>% mutate(nExp.ct= as.vector(cut(nExp.ct*100,c(0,5,25,50,75,100),label=c(0,25,50,75,100)))) %>% mutate(nExp.ct=ifelse(is.na(nExp.ct),"0",nExp.ct)) %>% mutate(sv=ifelse(is.na(sv),-1.75,sv)) %>% mutate(nExp.ct=as.numeric(nExp.ct)) %>% mutate(devTime=recode(devTime,"HS980_notrans_CM310"="before","HS980_trans"="after"))
  
  plot.results$BA.exp.dot <- temp.input %>% mutate(cell.cluster=paste(cell.cluster,devTime,sep=":"))%>% mutate(gene=factor(gene,(temp.mk.sel$gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,paste(rep(temp.lm.od,each=2),rep(c("before","after"),length(temp.lm.od)),sep=":"),ordered = T))%>% mutate(gene.cluster=factor(gene.cluster,temp.lm.od,ordered = T)) %>%  ggplot(aes(y=cell.cluster, x = gene, color =sv, size = nExp.ct)) + geom_point() + scale_color_gradient2(low="lightgrey",mid="#9595DF",high="blue") + cowplot::theme_cowplot() +theme(axis.text.x = element_text(angle =90, vjust = 1, hjust=1)) +ylab('')  +scale_size(name="",breaks = c(0,25,50,75,100),range = c(0.01, 4))+xlab("")+theme(axis.text.y = element_text(size = 9))+theme(axis.text.x = element_text(size = 7))+ theme(legend.position="top")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())
  plot.results$BA.exp.dot
  
  
  #temp.sel.genes <-c("INS","G6PC2","HOPX","UCN3","IAPP","SIX3","MAFA","CPE","FXYD2","MAFA","SLC30A8","ASCL1","CHGB","SCG2")
  #temp.sel.genes <-c("INS","SST","WNT4","CDH1","IAPP","MAFA","SLC30A8","ASCL1","CHGB","SCG2")
  #temp.sel.genes <- c("HOPX","UCN","IAPP","CPE","SIX3","BACE2","MAFA","FXYD2")
  temp.sel.genes <- c("INS", "IAPP","SIX3","MAFA","FXYD2","G6PC2", "SCGN", "PCSK1N", "HADH", "FXYD2", "RBP4")

  temp.M <- data.BA.trans.umap %>% filter(EML %in% c("early_beta","late_beta"))  %>% mutate(devTime=recode(devTime,"HS980_notrans_CM310"="before","HS980_trans"="after")) %>% mutate(SID=EML) 
  temp.norm <- BA.lognormExp.mBN[temp.sel.genes,temp.M$cell ]
  temp.exp <- temp.norm [temp.sel.genes,] %>% as.data.frame()%>% tibble::rownames_to_column("gene") %>% gather(cell,logExp,-gene) %>% mutate(Exp=expm1(logExp)) %>% inner_join(temp.M,by="cell") %>% group_by(gene,SID,devTime) %>% summarise(mean=mean(Exp),sd=sd(Exp),n=n_distinct(cell))  %>% mutate(sem=sd/(n^0.5))  %>% mutate(devTime=factor(devTime,c("before","after"),ordered = T))
  
  temp.plot <- list()
  for (g in temp.sel.genes) {
    if (g %in% c("INS","IAPP")) {
      temp.plot[[g]] <- temp.exp %>% filter(gene %in% g)%>%  ggplot(mapping=aes(x=SID,y=mean,fill=devTime))+geom_bar(stat="identity",position="dodge")+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),size=0.25,width=0.5,position=position_dodge(1))+theme_classic()+xlab("")+ylab("Normalized expression")+ggtitle(g)+theme(plot.title = element_text(hjust=0.5))+ theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=pj.col.set)+NoLegend()+ylab("")+scale_y_log10()
    }else{
      temp.plot[[g]] <- temp.exp %>% filter(gene %in% g)%>%  ggplot(mapping=aes(x=SID,y=mean,fill=devTime))+geom_bar(stat="identity",position="dodge")+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),size=0.25,width=0.5,position=position_dodge(1))+theme_classic()+xlab("")+ylab("Normalized expression")+ggtitle(g)+theme(plot.title = element_text(hjust=0.5))+ theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=pj.col.set)+NoLegend()+ylab("")
    }
    
  }
  expm1(BA.lognormExp.mBN[c("INS", "IAPP","MAFA","FXYD2","G6PC2", "SCGN", "PCSK1N", "HADH", "RBP4"),temp.M$cell ]) %>% tibble::rownames_to_column("gene")%>%  write.table("tmp_data/temp.source.data.fig6E.tsv",col.names = T,row.names = F,quote = F,sep="\t")
  plot.results$DEG.bar.trans <- temp.plot
  cowplot::plot_grid(plotlist = temp.plot)
  #cowplot::plot_grid(plotlist = temp.plot$before,nrow=1)
  #cowplot::plot_grid(plotlist = temp.plot$after,nrow=1)
}


#' check the mature signal
plot.results$trans.mature.sig <- ggplot(data=data.BA.trans.umap %>% filter(EML %in% c("late_beta")),mapping=aes(x=devTime,y=matureSig,fill=devTime)) + geom_violin() +  theme_cowplot() + geom_signif(comparisons = list(c("HS980_notrans_CM310", "HS980_trans")), map_signif_level = TRUE)+scale_fill_manual(values = pj.col.set) + xlab("")
plot.results$trans.mature.sig

#' plot the bar
temp.umap <- data.BA.trans.umap %>% mutate(cluster_EML=EML) %>% select(cell,EML,cluster_EML,pj,umap_1,umap_2) %>% rename(UMAP_1=umap_1,UMAP_2=umap_2) %>% mutate(EML=cluster_EML) %>% mutate(cluster_EML=recode(cluster_EML,"late_beta"="beta","early_beta"="beta")) %>% mutate(pj=recode(pj,"HS980_notrans_CM310"="before","HS980_trans"="after"))
temp.input1 <- temp.umap %>% group_by(pj,cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% group_by(pj) %>% mutate(prop=nCell/sum(nCell))  %>% mutate(pj=factor(pj,c("before","after"),ordered = T))  %>% mutate(prop=100*prop)
temp.input2 <- temp.umap %>% group_by(pj,EML,cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% group_by(pj) %>% mutate(prop=nCell/sum(nCell))  %>% mutate(pj=factor(pj,c("before","after"),ordered = T))  %>% mutate(prop=100*prop)
plot.results$HS980.BA.bar.left <- ggplot()+geom_bar(data=temp.input1  %>% filter(cluster_EML %in% c("beta","alpha","SCEC"))  %>% mutate(cluster_EML=factor(cluster_EML,c("beta","alpha","SCEC"),ordered = T)) , mapping=aes(x=cluster_EML,y=prop,fill=pj),stat="identity",position="dodge",alpha=0.75,width=0.75,col="black")+geom_bar(data=temp.input2 %>% filter(EML %in% c("early_beta","alpha","SCEC")) %>% mutate(cluster_EML=factor(cluster_EML,c("beta","alpha","SCEC"),ordered = T)),mapping=aes(x=cluster_EML,y=prop,fill=pj),stat="identity",position="dodge",width=0.75,col="grey33")+scale_fill_manual(values=pj.col.set)+xlab("")+theme_classic()+ylim(0,80)+ylab("")
plot.results$HS980.BA.bar.left

plot.results$HS980.BA.bar.right <- ggplot()+geom_bar(data=temp.input1 %>% filter(cluster_EML %in% c("delta","polyhormonal","prolif","psc"))  %>% mutate(cluster_EML=factor(cluster_EML,c("delta","polyhormonal","prolif","psc"),ordered = T)) , mapping=aes(x=cluster_EML,y=prop,fill=pj),stat="identity",position="dodge",width=0.75,col="black")+scale_fill_manual(values=pj.col.set)+xlab("")+theme_classic()+ylim(0,7.5)+ylab("")

temp.umap %>% select(cell,EML,pj,UMAP_1,UMAP_2) %>% rename(annotation=EML) %>% write.table("tmp_data/temp.source.data.fig6A_C.tsv",col.names = T,row.names = F,quote = F,sep="\t")

# check activate /quisent
if (check) {
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))[,meta.filter$cell]
  BA.lognormExp.mBN <- readRDS(file=paste0("tmp_data/",TD,"/before.after.trans.fastMNN.lognormExp.mBN.rds"))
  
  
  sa="HS980_trans"
  temp.M <- data.BA.trans.umap %>% filter(devTime %in% sa) %>% select(cell:mt.perc)
  temp.sel.expG <- rownames(BA.lognormExp.mBN)
  
  data.temp <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.temp@assays$RNA$data <- as.matrix(BA.lognormExp.mBN[temp.sel.expG,rownames(data.temp@meta.data)])
  data.temp <- data.temp  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:25,verbose=F)
  data.temp@reductions$umap@cell.embeddings[,1] <- (data.BA.trans.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"umap_1"]
  data.temp@reductions$umap@cell.embeddings[,2] <- (data.BA.trans.umap  %>% as.data.frame()%>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"umap_2"]
  Idents(data.temp) <- factor(data.temp@meta.data$EML)
  data.temp@meta.data$SS <- "SS"
  temp.ft.genes <- c("DES", "GFAP", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "FN1")
  temp.plot <- list()
  for (g in temp.ft.genes ) {
    if (!g %in% rownames(data.temp@assays$RNA$counts)) {
      temp.plot[[g]] <- DimPlot(data.temp,group.by = "SS",cols = "lightgrey",pt.size = 0.1)+ggtitle(paste(g))+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoAxes()+NoLegend()
      ggplot()+theme_void()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="plain"))
    }else if ( sum(data.temp@assays$RNA$counts[g,])==0) {
      temp.plot[[g]] <- FeaturePlot(data.temp,g,pt.size = 0.001,cols =c("lightgrey", "lightgrey"))+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoLegend()+NoAxes()+xlim(-5,13)+ylim(-7.5,6)
    }else{
      temp.plot[[g]] <- FeaturePlot(data.temp,g,pt.size = 0.001)+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="plain"))+NoLegend()+NoAxes()+xlim(-5,13)+ylim(-7.5,6)
      
    }
  }
  cowplot::plot_grid(plotlist = temp.plot)
  plot.results$stella.qa.plot <- temp.plot  
  temp <- counts.filter[temp.ft.genes ,data.BA.trans.umap %>% filter(devTime %in% sa & EML=="psc") %>% pull(cell)] %>% tibble::rownames_to_column("gene") %>% gather(cell,ct,-gene) %>% tbl_df() %>% mutate(ep=ifelse(ct >0,"Exp","NE")) %>% group_by(gene,ep) %>% summarise(nCell=n_distinct(cell)) %>%group_by(gene) %>% mutate(prop=nCell/sum(nCell)) %>% select(gene,ep,prop) %>% spread(ep,prop) %>% replace(.,is.na(.),0) %>% gather(ep,prop,-gene) %>% mutate(gene=factor(gene,temp.ft.genes,ordered = T))
  plot.results$stella.qa.ep.plot <- temp %>% filter(ep=="Exp") %>% ggplot()+geom_bar(mapping=aes(x=gene,y=prop*100),stat="identity",width=0.6,fill="grey66")+xlab("")+ylab("Proportion of stella cells with expression")+theme_classic()
}
if (check) {
  prolif.mk <- c("MKI67","CDK1","TOP2A","CCNB2","CCNA2","PBK")
  data.ob <- readRDS(paste0("tmp_data/",TD,"/before.after.trans.data.ob.rds"))
  data.temp <- JoinLayers(data.ob)
  data.temp <- AddModuleScore(data.temp, features = list(prolifSig=prolif.mk ), name = "prolifSig")
  temp <- data.frame(cell=rownames(data.temp@meta.data),prolifSig=data.temp@meta.data$prolifSig1) %>% tbl_df() %>% inner_join(data.BA.trans.umap ,by="cell")
  temp %>% filter(pj=="HS980_trans" & EML=="psc") %>% filter(prolifSig >0)
}
