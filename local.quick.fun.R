# 
EML.od <- c("early_beta","late_beta","alpha","delta","SCEC","polyhormonal","prolif")
mk.od <-  c("early_beta","late_beta","beta","alpha","delta","SCEC","polyhormonal","prolif")

pj.od <- c("H1","Veres_2019","Augsor_2022","BalBoa_2022")
pj.col.set <- c(H1="royalblue",Veres_2019="#00B6EB",Augsor_2022="#00ced1", BalBoa_2022="#00C08B","HS980_notrans_CM310"="#F8766D","before"="#F8766D","HS980"="#F8766D" ,HS980_trans="#A58AFF","after"="#A58AFF")


trans.lineage.col.set <- c("Alpha"="#eead0e","Early_Beta"="#ff7256","Mature_Beta"="#8b3e2f","Stellate_like"="#00ced1", "EC_like"="#228b22","Exocrine_progenitor"="#00B9E3","Polyhormonal"="#838b8b","Delta"="#d02090","Proliferating_Endocrine"="#ee00ee")

#BAT.col.set <- c(BT="#00BFC4",AT="#F8766D")
#BAT.EML.set.od <- c("Mature_Beta","Early_Beta","EC_like","Alpha","Polyhormonal","Stellate_like","Delta", "Exocrine_progenitor","Proliferating_Endocrine")
#UpDown.col.set  <- c("DownRe"="#00BFC4","UpRe"="#F8766D")

mature.score.mk <- c("INS", "G6PC2", "HOPX", "UCN3", "IAPP", "CPE", "SIX3", "BACE2", "MAFA", "FXYD2")

EML.lineage.col.set <- c(alpha="#F8766D",delta="#B79F00",early_beta="#53B400",late_beta="#00B6EB",beta="#2AB576",polyhormonal="#DE8C00",prolif="#A58AFF",SCEC="#FB61D7",exo="#C2B59B","psc"="#838b8b")

EML.lineage.col.extra.set <- c("Endo_Prog"="#00BA38","epsilon"="#CC4678","gamma"="#7E03A8","acinar"="#F0F921", "ductal"="#094334", "other"="#00BA38","mesenchymal"="#A4C9DD")



heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)




human.mk.list <- list(alpha=c("GCG","TTR","PCSK2"), beta=c("INS","ENO1","EDARADD"),delta=c("SST","HHEX","LEPR"),S4P=c("FOXJ1"),epsilon=c("GHRL","NPY1R","LINC00852","VTN"),gamma=c("PPY","GPC5-AS1","THSD7A"),co=c("IAPP","ADCYAP1","NPTX2"), ductal=c("KRT19","SPP1","HNF1B"),acinar=c("PRSS1","PNLIP"),endothelial=c("GNG11","SPARCL","SNHPG7","SPARCL1"),PSC=c("CALD1","VIM","LGALS1"), endocrine=c("CHGB","SCG5","CPE"), CC=c("TOP2A","CDK1"),CCI=c("CDKN1C"),alpha_prec=c("INS","GCG","ARX"),SCEC_p=c("CHGA","TPH1","PCSK1","NKX6.1","TPH1","DDC","SLC18A1","LMX1A","ADRΑ2A","FEV","TAC1","CXCL14"), SCEC_N=c("G6PC2","NPTX2","ISL1","PDX1")) #pancreatic stellate cells
human.mk.list.short <- human.mk.list %>% lapply(function(x) {head(x,3)})

main.mk <- list(beta="INS",alpha="GCG",endo=c("CHCG","TPH1"),late_beta=c("IAPP","HOPX"),early_beta=c("ASCL1","HADH"),late_beta=c("IGF2"),prolif=c("MKI67","CDK1","TOP2A"),delta="SST",psc=c("LGALS1","VIM"))

prolif.mk <- c("MKI67","CDK1","TOP2A","CCNB2","CCNA2","PBK")
#c(beta="INS",alpha="GCG",early_beta="HADH",mature_beta="BACE2",beta="PCDH7",polyhormonal="ASCL1")

beta.stage.mk <- list(up=c("IAPP","HOPX","NEFM","SIX2"),down=c("IGF2","LDHA"))
bridge.mk <- c("DPYSL3","GAP43","MAPT","NEUROG3","RIMS3","STMN4","LICAM")


FunDEGCal <- function(temp.compair,temp.cell.list,data.deg.temp) {
  
  s1=unlist(strsplit(temp.compair,split="_vs_"))[1]
  s2=unlist(strsplit(temp.compair,split="_vs_"))[2]
  c1=cell.list[[s1]]
  c2=cell.list[[s2]]
  data.deg.temp= subset(data.deg.temp,cell=c(c1,c2))
  Idents(data.deg.temp)=as.factor(c(rep("g1",length(c1)),rep("g2",length(c2))) %>% setNames(c(c1,c2)))[rownames(data.deg.temp@meta.data)]
  
  G1G2.DEG <- suppressMessages(FindMarkers(data.deg.temp,ident.1="g1",ident.2="g2",assay="RNA",verbose=F,logfc.threshold=0.1) %>% tibble::rownames_to_column(var="gene")) %>% tbl_df()
  G1.sig.up <- G1G2.DEG %>% filter(avg_log2FC >0.25  &  p_val_adj <0.05)
  G2.sig.up <- G1G2.DEG %>% filter(avg_log2FC < -0.25 & p_val_adj <0.05)
  
  
  G1G2.DEG.sig <- G1.sig.up %>% bind_rows(G2.sig.up)
  
  temp.up.ID <- G1.sig.up$gene
  temp.down.ID<- G2.sig.up$gene
  
  temp.out <- list()
  
  temp.out[["DEG.all.result"]] <- G1G2.DEG
  temp.out[["DEG.result"]] <- G1G2.DEG.sig
  temp.out[["DEG.result.up"]] <-  G1.sig.up
  temp.out[["DEG.result.down"]] <-  G2.sig.up
  if (nrow(G1.sig.up) > 0) {
    #temp.out[["DEG.up.GO.result"]] <- topGO_enrichment2(temp.up.ID,ALL_gene,0.05,geneID2GO)
    temp.out[["DEG.up.GO.result"]] <- NA
  }else{
    DEG.results[[temp.compair]] [["DEG.up.GO.result"]] <- NA
  }
  if (nrow(G2.sig.up) > 0) {
    #temp.out[["DEG.down.GO.result"]] <- topGO_enrichment2(temp.down.ID,ALL_gene,0.05,geneID2GO)
    temp.out[["DEG.down.GO.result"]] <- NA
  }else{
    temp.out[["DEG.down.GO.result"]] <- NA
  }
  return(temp.out)
}
