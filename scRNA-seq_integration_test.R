options(stringsAsFactors = F)
source("scripts/scoreCluster.R")
library(Seurat)
source("scripts/FindInfosGenes.R")
outputDir <- file.path("analyses","test_scRNAseq_integration")
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")


#jeu de res
res<-read.csv2("analyses/main/top1000Locis_Group.Sex__model1BasicsLocisFilteringAndSamples_F_annotated_280320.csv",row.names = 1)


#scRNA-seq Data
samples<-readRDS("../../../Alexandre_SC/analyses/test_batch_integration/CTRLandLGA_SCTransform_Integrated.rds")
samples
# 
# head(samples[["SCT"]])[,1:10]
# tail(samples[["SCT"]])[,1:10]
# samples[["SCT"]] #9k cells donc LGA+CTRL
# head(samples@meta.data)
# tail(samples@meta.data)
#annoter cluster
#unique(samples@meta.data$seurat_clusters)
# # 0, 1,2,3 : HSC; 4 : EMP, 5 : GMP, 6 : LyP, 7 : MkP, 8 :pro T, 9 : LMPP, 10 : pre B, 11 : LT-HSC, 12 : Neu, 13 : Ly B, 14 Ly-ETS1, 15 : DC
new.cluster.ids <- c("HSC-SELL2", "HSC-AVP", "HSC-Er", "HSC-SELL1", "EMP", "GMP",
                     "LyP", "MkP", "proT","LMPP","preB","LT-HSC","Neu","LyB","Ly-ETS1","DC")
names(new.cluster.ids) <- levels(samples)
samples <- RenameIdents(samples, new.cluster.ids)
# DimPlot(samples, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# 
# #genes expr in CTRL
# CTRL<-subset(samples, orig.ident=="freshCD34_CBP547_CTRL")
# head(CTRL@meta.data)
# DimPlot(CTRL, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# levels(CTRL)
# CTRL_expr_by_pop <- data.frame(row.names = rownames(CTRL[["SCT"]]@data))
# 
# for (pop in levels(CTRL) ){
#   CTRL_expr_by_pop[,pop] <- as.vector(rowMeans(x = as.matrix(CTRL[["SCT"]]@data[,Idents(CTRL)==pop])))
# }
# #bulk expr in CTRL : 
# CTRL_expr_by_pop[,"bulk"]<-rowMeans(as.matrix(CTRL[["SCT"]]@data))
# head(CTRL_expr_by_pop)
# 
# #proportion of cells in each cluster : 
# for (pop in levels(CTRL) ){
#   CTRL_expr_by_pop["proportion",pop]<-sum(Idents(CTRL)==pop)/ncol(CTRL)
# }
# CTRL_expr_by_pop["proportion",]
# sum(CTRL_expr_by_pop["proportion",],na.rm = T)
# 
# #ya til des des genes ds matrix  pas du tout expr ?
# sum(rowSums(CTRL_expr_by_pop,na.rm = T)==0) #75, on les flag
# CTRL_expr_by_pop[,"expr"]<-F
# CTRL_expr_by_pop[rowSums(CTRL_expr_by_pop,na.rm = T)!=0,"expr"]<-T
# sum(CTRL_expr_by_pop$expr) #7518
# 
# #genes expr in LGA
# LGA<-subset(samples, orig.ident=="freshCD34_CBP552_LGA")
# head(LGA@meta.data)
# DimPlot(LGA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# levels(LGA)
# LGA_expr_by_pop <- data.frame(row.names = rownames(LGA[["SCT"]]@data))
# 
# for (pop in levels(LGA) ){
#   LGA_expr_by_pop[,pop] <- as.vector(rowMeans(x = as.matrix(LGA[["SCT"]]@data[,Idents(LGA)==pop])))
# }
# #bulk expr in LGA : 
# LGA_expr_by_pop[,"bulk"]<-rowMeans(as.matrix(LGA[["SCT"]]@data))
# head(LGA_expr_by_pop)
# 
# #proportion of cells in each cluster : 
# for (pop in levels(LGA) ){
#   LGA_expr_by_pop["proportion",pop]<-sum(Idents(LGA)==pop)/ncol(LGA)
# }
# LGA_expr_by_pop["proportion",]
# sum(LGA_expr_by_pop["proportion",],na.rm = T)
# 
# #ya til des des genes ds matrix  pas du tout expr ?
# sum(rowSums(LGA_expr_by_pop,na.rm = T)==0) #1886, on les flag
# LGA_expr_by_pop[,"expr"]<-F
# LGA_expr_by_pop[rowSums(LGA_expr_by_pop,na.rm = T)!=0,"expr"]<-T
# sum(LGA_expr_by_pop$expr) #5749
# 
# write.csv2(CTRL_expr_by_pop,file.path(outputDir,"geneExprInCTRL.csv"),row.names = T)
# write.csv2(LGA_expr_by_pop,file.path(outputDir,"geneExprInLGA.csv"),row.names = T)

#ajout au res
#genes res
#count locis par genes

resLocisToGenes<-function(resLocis){
  genes<-na.omit(unique(res$gene))
  df<-data.frame(row.names =genes,nLocis=table(res$gene)[genes])
  colnames(df)<-c("gene","nLocis")
  return(df)
}
resGenes<-resLocisToGenes(res)
head(resGenes)


#genes Expr
CTRL_expr_by_pop<-read.csv2(file.path(outputDir,"geneExprInCTRL.csv"),row.names = 1)
LGA_expr_by_pop<-read.csv2(file.path(outputDir,"geneExprInLGA.csv"),row.names = 1)

listExprMat<-list(CTRL=CTRL_expr_by_pop,LGA=LGA_expr_by_pop)

findGenesExpr<-function(resGenes,listExprMat,retourne="all"){
  #retourne df avec expr des genes dans chaque pop (colonnes des mat)
  allGenes<-unique(as.vector(sapply(listExprMat, rownames)))
  genesInCommon<-intersect(rownames(resGenes),allGenes)
  
  print(paste(length(genesInCommon),"genes en commun")) 
  popNames<-paste(colnames(listExprMat[[1]]),names(listExprMat)[1],sep = "_")
  
  resGenesExpr<-data.frame(row.names = genesInCommon,
                           
                           listExprMat[[1]][genesInCommon,])
  cols<-c(popNames)
  colnames(resGenesExpr)<- cols
  
  for(i in 2:length(listExprMat)){
    resGenesExpr<-data.frame(row.names = genesInCommon,
                             resGenesExpr[genesInCommon,],
                             listExprMat[[i]][genesInCommon,])
    popNames<-paste(colnames(listExprMat[[i]]),names(listExprMat)[i],sep = "_")
    cols<-c(cols,popNames)
    colnames(resGenesExpr)<- cols

    
  }
  #mettre les cols d'une meme pop a coté
  
  cols<-c() #colone "nLocis" a conserver au debut
  decalcols<-1:(length(listExprMat)-1)
  ncols.ech<-ncol(listExprMat[[1]]) #nb de colonne par ech pour faire le decallage
  
  for(i in 1:ncols.ech){
    pos1st<-i
    newCols<-c(pos1st,pos1st+ncols.ech*decalcols)
    cols<-c(cols,newCols)
    
  }
  resGenesExpr<-resGenesExpr[,cols]
  
  if(retourne=="all"){
    df<-merge.data.frame(resGenes,resGenesExpr,by="row.names",all.x = T,all.y = T)
    rownames(df)<-df$gene
    return(df)
  }else{
    return(data.frame(row.names = genesInCommon,nLocis=resGenes[genesInCommon,"nLocis"],resGenesExpr[genesInCommon,]))
  }
                           
  
  
}

resGenes<-findGenesExpr(resGenes,listExprMat = listExprMat,retourne = "expr")
head(resGenes)


#markers of clusters
# 0, 1,2,3 : HSC; 4 : EMP, 5 : GMP, 6 : LyP, 7 : MkP, 8 :pro T, 9 : LMPP, 10 : pre B, 11 : LT-HSC, 12 : Neu, 13 : Ly B, 14 Ly-ETS1, 15 : DC

markers<-read.csv2("../../../Alexandre_SC/analyses/test_batch_integration/all.markers_SCTransform_CTRLandLGA_integrated.csv",row.names = 1)
head(markers)
for(i in 1:length(new.cluster.ids)){
  markers$cluster[markers$cluster==as.numeric(names(new.cluster.ids)[i])]<-new.cluster.ids[i]
}

scoresCluster<-scoreMarquageCluster(markers,samples,seuil = "intraClusterFixe",filtreMin = 2)
head(scoresCluster)

addMarqueursClusters<-function(scores,resGenes){
  genesMarqueurs<-intersect(rownames(scoresCluster),rownames(resGenes))
  print(paste(length(genesMarqueurs),"marqueurs de cluster"))
  resGenes[genesMarqueurs,"marqueur_de"]<-scoresCluster[genesMarqueurs,"dans_clusters_"]
  clusters<-colnames(scoresCluster)[str_detect(colnames(scoresCluster),"^cluster")]
  
  for(gene in genesMarqueurs){
    resGenes[gene,"scores"]<-paste(scores[gene,clusters][scores[gene,clusters]>0],collapse = "/")
  }
  return(resGenes)
}

resGenes<-addMarqueursClusters(scoresCluster,resGenes)
head(resGenes,100)

#markers known

findMarqueursPops<-function(genes,df=NULL){
  source("scripts/ManagePopMarqueurs.R")
  if(is.null(df)){
    df<-data.frame()
  }
  
  for (gene in genes){
    pops<-GetPopDe(gene)
    if(!is.null(pops)){
      print(paste(gene,"marqueurs de ",paste(pops,collapse = ", ")))
      df[gene,"MarqueurPops"]<-paste(pops,collapse = ";")
    }else{
      df[gene,"MarqueurPops"]<-NA
    }
    
  }
  print(paste(sum(!is.na(df$MarqueurPops)),"genes sont des marqueurs connus de sous-population"))
  return(df)
}

resGenes<-findMarqueursPops(rownames(resGenes),df = resGenes)
head(resGenes)
sum(!is.na(resGenes$MarqueurPops)) #25

#TF et fct
genesToAnnot<-rownames(resGenes)[(resGenes$expr_CTRL|resGenes$expr_LGA)&(!is.na(resGenes$marqueur_de)|!is.na(resGenes$MarqueurPops))]
length(genesToAnnot) #641

resGenes<-find_fonction(genesToAnnot,df = resGenes,save=F)
head(resGenes)



resGenes<-findIfTF(rownames(resGenes),resGenes,save=F,large = T)


head(resGenes)
resGenes<-findCbIn(rownames(resGenes),
                     listeKeywords = list(HSPC="hematopo|myeloid|lymphoid|HSC|HSPC",
                                          lineage="lineage decision|differentiation|cell fate",
                                          stress="stress",
                                          signaling="kinase|signaling|pathway"),
                     genes_infos=resGenes)
head(resGenes)
head(resGenes)
write.csv2(resGenes,file.path(outputDir,"resGenesIntegration.csv"),row.names = T)

#locis : 

for(gene in genesInCommon){
  pos<-which(res$gene==gene)
  #ajout expr in bulk
  res[pos,"bulk_CTRL"]<-genesResExpr[gene,"bulk_CTRL"]
  res[pos,"bulk_LGA"]<-genesResExpr[gene,"bulk_LGA"]
  
  #ajout expr in pop
  pops<-colnames(genesResExpr)[2:(which((colnames(genesResExpr)=="exprIn_CTRL"))-1)]
  ExprIn<-pops[genesResExpr[gene,pops]>0]
  ExprInOrdered<-ExprIn[order(genesResExpr[gene,ExprIn],decreasing = T)]
  if(length(ExprInOrdered)>4){
    ExprInOrdered<-ExprInOrdered[1:4]
  }
  res[pos,"ExprIn"]<-paste(ExprInOrdered,collapse = "/")
  res[pos,"log2Expr"]<-paste(round(genesResExpr[gene,ExprInOrdered],1),collapse = "/")
  
  #ajout top DEmarkers
  
  res[pos,"marqueur_de"]<-genesResExpr[gene,"marqueur_de"]
  res[pos,"scores"]<-genesResExpr[gene,"scores"]
  res[pos,"marqueur_connu_de"]<-genesResExpr[gene,"MarqueurPops"]
  res[pos,"TF"]<-genesResExpr[gene,"TF"]
  res[pos,"HSPC"]<-genesResExpr[gene,"HSPC"]
  res[pos,"lineage"]<-genesResExpr[gene,"lineage"]
  res[pos,"stress"]<-genesResExpr[gene,"stress"]
  res[pos,"fonction"]<-genesResExpr[gene,"Fonction"]
}
head(res[!(is.na(res$ExprIn)),],100)
write.csv2(res,file=file.path(outputDir,"res_with_scRNA-seq_integration.csv"),row.names = T)


# #! a faire : 
# #pour scorer les markers selon leur proximité phenotypique avec les autre markers, need bioprocess infos
# genes<-unique(DEmarkers$gene)
# genesInfos<-find_description(genes)
# genesInfos<-find_GoBioProcess(rownames(scoresCluster),df = genesInfos,save = T)
# #que l'on nettoir (pour eviter redondance de terme)
# genesInfos<-cleanBioProcess(genes, genesInfos)
# 
# clusters<-names(scoresCluster)[str_detect(names(scoresCluster),"cluster[0-9]")]
# linksScore<-data.frame(row.names =rownames(scoresCluster) )
# 
# 
# for (cluster in clusters){
#   numCluster<-as.numeric(strsplit(cluster,"ster")[[1]][2])
#   
#   genes<-rownames(scoresCluster)[scoresCluster[,cluster]!=0]
#   
#   for (gene in genes){
#     DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster, "score"]<-scoresCluster[gene,cluster]
#     DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster, "dans_clusters"]<-scoresCluster[gene,"dans_clusters_"]
#     #trouve genes links between genes
#     linkedGenes<-TrouveGenesLinks(gene,genes,genesInfos)
#     DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster,"LinkedGenes"]<-paste(linkedGenes,collapse = ";")
#     if(length(linkedGenes)<5){
#       score<-length(linkedGenes)
#     }else{
#       score<-5
#     }
#     
#     DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster,"PhenotLinksScore"]<-score
#     linksScore[gene,cluster]<-score
#   }
#   print(numCluster)
#   
# }
# 
# DEmarkers<-subset(DEmarkers,!is.na(score))
# 
# DEmarkers<-subset(DEmarkers,score>=scoreMin)
# head(DEmarkers)
# linksScore<-as.matrix(linksScore)
# linksScore[is.na(linksScore)]<-0
# 
# ##annotation genes
# 
# #genesinfos<-annotation_auto_genes(genes = unique(DEmarkers$gene))
# #ùettre Fonction, description, TF, tissue match
# genes<-unique(DEmarkers$gene)
# genesInfos<-subset(genesInfos,rownames(genesInfos)%in%genes)
# 
# 
# genesInfos<-find_fonction(genes,genesInfos,save=T)
# genesInfos<-findIfTF(genes,genesInfos,save=T,large = T)
# 
# 
# genesInfos<-findCbIn(genes,
#                      listeKeywords = list(neuron="neur(o|a)",astrocyte="astro(c|g)|macrogli",oligodendrocyte="oligodendrocyt",progenitor="progen|precur",gliale="glial"),genes_infos=genesInfos)
# 
# 
# View(genesInfos)
# 
# for(gene in rownames(genesInfos)){
#   for (annotation in annotationsDesGenes){
#     
#     if(gene %in% DEmarkers$gene){
#       DEmarkers[DEmarkers$gene==gene,annotation]<-genesInfos[gene,annotation]
#     }
#     
#     
#   }
#   
# }
# View(DEmarkers)
# #top10 par cluster selon leur scoreCLuster+scorePhenotLinks
# topN<-DEmarkers %>% group_by(cluster) %>%arrange(desc(score)) %>% top_n(n = NbTopMarkerDansHeatmap, wt = score + PhenotLinksScore)
# topN<-topN %>% group_by(cluster) %>% top_n(n = NbTopMarkerDansHeatmap, wt = avg_logFC)%>%arrange(cluster)
# 
# DEmarkersTop<-data.frame(topN)
# View(DEmarkersTop)
# 
# #enlever col pas informative
# names(DEmarkersTop)
# colToRm<-c("astrocyte","p_val")
# cols<-setdiff(names(DEmarkersTop),colToRm)
# cols<-cols[c(6,1:5,7,10,9,11:17,8)]
# colsSansAnnot<-cols[1:7]
# #sauvegarder
# write.csv2(DEmarkersTop[cols], file =file.path(outputDir,paste(sample_desc,"TopMarkersAvecAnnotation_resol",resol,".csv",sep="_")))
# write.csv2(DEmarkersTop[colsSansAnnot], file =file.path(outputDir,paste(sample_desc,"TopMarkers_resol",resol,".csv",sep="_")))

