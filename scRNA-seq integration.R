
options(stringsAsFactors = F)
source("scripts/scoreCluster.R")

source("scripts/FindInfosGenes.R")
outputDir <- file.path("analyses","test_scRNAseq_integration")
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")


#jeu de res
res<-read.csv2("analyses/main/top1000Locis_Group.Sex__model1BasicsLocisFilteringAndSamples_F_annotated_280320.csv",row.names = 1)
#genes
genes<-na.omit(unique(res$gene))
head(genes)

#scRNA-seq Data
samples<-readRDS("../../../Alexandre_SC/analyses/test_batch_integration/CTRLandLGA_SCTransform_Integrated.rds")
samples

head(samples[["SCT"]])[,1:10]
tail(samples[["SCT"]])[,1:10]
samples[["SCT"]] #9k cells donc LGA+CTRL
head(samples@meta.data)
tail(samples@meta.data)
#annoter cluster 
unique(samples@meta.data$seurat_clusters)
# 0, 1,2,3 : HSC; 4 : EMP, 5 : GMP, 6 : LyP, 7 : MkP, 8 :pro T, 9 : LMPP, 10 : pre B, 11 : LT-HSC, 12 : Neu, 13 : Ly B, 14 Ly-ETS1, 15 : DC
new.cluster.ids <- c("HSC-SELL2", "HSC-AVP", "HSC-Er", "HSC-SELL1", "EMP", "GMP", 
                     "LyP", "MkP", "proT","LMPP","preB","LT-HSC","Neu","LyB","Ly-ETS1","DC")
names(new.cluster.ids) <- levels(samples)
samples <- RenameIdents(samples, new.cluster.ids)
DimPlot(samples, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#genes expr in CTRL
CTRL<-subset(samples, orig.ident=="freshCD34_CBP547_CTRL")
head(CTRL@meta.data)
DimPlot(CTRL, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
levels(CTRL)
CTRL_expr_by_pop <- data.frame(row.names = rownames(CTRL[["SCT"]]@data))

for (pop in levels(CTRL) ){
  CTRL_expr_by_pop[,pop] <- as.vector(rowMeans(x = as.matrix(CTRL[["SCT"]]@data[,Idents(CTRL)==pop])))
}
#bulk expr in CTRL : 
CTRL_expr_by_pop[,"bulk"]<-rowMeans(as.matrix(CTRL[["SCT"]]@data))
head(CTRL_expr_by_pop)

#proportion of cells in each cluster : 
for (pop in levels(CTRL) ){
  CTRL_expr_by_pop["proportion",pop]<-sum(Idents(CTRL)==pop)/ncol(CTRL)
}
CTRL_expr_by_pop["proportion",]
sum(CTRL_expr_by_pop["proportion",],na.rm = T)


#genes expr in LGA
LGA<-subset(samples, orig.ident=="freshCD34_CBP552_LGA")
head(LGA@meta.data)
DimPlot(LGA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
levels(LGA)
LGA_expr_by_pop <- data.frame(row.names = rownames(LGA[["SCT"]]@data))

for (pop in levels(LGA) ){
  LGA_expr_by_pop[,pop] <- as.vector(rowMeans(x = as.matrix(LGA[["SCT"]]@data[,Idents(LGA)==pop])))
}
#bulk expr in LGA : 
LGA_expr_by_pop[,"bulk"]<-rowMeans(as.matrix(LGA[["SCT"]]@data))
head(LGA_expr_by_pop)

#proportion of cells in each cluster : 
for (pop in levels(LGA) ){
  LGA_expr_by_pop["proportion",pop]<-sum(Idents(LGA)==pop)/ncol(LGA)
}
LGA_expr_by_pop["proportion",]
sum(LGA_expr_by_pop["proportion",],na.rm = T)


#ajout au res
#genes
genesInCommon<-genes[genes%in%union(rownames(CTRL_expr_by_pop),rownames(LGA_expr_by_pop))]
genesRes<-data.frame(row.names = genesInCommon, 
                    CTRL_expr_by_pop[genesInCommon,], 
                    LGA_expr_by_pop[genesInCommon,]
                    
                    )
colnames(genesRes)<- c("HSC.SELL2_CTRL","HSC.AVP_CTRL","HSC.Er_CTRL","HSC.SELL1_CTRL","EMP_CTRL",
                       "GMP_CTRL","LyP_CTRL", "MkP_CTRL", "proT_CTRL", "LMPP_CTRL" ,"preB_CTRL","LT.HSC_CTRL",
                       "Neu_CTRL","LyB_CTRL","Ly.ETS1_CTRL","DC_CTRL","bulk_CTRL",
                       "HSC.SELL2_LGA","HSC.AVP_LGA","HSC.Er_LGA","HSC.SELL1_LGA","EMP_LGA",
                       "GMP_LGA","LyP_LGA", "MkP_LGA", "proT_LGA", "LMPP_LGA" ,"preB_LGA","LT.HSC_LGA",
                       "Neu_LGA","LyB_LGA","Ly.ETS1_LGA","DC_LGA","bulk_LGA")
ncol(genesRes)
cols<-c()
for(i in 1:(ncol(genesRes)/2)){
  cols<-c(cols,i,i+17)
}


genesRes<-genesRes[,cols]
head(genesRes)


#locis

for(gene in genesInCommon){
  pos<-which(res$gene==gene)
  res[pos,"bulk_CTRL"]<-genesRes[gene,"bulk_CTRL"]
  res[pos,"bulk_LGA"]<-genesRes[gene,"bulk_LGA"]
  
  ExprIn<-colnames(genesRes)[genesRes[gene,]>0]
  ExprInOrdered<-ExprIn[order(genesRes[gene,ExprIn],decreasing = T)]
  if(length(ExprInOrdered)>4){
    ExprInOrdered<-ExprInOrdered[1:4]
  }
  
  res[pos,"ExprIn"]<-paste(ExprInOrdered,collapse = "/")
  res[pos,"log2Expr"]<-paste(round(genesRes[gene,ExprInOrdered],1),collapse = "/")
}
head(res[!(is.na(res$ExprIn)),],100)
write.csv2(res,file=file.path(outputDir,"res_with_scRNA-seq_integration.csv"),row.names = T)


#markers of clusters
markers<-read.csv2("../../../Alexandre_SC/analyses/test_batch_integration/all.markers_SCTransform_CTRLandLGA_integrated.csv")
head(markers)
for(i in 1:length(new.cluster.ids)){
  markers$cluster[markers$cluster==as.numeric(names(new.cluster.ids)[i])]<-new.cluster.ids[i]
}

  
scoresCluster<-scoreMarquageCluster(markers,samples,seuil = "intraClusterFixe",filtreMin = 2)
#pour scorer les markers selon leur proximité phenotypique avec les autre markers, need bioprocess infos
genes<-unique(DEmarkers$gene)
genesInfos<-find_description(genes)
genesInfos<-find_GoBioProcess(rownames(scoresCluster),df = genesInfos,save = T)
#que l'on nettoir (pour eviter redondance de terme)
genesInfos<-cleanBioProcess(genes, genesInfos)

clusters<-names(scoresCluster)[str_detect(names(scoresCluster),"cluster[0-9]")]
linksScore<-data.frame(row.names =rownames(scoresCluster) )


for (cluster in clusters){
  numCluster<-as.numeric(strsplit(cluster,"ster")[[1]][2])
  
  genes<-rownames(scoresCluster)[scoresCluster[,cluster]!=0]
  
  for (gene in genes){
    DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster, "score"]<-scoresCluster[gene,cluster]
    DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster, "dans_clusters"]<-scoresCluster[gene,"dans_clusters_"]
    #trouve genes links between genes
    linkedGenes<-TrouveGenesLinks(gene,genes,genesInfos)
    DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster,"LinkedGenes"]<-paste(linkedGenes,collapse = ";")
    if(length(linkedGenes)<5){
      score<-length(linkedGenes)
    }else{
      score<-5
    }
    
    DEmarkers[DEmarkers$gene == gene & DEmarkers$cluster== numCluster,"PhenotLinksScore"]<-score
    linksScore[gene,cluster]<-score
  }
  print(numCluster)
  
}

DEmarkers<-subset(DEmarkers,!is.na(score))

DEmarkers<-subset(DEmarkers,score>=scoreMin)
head(DEmarkers)
linksScore<-as.matrix(linksScore)
linksScore[is.na(linksScore)]<-0

##annotation genes

#genesinfos<-annotation_auto_genes(genes = unique(DEmarkers$gene))
#ùettre Fonction, description, TF, tissue match
genes<-unique(DEmarkers$gene)
genesInfos<-subset(genesInfos,rownames(genesInfos)%in%genes)


genesInfos<-find_fonction(genes,genesInfos,save=T)
genesInfos<-findIfTF(genes,genesInfos,save=T,large = T)


genesInfos<-findCbIn(genes,
                     listeKeywords = list(neuron="neur(o|a)",astrocyte="astro(c|g)|macrogli",oligodendrocyte="oligodendrocyt",progenitor="progen|precur",gliale="glial"),genes_infos=genesInfos)


View(genesInfos)

for(gene in rownames(genesInfos)){
  for (annotation in annotationsDesGenes){
    
    if(gene %in% DEmarkers$gene){
      DEmarkers[DEmarkers$gene==gene,annotation]<-genesInfos[gene,annotation]
    }
    
    
  }
  
}
View(DEmarkers)
#top10 par cluster selon leur scoreCLuster+scorePhenotLinks
topN<-DEmarkers %>% group_by(cluster) %>%arrange(desc(score)) %>% top_n(n = NbTopMarkerDansHeatmap, wt = score + PhenotLinksScore)
topN<-topN %>% group_by(cluster) %>% top_n(n = NbTopMarkerDansHeatmap, wt = avg_logFC)%>%arrange(cluster)

DEmarkersTop<-data.frame(topN)
View(DEmarkersTop)

#enlever col pas informative
names(DEmarkersTop)
colToRm<-c("astrocyte","p_val")
cols<-setdiff(names(DEmarkersTop),colToRm)
cols<-cols[c(6,1:5,7,10,9,11:17,8)]
colsSansAnnot<-cols[1:7]
#sauvegarder
write.csv2(DEmarkersTop[cols], file =file.path(outputDir,paste(sample_desc,"TopMarkersAvecAnnotation_resol",resol,".csv",sep="_")))
write.csv2(DEmarkersTop[colsSansAnnot], file =file.path(outputDir,paste(sample_desc,"TopMarkers_resol",resol,".csv",sep="_")))

topScores<-scoresCluster[unique(topN$gene),clusters]
#rownames(topScores)<-paste(rownames(topScores),genesinfos[rownames(topScores),"MarqueurPops"],sep = ",")

ngenes<-nrow(topScores)
dfAnnot<-data.frame(row.names =rownames(topScores),
                    Expr_TopCluster=apply(norm.data[rownames(topScores),],1,max),
                    TF=as.numeric(genesInfos[rownames(topScores),"TF"]),
                    neuron_rel=genesInfos[rownames(topScores),"neuron"])

#juste score Cluster
pm<-pheatmap(topScores,fontsize = 9,
             
             annotation_row = dfAnnot
             
)