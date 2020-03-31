
resLocisToGenes<-function(resLocis){
  genes<-na.omit(unique(resLocis$gene))
  df<-data.frame(row.names =genes,nLocis=table(resLocis$gene)[genes])
  colnames(df)<-c("gene","nLocis")
  return(df)
}


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
    return(df[,-c(1,2)])
  }else{
    return(data.frame(row.names = genesInCommon,nLocis=resGenes[genesInCommon,"nLocis"],resGenesExpr[genesInCommon,]))
  }
                           
  
  
}



addMarqueursClusters<-function(scores,resGenes){
  genesMarqueurs<-intersect(rownames(scores),rownames(resGenes))
  print(paste(length(genesMarqueurs),"marqueurs de cluster"))
  resGenes[genesMarqueurs,"marqueur_de"]<-scores[genesMarqueurs,"dans_clusters_"]
  clusters<-colnames(scores)[str_detect(colnames(scores),"^cluster")]
  
  for(gene in genesMarqueurs){
    resGenes[gene,"scores"]<-paste(scores[gene,clusters][scores[gene,clusters]>0],collapse = "/")
  }
  return(resGenes)
}


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


#locis : 

mergeCols<-function(df,mergeColsName,colsToMerge=NULL,top=4,filter=0,abs=TRUE,roundNum=1,conserveCols=FALSE){
  if(is.null(colsToMerge)){
    colsToMerge<-colnames(df)
  }
  #use top=4 pour garder que les 4 plus grosses valeurs des cols tomerge (numeriques)
  if(conserveCols){
    newDf<-data.frame(row.names = rownames(df),df[,])
  }else{
    newDf<-data.frame(row.names = rownames(df),df[,!(colnames(df)%in%colsToMerge)])
  }
  
  if(is.numeric(top)){
    for(line in rownames(df)){
      if(abs){
        colsToMergeF<-colsToMerge[which(sapply(df[line,colsToMerge],abs)>filter)]
        if(length(colsToMergeF)>0){
          colsToMergeOrdered<-colsToMergeF[order(sapply(df[line,colsToMergeF],abs),decreasing = T)]
        }else{
          colsToMergeOrdered<-NA
        }
        
        
      }else{
        colsToMergeF<-na.omit(colsToMerge[df[line,colsToMerge]>filter])
        if(length(colsToMergeF)>0){
          colsToMergeOrdered<-colsToMergeF[order(df[line,colsToMergeF],decreasing = T)]
        }else{
          colsToMergeOrdered<-NA
        }
        
        
      }
      
      

      if(length(colsToMergeOrdered)>abs(top)){
        if(top>0){
          colsToMergeOrdered<-colsToMergeOrdered[1:top]
        }else if(top<0){
          colsToMergeOrdered<-colsToMergeOrdered[(length(colsToMergeOrdered)+(top-1)):length(colsToMergeOrdered)]
        }
        
        
      }
      subColNames<-paste0("top",top,mergeColsName)
      newDf[line,subColNames]<-paste(colsToMergeOrdered,collapse = "/")
      if(all(is.na(colsToMergeOrdered))){
        valeurs<-NA
      }else{
        valeurs<-paste(round(df[line,colsToMergeOrdered],roundNum),collapse = "/")
      }
      newDf[line,paste0(subColNames,"Value")]<-valeurs
    }
    return(newDf)
    
  }else{
    
    newDf[,paste0(mergeColsName,"Cols")]<-paste(colsToMerge,collapse = "/")
    newDf[,mergeColsName]<-apply(df[,colsToMerge],1,paste,collapse = "/")
    return(newDf)
    
  }
}


#annot locis 
annotLocis<-function(resLocis,resGenes,annots=NULL,genes=NULL,filter=T){
  if(is.null(annots)){
    annots<-colnames(resGenes)
  }
  if(is.null(genes)){
    genes<-rownames(resGenes)
  }
  poss<-c()
  for(gene in genes){
    pos<-which(resLocis$gene==gene)
    poss<-union(poss,pos)
    for(annot in annots){
      resLocis[pos,annot]<-resGenes[gene,annot]
      
      
    }
   
  }
  if(filter){
    return(resLocis[poss,])
  }else{
    return(resLocis)
  }
  

}



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

