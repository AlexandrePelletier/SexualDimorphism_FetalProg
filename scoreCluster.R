# #mis en place de la fct score cluster qui fixe un seuil de pval en fct du nombre de cellules 
# #determiné avec une droite de régression nb de cellules ~ signif moyen des top10 marqueurs
# 
# sample<-readRDS("TestCycleCell/freshCTRL_SCTransform_norm_Mito_CycleCellDifference.rds")
# 
# #DEFINITION DU SEUIL DE PVAL
# #check si nb de cellule influence bien la pvaladj
# #en reduisant le nb de cellule
# #a partir CTRL
# plot(density(log10(DEmarkersWilco$p_val_adj)))
# for (i in 2:10){
# 
#   downsample<-sample[,sample(colnames(sample),round(ncol(sample)/i))]
#   markers<-FindAllMarkers(downsample, min.pct = 0.25,only.pos = T, logfc.threshold = 0.25)
#   lines(density(log10(markers$p_val_adj)), col = i)
#   
# }
# #oui, plus tu diminue en nb de cellule, plus tu diminue la significativité
# #il faut qu'il y ai  le meme nombre ratio de marker qui passe le seuil qlq soit le nb de cellules
# #pour cela, on prend le quantile qui correspond au seuil 10^50 definit avant
# abline(v = log10(quantile(DEmarkersWilco$p_val_adj,0.15))) #q15 semble etre un bon seuil pour scorer pval_adj
# #verifions que tout les markers soit signif avec ce seuil qlq soit le nb de cellule
# print(summary(DEmarkersWilco[DEmarkersWilco$p_val_adj<quantile(DEmarkersWilco$p_val_adj,0.15),]))
# for (i in 2:5){
#   
#   downsample<-sample[,sample(colnames(sample),round(ncol(sample)/i))]
#   markers<-FindAllMarkers(downsample, min.pct = 0.25,only.pos = T, logfc.threshold = 0.25)
#   print(summary(markers[markers$p_val_adj<quantile(markers$p_val_adj,0.15),]))
#   
# }
# # pour faire un seuil fixe, on fait un droit de regression entre nb de cellule (x), et significativité pval(y)
# x<-c(ncol(sample))
# y<-c(quantile(DEmarkersWilco$p_val_adj,0.15))
# plot(density(log10(DEmarkersWilco$p_val_adj)),xlim=c(-100,10))
# abline(v=log10(quantile(DEmarkersWilco$p_val_adj,0.15)),col=1)
# for (i in 1:6+0.5){
#   nbcell<-round(ncol(sample)/i)
#   x<-c(x,nbcell)
#   downsample<-sample[,sample(colnames(sample),nbcell)]
#   downsample <- RunPCA(object = downsample, features = VariableFeatures(object = downsample))
#   downsample <- FindNeighbors(object = downsample, dims = c(1:30)) 
#   downsample <- FindClusters(object = downsample, resolution = 0.6) 
#   
#   markers<-FindAllMarkers(downsample, min.pct = 0.25,only.pos = T, logfc.threshold = 0.25)
#   q<-quantile(markers$p_val_adj,0.15)
#   y<-c(y,q)
#   lines(density(log10(markers$p_val_adj)), col = i)
#   abline(v=log10(q),col=i)
#   
# }
# 
# plot(x,log10(y))
# scatter.smooth(x, log10(y), main="q15(pval_adj) ~ nbCell") # on voit une droite mais pas assez de point. ON refait avec plus de point
# q_lm<-lm(log10(y)~x)
# abline(q_lm)
# summary(q_lm)
# "Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.040e+01  2.363e+00  -4.401 0.000604 ***
#   x           -8.824e-03  9.595e-04  -9.196 2.61e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.788 on 14 degrees of freedom
# Multiple R-squared:  0.858,	Adjusted R-squared:  0.8478 
# F-statistic: 84.57 on 1 and 14 DF,  p-value: 2.61e-07"
# 
deterSeuilPval<-function(nbCell){
  a<-(-8.824e-03)
  b<-(-10.4)
  return(10^(a*nbCell + b))
}
# 
# #DEFINITION DU SEUIL Delta PCT de detection
# #seuil interCLuster (intraEch)
# dt.pct<-abs(DEmarkersWilco$pct.1-DEmarkersWilco$pct.2)
# plot(density(dt.pct))
# # on prend le quantile qui correspond au seuil 10^50 definit avant
# abline(v = quantile(dt.pct,0.9)) #q90 (0.42) semble etre un bon seuil pour scorer dt.pct
# 
# #check si pct0 influence bien la dt.pct
# #en reduisant le nb genes avec bcp de zeros : enlever genes avec 3000:6000 0
# mat<-as.matrix(sample[["RNA"]]@counts)
# Xpct0<-c(sum(mat==0)/length(mat))
# Ydt.pct<-c(quantile(dt.pct,0.9))
# for (i in 1:10){
#   maxNbZeros<-(6000/i)
#   keep<-rowSums(mat==0)<maxNbZeros
#   print(sum(keep))
#   downsample<-sample[keep,]
#   submat<-as.matrix(downsample[["RNA"]]@counts)
#   pct0<-sum(submat==0)/length(submat)
#   Xpct0<-c(Xpct0,pct0)
#   print(pct0)
#   downsample <- RunPCA(object = downsample, features = VariableFeatures(object = downsample))
#   downsample <- FindNeighbors(object = downsample, dims = c(1:30)) 
#   downsample <- FindClusters(object = downsample, resolution = 0.6) 
#   
#   
#   markers<-FindAllMarkers(downsample, min.pct = 0.25,only.pos = T, logfc.threshold = 0.25)
#   
#   dt.pct<-abs(markers$pct.1-markers$pct.2)
#   lines(density(dt.pct), col = i)
#   print(quantile(dt.pct,0.9))
#   Ydt.pct<-c(Ydt.pct,quantile(dt.pct,0.9))
#   abline(v = quantile(dt.pct,0.9),col=i)
# }
# 
# plot(Xpct0,Ydt.pct^2) #need more points entre 0.1 et 0.8 de pct0
# 
# for (i in 1:4){
#   maxNbZeros<-(6400-3000/i)
#   keep<-rowSums(mat==0)<maxNbZeros
#   print(sum(keep))
#   downsample<-sample[keep,]
#   submat<-as.matrix(downsample[["RNA"]]@counts)
#   pct0<-sum(submat==0)/length(submat)
#   Xpct0<-c(Xpct0,pct0)
#   print(pct0)
#   downsample <- RunPCA(object = downsample, features = VariableFeatures(object = downsample))
#   downsample <- FindNeighbors(object = downsample, dims = c(1:30)) 
#   downsample <- FindClusters(object = downsample, resolution = 0.6) 
#   
#   
#   markers<-FindAllMarkers(downsample, min.pct = 0.25,only.pos = T, logfc.threshold = 0.25)
#   
#   dt.pct<-abs(markers$pct.1-markers$pct.2)
#   lines(density(dt.pct), col = i)
#   print(quantile(dt.pct,0.9))
#   Ydt.pct<-c(Ydt.pct,quantile(dt.pct,0.9))
#   abline(v = quantile(dt.pct,0.9),col=i)
# }
# plot(sqrt(Xpct0),Ydt.pct)
# #oui, plus tu diminue en pct de zeros, plus tu diminue q90 du dt.pct
# #il faut qu'il y ai  le meme nombre ratio de marker qui passe le seuil qlq soit le pct0 (complexite/profondeur de la librairie)
# 
# 
# scatter.smooth(sqrt(Xpct0),Ydt.pct, main="q10(dt.pct) ~ pct0") # on voit une droite mais pas assez de point. ON refait avec plus de point
# q_lm<-lm(Ydt.pct~sqrt(Xpct0))
# abline(q_lm)
# summary(q_lm)
# "Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.03232    0.01481   2.182    0.048 *  
# sqrt(Xpct0)  0.44754    0.02773  16.138 5.58e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03166 on 13 degrees of freedom
# Multiple R-squared:  0.9525,	Adjusted R-squared:  0.9488 
# F-statistic: 260.4 on 1 and 13 DF,  p-value: 5.581e-10"
# 
deterSeuilDtPct<-function(pct0){
  a<-0.44754
  b<-0.03232
  return(a*sqrt(pct0)+ b)
}
# 
# #DEFINITION DU SEUIL logFC 
# #seuil interCLuster (intraEch)
# 
# plot(density(DEmarkersWilco$avg_logFC))
# 
# # on prend le quantile qui prend que les meilleurs FC >0.6
# abline(v = quantile(DEmarkersWilco$avg_logFC,0.85)) #q80 (0.68) semble etre un bon seuil pour scorer logFC
# 
# #check si pct0 influence le logFC
# #en reduisant le nb genes avec bcp de zeros : enlever genes avec 3000:6000 0
# mat<-as.matrix(sample[["RNA"]]@counts)
# Xpct0<-c(sum(mat==0)/length(mat))
# YlogFC<-c(quantile(DEmarkersWilco$avg_logFC,0.8))
# for (i in 1:10){
#   maxNbZeros<-(6000/i)
#   keep<-rowSums(mat==0)<maxNbZeros
#   print(sum(keep))
#   downsample<-sample[keep,]
#   submat<-as.matrix(downsample[["RNA"]]@counts)
#   pct0<-sum(submat==0)/length(submat)
#   Xpct0<-c(Xpct0,pct0)
#   print(pct0)
#   downsample <- RunPCA(object = downsample, features = VariableFeatures(object = downsample))
#   downsample <- FindNeighbors(object = downsample, dims = c(1:30)) 
#   downsample <- FindClusters(object = downsample, resolution = 0.6) 
#   
#   
#   markers<-FindAllMarkers(downsample, min.pct = 0.25,only.pos = T, logfc.threshold = 0.25)
#   
#   
#   lines(density(markers$avg_logFC), col = i)
#   print(quantile(markers$avg_logFC,0.8))
#   YlogFC<-c(YlogFC,quantile(markers$avg_logFC,0.8))
#   abline(v = quantile(markers$avg_logFC,0.8),col=i)
# }
# 
# plot(Xpct0,YlogFC) #no influence donc on fixe le seuil à 0.68
# 

scoreMarquageCluster<-function(markersDE,ObjetSeurat,clusters="all",seuil="fixe", bonusExclu=T,filtreMin=2,Dfimportances=NULL){
  library(stringr)
  if (seuil=="fixe"){
    seuilPval=deterSeuilPval(ncol(ObjetSeurat))
    mat<-as.matrix(ObjetSeurat[["RNA"]]@counts)
    seuilDeltaPourcentageDexpr=deterSeuilDtPct(sum(mat==0)/length(mat))
    seuilLogFC=0.68
    print(paste("seuil Pval :",seuilPval))
    AjusteSeuilAuCluster<-F
  }else if(seuil=="interCluster"){
    seuilPval=quantile(markersDE$p_val_adj,0.15)
    seuilDeltaPourcentageDexpr=quantile(abs(markersDE$pct.1-markersDE$pct.2),0.9)
    seuilLogFC=quantile(markersDE$avg_logFC,0.85)
    print(paste("seuil Pval :",seuilPval))
    AjusteSeuilAuCluster<-F
  }else if (seuil=="intraClusterFixe"){
    #pour ne pas pénaliser les resolutions plus grande, qui ont donc plus de cluster et moins de cellules a l'interieur
    #seuilPval=deterSeuilPval(ncol(ObjetSeurat)) #!a mettre dans future version
    AjusteSeuilAuCluster<-T
    mat<-as.matrix(ObjetSeurat[["RNA"]]@counts)
    seuilDeltaPourcentageDexpr=deterSeuilDtPct(sum(mat==0)/length(mat))
    seuilLogFC=0.68
    
  }else if(seuil =="intraCluster"){
    #! a faire
  }
  
  genes_deja_passes<-c()
  BonusExcluTot<-1 #si marque qu'un cluster
  BonusExcluPartiel<-0.5 #si marque qu'un cluster avec score>1
  
  if(clusters=="all"){
    if (is.numeric(markersDE$cluster)){
      clusters<-0:tail(markersDE$cluster,1)
    }else if (is.character(markersDE$cluster)){
      if(all(str_length(markersDE$cluster)<3)){
        clusters<-as.numeric(levels(as.factor(markersDE$cluster)))
      }else{
        clusters<-unique(markersDE$cluster)
      }
      
    }else{
      clusters<-1:tail(markersDE$cluster,1)-1
    }
    
  }
  
  
  if(is.null(Dfimportances)){
    importanceMarkers<-data.frame(row.names=unique(markersDE$gene))
  }else {
    importanceMarkers<-Dfimportances
  }
  print(paste("seuil dt.pct :",seuilDeltaPourcentageDexpr))
  print(paste("seuil logFC :",seuilLogFC))
  #recuperer un dataframe contenant tout les gènes (nbdecluster*15) annoté
  #ajouter une colonne par cluster avec l'importance du gène en tant que marqueur de ce cluster 
  for (num_cluster in clusters ){
    
    if (AjusteSeuilAuCluster){
      ncells<-sum(Idents(ObjetSeurat)==num_cluster)
      seuilPval<-deterSeuilPval(ncells) #! verifier que seuil pas trop light,
      #et faire plutot une droite regressions nbCluster ~ pval, utiliser les coefs pour normaliser valeurs
      print(paste("seuil Pval pour le cluster", num_cluster,":",seuilPval))
    }
    
    for (gene in unique(markersDE$gene)){  
      importance_gene<-0
      
      if (gene %in% subset(markersDE,cluster==num_cluster)$gene){
        
        importance_gene<-1
        
        if (gene %in% genes_deja_passes){
          
          importanceMarkers[gene,"dans_clusters_"]<-paste(importanceMarkers[gene,"dans_clusters_"],num_cluster,sep = ", ")
          
        }else{
          importanceMarkers[gene,"dans_clusters_"]<-num_cluster
        }
        
        #recuperer le numero de ligne du resultats de signif marker
        ligne<-which(markersDE$cluster==num_cluster & markersDE$gene ==gene)
        
        #on recupere pct1, pct2, pval et logFC
        
        pct1<-markersDE[ligne,"pct.1"]
        pct2<-markersDE[ligne,"pct.2"]
        pval<-markersDE[ligne,"p_val_adj"]
        logFC<-markersDE[ligne,"avg_logFC"]
        
        
        #calculer l'importance de ce gène en tant que marqueur du cluster 
        #si pct.1-pct.2 >0.4  pval>10-50  logFC>1 --> +1 a chaque fois
        
        if(abs(pct1-pct2)>seuilDeltaPourcentageDexpr){
          importance_gene<-importance_gene+1
        }
        if (pval<seuilPval){
          importance_gene<-importance_gene+1
        }
        if(abs(logFC)>seuilLogFC){
          importance_gene<-importance_gene+1
        }
        
        #on met dans le score d'importance du gène pour ce cluster
        importanceMarkers[gene,paste("cluster",num_cluster,sep="")]<-importance_gene
        
        genes_deja_passes<-c(genes_deja_passes,gene)
        
      }else{
        importanceMarkers[gene,paste("cluster",num_cluster,sep="")]<-importance_gene
      }
      
    }
  }
  #on met un bonus si marque qu'un cluster
  if (bonusExclu){
    for (num_cluster in clusters){
      cluster<-paste("cluster",num_cluster,sep="")
      
      GenesRelevantDansCluster<-rownames(importanceMarkers)[importanceMarkers[,cluster] >= filtreMin] 
      if(length(GenesRelevantDansCluster)>0){
        score_markers<-importanceMarkers[GenesRelevantDansCluster,sapply(importanceMarkers,is.numeric)]
        
        marqueursExclu<-rownames(score_markers)[which(rowSums(score_markers!=0)==1)]
        
        importanceMarkers[marqueursExclu,cluster]<-importanceMarkers[marqueursExclu,cluster]+BonusExcluTot
        
        marqueursRelevanceExclu<-setdiff(rownames(score_markers)[rowSums(score_markers>1)==1],marqueursExclu)
        importanceMarkers[marqueursRelevanceExclu,cluster]<-importanceMarkers[marqueursRelevanceExclu,cluster]+BonusExcluPartiel
        
        
      }
      
      
    }
  }
  GenesRelevantDansClusters<-c()
  for (num_cluster in clusters){
    cluster<-paste("cluster",num_cluster,sep="")
    
    GenesRelevantDansCluster<-rownames(importanceMarkers)[importanceMarkers[,cluster] >= filtreMin]
    
    GenesRelevantDansClusters<-union(GenesRelevantDansClusters,GenesRelevantDansCluster)
  }
  
  #on filtre pour recuperer que le genes qui ont un score >= filtreMin
  importanceMarkers<-importanceMarkers[GenesRelevantDansClusters,]
  return(importanceMarkers)
}
