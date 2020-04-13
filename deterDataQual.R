

topVarFeatures<-function(mat,topSD){
  sds<-matrixStats::rowSds(mat,na.rm = T)
  return(rownames(mat)[order(sds,decreasing = T)][1:topSD])
}

deterQual<-function(matCpG, maxMethyl =5, minPct0=0.35){
  #parmi les locis avec pct0>, quelle pct de locis avec vrais zeros ?
  #1) remplace NA par 0
  matCpG[is.na(matCpG)]<-0
  #2) recupère mat avec plus de 35% de 0 dans locis
  matCpG<-matCpG[((rowSums(matCpG==0)/ncol(matCpG))>minPct0),]
  pctLocisAvecVraisZeros<-sum(rowSums(matCpG>0&matCpG<maxMethyl)>0)/nrow(matCpG)
  print(paste("pct Locis avec Vrais zeros =",round(pctLocisAvecVraisZeros,3)))
  return(pctLocisAvecVraisZeros)
}

deterQual2<-function(mat,df.covars,covarNum="Group_Complexity",pcTestes=5){
  #quelle corralation avec library complexity?
  res<-list(p=1, r2=0,PC=0, pctPC=0)
  
  mat[is.na(mat)]<-0
  pca<-prcomp(t(mat))
  pc<-pca$x
  library<-as.numeric(df.covars[rownames(pc),covarNum])
   #explique que 3.7% !
  #library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
  for (i in 1:pcTestes){
    
    resLm<-lm(pc[,i]~library)
    
    if(anova(resLm)$Pr[1]<res$p & summary(resLm)$adj.r.squared>res$r2){
      
      res$p<-anova(resLm)$Pr[1]
      res$r2<-summary(resLm)$adj.r.squared
      res$PC<-i
      res$pctPC<-round(pca$sd[i]^2/sum(pca$sdev^2),3)
      print(paste("PC",i," (",res$pctPC*100,"% de la variance","a R2 avec",covarNum,"=",round(summary(resLm)$adj.r.squared,2),"et pval = 10^",log10(anova(resLm)$Pr[1])))
      
    }
    
  }
  
  

  
  return(res)
  
}




deterQual3<-function(mat,df.covars,topSD=10000,covar){
  if(is.numeric(topSD)){
    locis<-topVarFeatures(mat,topSD)
    mat<-mat[locis,]
  }
  mat<-na.omit(mat)
  dis <- dist(t(mat))
  dis_mat<-as.matrix(dis)
  distsGroupe<-c()
  pos<-0
  for(lvl in levels(as.factor(df.covars[,covar]))){
    samplesGroupe<-na.omit(rownames(df.covars)[df.covars[,covar]==lvl])
    #print(samplesGroupe)
    distsGroupe<-c(distsGroupe,rep(0,(length(samplesGroupe)*(length(samplesGroupe)-1)/2)))
    
    for(i in 1:(length(samplesGroupe)-1)){
      
      for(j in (i+1):length(samplesGroupe)){
        pos<-pos+1
        #print(paste(samplesGroupe[i],samplesGroupe[j]))
        distsGroupe[pos]<-dis_mat[samplesGroupe[i],samplesGroupe[j]]
        
      }
      
    }
    
  }
  
  return(mean(distsGroupe))
  
}



