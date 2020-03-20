

topVarFeatures<-function(mat,topSD){
  sds<-matrixStats::rowSds(mat,na.rm = T)
  return(rownames(mat)[order(sds,decreasing = T)][1:topSD])
}

deterQual<-function(matCpG, maxMethyl =5, minPct0=0.7){
  #parmi les locis avec pct0>, quelle pct de locis avec vrais zeros ?
  #1) remplace NA par 0
  matCpG[is.na(matCpG)]<-0
  #2) recupère mat avec plus de 70% de CpG
  matCpG<-matCpG[((rowSums(matCpG==0)/ncol(matCpG))>minPct0),]
  
  return(sum(rowSums(matCpG>0&matCpG<maxMethyl)>0)/nrow(matCpG))
}


deterQual2<-function(mat,df.covars,topSD=10000,covar){
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



