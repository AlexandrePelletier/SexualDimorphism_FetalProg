

pctPC<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)*100
  names(pct.varPCs)<-rngPCs
  return( pct.varPCs)
}

correl<-function(x,y,ret="all",verbose=T){
  if(is.numeric(x)){
    if(verbose){
      print("linear modeling ")
    }
    
    res<-lm(x~y)
    if(verbose){
      print(summary(res))
    }
  if(ret=="r2"){
    
    return(summary(res)$adj.r.squared)
  }else if(ret=="pval"){
    return(anova(res)$Pr[1])
    
  }else{
    return(
      list(p=anova(res)$Pr[1],r2=summary(res)$adj.r.squared)
      )
  }
    
  }else if(all(sapply(list(x,y),is.factor))){
    if(verbose){
      print("Chi-squared independance test")
    }
    
    tableF<-table(x,y)
    test<-chisq.test(tableF)
    if(verbose){
      print(test)
    }
    
    return(test$p.value)
  }
  
}


trunkName<-function(names,maxLeng=4,n.mot=NULL){
  library(stringr)
  if(is.null(n.mot)){
    return(str_sub(names,1,maxLeng))
  }else if(is.numeric(n.mot)){
    liste_mots<-strsplit(names," ")
    
    vec2Mots<-rep("",length(names))
    for (j in 1:length(names)){
      mots<-liste_mots[[j]]
      mots<-mots[str_length(mots)>2]
      vecMots<-c()
      for(i in 1:n.mot){
        if(!is.na(mots[i])){
          vecMots<-c(vecMots,str_sub(mots[i],1,maxLeng)) 
        }
        
      }
      vec2Mots[j]<-paste(vecMots,collapse = ". ")
      
    }
    return(vec2Mots)
  }
  
  
  
  
  
  return(vecMots)
}