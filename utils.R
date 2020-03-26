

pctPC<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)*100
  names(pct.varPCs)<-rngPCs
  return( pct.varPCs)
}

correl<-function(x,y,ret="pval",verbose=T){
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


# trunkName<-function(names,maxLeng=5,n.Mot=3){
#   library(stringr)
#   str_split()
#   str_sub(names,)
# }