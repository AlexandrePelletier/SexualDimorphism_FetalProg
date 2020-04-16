

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

resLocisToGenes<-function(resLocis, withNCpGTot=FALSE){
  genes<-na.omit(unique(resLocis$gene))
  
  if(withNCpGTot){
    annot<-read.csv2("../../ref/annotation_CpG_HELP_ALL_070420.csv",row.names = 1)
    df<-data.frame(row.names =genes,nCpG=table(resLocis$gene)[genes],nCpGTot=table(annot$gene)[genes])
    colnames(df)<-c("gene","nCpG","Asuppr","nCpGtot")
    df<-df[,names(df)!="Asuppr"]
  }else{
    df<-data.frame(row.names =genes,nCpG=table(resLocis$gene)[genes])
    colnames(df)<-c("gene","nCpG")
  }
  
  return(df)
}



makeGeneList<-function(genes,res,tradInENTREZID=F,score="FC",aggregFUN=mean){
  if(tradInENTREZID){
    library(clusterProfiler)
    gene.df <- bitr(genes, fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)
    genes<-gene.df$ENTREZID
    geneList<-rep(0,length(genes))
    names(geneList)<-as.character(genes)
  }else{
    geneList<-rep(0,length(genes))
    names(geneList)<-genes
  }
  
    for(i in 1:length(genes)){
      
      if(tradInENTREZID){
        gene<-gene.df$SYMBOL[gene.df$ENTREZID==genes[i]]
      }else{
        gene<-genes[i]
      }
    
      geneList[i]<-aggregFUN(res[res$FC>0&res$gene==gene,score])
      
    }
    return(sort(geneList,decreasing = T))
}


tr<-function(ids_sepBySlash,tradEntrezInSymbol=F){
  library(clusterProfiler)
  
  IDs<-c(strsplit(ids_sepBySlash,"/")[[1]])
  if(tradEntrezInSymbol){
    return(bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL)
  }else{
    return(IDs)
  }
  
}
