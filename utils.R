
##utils everywhere
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

paste1<-function(vec){
  return(paste(vec,collapse = "/"))
}

##function explor data before modeling
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


##fct for res modeling
checkSignifPval<-function(res,pvalsPermut,cutoffStopIter=0.05,cutoffSigCol="q5",flag=T,colPval='pval'){
  
  res<-res[order(res[,colPval]),]
  locisToCheck<-rownames(res)
  pvalsPermut<-pvalsPermut[locisToCheck,]
  pctNonSig<-0
  i<-1
  nNonSig<-1
  locisPassPas<-c()
  while(pctNonSig<cutoffStopIter){
    if(i<=nrow(res)){
      if(res[i,colPval]>pvalsPermut[i,cutoffSigCol]){
        nNonSig<-nNonSig+1
        pctNonSig<-nNonSig/i
        
        if(flag){
          res[i,"PassPermut"]<-F
        }else{
          locisPassPas<-c(locisPassPas,i)
        }
        
      }else{
        if(flag){
          res[i,"PassPermut"]<-T
        }
      }
      
      i<-i+1
    }else {
      print(paste("le set de locis Sig est fiable, le pct de locis ne passant pas le cutoffPermut=",pctNonSig))
      if(flag){
        return(res)
      }else{
        return(locisToCheck[locisPassPas])
      }
      
    }
      
  }
  print(paste("le set de locis Sig est fiable jusqu'au locis ",i,"(pval=",res[i,colPval],")"))
  if(flag){
    return(res)
  }else{
    return(locisToCheck[locisPassPas])
  }
  
    
}


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

#locis > genes
aggregLocis<-function(gene,resLocis,colname,FUN=min,returnLocis=FALSE,colGene="gene"){
  if(returnLocis){
    return(rownames(resLocis)[which((resLocis[,colGene]==gene)&(resLocis[,colname]==FUN(resLocis[which(resLocis[,colGene]==gene),colname])))])
  }
  return(FUN(resLocis[which(resLocis[,colGene]==gene),colname]))
  
}


resLocisToGenes<-function(resLocis,withColLocis=TRUE, withNCpGTot=FALSE,aggregPval.FUN=NULL,aggregFC.FUN=NULL,colID="id",colGene="gene",colPval="pval",colFC="FC"){
  genes<-na.omit(unique(resLocis[,colGene]))
  df<-data.frame(row.names =genes,
                 nCpG=table(resLocis$gene)[genes])
  colnames(df)<-c("gene","nCpG")
  if(withColLocis==TRUE|withColLocis=="all"){
    df$locis<-sapply(genes,aggregLocis,resLocis,colID,paste1)
  }
  
  if(withNCpGTot){
    annot<-read.csv2("../../ref/annotation_CpG_HELP_ALL_070420.csv",row.names = 1)
    df<-data.frame(row.names =genes,df,nCpGTot=table(annot$gene)[genes])
    colnames(df)[ncol(df)-1:ncol(df)]<-c("Asuppr","nCpGtot")
    df<-df[,names(df)!="Asuppr"]
  }else{
    
  }
  
  if(!is.null(aggregPval.FUN)){
    print(paste("aggregation des", colPval))
    df[,colPval]<-sapply(genes,aggregLocis,resLocis,colPval,aggregPval.FUN)
    if(withColLocis=="all"){
      df[,paste0(colPval,".id")]<-sapply(genes,aggregLocis,resLocis,colPval,aggregPval.FUN,returnLocis=T)
    }
    
  }
  
  
  if(!is.null(aggregFC.FUN)){
    print(paste("aggregqtion des", colFC))
    df[,colFC]<-sapply(genes,aggregLocis,resLocis,colFC,aggregFC.FUN)
    if(withColLocis=="all"){
      df[,paste0(colPval,".id")]<-sapply(genes,aggregLocis,resLocis,colFC,aggregFC.FUN,returnLocis=T)
    }
  }
  
  
  return(df)
}

#genes>locis
annotLocis<-function(resLocis,resGenes,annots=NULL,genes=NULL,rmLocisSansGenes=T){
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
  if(rmLocisSansGenes){
    locisWithGenes<-rownames(resLocis)[poss]
    nLocis<-sum(!(rownames(resLocis)%in%locisWithGenes))
    print(paste(nLocis,"ont ete enleves car non reliés a un gène"))
    return(resLocis[poss,])
  }else{
    return(resLocis)
  }
  
  
}


# after modeling (GSEA..)
makeGeneList<-function(res,score="FC",returnRank=T,withResLocis=F,aggregLocisFUN=mean){
  if(withResLocis){
    genes<-na.omit(unique(res$gene))
  }else{
    genes<-rownames(res)
  }
  
  if(withResLocis){
    geneList<-rep(0,length(genes))
    for(i in 1:length(genes)){
      
      gene<-genes[i]
      
      geneList[i]<-aggregLocisFUN(res[res$gene==gene,score])
    }
    
    
    if(returnRank){
      geneList<-rank(geneList)
      
    }
    names(geneList)<-genes
    
  }else{
    if(returnRank){
      geneList<-rank(res[,score])
      
      
    }else{
      geneList<-res[,score]
    }
    names(geneList)<-genes
  }
  return(sort(geneList,decreasing = T))
}


tr<-function(ids_sepBySlash,retourne="all",sep="/",tradEntrezInSymbol=FALSE,uniqu=TRUE){
  IDs<-as.vector(strsplit(ids_sepBySlash,sep)[[1]])
  if(retourne=="all"){
    ret<-1:length(IDs)
  }else{
    ret<-retourne
  }
  if(tradEntrezInSymbol){
    require(clusterProfiler)
    library(org.Hs.eg.db)
    if(retourne=="all"){
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL)
    }
    else{
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL[ret])
    }
    
  }else{
    if(uniqu){
      return(unique(IDs[ret]))
    }else{
      return(IDs[ret])
    }
    
    
  }
  
}


#ensemble annot
putInChrRegFormat<-function(df,colChr="chr",colStart="start",colEnd="start",chrWithchr=TRUE){
  return(apply(df,1,function(x){
    if(chrWithchr){
      chr<-as.numeric(strsplit(x[colChr],"r")[[1]][2])
    }else{
      chr<-as.numeric(x[colChr])
    }
    
    start<-as.numeric(x[colStart])
    end<-as.numeric(x[colEnd])
    return(paste(chr,start,end,sep= ":"))
  }
  ))
  
}

regInReg<-function(ChrRangeVec1,ChrRangeVec2){
  return((ChrRangeVec1[1]==ChrRangeVec2[1]&ChrRangeVec1[2]>ChrRangeVec2[2]&ChrRangeVec1[3]<ChrRangeVec2[3]))
  
}
matchReg<-function(CpGsPos_InChrRegFormat,ChrReg){
  listCpG<-strsplit(CpGsPos_InChrRegFormat,":") 
  return(sapply(listCpG,regInReg,strsplit(ChrReg,":")[[1]]))
}
annotCpGResWithBMRes<-function(resCpG,resBM=NULL,chrPosCol="chrRegion",
                               attrs=c('chromosome_name','chromosome_start','chromosome_end','feature_type_name','regulatory_stable_id'),
                               annot='feature_type_name',
                               ensembl=NULL){
  library(stringr)
  if(is.null(resBM)){
    library(biomaRt)
    
    ensembl<-useEnsembl(biomart="ENSEMBL_MART_FUNCGEN", dataset="hsapiens_regulatory_feature", GRCh=37,version = 95)
    print(ensembl)
    resBM <- getBM(attributes=attrs,
                   filters = 'chromosomal_region',
                   values = resCpG$chrRegion,
                   mart = ensembl)
    
  }
  if("chrReg"%in%names(resBM)){
    print("chrReg format already here")
    
  }else{
    resBM$chrReg<-putInChrRegFormat(df=resBM,
                                    colChr = "chromosome_name",colStart ="chromosome_start",
                                    colEnd = "chromosome_end",chrWithchr = F)
  }
  resCpG[,annot]<-""
  features<-as.vector(unique(resBM[,annot]))
  resCpG[,features]<-""
  for(ChrReg in unique(resBM$chrReg)){
    
    locis<-rownames(resCpG)[matchReg(resCpG[,chrPosCol],ChrReg)]
    features<-resBM[resBM$chrReg==ChrReg&!duplicated(resBM$chrReg),annot]
    resCpG[locis,annot]<-paste0(resCpG[locis,annot],paste0(resBM[resBM$chrReg==ChrReg&!duplicated(resBM$chrReg),annot],sep="/"))
    resCpG[locis,features]<-paste(ChrReg,resBM[resBM$chrReg==ChrReg&!duplicated(resBM$chrReg),'regulatory_stable_id'],sep = "/")
    
  }
  
  resCpG[,annot]<-str_sub(resCpG[,annot],1,(str_length(resCpG[,annot])-1))
  return(resCpG)
}




