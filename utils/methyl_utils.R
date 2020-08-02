

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
  print(paste("le set de locis Sig est fiable jusqu'au locis num",i,"(pval=",res[i,colPval],")"))
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
    if(!(colID%in%colnames(resLocis))){
      resLocis[,colID]<-rownames(resLocis)
    }
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


extr<-function(x){
  r<-range(x)
  return(r[which(abs(r)==max(abs(r)))])
}


plotMeth<-function(cpgs,
                   factor="Group_Sex",
                   levels=c("C_F","C_M","L_F","L_M"),
                   methyl_df=NULL,
                   batch=NULL,
                   plot="boxplot",
                   group="locisID",
                   wrap=FALSE){
  library(ggplot2)
  library(data.table)
  if(is.null(methyl_df)){
    methyl_df<-fread("../../ref/2020-05-25_methyl_data_before_limma.csv")
  }
  if(is.null(batch)){
    batch<-fread("../../ref/cleaned_batch_CD34_library_date_220620.csv")
    
    keep<-batch[[factor]]%in%levels
    batch<-batch[as.vector(keep)]
  }
  
  if(is.data.table(cpgs)){
    dt<-copy(cpgs)
    cpgs<-dt$locisID
  }else{
    dt<-NA
  }
  
  if(is.numeric(cpgs)){

        cpgs_data<-methyl_df[locisID%in%cpgs]
        #need transformer methyl_df en df with sample,locisID,meth,group..
        samples<-colnames(cpgs_data)[colnames(cpgs_data)!="locisID"]
        
        cpgs_data2<-data.table(expand.grid(sample=samples,
                               locisID=cpgs_data$locisID))
        cpgs_score<-Reduce(rbind,lapply(samples, function(sampleID){
          
          return(cpgs_data[,unmeth:=.SD,.SDcols=sampleID][,sample:=sampleID][,.(locisID,sample,unmeth)])
        }))
        
        cpgs_data2<-merge(cpgs_data2,cpgs_score,by=c("locisID","sample"))
        
        
        cpgs_data2<-cpgs_data2[!is.na(unmeth)]
        cpgs_data_batch<-merge(cpgs_data2,batch[match(cpgs_data2$sample,sample)][!is.na(sample)],by="sample",allow.cartesian=TRUE)
        
        if(any(!is.na(dt))){
          cpgs_data_batch<-merge(cpgs_data_batch,dt,by="locisID")
          
        }
        
        ord<-as.character(unique(lapply(cpgs_data_batch[,.SD,.SDcols=group],sort)[[1]]))
        
        cpgs_data_batch[,(group):=lapply(.SD,function(x)factor(as.character(x),levels=ord)),.SDcols=group]
        
        
        if(plot=="boxplot"){
          if(wrap){
            return(ggplot(cpgs_data_batch)+geom_boxplot(aes_string(group,"unmeth"))+facet_wrap(factor))
          }
          return(ggplot(cpgs_data_batch)+geom_boxplot(aes_string(factor,"unmeth",fill=group))+scale_color_manual(breaks=as.character(ord)))
          
        }else if(plot=="jitter"){
          if(wrap){
            return(ggplot(cpgs_data_batch)+geom_jitter(aes_string(group,"unmeth"),width = 0.25)+facet_wrap(factor))
          }else{
            return(ggplot(cpgs_data_batch)+geom_jitter(aes_string(factor,"unmeth",color=group),width = 0.25))
          }
          
          
          
        }
      
    
    }else{
    print("enter numeric value")
  }
}


RunMethAnalysis<-function(methyl_df,batch_filtered,formule,cpg.regs_ref,compas_df,sumToGene=FALSE,verbose=TRUE){
  library(limma)
  
  if("data.table"%in%class(methyl_df)){
    print("transforming methyl_data in dataframe")
    methyl_df<-data.frame(methyl_df,row.names = methyl_df$locisID)
    print(head(methyl_df))
  }
  
  if("data.table"%in%class(batch_filtered)){
    print("transforming batch metadata in dataframe")
    batch_F<-data.frame(batch_filtered)
    print(head(batch_filtered))
  }
  
  design<-model.matrix(formule,data = batch_filtered)

  fit <- lmFit(methyl_df[,batch_filtered$sample], design)  
  
  cont.matrix <- makeContrasts(contrasts = compas_df$compa,
                               levels=design)
  print(colnames(cont.matrix))
  print(compas_df$abbrev_compa)
  colnames(cont.matrix)<-compas_df$abbrev_compa
  
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2) 
  if(nrow(compas_df)>1){
    
    res_list<-lapply(compas_df$abbrev_compa,function(compa){
      res<-topTable(fit2,coef=compa,n =Inf)
      
      res<-CalcGeneScore(res,cpg.regs_ref,sumToGene=sumToGene,verbose = verbose)
      return(res[order(-GeneScore)])
    })
    names(res_list)<-compas_df$abbrev_compa
    return(res_list)
    
  }else{
    res<-topTable(fit2,coef=compas_df$abbrev_compa,n =Inf)
    res<-CalcGeneScore(res,cpg.regs_ref,sumToGene=sumToGene,verbose=verbose)
    return(res[order(-GeneScore)])
  }
  
  
  
}

#CALCUL GENESCORE

CalcCpGWeights<-function(cpgs_genes){
  
  cpgs_genes[,TypeScore:=sapply(type,function(x){
    if(x%in%c(4,6)){
      score<-1
    }else if(x==5){
      score<-0.75
    }else if(x%in%1:3){
      score<-0.5
    }else{
      score<-0
    }
    return(score)
  })]
  
  
  #   EnsRegScore{0,0.25,0.5,0.75,1}
  cpgs_genes[,EnsRegScore:=sapply(feature_type_name,function(x){
    vecX<-strsplit(as.character(x),"/")[[1]]
    if(any(c("CTCF Binding Site","Promoter","Enhancer")%in%vecX)){
      score<-0.5
    }else if(any(c("Open chromatin","Promoter Flanking Region")%in%vecX)){
      score<-0.25
    }else{
      score<-0
    }
    if("TF binding site"%in%vecX){
      score<-score+0.5
    }
    return(score)
  })]
  
  # Regulatory weight:
  cpgs_genes[,RegWeight:=(0.5+1.5*((TypeScore+EnsRegScore)/2))]
  
  # > 3) Links Weight[0.1-1] = 0.1+0.9*max(LinkScore [0-1])
  #   - for CpG associated to a gene based on TSS proximity  
  cpgs_genes[in_eQTR==FALSE,LinkScore:=sapply(abs(tss_dist),function(x){
    if(x<1000){
      score<-1
    }else if(x<20000){
      score<-0.5+0.5*sqrt(1000/x)
    }else{
      score<-0.5*sqrt(20000/x)
      
    }
    return(score)
  })]
  
  cpgs_genes[in_eQTR==TRUE,LinkScore:=(0.5+(avg.mlog10.pv.eQTLs/50))/1.5] 
  
  cpgs_genes[in_wb_eQTR==TRUE,in_meta_eQTR:=FALSE]
  cpgs_genes[in_meta_eQTR==TRUE,in_wb_eQTR:=FALSE]
  
  cpgs_genes[,inBoth_eQTR:=any(in_meta_eQTR==TRUE)&any(in_wb_eQTR==TRUE),by=c("locisID","gene")]
  
  cpgs_genes[inBoth_eQTR==TRUE,LinksWeight:=1,by=c("locisID","gene")]
  
  cpgs_genes[inBoth_eQTR==FALSE,LinksWeight:=max(LinkScore),by=c("locisID","gene")]
  
  cpgs_genes[inBoth_eQTR==TRUE,LinksWeight:=1]
  
  return(cpgs_genes)
}

CalcCpGScore<-function(res,cpg_genes=NULL,verbose=TRUE){
  if(!is.null(cpg_genes)&
     !all(c("LinksWeight","RegWeight")%in%colnames(res))){
    if(verbose){
      print("merging limma res and cpgs annotations...")
    }
    
    if(all(c("FC","pval")%in%colnames(res))){
      res$logFC<-res$FC
      res$P.Value<-res$pval
    }else if("meth.change"%in%colnames(res)){
      res$logFC<-res$meth.change
      
    }
    if(!("locisID"%in%colnames(res))){
      res$locisID<-rownames(res)
    }
    res<-data.table(locisID=as.numeric(res$locisID),
                    meth.change=res$logFC,
                    pval=res$P.Value)[order(locisID)]
    
    res<-merge(res,cpg_genes,by="locisID")
    
  }else if(is.null(cpg_genes)&
           !all(c("LinksWeight","RegWeight")%in%colnames(res))){
    library(data.table)
    cpg_genes<-fread("../../ref/2020-06-29_All_CpG-Gene_links.csv")
    res<-merge(res,cpg_genes,by="locisID")
  }
  
  if(all(c("LinksWeight","RegWeight")%in%colnames(res))){
    res[,CpGScore:=(-log10(pval)/4*meth.change)*RegWeight*LinksWeight]
  }else{
    print("need calculate linksWeight and RegWeight of CpGs before.")
  }
  
  
  return(res)
  
}

CalcGeneScore<-function(res,cpg.regs_ref=NULL,pvalSig=10^-3,sumToGene=FALSE,test=FALSE,calcCpGScore=FALSE,verbose=TRUE){
  
  if((!"CpGScore"%in%colnames(res))|calcCpGScore==TRUE){
    if(verbose){
      print("calculating CpGScore...")
    }
    
    res<-CalcCpGScore(res,cpg.regs_ref,verbose = verbose)
  }
  
  if(verbose){
    print("calculating GeneScore...")
    print("1) add column number of CpG by Gene")
  }
  
  
  res[,nCpG.Gene:=.N,by=.(gene)]
  
  res[,nCpGSig.Gene:=sum(pval<pvalSig),by=.(gene)]
  if(verbose){
    print("2) the nCpGWeight : (1/sum(1/(abs(CpGScore)+1)))^(1/4)")
  }
  
  res[,nCpGWeight:=(1/sum(1/(abs(CpGScore)+1)))^(1/4),by="gene"]
  
  if(verbose){
    print("3) the GeneScore : sum(CpGScore)*nCpGWeight")
  }
  res[,GeneScore:=sum(CpGScore)*nCpGWeight,by="gene"]
  
  if(test==TRUE){
    library(ggplot2)
    library(patchwork)
    print("plot GeneScore ~ nCpG")
    p1<-ggplot(unique(res,by="gene"))+
      geom_boxplot(aes(x = as.factor(nCpG.Gene),y =nCpGWeight )) 
    
    p2<-ggplot(unique(res,by="gene"))+
      geom_boxplot(aes(x = as.factor(nCpG.Gene),y =GeneScore )) 
    p_all<-p1+p2
    print(p_all)
  }
  
  if(sumToGene){
    return(unique(res[order(-GeneScore,pval)],by="gene"))
  }else{
    return(res[order(-GeneScore,pval)])
  }
}

#deter genescore cutoff x / deter EpigenAffGene (EAG) : 
# 1000 permuts of all samples CTRL / LGA => LIMMA => pval and FC => df GeneScore => save pct genescore permuts >= geneScore obs
#=> gene score cutoff = Accept Gene in EAF if appeared < x%  => genesF and genesM spe
#make a permutFunction

prepBatchDf<-function(batch,varNumToModel,varFacToModel){
  
  varToModel<-c(varNumToModel,varFacToModel)
  
  sample_F<-batch$sample[!(apply(is.na(batch[,..varToModel]),1,any))]
  
  print(paste(length(sample_F),"after filtration for NA "))
  
  batch<-batch[sample%in%sample_F,c("sample",..varToModel),]
  
  batch[,(varNumToModel):=lapply(.SD,as.numeric),.SDcols=varNumToModel]
  batch[,(varFacToModel):=lapply(.SD,as.factor),.SDcols=varFacToModel]
  
  return(batch)
}
permutGeneScore<-function(res_to_permut,methyl_df,cpg.regs_ref,batch,var_to_permut,n_perm=1000,seed=1234,
                          varNumToModel=c("Mat.Age"),varFacToModel=c("Group_Sex",'batch',"latino","Group_Complexity_Fac"),
                          formule= ~0 + Group_Sex  + batch  + latino + Mat.Age + Group_Complexity_Fac,
                          verbose=FALSE
){
  library(data.table)
  library(stringr)
  
  if(!is.list(res_to_permut)){
    print("need results in a  list named by comparison")
  }
  batch<-prepBatchDf(batch,varNumToModel,varFacToModel)
  
  compas_df<-data.frame(compa=unlist(lapply(strsplit(names(res_to_permut),"-"),function(x)paste(paste0(var_to_permut,x),collapse = "-"))),
                        abbrev_compa=names(res_to_permut))
  print(compas_df)
  
  
  res_all<-data.table(gene=sort(unique(res_to_permut[[1]]$gene)))
  
  
  batch_sim<-copy(batch)
  print(paste("set seed to",seed))
  set.seed(seed)
  
  print(paste("permutation de",var_to_permut,n_perm,"fois..."))
  for(i in 1:n_perm){
    print(paste(i,"/",n_perm))
    
    batch_sim[,Group_Sex:=sample(Group_Sex)]
    res_list<-RunMethAnalysis(methyl_df = methyl_df,
                              batch_F = batch_sim,
                              formule = formule,
                              compas_df =compas_df,
                              cpg.regs_ref = cpg.regs_ref,
                              sumToGene = T,
                              verbose=verbose)
    
    
    res_list<-lapply(names(res_list),function(compa){
      res<-res_list[[compa]]
      res<-res[order(gene)][,.(gene,GeneScore)]
      col<-paste0(str_sub(compa,str_length(compa)),i)
      return(res[,(col):=as.integer(GeneScore)][,-"GeneScore"])
    })
    
    res_list<-Reduce(merge,res_list)
    res_all<-merge(res_all,res_list,by="gene")
    
  }
  
  res_to_permut<-lapply(1:length(res_to_permut),function(i){
    
    res<-unique(res_to_permut[[i]],by="gene")
    compa<-names(res_to_permut)[i]
    cols<-paste0(str_sub(compa,str_length(compa)),1:n_perm)
    col<-paste0("pval",n_perm,"perm")
    res<-res[order(gene)][,(col):=rowSums(..res_all[order(gene)][,.SD,.SD=cols]>GeneScore)/(..n_perm)]
    cols<-c("gene",col)
    return(merge(res_to_permut[[i]],res[,.SD,.SD=cols],by="gene"))
  })
  names(res_to_permut)<-compas_df$abbrev_compa
  return(res_to_permut)
}


getGenesKEGGPathw<-function(pathID){
  library(KEGGREST)
  g<-keggGet(pathID)[[1]]$GENE
  g<-g[1:length(g)%%2==0]
  return(as.vector(sapply(g,function(x)strsplit(x,";")[[1]][1])))
}




