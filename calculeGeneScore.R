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

CalcCpGScore<-function(res,cpg_genes=NULL){
  if(!is.null(cpg_genes)){
    print("merging limma res and cpgs annotations...")
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
    
  }
  
  if(all(c("LinksWeight","RegWeight")%in%colnames(res))){
    res[,CpGScore:=((-log10(pval)/4)*meth.change)*RegWeight*LinksWeight]
  }else{
    print("need calculate linksWeight and RegWeight of CpGs before.")
  }
  
  
  return(res)
  
}

CalcGeneScore<-function(res,cpg.regs_ref=NULL,pvalSig=10^-3,sumToGene=FALSE,test=FALSE){
  
  if(!"CpGScore"%in%colnames(res)){
    print("calculating CpGScore...")
    res<-CalcCpGScore(res,cpg.regs_ref)
  }
  
  print("calculating GeneScore...")
  print("1) add column number of CpG by Gene")
  res[,nCpG.Gene:=.N,by=.(gene)]
  
  res[,nCpGSig.Gene:=sum(pval<pvalSig),by=.(gene)]
  
  print("2) the nCpGWeight : (1/sum(1/(abs(CpGScore)+1)))^(1/4)")
  res[,nCpGWeight:=(1/sum(1/(abs(CpGScore)+1)))^0.2,by="gene"]
  
  print("3) the GeneScore : sum(CpGScore)*nCpGWeight")
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


