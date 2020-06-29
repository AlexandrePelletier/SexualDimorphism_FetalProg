
#2020-06-25_Def_GeneScore_WithSim_Data
script<-"2020-06-25_Def_GeneScore_WithSim_Data"

#see behaviour of my CpGScore/GeneScore with simulated

options(stringsAsFactors=F)
set.seed(12345)
library(data.table)
library(stringr)
# I) test my CpG Score 
#   1) create Simulated CpG
#for factors : all interesting factors, for numeric 5 : min q25 median q75 max
fc_pos<- c(0,15,30,60)
fc_pos_neg<-sort(c(-fc_pos,fc_pos,0))

pval <- c(1,0.1,0.001,10^-5,10^-8)

type <-as.factor(c("heterochr","transcribed","CRE","cre"))

ensReg1 <-c("CTCF Binding Site","Promoter", "Enhancer","Open chromatin","Promoter Flanking Region")
ensReg1
ensReg2<-c("",ensReg1,apply(expand.grid(ensReg1,c("TF binding site")),1,paste,collapse="/"),"TF binding site")
ensReg2<-as.factor(ensReg2)

in_eQTR<-c(FALSE,TRUE)
in_wb_eQTR<-c(FALSE,TRUE)
in_meta_eQTR<-c(FALSE,TRUE)

distTSS<-c(0,2000,10000,50000,10^6)
avg.mlog10.pv.eQTLs<-c(10,30,50)

sim_cpg_tss <- data.table(expand.grid("meth.change"=fc_pos_neg,"pval"= pval,
                                      "type"= type,"feature_type_name"=ensReg2,
                                      "distTSS"=distTSS,
                                      "in_eQTR"=F,"in_wb_eQTR"=F,
                                      "in_meta_eQTR"=F,"avg.mlog10.pv.eQTLs"=NA))
                                  


sim_cpg_eqtr <- data.table(expand.grid("meth.change"=fc_pos_neg,"pval"= pval,
                                       "type"= type,"feature_type_name"=ensReg2,
                                       "distTSS"=100000,
                                       "in_eQTR"=T,"in_wb_eQTR"=in_wb_eQTR,
                                   "in_meta_eQTR"=in_meta_eQTR,"avg.mlog10.pv.eQTLs"=avg.mlog10.pv.eQTLs))



sim_cpg<-rbind(sim_cpg_tss,sim_cpg_eqtr)
sim_cpg#36720 fake CpG



#calc cpg score :

sim_cpg[,TypeScore:=sapply(type,function(x){
  if(x=="CRE"){
    score<-1
  }else if(x=="cre"){
    score<-0.75
  }else if(x=="transcribed"){
    score<-0.5
  }else{
    score<-0
  }
  return(score)
})]


#   EnsRegScore{0,0.25,0.5,0.75,1}
sim_cpg[,EnsRegScore:=sapply(feature_type_name,function(x){
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
sim_cpg[,RegWeight:=(0.5+1.5*((TypeScore+EnsRegScore)/2))]

# > 3) Links Weight[0.1-1] = 0.1+0.9*max(LinkScore [0-1])
#   - for CpG associated to a gene based on TSS proximity  
sim_cpg[in_eQTR==FALSE,LinkScore:=sapply(abs(distTSS),function(x){
  if(x<1000){
    score<-1
  }else if(x<20000){
    score<-0.5+0.5*sqrt(1000/x)
  }else{
    score<-0.5*sqrt(20000/x)
    
  }
  return(score)
})]

sim_cpg[in_eQTR==TRUE,LinkScore:=(0.5+(avg.mlog10.pv.eQTLs/50))/1.5] 

sim_cpg[,inBoth_eQTL:=in_meta_eQTR==TRUE&in_wb_eQTR==TRUE]

sim_cpg[inBoth_eQTL==TRUE,LinksWeight:=1]

sim_cpg[inBoth_eQTL==FALSE,LinksWeight:=0.1+0.9*LinkScore]

sim_cpg[,CpGScore:=((-log10(pval)/4)*meth.change)*RegWeight*LinksWeight]

sim_cpg<-sim_cpg[order(-CpGScore)]
summary(sim_cpg$CpGScore)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   7.312  23.103  30.938 240.000 


sim_cpg[CpGScore>200] #en top cpg rank bien 

#en q75: comment avoir un score plus grand que 30 ?
head(sim_cpg[CpGScore>30][order(CpGScore)],60)
#ya til des cpgscore >30 mais pval / fc nul ?
sim_cpg[CpGScore>30&pval>0.01] #non
sim_cpg[CpGScore>30&meth.change<20] #oui ! mais parce que pval ++ signif, en pratique c'est pas trop possible

#RegWeight ya til des cpgscore >30 mais heterochromatine? 
sim_cpg[CpGScore>30&type=="heterochr"] #oui mais car meth change (min 30) ++ pval (min 10^-5) ++, et feature_type
levels(sim_cpg[CpGScore>30&type=="heterochr"]$feature_type_name) #il y a tout

#MinksWeight ya til des cpgscore >30 mais cpg très loin ? 
summary(sim_cpg[in_eQTR==F&CpGScore>30&distTSS==1000000 ]$CpGScore) 
#oui, donc enlevele 0.1+ : 
sim_cpg[inBoth_eQTL==FALSE,LinksWeight:=LinkScore]
summary(sim_cpg$LinksWeight)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07071 0.46667 0.73333 0.73521 1.00000 1.00000 

sim_cpg[,CpGScore:=((-log10(pval)/4)*meth.change)*RegWeight*LinksWeight]
#et maintenant ?
summary(sim_cpg[in_eQTR==F&CpGScore>30&distTSS==1000000 ]$CpGScore) #il y en a plus !
max(sim_cpg[in_eQTR==F&CpGScore>10&distTSS==1000000 ]$CpGScore) #=17

#Gene Score
sim_cpg[,locisID:=1:.N]
#gene with 1 cpg
sim_gene1<-copy(sim_cpg)
sim_gene1[,gene:=paste0("gene",1:.N,"_1")]
sim_gene1

#gene with 2, 4, 8, 16, 32, 100, 500
#with 2 :
MakePair<-function(ids,n.cpg,seed=1234){
  set.seed(seed)
  pairing<-rep(NA,length(ids))
  ids_toy<-ids
  while(length(ids_toy)>=n.cpg){
    pair<-sample(ids_toy,n.cpg)
    pos_pair<-which(ids%in%pair)
    pairing[pos_pair]<-paste(pair,collapse = "_")
    
    
    ids_toy<-ids_toy[!(ids_toy%in%pair)]
    
  }
  if(any(is.na(pairing))){
    pos_uniq<-which(is.na(pairing))
    pairing[pos_uniq]<-paste0(ids,"_")
  }
  
  return(pairing)
}


sim_gene2<-copy(sim_cpg)
sim_gene2[,gene:=MakePair(locisID,2)]

sim_gene2

#make 20 different combination

sim_gene2<-lapply(1:20, function(x){
  sim_gene<-copy(sim_cpg)
  return(sim_gene[,gene:=MakePair(locisID,2,seed = x)])
})
sim_gene2
sim_gene2<-Reduce(rbind,sim_gene2)

#gene with 4, 8, 16, 32, 100, 500
sim_gene5<-lapply(1:20, function(x){
  sim_gene<-copy(sim_cpg)
  return(sim_gene[,gene:=MakePair(locisID,5,seed = x)])
})

sim_gene5<-Reduce(rbind,sim_gene5)

sim_gene15<-lapply(1:20, function(x){
  sim_gene<-copy(sim_cpg)
  return(sim_gene[,gene:=MakePair(locisID,15,seed = x)])
})

sim_gene15<-Reduce(rbind,sim_gene15)



sim_gene50<-lapply(1:20, function(x){
  sim_gene<-copy(sim_cpg)
  return(sim_gene[,gene:=MakePair(locisID,50,seed = x)])
})

sim_gene50<-Reduce(rbind,sim_gene50)


sim_gene200<-lapply(1:20, function(x){
  sim_gene<-copy(sim_cpg)
  return(sim_gene[,gene:=MakePair(locisID,200,seed = x)])
})

sim_gene200<-Reduce(rbind,sim_gene200)

sim_genes<-Reduce(rbind,list(sim_gene1,
                       sim_gene2,
                       sim_gene5,
                       sim_gene15,
                       sim_gene50,
                       sim_gene200
                       ))
sim_genes

dir.create("analyses/sim_gene_score")
fwrite(sim_genes,"analyses/sim_gene_score/sim_data.txt")

#2020-06-29 
#calculate GeneScore
library(data.table)
sim_genes<-fread("analyses/sim_gene_score/sim_data.txt")
colnames(sim_genes)

CalcGeneScore<-function(res,cpg.regs_ref=NULL,pvalSig=10^-3,sumToGene=F){
 
  
  if(!is.null(cpg.regs_ref)){
    print("merging limma res and cpgs annotations...")
    res<-data.table(locisID=as.numeric(rownames(res)),
                    meth.change=res$logFC,
                    pval=res$P.Value)[order(locisID)]
    
    res<-merge(res,cpg.regs_ref,by="locisID")
    
  }
  
  if(!"CpGScore"%in%colnames(res)){
    print("calculating CpGScore...")
    res[,CpGScore:=((-log10(pval)/4)*meth.change)*RegWeight*LinksWeight]
  }
  
  print("calculating GeneScore...")
  print("1) the nCpGWeight : (1/sum(1/(abs(CpGScore)+1)))^(1/3)")
  res[,nCpGWeight:=(1/sum(1/(abs(CpGScore)+1)))^(1/3),by="gene"]
  print("2) the GeneScore : sum(CpGScore)*nCpGWeight")
  res[,GeneScore:=sum(CpGScore)*nCpGWeight,by="gene"]
  
  #bonus annot
  print("add column number of CpG by Gene")
  res[,nCpG.Gene:=.N,by=.(gene)]
  
  res[,nCpGSig.Gene:=sum(pval<pvalSig),by=.(gene)]
  
  if(sumToGene){
    return(unique(res[order(-GeneScore,pval)],by="gene"))
  }else{
    return(res[order(-GeneScore,pval)])
  }
}



sim_genes<-CalcGeneScore(sim_genes,)
sim_genes<-sim_genes[nCpG.Gene%in%c(1,2,5,15,50,200)]

#compares genes hi : 
#in function of the nCpG
library(ggplot2)
ggplot(unique(sim_genes[GeneScore>100],by="gene"))+
  geom_boxplot(aes(x = as.factor(nCpG.Gene),y =GeneScore )) 
#♣ccl : help to much -gene with only 1-2 cpg


#need optimize ncpg weight : 
#objectif : best GeneScore for 5-15 nCpG related genes 
#tests :
CalcGeneScore<-function(res,cpg.regs_ref=NULL,pvalSig=10^-3,sumToGene=FALSE,test=FALSE){
  if(!is.null(cpg.regs_ref)){
    print("merging limma res and cpgs annotations...")
    res<-data.table(locisID=as.numeric(rownames(res)),
                    meth.change=res$logFC,
                    pval=res$P.Value)[order(locisID)]
    
    res<-merge(res,cpg.regs_ref,by="locisID")
    
  }
  if(!"CpGScore"%in%colnames(res)){
    print("calculating CpGScore...")
    res[,CpGScore:=((-log10(pval)/4)*meth.change)*RegWeight*LinksWeight]
  }
  print("calculating GeneScore...")
  print("1) add column number of CpG by Gene")
  res[,nCpG.Gene:=.N,by=.(gene)]
  
  res[,nCpGSig.Gene:=sum(pval<pvalSig),by=.(gene)]
  
  print("2) the nCpGWeight : (1/sum(1/(abs(CpGScore)+1)))^(1/4)")
  res[,nCpGWeight:=(1/sum(1/(abs(CpGScore)+1)))^0.2,by="gene"]
  
  print("3) the GeneScore : sum(CpGScore)*nCpGWeight")
  res[,GeneScore:=sum(CpGScore)*nCpGWeight,by="gene"]
  
  #bonus annot
  
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
sim_genes<-CalcGeneScore(sim_genes,test = T)
#•fixed parameter 
#genes with all cpg with hi score (>50)

sim_genes[,genesHi:=all(CpGScore>30),by="gene"]
sim_genes[genesHi==T] #only genes with 1 and 2
ggplot(unique(sim_genes[genesHi==T],by="GeneScore"))+
  geom_boxplot(aes(x = as.factor(nCpG.Gene),y =GeneScore )) 

#genes with 1 cpg au max and the other au min

sim_genes[,genes1Hi:=sum(CpGScore>30)==1,by="gene"]
sim_genes[genes1Hi==T]
ggplot(unique(sim_genes[genes1Hi==T],by="GeneScore"))+
  geom_boxplot(aes(x = as.factor(nCpG.Gene),y =GeneScore )) 

#genes au moins 0.75 des cpg with hi cpgScore

sim_genes[,genes0.75Hi:=sum(CpGScore>30)>=nCpG.Gene*3/4,by="gene"]
sim_genes[genes0.75Hi==T] 

ggplot(unique(sim_genes[genes0.75Hi==T],by="GeneScore"))+
  geom_boxplot(aes(x = as.factor(nCpG.Gene),y =GeneScore )) 
#good car gene with 5 cpg have a median genescore gene with only 2

#genes au moins 0.5 des cpg with medium hi cpgScore (10)
sim_genes[,genes0.5Med:=sum(CpGScore>10)>=nCpG.Gene/2,by="gene"]
sim_genes[genes0.5Med==T] 

ggplot(unique(sim_genes[genes0.5Med==T],by="GeneScore"))+
  geom_boxplot(aes(x = as.factor(nCpG.Gene),y =GeneScore )) 
#good car gene with 15 cpg have a median genescore > que gene ac 2 et 5 cpg

#fianlly : 
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

  cpgs_genes[,inBoth_eQTL:=any(in_meta_eQTR==TRUE)&any(in_wb_eQTR==TRUE),by=c("locisID","gene")]

  cpgs_genes[inBoth_eQTL==TRUE,LinksWeight:=1,by=c("locisID","gene")]
  
  cpgs_genes[inBoth_eQTL==FALSE,LinksWeight:=max(LinkScore),by=c("locisID","gene")]
  
  cpgs_genes[inBoth_eQTL==TRUE,LinksWeight:=1]
  
  return(cpgs_genes)
}

CalcCpGScore<-function(res,cpg_genes=NULL){
  if(!is.null(cpg_genes)){
    print("merging limma res and cpgs annotations...")
    res<-data.table(locisID=as.numeric(rownames(res)),
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

#save in