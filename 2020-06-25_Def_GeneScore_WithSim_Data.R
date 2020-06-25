
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
sim_cpg[CpGScore>30&distTSS>99999]


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

#make 10
#gene with 3