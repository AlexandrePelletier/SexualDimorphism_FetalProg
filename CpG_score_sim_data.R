

#see behaviour of my CpGScore/GeneScore with simulated
options(stringsAsFactors=F)
set.seed(12345)
library(data.table)
# I) test my CpG Score 
#   1) create Simulated CpG
fc <- c(5,10,15,20,25,30,40,50,70 )
fc<-sort(c(-fc,fc,0))
fc
pval <- c(0.5,0.01,0.001,0.0001,10^-5,5*10^-6,10^-6,5*10^-7,10^-7,10^-10)

type <- 0:6
ensReg1 <- c("CTCF Binding Site","Promoter", "Enhancer","Open chromatin","Promoter Flanking Region")


ensReg2<-c("",ensReg1,apply(expand.grid(ensReg1,c("TF binding site")),1,paste,collapse="/"),"TF binding site")
ensReg2
eQTLlinked<-c(FALSE,TRUE)
distTSS<-c(0,50,100,200,1000,2000,5000,10000,20000,50000,100000,200000,10^6,10^7,10^8)
eQTRWidth<-c(1,100,1000, 2000,5000,10000,20000,50000,100000)

sim_cpg <- data.table(expand.grid("meth.change"=fc,"pval"= pval,"type"= type,
                                  "feature_type_name"=ensReg2,"in_eQTR"=eQTLlinked,
                                  "distTSS"=distTSS,"eQTR_width"=eQTRWidth))
sim_cpg #4309200 fake CpG

#   2) calculate CpG Score
sim_cpg[,PvalWeight:=(-log10(pval))]

sim_cpg[,DMCScore:=PvalWeight*meth.change]


sim_cpg[,TypeScore:=sapply(type,function(x){
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

source("scripts/utils.R")
sim_cpg[,EnsRegScore:=sapply(as.character(feature_type_name),function(x){
  vecX<-strsplit(x,"/")[[1]]
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



sim_cpg[,RegWeight:=(0.5+1.5*((TypeScore+EnsRegScore)/2))]
sort(unique(sim_cpg$RegWeight))


sim_cpg[in_eQTR==F,LinkScore:=sapply(abs(distTSS),function(x){
  if(x<1000){
    score<-1
  }else if(x<20000){
    score<-3/(log10(x))
  }else if (x<100000){
    score<-(3/log10(x))^2
    
  }else{
    score<-(3/log10(x))^3
  }
  return(score)
})]
summary(sim_cpg$LinkScore)

sim_cpg[in_eQTR==T,LinkScore:=sapply(eQTR_width,function(width){
  if(width<5000){
    score<-1+0.2*log10(2)/log10(width+1)
  }else {
    score<-log10(5000)/log10(width+1)
  }
  return(score)
})]

sim_cpg[,LinksWeight:=0.1+0.9*LinkScore]
summary(sim_cpg$LinksWeight)

ggplot(unique(sim_cpg[in_eQTR==F],by=c('distTSS')),aes(log10(distTSS),LinksWeight))+geom_point()
sim_cpg[in_eQTR==F&distTSS>20000]

sim_cpg[,CpGScore:=DMCScore*RegWeight*LinksWeight]

sim_cpg


#   3) plot components to see its behaviour
#each component influence the score ? (use facet_grid,)
summary(sim_cpg$CpGScore)
library(ggplot2)

#1 par 1 
#pval
ggplot(unique(sim_cpg,by=c("pval","CpGScore")),aes(x=factor(round(-log10(pval),1)),y=CpGScore))+geom_boxplot()
ggplot(unique(sim_cpg[meth.change>20],by=c("pval","CpGScore")),aes(x=factor(round(-log10(pval),1)),y=CpGScore))+geom_boxplot()


#meth change
ggplot(unique(sim_cpg,by=c("meth.change","CpGScore")),aes(x=factor(meth.change),y=CpGScore))+geom_boxplot()
ggplot(unique(sim_cpg[pval<10^-3],by=c("meth.change","CpGScore")),aes(x=factor(meth.change),y=CpGScore))+geom_boxplot()

#typescore
ggplot(unique(sim_cpg[meth.change>0],by=c("type","CpGScore")),aes(x=factor(type),y=CpGScore))+geom_boxplot()
ggplot(unique(sim_cpg[meth.change>20 & pval<10^-3 ],by=c("type","CpGScore")),aes(x=factor(type),y=CpGScore))+geom_boxplot()

#ensembl_regscore
ggplot(unique(sim_cpg[meth.change>0],by=c("feature_type_name","CpGScore")),aes(x=factor(feature_type_name),y=CpGScore))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(unique(sim_cpg[meth.change>20 & pval<10^-3],by=c("feature_type_name","CpGScore")),aes(x=factor(feature_type_name),y=CpGScore))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#linkscore
ggplot(unique(sim_cpg[in_eQTR==F & meth.change>20 & pval<10^-3],by=c("distTSS","CpGScore")),aes(x=factor(distTSS),y=CpGScore))+geom_boxplot()

ggplot(unique(sim_cpg[in_eQTR==T & meth.change>20 & pval<10^-3],by=c("eQTR_width","CpGScore")),aes(x=factor(eQTR_width),y=CpGScore))+geom_boxplot()


#etude de cas :
# 1 cpg "highMet" mais  regul = 1 et distTSS = 10000
# 1 cpg "regulmax" mais methchange = 20 et pval = 10^-4 et distTSS = 10000
#1 cpg "linksMax" avec distTSS =1  mais methchange = 20 et pval = 5.10^-4 et regul = 1

#quelle est le classement pour highmet avec FC =40 et pval =10^-6 ?

sim_cpg[in_eQTR==F&meth.change==40 & pval==10^-6 & EnsRegScore==0.5 & TypeScore==0.5 & distTSS==10000]$CpGScore
#232.5
sim_cpg[in_eQTR==F&meth.change==20 & pval==10^-4 & EnsRegScore==1 & TypeScore==1 & distTSS==10000]$CpGScore
#124

sim_cpg[in_eQTR==F&meth.change==20 & pval==10^-4 & EnsRegScore==0.5 & TypeScore==0.5 & distTSS==0]$CpGScore
#100

#c'est bien que regulmax > linksMax par contre highMet trop important ?
sim_cpg[in_eQTR==F&meth.change==25 & pval==10^-6 & EnsRegScore==0 & TypeScore==0.5 & distTSS==10000]$CpGScore
#meme a fc=25 avec cette pvalue on est dessus de regulmax pval^-4 => trop d'importance à la pval ?
ggplot(sim_cpg[in_eQTR==F&meth.change==25 & pval<10^-4 & EnsRegScore==0.5 & TypeScore==0.5 & distTSS==10000],aes(log10(pval),CpGScore))+geom_point()


#compox * compo y
ggplot(unique(sim_cpg[meth.change>5],by=c("pval","CpGScore")),aes(x=factor(round(-log10(pval),1)),y=CpGScore))+geom_boxplot()+facet_wrap("meth.change")



#   4)adjust weight of each component

sim_cpg[in_eQTR==F,LinkScore:=sapply(abs(distTSS),function(x){
  if(x<1000){
    score<-1
  }else if(x<20000){
    score<-3/(log10(x))
  }else if (x<100000){
    score<-(3/log10(x))^2
    
  }else{
    score<-(3/log10(x))^3
  }
  return(score)
})]
summary(sim_cpg$LinkScore)

sim_cpg[in_eQTR==T,LinkScore:=sapply(eQTR_width,function(width){
  #! a rechanger en regardant prox methyl with highest locus/locis signif
  # ++ meQTL et metaeQTL
  if(width<5000){
    score<-1+0.2*log10(2)/log10(width+1)
  }else {
    score<-log10(5000)/log10(width+1)
  }
  return(score)
})]

sim_cpg[,LinksWeight:=0.1+0.9*LinkScore]

#   5) search validated-CpG

#   6) calculate CpG Score for it

#   7) plot with simulated
#it is in the top CpG ?


# II) test my Gene Score 
#    1) create simulated 1000 Gene * 1:500 nCpG
n_cpg <- 2
gene_data <- simul_data[sample(seq(nrow(simul_data)), size = n_cpg), ]

#   2) calculate GeneScore

#   3) plot components to see its behaviour

#   4) adjust weight of each component

#   5) search validated multiple CpG-genes

#   6) calculate Gene Score for it

#   7) plot with simulated 



# III) compare with DMRcat (par exemple)
#   Find DMR 
# annotate DMR, with ENSEMBL












