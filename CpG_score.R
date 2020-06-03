
#introducing my Methyl-Gene score
#Goal : we would like to know which Biological pathway/activity is affected by stress-induced epigenetic change.



#We have our list of differentially methylated CpG (DMC) in stressed condition (LGA)
#to make GSEA/overrepresentation we need a list of genes
#So we need to links the CpGs to a genes
#Normally we assume that DMCs are associated with the closest genes (based on its TSS localisation) but it's a too simple assumption

#We are looking to make a score for each DMC reflecting their chance to  be link with a gene expression regulation.

# Based on evidence we define 3 criteria for a DMC to be most probably associated with a gene epigenetic regulation :

# 1) the DMC have a high Fold Change and pvalue => DMC Score
 
# 2) the CpG is on a regulatory region => Regulatory weight 
    #- promoter or enhancer chromatin region define thx to CD34+ ChIP-seq data 
    #- regulatory region in ENSEMBL annotations (CTCF binding site, promoter, re)
    #- on/close to a TF motif

# 3) the Genomic region of the CpG is linked to the gene expression (eQTL study), 
# with a confidence Score of the links => LinksWeight

#2020-06-03 GeneScore v2

options(stringsAsFactors=F)
set.seed(12345)

library(data.table)



#~~ I) link CpGs to genes ~~
# 2 types of CpG-Gene links
# - the closest TSS from a CpG is associated with them
annot<-fread("../../ref/allCpG_annotation_genes_and_feature_regulatory_domain_070520.csv")
cpg.gene1<-annot[!(is.na(gene)|gene=="")&is.na(start.eQTR2)][,.(locisID,gene,distTSS)]
annot<-annot[,.(locisID,chr,pos,type,feature_type_name)]
cpg.gene1
annot

# - the CpG is in a region with eQTL (SNP associated with a change of expression of a gene)
cpg.gene2<-fread("../../ref/2020-06-01_CpG_Gene_links_based_on_whole_blood_eQTL.csv")
cpg.gene2<-cpg.gene2[,in.eQTR:=TRUE][,.(locisID,gene,tss_dist,in.eQTR,avg.mlog10.pv.eQTLs,eQTL_dist)]
cpg.gene2
intersect(cpg.gene1[gene=="SIRT1"]$locisID,cpg.gene2[gene=="SIRT1"]$locisID)
cpg.gene1[locisID==1187983]
cpg.gene2[locisID==1187983]
#homogeneize the 2 dt before merge
cpg.gene1<-cpg.gene1[,tss_dist:=distTSS][,-'distTSS']
cpg.gene2[,tss_dist:=tss_dist+1]
merge(cpg.gene1,cpg.gene2)#only 10k5/236k match cpg-gene => ~5%
unique(cpg.gene2,by="locisID")

cpgs.genes<-merge(cpg.gene1,cpg.gene2,all=T)
cpgs.genes
cpgs.genes[locisID==987]
cpgs.genes[is.na(in.eQTR),in.eQTR:=FALSE]
cpgs.genes

# It is based on blood eQTL get in : https://gtexportal.org/home/datasets 
#file path is : GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt

##IMPROV1 : i plan to integrate also the meta_analysis provide by GTEx of all tissues to get 
# the variant_gene_pairs conserved across tissue (to improve confidence of the link)

#add annot 
cpgs.genes<-merge(cpgs.genes,annot,all.x = T,by="locisID")
#filter for necessary Missing values
cpgs.genes[is.na(type)]#3
cpgs.genes<-cpgs.genes[!is.na(type)]

#Calculate regulatory and links weight : 
# > Regulatory weight [0.5-2] = 0.5 + 1.5 * (TypeScore{0,0.4,0.66,1} + EnsRegScore{0,0.25,0.5,0.75,1})/2 
#     TypeScore{0,0.4,0.66,1} 
cpgs.genes[,TypeScore:=sapply(type,function(x){
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
plot(density(cpgs.genes$TypeScore))

#   EnsRegScore{0,0.25,0.5,0.75,1}
cpgs.genes[,EnsRegScore:=sapply(feature_type_name,function(x){
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
plot(density(cpgs.genes$EnsRegScore))
# Regulatory weight:
cpgs.genes[,RegWeight:=(0.4+1.6*((TypeScore+EnsRegScore)/2))]

plot(density(cpgs.genes[!duplicated(locisID)]$RegWeight))



# > 3) Links Weight[0.5-1] = 0.5+0.5*max(LinkScore [0-1])
#   - for CpG associated to a gene based on TSS proximity  
cpgs.genes[in.eQTR==FALSE,LinkScore:=sapply(abs(tss_dist),function(x){
  if(x<1000){
    score<-1
  }else {
    score<-(3/log10(x))^2
  }
  return(score)
})]
plot(density(cpgs.genes[in.eQTR==FALSE]$LinkScore))
is<-sample(seq(nrow(cpgs.genes[in.eQTR==FALSE])),100000)
plot(cpgs.genes[is,][abs(tss_dist)<20000]$tss_dist,cpgs.genes[is,][abs(tss_dist)<20000]$LinkScore)
# - for CpG associated to a gene based on eQTL studies: 
#de base, linksscore =1

cpgs.genes[in.eQTR==TRUE,LinkScore:=1] 

#but, filter for low signif asso, and add bonus if highly signif eQTR and close to best eQTL
plot(density(na.omit(cpgs.genes$eQTL_dist)))
cpgs.genes[in.eQTR==TRUE,disteQTLScore:=sapply(abs(eQTL_dist),function(x){
  if(x<500){
    distScore<-1 
  }else {
    distScore<-1-0.5*x/2500
  }
  
  return(distScore)
})]
plot(density(na.omit(cpgs.genes$avg.mlog10.pv.eQTLs)))

cpgs.genes[in.eQTR==TRUE,CpG.eQTRScore:=avg.mlog10.pv.eQTLs*disteQTLScore] #integrate with RegScore : mean of pval of eQTL in the region

plot(density(na.omit(cpgs.genes$CpG.eQTRScore))) 
#remove low signif links (<4) and add bonus if 
abline(v=4)
summary(na.omit(cpgs.genes$CpG.eQTRScore))
cpgs.genes<-cpgs.genes[CpG.eQTRScore>4|in.eQTR==F]

#and add bonus ; +CpG.eQTRScore/220
cpgs.genes[in.eQTR==TRUE,LinkScore:=LinkScore+(CpG.eQTRScore/max(CpG.eQTRScore))]
plot(density(cpgs.genes$LinkScore))
cpgs.genes[LinkScore>1.5]
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.001   5.056   7.971  12.019  13.895 221.354 

#there is duplicated rows :
cpgs.genes<-unique(cpgs.genes)
cpgs.genes[,n.cpg.gene.links:=.N,by=c("locisID","gene")]
plot(density(cpgs.genes$n.cpg.gene.links))
cpgs.genes[n.cpg.gene.links>2]
cpgs.genes[n.cpg.gene.links>1]
#...we select the best linkScore by cpg-gene pair
cpgs.genes<-cpgs.genes[,isBestLinks:=LinkScore==max(LinkScore),by=c("locisID","gene")][isBestLinks==TRUE]
cpgs.genes<-unique(cpgs.genes[,-c("isBestLinks","n.cpg.gene.links")])

cpgs.genes[,n.cpg.gene.links:=.N,by=c("locisID","gene")]
plot(density(cpgs.genes$n.cpg.gene.links))
cpgs.genes[n.cpg.gene.links>2][order(gene,locisID,abs(eQTL_dist))]
#there is case of sampe LinkScore,we order by eQTL_dist :

cpgs.genes<-unique(cpgs.genes[order(gene,locisID,abs(eQTL_dist))],by=c("locisID","gene"))
#and finally calculate the linksWeight
cpgs.genes[,LinksWeight:=LinkScore]
plot(density(cpgs.genes$LinksWeight))
summary(cpgs.genes$LinksWeight)

cpgs.genes[round(LinkScore,3)==0.405]

#add utils infos :
cpgs.genes[,nCpG.gene.bef.filtr.:=.N,by=.(gene)]
#Save final cpg-gene links ref
cpgs.genes<-cpgs.genes[,-c("n.cpg.gene.links","disteQTLScore","CpG.eQTRScore","LinkScore")]
fwrite(cpgs.genes,"../../ref/2020-06-03_All_CpG-Gene_links.csv",sep=";")


#~~ II) Calculate The CpG Score ~~
#CpGScore = DMCScore [-inf:inf] * Regulatory weight[0.4-2] * LinksWeight [0.2-2]

##CPGSCORE FOR FEMALE
#differentially methylated datas :
resCpGF<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_FC.FL_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",
               dec = ",")
resCpGF<-resCpGF[,locisID:=V1][,pos:=start-1][,meth.change:=FC][,-c("V1","start")][,.(locisID,chr,pos,pval,meth.change)]
resCpGF
#merge with res
resCpGF.Genes<-merge(resCpGF,cpgs.genes,by=c("locisID","chr","pos"),all.x = T)
resCpGF.Genes

#there is CpG without gene associated, we remove them
sum(is.na(resCpGF.Genes$gene)) #315608
length(unique(resCpGF.Genes$locisID)) # 786277
resCpGF.Genes<-resCpGF.Genes[!is.na(gene)]
resCpGF.Genes

#in eQTL analysis there is CpG linked to several Gene : 
resCpGF.Genes[,nGene.CpG:=.N,by=.(locisID)]
hist(resCpGF.Genes$nGene.CpG,breaks = 50)
resCpGF.Genes[nGene.CpG>10]

resCpGF.Genes

# > DMCScore = meth.change * PvalWeight[0,1]
plot(density(resCpGF.Genes$meth.change))
#PvalWeight[0-1]
plot(density(log10(resCpGF.Genes$pval)))
resCpGF.Genes[,PvalWeight:=(-log10(pval)/3)] #
plot(density(resCpGF.Genes$PvalWeight))
summary(resCpGF.Genes$PvalWeight)
sum(resCpGF.Genes$PvalWeight==0)#146
#multiply the 2 score
resCpGF.Genes[,DMCScore:=PvalWeight*meth.change]
plot(density(resCpGF.Genes$DMCScore))

# > finally :
resCpGF.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]

plot(density(resCpGF.Genes$CpGScore))
resCpGF.Genes[CpGScore>100]
# VALIDATION : 
lines(density(resCpGF.Genes[pval<0.001]$CpGScore),col=2) #majority of significant CpG have a good score
lines(density(resCpGF.Genes[meth.change>20]$CpGScore),col=3) #high meth.Change CpG are pull down 
lines(density(resCpGF.Genes[pval<0.001& meth.change>20]$CpGScore),col=4) # but keep a good score if significant

nrow(resCpGF.Genes[pval>0.01 & CpGScore>40]) # 29 CpG have a high CpGScore wereas pvalue >1%
resCpGF.Genes[pval>0.01 & CpGScore>40][order(-CpGScore)] #because 1) in highly signif eQTR, high meth.change, pval close to 1%, in promoter, and close to gene

plot(density(resCpGF.Genes[pval>0.01 ]$CpGScore)) 
abline(v=quantile(resCpGF.Genes[pval>0.01 ]$CpGScore,0.999)) #99,9% non significant CpG have a CpG score < 103

nrow(resCpGF.Genes[meth.change<20 & CpGScore>40]) # no CpG with small meth.change have a CpGScore >40
plot(density(resCpGF.Genes[abs(meth.change)<20 ]$CpGScore)) 
abline(v=quantile(resCpGF.Genes[abs(meth.change)<20 ]$CpGScore,0.999))#99,9% of small meth.change CpG have a CpG score < 55

# some validation plots : 
# my score highlight mostly CpG with high pvalue and meth.change : 
library(ggplot2)
library(patchwork)
m <- ggplot(resCpGF.Genes, aes( y = CpGScore))

m1<-m +geom_point(aes(x=-log10(pval)))
md1<-m1 + stat_density_2d(aes(x=-log10(pval),fill = after_stat(level)), geom = "polygon")


m2<-m +geom_point(aes(x = meth.change))
md2<-m2 + stat_density_2d(aes(x = meth.change,fill = after_stat(level)), geom = "polygon")
 

m3<-m + geom_point(aes(x = DMCScore))
md3<-m3 + stat_density_2d(aes(x = DMCScore,fill = after_stat(level)), geom = "polygon")

m_all<-(md1 | md2) / 
  md3
m_all


# my score highlight in a second level :
#   - CpGs in regulatory region :
m4<-m + geom_boxplot(aes(x = as.factor(RegWeight)))
m4


#   - CpGs with a strong link with the gene :
m5<-m + geom_point(aes(x = LinksWeight))
md5<-m5 + stat_density_2d(aes(x = LinksWeight,fill = after_stat(level)), geom = "polygon")
md5


##CPGSCORE FOR MALE

resCpGM<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_MC.ML_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",dec = ",")
resCpGM<-resCpGM[,locisID:=V1][,pos:=start-1][,meth.change:=FC][,-c("V1","start")][,.(locisID,chr,pos,pval,meth.change)]
resCpGM
#merge with res
resCpGM.Genes<-merge(resCpGM,cpgs.genes,by=c("locisID","chr","pos"),all.x = T)
resCpGM.Genes

#there is CpG without gene associated, we remove them
sum(is.na(resCpGM.Genes$gene)) #315608
length(unique(resCpGM.Genes$locisID)) # 786277
resCpGM.Genes<-resCpGM.Genes[!is.na(gene)]
resCpGM.Genes

#in eQTL analysis there is CpG linked to several Gene : 
resCpGM.Genes[,nGene.CpG:=.N,by=.(locisID)]
hist(resCpGM.Genes$nGene.CpG,breaks = 50)
resCpGM.Genes[nGene.CpG>10]

resCpGM.Genes

# > DMCScore = meth.change * PvalWeight[0,1]
plot(density(resCpGM.Genes$meth.change))
#PvalWeight[0-1]
plot(density(log10(resCpGM.Genes$pval)))
resCpGM.Genes[,PvalWeight:=(-log10(pval)/3)] #
plot(density(resCpGM.Genes$PvalWeight))
summary(resCpGM.Genes$PvalWeight)
sum(resCpGM.Genes$PvalWeight==0)#146
#multiply the 2 score
resCpGM.Genes[,DMCScore:=PvalWeight*meth.change]
plot(density(resCpGM.Genes$DMCScore))

# > finally :
resCpGM.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]

plot(density(resCpGM.Genes$CpGScore))
resCpGM.Genes[CpGScore>100]



#~~ III) summarize the CpGs Scores to a Gene Score  ~~

##GENESCORE FOR FEMALE
# problem : genes are linked to several CpG...
library(ggplot2)
resCpGF.Genes[,nCpG.Gene:=.N,by=.(gene)]
resCpGF.Genes[,nCpGSig.Gene:=sum(pval<10^-3),by=.(gene)]
plot(density(resCpGF.Genes$CpGScore))
abline(v=20)
resCpGF.Genes[,nHiCpGScore.Gene:=sum(CpGScore>20),gene]

ggplot(resCpGF.Genes)+
  geom_bar(aes(x=nCpG.Gene))

ggplot(resCpGF.Genes)+
  geom_bar(aes(x=nCpGSig.Gene))

ggplot(resCpGF.Genes)+
  geom_bar(aes(x=nHiCpGScore.Gene))

# ...How to summarize the CpG score to a Gene score ??
#to consider :
#- a big gene is most susceptible to have a lot of CpG associated, so more easy for it to have a high CpGscore => distrib most important
#- more a gene have a cpg with a hi CpGScore, more we are confident for the epigentic regulation of the gene
#- we could take just the sum of cpg score link to the gene but some gene have few CpG, so less chance to have high gene score 

#=>  to recapitulate we woulk like a GeneScore which highlight :

# 1) gene with 1/few cpg but with Very Hi CpGScore
# 2) gene with a lot of CpG with hi CpGScore

#for that i decided to weight the sum(CpGScore) in function of the numberofCpG links to the gene with no meth.change.

#GeneScore = sum(CpGScore) * NoChangeCpGWeight
#NoChangeCpGWeight is based on the inverse of the CpgScore :

resCpGF.Genes[,NoChangeCpGScore:=sum(1/(abs(CpGScore)+1)),by="gene"] #bad cpg (~0 ) will have 1/1 
#so if a lot of bad cpg in a gene that make 1+1+1+1+1, so hi nochangeMethScore
plot(density(resCpGF.Genes$NoChangeCpGScore)) # more the number of cpg have, more the score increase
plot(resCpGF.Genes$nCpG.Gene,resCpGF.Genes$NoChangeCpGScore)

#get the inverse 
resCpGF.Genes[,NoChangeCpGScore.inv:=1/NoChangeCpGScore]
plot(density(resCpGF.Genes$NoChangeCpGScore.inv))
plot(resCpGF.Genes$nCpG.Gene,resCpGF.Genes$NoChangeCpGScore.inv)

#and "centralize" with sqrt to avoid too much influence from the extremities

resCpGF.Genes[,NoChangeCpGWeight:=sqrt(NoChangeCpGScore.inv)]
plot(density(resCpGF.Genes$NoChangeCpGWeight))
plot(resCpGF.Genes$nCpG.Gene,resCpGF.Genes$NoChangeCpGWeight) 

resGenesF<-unique(resCpGF.Genes[order(-CpGScore)],by="gene")
summary(resCpGF.Genes[nCpG.Gene==2]$NoChangeCpGWeight)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7071  0.7793  0.8759  1.0198  1.0527  3.9663 
summary(resCpGF.Genes[nCpG.Gene==1]$NoChangeCpGWeight)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.044   1.215   1.580   1.691   6.241 
summary(resCpGF.Genes[nCpG.Gene==7]$NoChangeCpGWeight)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3858  0.4689  0.5328  0.5646  0.6273  1.0891 

#a gene with 1 cpgSig/1cpgtot have a high NoChangeCpGWeight, wereas a gene with 1 cpgSig/10 or 100 cpgtot no
resGenesF[nCpG.Gene==1&nHiCpGScore.Gene==1]$NoChangeCpGWeight
resGenesF[nCpG.Gene==10&nHiCpGScore.Gene==1]$NoChangeCpGWeight #weight ~ 12fois moins

summary(resGenesF$NoChangeCpGWeight)#min =0.04
summary(resGenesF[nCpG.Gene==100]$NoChangeCpGWeight)
summary(resCpGF.Genes[nCpG.Gene>100&nHiCpGScore.Gene>10]$NoChangeCpGWeight)
#to avoid completely exclude of the competition genes with ++cpg, i need to additionate to this wieght a fix coefficient : 
#determine fix coef  (0.1,0.2, 0.3, 0.4 ?); [RERUN with simulated data]
resCpGF.Genes[,GeneScore:=sum(CpGScore)*(0.4+NoChangeCpGWeight),by="gene"]
resGenesF<-unique(resCpGF.Genes[order(-CpGScore)],by="gene")
resGenesF[,ratio.HiCpGScore.Gene:=nHiCpGScore.Gene/nCpG.Gene]
plot(density(resGenesF$GeneScore))
plot(resGenesF$nCpG.Gene,resGenesF$GeneScore)
plot(resGenesF$nHiCpGScore.Gene,resGenesF$GeneScore)
plot(resGenesF$nHiCpGScore.Gene/resGenesF$nCpG.Gene,resGenesF$GeneScore) 
cor(resGenesF$nHiCpGScore.Gene/resGenesF$nCpG.Gene,resGenesF$GeneScore)
plot(resGenesF$NoChangeCpGWeight,resGenesF$GeneScore)

#see gene point en haut a droite
resGenesF[GeneScore>200&nCpG.Gene>300] #gene with numerous hiCpGScore so normal to have a hi GeneScore
resCpGF.Genes[gene=="ZFP57"][order(-CpGScore)]
resCpGF.Genes[locisID==724945]
#intersting gene but too intersting than point en haut a gauche ?
resGenesF[GeneScore>300&nCpG.Gene<10]
#a good normal gene : 16 cpg, with 5 high cpgscore
resGenesF
resGenesF[nCpG.Gene==16&nHiCpGScore.Gene==5]
#LOC644554 have hi Genescore wereas just 1/2 hiCpgScoez 
#after several try, 0.4 is good compromise to put in the same plan small and big gene in term of cpg
resCpGF.Genes[,GeneScore:=sum(CpGScore)*(0.4+NoChangeCpGWeight),by="gene"]
resGenesF<-unique(resCpGF.Genes[order(-CpGScore)],by="gene")

plot(density(resGenesF$GeneScore))
plot(resGenesF[nCpG.Gene<50]$nCpG.Gene,resGenesF[nCpG.Gene<50]$GeneScore)
plot(as.factor(resGenesF[nCpG.Gene<25]$nCpG.Gene),resGenesF[nCpG.Gene<25]$GeneScore)

fwrite(resCpGF.Genes[,-c("GeneScore.v1","NoChangeCpGScore.inv","NoChangeCpGScore")],"analyses/withoutIUGR/2020-06-03_CF.LF_CpG_Gene_scorev3.2.csv",sep=";")


##GENESCORE FOR MALE
resCpGM.Genes[,nCpG.Gene:=.N,by=.(gene)]
resCpGM.Genes[,nCpGSig.Gene:=sum(pval<10^-3),by=.(gene)]
plot(density(resCpGM.Genes$CpGScore))
abline(v=20)
resCpGM.Genes[,nHiCpGScore.Gene:=sum(CpGScore>20),gene]
resCpGM.Genes[CpGScore>100] # cpg of METTL26 intersting : in +++ signif EQTR and close to tss

resCpGM.Genes[,NoChangeCpGScore:=sum(1/(abs(CpGScore)+1)),by="gene"]
resCpGM.Genes[,NoChangeCpGScore.inv:=1/NoChangeCpGScore]
resCpGM.Genes[,NoChangeCpGWeight:=sqrt(NoChangeCpGScore.inv)]

resGenesM<-unique(resCpGM.Genes[order(-CpGScore)],by="gene")
plot(density(resGenesM$NoChangeCpGWeight))
plot(resGenesM$nCpG.Gene,resGenesM$NoChangeCpGWeight) 

resCpGM.Genes[,GeneScore:=sum(CpGScore)*(0.4+NoChangeCpGWeight),by="gene"]

resGenesM<-unique(resCpGM.Genes[order(-CpGScore)],by="gene")

plot(density(resGenesM$GeneScore))
plot(resGenesM$nCpG.Gene,resGenesM$GeneScore)
plot(resGenesM$nHiCpGScore.Gene,resGenesM$GeneScore)
plot(resGenesM$NoChangeCpGWeight,resGenesM$GeneScore)

#see gene point en haut a droite
resGenesM[GeneScore>400&nCpG.Gene>250] #gene with numerous hiCpGScore so normal to have a hi GeneScore
resCpGM.Genes[gene=="LRRC37A2"][order(-CpGScore)] #interesting because lot of cpg in eQTR highly signif
resCpGM.Genes[locisID==1824925] #interesting because links by dist tss at "KANSL1-AS1" (-78 ) and eQTR links allow link to "KANSL1"

plot(as.factor(resGenesM[nCpG.Gene<50]$nCpG.Gene),resGenesM[nCpG.Gene<50]$GeneScore)

fwrite(resCpGM.Genes[,-c("GeneScore.v1","NoChangeCpGScore.inv","NoChangeCpGScore")],"analyses/withoutIUGR/2020-06-03_CM.LM_CpG_Gene_scorev3.2.csv",sep=";")



#~~ OVER-REPRESENTATION TEST ~~
#For Female
library(clusterProfiler)
library(org.Hs.eg.db)
source("scripts/utils.R")
plot(density(resGenesF$GeneScore))
abline(v=170)
genesGS_F<-resGenesF[GeneScore>170]$gene

length(genesGS_F) #1853

rkGS_F <- enrichKEGG(gene         = bitr(genesGS_F,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)


dfkGS_F<-as.data.frame(rkGS_F)

nrow(dfkGS_F) #33
dfkGS_F
dotplot(rkGS_F,showCategory=30)

dfkGS_F<-data.table(dfkGS_F)
dfkGS_F[,GeneScore.avg:=mean(resGenesF$GeneScore[resGenesF$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGS_F
dotplot(rkGS_F,x=dfkGS$GeneScore.avg,showCategory=33)
emapplot(rkGS_F)

enrichplot::upsetplot(rkGS_F,n=25)

#[end 03/06/20 AM updated scripts]

#SEE GENE
library(stringr)
library(clusterProfiler)
library(pathview)
#numeric vector
geneList<-resGenes2$GeneScore3
#named vector
names(geneList)<-resGenes2$gene
#sorted vector
geneList<-sort(geneList,decreasing = T)
genes.df<-bitr(names(geneList),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneList[genes.df$SYMBOL]
geneList.Entrez<-genes.df$GeneScore
names(geneList.Entrez)<-genes.df$ENTREZID

pathw <- pathview(gene.data  = geneList.Entrez,
                  pathway.id = "hsa04930", 
                  species    = "hsa",
                  limit      = list(gene=max(abs(geneList.Entrez)), cpd=1),
                  kegg.dir ="analyses/withoutIUGR/pathwaysMap" 
)
  #for Male
#[CORRECTED] 03/06/20
cpgs.genes.score<-resCpGF.Genes[,-c("pval","meth.change","DMCScore",
                                    "PvalWeight","GeneScore","GeneScore2","GeneScore9","nCpGSig.Gene","nHiCpGScore.Gene")]

resCpGM<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_MC.ML_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",dec = ",")
resCpGM<-resCpGM[,locisID:=V1][,pos:=start-1][,meth.change:=FC][,-c("V1","start")][,.(locisID,chr,pos,pval,meth.change)]
resCpGM
resCpGM.Genes<-merge(resCpGM,cpgs.genes.score,by=c("locisID","chr","pos"))
resCpGM.Genes
resCpGM.Genes[,nCpGSig.Gene:=sum(pval<10^-3),gene]
resCpGM.Genes[,nHiCpGScore.Gene:=sum(CpGScore>10),gene]
#[end Correction] 03/06/20
# > 1) DMCScore = meth.change * PvalWeight[0,1]
plot(density(resCpGM.Genes$meth.change))
#PvalWeight[0-1]
plot(density(log10(resCpGM.Genes$pval)))
resCpGM.Genes[,PvalWeight:=(-log10(pval)/3)] #
plot(density(resCpGM.Genes$PvalWeight))
#multiply the 2 score
resCpGM.Genes[,DMCScore:=PvalWeight*meth.change]
plot(density(resCpGM.Genes$DMCScore))

resCpGM.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]
plot(density(resCpGM.Genes$CpGScore))
#geneSCORE : the sum over th ncpgtot
resCpGM.Genes[,GeneScore:=sum(CpGScore),by="gene"]
resGenes<-unique(resCpGM.Genes,by="gene")[order(-GeneScore)]
resGenes
resGenes$gene[1:100]
plot(density(resGenes$GeneScore))
abline(v=100)
genesGS<-resGenes[abs(GeneScore)>300]$gene
length(genesGS) #203
rkGSM <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

dfkGSM<-as.data.frame(rkGSM)
nrow(dfkGSM) #14
dfkGSM
dotplot(rkGSM,showCategory=30)
dfkGSM<-data.table(dfkGSM)
dfkGSM[,GeneScore.avg:=mean(resGenes2$GeneScore3[resGenes2$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGSM
dotplot(rkGSM,x=dfkGSM$GeneScore.avg,showCategory=38)
emapplot(rkGSM)

#diff with female
setdiff(dfkGS$Description,dfkGSM$Description)
setdiff(dfkGSM$Description,dfkGS$Description)

#with GeneScore9 
library(clusterProfiler)
library(org.Hs.eg.db)
plot(density(resGenes2$GeneScore9))
abline(v=170)
genesGS<-resGenes2[GeneScore9>170]$gene

length(genesGS) #1592

rkGS <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)


dfkGS<-as.data.frame(rkGS)

nrow(dfkGS) #25
dfkGS
dotplot(rkGS,showCategory=30)

dfkGS<-data.table(dfkGS)
dfkGS[,GeneScore.avg:=mean(resGenes2$GeneScore9[resGenes2$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGS
dotplot(rkGS,x=dfkGS$GeneScore.avg,showCategory=25)
emapplot(rkGS)

enrichplot::upsetplot(rkGS,n=25)
#SEE GENE
library(stringr)
library(clusterProfiler)
library(pathview)
#numeric vector
geneList<-resGenes2$GeneScore9
#named vector
names(geneList)<-resGenes2$gene
#sorted vector
geneList<-sort(geneList,decreasing = T)
genes.df<-bitr(names(geneList),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneList[genes.df$SYMBOL]
geneList.Entrez<-genes.df$GeneScore
names(geneList.Entrez)<-genes.df$ENTREZID

pathw <- pathview(gene.data  = geneList.Entrez,
                  pathway.id = "hsa04930", 
                  species    = "hsa",
                  limit      = list(gene=max(abs(geneList.Entrez)), cpd=1),
                  kegg.dir ="analyses/withoutIUGR/pathwaysMap" 
)
#for Male
cpgs.genes.score<-resCpG.Genes[,-c("pval","meth.change","DMCScore","PvalWeight","GeneScore","GeneScore2","GeneScore3")]

resCpGM<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_MC.ML_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",dec = ",")
resCpGM<-resCpGM[,locisID:=V1][,pos:=start-1][,meth.change:=FC][,-c("V1","start")][,.(locisID,chr,pos,pval,meth.change)]
resCpGM
resCpGM.Genes<-merge(resCpGM,cpgs.genes.score,by=c("locisID","chr","pos"))
# > 1) DMCScore = meth.change * PvalWeight[0,1]
plot(density(resCpGM.Genes$meth.change))
#PvalWeight[0-1]
plot(density(log10(resCpGM.Genes$pval)))
resCpGM.Genes[,PvalWeight:=(-log10(pval)/3)] #
plot(density(resCpGM.Genes$PvalWeight))
#multiply the 2 score
resCpGM.Genes[,DMCScore:=PvalWeight*meth.change]
plot(density(resCpGM.Genes$DMCScore))

resCpGM.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]
plot(density(resCpGM.Genes$CpGScore))
#geneSCORE9
resCpGM.Genes[,GeneScore.inv:=sum(1/(abs(CpGScore)+1)),by="gene"]
resCpGM.Genes[,nCpGWeight3:=1/GeneScore.inv]
resCpGM.Genes[,nCpGWeight.norm3.2:=sqrt(nCpGWeight3)]
resCpGM.Genes[,GeneScore9:=sum(CpGScore)*(0.4+nCpGWeight.norm3.2),by="gene"]
resGenes<-unique(resCpGM.Genes,by="gene")
resGenes
resGenes[order(-GeneScore9)]$gene[1:100]

#Male WIth the same cutof
plot(density(resGenes$GeneScore9))
abline(v=170)
genesGS<-resGenes[abs(GeneScore9)>170]$gene
length(genesGS) #214
rkGSM <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05)

dfkGSM<-as.data.frame(rkGSM)

nrow(dfkGSM) #17
dfkGSM
dotplot(rkGSM,showCategory=30)
dfkGSM<-data.table(dfkGSM)
dfkGSM[,GeneScore.avg:=mean(resGenes2$GeneScore3[resGenes2$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGSM
dotplot(rkGSM,x=dfkGSM$GeneScore.avg,showCategory=38)
emapplot(rkGSM)

#diff with female
setdiff(dfkGS$Description,dfkGSM$Description)
# [1] "Transcriptional misregulation in cancer"            
# [2] "Basal cell carcinoma"                               
# [3] "Rap1 signaling pathway"                             
# [4] "Human papillomavirus infection"                     
# [5] "Gastric cancer"                                     
# [6] "Maturity onset diabetes of the young"               
# [7] "MAPK signaling pathway"                             
# [8] "Parathyroid hormone synthesis, secretion and action"
# [9] "Breast cancer"                                      
# [10] "Regulation of actin cytoskeleton"                   
# [11] "Proteoglycans in cancer"                            
# [12] "Notch signaling pathway"                            
# [13] "TGF-beta signaling pathway"                         
# [14] "Prostate cancer"                                    
# [15] "Hedgehog signaling pathway"                         
# [16] "Thyroid cancer"                                     
# [17] "Melanogenesis"                                      
# [18] "Phospholipase D signaling pathway"                  
# [19] "Cortisol synthesis and secretion"    
setdiff(dfkGSM$Description,dfkGS$Description)
# [1] "Wnt signaling pathway"                               
# [2] "Neurotrophin signaling pathway"                      
# [3] "AGE-RAGE signaling pathway in diabetic complications"
# [4] "HIF-1 signaling pathway"                             
# [5] "Oxytocin signaling pathway"                          
# [6] "Chronic myeloid leukemia"                            
# [7] "Cellular senescence"                                 
# [8] "Yersinia infection"                                  
# [9] "B cell receptor signaling pathway"                   
# [10] "C-type lectin receptor signaling pathway"            
# [11] "GnRH secretion"           

#male with the same nb of gene than female (1592)

genesGS<-resGenes[order(-GeneScore9)]$gene[1:1592]
length(genesGS) #1592
plot(density(resGenes$GeneScore9))
abline(v=resGenes[order(-GeneScore9)]$GeneScore9[1592])#88
rkGSM <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05)

dfkGSM<-as.data.frame(rkGSM)

nrow(dfkGSM) #55
dfkGSM
dotplot(rkGSM,showCategory=30)
dfkGSM<-data.table(dfkGSM)
dfkGSM[,GeneScore.avg:=mean(resGenes$GeneScore9[resGenes$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGSM
dotplot(rkGSM,x=dfkGSM$GeneScore.avg,showCategory=55)
emapplot(rkGSM,showCategory=55)

#diff with female
setdiff(dfkGS$Description,dfkGSM$Description)
# [1] "Parathyroid hormone synthesis, secretion and action"
# [2] "Regulation of actin cytoskeleton"                   
# [3] "Notch signaling pathway"                            
# [4] "Hedgehog signaling pathway"                         
# [5] "Thyroid cancer"                                     
# [6] "Cortisol synthesis and secretion"    


setdiff(dfkGSM$Description,dfkGS$Description)

# [1] "Wnt signaling pathway"                               
# [2] "mTOR signaling pathway"                              
# [3] "Acute myeloid leukemia"                              
# [4] "Spinocerebellar ataxia"                              
# [5] "AMPK signaling pathway"                              
# [6] "Non-small cell lung cancer"                          
# [7] "Neurotrophin signaling pathway"                      
# [8] "Choline metabolism in cancer"                        
# [9] "Endometrial cancer"                                  
# [10] "Glioma"                                              
# [11] "Chronic myeloid leukemia"                            
# [12] "Fc gamma R-mediated phagocytosis"                    
# [13] "Circadian rhythm"                                    
# [14] "Central carbon metabolism in cancer"                 
# [15] "Colorectal cancer"                                   
# [16] "Cellular senescence"                                 
# [17] "Phosphatidylinositol signaling system"               
# [18] "Growth hormone synthesis, secretion and action"      
# [19] "Fc epsilon RI signaling pathway"                     
# [20] "Platelet activation"                                 
# [21] "Prolactin signaling pathway"                         
# [22] "EGFR tyrosine kinase inhibitor resistance"           
# [23] "Longevity regulating pathway - multiple species"     
# [24] "cAMP signaling pathway"                              
# [25] "Endocrine resistance"                                
# [26] "Adherens junction"                                   
# [27] "Longevity regulating pathway"                        
# [28] "Autophagy - animal"                                  
# [29] "Alzheimer disease"                                   
# [30] "HIF-1 signaling pathway"                             
# [31] "AGE-RAGE signaling pathway in diabetic complications"
# [32] "Bacterial invasion of epithelial cells"              
# [33] "Insulin signaling pathway"                           
# [34] "Adrenergic signaling in cardiomyocytes"              
# [35] "Ras signaling pathway"                               
# [36] "Thyroid hormone signaling pathway"   

#compa male female
candidat_genes.list<-list(female=bitr(resGenes2[GeneScore9>170]$gene,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                          male=bitr(resGenes[GeneScore9>170]$gene,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID)
res.compa<-compareCluster(candidat_genes.list,
                          fun = "enrichKEGG",
                          organism="hsa",
                          pvalueCutoff=0.05) 
dotplot(res.compa,showCategory =30,size='count')

emapplot(res.compa,pie="count",showCategory =30 )
#save dts
fwrite(resCpGF.Genes[,GeneScore.v1:=GeneScore]
       [,GeneScore.v2:=GeneScore3]
       [,GeneScore.v3:=GeneScore9]
       [,.(locisID,chr,pos,pval,PvalWeight,meth.change,DMCScore,nGene.CpG,
                gene,tss_dist,in.eQTR,avg.mlog10.pv.eQTLs,eQTL_dist,LinkScore,LinksWeight,
                type,TypeScore,feature_type_name,EnsRegScore,RegWeight,CpGScore,
                nCpG.Gene,nCpGSig.Gene,nHiCpGScore.Gene,GeneScore.v1,GeneScore.v2,GeneScore.v3)],
       file = "analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev3.csv",
       sep = ";")

resCpGM.Genes[,coef:=0.5+0.5*sqrt(1/nCpG.Gene)]

resCpGM.Genes[,GeneScore:=coef*CpGScore[which.max(abs(CpGScore))]+(1-coef)*median(CpGScore),by="gene"]


resCpGM.Genes[,GeneScore3:=sum(CpGScore),by="gene"]

fwrite(resCpGM.Genes[,GeneScore.v1:=GeneScore]
       [,GeneScore.v2:=GeneScore3]
       [,GeneScore.v3:=GeneScore9]
       [,.(locisID,chr,pos,pval,PvalWeight,meth.change,DMCScore,nGene.CpG,
           gene,tss_dist,in.eQTR,avg.mlog10.pv.eQTLs,eQTL_dist,LinkScore,LinksWeight,
           type,TypeScore,feature_type_name,EnsRegScore,RegWeight,CpGScore,
           ncpg.gene,nCpGSig.Gene,,nHiCpGScore.Gene,GeneScore.v1,GeneScore.v2,GeneScore.v3)],
       file = "analyses/withoutIUGR/2020-06-02_CM.LM_CpG_Gene_scorev3.csv",
       sep = ";")

#2020-06-02 : search pathways of interest
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
source("scripts/utils.R")

resCpG.Genes<-fread("analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev3.csv")
resCpGM.Genes<-fread("analyses/withoutIUGR/2020-06-02_CM.LM_CpG_Gene_scorev3.csv")
resCpG.Genes[,ncpg.gene:=.N,by="gene"]
resCpGM.Genes[,ncpg.gene:=.N,by="gene"]
fwrite(resCpG.Genes,"analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev3.csv",sep=";")
fwrite(resCpGM.Genes,"analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev3.csv",sep=";")
#see hippo and wnt pathway in female/male
#female
resGenesF<-unique(resCpG.Genes[order(pval)],by="gene")
genesGSF<-resGenesF[GeneScore.v3>170]$gene
length(genesGSF) #1592
plot(density(resGenesM$GeneScore.v3))
rkGSF <- enrichKEGG(gene         = bitr(genesGSF,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05)
dfkGSF<-as.data.frame(rkGSF)

nrow(dfkGSF) #25
dfkGSF
dotplot(rkGSF,showCategory=30)

dfkGSF<-data.table(dfkGSF)
dfkGSF[,GeneScore.avg:=mean(resGenesF$GeneScore.v3[resGenesF$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGSF
dotplot(rkGSF,x=dfkGSF$GeneScore.avg,showCategory=25)
emapplot(rkGSF)

enrichplot::upsetplot(rkGS,n=25)

resGenesF[gene%in%c("SIRT1","TP53","FOXO3","SIRT3", "SOD2", "IDH1","ACACA","ACACB","MTOR",
                    "RPTOR","MLST8","PRAS40","DEPTOR",
                    "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG1", "PRKAG2", "PRKAG3")][order(-GeneScore.v3)]

#male
resGenesM<-unique(resCpGM.Genes[order(pval)],by="gene")
genesGSM<-resGenesM[GeneScore.v3>170]$gene
length(genesGS) #214
plot(density(resGenesM$GeneScore.v3))
abline(v=resGenesM[order(-GeneScore9)]$GeneScore9[1592])#88
rkGSM <- enrichKEGG(gene         = bitr(genesGSM,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05)


#SEE GENE
library(stringr)
library(clusterProfiler)
library(pathview)
#numeric vector
geneList<-resGenes2$GeneScore9
#named vector
names(geneList)<-resGenes2$gene
#sorted vector
geneList<-sort(geneList,decreasing = T)
genes.df<-bitr(names(geneList),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneList[genes.df$SYMBOL]
geneList.Entrez<-genes.df$GeneScore
names(geneList.Entrez)<-genes.df$ENTREZID

pathw <- pathview(gene.data  = geneList.Entrez,
                  pathway.id = "hsa04930", 
                  species    = "hsa",
                  limit      = list(gene=max(abs(geneList.Entrez)), cpd=1),
                  kegg.dir ="analyses/withoutIUGR/pathwaysMap" 
)

#2020-06-02 : 1) correct some details, 2) sea if i can improve the ncpgweight
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
source("scripts/utils.R")
#☻1) correct some details
#cpgsig count :
resCpGF.Genes<-fread("analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev3.csv")
resCpGF.Genes
resCpGF.Genes[,nCpGSig.Gene:=sum(pval<10^-3),gene]
resCpGF.Genes[,nHiCpGScore.Gene:=sum(CpGScore>10),gene]
resCpGF.Genes[gene=="MTOR"]
resCpGF.Genes[gene=="FOXO3"][order(pval)]
#add true cpgsig count for male
resCpGM.Genes<-fread("analyses/withoutIUGR/2020-06-02_CM.LM_CpG_Gene_scorev3.csv")
resCpGM.Genes[,nCpGSig.Gene:=sum(pval<10^-3),gene]
resCpGM.Genes[gene=="MTOR"]
resCpGM.Genes[gene=="FOXO3"][order(pval)]
resCpGM.Genes[,nHiCpGScore.Gene:=sum(CpGScore>10),gene]


#i correct in the script plus haut ([corrected])
#save :
fwrite(resCpGF.Genes,"analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev3.csv",sep=";")
fwrite(resCpGM.Genes,"analyses/withoutIUGR/2020-06-02_CM.LM_CpG_Gene_scorev3.csv",sep=";")

#IMPROV  : test the reliability and robustness of my score;
#ideas : 
# - change number of genes / geneScore threshold => Are the pathway enrichment robust ? or highly variable
# - test with simulated data (random permutation..)
# - test with disease signature methylation data
