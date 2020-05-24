
#introducing my CpG score
#Goal : we would like to know which Biological pathway/activity is affected by stress-induced epigenetic change.

#We have our list of differentially methylated CpG (DMC) in stressed condition (LGA)
#to make GSEA/overrepresentation we need a list of genes
#So we need to links the CpGs to a genes
#Normally we assume that DMCs are associated with the closest genes (based on its TSS localisation) but it's a too simple assumption

#We are looking to make a score allowing to sort DMC based on their chance to really impact the gene expression.

# Based on evidence we define 3 criteria for a DMC to be most probably associated with a gene epigenetic regulation :

# 1) the DMC have a high Fold Change and pvalue => DMC Score
 

# 2) the CpG is on a regulatory region => Regulatory weight 
    #- promoter or enhancer chromatin region define thx to CD34+ ChIP-seq data 
    #- regulatory region in ENSEMBL annotations (CTCF binding site, promoter, re)
    #- on/close to a TF motif

# 3) the Genomic region of the CpG is linked to the gene expression (eQTL study), 
# with a confidence Score of the links => LinksWeight

options(stringsAsFactors=F)
set.seed(12345)

library(data.table)


#differentially methylated datas :
resCpG<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_FC.FL_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",
              dec = ",")
resCpG<-resCpG[,locisID:=V1][,pos:=start-1][,-c("V1","start")][,.(locisID,chr,pos,pval,FC)]
resCpG


#~~ I) link CpGs to genes ~~

annot<-fread("../../ref/allCpG_annotation_genes_and_feature_regulatory_domain_070520.csv")
annot<-annot[!(is.na(gene)|gene=="")]
annot
# 2 types of CpG-Gene links
# - the closest TSS from a CpG is associated with them
annot[is.na(start.eQTR2)]
max(annot[is.na(start.eQTR2)]$distTSS)

# - the CpG is in a region with eQTL (SNP associated with a change of expression of a gene)
resCpG.Genes<-merge(resCpG,annot,by=c("locisID","chr","pos"),all.x = T)
resCpG.Genes

#there is CpG without gene associated, we remove them
sum(is.na(resCpG.Genes$gene)) #339390
length(unique(resCpG.Genes$locisID)) # 786277
resCpG.Genes<-resCpG.Genes[!is.na(gene)]
resCpG.Genes

#in eQTL analysis there is CpG with several CpG, to avoid noise we select just 1 gene by CpG : the closest gene
resCpG.Genes[is.na(start.eQTR2),closest.Genes:=T,by=.(locisID)]
resCpG.Genes[!is.na(start.eQTR2),closest.Genes:=distTSS==min(distTSS),by=.(locisID)]
resCpG.Genes<-resCpG.Genes[closest.Genes==T]
resCpG.Genes

#~~ II) Calculate The CpG Score ~~
#CpGScore = DMCScore [-inf:inf] * Regulatory weight[0.5-2] * LinksWeight [0.5-1]

# > 1) DMCScore = FC * PvalWeight[0.1,1]
plot(density(resCpG.Genes$FC))
#PvalWeight[0-1]
plot(density(log10(resCpG.Genes$pval)))
resCpG.Genes[,PvalWeight:=(-log10(pval)/8)]
plot(density(resCpG.Genes$PvalWeight))

#multiply the 2 score
resCpG.Genes[,DMCScore:=PvalWeight*FC]
plot(density(resCpG.Genes$DMCScore))

# > 2) Regulatory weight [0.5-2] = 0.5 + 1.5 * (TypeScore{0,0.4,0.66,1} + EnsRegScore{0,0.25,0.5,0.75,1})/2 
#     TypeScore{0,0.4,0.66,1} 
resCpG.Genes[,TypeScore:=sapply(type,function(x){
  if(x%in%c(4,6)){
    score<-1
  }else if(x==5){
    score<-0.66
  }else if(x%in%1:3){
    score<-0.4
  }else{
    score<-0
  }
  return(score)
})]
plot(density(resCpG.Genes$TypeScore))

#   EnsRegScore{0,0.25,0.5,0.75,1}
source("scripts/utils.R")
resCpG.Genes[,EnsRegScore:=sapply(feature_type_name,function(x){
  vecX<-tr(x)
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
plot(density(resCpG.Genes$EnsRegScore))
# Regulatory weight:
resCpG.Genes[,RegWeight:=(0.5+1.5*((TypeScore+EnsRegScore)/2))]
plot(density(resCpG.Genes[!duplicated(locisID)]$RegWeight))



# > 3) Links Weight[0.5-1] = 0.5 + 0.5 * (TypeScore{0,0.4,0.66,1} + EnsRegScore{0,0.25,0.5,0.75,1}) /2
#   - for CpG associated to a gene based on TSS proximity  
resCpG.Genes[is.na(start.eQTR2),LinkScore:=sapply(abs(distTSS),function(x){
  if(x<1000){
    score<-1
  }else if(x<20000){
    score<-(3/(log10(x)))^0.5
  }else if (x<100000){
    score<-3/log10(x)
    
  }else{
    score<-(3/log10(x))^2
  }
  return(score)
})]
plot(density(resCpG.Genes[is.na(start.eQTR2)]$LinkScore))

# - for CpG associated to a gene based on eQTL studies: 
#   *need eQTR : most precise eQTL region
eQTRs<-fread("../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv")
eQTRs<-eQTRs[,.(chr,start.eQTR2,end.eQTR2,start.eQTR,end.eQTR,gene)]
eQTRs
resCpG.Genes<-merge(resCpG.Genes,eQTRs,all.x = TRUE,by=c("chr","start.eQTR2", "end.eQTR2","gene"))[order(pval)]
resCpG.Genes

#   * I made a function to calculate confidence of the links CpG-Gene for eQTL linked Gene : 
calcLinkScore<-function(pos,start,end,signif=NULL,maxSignif=0.001475078){
  
  # i) linkScore ~ CpG distance of the eQTL region (eQTR)
  if(pos > (start-50)& pos < (end+50)){
    score<-1
  }else {
    x<-min(abs(c(start-pos,end-pos)))
    score<-(log10(50)/(log10(x)))^0.5
  }
  
  if(!is.null(signif)){
    print('test')
    #si min pval >*1 si max pval > *0.75
    score<-score*(0.75+0.25*(signif/maxSignif))
  }
  
  
  # ii) ajust score if eQTL region too large (because too much uncertainty) 
  #    si >200, ajuste score sinon descend progressivement
  largeur<-end-start+1
  if(largeur<200){
    score<-score*1
  }else if(largeur<10000){
    score<-score*(0.5+0.5*(log10(200)/log10(largeur+1))^0.5)
  }else if(largeur<100000){
    score<-score*(0.5+0.5*(log10(200)/log10(largeur+1)))
  }else{
    score<-score*(0.5+0.5*(log10(200)/log10(largeur+1))^2)
  }
  
  
  return(score)
  
}
  #   *calculate linkScore for eQTL-based CpG-Gene pairs :
resCpG.Genes[!is.na(start.eQTR2),LinkScore:=calcLinkScore(pos[1],start.eQTR[1],end.eQTR[1]),by=.(locisID,gene,start.eQTR2)]

lines(density(resCpG.Genes[!is.na(start.eQTR2)]$LinkScore),col=2)

resCpG.Genes[,LinksWeight:=0.5+0.5*max(LinkScore),by=.(locisID,gene)]

#because eQTR can overlap :
resCpG.Genes<-resCpG.Genes[,nEQTR.CpgGeneLink:=.N,by=c("locisID","gene")]
plot(density(resCpG.Genes$nEQTR.CpgGeneLink))

#we select the best linkScore by cpg-gene pair
resCpG.Genes<-resCpG.Genes[,isBestLinks:=LinkScore==max(LinkScore),by=c("locisID","gene")][isBestLinks==TRUE]
resCpG.Genes<-unique(resCpG.Genes)
plot(density(resCpG.Genes$LinksWeight))

#Note : i would like to be more precise in the calculation of this links score


# > finally :
resCpG.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]
plot(density(resCpG.Genes$CpGScore))

lines(density(resCpG.Genes[pval<0.001]$CpGScore),col=2)
lines(density(resCpG.Genes[FC>20]$CpGScore),col=3)
lines(density(resCpG.Genes[pval<0.001& FC>20]$CpGScore),col=4)

nrow(resCpG.Genes[pval>0.05 & CpGScore>10])
plot(density(resCpG.Genes[pval>0.01 ]$CpGScore))

nrow(resCpG.Genes[FC<15 & CpGScore>10])
plot(density(resCpG.Genes[abs(FC)<20 ]$CpGScore))
resCpG.Genes[abs(FC)<20 & CpGScore>10]


#~~ III) summarize the CpGs Scores to a Gene Score  ~~
# problem : genes are linked to several CpG...
library(ggplot2)
resCpG.Genes[,nCpG.Gene:=.N,by=.(gene)]
resCpG.Genes[CpGScore>10,nCpGSig.Gene:=.N,by=.(gene)]

resCpG.Genes
ggplot(resCpG.Genes)+
  geom_bar(aes(x=nCpG.Gene))

ggplot(resCpG.Genes[CpGScore>10])+
  geom_bar(aes(x=nCpGSig.Gene))


# ...How to summarize the CpG score to a gene score ??
#to consider :
#- a big gene is most susceptible to have a lot of CpG, so more easy for it to have a high CpGscore => distrib more important
#- some CpGs in the same regulatory region, they are markers for the same regulation

#=>  the GeneScore need to integrate :
# 1) the CpGScore max 
# 2) the number of CpG associated with the Gene 
# 3) the distribution of their CpGscores associated with it,

#GeneScore[-inf:inf] = c*max(CpG) + (1-c)*median(CpGScore)
#with c[0.5-1] ~ nCpG associated to the gene. 
# for a gene with a lot of CpG, the CpG score distribution is more important => c=0.5, the maximum and the median have the same weight
# for a gene with 2 CpG => c ~= 1, the max(CpGScore) is more important
#  c= 0.5+0.5*sqrt(1/nCpG)
resCpG.Genes[,coef:=0.5+0.5*sqrt(1/nCpG.Gene)]


resCpG.Genes[,GeneScore:=coef*extr(CpGScore)[1]+(1-coef)*median(CpGScore),by="gene"]

resGenes<-unique(resCpG.Genes,by="gene")[order(-GeneScore)]
resGenes
resGenes$gene[1:100]
plot(density(resGenes$nCpG.Gene))
plot(density(resGenes$coef))
plot(density(resGenes$GeneScore))
abline(v=15)

#more precise : GeneScore2[-inf:inf]= moy ponder ci*CpgScore/sum(ci) avec ci~1/nCpG

moyPond<-function(CpGScores){
  CpGScores<-CpGScores[order(abs(CpGScores),decreasing = T)]
  coefs<-sapply(1:length(CpGScores), function(x)1/(x))
  
  return(sum(CpGScores*coefs)/sqrt(sum(coefs)))
}

resCpG.Genes[,GeneScore2:=moyPond(CpGScore),by="gene"]

resGenes2<-unique(resCpG.Genes,by="gene")[order(-GeneScore2)]

resGenes2
resGenes2$gene[1:100]

plot(density(resGenes2$GeneScore2))
abline(v=20)

#validation : 

plot(-log10(resGenes$pval),resGenes$GeneScore2)
abline(h=2)
plot(resGenes$FC,resGenes$GeneScore2)
abline(h=10)
plot(resGenes$DMC,resGenes$GeneScore2)
plot(resGenes$LinksWeight,resGenes$GeneScore)
plot(as.factor(resGenes$RegWeight),resGenes$GeneScore)


#idea to improve the summarize : summarize CpGScore by regulatory region and sum this scores to have a GeneScore


library(clusterProfiler)
library(org.Hs.eg.db)

#OVER-REPRESENTATION TEST
#genes candidat : GeneScore >10
plot(density(resGenes$GeneScore))
abline(v=15)
genes15<-resGenes[GeneScore>15]$gene
length(genes15) #1700

rk15 <- enrichKEGG(gene         = bitr(genes15,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1)

dfk15<-as.data.frame(rk15)
nrow(dfk15)
dfk15
dotplot(rk15,showCategory=30)

source("scripts/utils.R")
library(data.table)
dfk15<-data.table(dfk15)
dfk15[,GeneScore.avg:=mean(resGenes$GeneScore[resGenes$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfk15


dotplot(rk15,x=dfk15$GeneScore.avg)

saveRDS(rk15,"analyses/withoutIUGR/2020-05-24_CF.LF_clusterProf_object_GSEA_GeneScore1.rds")
write.csv2(as.data.frame(setReadable(rk15,org.Hs.eg.db,keyType = "ENTREZID")),"analyses/withoutIUGR/2020-05-24_CF.LF_res_GSEA_GeneScore1.csv")

#GeneScore2



#pvalue


#comparaison 





#cnetplot 
# - make the gene list (for foldchange) :
#numeric vector
geneList<-resGenes$GeneScore
#named vector
names(geneList)<-resGenes$gene
#sorted vector
geneList<-sort(geneList,decreasing = T)
genes.df<-bitr(names(geneList),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneList[genes.df$SYMBOL]
geneList.Entrez<-genes.df$GeneScore
names(geneList.Entrez)<-genes.df$ENTREZID

cnetplot(setReadable(rk15,keyType = "ENTREZID",OrgDb = org.Hs.eg.db),foldChange =geneList ,readable=T)




## ANALYSE MALE 
#compare with FC of the min pvalue locis by gene
resCpG.Genes[,is.min.pval:=pval==min(pval),by=gene]
resGenes2<-resCpG.Genes[is.min.pval==TRUE,]
resGenes2<-unique(resGenes2,by='gene')
geneList2<-resGenes2$FC
names(geneList2)<-resGenes2$gene
geneList2<-sort(geneList2,decreasing = T)
genes.df<-bitr(names(geneList2),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$FC<-geneList2[genes.df$SYMBOL]
geneList.Entrez<-genes.df$FC
names(geneList.Entrez)<-genes.df$ENTREZID
resKEGG2<- gseKEGG(geneList     = geneList.Entrez,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05,
                  verbose = FALSE)

#compare 
nrow(as.data.frame(resKEGG)) #36
nrow(as.data.frame(resKEGG2)) #82
setdiff(as.data.frame(resKEGG)$Description,as.data.frame(resKEGG2)$Description)
setdiff(as.data.frame(resKEGG2)$Description,as.data.frame(resKEGG)$Description)

library(patchwork)
p1<-dotplot(resKEGG, showCategory=30)+ggtitle("with my CpG Score")

p2<-dotplot(resKEGG2, showCategory=30)+ggtitle("with the Fold Change")

p1+p2

emapplot(resKEGG2)




#les males
#differentially methylated datas :
resCpGM<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_MC.ML_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",
              dec = ",")
resCpGM<-resCpGM[,locisID:=V1][,pos:=start-1][,-c("V1","start")][,.(locisID,chr,pos,pval,FC)]
resCpGM


#~~ I) link CpGs to genes ~~

annot<-fread("../../ref/allCpG_annotation_genes_and_feature_regulatory_domain_070520.csv")

annot<-annot[!(is.na(gene)|gene=="")]
annot # 2 types of links
annot[!is.na(start.eQTR2)]
resCpGM.Genes<-merge(resCpGM,annot,by=c("locisID","chr","pos"),all.x = T)
resCpGM.Genes
#there is CpG without gene associated, we remove them
sum(is.na(resCpGM.Genes$gene)) #339390
length(unique(resCpGM.Genes$locisID)) # 786277
resCpGM.Genes<-resCpGM.Genes[!is.na(gene)]
resCpGM.Genes

#in eQTL analysis there is CpG with several CpG, to avoid noise we select just 1 gene by CpG : the closest gene 
resCpGM.Genes[is.na(start.eQTR2),closest.Genes:=T,by=.(locisID)]
resCpGM.Genes[!is.na(start.eQTR2),closest.Genes:=distTSS==min(distTSS),by=.(locisID)]
resCpGM.Genes<-resCpGM.Genes[closest.Genes==T]
resCpGM.Genes

#~~ II) Calculate The CpG Score ~~
#CpGScore = DMCScore [-inf:inf] * Regulatory weight[0.5-1] * LinksWeight [0.5-1]

# #>0) filter no significant locis : pval> 0.001 and FC <20
# plot(density(resCpGM.Genes$FC))
# resCpGM.Genes<-resCpGM.Genes[pval>0.001,FC:=0]
# resCpGM.Genes<-resCpGM.Genes[FC<10,FC:=0]
# resCpGM.Genes[FC>0]

# > 1) DMCScore = FC * PvalWeight[0.1,1]
plot(density(resCpGM.Genes$FC))
#PvalWeight[0-1]
min(resCpGM.Genes$pval)
resCpGM.Genes[,PvalWeight:=(-log10(pval)/8)]
plot(density(resCpGM.Genes$PvalWeight))
min(resCpGM.Genes$PvalWeight)

resCpGM.Genes[,PvalWeight2:=1-pval]
plot(density(resCpGM.Genes$PvalWeight2))

#multiply the 2 score
resCpGM.Genes[,DMCScore:=PvalWeight*FC]
plot(density(resCpGM.Genes$DMCScore))

# > 2) Regulatory weight [0.5-2] = 0.5 + 1.5 * (TypeScore{0,0.4,0.66,1} + EnsRegScore{0,0.25,0.5,0.75,1}) 
#     TypeScore{0,0.2,0.33,0.5}
resCpGM.Genes[,TypeScore:=sapply(type,function(x){
  if(x%in%c(4,6)){
    score<-1
  }else if(x==5){
    score<-0.66
  }else if(x%in%1:3){
    score<-0.4
  }else{
    score<-0
  }
  return(score)
})]
plot(density(resCpGM.Genes$TypeScore))
#   EnsRegScore{0,0.25,0.5,0.75,1}
source("scripts/utils.R")
resCpGM.Genes[,EnsRegScore:=sapply(feature_type_name,function(x){
  vecX<-tr(x)
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
plot(density(resCpGM.Genes$EnsRegScore))
# Regulatory weight:
resCpGM.Genes[,RegWeight:=(0.5+1.5*((TypeScore+EnsRegScore)/2))]
plot(density(resCpGM.Genes[!duplicated(locisID)]$RegWeight))



# > 3) Links Weight[0.5-1] = 0.5 + 0.5 * (TypeScore{0,0.4,0.66,1} + EnsRegScore{0,0.25,0.5,0.75,1}) /2
#   - for CpG associated to a gene based on TSS proximity  
resCpGM.Genes[is.na(start.eQTR2),LinkScore:=sapply(abs(distTSS),function(x){
  if(x<1000){
    score<-1
  }else if(x<20000){
    score<-(3/(log10(x)))^0.5
  }else if (x<100000){
    score<-3/log10(x)
    
  }else{
    score<-(3/log10(x))^2
  }
  return(score)
})]
plot(density(resCpGM.Genes[is.na(start.eQTR2)]$LinkScore))

# - for CpG associated to a gene based on eQTL studies: 
#   *need eQTR : most precise eQTL region
eQTRs<-fread("../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv")
eQTRs<-eQTRs[,.(chr,start.eQTR2,end.eQTR2,start.eQTR,end.eQTR)]
eQTRs
resCpGM.Genes<-merge(resCpGM.Genes,eQTRs,all.x = TRUE,by=c("chr","start.eQTR2", "end.eQTR2"))[order(pval)]
resCpGM.Genes

#   * I made a function to calculate confidence of the links CpG-Gene for eQTL linked Gene : 
calcLinkScore<-function(pos,start,end,signif=NULL,maxSignif=0.001475078){
  
  # i) linkScore ~ CpG distance of the eQTL region (eQTR)
  if(pos > (start-50)& pos < (end+50)){
    score<-1
  }else {
    x<-min(abs(c(start-pos,end-pos)))
    score<-(log10(50)/(log10(x)))^0.5
  }
  
  if(!is.null(signif)){
    print('test')
    #si min pval >*1 si max pval > *0.75
    score<-score*(0.75+0.25*(signif/maxSignif))
  }
  
  
  # ii) ajust score if eQTL region too large (because too much uncertainty) 
  #    si >200, ajuste score sinon descend progressivement
  largeur<-end-start+1
  if(largeur<200){
    score<-score*1
  }else if(largeur<10000){
    score<-score*(0.5+0.5*(log10(200)/log10(largeur+1))^0.5)
  }else if(largeur<100000){
    score<-score*(0.5+0.5*(log10(200)/log10(largeur+1)))
  }else{
    score<-score*(0.5+0.5*(log10(200)/log10(largeur+1))^2)
  }
  
  
  return(score)
  
}
#   *calculate linkScore for eQTL-based CpG-Gene pairs :
resCpGM.Genes[!is.na(start.eQTR2),LinkScore:=calcLinkScore(pos[1],start.eQTR[1],end.eQTR[1]),by=.(locisID,gene,start.eQTR2)]
lines(density(resCpGM.Genes[!is.na(start.eQTR2)]$LinkScore),col=2)

resCpGM.Genes[,LinksWeight:=0.5+0.5*max(LinkScore),by=.(locisID,gene)]
resCpGM.Genes<-resCpGM.Genes[,isBestLinks:=LinkScore==max(LinkScore),by=c("locisID","gene")][isBestLinks==TRUE]
resCpGM.Genes<-unique(resCpGM.Genes)
plot(density(resCpGM.Genes$LinksWeight))

#Note : i would like to be more precise in the calculation of this links score


# > finally :
resCpGM.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]
plot(density(resCpGM.Genes$CpGScore))
sum(resCpGM.Genes$CpGScore<0)



#~~ summarize the CpGs Scores to a Gene Score  ~~
# problem : genes are linked to several CpG...
library(ggplot2)
resCpGM.Genes[,nCpG.Gene:=.N,by=.(gene)]
resCpGM.Genes[CpGScore!=0,nCpGSig.Gene:=.N,by=.(gene)]

resCpGM.Genes
ggplot(resCpGM.Genes[CpGScore!=0])+
  geom_bar(aes(x=nCpGSig.Gene))


resCpGM.Genes[CpGScore!=0,nGene.CpGSig:=.N,by=.(locisID)]
ggplot(resCpGM.Genes[CpGScore!=0])+
  geom_bar(aes(x=nGene.CpGSig))


# ...How to summarize the CpG score to a gene score ??
#to consider :
#- a big gene is most susceptible to have a lot of CpG, so more easy for it to have a high CpGscore => distrib more important
#- some CpGs in the same regulatory region, they are markers for the same regulation

#=>  the GeneScore need to integrate 1) the CpGScore max and 2) the distribution of the score, 3) the number of CpG associated with the score

#GeneScore[-inf:inf] = c*max(CpG) + (1-c)*q75(CpGScore)
#with c[0.5-1] ~ nCpG associated to the gene. 
# c=0.5 for a gene with a lot of CpG  => bigger coef for the q75Score because distrib most important than with few cpG 
# c=1 for a gene with 1 CpG => just the max(CpG) is important
#  c= 0.5+0.5*sqrt(1/nCpG)
resCpGM.Genes[,coef:=0.5+0.5*sqrt(1/nCpG.Gene)]


resCpGM.Genes[,GeneScore:=coef*extr(CpGScore)[1]+(1-coef)*median(CpGScore),by="gene"]

resGenesM<-unique(resCpGM.Genes,by="gene")[order(-GeneScore)]
resGenesM
resGenesM$gene[1:100]
plot(density(resGenesM$nCpG.Gene))
plot(density(resGenesM$coef))
plot(density(resGenesM$GeneScore))
sum(resGenesM$GeneScore<0)

#more precise : GeneScore2[-inf:inf]= moy ponder ci*CpgScore/sum(ci) avec ci~1/nCpG

moyPond<-function(CpGScores){
  CpGScores<-CpGScores[order(abs(CpGScores),decreasing = T)]
  coefs<-sapply(1:length(CpGScores), function(x)1/(x))
  
  return(sum(CpGScores*coefs)/sqrt(sum(coefs)))
}

resCpG.Genes[,GeneScore2:=moyPond(CpGScore),by="gene"]
resCpG.Genes
resGenesM<-unique(resCpG.Genes,by="gene")[order(-GeneScore2)]

resGenesM
resGenesM$gene[1:100]

plot(density(resGenesM$GeneScore2))
abline(v=10)

plot(-log10(resGenesM$pval),resGenesM$GeneScore2)
abline(h=2)
plot(resGenesM$FC,resGenesM$GeneScore2)
abline(h=10)
plot(resGenesM$DMC,resGenesM$GeneScore2)

plot(resGenesM$LinksWeight,resGenesM$GeneScore)

plot(as.factor(resGenesM$RegWeight),resGenesM$GeneScore)

#But don't take into account if the gene is multi-regulated.
#idea to improve : summarize CpGScore by regulatory region and sum this scores to have a GeneScore


library(clusterProfiler)
library(org.Hs.eg.db)

#OVER-REPRESENTATION TEST
#genes candidat : GeneScore >10


resGenesM<-unique(resCpGM.Genes,by="gene")[order(-GeneScore)]
resGenesM
plot(density(resGenesM$GeneScore))
abline(v=15)
genes15<-resGenesM[abs(GeneScore)>5]$gene
length(genes15) #251

rk15M <- enrichKEGG(gene         = bitr(genes15,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

dfk15M<-as.data.frame(rk15M)
nrow(dfk15M)
dfk15M
dotplot(rk15M,showCategory=30)


source("scripts/utils.R")
library(data.table)
dfk15M<-data.table(dfk15M)

dfk15M[,GeneScore.avg:=mean(resGenesM$GeneScore[resGenesM$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfk15M


dotplot(rk15M,x=dfk15M$GeneScore.avg[1:10],showCategory=10)
