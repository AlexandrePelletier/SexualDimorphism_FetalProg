
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

options(stringsAsFactors=F)
set.seed(12345)

library(data.table)


#differentially methylated datas :
resCpG<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_FC.FL_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",
              dec = ",")
resCpG<-resCpG[,locisID:=V1][,pos:=start-1][,meth.change:=FC][,-c("V1","start")][,.(locisID,chr,pos,pval,meth.change)]
resCpG


#~~ I) link CpGs to genes ~~

annot<-fread("../../ref/allCpG_annotation_genes_and_feature_regulatory_domain_070520.csv")
annot<-annot[!(is.na(gene)|gene=="")]
annot
# 2 types of CpG-Gene links
# - the closest TSS from a CpG is associated with them
annot[is.na(start.eQTR2)]

# - the CpG is in a region with eQTL (SNP associated with a change of expression of a gene)
annot[!is.na(start.eQTR2)]

#IMPROV1 : here is based on blood eQTL get in : https://gtexportal.org/home/datasets 
#file path is : GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt
#=> i plan to integrate also the meta_analysis provide by GTEx of all tissues to get 
# the variant_gene_pairs conserved across tissue (to improve confidence of the link)

resCpG.Genes<-merge(resCpG,annot,by=c("locisID","chr","pos"),all.x = T)
resCpG.Genes

#there is CpG without gene associated, we remove them
sum(is.na(resCpG.Genes$gene)) #339390
length(unique(resCpG.Genes$locisID)) # 786277
resCpG.Genes<-resCpG.Genes[!is.na(gene)]
resCpG.Genes

#in eQTL analysis there is CpG linked to several Gene, to avoid background noise
#we select just 1 gene by CpG : the closest gene
resCpG.Genes[is.na(start.eQTR2),closest.Genes:=T,by=.(locisID)]
resCpG.Genes[!is.na(start.eQTR2),closest.Genes:=distTSS==min(distTSS),by=.(locisID)]
resCpG.Genes<-resCpG.Genes[closest.Genes==T]
resCpG.Genes

#IMPROV2 : here i could choose the best links according to better metrics than distTSS (number and signifiance of eQTL-Gene link in the region,...)


#~~ II) Calculate The CpG Score ~~
#CpGScore = DMCScore [-inf:inf] * Regulatory weight[0.5-2] * LinksWeight [0.5-1]

# > 1) DMCScore = meth.change * PvalWeight[0,1]
plot(density(resCpG.Genes$meth.change))
#PvalWeight[0-1]
plot(density(log10(resCpG.Genes$pval)))
resCpG.Genes[,PvalWeight:=(-log10(pval)/8)] #
plot(density(resCpG.Genes$PvalWeight))

#multiply the 2 score
resCpG.Genes[,DMCScore:=PvalWeight*meth.change]
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
calcLinkScore<-function(pos,start,end){
  
  # i) linkScore ~ CpG distance of the eQTL region (eQTR)
  if(pos > (start-50)& pos < (end+50)){
    score<-1
  }else {
    x<-min(abs(c(start-pos,end-pos)))
    score<-(log10(50)/(log10(x)))^0.5
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

#because eQTR can overlap...
resCpG.Genes<-resCpG.Genes[,nEQTR.CpgGeneLink:=.N,by=c("locisID","gene")]
plot(density(resCpG.Genes$nEQTR.CpgGeneLink))

#...we select the best linkScore by cpg-gene pair
resCpG.Genes<-resCpG.Genes[,isBestLinks:=LinkScore==max(LinkScore),by=c("locisID","gene")][isBestLinks==TRUE]
resCpG.Genes<-unique(resCpG.Genes)
plot(density(resCpG.Genes$LinksWeight))

#IMPROV3 : maybe i could be more precise in the calculation of this links score 
#(take into account the signifiance of the link + the number of eQTL, and also the proximity of the CpG with the eQTL)


# > finally :
resCpG.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]
plot(density(resCpG.Genes$CpGScore))



# VALIDATION : 
lines(density(resCpG.Genes[pval<0.001]$CpGScore),col=2)
lines(density(resCpG.Genes[meth.change>20]$CpGScore),col=3)
lines(density(resCpG.Genes[pval<0.001& meth.change>20]$CpGScore),col=4) #

nrow(resCpG.Genes[pval>0.05 & CpGScore>10])
plot(density(resCpG.Genes[pval>0.01 ]$CpGScore))

nrow(resCpG.Genes[meth.change<15 & CpGScore>10])
plot(density(resCpG.Genes[abs(meth.change)<20 ]$CpGScore))

# some validation plots : 
# my score highlight mostly CpG with high pvalue and meth.change : 
plot(-log10(resCpG.Genes$pval),resCpG.Genes$CpGScore) 
plot(resCpG.Genes$meth.change,resCpG.Genes$CpGScore)
plot(resCpG.Genes$DMC,resCpG.Genes$CpGScore)
# my score highlight in a second level :
#   - CpGs in regulatory region :
plot(as.factor(resCpG.Genes$RegWeight),resCpG.Genes$CpGScore)

#   - CpGs where the link to the gene is strong :
plot(resCpG.Genes$LinksWeight,resCpG.Genes$CpGScore)

# QU1 : Is there better way to "statistically" validate my score ?



#~~ III) summarize the CpGs Scores to a Gene Score  ~~

# problem : genes are linked to several CpG...
library(ggplot2)
resCpG.Genes[,nCpG.Gene:=.N,by=.(gene)]
resCpG.Genes[CpGScore>10,nCpGSig.Gene:=.N,by=.(gene)]

ggplot(resCpG.Genes)+
  geom_bar(aes(x=nCpG.Gene))

ggplot(resCpG.Genes[CpGScore>10])+
  geom_bar(aes(x=nCpGSig.Gene))

# ...How to summarize the CpG score to a Gene score ??
#to consider :
#- a big gene is most susceptible to have a lot of CpG associated, so more easy for it to have a high CpGscore => distrib more important
#- some CpGs in the same regulatory region, they are markers for the same regulation

#=>  the GeneScore need to integrate :
# 1) the max CpGScore linked to the gene : max(CpGscore)
# 2) the number of CpG associated with the Gene  : nCpG
# 3) the distribution of their CpGscores associated with it : median(CpGscore)

#GeneScore[-inf:inf] = c*max(CpG) + (1-c)*median(CpGScore)
#with c[0.5-1] ~ nCpG associated to the gene. 
# for a gene with a lot of CpG, the CpG score distribution is more important => c=0.5, the maximum and the median have the same weight
# for a gene with 2 CpG => c ~= 1, the max(CpGScore) is more important
#  c= 0.5+0.5*sqrt(1/nCpG)
resCpG.Genes[,coef:=0.5+0.5*sqrt(1/nCpG.Gene)]
plot(density(resGenes$coef))
plot(unique(resCpG.Genes$coef),unique(resCpG.Genes$nCpG.Gene),log="y")

resCpG.Genes[,GeneScore:=coef*extr(CpGScore)[1]+(1-coef)*median(CpGScore),by="gene"]


plot(density(resCpG.Genes[,medCpGScore:=median(CpGScore),by="gene"]$medCpGScore))

resGenes<-unique(resCpG.Genes,by="gene")[order(-GeneScore)]
resGenes
resGenes$gene[1:100]
plot(density(resGenes$nCpG.Gene))

plot(density(resGenes$GeneScore))
abline(v=15)



# VALIDATION : 

# some validation plots : 
# the score highlight Gene with high significant  CpG 
plot(-log10(resGenes$pval),resGenes$GeneScore)

# the score highligh also high methylation change, and in the same directory (hypo : negative value, hyper : postive value)
plot(resGenes$meth.change,resGenes$GeneScore)
plot(resGenes$DMC,resGenes$GeneScore) 
#We see for some gene opposite score :  meth.change negative  but a GeneScore positive, and inversely
#Is possible if there are 2 signif CpG with opposite methylation in a gene 

#Qu ) How integrate Opposite value ?


#IMPROV 4 :  take into account the length of the gene to normalize the gene-score


library(clusterProfiler)
library(org.Hs.eg.db)


#OVER-REPRESENTATION TEST
#genes candidat : GeneScore >15
plot(density(resGenes$GeneScore))
abline(v=15)
genesGS<-resGenes[GeneScore>15]$gene

length(genesGS) #1726

rkGS <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1)


dfkGS<-as.data.frame(rkGS)
nrow(dfkGS) #10
dfkGS
dotplot(rkGS,showCategory=30)

source("scripts/utils.R")
library(data.table)
dfkGS<-data.table(dfkGS)
dfkGS[,GeneScore.avg:=mean(resGenes$GeneScore[resGenes$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGS


dotplot(rkGS,x=dfkGS$GeneScore.avg)
emapplot(rkGS)

saveRDS(rkGS,"analyses/withoutIUGR/2020-05-24_CF.LF_clusterProfObject_ORA_GeneScore1_thr15.rds")
write.csv2(as.data.frame(setReadable(rkGS,org.Hs.eg.db,keyType = "ENTREZID")),"analyses/withoutIUGR/2020-05-24_CF.LF_res_ORA_GeneScore1_thr15.csv")

# COMPARISON WITH A PVALUE-BASED GENE SELECTION
#with 1726 genes
resCpG.GenesF<-resCpG.Genes[meth.change>0][order(pval)]
genesPval<-unique(resCpG.GenesF$gene)[1:1726]

intersect(genesPval,genesGS)
setdiff(genesGS,genesPval)
setdiff(genesPval,genesGS)


rkPval <- enrichKEGG(gene         = bitr(genesPval,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1)
dfkPval<-as.data.frame(rkPval)
nrow(dfkPval) #5
dfkPval
dotplot(rkPval,showCategory=30)

library(patchwork)
p1<-dotplot(rkGS, showCategory=30)+ggtitle("Methyl-Gene Score based")+theme(plot.title = element_text(size = 11, face = "bold"))

p2<-dotplot(rkPval, showCategory=30)+ggtitle("Pvalue based")+theme(plot.title = element_text(size = 11, face = "bold"))

p1+p2+plot_annotation("Most methylation-affecting candidat genes","top 1726 genes. KEGG pathways Over-representation test, p.adj_cutoff = 0.1")


candidat_genes.list<-list(GeneScore=bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                          Pvalue=bitr(genesPval,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID)

res.compa<-compareCluster(candidat_genes.list,
                          fun = "enrichKEGG",
                          organism="hsa",
                          pvalueCutoff=0.1) #can compare in other ontology (GO, Reactome, DOSE..)
nrow(as.data.frame(res.compa)) #15
dotplot(res.compa,showCategory =15)

emapplot(res.compa,pie="count",showCategory =15 )

#IMPROV 5 : test the reliability and robustness of my score;
#ideas : 
# - change number of genes / geneScore threshold => Are the pathway enrichment robust ? or highly variable
# - test with simulated data (random permutation..)
# - test with disease signature methylation data

#QU2 : Have you any idea what relevant test could i perform to validate the interest/reliability of my geneScore  ?



#ANNEXE
#OTHER TEST OF CpG-score SUMMARIZING : 
#more precise than just the median , my first idea was symply with a ponderate mean:
#Genescore = sum(c_i*CpgScore)/sum(c_i) avec c_1 = 1 for the max(CpGScore) , then c_i = 1/rank(CpGScore_i)
#but the geneScore for gene with numerous CpG was too much push down (because lot of CpG with CpG score=0)

#so i decided to alleviate the weigth of the number of CpG like this :
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







