
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

#2020-06-01 CpGScore v2

options(stringsAsFactors=F)
set.seed(12345)

library(data.table)


#differentially methylated datas :
resCpGF<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_FC.FL_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",
              dec = ",")
resCpGF<-resCpGF[,locisID:=V1][,pos:=start-1][,meth.change:=FC][,-c("V1","start")][,.(locisID,chr,pos,pval,meth.change)]
resCpGF


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
resCpGF.Genes[nGene.CpG>50]
#there is duplicated row, we conseve unique locisID-gene pairing 
resCpGF.Genes<-unique(resCpGF.Genes)
hist(resCpGF.Genes$nGene.CpG)
resCpGF.Genes[nGene.CpG>50]
resCpGF.Genes[locisID==8342]
#there is multiple row for the same eQTR, just the eQTL-dist change, we select the closest
resCpGF.Genes[in.eQTR==TRUE,isClosest:=eQTL_dist==min(eQTL_dist),by=.(locisID,gene)]
resCpGF.Genes[in.eQTR==FALSE,isClosest:=TRUE]

resCpGF.Genes<-resCpGF.Genes[isClosest==TRUE]
resCpGF.Genes[,nGene.CpG:=.N,by=.(locisID)]
hist(resCpGF.Genes$nGene.CpG)
resCpGF.Genes[nGene.CpG>10]

resCpGF.Genes


#~~ II) Calculate The CpG Score ~~
#CpGScore = DMCScore [-inf:inf] * Regulatory weight[0.5-2] * LinksWeight [0.5-1]

# > 1) DMCScore = meth.change * PvalWeight[0,1]
plot(density(resCpGF.Genes$meth.change))
#PvalWeight[0-1]
plot(density(log10(resCpGF.Genes$pval)))
resCpGF.Genes[,PvalWeight:=(-log10(pval)/3)] #
plot(density(resCpGF.Genes$PvalWeight))

#multiply the 2 score
resCpGF.Genes[,DMCScore:=PvalWeight*meth.change]
plot(density(resCpGF.Genes$DMCScore))

# > 2) Regulatory weight [0.5-2] = 0.5 + 1.5 * (TypeScore{0,0.4,0.66,1} + EnsRegScore{0,0.25,0.5,0.75,1})/2 
#     TypeScore{0,0.4,0.66,1} 
resCpGF.Genes[,TypeScore:=sapply(type,function(x){
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
plot(density(resCpGF.Genes$TypeScore))

#   EnsRegScore{0,0.25,0.5,0.75,1}
resCpGF.Genes[,EnsRegScore:=sapply(feature_type_name,function(x){
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
plot(density(resCpGF.Genes$EnsRegScore))
# Regulatory weight:
resCpGF.Genes[,RegWeight:=(0.4+1.6*((TypeScore+EnsRegScore)/2))]

plot(density(resCpGF.Genes[!duplicated(locisID)]$RegWeight))



# > 3) Links Weight[0.5-1] = 0.5+0.5*max(LinkScore [0-1])
#   - for CpG associated to a gene based on TSS proximity  
resCpGF.Genes[in.eQTR==FALSE,LinkScore:=sapply(abs(tss_dist),function(x){
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
plot(density(resCpGF.Genes[in.eQTR==FALSE]$LinkScore))

# - for CpG associated to a gene based on eQTL studies: 
#de base, linksscore =1

resCpGF.Genes[in.eQTR==TRUE,LinkScore:=1] 

#but, filter for low signif asso, and add bonus if highly signif eQTR and close to best eQTL
plot(density(na.omit(resCpGF.Genes$eQTL_dist)))
resCpGF.Genes[in.eQTR==TRUE,distScore:=sapply(abs(eQTL_dist),function(x){
  if(x<500){
    distScore<-1 
  }else {
    distScore<-1-0.5*x/2500
  }
  
  return(distScore)
})]
plot(density(na.omit(resCpGF.Genes$avg.mlog10.pv.eQTLs)))

plot(density(na.omit(resCpGF.Genes$avg.mlog10.pv.eQTLs)))
resCpGF.Genes[in.eQTR==TRUE,RegScore:=avg.mlog10.pv.eQTLs*distScore] #integrate with RegScore : mean of pval of eQTL in the region

plot(density(na.omit(resCpGF.Genes$RegScore))) 
#remove low signif links (<4) and add bonus if 
abline(v=4)
summary(na.omit(resCpGF.Genes$RegScore))
resCpGF.Genes<-resCpGF.Genes[RegScore>4|in.eQTR==F]
#and add bonus ; +1*regScore/220

resCpGF.Genes[in.eQTR==TRUE,LinkScore:=LinkScore+1*RegScore/max(RegScore)]
plot(density(resCpGF.Genes$LinkScore))
resCpGF.Genes[LinkScore>1.5]
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.001   5.056   7.971  12.019  13.895 221.354 

#because of some tss_distance dismatch between eQTRlinked and not, there is different linkScore
resCpGF.Genes[,nEQTR.CpgGeneLink:=.N,by=c("locisID","gene")]
plot(density(resCpGF.Genes$nEQTR.CpgGeneLink))
resCpGF.Genes[nEQTR.CpgGeneLink>1]
#...we select the best linkScore by cpg-gene pair
resCpGF.Genes<-resCpGF.Genes[,isBestLinks:=LinkScore==max(LinkScore),by=c("locisID","gene")][isBestLinks==TRUE]
resCpGF.Genes<-unique(resCpGF.Genes[,-c("isClosest","isBestLinks")])
#for same linksWeught
resCpGF.Genes<-unique(resCpGF.Genes,by=c("locisID","gene"))
#and finally calculate the linksWeight
resCpGF.Genes[,LinksWeight:=0.3+1.2*(LinkScore/2)]
plot(density(resCpGF.Genes$LinksWeight))





# > finally :
resCpGF.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]

plot(density(resCpGF.Genes$CpGScore))



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


#~~ III) summarize the CpGs Scores to a Gene Score  ~~

# problem : genes are linked to several CpG...
library(ggplot2)
resCpGF.Genes[,nCpG.Gene:=.N,by=.(gene)]
#resCpGF.Genes[pval<10^-3,nCpGSig.Gene:=.N,by=.(gene)]
#[corrected]
resCpGF.Genes[,nCpGSig.Gene:=sum(pval<10^-3),by=.(gene)]
resCpGF.Genes[gene=="FOXO3"][order(-CpGScore)]
plot(density(resCpGF.Genes$CpGScore))
abline(v=10)
resCpGF.Genes[,nHiCpGScore.Gene:=sum(CpGScore>10),gene]
#[end correction]
ggplot(resCpGF.Genes)+
  geom_bar(aes(x=nCpG.Gene))

ggplot(resCpGF.Genes)+
  geom_bar(aes(x=nCpGSig.Gene))

# ...How to summarize the CpG score to a Gene score ??
#to consider :
#- a big gene is most susceptible to have a lot of CpG associated, so more easy for it to have a high CpGscore => distrib more important
#- some CpGs in the same regulatory region, they are markers for the same regulation

#=>  the GeneScore need to integrate :
# 1) the max CpGScore linked to the gene : max(CpGscore)
# 2) the number of CpG associated with the Gene  : nCpG
# 3) the distribution of their CpGscores associated with it : median(CpGscore)

#GeneScore[-200:200] = c*max(CpG) + (1-c)*median(CpGScore)
#with c[0.5-1] ~ nCpG associated to the gene. 
# for a gene with a lot of CpG, the CpG score distribution is more important => c=0.5, the maximum and the median have the same weight
# for a gene with 2 CpG => c ~= 1, the max(CpGScore) is more important
#  c= 0.5+0.5*sqrt(1/nCpG)
resCpGF.Genes[,coef:=0.5+0.5*sqrt(1/nCpG.Gene)]
plot(unique(resCpGF.Genes$coef),unique(resCpGF.Genes$nCpG.Gene))

resCpGF.Genes[,GeneScore:=coef*CpGScore[which.max(abs(CpGScore))]+(1-coef)*median(CpGScore),by="gene"]

fwrite(resCpGF.Genes,"analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev2.csv",sep=";")

resGenes<-unique(resCpGF.Genes[order(pval)],by="gene")[order(-GeneScore)]

plot(density(resGenes$GeneScore))
abline(v=40)

# VALIDATION : 

# some validation plots : 
# the score highlight Gene with high significant  CpG 
plot(-log10(resGenes$pval),resGenes$GeneScore)

# the score highligh also high methylation change, and in the same directory (hypo : negative value, hyper : postive value)
plot(resGenes$meth.change,resGenes$GeneScore)
plot(resGenes$DMC,resGenes$GeneScore) 
#We see for some gene opposite score :  meth.change negative  but a GeneScore positive, and inversely
#Possible if 2 signif CpG with opposite methylation change in a gene 
resGenes[meth.change>20&GeneScore<(-5)]
resCpGF.Genes[gene=="L2HGDH"] #1CpG avec score = 8 et un aute avec score =-8
median(resCpGF.Genes[gene=="L2HGDH"]$CpGScore)
#QU2 : any idea to better integrate Opposite methylation change ?

#IMPROV4 :  take into account the length of the gene to normalize the gene-score




#OVER-REPRESENTATION TEST
#genes candidat : GeneScore >15
library(clusterProfiler)
library(org.Hs.eg.db)
plot(density(resGenes$GeneScore))
abline(v=40)
genesGS<-resGenes[GeneScore>40]$gene

length(genesGS) #991

rkGS <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1)


dfkGS<-as.data.frame(rkGS)
nrow(dfkGS) #7
dfkGS
dotplot(rkGS,showCategory=30)

dfkGS<-data.table(dfkGS)
dfkGS[,GeneScore.avg:=mean(resGenes$GeneScore[resGenes$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGS
dotplot(rkGS,x=dfkGS$GeneScore.avg)
emapplot(rkGS)

#IN PATHWAY MAP

#need genelist.entrez
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

rkGS<-readRDS("analyses/withoutIUGR/2020-05-24_CF.LF_clusterProfObject_ORA_GeneScore1_thr15.rds")

library("pathview")
as.data.frame(rkGS)

#Longevity regulating pathway, multiple species :
dir.create("analyses/withoutIUGR/pathwaysMap",recursive = T)
pathw <- pathview(gene.data  = geneList.Entrez,
                  pathway.id = "hsa04213", 
                  species    = "hsa",
                  limit      = list(gene=max(abs(geneList)), cpd=1),
                  kegg.dir ="analyses/withoutIUGR/pathwaysMap" 
)




# COMPARISON WITH A PVALUE-BASED GENE SELECTION
#with 1726 genes
resCpGF.GenesF<-resCpGF.Genes[meth.change>0][order(pval)]
genesPval<-unique(resCpGF.GenesF$gene)[1:1697]

intersect(genesPval,genesGS) #1164 common genes
setdiff(genesGS,genesPval) #533 genes diff
setdiff(genesPval,genesGS)


rkPval <- enrichKEGG(gene         = bitr(genesPval,fromType = "SYMBOL",toType = "ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1)

dfkPval<-as.data.frame(rkPval)
nrow(dfkPval) #4
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
nrow(as.data.frame(res.compa)) #11
dotplot(res.compa,showCategory =11)

emapplot(res.compa,pie="count",showCategory =11 )

#IMPROV 5 : test the reliability and robustness of my score;
#ideas : 
# - change number of genes / geneScore threshold => Are the pathway enrichment robust ? or highly variable
# - test with simulated data (random permutation..)
# - test with disease signature methylation data

#QU3 : Have you any idea what relevant test could i perform to validate the interest/reliability of my geneScore  ?



#NEW GENESCORE pour VALORISER GENE AVEC ++ High CpGScore
#new geneSCORE
resCpGF.Genes[,ismaxCpGScore:=CpGScore==max(CpGScore),by=.(locisID,gene)]
resCpGF.Genes<-resCpGF.Genes[ismaxCpGScore==T][order(gene,-CpGScore)]
resCpGF.Genes
resCpGF.Genes[,GeneScore2:=CpGScore[which.max(abs(CpGScore))],by="gene"]
resCpGF.Genes[nCpG.Gene>1,GeneScore2:=GeneScore2+mean(CpGScore[2:.N]),by="gene"]
resGenes2<-unique(resCpGF.Genes,by="gene")[order(-GeneScore2)]

resGenes2
resGenes2$gene[1:100]

plot(density(resGenes2$GeneScore2))
abline(v=60)
#OVER-REPRESENTATION TEST
#genes candidat : GeneScore2 >25
library(clusterProfiler)
library(org.Hs.eg.db)
plot(density(resGenes2$GeneScore2))
abline(v=60)
genesGS<-resGenes2[GeneScore2>60]$gene

length(genesGS) #1841

rkGS <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1)


dfkGS<-as.data.frame(rkGS)
nrow(dfkGS) #18
dfkGS
dotplot(rkGS,showCategory=30)

dfkGS<-data.table(dfkGS)
dfkGS[,GeneScore.avg:=mean(resGenes$GeneScore[resGenes$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGS
dotplot(rkGS,x=dfkGS$GeneScore.avg,showCategory=18)
emapplot(rkGS)


#genes candidat : top100
genesGS<-resGenes2$gene[1:100]
plot(density(resGenes2$GeneScore2))
abline(v=min(resGenes2[gene%in%genesGS]$GeneScore2))


rkGS <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1)


dfkGS<-as.data.frame(rkGS)
nrow(dfkGS) #3
dfkGS
dotplot(rkGS,showCategory=30)

dfkGS<-data.table(dfkGS)
dfkGS[,GeneScore.avg:=mean(resGenes$GeneScore[resGenes$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGS
dotplot(rkGS,x=dfkGS$GeneScore.avg,showCategory=17)
emapplot(rkGS)





#OTHER TEST OF CpG-score SUMMARIZING : 
#so i decided to alleviate the weigth of the number of CpG like this :


#geneSCORE3 : the sum

resCpGF.Genes[,GeneScore3:=sum(CpGScore),by="gene"]
plot(density(na.omit(log10(resCpGF.Genes$CpGScore))))
resCpGF.Genes[CpGScore==0]

resCpGF.Genes[,GeneScore.inv:=sum(1/(abs(CpGScore)+1)),by="gene"]
resCpGF.Genes[,sqrt.inv.geneScore:=sum(sqrt(1/(abs(CpGScore)+1))),by="gene"]
resCpGF.Genes[,nCpGWeight:=1/log1p(GeneScore.inv)]

resCpGF.Genes[,nCpGWeight2:=1+(1/sqrt.inv.geneScore)]

resCpGF.Genes[,nCpGWeight3:=1/GeneScore.inv]

resCpGF.Genes[nCpGWeight>10 & nCpG.Gene>1]
resGenes2<-unique(resCpGF.Genes[order(-CpGScore)],by="gene")

plot(resGenes2$nCpG.Gene,resGenes2$nCpGWeight)

plot(resGenes2$nCpG.Gene,resGenes2$nCpGWeight)
plot(resGenes2$nCpG.Gene,resGenes2$nCpGWeight2)

summary(resGenes2$nCpG.Gene)
summary(resGenes2$nCpGWeight3)
plot(resGenes2$nCpG.Gene,resGenes2$nCpGWeight3)

resCpGF.Genes[,nCpGWeight.norm:=sqrt(log1p(nCpGWeight)),by="gene"]


resCpGF.Genes[,nCpGWeight.norm2:=sqrt(log1p(GeneScore.inv)),by="gene"]

resCpGF.Genes[,nCpGWeight.norm3:=sqrt(log1p(nCpGWeight3)),by="gene"]
plot(resGenes2$nCpG.Gene,resGenes2$nCpGWeight.norm3)
resCpGF.Genes[,nCpGWeight.norm3.2:=sqrt(nCpGWeight3)]

abline(v=20) #need def opti ncpg.sig pour avoir confiance dans le gene
plot(resGenes2$nCpG.Gene,resGenes2$nCpGSig.Gene)
plot(resGenes2[nCpG.Gene<0]$nCpG.Gene,resGenes2[nCpG.Gene<50]$nCpGSig.Gene)

#=> 7 CpG link to a gene is optimal to have confidence of the meth.change impact
plot(resGenes2$nCpG.Gene,resGenes2$nCpGWeight.norm3)
abline(v=7)
resGenes2[nCpGWeight.norm3>1.6][order(-nCpGWeight.norm3)]
summary(resGenes2[nCpG.Gene==7]$nCpGWeight.norm3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3746  0.4519  0.5021  0.5213  0.5711  0.8583 

plot(resGenes2[nCpG.Gene<50]$nCpG.Gene,resGenes2[nCpG.Gene<50]$nCpGWeight.norm3.2)
summary(resGenes2$nCpGWeight.norm3.2)
abline(v=7) #need def opti ncpg.sig pour avoir confiance dans le gene

resGenes2[nCpGWeight.norm3.2>5][order(-nCpGWeight.norm3)]
summary(resGenes2[nCpG.Gene==7]$nCpGWeight.norm3.2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3881  0.4760  0.5354  0.5634  0.6210  1.0435 

resCpGF.Genes[,GeneScore4:=sum(CpGScore)*nCpGWeight.norm,by="gene"]
resCpGF.Genes[,GeneScore5:=sum(CpGScore)*nCpGWeight.norm2,by="gene"]

resCpGF.Genes[,GeneScore6:=sum(CpGScore)*nCpGWeight2,by="gene"]
resCpGF.Genes[,GeneScore7:=sum(CpGScore)*nCpGWeight.norm3,by="gene"]
resCpGF.Genes[,GeneScore8:=sum(CpGScore)*(0.2+nCpGWeight.norm3),by="gene"]
resCpGF.Genes[,GeneScore9:=sum(CpGScore)*(0.4+nCpGWeight.norm3.2),by="gene"]

#determine fix coef  (0.1,0.2, 0.3, 0.4 ?); 
resGenes2<-unique(resCpGF.Genes[order(-CpGScore)],by="gene")
plot(resGenes2$nHiCpGScore.Gene,resGenes2$GeneScore9)
plot(resGenes2$nCpGWeight.norm3.2,resGenes2$GeneScore9)
resGenes2[GeneScore9>380&nHiCpGScore.Gene<5]
summary(resGenes2$nCpG.Gene)
plot(density(resGenes2$CpGScore))
abline(v=10)
#a good normal gene : 16 cpg, with 5 high cpgscore
resGenes2
resGenes2[nCpG.Gene==16&nHiCpGScore.Gene==5]


plot(density(resGenes2$GeneScore3))
plot(density(resGenes2$GeneScore4))
plot(density(resGenes2$GeneScore5))
plot(density(resGenes2$GeneScore6))
plot(density(resGenes2$GeneScore7))
plot(density(resGenes2$GeneScore8))
plot(density(resGenes2$GeneScore9))
#♦cutoff : 
abline(v=100)
length(resGenes2[GeneScore9>100]$gene) #1706

resGenes2[,cpg.sig.ratio:=nCpGSig.Gene/nCpG.Gene]
plot(resGenes2$nCpG.Gene,resGenes2$nCpGWeight.norm)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore3)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore4)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore5)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore6)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore7)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore8)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore9)
plot(resGenes2$nCpG.Gene,resGenes2$GeneScore9,col=(as.numeric(resGenes2$nCpGSig.Gene>3)+1))


plot(resGenes2$nCpG.Gene,resGenes2$GeneScore8,col=as.numeric(resGenes2$locisID==1072294))


plot(resGenes2[nCpG.Gene<100]$nCpG.Gene,resGenes2[nCpG.Gene<100]$GeneScore3)
plot(resGenes2[nCpG.Gene<100]$nCpG.Gene,resGenes2[nCpG.Gene<100]$GeneScore4)
plot(resGenes2[nCpG.Gene<100]$nCpG.Gene,resGenes2[nCpG.Gene<100]$GeneScore6)
plot(resGenes2[nCpG.Gene<100]$nCpG.Gene,resGenes2[nCpG.Gene<100]$GeneScore7)
plot(resGenes2[nCpG.Gene<100]$nCpG.Gene,resGenes2[nCpG.Gene<100]$GeneScore8)
plot(resGenes2[nCpG.Gene<100]$nCpG.Gene,resGenes2[nCpG.Gene<100]$GeneScore9)

plot(as.factor(resGenes2[nCpG.Gene<25]$nCpG.Gene),resGenes2[nCpG.Gene<25]$GeneScore9)
plot(as.factor(resGenes2[nCpG.Gene<25]$nCpG.Gene),resGenes2[nCpG.Gene<25]$GeneScore9)

abline(h=300)
abline(v=20)

resGenes2[nCpG.Gene==1&GeneScore9>=200]

resGenes2
genesDiff<-setdiff(resGenes2[order(-GeneScore4)]$gene[1:1000],resGenes2[order(-GeneScore3)]$gene[1:1000])
length(genesDiff) #103
resGenes2[gene%in% genesDiff & rank(GeneScore4)-rank(GeneScore3)>quantile(rank(GeneScore4)-rank(GeneScore3),0.75)] #> 243 rank de diff pour les 25% meilleurs decallage de gene
#remonte bien SETD9, qui a un bon methchange mais relativement peu de cpg


genesDiff2<-setdiff(resGenes2[order(-GeneScore5)]$gene[1:1000],resGenes2[order(-GeneScore3)]$gene[1:1000])
length(genesDiff2)
q75.2<-quantile(rank(resGenes2$GeneScore5)-rank(resGenes2$GeneScore3),0.75) #> 382 rank de diff pour les 25% meilleurs decallage de gene
resGenes2[gene%in% genesDiff2 & rank(GeneScore5)-rank(GeneScore3)>q75.2] #remonte si 1cpg+++signif, si ++ cpgsig /gene
#mais n'avantage pas assez les genes avec peu de cpg

genesDiff3<-setdiff(resGenes2[order(-GeneScore6)]$gene[1:1000],resGenes2[order(-GeneScore3)]$gene[1:1000])
length(genesDiff3) #18
q75.2<-quantile(rank(resGenes2$GeneScore5)-rank(resGenes2$GeneScore3),0.75) #> 382 rank de diff pour les 25% meilleurs decallage de gene
resGenes2[gene%in% genesDiff2 & rank(GeneScore5)-rank(GeneScore3)>q75.2] 

genesDiff4<-setdiff(resGenes2[order(-GeneScore9)]$gene[1:1000],resGenes2[order(-GeneScore3)]$gene[1:1000])
length(genesDiff4) #151/1000
q75.3<-quantile(rank(resGenes2$GeneScore9)-rank(resGenes2$GeneScore3),0.75) #> 322 rank de diff pour les 25% meilleurs decallage de gene
subres<-resGenes2[gene%in% genesDiff4 & rank(GeneScore9)-rank(GeneScore3)>q75.3] 
summary(subres[,.(pval,meth.change,nCpG.Gene,nCpGSig.Gene)])
# pval            meth.change      nCpG.Gene      nCpGSig.Gene   
# Min.   :7.010e-07   Min.   :27.50   Min.   : 1.00   Min.   : 1.000  
# 1st Qu.:1.936e-04   1st Qu.:35.18   1st Qu.:12.50   1st Qu.: 1.000  
# Median :5.598e-04   Median :39.37   Median :25.00   Median : 1.000  
# Mean   :1.233e-03   Mean   :40.08   Mean   :21.08   Mean   : 1.916  
# 3rd Qu.:1.645e-03   3rd Qu.:44.82   3rd Qu.:29.00   3rd Qu.: 2.000  
# Max.   :1.097e-02   Max.   :59.80   Max.   :43.00   Max.   :20.000
#c'est divers ==> le score ne privéligie pas un certain type de gene (e.g, gene a 1 cpg)


plot(resGenes2[nCpG.Gene<20]$nCpG.Gene,resGenes2[nCpG.Gene<20]$GeneScore3,col=3)
points(resGenes2[nCpG.Gene<20]$nCpG.Gene,resGenes2[nCpG.Gene<20]$GeneScore9,col=2)
points(resGenes2[nCpG.Gene<20&nCpGSig.Gene>5]$nCpG.Gene,resGenes2[nCpG.Gene<20&nCpGSig.Gene>5]$GeneScore9,col=4)


#OVER-REPRESENTATION TEST
#with GeneScore3 
library(clusterProfiler)
library(org.Hs.eg.db)
plot(density(resGenes2$GeneScore3))
abline(v=300)
genesGS<-resGenes2[GeneScore3>300]$gene

length(genesGS) #1212

rkGS <- enrichKEGG(gene         = bitr(genesGS,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)


dfkGS<-as.data.frame(rkGS)

nrow(dfkGS) #25
dfkGS
dotplot(rkGS,showCategory=30)

dfkGS<-data.table(dfkGS)
dfkGS[,GeneScore.avg:=mean(resGenes2$GeneScore3[resGenes2$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dfkGS
dotplot(rkGS,x=dfkGS$GeneScore.avg,showCategory=25)
emapplot(rkGS)
upsetplot(rkGS)
enrichplot::upsetplot(rkGS,n=25)
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
