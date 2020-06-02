
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
resCpG<-fread("analyses/withoutIUGR/2020-04-16_res_locis_in_FC.FL_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",
              dec = ",")
resCpG<-resCpG[,locisID:=V1][,pos:=start-1][,meth.change:=FC][,-c("V1","start")][,.(locisID,chr,pos,pval,meth.change)]
resCpG


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
resCpG.Genes<-merge(resCpG,cpgs.genes,by=c("locisID","chr","pos"),all.x = T)
resCpG.Genes

#there is CpG without gene associated, we remove them
sum(is.na(resCpG.Genes$gene)) #315608
length(unique(resCpG.Genes$locisID)) # 786277
resCpG.Genes<-resCpG.Genes[!is.na(gene)]
resCpG.Genes

#in eQTL analysis there is CpG linked to several Gene : 
resCpG.Genes[,nGene.CpG:=.N,by=.(locisID)]
hist(resCpG.Genes$nGene.CpG,breaks = 50)
resCpG.Genes[nGene.CpG>50]
#there is duplicated row, we conseve unique locisID-gene pairing 
resCpG.Genes<-unique(resCpG.Genes)
hist(resCpG.Genes$nGene.CpG)
resCpG.Genes[nGene.CpG>50]
resCpG.Genes[locisID==8342]
#there is multiple row for the same eQTR, just the eQTL-dist change, we select the closest
resCpG.Genes[in.eQTR==TRUE,isClosest:=eQTL_dist==min(eQTL_dist),by=.(locisID,gene)]
resCpG.Genes[in.eQTR==FALSE,isClosest:=TRUE]

resCpG.Genes<-resCpG.Genes[isClosest==TRUE]
resCpG.Genes[,nGene.CpG:=.N,by=.(locisID)]
hist(resCpG.Genes$nGene.CpG)
resCpG.Genes[nGene.CpG>10]

resCpG.Genes


#~~ II) Calculate The CpG Score ~~
#CpGScore = DMCScore [-inf:inf] * Regulatory weight[0.5-2] * LinksWeight [0.5-1]

# > 1) DMCScore = meth.change * PvalWeight[0,1]
plot(density(resCpG.Genes$meth.change))
#PvalWeight[0-1]
plot(density(log10(resCpG.Genes$pval)))
resCpG.Genes[,PvalWeight:=(-log10(pval)/3)] #
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
    score<-0.75
  }else if(x%in%1:3){
    score<-0.5
  }else{
    score<-0
  }
  return(score)
})]
plot(density(resCpG.Genes$TypeScore))

#   EnsRegScore{0,0.25,0.5,0.75,1}
resCpG.Genes[,EnsRegScore:=sapply(feature_type_name,function(x){
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
plot(density(resCpG.Genes$EnsRegScore))
# Regulatory weight:
resCpG.Genes[,RegWeight:=(0.4+1.6*((TypeScore+EnsRegScore)/2))]

plot(density(resCpG.Genes[!duplicated(locisID)]$RegWeight))



# > 3) Links Weight[0.5-1] = 0.5+0.5*max(LinkScore [0-1])
#   - for CpG associated to a gene based on TSS proximity  
resCpG.Genes[in.eQTR==FALSE,LinkScore:=sapply(abs(tss_dist),function(x){
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
plot(density(resCpG.Genes[in.eQTR==FALSE]$LinkScore))

# - for CpG associated to a gene based on eQTL studies: 
#de base, linksscore =1

resCpG.Genes[in.eQTR==TRUE,LinkScore:=1] 

#but, filter for low signif asso, and add bonus if highly signif eQTR and close to best eQTL
plot(density(na.omit(resCpG.Genes$eQTL_dist)))
resCpG.Genes[in.eQTR==TRUE,distScore:=sapply(abs(eQTL_dist),function(x){
  if(x<500){
    distScore<-1 
  }else {
    distScore<-1-0.5*x/2500
  }
  
  return(distScore)
})]
plot(density(na.omit(resCpG.Genes$avg.mlog10.pv.eQTLs)))

plot(density(na.omit(resCpG.Genes$avg.mlog10.pv.eQTLs)))
resCpG.Genes[in.eQTR==TRUE,RegScore:=avg.mlog10.pv.eQTLs*distScore] #integrate with RegScore : mean of pval of eQTL in the region

plot(density(na.omit(resCpG.Genes$RegScore))) 
#remove low signif links (<4) and add bonus if 
abline(v=4)
summary(na.omit(resCpG.Genes$RegScore))
resCpG.Genes<-resCpG.Genes[RegScore>4|in.eQTR==F]
#and add bonus ; +1*regScore/220

resCpG.Genes[in.eQTR==TRUE,LinkScore:=LinkScore+1*RegScore/max(RegScore)]
plot(density(resCpG.Genes$LinkScore))
resCpG.Genes[LinkScore>1.5]
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.001   5.056   7.971  12.019  13.895 221.354 

#because of some tss_distance dismatch between eQTRlinked and not, there is different linkScore
resCpG.Genes[,nEQTR.CpgGeneLink:=.N,by=c("locisID","gene")]
plot(density(resCpG.Genes$nEQTR.CpgGeneLink))
resCpG.Genes[nEQTR.CpgGeneLink>1]
#...we select the best linkScore by cpg-gene pair
resCpG.Genes<-resCpG.Genes[,isBestLinks:=LinkScore==max(LinkScore),by=c("locisID","gene")][isBestLinks==TRUE]
resCpG.Genes<-unique(resCpG.Genes[,-c("isClosest","isBestLinks")])
#for same linksWeught
resCpG.Genes<-unique(resCpG.Genes,by=c("locisID","gene"))
#and finally calculate the linksWeight
resCpG.Genes[,LinksWeight:=0.3+1.2*(LinkScore/2)]
plot(density(resCpG.Genes$LinksWeight))





# > finally :
resCpG.Genes[,CpGScore:=DMCScore*RegWeight*LinksWeight]

plot(density(resCpG.Genes$CpGScore))



# VALIDATION : 
lines(density(resCpG.Genes[pval<0.001]$CpGScore),col=2) #majority of significant CpG have a good score
lines(density(resCpG.Genes[meth.change>20]$CpGScore),col=3) #high meth.Change CpG are pull down 
lines(density(resCpG.Genes[pval<0.001& meth.change>20]$CpGScore),col=4) # but keep a good score if significant

nrow(resCpG.Genes[pval>0.01 & CpGScore>40]) # 29 CpG have a high CpGScore wereas pvalue >1%
resCpG.Genes[pval>0.01 & CpGScore>40][order(-CpGScore)] #because 1) in highly signif eQTR, high meth.change, pval close to 1%, in promoter, and close to gene

plot(density(resCpG.Genes[pval>0.01 ]$CpGScore)) 
abline(v=quantile(resCpG.Genes[pval>0.01 ]$CpGScore,0.999)) #99,9% non significant CpG have a CpG score < 103

nrow(resCpG.Genes[meth.change<20 & CpGScore>40]) # no CpG with small meth.change have a CpGScore >40
plot(density(resCpG.Genes[abs(meth.change)<20 ]$CpGScore)) 
abline(v=quantile(resCpG.Genes[abs(meth.change)<20 ]$CpGScore,0.999))#99,9% of small meth.change CpG have a CpG score < 55

# some validation plots : 
# my score highlight mostly CpG with high pvalue and meth.change : 
library(ggplot2)
library(patchwork)
m <- ggplot(resCpG.Genes, aes( y = CpGScore))

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
resCpG.Genes[,nCpG.Gene:=.N,by=.(gene)]
resCpG.Genes[CpGScore>40,nCpGSig.Gene:=.N,by=.(gene)]

ggplot(resCpG.Genes)+
  geom_bar(aes(x=nCpG.Gene))

ggplot(resCpG.Genes[CpGScore>40])+
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
resCpG.Genes[,coef:=0.5+0.5*sqrt(1/nCpG.Gene)]
plot(unique(resCpG.Genes$coef),unique(resCpG.Genes$nCpG.Gene))

resCpG.Genes[,GeneScore:=coef*CpGScore[which.max(abs(CpGScore))]+(1-coef)*median(CpGScore),by="gene"]

fwrite(resCpG.Genes,"analyses/withoutIUGR/2020-06-02_CF.LF_CpG_Gene_scorev2.csv",sep=";")

resGenes<-unique(resCpG.Genes[order(pval)],by="gene")[order(-GeneScore)]

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
resCpG.Genes[gene=="L2HGDH"] #1CpG avec score = 8 et un aute avec score =-8
median(resCpG.Genes[gene=="L2HGDH"]$CpGScore)
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
resCpG.GenesF<-resCpG.Genes[meth.change>0][order(pval)]
genesPval<-unique(resCpG.GenesF$gene)[1:1697]

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
resCpG.Genes[,ismaxCpGScore:=CpGScore==max(CpGScore),by=.(locisID,gene)]
resCpG.Genes<-resCpG.Genes[ismaxCpGScore==T][order(gene,-CpGScore)]
resCpG.Genes
resCpG.Genes[,GeneScore2:=CpGScore[which.max(abs(CpGScore))],by="gene"]
resCpG.Genes[nCpG.Gene>1,GeneScore2:=GeneScore2+mean(CpGScore[2:.N]),by="gene"]
resGenes2<-unique(resCpG.Genes,by="gene")[order(-GeneScore2)]

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

#numeric vector
geneList<-resGenes2$GeneScore2
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
                  pathway.id = "hsa00640", 
                  species    = "hsa",
                  limit      = list(gene=max(abs(geneList)), cpd=1),
                  kegg.dir ="analyses/withoutIUGR/pathwaysMap" 
)


#SEE GENE
library(stringr)
library(clusterProfiler)
resCpG.Genes<-fread("analyses/withoutIUGR/2020-05-26_CF.LF_CpG_Gene_score.csv",sep=";")
resCpG.Genes[gene=="SOD2"]


resCpG.Genes[str_detect(gene,"IRS")]
resGenes<-unique(resCpG.Genes,by="gene")[order(-GeneScore)]


#OTHER TEST OF CpG-score SUMMARIZING : 
#so i decided to alleviate the weigth of the number of CpG like this :


#geneSCORE3 : the sum

resCpG.Genes[,GeneScore3:=sum(CpGScore),by="gene"]

resGenes2<-unique(resCpG.Genes,by="gene")[order(-GeneScore3)]

resGenes2
resGenes2$gene[1:100]


plot(density(resGenes2$GeneScore3))


abline(v=60)
#OVER-REPRESENTATION TEST
#genes candidat : GeneScore2 >25
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
dotplot(rkGS,x=dfkGS$GeneScore.avg,showCategory=38)
emapplot(rkGS)

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

#diff with femame
setdiff(dfkGS$Description,dfkGSM$Description)


setdiff(dfkGSM$Description,dfkGS$Description)

