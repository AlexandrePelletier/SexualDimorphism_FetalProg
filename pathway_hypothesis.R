
#DISCOVER PATHWAY / BIOLOGICAL PROCESS AFFECTED IN FEMALE :
#based on meth F only, what are the pathways/BP/msigdb/gwas ?

library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
library(ggrepel)
source("scripts/utils/methyl_utils.R")
set.seed(1234)


resF<-fread("analyses/withoutIUGR/2020-07-03_resCF_LF.csv")
resM<-fread("analyses/withoutIUGR/2020-07-03_resCM_LM.csv")


resF<-CalcGeneScore(resF,calcCpGScore = T)

ggplot(unique(resF[order(-GeneScore,pval)],by="gene"))+geom_point(aes(nCpG.Gene,GeneScore))

resM<-CalcGeneScore(resM,calcCpGScore = T)
ggplot(unique(resM[order(-GeneScore,pval)],by="gene"))+geom_point(aes(nCpG.Gene,GeneScore))


fwrite(resF,"analyses/withoutIUGR/2020-07-30_resCF_LF.csv")
fwrite(resM,"analyses/withoutIUGR/2020-07-30_resCM_LM.csv")

resF[,compa:="ctrlF_vs_lgaF"]
resM[,compa:="ctrlM_vs_lgaM"]
resA<-rbind(resF,resM)
resA[,sex:=ifelse(compa=="ctrlF_vs_lgaF","female","male")]
# Methylome change in LGA (compared to Control) in female and in male

p<-ggplot(unique(resA,by=c("locisID","compa")))+
  geom_point(aes(x=meth.change,y=-log10(pval),col=abs(meth.change)>20&pval<10^-3))+
  facet_wrap("compa")+
  theme_classic()

p + scale_color_manual(values = c("grey2","red")) + theme_minimal() +theme(legend.position = "none")

#Gene Score :
plot(density(unique(resF,by="gene")$GeneScore))
plot(density(unique(resM,by="gene")$GeneScore))

abline(v=90)
ggplot(unique(resA[order(pval)],by=c("gene","sex")),aes(x = sex,y=GeneScore,col=sex))+
  geom_boxplot()+ theme_minimal() +theme(legend.position = "none")

 # p<-p+ geom_label_repel(aes(label = ifelse(GeneScore>100&pval<10^-5,gene,"")),
 #                   box.padding   = 0.35, 
 #                   point.padding = 0.5,
 #                   segment.color = 'grey50')

p + scale_color_manual(values = c("grey","red")) + theme_minimal()

unique(resF,by="gene")[GeneScore>90&pval1000perm<0.01] #1554
unique(resM,by="gene")[GeneScore>90&pval1000perm<0.01] #120

#make genelist
#for female
geneListF<-unique(resF,by="gene")$GeneScore
names(geneListF)<-unique(resF,by="gene")$gene
geneListF<-sort(geneListF,decreasing = T)
head(geneListF,20)
genes.df<-bitr(names(geneListF),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneListF[genes.df$SYMBOL]
geneListF.Entrez<-genes.df$GeneScore
names(geneListF.Entrez)<-genes.df$ENTREZID
head(geneListF.Entrez)

#for male
geneListM<-unique(resM,by="gene")$GeneScore
names(geneListM)<-unique(resM,by="gene")$gene
geneListM<-sort(geneListM,decreasing = T)
head(geneListM,20)
genes.df<-bitr(names(geneListM),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneListM[genes.df$SYMBOL]
geneListM.Entrez<-genes.df$GeneScore
names(geneListM.Entrez)<-genes.df$ENTREZID
head(geneListM.Entrez)


#genes epigenetically impacted in LGA = queue de distrib GeneScore et pval1000perm<0.01

genesF<-unique(resF,by="gene")[GeneScore>95&pval1000perm<=0.01]$gene
length(genesF)#1336

genesM<-unique(resM,by="gene")[GeneScore>95&pval1000perm<=0.01]$gene
length(genesM)#93

genesF2<-unique(resF,by="gene")[order(-GeneScore)][pval1000perm<=0.005]$gene[1:500]
genesM2<-unique(resM,by="gene")[order(-GeneScore)][pval1000perm<=0.005]$gene[1:500]
length(genesM2)
#KEGG  1
candidat_genes.list<-list(female=genesF,male=genesM)

resFM_KEGG <- compareCluster(gene         = lapply(candidat_genes.list,function(genes)bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID),
                              fun = "enrichKEGG",
                              organism="hsa",
                             minGSSize    = 50,
                             maxGSSize    = 400,
                              pvalueCutoff=0.1)

nrow(data.frame(resFM_KEGG)) #54
dotplot(resFM_KEGG,showCategory=58)

emapplot(resFM_KEGG,pie="count",showCategory = 58)
resF_KEGG <- enrichKEGG(gene         = bitr(genesF,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)
nrow(data.frame(resF_KEGG)) #33
dotplot(resF_KEGG,showCategory=15)
emapplot(resF_KEGG,showCategory=26)

#KEGG  2
candidat_genes.list<-list(female_GS60=genesF,female_top500=genesF2,male_GS60=genesM,male_top500=genesM2[!is.na(genesM2)])

resFM_KEGG <- compareCluster(gene         = lapply(candidat_genes.list,function(genes)bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID),
                             fun = "enrichKEGG",
                             organism="hsa",
                             pvalueCutoff=0.1)

nrow(data.frame(resFM_KEGG)) #96
dotplot(resFM_KEGG,showCategory=52)

emapplot(resFM_KEGG,pie="count",showCategory = 30)



resF_KEGG <- enrichKEGG(gene         = bitr(genesF,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.2)
nrow(data.frame(resF_KEGG)) #62
dkF<-as.data.frame(resF_KEGG)
nrow(dkF) #13
dkF<-data.table(dkF)
dkF[,GeneScore.avg:=mean(unique(resF,by="gene")$GeneScore[unique(resF,by="gene")$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dotplot(resF_KEGG,x=dkF$GeneScore,showCategory=62)



dotplot(resF_KEGG,showCategory=62)
emapplot(resF_KEGG,showCategory=62)

pathwayofInterest<-c("Maturity onset diabetes of the young",
                     "Type II diabetes mellitus",
                     "Cushing syndrome",
                     "Gastric cancer",
                     "TGF-beta signaling pathway",
                     "Transcriptional misregulation in cancer",
                     "MAPK signaling pathway")



#try with metascape : 
fwrite(resF[GeneScore>115&pval1000perm<0.01,.(gene,GeneScore)],"analyses/withoutIUGR/genesLGAF_GeneScore115_pvalPerm0.01.csv",sep = ",")
fwrite(resM[GeneScore>115&pval1000perm<0.01,.(gene,GeneScore)],"analyses/withoutIUGR/genesLGAM_GeneScore115_pvalPerm0.01.csv",sep = ",")
#in male, all genes in the queues
plot(density(unique(resM,by="gene")$GeneScore))
abline(v=50)
genesM2<-unique(resM[GeneScore>50&pval1000perm<0.01]$gene)
length(genesM2)#863
fwrite(resM[GeneScore>50&pval1000perm<0.01,.(gene,GeneScore)],"analyses/withoutIUGR/genesLGAM_GeneScore50_pvalPerm0.01.csv",sep = ",")

pathwayofInterest<-c(pathwayofInterest,
                     "HIF-1 signaling pathway")
pathwayofInterest<-c(pathwayofInterest,
                     "Hepatocellular carcinoma",
                     "Insulin signaling pathway")
library("pathview")
#pathview of interesting pathway
#need pathway id
pathIDofInterest<-data.table(Reduce(rbind,list(data.frame(resF_KEGG),data.frame(resM_KEGG),data.frame(resM2_KEGG))))[Description%in%pathwayofInterest]$ID

pathF<-lapply(pathIDofInterest, function(id)return(pathview(gene.data  = geneListF.Entrez,
                                                     pathway.id = id, 
                                                     species    = "hsa",
                                                     limit      = list(gene=max(abs(geneListF)), cpd=1))))

pathM<-lapply(pathIDofInterest, function(id)return(pathview(gene.data  = geneListM.Entrez,
                                                            pathway.id = id, 
                                                            species    = "hsa",
                                                            limit      = list(gene=max(abs(geneListF)), cpd=1))))



#GSEA 
resF_KEGG.GSEA<- gseKEGG(geneList     = rank(geneListF.Entrez),
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 0.001,
                        verbose = FALSE)
nrow(as.data.frame(resF_KEGG.GSEA))#77

dotplot(resF_KEGG.GSEA,showCategory=77)
dkF<-as.data.frame(resF_KEGG.GSEA)
nrow(dkF) #13
dkF
dkF<-data.table(dkF)
dkF[,GeneScore.avg:=mean(unique(resF,by="gene")$GeneScore[unique(resF,by="gene")$gene %in% tr(core_enrichment,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dkF
dotplot(resF_KEGG.GSEA,x=dkF$GeneScore,showCategory=77)

resF_KEGG.GSEA
resF_KEGG.GSEA@result$Description
resF_KEGG.GSEA@result[resF_KEGG.GSEA@result$Description=="Longevity regulating pathway",]
emapplot(resF_KEGG.GSEA,showCategory=77,layout = "kk")

resM_KEGG.GSEA<- gseKEGG(geneList     = rank(geneListM.Entrez),
                         organism     = 'hsa', 
                         minGSSize    = 50,
                         pvalueCutoff = 0.001,
                         verbose = FALSE)
nrow(as.data.frame(resM_KEGG.GSEA))#42
dotplot(resM_KEGG.GSEA,showCategory=42)

emapplot(resM_KEGG.GSEA,showCategory=42)

resF_KEGG.GSEA@result[which(resF_KEGG.GSEA@result$Description=="Alzheimer disease"),]

resM_KEGG.GSEA@result[which(resM_KEGG.GSEA@result$Description=="Alzheimer disease"),]
resF_KEGG.GSEA@result

#compared with classical analysis (without genescore)

#for female
#i choose here to rank gene based on the Fold change of the most significant CpG in the gene
resF[,is.min:=pval==min(pval),by=gene]
resF.genes<-unique(resF[is.min==TRUE],by='gene')
#numeric vector
geneListF0<-resF.genes$meth.change
#named vector
names(geneListF0)<-resF.genes$gene
#sorted vector
geneListF0<-sort(geneListF0,decreasing = T)

head(geneListF0,20)


genes.df<-bitr(names(geneListF0),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$FC<-geneListF0[genes.df$SYMBOL]
geneListF0.Entrez<-genes.df$FC
names(geneListF0.Entrez)<-genes.df$ENTREZID
head(geneListF0.Entrez)


resF0_KEGG.GSEA<- gseKEGG(geneList     = rank(geneListF0.Entrez),
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 0.002,
                        verbose = FALSE)
nrow(as.data.frame(resF0_KEGG.GSEA)) #11

dotplot(resF0_KEGG.GSEA,showCategory=11)

emapplot(resF0_KEGG.GSEA,showCategory=11)


#for male
resM[,is.min:=pval==min(pval),by=gene]
resM.genes<-unique(resM[is.min==TRUE],by='gene')
#numeric vector
geneListM0<-resM.genes$meth.change
#named vector
names(geneListM0)<-resM.genes$gene
#sorted vector
geneListM0<-sort(geneListM0,decreasing = T)

head(geneListM0,20)


genes.df<-bitr(names(geneListM0),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$FC<-geneListM0[genes.df$SYMBOL]
geneListM0.Entrez<-genes.df$FC
names(geneListM0.Entrez)<-genes.df$ENTREZID
head(geneListM0.Entrez)


resM0_KEGG.GSEA<- gseKEGG(geneList     = rank(geneListM0.Entrez),
                          organism     = 'hsa', 
                          minGSSize    = 50,
                          pvalueCutoff = 0.002,
                          verbose = FALSE)
nrow(as.data.frame(resM0_KEGG.GSEA)) #38
dotplot(resM0_KEGG.GSEA,showCategory=38)
pM0<-resM0_KEGG.GSEA@result$Description
pM1<-resM_KEGG.GSEA@result$Description
length(intersect(pM0,pM1))
length(pM1)
emapplot(resM0_KEGG.GSEA,showCategory=11)

data.frame(nP)



#genes longevity
longevityPaths<-c("AMPK signaling pathway","mTOR signaling pathway","Longevity regulating pathway",
                  "PI3K-Akt signaling pathway","Wnt signaling pathway","FoxO signaling pathway","Insulin signaling pathway",
                  "Ras signaling pathway","Rap1 signaling pathway","Regulation of lipolysis in adipocytes")
longevityPathsID<-unique(data.table(Reduce(rbind,list(data.table(as.data.frame(resF_KEGG.GSEA))[,.(ID,Description)],
                                                      data.table(as.data.frame(resFM_KEGG))[,.(ID,Description)])))[Description%in%longevityPaths]$ID)

longevityGenes_list<-lapply(longevityPathsID,getGenesKEGGPathw)

names(longevityGenes_list)<-longevityPaths
longevityGenes<-unique(unlist(longevityGenes_list))

resA[,longevity:=gene%in%longevityGenes]
resA[,random:=F]
resA[,random:=gene%in%sample(unique(gene),length(longevityGenes))]
resA[,GeneSet:=NA]
resA[longevity==T,GeneSet:="longevity_related"]
resA[random==T,GeneSet:="random1"]
resA[gene%in%sample(unique(gene),length(longevityGenes)),GeneSet:="random2"]
resA[gene%in%sample(unique(gene),length(longevityGenes)),GeneSet:="random3"]
resA[gene%in%sample(unique(gene),length(longevityGenes)),GeneSet:="random4"]
resA[gene%in%sample(unique(gene),length(longevityGenes)),GeneSet:="random5"]

ggplot(unique(resA[!is.na(GeneSet)],by=c("sex","gene")))+geom_boxplot(aes(x=sex,y=GeneScore,fill=GeneSet))

#or
ggplot(unique(resA,by=c("sex","gene")))+
  geom_boxplot(aes(x=sex,y=GeneScore,fill=gene%in%longevityGenes),outlier.shape = NA)+
  coord_cartesian(ylim = c(-100,200))+
  theme_minimal()
length(longevityGenes)

unique(resF,by="gene")[!(gene%in%longevityGenes)]

t.test(unique(resF,by="gene")[gene%in%longevityGenes]$GeneScore,unique(resF,by="gene")[!(gene%in%longevityGenes)]$GeneScore)

t.test(unique(resM,by="gene")[gene%in%longevityGenes]$GeneScore,unique(resM,by="gene")[!(gene%in%longevityGenes)]$GeneScore)


#Pathway colored with GeneScore:
library(pathview)

pathF<-lapply(longevityPathsID, function(id)return(pathview(gene.data  = geneListF.Entrez,
                                                    pathway.id = id, 
                                                    species    = "hsa",
                                                    limit      = list(gene=max(abs(geneListF)), cpd=1))))

pathM<-lapply(longevityPathsID, function(id)return(pathview(gene.data  = geneListM.Entrez,
                                                    pathway.id = id, 
                                                    species    = "hsa",
                                                    limit      = list(gene=max(abs(geneListF)), cpd=1))))


#2020-07-30 figure pour pres
#DEG ARE ENRICHED IN GENESCORE
library(data.table)
library(ggplot2)
library(patchwork)
set.seed(1234)
options(stringsAsFactors = F)
source("scripts/utils/methyl_utils.R")

resF<-fread("analyses/withoutIUGR/2020-07-30_resCF_LF.csv")
resF<-CalcGeneScore(resF,calcCpGScore = T)


res_DEG_meth<-fread("../../../Alexandre_SC/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")

res_DEG_meth<-res_DEG_meth[,.SD,.SDcols=-"GeneScore"]

res_DEG_meth<-merge(res_DEG_meth[,1:7],resF,by=c("gene"),all.x=T)
res_DEG_meth[is.na(padj),padj:=1]


p1<-ggplot(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene"))+
  geom_boxplot(aes(x = padj<0.05,y=meth.change,col = padj<0.05),outlier.shape = NA)+
  scale_color_manual(values = c("grey2","red"))+
  coord_cartesian(ylim = c(-10,75))+theme_minimal()+theme(legend.position = "none")

p2<-ggplot(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene"))+
  geom_boxplot(aes(x = padj<0.05,y=-log10(pval),col = padj<0.05),outlier.shape = NA)+
  scale_color_manual(values = c("grey2","red"))+
  coord_cartesian(ylim = c(0,5))+theme_minimal()+theme(legend.position = "none")

p3<-ggplot(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene"))+
  geom_boxplot(aes(x = padj<0.05,y=GeneScore,col = padj<0.05),outlier.shape = NA)+
  scale_color_manual(values = c("grey2","red"))+
  coord_cartesian(ylim = c(-40,175))+theme_minimal()+theme(legend.position = "none")

p_all<-p1+p2+p3

p_all
t.test(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj<0.05]$GeneScore,unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj>0.05]$GeneScore)

t.test(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj<0.05]$meth.change,unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj>0.05]$meth.change)

t.test(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj<0.05]$pval,unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj>0.05]$pval)

#DEG are related to longevity pathways
stem_DEG<-c("FOSB","JUN","ID2","KLF2","EGR1","KLF4","FOS","HES1","JUNB","ZFP36","HES4","IER2")


ggplot(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene"),aes(x = log2FoldChange,y=-log10(pvalue),col=padj<0.05))+
  geom_point()+
  geom_label_repel(aes(label = ifelse(padj<=0.05&gene%in%stem_DEG,gene,"")),col="deepskyblue2",
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()

#longevity/nutr sens
ggplot(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene"),aes(x = log2FoldChange,y=-log10(pvalue),col=padj<0.05))+
  geom_point()+
  geom_label_repel(aes(label = ifelse(padj<=0.05&gene%in%longevityGenes,gene,"")),col="darkorange1",
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()

#hormonal
resF_KEGG.GSEA@result$Description
hormonalPaths<-c("GnRH secretion","Morphine addiction","Adrenergic signaling in cardiomyocytes","Oxytocin signaling pathway",
                  "Prolactin signaling pathway","Aldosterone synthesis and secretion","Cholinergic synapse","Parathyroid hormone synthesis, secretion and action",
                  "Endocrine resistance","Growth hormone synthesis, secretion and action","Insulin secretion","Cortisol synthesis and secretion")
hormonalPathsID<-unique(data.table(Reduce(rbind,list(data.table(as.data.frame(resF_KEGG.GSEA))[,.(ID,Description)],
                                                      data.table(as.data.frame(resFM_KEGG))[,.(ID,Description)])))[Description%in%hormonalPaths]$ID)
#add Estrogen signaling pathway", 
hormonalPathsID<-c(hormonalPathsID,"hsa04915")
hormonalPaths<-c(hormonalPaths,"Estrogen signaling pathway")

hormonalGenes_list<-lapply(hormonalPathsID,getGenesKEGGPathw)
names(hormonalGenes_list)<-hormonalPathsID
hormonalGenes<-unique(unlist(hormonalGenes_list))

ggplot(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene"),aes(x = log2FoldChange,y=-log10(pvalue),col=padj<0.05))+
  geom_point()+
  geom_label_repel(aes(label = ifelse(padj<=0.05&gene%in%hormonalGenes,gene,"")),col="chartreuse4",
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()


#DEG - GeneScore
plot(density(unique(resF,by="gene")$GeneScore))
unique(resF,by="gene")[GeneScore>50&pval1000perm<0.05]
ggplot(unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj<0.05],aes(x = log2FoldChange,y=GeneScore))+
  geom_point(aes(col=GeneScore>50&pval1000perm<0.05))+
  geom_label_repel(aes(label = ifelse(GeneScore>50&pval1000perm<0.05&gene%in%c(longevityGenes,stem_DEG,hormonalGenes),gene,"")),col="black",
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()



#KEGG on this DEG
DEG<-unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj<=0.05]$gene
length(DEG)
unique(res_DEG_meth[order(-GeneScore,pval)],by="gene")[padj>0.05]
resDEG_KEGG <- enrichKEGG(gene         = bitr(DEG,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,minGSSize = 30
                         )

nrow(data.frame(resDEG_KEGG)) #17
dotplot(resDEG_KEGG,showCategory=32)
emapplot(resDEG_KEGG,showCategory=32)





#ANNEXE /Test  [a trier]
#but a lot of genes in LGAF > more stringeant to detect Highly affected pathway
plot(density(unique(resF,by="gene")$GeneScore))
abline(v=150)
genesF2<-unique(resF[GeneScore>150&pval1000perm<0.01]$gene)
length(genesF2)#394
resF2_KEGG <- enrichKEGG(gene         = bitr(genesF2,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.2)
nrow(data.frame(resF2_KEGG)) #9
dotplot(resF2_KEGG,showCategory=9)
emapplot(resF2_KEGG,showCategory=9)

genesM3<-unique(resM[GeneScore>150&pval1000perm<0.01]$gene)
length(genesM3)#14
resM3_KEGG <- enrichKEGG(gene         = bitr(genesM3,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.2)
nrow(data.frame(resM3_KEGG)) #12
dotplot(resM3_KEGG,showCategory=9)
emapplot(resM3_KEGG,showCategory=9)
#the same pathways


#see also longevity/IGF1 pathway :
pathview(gene.data  = geneListF.Entrez,
         pathway.id = "hsa04211", 
         species    = "hsa",
         limit      = list(gene=max(abs(geneListF)), cpd=1))
pathview(gene.data  = geneListM.Entrez,
         pathway.id = "hsa04211", 
         species    = "hsa",
         limit      = list(gene=max(abs(geneListF)), cpd=1))


#MAPKAPK3 interessant car only methylated in LGAF
plotMeth(resF[gene=="MAPKAPK3"]$locisID[1],plot = "jitter")

#Phospho HSPB1, also only methylated in LGAF
plotMeth(resF[gene=="HSPB1"&CpGScore>20]$locisID)
#Phospho ZFP36, KD in LGA



#gene only in F, 

genesOF<-genesF[!genesF%in%genesM2]
length(genesOF) #967
resOF_KEGG <- enrichKEGG(gene         = bitr(genesOF,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05)
nrow(data.frame(resOF_KEGG)) #18
dotplot(resOF_KEGG,showCategory=18)
emapplot(resOF_KEGG,showCategory=18)



genesOF2<-genesF2[!(genesF2%in%genesM2)]
length(genesOF2) #339
resOF2_KEGG <- enrichKEGG(gene         = bitr(genesOF2,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05)
nrow(data.frame(resOF2_KEGG)) #2
dotplot(resOF2_KEGG,showCategory=34)
emapplot(resOF2_KEGG,showCategory=34)
#MODY with PDX1
tr(data.frame(resOF2_KEGG)$geneID,tradEntrezInSymbol = T)
plotMeth(resF[gene=="PDX1"&CpGScore>20]$locisID)

#gene only  in M
genesOM<-genesM[!genesM%in%genesF]
length(genesOM) #32
resOM_KEGG <- enrichKEGG(gene         = bitr(genesOM,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05)

nrow(data.frame(resOM_KEGG)) #1
dotplot(resOM_KEGG,showCategory=1)
emapplot(resOM_KEGG,showCategory=18)



genesOM2<-genesM2[!(genesM2%in%genesF)]
length(genesOM2) #713
resOM2_KEGG <- enrichKEGG(gene         = bitr(genesOM2,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05)
nrow(data.frame(resOM2_KEGG)) #5
dotplot(resOM2_KEGG,showCategory=34)
emapplot(resOM2_KEGG,showCategory=34)
#hippo 
hippoM_only<-tr(data.frame(resOM2_KEGG)$geneID[4],tradEntrezInSymbol = T)
resM[gene%in%hippoM_only&CpGScore>30]
plotMeth(resM[gene=="LLGL2"&CpGScore>20]$locisID)


#2020-07-23 : LGAF vs LGAM
options(stringsAsFactors = F)
source("scripts/utils/methyl_utils.R")
cpg.regs_ref<-fread("../../ref/2020-06-29_All_CpG-Gene_links.csv")
batchF<-prepBatchDf(batch,
                    varNumToModel = c("Mat.Age"),
                    varFacToModel=c("Group_Sex",'batch',"latino","Group_Complexity_Fac"))



res<-RunMethAnalysis(methyl_df = methyl_df,
                     batch_filtered = batchF,
                     formule = ~0 + Group_Sex  + batch  + latino + Mat.Age + Group_Complexity_Fac,
                     compas_df = data.frame(compa=c("Group_SexL_F-Group_SexL_M","Group_SexC_F-Group_SexC_M"),
                                            abbrev_compa=c("L_F-L_M","C_F-C_M")),
                     cpg.regs_ref =cpg.regs_ref )

resC<-res$`C_F-C_M`
resL<-res$`L_F-L_M`
rm(res)

library(ggplot2)
library(patchwork)

p1<-ggplot(unique(resC,by="gene"))+geom_point(aes(x=GeneScore,y=-log10(pval)))
p2<-ggplot(unique(resL,by="gene"))+geom_point(aes(x=GeneScore,y=-log10(pval)))
p1+p2

#KEGG

# genescore pos = ?
unique(resL,by="gene")[order(-GeneScore)] #NFIC
plotMeth(resL[gene=="WASH7P"&abs(CpGScore)>20]$locisID)
plotMeth(resL[gene=="KANSL1-AS1"&abs(CpGScore)>20]$locisID) #GeneScore pos = hypermeth chez M
plotMeth(resL[gene=="PLEC"&abs(CpGScore)>20]$locisID) #GeneScore neg = hypermeth chez F
plotMeth(resL[gene=="GNAS"&abs(CpGScore)>20]$locisID)

plotMeth(resF[gene=="SOCS3"&abs(CpGScore)>20]$locisID)

p1<-ggplot(unique(resL,by="gene"))+geom_point(aes(x=meth.change,y=-log10(pval)))+coord_cartesian(ylim = c(0,7))+ggtitle("LGAF vs LGAM")
p2<-ggplot(unique(resC,by="gene"))+geom_point(aes(x=meth.change,y=-log10(pval)))+coord_cartesian(ylim = c(0,7))+ggtitle("ctrlF vs ctrlM")
p1+p2
# hypermet chez LGA femelle, 
plot(density(unique(resL,by="gene")$GeneScore))
abline(v=-60)
genesLF<-unique(resL[GeneScore<(-60)]$gene)
length(genesLF)#521
resLF_KEGG <- enrichKEGG(gene         = bitr(genesLF,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

nrow(data.frame(resLF_KEGG)) #32
dotplot(resLF_KEGG,showCategory=32)
emapplot(resLF_KEGG,showCategory=32)

# hypermet chez LGA male, 
plot(density(unique(resL,by="gene")$GeneScore))
abline(v=60)
genesLM<-unique(resL[GeneScore>60]$gene)
length(genesLM)#100
resLM_KEGG <- enrichKEGG(gene         = bitr(genesLM,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

nrow(data.frame(resLM_KEGG)) #3
dotplot(resLM_KEGG,showCategory=32)
emapplot(resLM_KEGG,showCategory=32)

#interesting : hippo signal and signalling pathway regulating pluripotency appear hyperMet chez femelle and Male, 
#1) these pathways are not in ctrl ?
ggplot(unique(resC,by="gene"))+geom_point(aes(x=meth.change,y=-log10(pval)))


plot(density(unique(resC,by="gene")$GeneScore))
abline(v=60)
genesCM<-unique(resC[GeneScore>60]$gene)
length(genesCM)#879

resCM_KEGG <- enrichKEGG(gene         = bitr(genesCM,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

nrow(data.frame(resCM_KEGG)) #17
dotplot(resCM_KEGG,showCategory=32)
emapplot(resCM_KEGG,showCategory=32)

#2) which gene in these pathway are hypermet chez femelle, chez male ?

pathwayofInterest2<-c("Signaling pathways regulating pluripotency of stem cells",
                     "Hippo signaling pathway")

#pathview of interesting pathway
library("pathview")

#for LGA
geneListL<-unique(resL,by="gene")$GeneScore
names(geneListL)<-unique(resL,by="gene")$gene
geneListL<-sort(geneListL,decreasing = T)
head(geneListL,20)
genes.df<-bitr(names(geneListL),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneListL[genes.df$SYMBOL]
geneListL.Entrez<-genes.df$GeneScore
names(geneListL.Entrez)<-genes.df$ENTREZID
head(geneListL.Entrez)



#need pathway id
pathIDofInterest2<-unique(data.table(Reduce(rbind,list(data.frame(resLF_KEGG),data.frame(resCM_KEGG))))[Description%in%pathwayofInterest2]$ID)

pathL<-lapply(pathIDofInterest2, function(id)return(pathview(gene.data  = geneListL.Entrez,
                                                            pathway.id = id, 
                                                            species    = "hsa",
                                                            limit      = list(gene=max(abs(geneListL)), cpd=1))))

plotMeth(resL[gene=="FGFR1"&abs(CpGScore)>20]$locisID)
plotMeth(resL[gene=="IGF1R"&abs(CpGScore)>20]$locisID)
plotMeth(resL[gene=="SMAD3"&abs(CpGScore)>20]$locisID)
plotMeth(resL[gene=="SMAD1"&abs(CpGScore)>20]$locisID)
#make genelist
#for ctrl
geneListC<-unique(resC,by="gene")$GeneScore
names(geneListC)<-unique(resC,by="gene")$gene
geneListC<-sort(geneListC,decreasing = T)
head(geneListC,20)
genes.df<-bitr(names(geneListC),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneListC[genes.df$SYMBOL]
geneListC.Entrez<-genes.df$GeneScore
names(geneListC.Entrez)<-genes.df$ENTREZID
head(geneListC.Entrez)
pathC<-lapply(pathIDofInterest2, function(id)return(pathview(gene.data  = geneListC.Entrez,
                                                            pathway.id = id, 
                                                            species    = "hsa",
                                                            limit      = list(gene=max(abs(geneListL)), cpd=1))))


#Validate longevity pathways /nutrient sensing pathways upregulate in LGAF comapred to ctrlF: 


#1) create/find a gene list of genes involved nutrient sensing pathway :
#a quoi ressemble une custom ?
custom_exampl<-read.gmt("ref/msigdb/c2.cp.kegg.v7.1.symbols.gmt") 
head(custom_exampl)
custom_exampl[str_detect(custom_exampl$term,"LONGEVITY")] #no longevity in this db

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")
getGenesKEGGPathw<-function(pathID){
  library(KEGGREST)
  g<-keggGet(pathID)[[1]]$GENE
  g<-g[1:length(g)%%2==0]
  return(as.vector(sapply(g,function(x)strsplit(x,";")[[1]][1])))
}

longevity<-getGenesKEGGPathw("hsa04211")

resA<-merge(resF[,compa:="ctrlFvslgaF"],resC[,compa:="ctrlFvsctrlM"],all=T)

ggplot(resA[,longevityPathway:=gene%in%longevity])+geom_boxplot(aes(x=compa,y=GeneScore,fill=longevityPathway))


length(longevity)

insulin<-getGenesKEGGPathw("hsa04910")

diabetesII<-getGenesKEGGPathw("hsa04930")

cushing<-getGenesKEGGPathw("hsa04934")
  
insulinR<-getGenesKEGGPathw("hsa04931")

nutrientSens<-fread("../List_geneofinterest_nutrient_070720.csv",header = F)$V1
ggplot(resA[,nutrient_related:=gene%in%nutrientSens])+geom_boxplot(aes(x=compa,y=GeneScore,fill=nutrient_related))




#but lack ctrlF vs LGAM
resCFLM<-RunMethAnalysis(methyl_df = methyl_df,
                     batch_filtered = batchF,
                     formule = ~0 + Group_Sex  + batch  + latino + Mat.Age + Group_Complexity_Fac,
                     compas_df = data.frame(compa="Group_SexC_F-Group_SexL_M",
                                            abbrev_compa="C_F-L_M"),
                     cpg.regs_ref =cpg.regs_ref )

resA<-merge(resA,resCFLM[,compa:="ctrlFvslgaM"],all=T)
resA[,longevityPathway:=gene%in%longevity]
resA[,nutrient_related:=gene%in%nutrientSens]

ggplot(resA)+geom_boxplot(aes(x=compa,y=GeneScore,fill=longevityPathway))

ggplot(resA)+geom_boxplot(aes(x=compa,y=GeneScore,fill=nutrient_related))

ggplot(resA[gene%in%insulin])+geom_boxplot(aes(x=compa,y=GeneScore))

ggplot(resA[gene%in%diabetesII])+geom_boxplot(aes(x=compa,y=GeneScore))

ggplot(resA[gene%in%cushing])+geom_boxplot(aes(x=compa,y=GeneScore))
ggplot(resA[gene%in%insulinR])+geom_boxplot(aes(x=compa,y=GeneScore))


#dtGeneScore : is there gene specifically meth in LGAF ?
resG<-unique(resA,by=c('compa','gene'))[,dtGS_lgaF:=GeneScore[compa=="ctrlFvslgaF"]-GeneScore[compa=="ctrlFvslgaM"],by="gene"]
resG
plot(density(resG$dtGS_lgaF))
resG[dtGS_lgaF>100]$gene


#pathway spé femelle :
rk_dtgsF<-enrichKEGG(gene         = bitr(resG[dtGS_lgaF>100]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05)

nrow(data.frame(rk_dtgsF)) #1
dotplot(rk_dtgsF,showCategory=32)
emapplot(rk_dtgsF,showCategory=32)
library(pathview)


path<-pathview(gene.data  = geneListF.Entrez,
                             pathway.id ="hsa04928" , 
                             species    = "hsa",
                             limit      = list(gene=max(abs(geneListF)), cpd=1))
path<-pathview(gene.data  = geneListM.Entrez,
               pathway.id ="hsa04928" , 
               species    = "hsa",
               limit      = list(gene=max(abs(geneListF)), cpd=1))

ggplot(resA[gene%in%getGenesKEGGPathw("hsa04060")])+geom_boxplot(aes(x=compa,y=GeneScore))

ggplot(resA[gene%in%tr("1385/2778/9935/4929/860/6256",tradEntrezInSymbol = T)])+geom_boxplot(aes(x=compa,y=GeneScore))
tr("1385/2778/9935/4929/860/6256",tradEntrezInSymbol = T)

#STOP HRE
#pathway spe male :

rk_dtgsM<-enrichKEGG(gene         = bitr(resG[dtGS_lgaF>-100]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05)
nrow(data.frame(rk_dtgsM)) #1
dotplot(rk_dtgsM,showCategory=32)
emapplot(rk_dtgsM,showCategory=32)


#TO DO +++ : is there change really significative ? permut pour enrichKEGG, ttest pour ggplot 

#GSEA 


#on this pathway could we see gene expression dysregulation in LGAF ?

#=> pathways pertubed in LGAF, what are the hypothesis for this pertubance, and what are the possible biological consequences?

#how validate this hypothesis ?
