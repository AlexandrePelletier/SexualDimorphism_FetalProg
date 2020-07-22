
#DISCOVER PATHWAY / BIOLOGICAL PROCESS AFFECTED IN FEMALE :
#based on meth F only, what are the pathways/BP/msigdb/gwas ?
resF<-fread("analyses/withoutIUGR/2020-07-03_resCF_LF.csv")
resM<-fread("analyses/withoutIUGR/2020-07-03_resCM_LM.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

#Female : 
#genes epigenetically impacted in LGA = queue de distrib GeneScore et pval1000perm<0.01
plot(density(unique(resF,by="gene")$GeneScore))
abline(v=115)
genesF<-unique(resF[GeneScore>115&pval1000perm<0.01]$gene)
length(genesF)#1117




#try with metascape : 
fwrite(resF[GeneScore>115&pval1000perm<0.01,.(gene,GeneScore)],"analyses/withoutIUGR/genesLGAF_GeneScore115_pvalPerm0.01.csv",sep = ",")

#KEGG 

resF_KEGG <- enrichKEGG(gene         = bitr(genesF,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)
nrow(data.frame(resF_KEGG)) #26
dotplot(resF_KEGG,showCategory=26)
emapplot(resF_KEGG,showCategory=26)

pathwayofInterest<-c("Maturity onset diabetes of the young",
                     "Type II diabetes mellitus",
                     "Cushing syndrome",
                     "Gastric cancer",
                     "TGF-beta signaling pathway",
                     "Transcriptional misregulation in cancer",
                     "MAPK signaling pathway")


#in male, same filter

genesM<-unique(resM[GeneScore>115&pval1000perm<0.01]$gene)
length(genesM)#73
#try with metascape : 
fwrite(resM[GeneScore>115&pval1000perm<0.01,.(gene,GeneScore)],"analyses/withoutIUGR/genesLGAM_GeneScore115_pvalPerm0.01.csv",sep = ",")

resM_KEGG <- enrichKEGG(gene         = bitr(genesM,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.2)
nrow(data.frame(resM_KEGG)) #2
dotplot(resM_KEGG,showCategory=26)
emapplot(resM_KEGG,showCategory=26)


pathwayofInterest<-c(pathwayofInterest,
                     "HIF-1 signaling pathway")
#in male, all genes in the queues

plot(density(unique(resM,by="gene")$GeneScore))
abline(v=50)
genesM2<-unique(resM[GeneScore>50&pval1000perm<0.01]$gene)
length(genesM2)#863

#try with metascape : 
fwrite(resM[GeneScore>50&pval1000perm<0.01,.(gene,GeneScore)],"analyses/withoutIUGR/genesLGAM_GeneScore50_pvalPerm0.01.csv",sep = ",")


resM2_KEGG <- enrichKEGG(gene         = bitr(genesM2,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)
nrow(data.frame(resM2_KEGG)) #14
dotplot(resM2_KEGG,showCategory=26)
emapplot(resM2_KEGG,showCategory=26)


pathwayofInterest<-c(pathwayofInterest,
                     "Hepatocellular carcinoma",
                     "Insulin signaling pathway")

#pathview of interesting pathway
library("pathview")


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

#on this pathway could we see gene expression dysregulation in LGAF ?

#=> pathways pertubed in LGAF, what are the hypothesis for this pertubance, and what are the possible biological consequences?

#how validate this hypothesis ?
