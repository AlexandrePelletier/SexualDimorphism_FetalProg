
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
custom_exampl<-read.gmt("../../ref/msigdb/c2.cp.kegg.v7.1.symbols.gmt") 
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
