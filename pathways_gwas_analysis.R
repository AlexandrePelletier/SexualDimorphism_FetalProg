#2020-06-19 
#PATHWAYS ANALYSIS
#Based on GeneScore, i want perform pathway and gwas analysis to determine master gene regulator for future biological validation.

#let's first focus on Female Genes : 
library(data.table)
library(stringr)
library(clusterProfiler)
library(pathview)
library(ggplot2)
resF<-fread("analyses/withoutIUGR/2020-06-03_CF.LF_CpG_Gene_scorev3.2.csv")
source("scripts/utils.R")

#lets defined threshold ~ best compromise between keeping genes with hi Pvalue and genes with weak.
resF[,noSigGene:=all(pval>10e-3),by="gene"]
resF[,SigGene:=any(pval<10e-5),by="gene"]
plot(density(resF$GeneScore))
lines(density(unique(resF[noSigGene==TRUE],by="gene")$GeneScore),col=2)

lines(density(unique(resF[SigGene==TRUE],by="gene")$GeneScore),col=3)
#5there is pva >10^-3 but with hi g
abline(v=200)
gF<-unique(resF[GeneScore>200]$gene)
length(gF) #1123

kF <- enrichKEGG(gene         = bitr(gF,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05)
dkF<-as.data.frame(kF)
nrow(dkF) #13
dkF
dkF<-data.table(dkF)
dkF[,GeneScore.avg:=mean(unique(resF,by="gene")$GeneScore[unique(resF,by="gene")$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dkF
dotplot(kF,x=dkF$GeneScore,showCategory=13)
emapplot(kF)
enrichplot::upsetplot(kF,n=13)
#3 things :  
#1) MAPK have 13 genes exclusif, 
#2) 6 genes are shared between HPI, plurisigpathway, cush syndr, gastric cancer, hippo sig path, and basal cell carcinoma
#3) sig pathway and MODY share 2 exclu genes.
#what are this genes and plot pathview with genescore ?
dkF[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "/"),by="ID"]
#1)
genes1<-setdiff(tr(dkF[ID=="hsa04010"]$genes),unique(unlist(lapply(dkF[ID!="hsa04010"]$genes,tr))))
genes1
# [1] "MAPKAPK3" "DAXX"     "RPS6KA2"  "DUSP4"    "TRAF2"    "MAP3K8"   "IGF2"     "MAP3K11" 
# [9] "CACNB3"   "MAP2K5"   "TAOK2"    "NFATC3"   "NFATC1"   "MKNK2" 
#first make geneList of genescore
geneList<-unique(resF,by="gene")$GeneScore
names(geneList)<-unique(resF,by="gene")$gene
geneList<-sort(geneList,decreasing = T)
genes.df<-bitr(names(geneList),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneList[genes.df$SYMBOL]
geneList.Entrez<-genes.df$GeneScore
names(geneList.Entrez)<-genes.df$ENTREZID
pathview(gene.data  = geneList.Entrez,
         pathway.id = "hsa04010", 
         species    = "hsa",
         limit      = list(gene=max(abs(geneList.Entrez)), cpd=1),
         kegg.dir ="analyses/withoutIUGR/pathwaysMap")


#2020-06-30 ; trouver sex spe signif pathways /genes 
#I) by path enrichment of sex spé sig epigen affected genes 
#a) first, deter genescore cutoff x / deter EpigenAffGene (EAG) : 
source("scripts/2020-06-25_deterSexSpeMethGenes.R")
resM<-read.csv2("analyses/withoutIUGR/2020-04-16_res_locis_in_MC.ML_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",row.names = 1)
resM$locisID<-rownames(resM)
head(resM)
resF<-read.csv2("analyses/withoutIUGR/2020-04-16_res_locis_in_FC.FL_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",row.names = 1)
resF$locisID<-rownames(resF)
res_list<-list("C_M-L_M"=resM,"C_F-L_F"=resF)
names(res_list)
res_list<-deterEpigenAffGene(res_list,var_to_permut="Group_Sex",permut=1000,pvalThr=0.01,)


#b) deter sex spé sig genes   : 
#=> 100 permutation of all samples CTRL / LGA => LIMMA => pval and FC => gene score => genescore cutoff of x (for each save nb of genes > x)
#=> setdiff(genesF, genesM) and inversely => save genes spéF and M => all genes obs dans <5 permut => sex spé sig genes

#b) pathway enrichment of this sex spé genes

#II) by path enrichment of of sig affected genes
