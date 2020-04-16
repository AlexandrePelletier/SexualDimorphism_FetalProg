#determiner comment utiliser GSEA pour nos données de methylation

#config et library
options(stringsAsFactors=F)
set.seed(12345)
library(limma)
source("scripts/scRNA-seq integration.R")
#output dir 
script_name <- "GSEA_test"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
date<-Sys.Date()
output <- file.path(outputDir,date)

#selon gsea-msigDB, pour que GSEA donne des résultats relavant pour des datas autres que RNA-seq, il est préférable de 
#1) ranker les genes "in order of most (largest value) to least (smallest value) "of interest" 
#2) si metric pas relevant biologiquement, vaux mieux juste mettre le rank en metrique

#> donc, nous chercons a faire un score permettant de trier les CpG selon leur chance d'etre rellement impactant sur l'expression du gène relié.
#pour cela, un cpg est le plus susceptible d'influencer un gène si :
#1) il a un FC elevé : 
#2) la difference avec les CTRL est très significative 
#3) il est sur un prom ou un enh(4, 6)
#4) il est sur un zone définit comme régulatrice dans ENSEMBL(ou UCSC ?)
#5) il est a proximité d'un motifTF
#6) son changement de niveau de méhylation est  cross corrélé avec d'autres CpG 


#=> ça nous donne un score par CpG (CpGScore) susceptible d'influencé tel ou tel gène.
#Aves ces scores, nous créons ensuite un score par gène :
#1) en faisant la somme des CpGscores lié au gène* normalisé par le nb totale de CpG lié à ce gène
#2) pour valoriser les gènes réellement exprimé dans les HSPC** (i.e dont la methylation est reelement susceptible d'impacter son expression),
#bonus1) il est dans une DMR
#bonus2) si le genes expr est tres différents entre les pop ?

#ce score est ensuite amoindrit si le gène n'est pas du tout retrouvé exprimés dans nos datas RNA-seq de référence**

#*si améliore par la suite le lien CpG-gene, en liant les CpG a plusieurs gènes par esemble (pas seulement le plus proche), on pourra facielement donné
#**  datas RNA-seq de HSPC publique, et nos datas (sc)RNAseq de  CTRL et LGA.

#le FC, Pval, si sur prom ou sur enh,  la distance du TSS, 
#with clusterProf, 

fit2<-readRDS("analyses/main/fit2_model13.rds")
colnames(fit2$coefficients)
res<- topTable(fit2,coef="FC.FL",n =Inf)
locis<-rownames(res)
res<-data.frame(row.names = locis,
                FC=res$logFC,
                pval=res$P.Value,
                pval.adj=res$adj.P.Val,
                annot[locis,c("chr","start","posAvant","gene","ENTREZID","type")])
#on enleve les locis non rattaché a des genes
sum(is.na(res$gene)) #5971
res<-res[!is.na(res$gene),]

#get ENSEMBL info : regulatory domain, et motifTF
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ensembldb")
# browseVignettes("ensembldb")
# #need db 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("EnsDb.Hsapiens.v86")
# 
# library(EnsDb.Hsapiens.v86)
# ## Making a "short cut"
# edb <- EnsDb.Hsapiens.v86
# ## print some informations for this package
# edb
# #get cross correl weigth.

#dont find in the docs how query chr and start and get annot regulatory,so with biomart :

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(biomaRt) 
#there is 1) Human Regulatory Features (hsapiens_regulatory_feature) and 2) Human Binding Motifs dataset (hsapiens_motif_feature) interesting (and 3: Human Regulatory Evidence (hsapiens_peak))
#1)  Human Regulatory Features
ensemblReg = useEnsembl(biomart="ENSEMBL_MART_FUNCGEN", dataset="hsapiens_regulatory_feature", GRCh=37,version = 95)

listFilters(ensemblReg) #chromosomal_region  e.g 1:100:10000000, 1:100000:200000
listAttributes(ensemblReg)
attrs<-c('chromosome_name','chromosome_start','chromosome_end','feature_type_name','feature_type_description','epigenome_name',
         'epigenome_description','activity','regulatory_stable_id')

#get info for 1 pos, par exemple le 1er locis :
head(res)
#               FC         pval   pval.adj   chr     start posAvant    gene ENTREZID type
# 844706  48.33741 2.379390e-08 0.01870860  chr7  36193351      515   EEPD1    80820    6
#need put in format 7:36193351:36193351
putCpGInChrRegFormat<-function(resCpG,colChr="chr",colStart="start"){
  return(apply(resCpG,1,function(x){
    chr<-as.numeric(strsplit(x[colChr],"r")[[1]][2])
    start<-as.numeric(x[colStart])
    return(paste(chr,start,start,sep= ":"))
    }
    ))
  
}
res$chrRegion<-putCpGInChrRegFormat(res)
head(res)
#so for the first locis
chrPos<-res$chrRegion[1]
result_BM <- getBM(attributes=attrs,
                   filters = 'chromosomal_region',
                   values = chrPos,
                   mart = ensemblReg)
result_BM
                                
tail(result_BM,100) #2 regions found with 2 type of features : promoters and CTCF binding site. seen in a lot of tissue/cells (epigenome_name)

#with a feature called enhancer (4) :
chrPos<-res$chrRegion[2]
result_BM <- getBM(attributes=attrs,
                   filters = 'chromosomal_region',
                   values = chrPos,
                   mart = ensemblReg)
result_BM #called promoter too, instead of enhancer, BUT, in activity col, all are denoted "active"

# to simply resutls dont need epigenome name in attr :
attrs<-c('chromosome_name','chromosome_start','chromosome_end','feature_type_name','activity','regulatory_stable_id')

#with much more locis
chrPos<-res$chrRegion[1:20]
result_BM <- getBM(attributes=attrs,
                   filters = 'chromosomal_region',
                   values = chrPos,
                   mart = ensemblReg,)
result_BM
#merge with our res
#pb, we loose the real pos of the CpG, we need need check wich cpg put chrname start 
#we would like a fnction 

res20<-head(res,20)
res20$locis<-rownames(res20)
res20.merge<-merge(res20,by=)
#with ecidence db


#with motif db

#conf annot with gene db

#get scRNAseq ecpr, need before to put res by genes :
head(res)
resGenes<-resLocisToGenes(res,withNCpGTot = T)
head(resGenes) #get also the number of CpG detected in at least 1 sample HPA2c library(nCpG), and the number of CpG tot detected by Msp1c ref

#then we can add scRNAseq expr : 
CTRL_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInCTRL.csv",row.names = 1)
LGA_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInLGA.csv",row.names = 1)
listExprMat<-list(CTRL=CTRL_expr_by_pop,LGA=LGA_expr_by_pop)
#add genes Expr in bulk scrnaseq and in subpop
resGenes<-find.scGeneExpr(resGenes,listExprMat)
head(resGenes)#keep only cols nCpG and Bulk_expr in ctrl and LGA
#add HSPC pub expr in resGENEs
resGenes$HSPC1<-annot[match(rownames(resGenes),annot$gene),"HSPC1"]
resGenes$HSPC2<-annot[match(rownames(resGenes),annot$gene),"HSPC2"]
head(resGenes)




                
                


