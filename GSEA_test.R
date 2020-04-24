#determiner comment utiliser GSEA pour nos données de methylation

#config et library
options(stringsAsFactors=F)
set.seed(12345)
library(limma)
source("scripts/scRNA-seq integration.R")
source("scripts/utils.R")
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

annot<-read.csv2("../../ref/annotation_CpG_HELP_ALL_070420.csv",row.names = 1)
head(annot)

#annot res :
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

res$chrRegion<-putInChrRegFormat(res)
head(res)
write.csv2(res,paste0(output,"_resCF.LF_AllLocis_with_AnnotNecessForCpGScore.csv"))
res<-read.csv2("analyses/GSEA_test/2020-04-16_resCF.LF_AllLocis_with_AnnotNecessForCpGScore.csv",row.names = 1)
head(res)
nrow(res)
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
                   mart = ensemblReg)
result_BM
result_BM$chrReg<-putInChrRegFormat(df=result_BM,
                  colChr = "chromosome_name",colStart ="chromosome_start",
                  colEnd = "chromosome_end",chrWithchr = F)

#merge with our res
#pb, we loose the real pos of the CpG, we need need check wich cpg put chrname start 
#we would like a fnction, for each chrRange in resBM, get locis in this region, and annotate res in these locis
  
res20<-head(res,20)
res20<-annotCpGResWithBMRes(res20,result_BM)
head(res20)

#with all res
#res<-annotCpGResWithBMRes(res) #take too much time, need to make step by step, chr by chr and split by 1k locis
table(res$chr)
#for efficiency, need only chrRegion and chr 
res2<-res[,c("type","chrRegion")]
head(res2)

for(chr in"chr1"){
  print(chr)
  n<-ceiling(table(res$chr)["chr1"]/1000)
  print(paste("partionner en ",n,"fois"))
  q2<-0
  for(i in 1){
    f<-round(i/n,2)
    q1<-quantile(res$start,f,na.rm = T)
    locis<-rownames(res)[res$chr==chr&res$start<=q1&res$start>q2]
    print(paste(chr,"part",i,"annotation de ",length(locis),"locis.."))
    subRes<-annotCpGResWithBMRes(res2[locis,])
    q2<-q1
    write.csv2(subRes,file = paste0(output,"annotEnsembl",chr,"part",i,".csv"),na = "")
    print("")
    print("..annotes et enregistrer")
  }
  
}
head(subRes,10)
subRes$RangeProm<-sapply(subRes$Promoter,tr,1) #work

for(chr in paste0("chr",9:21)){
  print(chr)
  n<-ceiling(table(res$chr)[chr]/5000)
  print(paste("partionner",n,"fois"))
  q2<-0
  for(i in 1:n){
    f<-round(i/n,2)
    q1<-quantile(res[res$chr==chr,"start"],f)
    locis<-rownames(res)[res$chr==chr&res$start<=q1&res$start>q2]
    if(length(locis)>0){
      print(paste(chr,"part",i,", annotation de",length(locis),"locis.."))
      subRes<-annotCpGResWithBMRes(res2[locis,])
      q2<-q1
      write.csv2(subRes,file = paste0(output,"annotEnsembl",chr,"part",i,".csv"))
      print("..annotes et enregistrer")
      
    }
    
  }
  
}

#add this annot in ref/annot
annot<-read.csv2("../../ref/annotation_CpG_HELP_ALL_070420.csv",row.names = 1)
head(annot)
chrs<-paste0("chr",1:21)
for(chr in chrs){
  print(chr)
  for(file in list.files("analyses/GSEA_test/",pattern = paste0(chr,"part"))){
    print(strsplit(file,chr)[[1]][2])
    ensAnno<-read.csv2(paste0("analyses/GSEA_test/",file),row.names = 1)
    locis<-rownames(ensAnno)
    annot[locis,"feature_type_name"]<-ensAnno$feature_type_name
    annot[locis,"CTCF.Binding.Site"]<-ensAnno$CTCF.Binding.Site
    annot[locis,"Promoter"]<-ensAnno$Promoter
    annot[locis,"Promoter.Flanking.Region"]<-ensAnno$Promoter.Flanking.Region
    annot[locis,"Open.chromatin"]<-ensAnno$Open.chromatin
    annot[locis,"Enhancer"]<-ensAnno$Enhancer
    annot[locis,"TF.binding.site"]<-ensAnno$TF.binding.site
  }
}

head(annot[which(annot$feature_type_name!=""),])
annot[which(annot$feature_type_name==""),]<-NA
head(annot)

#pendant qu'on y est, on ajoute annot bulk CTRL et LGA
CTRL_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInCTRL.csv",row.names = 1)
head(CTRL_expr_by_pop,20)
LGA_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInLGA.csv",row.names = 1)
for(gene in na.omit(unique(annot$gene))){
  if(gene%in%rownames(CTRL_expr_by_pop)){
    annot[which(annot$gene==gene),"CTRL_sc"]<-CTRL_expr_by_pop[gene,"bulk"]
  }
  if(gene%in%rownames(LGA_expr_by_pop)){
    annot[which(annot$gene==gene),"LGA_sc"]<-LGA_expr_by_pop[gene,"bulk"]
  }
  
}
head(annot[which(!is.na(annot$CTRL_sc)),])
write.csv2(annot,file = "../../ref/annotation_CpG_HELP_ALL_070420.csv",row.names = T)



#number of ensFeatures by type
#0) distrib type in 780k locis :
table(annot[rownames(data_F),"type"])
#distrib ensFeat by type :
allFeats<-list()
allFeatsCounts<-data.frame()

for(type in 0:6){
  
  #1) count by type
  allFeats[[as.character(type)]]<-unlist(lapply(res[which(res$type==type),"feature_type_name"],tr))
  print("ok")
  counts<-table(allFeats[[as.character(type)]])
  print("ok2")
  allFeatsCounts[names(counts),as.character(type)]<-counts
  allFeatsCounts["None",as.character(type)]<-sum(is.na(res[which(res$type==type),"feature_type_name"]))
  
}
allFeatsCounts

allFeatsCountsDf<-data.frame()
i<-1
for(type in colnames(allFeatsCounts)){
  for(ensFeat in rownames(allFeatsCounts)){
    allFeatsCountsDf[i,c("type","ensembFeat")]<-c(type,ensFeat)
    allFeatsCountsDf[i,"count"]<-allFeatsCounts[ensFeat,type]
    i<-i+1
  }
}
head(allFeatsCountsDf,10)
#barplot
library(ggplot2)
ggplot(allFeatsCountsDf)+
  geom_bar(aes(x = type, y = count,fill=ensembFeat), stat = 'identity')


#add in res
res<-data.frame(row.names = rownames(res),res,annot[rownames(res),12:18])
head(res)
#check if true
#type 4 and 6 enrichit en prom/enh features ?
m<-which(res$type%in%0:6 &res$feature_type_name!="")
length(m)#396k/780k
table(res$type[m],res$feature_type_name[m])

plot(factor(res$type[m]),factor(res$feature_type_name[m]))

#make df count by ensemblfeatures and type
library(dplyr)
featuresCounts<-res%>%count(type,feature_type_name,sort = T)
head(featuresCounts,20)
featuresCounts[featuresCounts$type==4,]#2k enh et 3k6 TF /150k
sum(res$type==4)
featuresCounts[featuresCounts$type==6,]#1k enh et no TF /200k
sum(res$type==6)
featuresCounts[featuresCounts$type==5,]#2k3 enh et noTF  /51k
sum(res$type==5)

featuresCounts[featuresCounts$type==0,]#8k enh et 1k TF /164k
sum(res$type==0)
featuresCounts[featuresCounts$type==1,]#3k enh et 500 TF /73k
sum(res$type==1)
featuresCounts[featuresCounts$type==2,]#3k5 enh et 500 TF /66k
sum(res$type==2)
featuresCounts[featuresCounts$type==3,]#3k enh et 500 TF /70k
sum(res$type==3)



#now, we can create our CpGscore : 
#we want FC=pvalue > 4/6 > regulatory element > TF > wieght of correlation > RNAseq reference
#so principal weight is FC[0:1] + pval[0:1] (max=2) (the DMC weigtht)
#We multiply this weight by the "regulatory weight" [0-1]
# +1 if (4,6), 0.5 if 5, 0.1 otherwise
# +0.7 if(CTCF,enh,prom,)+ 0.3 if prom flanking, openchrine, #+0.5 if(TFbind site), or 0 otherwise
#rank this and norm the rank to obtain a score between 0 and 1
 
##DMC weigtht
#(FC[0:1] : toavoid overfitting, we use the rank of abs(FC) as metriks normalized
res$FCScore<-rank(abs(res$FC))/nrow(res)
head(res)
tail(res)
#pval[0:1] : the same (with rank), to have the same weight
res$PvalScore<-rank(1/res$pval)/nrow(res)
res$DMCWeight<-res$PvalScore+res$FCScore

##regulatory weight
# +1 if (4,6), 0.5 if 5, 0.1 otherwise
res$TypeScore<-sapply(res$type,function(x){
  if(x%in%c(4,6)){
    score<-1
  }else if(x==5){
    score<-0.5
  }else{
    score<-0.1
  }
  return(score)
})

# +0.7 if(CTCF,enh,prom,)+ 0.3 if prom flanking, openchrine, #+0.5 if(TFbind site), or 0 otherwise
res$EnsRegScore<-sapply(res$feature_type_name,function(x){
  vecX<-tr(x)
  if(any(c("CTCF Binding Site","Promoter","Enhancer")%in%vecX)){
    score<-0.7
  }else if(any(c("Open chromatin","Promoter Flanking Region")%in%vecX)){
    score<-0.3
  }else{
    score<-0
  }
  if("TF binding site"%in%vecX){
    score<-score+0.5
  }
  return(score)
})

head(res)
tail(res)
#rank this and norm the rank to obtain a weight between 0 and 1
res$RegWeight<-rank(res$TypeScore+res$EnsRegScore)/max(rank(res$TypeScore+res$EnsRegScore))

##CpG score: DMCweight*RegWeight
res$CpGScore<-res$DMCWeight*res$RegWeight

#now we can calculate de GenesScores by sum(CpGScore in a gene)/nCpGtot in the genes and add a weight of this new score depending of the gene expr: 

#first, det the res by CpG before to put res by genes :
head(res)
resGenes<-resLocisToGenes(res,withNCpGTot = T)
head(resGenes) #get also the number CpG detected in at least 1 sample HPA2c library(nCpG), compared the number of CpG tot detected by Msp1c ref
for(gene in rownames(resGenes)){
  resGenes[gene,"GeneScore"]<-sum(res$CpGScore[which(res$gene==gene)])
}
head(resGenes)
#GSEA without weitghed for non expr genes :
library(clusterProfiler)
KeggGS <- read.gmt("../../ref/c2.cp.kegg.v7.1.symbols.gmt") 
geneList<-makeGeneList(resGenes,score = "GeneScore",returnRank = T)
head(geneList,50)

kk<-GSEA(geneList,TERM2GENE=KeggGS,nPerm = 1000,pvalueCutoff = 0.1)
head(kk,20) #done
emapplot(kk)

#normalized with the number of CpG in data_F by genes
for(gene in rownames(resGenes)){
  resGenes[gene,"GeneScore2"]<-sum(res$CpGScore[which(res$gene==gene)])/resGenes[gene,"nCpG"]
}
head(resGenes)
#GSEA without weitghed for non expr genes :

geneList<-makeGeneList(resGenes,score = "GeneScore2",returnRank = T)
head(geneList,50)

kk2<-GSEA(geneList,TERM2GENE=KeggGS,nPerm = 1000,pvalueCutoff = 0.1)
head(kk2,20) #done
emapplot(kk2)


#normalized with the number of CpG detected by msp1 by genes
for(gene in rownames(resGenes)){
  resGenes[gene,"GeneScore3"]<-sum(res$CpGScore[which(res$gene==gene)])/resGenes[gene,"nCpGtot"]
}
head(resGenes)
#GSEA without weitghed for non expr genes :

geneList<-makeGeneList(resGenes,score = "GeneScore3",returnRank = T)
head(geneList,50)

kk3<-GSEA(geneList,TERM2GENE=KeggGS,nPerm = 1000,pvalueCutoff = 0.1)
head(kk3,20) #done
emapplot(kk3)

saveRDS(list(resLocis=res,resGenes=resGenes),paste0(output,"resCLF.rds"))

#* confWeight [0-1] : TSS, ewas links, cross correl
#rank normalisé de sum de : 
#tss +1 : 
  #pour promcpg, <1000 => max point, descend progress au dela (0.9 <2000, 0.8<3000...)
  #pour les autres : <3000 max point, et descend progressivement apres

#ewas link : si locis +/-xpb rattaché à expression genes => +1

#cross correl : si cpg bouge pareillement a d'autres cpg a travers ech sans etre corrélé a group_complecity => poids 0>1


#need better norm than just divide by number of CpG tot, that too influencing
#moyenne pondéré /2: pour eviter inflaté genes avec beaucoup de locis, on veut ponderé le score des CpG :
#un gene est considéré comme tres impacté si il a 3 CpG important d'activer : (par exmple 1 sur enh, prom, genebody)
#donc, on notre moyenne pondéré favorisera surout les 3 topCpG score par gene : 
  #ScoreGene=CpG1+CpG2+CpG3 + 0.5*CpG4 + 0.2*CpG5 ... / (sommeCoeffTot)

#mais nb de CpG important dépend de cb CpG ont été rattaché au gene à la base : 
#donc plutot que 3 disons qu'il y a x CpG important par gene, pour trouver x, nous devons prendre en compte : 
#1) nb de CpG par gene : si peu de cpg rattaché, il ya moins de chance d'y avoir des cpg importants quand qd il yen a bcp
#2) taile du genes : plus il est grand plus on s'attend a qu'il est plus de cpg pour le reguler

#1) x=3 pour gene avec nCpG dans lespace interquartile , 2/4 si <q25/>q75, 1/5 si <q5/>q95
#2) +1/-1 si taille gene <q25/>q75

#ebayes/methylGSA ?


#we can add RNAseq expr to improve the score: 
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




                
                


