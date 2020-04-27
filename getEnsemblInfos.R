#get ENSEMBL info : regulatory domain, et motifTF
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


#for the 780k locis  after filtration, 
#to start get 1 results from 1 compa :
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


