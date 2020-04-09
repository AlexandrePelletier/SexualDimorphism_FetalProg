
#config et library
options(stringsAsFactors=F)
library(data.table)
library(stringr)
library(ggplot2)
set.seed(12345)
#output dir 
script_name <- "newComplexityVar_test"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
date<-Sys.Date()
output <- file.path(outputDir,date)

#get genes with no expr modif
library(Seurat)
samples.integrated<-readRDS("../../../../ULYSSE/analysis/CBP547_and_CBP552_integrated_Seurat_object.rds")
head(samples.integrated@meta.data)
tail(samples.integrated@meta.data)
Idents(samples.integrated)<-"orig.ident"
avg.expr<-log2(AverageExpression(samples.integrated, verbose = FALSE)$RNA)
#avg.expr$gene <- rownames(avg.expr)
head(avg.expr)
plot(density(avg.expr$freshCD34_CBP547_CTRL))
#remove genes with log2expr <0.1 in at least one sample
avg.expr<-na.omit(avg.expr[rowSums(avg.expr<0.5|avg.expr==Inf)==0,])
nrow(avg.expr) #931
plot(avg.expr$freshCD34_CBP547_CTRL,avg.expr$LGA_3k)
#keep only genes +/- sd from the regression
summary(lm(avg.expr$freshCD34_CBP547_CTRL~avg.expr$LGA_3k))
#b=3.9647; a=0.6555
calcExprCtrl<-function(x,nSd=1){
  a=0.6555
  sd=0.0161
  b=3.9647
  amin<-a-sd*nSd
  amax<-a+sd*nSd
  min<-amin*x+b
  max<-amax*x+b
  return(c(min,max))
}
calcExprCtrl(1)
genesStables<-c()
for(gene in rownames(avg.expr)){
  exprCTRLrange<-calcExprCtrl(avg.expr[gene,"LGA_3k"])
  if(avg.expr[gene,"freshCD34_CBP547_CTRL"]>exprCTRLrange[1]&avg.expr[gene,"freshCD34_CBP547_CTRL"]<exprCTRLrange[2]){
    genesStables<-c(genesStables,gene)
  }
}
genesStables#18 genes only
#+/-2sd
genesStables2<-c()
for(gene in rownames(avg.expr)){
  exprCTRLrange<-calcExprCtrl(avg.expr[gene,"LGA_3k"],nSd = 2)
  if(avg.expr[gene,"freshCD34_CBP547_CTRL"]>exprCTRLrange[1]&avg.expr[gene,"freshCD34_CBP547_CTRL"]<exprCTRLrange[2]){
    genesStables2<-c(genesStables2,gene)
  }
}

length(genesStables2) #33
#+/-3sd
genesStables3<-c()
for(gene in rownames(avg.expr)){
  exprCTRLrange<-calcExprCtrl(avg.expr[gene,"LGA_3k"],nSd = 3)
  if(avg.expr[gene,"freshCD34_CBP547_CTRL"]>exprCTRLrange[1]&avg.expr[gene,"freshCD34_CBP547_CTRL"]<exprCTRLrange[2]){
    genesStables3<-c(genesStables3,gene)
  }
}

length(genesStables3) #43 genes, 
plot(avg.expr[genesStables3,"freshCD34_CBP547_CTRL"],avg.expr[genesStables3,"LGA_3k"])

#+/-4sd
genesStables4<-c()
for(gene in rownames(avg.expr)){
  exprCTRLrange<-calcExprCtrl(avg.expr[gene,"LGA_3k"],nSd = 4)
  if(avg.expr[gene,"freshCD34_CBP547_CTRL"]>exprCTRLrange[1]&avg.expr[gene,"freshCD34_CBP547_CTRL"]<exprCTRLrange[2]){
    genesStables4<-c(genesStables4,gene)
  }
}

length(genesStables4) #57 genes, stop here
plot(avg.expr[genesStables4,"freshCD34_CBP547_CTRL"],avg.expr[genesStables4,"LGA_3k"])


#now lets create our score based on this genes
#how many CpG on this genes?
library(data.table)
library(stringr)
data_all<-fread("../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)
data_all<-data.frame(data_all)
rownames(data_all)<-data_all$id
samples<-names(data_all)[str_detect(names(data_all),"CBP")]
mat<-as.matrix(data_all[,samples])

data_F<-data_all[data_all$msp1c>10^-7&
                   rowSums(is.na(mat))==0,]
#data_F<-data_F[rowSums(data_F[,samples]>10)>3,]
rm(data_all,mat)
locisF<-rownames(data_F)
annot<-read.csv2("../../ref/annotation_CpG_HELP_ALL_070420.csv",row.names = 1)

annot<-annot[locisF,]
dim(annot)  #1392170
head(annot)
locis<-rownames(annot)[annot$gene%in% genesStables4 & annot$type%in%c(4,6) & abs(annot$posAvant)<1000]
length(locis) #542
annot[locis,]
data_F[locis,samples]#some locis are 90 or 0, and other are more sparse => locis tech = locis 90 or 0
#filter for locis with possibly true zeros : with score <30
locis2<-locis[rowSums(data_F[locis,samples]<30&data_F[locis,samples]>0)==0]
length(locis2) #81
data_F[locis2,samples]
plot(density(as.matrix(data_F[locis,samples]))) #2 peaks
batch<-read.csv2("../../ref/batch_CD34_library_date_032520.csv",row.names = 1)
head(batch)

batch[samples,"newComplexityScore"]<-colSums(data_F[locis,samples]!=0)/max(colSums(data_F[locis,samples]!=0))
batch[samples,"newComplexityScore2"]<-colSums(data_F[locis2,samples]!=0)/max(colSums(data_F[locis2,samples]!=0))

correl(batch$newComplexityScore,batch$newComplexityScore2) #highly correl
plot(as.factor(batch$Group_Complexity_Fac),batch$newComplexityScore)
source("scripts/visual_function.R")
source("scripts/utils.R")
correl(batch$newComplexityScore,batch$Group) #no signif
correl(batch$newComplexityScore2,batch$Group) #no signif

correl(batch$newComplexityScore,batch$PI)
correl(batch$newComplexityScore,batch$Weight.at.term..lbs.)
correl(batch$newComplexityScore,batch$Weight..g.)
correl(batch$newComplexityScore2,batch$Weight..g.)
correl(batch$newComplexityScore,batch$Mat.Age) #r2 0.04
correl(batch$newComplexityScore2,batch$Mat.Age) 
 
write.csv2(batch,"../../ref/batch_CD34_library_date_090420.csv",row.names = T)
