#config et library
options(stringsAsFactors=F)
set.seed(12345)
library(data.table)
library(stringr)
library(limma)

source("scripts/utils.R")
#output dir 
script_name <- "estim_pval_cutoff_with_permut"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
date<-Sys.Date()
output <- file.path(outputDir,date)

#data
data_all<-fread("../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)

data_all<-data.frame(data_all)
rownames(data_all)<-data_all$id
head(data_all)

dim(data_all) #1709224     132
samples<-names(data_all)[str_detect(names(data_all),"CBP")]


batch<-read.csv2("../../ref/batch_CD34_library_date_090420.csv",header=T,row.names = 1)

#locisF
mat<-as.matrix(data_all[,samples])
data_F<-data_all[data_all$msp1c>quantile(data_all$msp1c,0.125)&
                   rowSums(is.na(mat))==0,]

#puis on retire les locis full methylated
data_F<-data_F[rowSums(data_F[,samples]>10)>4,]

#plus conf Score, nbMethylNonzeros dans pct0 elevé :
data_F<-data_F[data_F$confidenceScore>quantile(data_F$confidenceScore,0.2),]
nrow(data_F) #791613

data_F<-data_F[!(data_F$pct0>0.7&data_F$nbMethylNonZeros==0),]
nrow(data_F) #786277
locisF<-rownames(data_F)