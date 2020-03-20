#determiner les CpG différentiellement méthylé entre les sexes a 'linterieur d'un groupe

#config et library
options(stringsAsFactors=F)
library(data.table)
library(stringr)
library(plyr)
library(calibrate)
library(pheatmap)
library(limma)
library('clusterProfiler')
library(org.Hs.eg.db)
set.seed(12345)
#source("scripts/deterQual.R")

#output dir
script_name <- "main"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")

#data
data_all<-fread("../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)

data_all<-data.frame(data_all)
rownames(data_all)<-data_all$id
head(data_all)
dim(data_all) #1709224     127
samples<-names(data_all)[str_detect(names(data_all),"CBP")]

batch<-read.csv2("../../ref/20-02-17_batch_119samples_IUGRLGACTRL_noNA_withHpaC.csv",header=T,row.names = 1)


#data exploration
plot(density(na.omit(as.matrix(data_all[,samples]))))
plot(density(data_all$mean))
plot(density(data_all$sd))
plot(density(data_all$pct0))
#etabli selon le script deter_QC_filter,  msp1c<q10%, nbMethylNon0 ==0 sauf sipct0<0.7
locis<-na.omit(rownames(data_all)[data_all$msp1c>quantile(data_all$msp1c,0.08)&(
  (data_all$nbMethylNonZeros>0)|(data_all$pct0<0.7))])
