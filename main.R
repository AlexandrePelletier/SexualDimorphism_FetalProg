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
source("scripts/deterDataQual.R")

#output dir
script_name <- "main"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")

#data
data_all<-fread("../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)

data_all<-data.frame(data_all)
rownames(data_all)<-data_all$id
head(data_all)
dim(data_all) #1709224     132
samples<-names(data_all)[str_detect(names(data_all),"CBP")]

batch<-read.csv2("../../ref/20-02-17_batch_119samples_IUGRLGACTRL_noNA_withHpaC.csv",header=T,row.names = 1)


#data exploration
plot(density(na.omit(as.matrix(data_all[,samples]))))
plot(density(data_all$mean))
plot(density(data_all$sd))
plot(density(data_all$pct0))


#Filtration des locis peu fiable
locis<-read.table("locisPassantQC.txt",row.names = 1,sep = ",",colClasses = "character")
head(locis)
locis<-locis$x
head(locis)

#gain en qualité
#avant
mat<-as.matrix(data_all[,samples])
deterQual(mat) #41% des locis avec Vrais zeros
deterQual(mat[locis,]) #62% de locis avec vrais zeros

deterQual2(mat,batch) #PC 1  ( 16.9 % de la variance) a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3801368616723


deterQual2(mat[locis,],batch) # PC 1  ( 18.7 % de la variance) a R2 avec Library_Complexity = 0.75 et pval = 10^ -36.5075044149743
deterQual2(mat[!(rownames(mat)%in%locis),],batch)
# [1] "PC 1  ( 7 % de la variance) a R2 avec Library_Complexity = 0.39 et pval = 10^ -13.8555937486684"
# [1] "PC 2  ( 3.3 % de la variance) a R2 avec Library_Complexity = 0.48 et pval = 10^ -17.713854390605"


deterQual2(mat[sample(rownames(mat),length(locis)),],batch)
#PC 1  ( 16.9 % de la variance) a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.4019834323515"
deterQual2(mat[sample(rownames(mat),length(locis)),],batch)
#PC 1  ( 16.9 % de la variance a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3775925094346

#petite amelioration comparé aux hasard mais on peut surement ameliorer filtration pour enlever plus de dépendance
  
#visualisation de nouvelle distribution
data_F<-data_all[locis,]
head(data_F)
sum(is.na(data_F[,samples]))

dim(data_F) #1097732     132

plot(density(as.matrix(data_F[,samples])))
plot(density(data_all$mean))
lines(density(data_F$mean),col=2)
plot(density(data_all$sd))
lines(density(data_F$sd),col=2)
plot(density(data_all$pct0))
lines(density(data_F$pct0),col=2)

hist(data_F$RankConfidenceScore,breaks = 100)
#comparé aux hasard : 
hist(data_all[sample(rownames(mat),length(locis)),"RankConfidenceScore"],breaks = 100)
hist(data_all[sample(rownames(mat),length(locis)),"RankConfidenceScore"],breaks = 100)

#heatmap 
#avant filtrage


topVarLocis1<-rownames(data_all)[order(data_all$sd,decreasing = T)[1:10000]]
head(data_all[topVarLocis,"sd"])
subMat<-mat[topVarLocis,]
p<-pheatmap(subMat,show_rownames = F,
         annotation_col = data.frame(row.names =colnames(subMat),
                                     groupe=batch[colnames(subMat),"Group_name"],
                                     sex=batch[colnames(subMat),"Gender"],
                                     batch=batch[colnames(subMat),"batch"],
                                     hpa2c=batch[colnames(subMat),"Library_Complexity"])
         
)
png(file.path(outputDir,paste("heatmapsAvantFiltration.png",sep="_")), width = 1200, height = 800)
p
dev.off()

#apres
topVarLocis2<-rownames(data_F)[order(data_F$sd,decreasing = T)[1:10000]]
sum(topVarLocis2 %in% topVarLocis1)/10000 #83% des topSD locis conservés
head(data_all[topVarLocis2,"sd"])
subMat<-mat[topVarLocis2,]
p<-pheatmap(subMat,show_rownames = F,
            annotation_col = data.frame(row.names =colnames(subMat),
                                        groupe=batch[colnames(subMat),"Group_name"],
                                        sex=batch[colnames(subMat),"Gender"],
                                        batch=batch[colnames(subMat),"batch"],
                                        hpa2c=batch[colnames(subMat),"Library_Complexity"])
            
)
png(file.path(outputDir,paste("heatmapsApresFiltration.png",sep="_")), width = 1200, height = 800)
p
dev.off()
#encore une clusterisation selon batch et Library Complexity(pas etonnant vu que 83% sont les meme locis)


#samples filtration
#
