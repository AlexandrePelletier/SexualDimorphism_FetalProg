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
source("scripts/visual_function.R")
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


#FILTRATION DES LOCIS
#see deterSeuilQC.R pour la détermination des seuils de QC
locis<-read.table("locisPassantQC.txt",row.names = 1,sep = ",",colClasses = "character")
head(locis)
locis<-locis$x
head(locis)
length(locis)

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
  
#FILTRATION
data_F<-data_all[locis,]
head(data_F)

#visualisation de nouvelle distribution
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


#FILTRATION DES SAMPLES
# en fct pct0

# batch[samples,"pct0"]<-colSums(data_all[,samples]==0,na.rm = T)/nrow(data_all)
# batch[samples,"pct0ApresF"]<-colSums(data_F[,samples]==0,na.rm = T)/nrow(data_F)
# write.csv2(batch,"../../ref/20-02-17_batch_119samples_IUGRLGACTRL_noNA_withHpaC.csv",row.names = T)


#plot les pct0 avant filtration, colorés par le groupe
y<-batch$pct0
o<-order(y)
plot(x=1:nrow(batch),y[o],col=(batch$Group[o]+1))
plot(x=1:nrow(batch),y[o],col=(batch$batch[o]))
abline(h=0.845)
samplesToRm1<-rownames(batch)[batch$pct0>0.845]
length(samplesToRm1)

#apres F
y<-batch$pct0ApresF
o<-order(y)
plot(x=1:nrow(batch),y[o],col=(batch$Group[o]+1))
plot(x=1:nrow(batch),y[o],col=(batch$batch[o]))
abline(h=0.775)
samplesToRm2<-rownames(batch)[batch$pct0ApresF>0.775]
length(samplesToRm2)

d<-setdiff(samplesToRm2,samplesToRm1)#CBP467
col<-as.numeric(rownames(batch)%in%d)+1
plot(x=1:nrow(batch),y[o],col=col[o])
y<-batch$pct0
o<-order(y)
plot(x=1:nrow(batch),y[o],col=col[o]) 

d<-setdiff(samplesToRm1,samplesToRm2)#CBP250
col<-as.numeric(rownames(batch)%in%d)+1
plot(x=1:nrow(batch),y[o],col=col[o])
y<-batch$pct0ApresF
o<-order(y)
plot(x=1:nrow(batch),y[o],col=col[o]) 

#ccl : la filtration ne change presque pas la qualité des ech, il fait juste  remonter CBP467 (devient un mauvais sample)
#et descendre CBP250. => on filtre les samples avec batch$pct0ApresF>0.775

data_F_S<-data_F[,!(names(data_F)%in%samplesToRm2)]
dim(data_F_S)
samples_F<-samples[samples%in%names(data_F_S)]
length(samples_F) #96

#gagne en qualité ?
deterQual(mat[locis,samples_F]) #62% > 65%
deterQual2(mat[locis,samples_F],batch) 
#PC 1  ( 10.2 % de la variance a R2 avec Library_Complexity = 0.74 et pval = 10^ -29.0080916905632
#au lieu de : 
# PC 1  ( 18.7 % de la variance) a R2 avec Library_Complexity = 0.75 et pval = 10^ -36.5075044149743

#VISUALISATION AVEC PCA
matNA<-mat
matNA[is.na(matNA)]<-0
sum(is.na(matNA))
pc1<-prcomp(t(matNA),center = T)
pc2<-prcomp(t(matNA[locis,]),center = T)
pc3<-prcomp(t(matNA[locis,samples_F]),center = T)

PCAlist<-list(pca_All=pc1,pca_F=pc2,pca_F_S=pc3)
rm(pc1,pc2,pc3)
saveRDS(PCAlist,file.path(outputDir,"PCAlist.rds"))

plotPCVarExplain(PCAlist[[3]],1:40,lineSeuilPct = 1)

#PCA
#PCAlist<-readRDS(file.path(outputDir,"PCAlist.rds"))
names(batch)
plotPCA(PCAlist[[3]],PCx=1,PCy=2,colorBatch="Group_Complexity_Fac",showSampleIDs=F)


#influence Covar sur PC
pc<-data.frame(pc3$x)
pcs<-1:30
names(batch)
var_fac<-names(batch)[c(2,4,8,10,11,12,13,15,16,19,24,27,29)]
var_num<-names(batch)[c(5,6,17,20,21,22,23,25,26,28,32)]

batch_num=batch[rownames(pc),var_num]
batch_fac=batch[rownames(pc),var_fac]
batch_num=t(batch_num)
batch_fac=t(batch_fac)
split_num=split(batch_num,rownames(batch_num))
split_fac=split(batch_fac,rownames(batch_fac))

pv_num=lapply(split_num,function(x){
  FAC1.p<-rep(0,length(pcs))
  for (i in pcs){
    FAC1<-as.numeric(x)
    FAC1<-lm(pc[,i]~FAC1)
    FAC1.p[i]<-anova(FAC1)$Pr[1]
  }
  return(FAC1.p)})

R_num=lapply(split_num,function(x){
  FAC1.r2<-rep(0,length(pcs))
  for (i in pcs){
    FAC1<-as.numeric(x)
    FAC1<-lm(pc[,i]~FAC1)
    FAC1.r2[i]<-summary(FAC1)$adj.r.squared
    
  }
  return(FAC1.r2)})

pv_fac=lapply(split_fac,function(x){
  FAC1.p<-rep(0,length(pcs))
  for (i in pcs){
    FAC1<-as.factor(x)
    FAC1<-lm(pc[,i]~FAC1)
    FAC1.p[i]<- anova(FAC1)$Pr[1]
  }
  return(FAC1.p)})

R_fac=lapply(split_fac,function(x){
  FAC1.r2<-rep(0,length(pcs))
  for (i in pcs){
    FAC1<-as.factor(x)
    FAC1<-lm(pc[,i]~FAC1)
    FAC1.r2[i]<-summary(FAC1)$adj.r.squared
    
  }
  return(FAC1.r2)})
pvals.num<-do.call(rbind,pv_num)
pvals.fac<-do.call(rbind,pv_fac)
final_pv=rbind(pvals.num,pvals.fac)
pv2=data.matrix(final_pv)
pv2[which(pv2>0.05)]<-1 ####here I basicaly put them to 1 if less than 0.05
logpvals.raw<--log10(pv2)

pheatmap(logpvals.raw[,1:10],cluster_rows = F,cluster_cols = F,labels_col= paste("PC",1:10),display_numbers = T)
pheatmap(logpvals.raw[,11:30],cluster_rows = F,cluster_cols = F,labels_col= paste("PC",11:30),display_numbers = T)

#pheatmaps des R2 aussi 

