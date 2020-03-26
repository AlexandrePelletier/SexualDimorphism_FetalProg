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
source("scripts/utils.R")
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

# #add date(v2), sequencing, dnaextraction, 
# batch2<-read.csv2("../../ref/batch_CD34_library_date_032520.csv",header = T)
# head(batch2[,25:30],30)
# batch2$V2
# batch2$DNA.extraction
# batch2$batch
# batch2$wasp.name
# batch2$V2
# batch2[batch2$batch==2,25:30] #ya 2 batch2 avec DNA bag
# #fusion des batch
# batch2<-batch2[batch2$colnames.data2..5.124.%in%samples,]
# d<-duplicated(batch2$colnames.data2..5.124.)
# batch2[d,]
# samplesDupli<-batch2$colnames.data2..5.124.[d]
# batch2[batch2$colnames.data2..5.124.%in%samplesDupli,]#les memes 
# batch2<-batch2[!d,]
# rownames(batch2)<-batch2$colnames.data2..5.124.
# batch$date<-batch2[samples,"V2"]
# batch$DNA.extraction<-batch2[samples,"DNA.extraction"]
# batch$sequencing<-batch2[samples,"Sequencing"]
# write.csv2(batch,"../../ref/batch_CD34_library_date_032520.csv",row.names = T)
head(batch)
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

#ccl : la filtration des locis ne change presque pas la qualité des ech, il fait juste  remonter CBP467 (devient un mauvais sample)
#et descendre CBP250. => on filtre les samples avec batch$pct0ApresF>0.775

data_F_S<-data_F[,!(names(data_F)%in%samplesToRm2)]
dim(data_F_S) # 1097732     109
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
# pc1<-prcomp(t(matNA),center = T)
# pc2<-prcomp(t(matNA[locis,]),center = T)
# pc3<-prcomp(t(matNA[locis,samples_F]),center = T)
# PCAlist<-list(pca_All=pc1,pca_F=pc2,pca_F_S=pc3)
# rm(pc1,pc2,pc3)
# saveRDS(PCAlist,file.path(outputDir,"PCAlist.rds"))
PCAlist<-readRDS(file.path(outputDir,"PCAlist.rds"))
pcaChoose<-"pca_All"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
names(batch)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Complexity_Fac",showSampleIDs=F)

#influence Covar sur PC
var_fac<-names(batch)[c(2,4,8,10,11,12,13,15,16,19,24,27,29)]
var_num<-names(batch)[c(5,6,17,20,21,22,23,25,26,28,32)]
varAdd<-c('Group_Complexity','GroupBatch_Complexity','GroupBatch_Complexity_Fac','Group_Complexity_Fac','pct0ApresF')
vardint<-c("Group","Group_Sex","Sex","Library_Complexity")
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd )

rowMeans(resPV[rownames(resPV)%in%vardint,])


pcaChoose<-"pca_F"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Complexity_Fac",showSampleIDs=F)
#influence Covar sur PC
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd )


pcaChoose<-"pca_F_S"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Complexity_Fac",showSampleIDs=F)
#influence Covar sur PC
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd )

#on voit ici que les pvalue augmente tous


#on veut recupèrer également les R2
pcaChoose<-"pca_F_S"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
resR2<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,res="r2",exclude = varAdd )
#ordonner covar importante en fct sum(R2*pctPC)

resR2[resR2<0]<-0
#enlever val si pval >0.01
resR2[resPV==1]<-0
pctPCs<-pctPC(PCAlist[[pcaChoose]],rngPCs =PCs1pct )

matR.pcPC<-t(resR2)*(pctPCs/100)
pctCovarExplPC<-colSums(matR.pcPC)*100
vars<-!(names(pctCovarExplPC)%in%varAdd)
barplot(sort(pctCovarExplPC[vars],decreasing = T),cex.names = 0.7)
varImportantes<-names(sort(pctCovarExplPC[vars],decreasing = T)[1:3])

#Library>mat.age>batch, (GDM > ethni > Group_Sex )

#on s'interesse a l'interaction entre sex et group, les var a modeliser sont donc :
varModelisable<-c(varImportantes,"Group","Sex","Group_Sex")

#LibraryC est une variable tech et bio,
#pour le prouver :
varModelisable
library<-as.numeric(batch[samples_F,"Library_Complexity"])
group<-as.factor(batch[samples_F,"Group_name"])
sex<-as.factor(batch[samples_F,"Gender"])
mat.age<-as.numeric(batch[samples_F,"Mat.Age"])
batches<-as.factor(batch[samples_F,"batch"])
group_sex<-as.factor(batch[samples_F,"Group_Sex"])

summary(lm(pc[,1]~library+library:group))#interaction Library:IUGR ameliore model pour expliquer PC1
summary(lm(pc[,1]~library:group)) #enlever library diminue pas model
summary(lm(pc[,2]~library:group))
summary(lm(pc[,2]~library:group+library))
summary(lm(pc[,2]~library))

summary(lm(pc[,1]~library+library:sex)) #no interaction library:sex
summary(lm(pc[,1]~library+mat.age+library:mat.age)) #l'interaction library:matage n'explique pas mieux le model
summary(lm(pc[,1]~library+library:batches)) #interaction Library:batch2 ameliore model pour expliquer PC1 que library seul
#et pour pc2 qui est expliqué par batch et library ?
summary(lm(pc[,2]~library+library:batches+batches)) #model non explicatif de pc2
summary(lm(pc[,2]~library+batches)) #library et batches2 signif mais pas intercept,

#donc l'interaction library:Group est importante, vaut mieux prendre ça que library seul

#est ce que l'interaction sex:group est importante ?
# verif que var group_sex == groupe:sex
#pc13 est ++ signif pour group_sex, l'est il aussi pour group:sex ?
summary(lm(pc[,13]~group:sex))#oui !
#les pc ne le sont pas par contre ? par exemple PC4 (ctl neg)

summary(lm(pc[,4]~group:sex))#effectivement

#pc13 est aussi expliqué par sex et group c'est bien indep ?
summary(lm(pc[,13]~sex+group)) #c'est bien l'interaction qui explique le plus le pc13
#ccl : en ajoutant des var dans le model, on perd en R2, donc 
#pour notre qu bio, on s'interesse a l'interaction groupe et sex => on met dans model group:sex

#le model serait donc : group:sex + group:library +mat.age+batch
#group library correspond à ça :
group_library<-batch[samples_F,"Group_Complexity"]

#pour que les coeff du model soit fiable, on doit eviter la colinearité entre les var du model : 
#il faut surtout que les coeff de group et sex soit bon, donc on check leur independance des autre vars


correl(mat.age,group_library,ret = "all")#correlé 
# [1] 0.002950581
# 
# $r2
# [1] 0.08134327 mais mat.age explique que 8% de libraryC

correl(group_library,batches,ret = "all") #no signif
correl(group_library,group:sex,ret = "all") #no signif


correl(mat.age,group_sex,ret = "all")#petite signif r2 0.07
#revient au meme que : 
correl(mat.age,group:sex,ret = "all")#petite signif r2 0.07n donc oui

correl(group_sex,batches,ret = "all") #no signif

#ccl : seulement colinearité entre mat.age et group mais group explique suelemnet 7% de mat.age donc on prend pas en compte

# LIMMA model
samples_F_F<-samples_F[!is.na(batch[samples_F,"Mat.Age"])]

library<-as.numeric(batch[samples_F_F,"Library_Complexity"])
group<-as.factor(batch[samples_F_F,"Group_name"])
sex<-as.factor(batch[samples_F_F,"Gender"])
mat.age<-as.numeric(batch[samples_F_F,"Mat.Age"])
batches<-as.factor(batch[samples_F_F,"batch"])
group_sex<-as.factor(batch[samples_F_F,"Group_Sex"])
group_library<-as.factor(batch[samples_F_F,"Group_Complexity"])

formule<-~0 + group:sex + group:library +mat.age+batches
design<-model.matrix(formule)
colnames(design)<-make.names(colnames(design))
 #loose just 1 ech
fit <- lmFit(data_F_S[,samples_F_F], design)  #Coefficients not estimable: groupL.sexM, why ???
#! continue ici


cont.matrix <- makeContrasts(C.I="(groupC.sexF+groupC.sexM)-(groupI.sexF+groupI.sexM)",
                             C.L = "(groupC.sexF+groupC.sexM)-(groupL.sexF+groupL.sexM)",
                             I.L = "(groupI.sexF+groupI.sexM)-(groupL.sexF+groupL.sexM)",
                             MC.ML="groupC.sexM-groupL.sexM",
                             MC.MI="groupC.sexM-groupI.sexM",
                             MI.ML="groupI.sexM-groupL.sexM",
                             FC.FL="groupC.sexF-groupL.sexF",
                             FC.FI="groupC.sexF-groupI.sexF",
                             ML.FL="groupL.sexM-groupL.sexF",
                             MI.FI="groupI.sexM-groupI.sexF",
                             MC.FC="groupC.sexM-groupC.sexF",
                             M.F="(groupC.sexF+groupL.sexF+groupI.sexF)-(groupI.sexM+groupC.sexM+groupL.sexM)",
                             
                             levels=design)
fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
top1000<-topTable(fit2, adjust="BH",number = 1000)
