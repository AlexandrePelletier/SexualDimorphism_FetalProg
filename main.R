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
library(VennDiagram)
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

batch<-read.csv2("../../ref/batch_CD34_library_date_032520.csv",header=T,row.names = 1)

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
write.csv2(batch,"../../ref/batch_CD34_library_date_032520.csv",row.names = T)

head(batch[,c("batch","DNA")],30)


#data exploration
plot(density(na.omit(as.matrix(data_all[,samples]))))
plot(density(data_all$mean))
plot(density(data_all$sd))
plot(density(data_all$pct0))


#FILTRATION DES LOCIS
mat<-as.matrix(data_all[,samples])
data_F<-data_all[data_all$msp1c>10^-7&
                   rowSums(is.na(mat))==0&
                   rowSums(data_all[,samples]>10,na.rm = T)>3,]
nrow(data_F) #1029401
locisF<-rownames(data_F)
#gain en qualité
#avant
deterQual(mat) #41% des locis avec Vrais zeros
deterQual(mat[locisF,]) #54% de locis avec vrais zeros

deterQual2(mat,batch) #PC 1  ( 16.9 % de la variance) a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3801368616723


deterQual2(mat[locisF,],batch) #PC 1  ( 17.8 % de la variance a R2 avec Library_Complexity = 0.74 et pval = 10^ -35.5622767364567
deterQual2(mat[!(rownames(mat)%in%locisF),],batch)
# "PC 1  ( 14.3 % de la variance a R2 avec Library_Complexity = 0.79 et pval = 10^ -40.6381624618702"

deterQual2(mat[sample(rownames(mat),length(locisF)),],batch)
#PC 1  ( 16.9 % de la variance) a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.4019834323515"
deterQual2(mat[sample(rownames(mat),length(locisF)),],batch)
#PC 1  ( 16.9 % de la variance a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3775925094346

#amelioration comparé aux hasard mais on peut surement ameliorer filtration pour enlever plus de dépendance

#visualisation de nouvelle distribution

dim(data_F) #1029401     132

plot(density(as.matrix(data_F[,samples])))
plot(density(data_all$mean))
lines(density(data_F$mean),col=2)
plot(density(data_all$sd))
lines(density(data_F$sd),col=2)
plot(density(data_all$pct0))
lines(density(data_F$pct0),col=2)

hist(data_F$RankConfidenceScore,breaks = 100)

#comparé aux hasard : 
hist(data_all[sample(rownames(mat),length(locisF)),"RankConfidenceScore"],breaks = 100)
hist(data_all[sample(rownames(mat),length(locisF)),"RankConfidenceScore"],breaks = 100)

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
sum(topVarLocis2 %in% topVarLocis1)/10000 #seulement 4.1% des topSD locis conservés !!!
head(data_all[topVarLocis2,"sd"])
subMat<-mat[topVarLocis2,]
p<-pheatmap(subMat,show_rownames = F,
            annotation_col = data.frame(row.names =colnames(subMat),
                                        groupe=batch[colnames(subMat),"Group_name"],
                                        sex=batch[colnames(subMat),"Gender"],
                                        mat.age=batch[colnames(subMat),"Mat.Age"],
                                        batch=batch[colnames(subMat),"batch"],
                                        seq=batch[colnames(subMat),"sequencing"],
                                        dnaExtr=batch[colnames(subMat),"DNA.extraction"],
                                        hpa2c=batch[colnames(subMat),"Library_Complexity"])
            
)
png(file.path(outputDir,paste("heatmapsApresFiltration.png",sep="_")), width = 1200, height = 800)
p
dev.off()
#encore une clusterisation selon batch et Library Complexity(pas etonnant vu que 80% sont les meme locis)


#FILTRATION DES SAMPLES
# en fct pct0

# batch[samples,"pct0"]<-colSums(data_all[,samples]==0,na.rm = T)/nrow(data_all)
 batch[samples,"pct0ApresF2"]<-colSums(data_F[,samples]==0,na.rm = T)/nrow(data_F)
 batch[samples,"pct0locisLo"]<-colSums(data_all[!(rownames(data_all)%in%locisF),samples]==0,na.rm = T)/sum(!(rownames(data_all)%in%locisF))
 # write.csv2(batch,"../../ref/20-02-17_batch_119samples_IUGRLGACTRL_noNA_withHpaC.csv",row.names = T)


#plot les pct0 avant filtration, colorés par le groupe
y<-batch$pct0
o<-order(y)
plot(x=1:nrow(batch),y[o],col=(batch$Group[o]+1))
plot(x=1:nrow(batch),y[o],col=(batch$batch[o]))
plot(x=1:nrow(batch),y[o],col=(batch$DNA.extraction[o]))
abline(h=0.845)
samplesToRm1<-rownames(batch)[batch$pct0>0.845]
length(samplesToRm1)

#Locis Lo
y<-batch$pct0locisLo
o<-order(y)
plot(x=1:nrow(batch),y[o],col=(batch$Group[o]+1))
plot(x=1:nrow(batch),y[o],col=(batch$batch[o]))
plot(x=1:nrow(batch),y[o],col=(batch$DNA.extraction[o]))
plot(x=1:nrow(batch),y[o],col=(batch$sequencing[o]))
abline(h=0.845)
samplesToRm1<-rownames(batch)[batch$pct0>0.845]
length(samplesToRm1)

#apres F
y<-batch$pct0ApresF2
o<-order(y)
plot(x=1:nrow(batch),y[o],col=(batch$Group[o]+1))
plot(x=1:nrow(batch),y[o],col=(batch$batch[o]))
plot(x=1:nrow(batch),y[o],col=(batch$Mat.Age_Fac[o]),main="pct0 dans les échantillons") 
plot(x=1:nrow(batch),y[o],col=(batch$sequencing[o]),main="pct0 dans les échantillons") 
plot(x=1:nrow(batch),y[o],col=(batch$DNA.extraction[o]),main="pct0 dans les échantillons")

plot(x=1:nrow(batch),y[o],col=(batch$sequencing),main="pct0 dans les échantillons")


abline(h=0.8)
samplesToRm2<-rownames(batch)[batch$pct0ApresF>0.8]
length(samplesToRm2) #21

d<-setdiff(samplesToRm2,samplesToRm1)#CBP467
col<-as.numeric(rownames(batch)%in%d)+1
plot(x=1:nrow(batch),y[o],col=col[o])
y<-batch$pct0
o<-order(y)
plot(x=1:nrow(batch),y[o],col=col[o]) 
batch[d,]


#ccl : la filtration des locis ne change presque pas la qualité des ech, il fait juste  remonter CBP467 (devient un mauvais sample)
# batch$pct0ApresF>0.8

data_F_S<-data_F[,!(names(data_F)%in%samplesToRm2)]
dim(data_F_S) # 1097732     111
samples_F<-samples[samples%in%names(data_F_S)]
length(samples_F) #98

#gagne en qualité ?
deterQual(mat[locisF,samples_F]) #54%> 56%
deterQual2(mat[locisF,samples_F],batch) 
#PC 1  ( 9.4 % de la variance a R2 avec Library_Complexity = 0.77 et pval = 10^ -32.0961766241104
#au lieu de : 
# PC 1  ( 16.9 % de la variance) a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.4019834323515"

#VISUALISATION AVEC PCA
# matNA<-mat
# matNA[is.na(matNA)]<-0
# sum(is.na(matNA))
# pc1<-prcomp(t(matNA),center = T)
# pc2<-prcomp(t(matNA[locisF,]),center = T)
# pc3<-prcomp(t(matNA[locisF,samples_F]),center = T)
# PCAlist[["pca_F"]]<-pc2
# PCAlist[["pca_F_S"]]<-pc3
# rm(pc1,pc2,pc3)
# saveRDS(PCAlist,file.path(outputDir,"PCAlist.rds"))
PCAlist<-readRDS(file.path(outputDir,"PCAlist.rds"))

#visual PCA
pcaChoose<-"pca_All"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
names(batch)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="sequencing",batch = batch,showSampleIDs=F)

plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group",batch = batch,showSampleIDs=F)

plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Sex",batch = batch,showSampleIDs=F)

#influence Covar sur PC
var_fac<-names(batch)[c(2,4,8,10,11,12,13,15,16,19,24,27,29,33,34,35)]
var_num<-names(batch)[c(5,6,17,20,21,22,23,25,26,28,32)]
varAdd<-c('GroupBatch_Complexity','GroupBatch_Complexity_Fac','pct0ApresF')
vardint<-c("Group","Group_Sex","Sex","Library_Complexity","Group_Complexity","Group_Complexity_Fac")
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd ) #date PC4 et 6

rowSums(-log10(resPV[rownames(resPV)%in%vardint,]))
# Library_Complexity              Group          Group_Sex                Sex 
# 40.975445           4.269812           5.193760           2.930691 

pcaChoose<-"pca_F"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Complexity_Fac",showSampleIDs=F)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="date",batch = batch,showSampleIDs=F)

#influence Covar sur PC
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd )
rowSums(-log10(resPV[rownames(resPV)%in%vardint,]))
# Group_Complexity   Library_Complexity                Group Group_Complexity_Fac            Group_Sex                  Sex 
# 38.669826            39.491955             7.206245            30.814766             6.855886             3.876376 

pcaChoose<-"pca_F_S"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Complexity_Fac",showSampleIDs=F)

length(LGA)
batch[LGA,]
#!conitnuer ici : regarder si LGA enlever etait les plus extrèmes et si oui, on les garde

plot(density(na.omit(as.numeric(batch[LGA,"PI"]))))
lines(density(na.omit(as.numeric(batch[CTRL,"PI"]))),col=3)
lines(density(na.omit(as.numeric(batch[LGA[LGA%in%samplesToRm2],"PI"]))),col=2) #enrich en bon LGA donc on suppr pas.

plot(density(batch[samplesToRm2,"PI"]))

#influence Covar sur PC
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd )
rowSums(-log10(resPV[rownames(resPV)%in%vardint,]))
# Group_Complexity   Library_Complexity                Group Group_Complexity_Fac            Group_Sex                  Sex 
# 31.515646            33.529202            10.371856            27.248603             8.230855             7.483574 
#on voit ici que les pvalue augmente pour les var d'int et diminue pour library comp


#on veut recupèrer également les R2
pcaChoose<-"pca_F_S"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)

resR2<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,res="r2",exclude = varAdd )
#ordonner covar importante en fct sum(R2*pctPC)

#enlever val si pval >0.01
resR2[resPV==1]<-0
pctPCs<-pctPC(PCAlist[[pcaChoose]],rngPCs =PCs1pct )

matR.pcPC<-t(resR2)*(pctPCs/100)
pctCovarExplPC<-colSums(matR.pcPC)*100
vars<-names(pctCovarExplPC)[!(names(pctCovarExplPC)%in%varAdd)]
varTrunk<-str_sub(vars,1,4)
o<-order(pctCovarExplPC[vars],decreasing = T)
barplot(pctCovarExplPC[vars][o],names.arg = varTrunk[o],cex.names = 0.9)
abline(h=1)

# varTrunk<-str_sub(names(pctCovarExplPC),1,4)
# o<-order(pctCovarExplPC,decreasing = T)
# barplot(pctCovarExplPC[o],names.arg = varTrunk[o],cex.names = 0.9)

varImportantes<-names(pctCovarExplPC[vars])[pctCovarExplPC[vars]>1]

#"Library_Complexity" "Mat.Age"            "batch"              "date"               "sequencing"

#on s'interesse a l'interaction entre sex et group, les var a modeliser sont donc :
varModelisable<-c(varImportantes,"Group","Sex","Group_Sex")

#LibraryC est une variable tech et bio ?
##est ce que library bonne var tech ? (=pas correlé a var bio : group, sex, matage, PI)
#0)
PI<-as.numeric(batch[samples_F,'PI'])
weight<-as.numeric(batch[samples_F,'Weight..g.'])
weightAtTerm<-as.numeric(batch[samples_F,"Weight.at.term..lbs."])
library<-as.numeric(batch[samples_F,"Library_Complexity"])
group_complexity<-as.numeric(batch[samples_F,"Group_Complexity"])

mat.age<-as.numeric(batch[samples_F,"Mat.Age"])


group_complexity_fac<-as.factor(batch[samples_F,"Group_Complexity_Fac"])
mat.age_fac<-as.factor(batch[samples_F,"Mat.Age_Fac"])
group<-as.factor(batch[samples_F,"Group_name"])
sex<-as.factor(batch[samples_F,"Gender"])
batches<-as.factor(batch[samples_F,"batch"])
group_sex<-as.factor(batch[samples_F,"Group_Sex"])
date<-as.factor(batch[samples_F,"date"])
sequencing<-as.factor(batch[samples_F,"sequencing"])

#1)correl ac bio ?
correl(library,PI,ret = "all")
correl(library,group,ret = "all") #no signif
correl(library,sex,ret = "all") #no signif
correl(library,weight,ret = "all")
correl(library,weightAtTerm,ret = "all")
correl(library,mat.age,ret = "all") #signif
# $p
# [1] 0.0008575249
# 
# $r2
# [1] 0.1036555
#donc library pas 100% independant/technique

#group_comp ?
correl(group_complexity,PI,ret = "all")
correl(group_complexity,group,ret = "all") #no signif
correl(group_complexity,sex,ret = "all") #no signif
correl(group_complexity,weight,ret = "all")
correl(group_complexity,weightAtTerm,ret = "all")
correl(group_complexity,mat.age,ret = "all") #encore signif mais moins

#group_comp_fac ?
correl(mat.age,group_complexity_fac,ret="all")# c'est pu signif ! donc c'est mieux


#2) PC1 explicable par group ?
pc<-PCAlist[[3]]$x

summary(lm(pc[,1]~group_complexity_fac+group+mat.age))
#0.75 !! alors que sans group > 0.67 donc on capture avec group de la var qu'on capturait pas avant

#est ce que l'interaction library:sequencing peut etre une bonne var techniqui ?
#pc1
summary(lm(pc[,1]~library)) #0.74
summary(lm(pc[,1]~library:sequencing)) #0.75
#pc2
summary(lm(pc[,2]~library)) #0.07, **
summary(lm(pc[,2]~sequencing)) #0.06
summary(lm(pc[,2]~library:sequencing)) #0.16, oui mais perd en signif
#s'ajoute au batch ?
summary(lm(pc[,2]~library:sequencing+batches))#oui, 0.25 et gagne en signif

#donc plutot que library, library:sequencing est une bonne var technique


#réponse : oui
summary(lm(pc[,2]~sequencing)) #0.06

summary(lm(pc[,1]~library:sequencing))

summary(lm(pc[,1]~library:sequencing + library:group)) #0.765


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


correl(group_complexity_fac,group,ret = "all") #indep
correl(group_complexity_fac,sex,ret = "all") #p 0.08
correl(group_complexity_fac,group_sex,ret = "all")

correl(mat.age,group_sex,ret = "all")#petite signif r2 0.1
correl(batches,group_sex,ret = "all")
correl(sequencing,group_sex,ret = "all")
#ccl : seulement colinearité entre mat.age et group mais group explique suelemnet 10% de mat.age.
batch$Mat.Age_Fac
hist(batch$Mat.Age,breaks =12 )
plot(batch$Mat.Age,batch$pct0ApresF2,col=(batch$Mat.Age_Fac)) #c'est ok
correl(mat.age_fac,group_sex) #tjr correlé donc on prend mat.age
table(mat.age_fac,group)


plot(batch$Mat.Age,batch$pct0ApresF2,col=(batch$Sex)+1)

LGA<-rownames(batch)[batch$Group_name=="L"]
IUGR<-rownames(batch)[batch$Group_name=="I"]
CTRL<-rownames(batch)[batch$Group_name=="C"]

M<-rownames(batch)[batch$Gender=="M"]
Fe<-rownames(batch)[batch$Gender=="F"]

correl(batch[LGA,"Mat.Age"],batch[LGA,"pct0ApresF2"])

correl(batch[CTRL,"Mat.Age"],batch[CTRL,"pct0ApresF2"])

correl(batch[IUGR,"Mat.Age"],batch[IUGR,"pct0ApresF2"])

correl(batch[M,"Mat.Age"],batch[M,"pct0ApresF2"])
correl(batch[Fe,"Mat.Age"],batch[Fe,"pct0ApresF2"])

LGAM<-rownames(batch)[batch$Group_name=="L"]
IUGRM<-rownames(batch)[batch$Group_name=="I"]
CTRLM<-rownames(batch)[batch$Group_name=="C"]

correl(batch$Mat.Age,batch$pct0ApresF2)

plot(batch[LGA,"Mat.Age"],batch[LGA,"pct0ApresF2"],col=batch$Mat.Age_Fac) 
plot(batch[CTRL,"Mat.Age"],batch[CTRL,"pct0ApresF2"],col=batch$Mat.Age_Fac) 
plot(batch[IUGR,"Mat.Age"],batch[IUGR,"pct0ApresF2"],col=batch$Mat.Age_Fac) 


# LIMMA model 
#A))) est ce que prendre group_complexity ameliorer inte  rpret bio des locis ?
varToModel<-varModelisable[c(6,3,4,7,10)]
samples_F_F<-samples_F[rowSums(is.na(batch[samples_F,varToModel]))==0] 
length(samples_F_F) #loose just 2 ech
sequencing<-as.factor(batch[samples_F_F,"sequencing"])
group<-as.factor(batch[samples_F_F,"Group_name"])
#sex<-as.factor(batch[samples_F_F,"Gender"])
mat.age<-as.numeric(batch[samples_F_F,"Mat.Age"])
mat.age_fac<-as.factor(batch[samples_F_F,"Mat.Age_Fac"])
batches<-as.factor(batch[samples_F_F,"batch"])
group_sex<-as.factor(batch[samples_F_F,"Group_Sex"])
group_sex<-revalue(group_sex,c("1"="FC","2"="MC","3"="FI","4"="MI","5"="FL","6"="ML"))
group_library<-as.numeric(batch[samples_F_F,"Group_Complexity"])
group_library_fac<-as.numeric(batch[samples_F_F,"Group_Complexity_Fac"])
group_mat.age_fac<-mat.age_fac:group
length(group_mat.age_fac)
formule1<-~0 + group_sex + sequencing +batches 
formule2<-~0 + group_sex + sequencing +mat.age+batches 
formule3<-~0 + group_sex + sequencing +mat.age+batches + group_library
formule4<-~0 + group_sex  +sequencing+ batches +mat.age+ group_library_fac  
formule5<-~0 + group_sex  + batches +mat.age+ group_library_fac  #mieux d'enlever sequencing car ~= group_library_fac et 21 lvl bcp ?
formule6<-~0 + group_sex  + batches + group_library_fac
models<-list(formule1,formule2,formule3,formule4,formule5)
model<-1
design<-model.matrix(models[[model]])
colnames(design)<-make.names(colnames(design))
fit <- lmFit(data_F_S[,samples_F_F], design)  
#! continue ici

cont.matrix <- makeContrasts(C.I="(group_sexFC+group_sexMC)-(group_sexFI+group_sexMI)",
                             C.L = "(group_sexFC+group_sexMC)-(group_sexFL+group_sexML)",
                             I.L = "(group_sexFI+group_sexMI)-(group_sexFL+group_sexML)",
                             MC.ML="group_sexMC-group_sexML",
                             MC.MI="group_sexMC-group_sexMI",
                             MI.ML="group_sexMI-group_sexML",
                             FC.FL="group_sexFC-group_sexFL",
                             FC.FI="group_sexFC-group_sexFI",
                             ML.FL="group_sexML-group_sexFL",
                             MI.FI="group_sexMI-group_sexFI",
                             MC.FC="group_sexMC-group_sexFC",
                             F.M="(group_sexFC+group_sexFL+group_sexFI)-(group_sexMI+group_sexMC+group_sexML)",
                             levels=design)
                             

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2) #pour formule1 et 2 : warning message :Zero sample variances detected, have been offset away from zero 
#concentrons nous sur le top1000
top1000<-topTable(fit2, adjust="BH",number = 1000)
head(top1000)
mean(top1000$adj.P.Val) #0.63

top12000<-topTable(fit2, adjust="BH",number = 12000)
max(top12000$P.Value)#0.01
plot(density(log10(fit2$p.value[fit2$p.value<0.05])))
abline(v=log10(max(top12000$P.Value))) #top12000 bien

#on enregistre top1000 locis
#locisListe<-list()
locisListe[[model]]<-rownames(top12000)
locis<-locisListe[[model]]
res<-data.frame(row.names = locis,fit2$p.value[locis,],annot[locis,c("chr","start","posAvant","gene","type","mean","Expr","Variable")],
                         data_F_S[locis,c("confidenceScore","complexity","msp1c","RankConfidenceScore")])

results <- decideTests(fit2)
sum(abs(results)) #7
colSums(abs(results))
# C.I   C.L   I.L MC.ML MC.MI MI.ML FC.FL FC.FI ML.FL MI.FI MC.FC   F.M 
# 0     0     1     0     0     2     0     0     4     0     0     0 

#enh/prom distrib
round(table(res$type)/length(res$type)*100,1) #tjr d'enrich en 4 et 6..
# 0    1    2    3    4    5    6 
# 30.6 11.0 10.2 11.5 14.8  7.1 14.3 


#distribListe<-list()
colSums(apply(res[,compas],2, function(x)return(x<0.001)))
distribListe[[model]]<-colSums(apply(res[,compas],2, function(x)return(x<0.001)))


#1) df FC :
compas<-names(res)[str_detect(names(res),"\\.")]
FCs<-data.frame()
pvals<-data.frame()
LocisParCompa<-list()
for(i in 1:length(compas)){
  print(compas[i])
  res2<- topTable(fit2,coef=compas[i],n =1000)
  if(any(res2$P.Value<0.01)){
    resF<-res2[(rownames(res2)%in%rownames(res)),]
    FCs[rownames(resF),compas[i]]<-resF$logFC
    pvals[rownames(resF),compas[i]]<-resF$P.Value
    print(paste(compas[i],"nb de locis ds top1000 :",nrow(resF)))
    LocisParCompa[[compas[i]]]<-rownames(resF)
    colors<-as.numeric(rownames(res2)%in%rownames(res))+1
    
    plot(res2$logFC,-log10(res2$P.Value),col=colors,main = compas[i])
    
    png(file.path(outputDir,paste(compas[i],"volcano_locisp0.1.png",sep="_")), width = 700, height = 500)
    plot(res2$logFC,-log10(res2$P.Value),col=colors,main = compas[i])
    dev.off()
    
  }else{
    print(paste("pas de CpG pour",compas[i]))
  }
  
  
}
head(FCs)
dim(FCs)#974  12
#FCliste<-list()
#Pvalliste<-list()
FCliste[[model]]<-FCs
Pvalliste[[model]]<-pvals

#topFC
locisSansMax<-c()
for(locus in rownames(res)){
  if(locus%in%rownames(FCs)){
    
    compasTop<-colnames(FCs)[which(sapply(FCs[locus,],abs)>20)]
    
    if(length(compasTop)==1){
      res[locus,"TopCompaSig"]<-compasTop
      res[locus,"topFC"]<-round(FCs[locus, compasTop],1)
      
    }else if(length(compasTop)>1){
      compasTopOrdered<- compasTop[order(sapply(FCs[locus, compasTop],abs),decreasing = T)]
      res[locus,"TopCompaSig"]<-compasTopOrdered[1]
      res[locus,"topFC"]<-round(FCs[locus, compasTopOrdered[1]],1)
      
      res[locus,"CompaSig"]<-paste(compasTopOrdered,collapse = "/")
      res[locus,"FC_CompaSig"]<-paste(round(FCs[locus, compasTopOrdered],0),collapse = "/")
    }else{
      #print(paste("locus : ",locus,"pas de max"))
      locisSansMax<-c(locisSansMax,locus)
    }
    
  }
  
  
  
}
length(locisSansMax)
                                    
head(res)
#resListe<-list()
resListe[[model]]<-res



write.csv2(res,file.path(outputDir,paste0("top1000Locis_Group.Sex__model",model,"BasicsLocisFilteringAndSamples_F_annotated_280320.csv")),row.names = T)

locisFC<-na.exclude(rownames(res)[res$topFC>20])
length(locisFC)#3269

hist(res$RankConfidenceScore,100)
abline(v=750000)
locisConf<-na.omit(rownames(res)[res$RankConfidenceScore>750000])
length(locisConf)#10692


#prox genes : 
plot(density(na.omit(res$posAvant))) #globale
abline(v=-50000)
abline(v=quantile(na.omit(res$posAvant),0.1),col=2) 
abline(v=quantile(na.omit(res$posAvant),0.9),col=2)
sum(res$posAvant>-20000&res$posAvant<20000,na.rm = T)/length(res$posAvant) #50% entre -20 et 20kb du TSS 

sum(res$posAvant>-2000&res$posAvant<2000,na.rm = T)/length(res$posAvant) #18% entre -2 et 2kb 

#all locis prox de genes
# randomProxAllData<-annot[sample(rownames(data_all),100000),]$posAvant
# sum(randomProxAllData>-20000&randomProxAllData<20000,na.rm = T)/length(randomProxAllData) #49% entre -20 et 20kb du TSS
# sum(randomProxAllData>-2000&randomProxAllData<2000,na.rm = T)/length(randomProxAllData) #17% entre -2 et 2kb

plot(density(na.omit(randomProxAllData[(randomProxAllData>(-50000) & randomProxAllData<(50000))])))
#locis signif :
lines(density(na.omit(res$posAvant[(res$posAvant>(-50000) & res$posAvant<(50000))])),col=2) #around 10kb of TSS

#random
hist(randomProxAllData[(randomProxAllData>(-5000) & randomProxAllData<(5000))],breaks = 100)
#signif
hist(res$posAvant[(res$posAvant>(-5000) & res$posAvant<(5000))],breaks = 100)

locisGenes10kb<-rownames(res)[res$posAvant>-10000&res$posAvant<10000]
length(locisGenes10kb)#4474
locisGenes2kb<-rownames(res)[res$posAvant>-2000&res$posAvant<2000]
length(locisGenes2kb)#2253


#CRE : 
locisCRE<-na.exclude(rownames(res)[res$type%in%c(4,6)])
length(locisCRE) #3494

# round(table(annot[sample(rownames(data_all),100000),]$type)/1000,1)
# # 0     1     2     3     4     5     6 
# # 0.304 0.108 0.101 0.118 0.146 0.072 0.146 
# round(table(annot[sample(rownames(data_F),100000),]$type)/1000,1)
# # 0    1    2    3    4    5    6 
# # 30.3 11.1 10.0 11.9 14.7  7.1 14.5 

#typeListe<-list()
typeListe[[model]]<-round(table(res$type)/length(res$type)*100,1)


#complexity : 
plot(density(data_all$complexity))
lines(density(res$complexity),col=2)
locisComplex<-rownames(res)[res$complexity>60]
length(locisComplex)#422

#vulcano #! a faire
# for(i in 1:length(compas)){
#   print(compas[i])
#   res2<- topTable(fit2,coef=compas[i],n =1000)
#   if(any(res2$P.Value<0.01)){
#     resF<-res2[(rownames(res2)%in%rownames(res)),]
#     
#     AllLocisP4F20[[compas[i]]]<-rownames(resFF)
#     
#     colors<-resFiltered[rownames(resF),"type"]+1
#     
#     
#     top20<-(rownames(resF)[order((length(resF$logFC)/rank(resF$logFC))+rank(resF$P.Value))])[1:20]
#     
#     
#     png(file.path(outputDir,paste(compas[i],"volcano.png",sep="_")), width = 700, height = 500)
#     plot(resF$logFC,-log10(resF$P.Value),col=colors,main = compas[i])
#     textxy(resF[top20,'logFC'],-log10(resF[top20,'P.Value']),resFiltered[top20,'gene'],offset= -.7)
#     dev.off()
#     
#   }else{
#     print(paste("pas de CpG pour",compas[i]))
#   }
#   
#   
# }


#locis Hi conf

locisHiConf<-intersect(intersect(intersect(intersect(locisComplex,locisConf),locisCRE),locisFC),locisGenes10kb) 
length(locisHiConf) #361
table(res[locisHiConf,"TopCompaSig"])
# C.I   C.L   F.M FC.FI FC.FL   I.L MC.FC MC.MI MC.ML MI.FI MI.ML ML.FL 
# 15    55    26     3    29    95    17     9    16     4    18    74 

#save filtr
#LocisF_liste<-list()
LocisF_liste[[model]]<-list(complex=locisComplex,conf=locisConf, CRE=locisCRE,
                            FC20=locisFC,Genes10kb=locisGenes10kb,Genes2kb=locisGenes2kb )
genes<-na.omit(res$gene)
#genesF_liste<-list()
genesF_liste[[model]]<-list(all=genes, complex=na.omit(res[locisComplex,'gene']),conf=na.omit(res[locisConf,'gene']), CRE=na.omit(res[locisCRE,'gene']),
                            FC20=na.omit(res[locisFC,'gene']),Genes10kb=na.omit(res[locisGenes10kb,'gene']),Genes2kb=na.omit(res[locisGenes2kb,'gene']) )
#locis
#locisAll.compa_list<-list()
#genesAll.compa_list<-list()
#LocisF.compa_liste<-list()
#genesF.compa_liste<-list()

# locisAll.compa_list[[model]]<-list()
# genesAll.compa_list[[model]]<-list()
LocisF.compa_liste[[model]]<-list()
genesF.compa_liste[[model]]<-list()
for(compa in compas){
  print(paste("save",compa))
  locis<-rownames(res)[res[,compa]<0.001]
  
  genes<-na.exclude(res[locis,"gene"])
  
  
  LocisF.compa_liste[[model]][[compa]]<-list(all=locis,complex=locis[locis%in%locisComplex],conf=locis[locis%in%locisConf], CRE=locis[locis%in%locisCRE],
                                             FC20=locis[locis%in%locisFC],Genes10kb=locis[locis%in%locisGenes10kb],Genes2kb=locis[locis%in%locisGenes2kb] )
  
  genesF.compa_liste[[model]][[compa]]<-list(all=genes,complex=genes[genes%in%genesF_liste[[model]]$complex],conf=genes[genes%in%genesF_liste[[model]]$conf], CRE=genes[genes%in%genesF_liste[[model]]$CRE],
                                             FC20=genes[genes%in%genesF_liste[[model]]$FC20],Genes10kb=genes[genes%in%genesF_liste[[model]]$Genes10kb],Genes2kb=genes[genes%in%genesF_liste[[model]]$Genes2kb])
  
}

saveRDS(list(locis=LocisF.compa_liste,genes=genesF.compa_liste),file = file.path(outputDir,'LocisetGenes_AllCompas_BasicsLocisFilteringAndSamples_F_280320.rds'))



#FOCUS sur LGAvs ctrl
#lga_CTRL
g<-genesF.compa_liste[[model]]$C.L$all
gM<-genesF.compa_liste[[model]]$MC.ML$all
gF<-genesF.compa_liste[[model]]$FC.FL$all


venn.diagram(
  x = list(g, gM, gF),
  category.names = c("C.L" , "CM.LM" , "CF.LF"),
  filename = file.path(outputDir,paste0('vennGenesCTRLvsLGA_model',model,'.png')),
  output=TRUE
)

intersect(g,gM)
intersect(g,gF)
intersect(gF,gM)


locisListe

#COMPA model
g1<-genesF.compa_liste[[1]]$C.L$all
gM1<-genesF.compa_liste[[1]]$MC.ML$all
gF1<-genesF.compa_liste[[1]]$FC.FL$all


g2<-genesF.compa_liste[[2]]$C.L$all
gM2<-genesF.compa_liste[[2]]$MC.ML$all
gF2<-genesF.compa_liste[[2]]$FC.FL$all

g3<-genesF.compa_liste[[3]]$C.L$all
gM3<-genesF.compa_liste[[3]]$MC.ML$all
gF3<-genesF.compa_liste[[3]]$FC.FL$all

venn.diagram(
  x = list(g1,g2, g3),
  category.names = c("noComplexity" , "Group_Complexity" , "Group_Complexity_Fac"),
  filename = file.path(outputDir,paste0('vennGenesCTRLvsLGA_compa_model1-3.png')),
  output=T
)

venn.diagram(
  x = list(gM1,gM2, gM3),
  category.names = c("noComplexity" , "Group_Complexity" , "Group_Complexity_Fac"),
  filename = file.path(outputDir,'vennGenesCTRLMvsLGAM_compa_model1-3.png'),
  output=T
)
venn.diagram(
  x = list(gF1,gF2, gF3),
  category.names = c("noComplexity" , "Group_Complexity" , "Group_Complexity_Fac"),
  filename = file.path(outputDir,'vennGenesCTRLFvsLGAF_compa_model1-3.png'),
  output=T
)

