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
source("scripts/FindInfosGenes.R")
source("scripts/scRNA-seq integration.R")
#output dir 
script_name <- "main"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
date<-Sys.Date()
output <- file.path(outputDir,date)

#PARAMS
filtres<-"locisF.msp1.NA.fullMethyl"

#data
data_all<-fread("../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)

data_all<-data.frame(data_all)
rownames(data_all)<-data_all$id
head(data_all)
dim(data_all) #1709224     132
samples<-names(data_all)[str_detect(names(data_all),"CBP")]

batch<-read.csv2("../../ref/batch_CD34_library_date_032520.csv",header=T,row.names = 1)


head(batch)
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
#write.csv2(batch,"../../ref/batch_CD34_library_date_032520.csv",row.names = T)

head(batch)


#data exploration
plot(density(na.omit(as.matrix(data_all[,samples]))))
plot(density(data_all$mean))
plot(density(data_all$sd))
plot(density(data_all$pct0))


#FILTRATION DES LOCIS
mat<-as.matrix(data_all[,samples])
data_F<-data_all[data_all$msp1c>10^-7&
                   rowSums(is.na(mat))==0&
                   ,]
data_F<-data_F[rowSums(data_F[,samples]>10)>3,]
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


#VISUALISATION DES SAMPLES
# en pct0

# batch[samples,"pct0"]<-colSums(data_all[,samples]==0,na.rm = T)/nrow(data_all)
 batch[samples,"pct0ApresF"]<-colSums(data_F[,samples]==0,na.rm = T)/nrow(data_F)
 #batch[samples,"pct0locisLo"]<-colSums(data_all[!(rownames(data_all)%in%locisF),samples]==0,na.rm = T)/sum(!(rownames(data_all)%in%locisF))
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


#apres F
y<-batch$pct0ApresF
o<-order(y)
plot(x=1:nrow(batch),y[o],col=(batch$Group[o]+1))
plot(x=1:nrow(batch),y[o],col=(batch$batch[o]))
plot(x=1:nrow(batch),y[o],col=(batch$Mat.Age_Fac[o]),main="pct0 dans les échantillons") 
plot(x=1:nrow(batch),y[o],col=(batch$sequencing[o]),main="pct0 dans les échantillons") 
plot(x=1:nrow(batch),y[o],col=(batch$DNA.extraction[o]),main="pct0 dans les échantillons")

plot(x=1:nrow(batch),y[o],col=(batch$sequencing),main="pct0 dans les échantillons")


abline(h=0.8)
#samplesToRm<-rownames(batch)[batch$pct0ApresF>0.8]


d<-setdiff(samplesToRm,samplesToRm1)#CBP467
col<-as.numeric(rownames(batch)%in%d)+1
plot(x=1:nrow(batch),y[o],col=col[o])
y<-batch$pct0
o<-order(y)
plot(x=1:nrow(batch),y[o],col=col[o]) 
batch[d,]


#ccl : la filtration des locis ne change presque pas la qualité des ech, il fait juste  remonter CBP467 (devient un mauvais sample)
# batch$pct0ApresF>0.8
#data_F_S<-data_F[,!(names(data_F)%in%samplesToRm)]

#gagne en qualité ?
deterQual(mat[locisF,samples]) #54%> 56%
deterQual2(mat[locisF,samples],batch) 
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

#visual PCA et INFLUENCE COVAR SUR PC
pcaChoose<-"pca_All"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
names(batch)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="sequencing",batch = batch,showSampleIDs=F)

plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group",batch = batch,showSampleIDs=F)

plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Sex",batch = batch,showSampleIDs=F)

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
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd )
rowSums(-log10(resPV[rownames(resPV)%in%vardint,]))
# Group_Complexity   Library_Complexity                Group Group_Complexity_Fac            Group_Sex                  Sex 
# 38.669826            39.491955             7.206245            30.814766             6.855886             3.876376 



pcaChoose<-"pca_F"
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


varImportantes<-names(pctCovarExplPC[vars])[pctCovarExplPC[vars]>1]
# [1] "Group_Complexity"     "Library_Complexity"  
# [3] "Mat.Age"              "batch"               
# [5] "date"                 "Group"               
# [7] "Group_Complexity_Fac" "Group_Sex"           
# [9] "latino"               "sequencing"   



###pk on prend pas samples filtres
# length(LGA)
# batch[LGA,]
# #!conitnuer ici : regarder si LGA enlever etait les plus extrèmes et si oui, on les garde
# 
# plot(density(na.omit(as.numeric(batch[LGA,"PI"]))))
# lines(density(na.omit(as.numeric(batch[CTRL,"PI"]))),col=3)
# lines(density(na.omit(as.numeric(batch[LGA[LGA%in%samplesToRm2],"PI"]))),col=2) #enrich en bon LGA donc on suppr pas.
# 


#on s'interesse a l'interaction entre sex et group, les var a modeliser sont donc :
varModelisable<-c(varImportantes)

#LibraryComplexity est une variable tech et bio ?
##est ce que library bonne var tech ? (=pas correlé a var bio : group, sex, matage, PI)
#0)
PI<-as.numeric(batch[samples,'PI'])
weight<-as.numeric(batch[samples,'Weight..g.'])
weightAtTerm<-as.numeric(batch[samples,"Weight.at.term..lbs."])
complexity<-as.numeric(batch[samples,"Library_Complexity"])
group_complexity<-as.numeric(batch[samples,"Group_Complexity"])

mat.age<-as.numeric(batch[samples,"Mat.Age"])


group_complexity_fac<-as.factor(batch[samples,"Group_Complexity_Fac"])
mat.age_fac<-as.factor(batch[samples,"Mat.Age_Fac"])
group<-as.factor(batch[samples,"Group_name"])
sex<-as.factor(batch[samples,"Gender"])
batches<-as.factor(batch[samples,"batch"])
group_sex<-as.factor(batch[samples,"Group_Sex"])
date<-as.factor(batch[samples,"date"])
sequencing<-as.factor(batch[samples,"sequencing"])
latino<-as.factor(batch[samples,"latino"])
#correl ac bio ?
correl(complexity,PI,ret = "all") #p= 0.01, r=0.05
correl(complexity,group,ret = "all") #pareil
correl(complexity,sex,ret = "all") #no signif
correl(complexity,weight,ret = "all")#p= 0.04, r=0.02
correl(complexity,weightAtTerm,ret = "all") #no
correl(complexity,mat.age,ret = "all") #signif++ p=0.001, r2= 0.08
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
correl(group_complexity,mat.age,ret = "all") #encore signif mais moins r2= 0.06

#group_comp_fac ?
correl(mat.age,group_complexity_fac,ret="all")# c'est pu signif ! donc c'est mieux


#correl entre groupe_sex et (group_complexity_fac,sequencing, mat.age, latino,batch)
correl(group_sex,sequencing,ret="all") #no

correl(mat.age,group_sex,ret="all") #yes, r2 0.09

correl(latino,group_sex,ret="all") #no
correl(batches,group_sex,ret="all") #no

correl(batches,group_sex,ret="all") #no

#mat.age correler avec group, le prendre dans model ou pas ?
pc<-PCAlist[[2]]$x
summary(lm(pc[,1]~group_complexity_fac)) #r2=0.63
summary(lm(pc[,1]~group_complexity_fac+group)) #r2 = 0.70
summary(lm(pc[,1]~group_complexity_fac+group+mat.age)) #r2 0.68, donc vaut mieux pas mettre mat.age
summary(lm(pc[,1]~group+mat.age))
summary(lm(pc[,1]~group_complexity_fac+group+mat.age+latino))  #0.79 !
summary(lm(pc[,1]~group_complexity_fac+group+latino))  #0.796 c'est mieux.
summary(lm(pc[,1]~group_complexity_fac+group+latino+sequencing)) #0.82 ! 
summary(lm(pc[,1]~group_complexity_fac+group+latino+sequencing+mat.age)) #0.83 !

#ccl mat.age peu nécessaire pour expliquer pc1, mais peut etre bon pour expliquer PC6
# on capture avec group de la var que complexity library n'expliquait seul

#le model opti serait donc :  model ~0 + group_sex + sequencing + batches + group_complexity_fac + latino + mat.age


### histoire de la colinearité entre group et mat.age
plot(batch$Mat.Age,batch$pct0ApresF,col=(batch$Group_Sex)+1) 

LGA<-rownames(batch)[batch$Group_name=="L"]
IUGR<-rownames(batch)[batch$Group_name=="I"]
CTRL<-rownames(batch)[batch$Group_name=="C"]

M<-rownames(batch)[batch$Gender=="M"]
Fe<-rownames(batch)[batch$Gender=="F"]

correl(batch[LGA,"Mat.Age"],batch[LGA,"pct0ApresF"])

correl(batch[CTRL,"Mat.Age"],batch[CTRL,"pct0ApresF"])

correl(batch[IUGR,"Mat.Age"],batch[IUGR,"pct0ApresF"])

correl(batch[M,"Mat.Age"],batch[M,"pct0ApresF"])
correl(batch[Fe,"Mat.Age"],batch[Fe,"pct0ApresF"])
LGAM<-rownames(batch)[batch$Group_name=="L"]
IUGRM<-rownames(batch)[batch$Group_name=="I"]
CTRLM<-rownames(batch)[batch$Group_name=="C"]
correl(batch$Mat.Age,batch$pct0ApresF2)
plot(batch[LGA,"Mat.Age"],batch[LGA,"pct0ApresF2"],col=batch$Mat.Age_Fac) 
plot(batch[CTRL,"Mat.Age"],batch[CTRL,"pct0ApresF2"],col=batch$Mat.Age_Fac) 
plot(batch[IUGR,"Mat.Age"],batch[IUGR,"pct0ApresF2"],col=batch$Mat.Age_Fac) 


# LIMMA CHOIX MODEL 
# model de ! group_sex, sequencing, batch, group_comp_fac , latino, 
#coix bas" sur l'enrichissement enh /prom et en distance du gene le plus proche
annot<-read.table("../annotation_CpG_HELP_ALL_011520.txt",sep="\t",header=T)
colnames(annot)<-c("chr","start","stop","id","V5","posAvant","ENSEMBL_ID","gene","type")
rownames(annot)<-annot$id
head(annot)


varToModel<-c("Group_Sex","sequencing",'batch',"Group_Complexity_Fac","Mat.Age","latino","GroupBatch_Complexity_Fac")
samples_F_F<-samples[rowSums(is.na(batch[samples,varToModel]))==0] 
length(samples_F_F) #107
sequencing<-factor(batch[samples_F_F,"sequencing"])
group<-factor(batch[samples_F_F,"Group_name"])
groupBatch_complexity_fac<-factor(batch[samples_F_F,"GroupBatch_Complexity_Fac"])
#sex<-as.factor(batch[samples_F_F,"Gender"])
mat.age<-as.numeric(batch[samples_F_F,"Mat.Age"])

batches<-factor(batch[samples_F_F,"batch"])
group_sex<-factor(batch[samples_F_F,"Group_Sex"])
group_sex<-revalue(group_sex,c("1"="FC","2"="MC","3"="FI","4"="MI","5"="FL","6"="ML"))
group_complexity<-as.numeric(batch[samples_F_F,"Group_Complexity"])
group_complexity_fac<-factor(batch[samples_F_F,"Group_Complexity_Fac"])
latino<-factor(batch[samples_F_F,"latino"])
#complexity<-as.numeric(batch[samples_F_F,"Library_Complexity"])

models<-list()
model<-6
formule<- ~0 + group_sex  + batches + group_complexity + latino + mat.age
models[[model]]<-formule
design<-model.matrix(models[[model]])
#resModels<-list()

colnames(design)<-make.names(colnames(design))
fit <- lmFit(data_F[,samples_F_F], design)  

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
fit2  <- eBayes(fit2) #warning message :Zero sample variances detected, have been offset away from zero 

results <- decideTests(fit2)

sum(abs(results)) #>7>1268 > 2 >15
colSums(abs(results))


# C.I   C.L   I.L MC.ML MC.MI MI.ML FC.FL FC.FI ML.FL MI.FI MC.FC   F.M 
# 0     3     1     1     0     0     0     0     1     1     0     0 
# C.I   C.L   I.L MC.ML MC.MI MI.ML FC.FL FC.FI ML.FL MI.FI MC.FC   F.M 
# 0    45  1196     0     0     9     4     0     7     7     0     0 

# C.I   C.L   I.L MC.ML MC.MI MI.ML FC.FL FC.FI ML.FL MI.FI MC.FC   F.M 
# 0     0     0     0     0     0     0     0     1     1     0     0 

# C.I   C.L   I.L MC.ML MC.MI MI.ML FC.FL FC.FI ML.FL MI.FI MC.FC   F.M 
# 0     0     4     0     2     0     0     0     3     6     0     0 


#bon model ? 1) enrichissment en enh et prom, 2) prox du gene
locisSig<-rownames(fit2$p.value)[apply(fit2$p.value<0.001,1,any)]
length(locisSig) #8617 > 18k > 26k > 12k > 15k
resSig<-data.frame(row.names = locisSig,fit2$p.value[locisSig,],annot[locisSig,c("chr","start","posAvant","gene","type")],
                   data_F[locisSig,c("confidenceScore","confidenceScoreNorm","complexity","msp1c","RankConfidenceScore")])


head(resSig,100)
#nb C-L 
sum(resSig$C.L<0.001)  #>3692>>869> 1079

#enrichissement en bon locis :
mean(resSig$RankConfidenceScore)/mean(na.omit(data_all$RankConfidenceScore)) #>1.35>1.39>1.28 > 1.31

#enh et prom

locisRall<-sample(rownames(data_all),100000)
locisRF<-sample(rownames(data_F),100000)

table(annot[locisRall,"type"])/length(locisRall)*100
# 0      1      2      3      4      5      6 
# 30.383 10.877 10.145 11.730 14.485  7.131 14.735 
table(annot[locisRF,"type"])/length(locisRF)*100
# 0      1      2      3      4      5      6 
# 25.488 10.194  9.289 10.180 17.424  7.017 19.974 

table(resSig$type)/length(resSig$type)*100
# 0         1         2         3         4         5         6 
# 25.089938 10.154346  9.144714  9.423233 18.324243  7.229894 20.262272 

# 0         1         2         3         4         5         6 
# 17.098901  7.780220  6.659341  6.967033 22.620879  6.098901 32.478022 

# 0         1         2         3         4         5         6 
# 13.877926  6.417870  5.076823  5.670715 23.893636  4.915897 39.840607 

# 0         1         2         3         4         5         6 
# 21.650452  9.097116  7.817034  8.474142 19.994880  6.716163 25.849121 

# 0         1         2         3         4         5         6 
# 19.981805  8.148678  6.888037  7.602833 20.865553  5.965300 30.092924 


#locis around5k of the gene
locis5k<-na.omit(rownames(annot)[annot$posAvant<5000&annot$posAvant>-5000])

sum(locisRall%in%locis5k)/length(locisRall) # 0.26

sum(locisRF%in%locis5k)/length(locisRF) #0.32

sum(locisSig%in%locis5k)/length(locisSig)  #0.34 > 0.46 > 0.54 > 0.39 > 0.43

hist(annot[locisRall[locisRall%in%locis5k],"posAvant"],breaks = 100)
hist(annot[locisRF[locisRF%in%locis5k],"posAvant"],breaks = 100)
hist(annot[locisSig[locisSig%in%locis5k],"posAvant"],breaks = 100)


#MODEL De confiance : 3 et 4
model<-2
formule<-~0 + group_sex  +batches +latino +mat.age  + group_complexity_fac +sequencing
models[[model]]<-formule
design<-model.matrix(models[[model]])
colnames(design)<-make.names(colnames(design))
fit <- lmFit(data_F[,samples_F_F], design)  

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
fit2  <- eBayes(fit2) #warning message :Zero sample variances detected, have been offset away from zero 

#concentrons nous sur le top1000 par compa
compas<-colnames(fit2$p.value)

#deter seuil pval pour avoir assez de locis signif
for(compa in compas){
  print(compa)
  top1000<-topTable(fit2,coef = compa, adjust="BH",number = 1000)
  #remove F si pval <0.001
  print(sum(top1000$P.Value<0.001))
 
  
}


FCs<-data.frame()
resParCompa<-list()
seuilPval<-0.001
for(compa in compas){
  print(compa)
  res<- topTable(fit2,coef=compa,n =Inf)
  
  if(any(res$P.Value<seuilPval)){
    resF<-res[res$P.Value<seuilPval,]
    print(paste(nrow(resF),"locis passe filtre pval"))
    FCs[rownames(resF),compa]<-resF$logFC
    locis<-rownames(resF)
    resA<-data.frame(row.names = locis,
                                         FC=resF$logFC,
                                         pval=resF$P.Value,
                                         pval.adj=resF$adj.P.Val,
                                         topCompa=compas[as.vector(apply(fit2$p.value[locis,compas],1,which.min))],annot[locis,c("chr","start","posAvant","gene","type")],
                                         data_F[locis,c("confidenceScore","complexity","msp1c","RankConfidenceScore")]
                                         
                                         )
    
    resParCompa[[compa]]<-resA
    
    #enh/prom distrib
    barplot(table(resA$type)/length(resA$type)*100,main = compa)
  }else{
    print(paste("pas de CpG pour",compa))
    resParCompa[[compa]]<-NA
  }
  
  
}

head(resParCompa)


#RES PAR COMPAS
#scRNA-seq annot
#genes Expr
CTRL_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInCTRL.csv",row.names = 1)
LGA_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInLGA.csv",row.names = 1)

listExprMat<-list(CTRL=CTRL_expr_by_pop,LGA=LGA_expr_by_pop)
markers<-read.csv2("../../../Alexandre_SC/analyses/test_batch_integration/all.markers_SCTransform_CTRLandLGA_integrated.csv",row.names = 1)
head(markers)
#reannot 
library(Seurat)
source("scripts/scoreCluster.R")
samples<-readRDS("../../../Alexandre_SC/analyses/test_batch_integration/CTRLandLGA_SCTransform_Integrated.rds")
new.cluster.ids <- c("HSC-SELL2", "HSC-AVP", "HSC-Er", "HSC-SELL1", "EMP", "GMP",
                     "LyP", "MkP", "proT","LMPP","preB","LT-HSC","Neu","LyB","Ly-ETS1","DC")
names(new.cluster.ids) <- levels(samples)
samples <- RenameIdents(samples, new.cluster.ids)

for(i in 1:length(new.cluster.ids)){
  markers$cluster[markers$cluster==as.numeric(names(new.cluster.ids)[i])]<-new.cluster.ids[i]
}
head(markers)
scoresCluster<-scoreMarquageCluster(markers,samples,seuil = "intraClusterFixe",filtreMin = 2)

resGenesParCompa<-list()
for(compa in compas){
  print(compa)
  #RES PAR GENES
  #make a df with rownames = genes, and  count of locis by genes 
  resGenes<-resLocisToGenes(resParCompa[[compa]])
  
  #add genes Expr in bulk scrnaseq and in subpop
  resGenes<-findGenesExpr(resGenes,listExprMat)
  #add markers of clusters
  resGenes<-addMarqueursClusters(scoresCluster,resGenes)
  #add canonical/published markers of cell types
  resGenes<-findMarqueursPops(rownames(resGenes),df = resGenes)
  
  #add fct and tags
  try(resGenes<-find_fonction(rownames(resGenes),df = resGenes,save=T))
  try(resGenes<-findIfTF(rownames(resGenes),resGenes,save=F,large = T))
  try(resGenes<-findCbIn(rownames(resGenes),
                         listeKeywords = list(HSPC="hematopo|myeloid|lymphoid|HSC|HSPC",
                                              lineage="lineage decision|differentiation|cell fate",
                                              stress="stress",
                                              signaling="kinase|signaling|pathway"),
                         genes_infos=resGenes))
  #save 
  resGenesParCompa[[compa]]<-resGenes
  
  
}


#write resLocis
FCm<-mergeCols(FCs,mergeColsName = "FC",top = 4,filter = 20,abs = T,roundNum = 0)
head(FCm)
colsToMerge<-colnames(resGenes)[2:(which(colnames(resGenes)=="bulk_CTRL")-1)]
for(compa in compas){
  
  print(compa)
  resGenes<-resGenesParCompa[[compa]]
  #garder seulement expr des 5 plus grande pop pour save ds csv
  resGenes.merge<-mergeCols(df = resGenes,colsToMerge = colsToMerge,mergeColsName = "Expr",filter = 0.1,top = 5,abs = F)
  
  #RES PAR LOCIS
  resLocis<-resParCompa[[compa]]
  resLocis<-data.frame(row.names = rownames(resLocis),resLocis,FCm[rownames(resLocis),])
  
  resLocis<-annotLocis(resLocis = resLocis,resGenes =resGenes.merge )
  
  #save
  resParCompa[[compa]]<-resLocis
  write.csv2(resLocis,paste(output,"res_locis_in",compa,"top",nrow(resLocis),"pval",seuilPval,filtres,"model",model,".csv",sep = "_"),row.names = T,na = "")
  
}
#SAVE ALL RES DU MODEL
#resModels<-list()
resModels[[model]]<-list(model=models[[model]],
                         locisSig=locisSig,
                         res=resSig,
                         seuilSig=seuilPval,
                         distrib.compas=colSums(apply(resSig[,compas],2, function(x)return(x<seuilPval))),
                         distrib.features=round(table(resSig$type)/length(resSig$type)*100,1),
                         res.compas=resParCompa,
                         resGenes.compas=resGenesParCompa
                         
                         )
                         
                         
resModels[[model]]$distrib.compas
saveRDS(resModels[[model]],file = paste(output,"LocisetGenes_AllCompas_pval",seuilPval,filtres,"model",model,'.rds'))


#LES FILTRATIONS DES LOCIS SIG  : !a revoir 
# #res generales
# locisSig<-c()
# for(compa in compas){
#   locisSig<-union(locisSig,rownames(resParCompa[[compa]]))
# }
# 
# length(locisSig) #25k
# resSig<-data.frame(row.names = locisSig,)
# #FC>30
# 
# locisFC<-list()
# for (compa in compas){
#   locisFC[[compa]]<-na.omit(rownames(resParCompa[[compa]])[resParCompa[[compa]]$FC>20]) 
#   
#   
# }
# 
# 
# locisConf<-na.omit(rownames(res)[res$RankConfidenceScore>750000])
# #dist from TSS
# locisGenes30kb<-rownames(res)[res$posAvant>-30000&res$posAvant<30000]
# locisGenes2kb<-rownames(res)[res$posAvant>-2000&res$posAvant<2000]
# #CRE (enhancer 4, prom 6)
# locisCRE<-na.exclude(rownames(res)[res$type%in%c(4,6)])
# #complexity >60 : 
# locisComplex<-rownames(res)[res$complexity>60]
# 
# #locis Expr in sc
# locisExpr<-annot[]
# 
# 
# #locis change in sc
# locisExpr
# 
# #locis Hi conf
# locisHiConf<-intersect(intersect(intersect(intersect(locisComplex,locisConf),locisCRE),locisFC),locisGenes30kb) 
# 

#SAVE ALL RES DU MODEL
resModels[[model]]<-list(model=models[[model]],
                         locisSig=rownames(top12000),
                         res=res,
                         seuilSig=max(top12000$P.Value),
                         distrib.compas=colSums(apply(res[,compas],2, function(x)return(x<max(top12000$P.Value)))),
                         distrib.features=round(table(res$type)/length(res$type)*100,1),
                         res.compas=resParCompa,
                         resGenes.compas=resGenesParCompa,
                         locisF=list(complex=locisComplex,conf=locisConf, CRE=locisCRE,
                                     FC20=locisFC,Genes30kb=locisGenes30kb,Genes2kb=locisGenes2kb,HiConf= locisHiConf),
                         distrib.hi.conf=table(res[locisHiConf,"topCompa"])
                         
                         
)
resModels[[model]]$model
resModels[[model]]$seuilSig
head(resModels[[model]]$locisSig)
head(resModels[[model]]$distrib.compas)
saveRDS(resModels[[model]],file = paste(output,"LocisetGenes_AllCompas_pval",seuilPval,filtres,"model",model,'.rds'))


##F.M
#vulcano 
CM.CF<-read.csv2("../Alexandre_Methyl/analyses/main/2020-04-01_res_locis_in_MC.FC_top_1224_pval_0.001_locisF.msp1.NA.fullMethyl_model_4_.csv")
C.L<-read.csv2("../Alexandre_Methyl/analyses/main/2020-04-01_res_locis_in_C.L_top_4878_pval_0.001_locisF.msp1.NA.fullMethyl_model_4_.csv")
CM.LM<-read.csv2("../Alexandre_Methyl/analyses/main/2020-04-01_res_locis_in_MC.ML_top_1758_pval_0.001_locisF.msp1.NA.fullMethyl_model_4_.csv")
CF.LF<-read.csv2("../Alexandre_Methyl/analyses/main/2020-04-01_res_locis_in_FC.FL_top_4495_pval_0.001_locisF.msp1.NA.fullMethyl_model_4_.csv")

##F.M in ctrl
library(calibrate)
names(CM.CF)
seuilFC<-30
seuilpval<-10^-4
colors<-as.numeric(CM.CF$pval<seuilpval & abs(CM.CF$FC)>seuilFC)+1
#plot(1:7,col=1:7)
#top30<-(rownames(CM.CF)[order((length(CM.CF$FC)/rank(CM.CF$FC))+rank(CM.CF$pval))])[1:30]
png(paste(output,"CM.CF","vulcano.png",sep="_"), width = 700, height = 500)
    
    plot(CM.CF$FC,-log10(CM.CF$pval),col=colors,main = "CM.CF",pch=18)
    
    abline(h=-log10(seuilpval))
    abline(v=seuilFC)
    abline(v=-seuilFC)
    #textxy(CM.CF[top30,'FC'],-log10(CM.CF[top30,'pval']),CM.CF[top30,'gene'],offset= -.7)
    
dev.off()
    
#C.L 
colors<-as.numeric(C.L$pval<seuilpval & abs(C.L$FC)>seuilFC)+1
png(paste(output,"C.L","vulcano.png",sep="_"), width = 700, height = 500)

plot(C.L$FC,-log10(C.L$pval),col=colors,main = "C.L",pch=18)

abline(h=-log10(seuilpval))
abline(v=seuilFC)
abline(v=-seuilFC)
#textxy(C.L[top30,'FC'],-log10(C.L[top30,'pval']),C.L[top30,'gene'],offset= -.7)

dev.off()

#CM.LM 
colors<-as.numeric(CM.LM$pval<seuilpval & abs(CM.LM$FC)>seuilFC)+1
png(paste(output,"CM.LM","vulcano.png",sep="_"), width = 700, height = 500)

plot(CM.LM$FC,-log10(CM.LM$pval),col=colors,main = "CM.LM",pch=18)

abline(h=-log10(seuilpval))
abline(v=seuilFC)
abline(v=-seuilFC)
#textxy(CM.LM[top30,'FC'],-log10(CM.LM[top30,'pval']),CM.LM[top30,'gene'],offset= -.7)

dev.off()

#CF.LF 
colors<-as.numeric(CF.LF$pval<seuilpval & abs(CF.LF$FC)>seuilFC)+1
png(paste(output,"CF.LF","vulcano.png",sep="_"), width = 700, height = 500)

plot(CF.LF$FC,-log10(CF.LF$pval),col=colors,main = "CF.LF",pch=18)

abline(h=-log10(seuilpval))
abline(v=seuilFC)
abline(v=-seuilFC)
#textxy(CM.LM[top30,'FC'],-log10(CM.LM[top30,'pval']),CM.LM[top30,'gene'],offset= -.7)

dev.off()

#very more locis in CF.LF than in CM.LM, 
#pathways analysis
#lets see pathway enrichment in the windows for gene of locis with FC>30 and pval<10-4
#M.F
library(clusterProfiler)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

genes<-CM.CF$gene[CM.CF$pval<seuilpval & abs(CM.CF$FC)>seuilFC]
length(genes)
genes
enrichKEGG(genes,organism = )




#VISUALISATION RES
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

#COMPA MODEL
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

