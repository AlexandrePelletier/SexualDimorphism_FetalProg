#determiner les CpG différentiellement méthylé entre les sexes a 'linterieur d'un groupe

#config et library
options(stringsAsFactors=F)
library(data.table)
library(stringr)
library(plyr)
library(calibrate)
library(pheatmap)
library(limma)
set.seed(12345)
source("scripts/deterDataQual.R")
source("scripts/visual_function.R")
source("scripts/utils.R")
source("scripts/FindInfosGenes.R")
source("scripts/scRNA-seq integration.R")
#output dir 
script_name <- "eQTL_integration"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
date<-Sys.Date()
output <- file.path(outputDir,date)

#PARAMS
filtres<-"locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros"

#data
data_all<-fread("../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)

data_all<-data.frame(data_all)
rownames(data_all)<-data_all$id
head(data_all)

dim(data_all) #1709224     132
samples<-names(data_all)[str_detect(names(data_all),"CBP")]


batch<-read.csv2("../../ref/batch_CD34_library_date_090420.csv",header=T,row.names = 1)


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

#seq Depth var,
plot(batch$Library_Complexity,log10(batch$SeqDepth)) #better to put in log
batch$SeqDepthLog<-log10(batch$SeqDepth)


head(batch)


#data exploration
plot(density(na.omit(as.matrix(data_all[,samples]))))
plot(density(data_all$mean))
plot(density(data_all$sd))
plot(density(data_all$pct0))


#FILTRATION DES LOCIS
source("scripts/deter_seuilQC.R")
names(data_all)
#d'abord en fct msp1c et et nbNA
deterSeuilQC(data_all,metrique = "msp1c",qTestes = 1:9/40) #exclu locis < q0.125
quantile(data_all$msp1c,0.125) #6.532433e-08

deterSeuilQC(data_all,metrique = "pctNA",qTestes = 0:5,test = "brut") #seuil exclu ; no NA

mat<-as.matrix(data_all[,samples])
data_F<-data_all[data_all$msp1c>quantile(data_all$msp1c,0.125)&
                   rowSums(is.na(mat))==0,]

#puis on retire les locis full methylated
data_F$nbNonFullMethyl<-rowSums(data_F[,samples]>10,na.rm = T)

deterSeuilQC(data_F,metrique = "nbNonFullMethyl",qTestes = 0:10,test = "brut",lowerThan = F) #exclu locis avec nbNonFUllmethyl<5

data_F<-data_F[rowSums(data_F[,samples]>10)>4,]
nrow(data_F) #989522
#plus conf Score, nbMethylNonzeros dans pct0 elevé :

names(data_F)
deterSeuilQC(data_F,metrique = "confidenceScore",qTestes = 1:9/10) #exclu locis < q0.2
data_F<-data_F[data_F$confidenceScore>quantile(data_F$confidenceScore,0.2),]
nrow(data_F) #791613

deterSeuilQC(data_F[data_F$pct0>0.7,],metrique = "nbMethylNonZeros",qTestes = 0:5,test = "brut") #exclu locis avec pct0>0.7 et nbMethylVraizers==0


data_F<-data_F[!(data_F$pct0>0.7&data_F$nbMethylNonZeros==0),]
nrow(data_F) #786277
locisF<-rownames(data_F)
#gain en qualité
#avant
deterQual(mat) #41% des locis avec Vrais zeros
deterQual(mat[locisF,]) #54%>64% de locis avec vrais zeros

deterQual2(mat,batch) #"PC 1  ( 16.9 % de la variance) a R2 avec Group_Complexity = 0.75 et pval = 10^ -36.1671884642798"


deterQual2(mat[locisF,],batch) #"PC 1  ( 19.4 % de la variance a R2 avec Group_Complexity = 0.72 et pval = 10^ -33.9025713104284"

deterQual2(mat[!(rownames(mat)%in%locisF),],batch)
# "PC 1  ( 10.6 % de la variance a R2 avec Group_Complexity = 0.73 et pval = 10^ -34.8942186739751"

deterQual2(mat[sample(rownames(mat),length(locisF)),],batch)
#"PC 1  ( 16.9 % de la variance a R2 avec Group_Complexity = 0.75 et pval = 10^ -36.2830174521662"
deterQual2(mat[sample(rownames(mat),length(locisF)),],batch)
#PC 1  ( 16.9 % de la variance a R2 avec Group_Complexity = 0.75 et pval = 10^ -36.0886977688336"

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
sum(topVarLocis2 %in% topVarLocis1)/10000 #seulement 14% des topSD locis conservés !!!
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
png(paste(output,"heatmapsApresFiltration.png",sep="_"), width = 1200, height = 800)
p
dev.off()
#encore une clusterisation selon batch et Library Complexity(pas etonnant vu que 80% sont les meme locis)


#VISUALISATION DES SAMPLES
# en pct0

# batch[samples,"pct0"]<-colSums(data_all[,samples]==0,na.rm = T)/nrow(data_all)
 batch[samples,"pct0ApresF"]<-colSums(data_F[,samples]==0,na.rm = T)/nrow(data_F)
 #batch[samples,"pct0locisLo"]<-colSums(data_all[!(rownames(data_all)%in%locisF),samples]==0,na.rm = T)/sum(!(rownames(data_all)%in%locisF))
 write.csv2(batch,"../../ref/batch_CD34_library_date_090420.csv",row.names = T)

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
table(batch$Group_name[batch$pct0ApresF>0.8])
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
PCAlist<-readRDS("analyses/main/PCAlist.rds")
# pc2<-prcomp(t(as.matrix(data_F[,samples])),center = T)
# PCAlist[["pca_F"]]<-pc2
# rm(pc2)
#visual PCA et INFLUENCE COVAR SUR PC
pcaChoose<-"pca_All"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
names(batch)
plot(density(batch$SeqDepth))
plot(batch$SeqDepth)
abline(h=quantile(batch$SeqDepth,0.91))

plot(batch$SeqDepth,batch$newComplexityScore2,col=as.numeric(batch$SeqDepth>2*10^7)+1)
batch$SeqDepthHi<-as.numeric(batch$SeqDepth>2*10^7)

plot(batch$SeqDepth,batch$newComplexityScore4,col=batch$Group+1)

plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="SeqDepthHi",batch = batch,showSampleIDs=F)

plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group",batch = batch,showSampleIDs=F)

plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Sex",batch = batch,showSampleIDs=F)

var_fac<-names(batch)[c(2,4,8,10,11,12,13,15,16,19,24,27,29,33,34,35,41)]
var_num<-names(batch)[c(5,6,17,20,21,22,23,25,26,28,32,37,39,42)]
varAdd<-c('GroupBatch_Complexity','GroupBatch_Complexity_Fac','pct0ApresF')
vardint<-c("Group","Group_Sex","Sex","Library_Complexity","Group_Complexity","Group_Complexity_Fac","newComplexityScore2","SeqDepthLog" )
resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd ) #date PC4 et 6




rowSums(-log10(resPV[rownames(resPV)%in%vardint,]))
# Group_Complexity   Library_Complexity                Group Group_Complexity_Fac            Group_Sex                  Sex 
# 40.745262            40.975445             4.269812            32.048345             5.193760             2.930691 

pcaChoose<-"pca_F"
PCs1pct<-plotPCVarExplain(PCAlist[[pcaChoose]],1:40,lineSeuilPct = 1)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="Group_Complexity_Fac",showSampleIDs=F)
plotPCA(PCAlist[[pcaChoose]],PCx=1,PCy=2,colorBatch="batch",batch = batch,showSampleIDs=F)

resPV<-plotCovarPCs(PCAlist[[pcaChoose]],PCs1pct,batch,var_num,var_fac,exclude = varAdd )
rowSums(-log10(resPV[rownames(resPV)%in%vardint,]))
# Group_Complexity   Library_Complexity                Group Group_Complexity_Fac            Group_Sex                  Sex 
# 38.053873            38.754115             7.637016            30.375926             5.377435             2.859952 


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
# [1] "Group_Complexity"     "Library_Complexity"   "Mat.Age"              "newComplexityScore2"  "SeqDepth"            
# [6] "SeqDepthLog"          "batch"                "date"                 "Group"                "Group_Complexity_Fac"
# [11] "latino"               "SeqDepthHi"           "sequencing"   



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
seqDepthLog<-as.numeric(batch[samples,"SeqDepthLog"])
seqDepth<-as.numeric(batch[samples,"SeqDepth"])
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
seqDepthHi<-as.factor(batch[samples,"SeqDepthHi"])
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
correl(group_complexity,group_sex,ret = "all") #no signif
correl(group_complexity,sex,ret = "all") #no signif
correl(group_complexity,weight,ret = "all")
correl(group_complexity,weightAtTerm,ret = "all")
correl(group_complexity,mat.age,ret = "all") #encore signif mais moins r2= 0.06

#group_comp_fac ?
correl(weightAtTerm,group_complexity_fac,ret="all")
correl(mat.age,group_complexity_fac,ret="all")# c'est pu signif ! donc c'est mieux


#correl entre groupe_sex et (group_complexity_fac,sequencing, mat.age, latino,batch)
correl(group_sex,sequencing,ret="all") #no

correl(mat.age,group_sex,ret="all") #yes, r2 0.09

correl(latino,group_sex,ret="all") #no
correl(batches,group_sex,ret="all") #no

correl(batches,group_sex,ret="all") #no

#mat.age correler avec group, le prendre dans model ou pas ?
pc<-PCAlist[[2]]$x
summary(lm(pc[,1]~group_complexity_fac)) #r2=0.707
summary(lm(pc[,1]~group_complexity_fac+group_sex)) #r2 = 0.78
summary(lm(pc[,1]~group_complexity_fac+group_sex+mat.age)) #r2 0.773 donc vaut mieux pas mettre mat.age
summary(lm(pc[,1]~group_complexity_fac+group_sex+latino)) #0.786 !
summary(lm(pc[,1]~group_complexity_fac+group_sex+latino+seqDepthLog)) #nameliore pas explication pc1
summary(lm(pc[,1]~group_complexity_fac+group_sex+latino+seqDepthHi)) #ameliore un peu : 0.79

#group_comp_fac et seqDepHi utile les 2 ?
correl(seqDepthHi,group_complexity_fac,ret="all") #pas independant, donc group_complexity_fac suffisant



correl(mat.age,group_complexity_fac)


summary(lm(pc[,1]~group_complexity_fac+group_sex+mat.age+latino))  #0.79 !
summary(lm(pc[,1]~group_complexity_fac+group+latino))  #0.796 c'est mieux.
summary(lm(pc[,1]~group_complexity_fac+group_sex+latino+seqDepthLog)) #0.82 ! 

summary(lm(pc[,1]~group_complexity_fac+group+latino+seqDepthLog+mat.age)) #0.83 !

#ccl mat.age peu nécessaire pour expliquer pc1, mais peut etre bon pour expliquer PC6
# on capture avec group de la var que complexity library n'expliquait seul

#le model opti serait donc :  model ~0 + group_sex + sequencing + batches + group_complexity_fac + latino + mat.age


### histoire de la colinearité entre group et mat.age
plot(batch$Mat.Age,batch$pct0ApresF,col=(batch$Group_Complexity_Fac)) 
# rownames(batch)[batch$Group_Complexity_Fac==4&batch$pct0ApresF>0.8]
# batch$Group_Complexity_Fac[batch$Group_Complexity_Fac==4&batch$pct0ApresF>0.8]<-1
batch2<-na.omit(batch[,c("Mat.Age","pct0ApresF")])

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

models<-list()
model<-13
names(batch)
varToModel<-c("Group_Sex",'batch',"Mat.Age","latino","Group_Complexity_Fac")
samples_F_F<-samples[rowSums(is.na(batch[samples,varToModel]))==0] 
length(samples_F_F) # 108 (sans sequencing)
table(batch[samples_F_F,"Group_name"])
# C  I  L 
# 34 38 36 
sequencing<-factor(batch[samples_F_F,"sequencing"])
#group<-factor(batch[samples_F_F,"Group_name"])
#groupBatch_complexity_fac<-factor(batch[samples_F_F,"GroupBatch_Complexity_Fac"])
#sex<-as.factor(batch[samples_F_F,"Gender"])
mat.age<-as.numeric(batch[samples_F_F,"Mat.Age"])

batches<-factor(batch[samples_F_F,"batch"])
group_sex<-factor(batch[samples_F_F,"Group_Sex"])
group_sex<-revalue(group_sex,c("1"="FC","2"="MC","3"="FI","4"="MI","5"="FL","6"="ML"))
#group_complexity<-as.numeric(batch[samples_F_F,"Group_Complexity"])
group_complexity_fac<-factor(batch[samples_F_F,"Group_Complexity_Fac"])
latino<-factor(batch[samples_F_F,"latino"])
#complexity<-as.numeric(batch[samples_F_F,"Library_Complexity"])


formule<- ~0 + group_sex  + batches  + latino + mat.age + group_complexity_fac
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
#saveRDS(fit2,"analyses/main/fit2_model13.rds")
results <- decideTests(fit2)

sum(abs(results)) #2039 > 3755 !
colSums(abs(results))

#4:
# C.I   C.L   I.L MC.ML MC.MI MI.ML FC.FL FC.FI ML.FL MI.FI MC.FC   F.M 
# 0    71  1955     0     0     0     4     0     7     2     0     0 

# C.I   C.L   I.L MC.ML MC.MI MI.ML FC.FL FC.FI ML.FL MI.FI MC.FC   F.M 
# 0   126  3602     0     1     0    17     0     8     1     0     0 

#13 : 
#bon model ? 1) enrichissment en enh et prom, 2) prox du gene
locisSig<-rownames(fit2$p.value)[apply(fit2$p.value<0.001,1,any)]
length(locisSig) # 27k (4) > 23k
resSig<-data.frame(row.names = locisSig,fit2$p.value[locisSig,],annot[locisSig,c("chr","start","posAvant","gene","type")],
                   data_F[locisSig,c("confidenceScore","confidenceScoreNorm","complexity","msp1c","RankConfidenceScore")])


head(resSig,100)
#nb C-L 
sum(resSig$C.L<0.001)  #5069>5008

#enrichissement en bon locis :
mean(resSig$RankConfidenceScore)/mean(na.omit(data_all$RankConfidenceScore)) #1.40>1.49 !

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


#4:
# 0         1         2         3         4         5         6 
# 13.877926  6.417870  5.076823  5.670715 23.893636  4.915897 39.840607 

#13 :
# 0         1         2         3         4         5         6 
# 9.529833  5.443812  3.930002  4.516129 25.738984  4.237824 46.341978 


#locis around5k of the gene
locis5k<-na.omit(rownames(annot)[annot$posAvant<5000&annot$posAvant>-5000])

sum(locisRall%in%locis5k)/length(locisRall) # 0.26

sum(locisRF%in%locis5k)/length(locisRF) #0.32

sum(locisSig%in%locis5k)/length(locisSig)  #0.54 > 0.61

hist(annot[locisRall[locisRall%in%locis5k],"posAvant"],breaks = 100)
hist(annot[locisRF[locisRF%in%locis5k],"posAvant"],breaks = 100)
hist(annot[locisSig[locisSig%in%locis5k],"posAvant"],breaks = 100)


#MODEL De confiance : 13

#concentrons nous sur le top1000 par compa
compas<-colnames(fit2$p.value)

#deter seuil pval pour avoir assez de locis signif
for(compa in compas){
  print(compa)
  top1000<-topTable(fit2,coef = compa, adjust="BH",number = 1000)
  #remove F si pval <0.001
  print(sum(top1000$P.Value<0.001))
 
  
}

resParCompa<-list()
#add quantile msp1c et conf plutot qua valeur brut (car illisible)
q9Msp1c<-quantile(data_F$msp1c,1:9/10)
q9ConfScore<-quantile(data_F$confidenceScore,1:9/10)

for(compa in compas){
  print(compa)
  res<- topTable(fit2,coef=compa,n =Inf)
  
  res<-data.table(locisID=as.numeric(rownames(res)),
                  FC=res$logFC,
                  pval=res$P.Value,
                  pval.adj=res$adj.P.Val,
                  complexity=data_F[rownames(res),"complexity"],
                  msp1Conf=sapply(data_F[rownames(res),"msp1c"],function(x){
                    return(sum(x>q9Msp1c))
                  }),
                  confScore=sapply(data_F[rownames(res),"confidenceScore"],function(x){
                    return(sum(x>q9ConfScore))
                  }))
  # save all locis pval and FC : 
  
  fwrite(res,paste(output,"res_locis_in",compa,"allLocis",filtres,"model",model,".csv",sep = "_"),sep=";")
  resParCompa[[compa]]<-res
  print("ok")
  
}
#now, annot res with genes :
#  annot1<-fread("../../ref/annotation_CpG_HELP_ALL_070420.csv")
#  annot1<-annot1[,distTSS:=posAvant][,-c("ENSEMBL_ID","posAvant", "ENTREZID","HSPC1","HSPC2",
#                                   "CTCF.Binding.Site","Promoter", "Promoter.Flanking.Region", "Open.chromatin","Enhancer",
#                                   "TF.binding.site", "CTRL_sc",
#                                   "LGA_sc", "eQTLScore")]
# annot1
# annot2<-fread("../../ref/allCpG_annotation_genes_links_with_blood_eQTL_and_ensembl_reg_dom.csv")
# annot2<-annot2[,pos:=pos-1][!is.na(start.eQTR2)][,gene:=geneLink][,-c("geneLink","length.eQTR2","nQTL.reg.gene","nGenes.reg","nGenes.cpg","nCpg.gene")]
# annot2
# annot<-merge(annot1,annot2,all=TRUE)
# annot
# fwrite(annot,"../../ref/allCpG_annotation_genes_and_feature_regulatory_domain_070520.csv",sep = ";")

annot<-fread("../../ref/allCpG_annotation_genes_and_feature_regulatory_domain_070520.csv")

for(compa in compas){
  res<-resParCompa[[compa]]
  res<-merge(res,annot,by="locisID",all.x = T)
  resParCompa[[compa]]<-res[order(pval)]
  
  #
  #fwrite(res,paste(output,"res_locis_in",compa,"allLocis",filtres,"model",model,".csv",sep = "_"),sep=";")
}
##SCORE CpG
#need eQTR pour confWeight Links CpG-Gene
eQTRs<-fread("../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv")
#add eQTR1 et pval et 
eQTRs<-eQTRs[,.(chr,start.eQTR2,end.eQTR2,start.eQTR,end.eQTR,avg.pval.thr)]
eQTRs
res<-merge(res,eQTRs,all.x = TRUE)[order(pval)]
res<-res[,avg.pval.thr:=abs(avg.pval.thr)]

#need also a function to calculate confidence of the links CpG-Gene for eQTL linked Gene : 
calcLinkScore<-function(pos,start,end,signif,maxSignif=0.001475078){
  pos<-pos[1]
  start<-start[1]
  end<-end[1]
  signif<-signif[1]
  
  if(pos > (start-50)& pos < (end+50)){
    score<-1
  }else {
    score<-log10(50)/(log10(1+abs(pos-min(abs(c(pos-(start-50),pos-(end+50)))))))
  }
  #si min pval >*1 si max pval > *0.75
  score<-score*(0.75+0.25*(signif/maxSignif))
  
  #ajutser pour les eQTR trop grands (donc trop incertain) : si >200, ajuste score sinon descend progressivement
  largeur<-end-start+1
  if(largeur>200){
    score<-score*(0.75+0.25*log10(200)/log10(largeur+1))
  }
  
  
  return(score)
  
}


for(compa in compas){
  res<-resParCompa[[compa]]
  #FILTER RES
  #without type and without gene
  nLocisAvant<-length(unique(res$locisID))
  res<-res[!is.na(type)][!is.na(gene)&gene!=""]
  nLocisApres<-length(unique(res$locisID))
  print(paste(round((nLocisAvant-nLocisApres)/nLocisAvant*100,1),"% des locis non associés a des gènes"))
  ##DMC weigtht
  #(FC[0:1] : toavoid overfitting, we use the rank of abs(FC) as metriks normalized
  res[,FCScore:=rank(abs(FC))/max(rank(abs(FC)))]

  #pval[0:1] : the same (with rank), to have the same weight
  res[,PvalScore:=rank(1/pval)/max(rank(1/pval))]
  
  res[,DMCWeight:=PvalScore+FCScore]
  
  ##regulatory weight : je veux un poids entre 0.5 et 2, 0.5 etant aucune annotation regulatrice,
  #et 2 etant toutes les annotations regulatrice (4/6, CTCF/prom/enh/, TF motif) 
  #=> on part d'un score de 0.5 et +0.5 a chaque annot
  
  # +1 if (4,6), +0.75 if 5, +0.5 otherwise
  res[,TypeScore:=sapply(type,function(x){
    if(x%in%c(4,6)){
      score<-1
    }else if(x==5){
      score<-0.75
    }else{
      score<-0.5
    }
    return(score)
  })]
  
  # +0.5 if(CTCF,enh,prom,)+ 0.25 if prom flanking, openchrine, #+0.5 if(TFbind site), or 0 otherwise
  res[,EnsRegScore:=sapply(feature_type_name,function(x){
    vecX<-tr(x)
    if(any(c("CTCF Binding Site","Promoter","Enhancer")%in%vecX)){
      score<-0.5
    }else if(any(c("Open chromatin","Promoter Flanking Region")%in%vecX)){
      score<-0.25
    }else{
      score<-0
    }
    if("TF binding site"%in%vecX){
      score<-score+0.5
    }
    return(score)
  })]
  
  # add the 2 score to make the regulWeight:
  res[,RegWeight:=TypeScore+EnsRegScore]
  
  #ConfWeight :
  # confWeight [0-1] : TSS, eQTL links, (cross correl ?)
  #rank normalisé de sum de : 
  #pour les genes pas asso par eQTR : tss +1 si <1000pb, descend progressivement sinon  : 
  res[is.na(start.eQTR2),LinkScore:=sapply(abs(distTSS),function(x){
    if(x<1000){
      score<-1
    }else {
      score<-3/log10(x)
    }
    return(score)
  })]
  
  #pour les genes asso par eQTR : 
  # 1 si dans eQTR1+/-50pb, baisse progress jusqu'a 0.8
  
  res[!is.na(start.eQTR2),LinkScore:=calcLinkScore(pos,start.eQTR,end.eQTR,avg.pval.thr),by=.(locisID,start.eQTR2,gene)]
  
  #confWeight : rank
  
  
  
  #fwrite(res,paste(output,"res_locis_in",compa,"allLocis",filtres,"model",model,".csv",sep = "_"),sep=";")
}




#FIGURES RES LIMMA
library(calibrate)
seuilFC<-30
seuilpval<-10^-4
for(compa in compas){
  print(compa)
  res<-resParCompa[[compa]]
  colors<-as.numeric(res$pval<seuilpval & abs(res$FC)>seuilFC)+1
  
  down<-sum(res$pval<seuilpval & res$FC<(-seuilFC))
  up<-sum(res$pval<seuilpval & res$FC>(+seuilFC)) 
  
  png(paste(output,compa,"vulcano_model",model,".png",sep="_"), width = 700, height = 700)
  
  plot(res$FC,-log10(res$pval),col=colors,main = compa,sub=paste0(nrow(res)," locis, with",up+down," locis sig. : ",down,"hypoM and",up," hyperM"),pch=18)
  abline(h=-log10(seuilpval))
  abline(v=seuilFC)
  abline(v=-seuilFC)
  #textxy(C.L[top30,'FC'],-log10(C.L[top30,'pval']),C.L[top30,'gene'],offset= -.7)
  
  dev.off()
}


#RES PAR genes /!\ a actualiser
#annot expr genes with CD34+ rnaseq data : 

#SCORE GENE


###PATHWAY ANALYSIS
# POUR model 4 in the windows for gene of locis with FC>30 and pval<10-4
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)



###DMR ANALYSIS : 
#quelle algo prendre ?
#from An evaluation of supervised methods for identifying differentially methylated regions in Illumina methylation arrays
#(https://academic.oup.com/bib/article/20/6/2224/5096828#191593271) : 
#comb-p showed the best sensitivity as well as good control of FP rate across all simulation scenarios
#comb-p identified substantially more TP DMRs than other methods, especially when effect sizes were small

#=> COMB-P
#command-line tool and a Python library [16].
#principe : corrected P-value at a CpG site will be smaller than the original P-value if the neighboring CpG sites also have comparatively small P-values
#input of comb-p is a .BED file with P-values and chromosome locations of the CpG sites
for (compa in compas){
  
  res<-topTable(fit2,coef = compa,n=Inf)
  head(res) 
  #need add chr start et stop
  res<-data.frame(row.names =rownames(res),annot[rownames(res),c("chr","start","stop")],pval=res$P.Value )
  #need rm NA : 
  dim(res)
  res2<-na.omit(res)
  head(res2)
  #need sort by genomic location
  res2<-arrange(res2,chr,start) #tjr pas bon, car commence apres chr1 il y a chr10
  #revalue levels chr
  res2$chrom
  tail(res2)
  dim(res2)
  
  #need chrom, start, end header
  colnames(res2)<-c("chrom","start","end","p-value")
  
  write.table(res2,file = paste0("analyses/DMR_with_comb_p/",compa,"model",model,"_allLocis_and_pval.bed"),row.names = F,sep = "\t",quote = F)
  
}


#in shell (dir : analyses/DMR_with_comb_p) :
# comb-p pipeline \
# -c 4 \          # p-values in 4th column
# --seed 1e-3 \   # require a p-value of 1e-3 to start a region
# --dist 200      # extend region if find another p-value within this dist
# -p FC.FL \
# --region-filter-p 0.1 \ # post-filter reported regions
# --anno hg38 \            # annotate with genome hg38 from UCSC
# FC.FL_allLocis_and_pval.bed                  # sorted BED file with pvals in 4th column


#with all compas : 
# for compa in 'C.L' 'MC.ML' 'FC.FL' 'ML.FL' 'MC.FC' 'F.M'; do comb-p pipeline -c 4 --seed 0.05 --dist 750 -p "model14""$compa" --region-filter-p 0.05 "$compa""model14_allLocis_and_pval.bed"; done 

#all compas get :
resDMRs<-list()

for(compa in compas){
  print(compa)
  resDMR<-read.table(paste0("analyses/DMR_with_comb_p/model13",compa,".regions-t.bed"),sep = "\t")
  colnames(resDMR)<-c("chrom","start","end","min_p","n_probes","z_p","z_sidak_p")
  print(nrow(resDMR)) #108 region
  #annotation :
  for(chrom in levels(as.factor(resDMR$chrom))){
    
    for (region in which(resDMR$chrom==chrom)){
      
      range<-c(resDMR[region,"start"],resDMR[region,"end"])
      genes<-annot$gene[which(annot$chr==chrom&annot$start>range[1]&annot$start<range[2])]
      
      resDMR[region,"gene"]<-paste(unique(genes),collapse = "/")
      resDMR[region,"nbLocis"]<-paste(table(genes),collapse = "/")
    }
  }
  write.csv2(resDMR,file = paste(output,compa,"resDMR_comb-p_model",model,".csv",sep = "_"))
  resDMRs[[compa]]<-resDMR
}


names(resDMRs)
for(compa in "F.M"){
  print(compa)
  resDMR<-read.table(paste0("analyses/DMR_with_comb_p/model13",compa,".regions-t.bed"),sep = "\t")
  colnames(resDMR)<-c("chrom","start","end","min_p","n_probes","z_p","z_sidak_p")
  print(nrow(resDMR)) #108 region
  #annotation :
  for(chrom in levels(as.factor(resDMR$chrom))){
    
    for (region in which(resDMR$chrom==chrom)){
      
      range<-c(resDMR[region,"start"],resDMR[region,"end"])
      genes<-annot$gene[which(annot$chr==chrom&annot$start>range[1]&annot$start<range[2])]
      
      resDMR[region,"gene"]<-paste(unique(genes),collapse = "/")
      resDMR[region,"nbLocis"]<-paste(table(genes),collapse = "/")
    }
  }
  write.csv2(resDMR,file = paste(output,compa,"resDMR_comb-p_model",model,".csv",sep = "_"))
  resDMRs[[compa]]<-resDMR
}
#du sens bio ? 
library(clusterProfiler)
resDMR<-resDMRs[["FC.FL"]]
gene.df<-bitr(unique(resDMR$gene),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene.df
kkDMR<-enrichKEGG(gene.df$ENTREZID,pvalueCutoff = 0.05, qvalueCutoff = 0.15)
head(kkDMR) #signaling pluripotency, gastic cancer et wnt sig

GO_DMR <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable = T)

head(GO_DMR) #transcription corepressor activity 




#BONUS : scRNA-seq annot

# #reannot 
# library(Seurat)
# source("scripts/scoreCluster.R")
# samples<-readRDS("../../../Alexandre_SC/analyses/test_batch_integration/CTRLandLGA_SCTransform_Integrated.rds")
# new.cluster.ids <- c("HSC-SELL2", "HSC-AVP", "HSC-Er", "HSC-SELL1", "EMP", "GMP",
#                      "LyP", "MkP", "proT","LMPP","preB","LT-HSC","Neu","LyB","Ly-ETS1","DC")
# names(new.cluster.ids) <- levels(samples)
# samples <- RenameIdents(samples, new.cluster.ids)
# 
# for(i in 1:length(new.cluster.ids)){
#   markers$cluster[markers$cluster==as.numeric(names(new.cluster.ids)[i])]<-new.cluster.ids[i]
# }
# head(markers)
# scoresCluster<-scoreMarquageCluster(markers,samples,seuil = "intraClusterFixe",filtreMin = 2)
# write.csv2(scoresCluster,"analyses/test_scRNAseq_integration/sc_scoreClustersHSPC.csv",row.names = T)

markers<-read.csv2("../../../Alexandre_SC/analyses/test_batch_integration/all.markers_SCTransform_CTRLandLGA_integrated.csv",row.names = 1)
head(markers)

scoresCluster<-read.csv2("analyses/test_scRNAseq_integration/sc_scoreClustersHSPC.csv",row.names = 1)
head(scoresCluster)

resGenesParCompa<-list()

#for genes Expr
CTRL_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInCTRL.csv",row.names = 1)
LGA_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInLGA.csv",row.names = 1)
listExprMat<-list(CTRL=CTRL_expr_by_pop,LGA=LGA_expr_by_pop)

for(compa in compas){
  print(compa)
  #RES PAR GENES
  #make a df with rownames = genes, and  count of locis by genes 
  resGenes<-resLocisToGenes(resParCompa[[compa]]) #!function a remplacer par un datatable
  
  #add genes Expr in bulk scrnaseq and in subpop
  resGenes<-find.scGeneExpr(resGenes,listExprMat)
  #add markers of clusters
  resGenes<-addMarqueursClusters(scoresCluster,resGenes)
  #add canonical/published markers of cell types
  resGenes<-findMarqueursPops(rownames(resGenes),df = resGenes)
  
  #save 
  resGenesParCompa[[compa]]<-resGenes
  
  
}

