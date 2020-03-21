#config et library
options(stringsAsFactors=F)
library(data.table)
library(stringr)
set.seed(12345)

script_name <- "deter_QC_filter"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")

#data
data_all<-fread("../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)
data_all=data.frame(data_all)
rownames(data_all)<-data_all$id
head(data_all)
dim(data_all) #1709224     127
samples<-names(data_all)[str_detect(names(data_all),"CBP")]

batch=read.csv2("../../ref/20-02-17_batch_119samples_IUGRLGACTRL_noNA_withHpaC.csv",header=T,row.names = 1)
#je veux enlever les locis de mauvaise qual selon un critére bien définit
#ccl des 2 premiers tests : methyl entre 0 et 20, 
#bon locis = pctdescore entre 0 et 20 élévé ou pct0 faible (<0.7) (ou sd > 25 (locis ++ unmethylé), mais pas utilisé pour le moment)
#qualitéData=pct0*pctEntre0et20

#seuil optimale = maximiser ses 3 metriques avec msp1c, confscore, nbVraisZeros, pct0 

#comment : chercher seuil bas (d'exclusion) et haut (d'inclusion)
#1)seuil bas : pour quantile de 0.1 a 0.9, tant que locis < quantile x n'augmente pas pctBonlocis augmenter

#1.1 plot coude : 
#a) fct deterQual
source('scripts/old/deterQual.R')

#b) test arrondi
mat<-as.matrix(data_all[,samples])
deterQual(mat) #0.00938


#c) plot coude de 0.01 a 0.20 de msp1c

deterSeuilQC<-function(dataCpG,metrique,qTestes,test="quantile"){
  quals<-rep(0,length(qTestes))
  print(paste("determine si",metrique,"peut permettre d'exclure des locis de faible confiance"))
  for (i in 1:length(qTestes)){
    q<-qTestes[i]
    print(paste("seuil",q))
    if(test=="quantile"){
      q<-quantile(dataCpG[,metrique],q,na.rm=T)
      
    }
    
    locisLo<-na.exclude(rownames(dataCpG)[dataCpG[,metrique]<=q])
    print(paste("analyses qualité des",length(locisLo),"locis avec",metrique,"<",round(q,2)))
    
    if(length(locisLo)>100000){
      locisLo<-sample(locisLo,100000)
    }
    dataCpGLo<-as.matrix(dataCpG[locisLo,samples])
    
    
    quals[i]<-deterQual(dataCpGLo,maxMethyl = 5)
  }
  
  print(plot(qTestes,quals))
  return(quals)
}
#on veut enlelver les locis faux zeros parmi les locis pct0>0.7
locis0<-rownames(data_all)[data_all$pct0>0.7]
#tester, msp1c, conf
names(data_all)

#msp1c
qTestes<-1:20/50
qualsMsp1c<-deterSeuilQC(data_all,"msp1c",qTestes) #0.08

#confScore
qTestes<-1:20/50
qualsConf<-deterSeuilQC(data_all,"confidenceScore",qTestes) #0.08 aussi

#mean
qTestes<-1:19/40
qualsMean<-deterSeuilQC(data_all,"mean",qTestes) #quals diminue jusqu'a 0.2 puis augmente ! donc on enleve locis avec mean <0.6 !
plot(density(data_all$mean))
abline(v=0.1)
#entre q0.1 et 0.2
locis0IntMean<-locis0[data_all[locis0,"mean"]>quantile(data_all$mean,0.1)&data_all[locis0,"mean"]<quantile(data_all$mean,0.2)]
length(locis0IntMean)/nrow(data_all)#10%
locis0ExtMean<-locis0[!(locis0%in%locis0IntMean)]
locis0LoMean<-locis0[data_all[locis0,"mean"]<quantile(data_all$mean,0.2)]
locis0HiMean<-locis0[data_all[locis0,"mean"]>quantile(data_all$mean,0.2)]

#locislomean bad qual locis ?
plot(density(log10(data_all$msp1c)))
lines(density(na.omit(log10(data_all[locis0,"msp1c"]))),col=2)
lines(density(na.omit(log10(data_all[locis0IntMean,"msp1c"]))),col=3)
lines(density(na.omit(log10(data_all[locis0ExtMean,"msp1c"]))),col=4)
lines(density(na.omit(log10(data_all[locis0LoMean,"msp1c"]))),col=5)
lines(density(na.omit(log10(data_all[locis0HiMean,"msp1c"]))),col=6)

plot(density(na.omit(log10(data_all$confidenceScore))))
lines(density(na.omit(log10(data_all[locis0,"confidenceScore"]))),col=2)
lines(density(na.omit(log10(data_all[locis0IntMean,"confidenceScore"]))),col=3)#bad locis
lines(density(na.omit(log10(data_all[locis0ExtMean,"confidenceScore"]))),col=4)#best locis
lines(density(na.omit(log10(data_all[locis0LoMean,"confidenceScore"]))),col=5)
lines(density(na.omit(log10(data_all[locis0HiMean,"confidenceScore"]))),col=6)
#ce serait effectivmeny les locis int les mauvais
#check avec library ~PC1
#all data
mat0<-mat
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8436.5 -3856.7  -757.4  2651.5 29657.6 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 21935.13    1240.34   17.68   <2e-16 ***
#   library     -1179.92      61.22  -19.27   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5377 on 117 degrees of freedom
# Multiple R-squared:  0.7605,	Adjusted R-squared:  0.7584 
# F-statistic: 371.4 on 1 and 117 DF,  p-value: < 2.2e-16

#juste pct0>0.7
mat0<-mat[locis0,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5086.2 -1353.8  -567.5   499.8 13327.5 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -5775.89     605.54  -9.538 2.67e-16 ***
#   library       310.69      29.89  10.395  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2625 on 117 degrees of freedom
# Multiple R-squared:  0.4801,	Adjusted R-squared:  0.4757 
# F-statistic:   108 on 1 and 117 DF,  p-value: < 2.2e-16

#juste locisExt
mat0<-mat[locis0ExtMean,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
summary(pca)
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -7156.5  -408.5   480.6  1010.7  4351.9 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5071.70     402.43   12.60   <2e-16 ***
#   library      -272.81      19.86  -13.73   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1745 on 117 degrees of freedom
# Multiple R-squared:  0.6172,	Adjusted R-squared:  0.6139 
# F-statistic: 188.6 on 1 and 117 DF,  p-value: < 2.2e-16


#juste locisInt
mat0<-mat[locis0IntMean,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
summary(pca)
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)

# Residuals:
# Min     1Q Median     3Q    Max 
# -488.8  -50.5  -40.8  -31.0 5527.4 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -49.757    118.860  -0.419    0.676
# library        2.677      5.867   0.456    0.649
# 
# Residual standard error: 515.3 on 117 degrees of freedom
# Multiple R-squared:  0.001776,	Adjusted R-squared:  -0.006756 
# F-statistic: 0.2081 on 1 and 117 DF,  p-value: 0.6491

#! tester si juste Mean seuil haut suffit.pour gagner en locis
#juste locisLo
mat0<-mat[locis0LoMean,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
summary(pca) #explique que 3.7% !
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)


#ccl : locis avec intMean sont bon, bizattre, ptete cachéds pc2

#juste locishi
mat0<-mat[locis0HiMean,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
summary(pca) #explique que 3.7% !
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4283.5 -1005.4  -394.8   335.0  6024.0 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -4925.73     373.25  -13.20   <2e-16 ***
#   library       264.96      18.42   14.38   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1618 on 117 degrees of freedom
# Multiple R-squared:  0.6387,	Adjusted R-squared:  0.6356 
# F-statistic: 206.8 on 1 and 117 DF,  p-value: < 2.2e-16
# #ccl : locis avec HiMean sont mauvais, on doit les enlever

#sd
qTestes<-1:19/20
qualsSD<-deterSeuilQC(data_all,"sd",qTestes) 
#quals diminue jusqu'a 0.4, plateau jusqu'a 0.7 , et redescend ensuite : 
plot(density(data_all$sd))
lines(density(data_all[locis0,"sd"]))
abline(v=quantile(data_all$sd,0.4))
abline(v=quantile(data_all$sd,0.7))


#on enleve sd <0.4 et > 0.7 ? ça garde seulement 30% 
#♥mais pas si inclusion des locis>70% des de pct0 !

locis0IntSD<-locis0[(data_all[locis0,"sd"]>10.33)&(data_all[locis0,"sd"]<21.01)]
locis0ExtSD<-locis0[!(locis0%in%locis0IntSD)]
locis0LoSD<-locis0[(data_all[locis0,"sd"]<10.33)]
locis0HiSD<-locis0[(data_all[locis0,"sd"]>21.01)]


length(locis0ExtSD)/nrow(data_all) #45%

plot(density(log10(data_all$msp1c)))

lines(density(na.omit(log10(data_all[locis0,"msp1c"]))),col=2)
lines(density(na.omit(log10(data_all[locis0IntSD,"msp1c"]))),col=3)
lines(density(na.omit(log10(data_all[locis0ExtSD,"msp1c"]))),col=4)
lines(density(na.omit(log10(data_all[locis0LoSD,"msp1c"]))),col=5)
lines(density(na.omit(log10(data_all[locis0HiSD,"msp1c"]))),col=6)

plot(density(na.omit(log10(data_all$confidenceScore))))
lines(density(na.omit(log10(data_all[locis0,"confidenceScore"]))),col=2)
lines(density(na.omit(log10(data_all[locis0IntSD,"confidenceScore"]))),col=3)
lines(density(na.omit(log10(data_all[locis0ExtSD,"confidenceScore"]))),col=4)
lines(density(na.omit(log10(data_all[locis0LoSD,"confidenceScore"]))),col=5)
lines(density(na.omit(log10(data_all[locis0HiSD,"confidenceScore"]))),col=6)
#a voir si c'est vraiment des bon locis avec dist, PC1~Library
#all data
mat0<-mat
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8436.5 -3856.7  -757.4  2651.5 29657.6 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 21935.13    1240.34   17.68   <2e-16 ***
#   library     -1179.92      61.22  -19.27   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5377 on 117 degrees of freedom
# Multiple R-squared:  0.7605,	Adjusted R-squared:  0.7584 
# F-statistic: 371.4 on 1 and 117 DF,  p-value: < 2.2e-16

#juste pct0>0.7
mat0<-mat[locis0,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5086.2 -1353.8  -567.5   499.8 13327.5 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -5775.89     605.54  -9.538 2.67e-16 ***
#   library       310.69      29.89  10.395  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2625 on 117 degrees of freedom
# Multiple R-squared:  0.4801,	Adjusted R-squared:  0.4757 
# F-statistic:   108 on 1 and 117 DF,  p-value: < 2.2e-16

#juste locisExt
mat0<-mat[locis0ExtSD,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -7156.5  -408.5   480.6  1010.7  4351.9 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5071.70     402.43   12.60   <2e-16 ***
#   library      -272.81      19.86  -13.73   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1745 on 117 degrees of freedom
# Multiple R-squared:  0.6172,	Adjusted R-squared:  0.6139 
# F-statistic: 188.6 on 1 and 117 DF,  p-value: < 2.2e-16


#juste locisInt
mat0<-mat[locis0IntSD,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2765.3  -949.8  -367.1   296.8 11771.4 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2987.65     478.40  -6.245 7.06e-09 ***
#   library       160.71      23.61   6.806 4.55e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2074 on 117 degrees of freedom
# Multiple R-squared:  0.2836,	Adjusted R-squared:  0.2775 
# F-statistic: 46.32 on 1 and 117 DF,  p-value: 4.554e-10

#! tester si juste sd seuil haut suffit.pour gagner en locis
#juste locisLo
mat0<-mat[locis0LoSD,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
summary(pca) #explique que 3.7% !
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
summary(resLm)

# Residuals:
#   Min     1Q Median     3Q    Max 
# -702.2 -248.3 -146.7    9.7 9158.3 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -769.92     226.49  -3.399 0.000924 ***
#   library        41.42      11.18   3.704 0.000325 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 981.9 on 117 degrees of freedom
# Multiple R-squared:  0.105,	Adjusted R-squared:  0.09733 
# F-statistic: 13.72 on 1 and 117 DF,  p-value: 0.0003248

#ccl : locis avec loSD sont bon, on doit les garder

#juste locisLo
mat0<-mat[locis0HiSD,]
mat0[is.na(mat0)]<-0
pca<-prcomp(t(mat0))
pc<-pca$x
summary(pca) #explique que 3.7% !
#library<-as.numeric(batch[rownames(pc),"Library_Complexity"])
resLm<-lm(pc[,1]~library)
anova(resLm)$Pr[1]

summary(resLm)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4283.5 -1005.4  -394.8   335.0  6024.0 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -4925.73     373.25  -13.20   <2e-16 ***
#   library       264.96      18.42   14.38   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1618 on 117 degrees of freedom
# Multiple R-squared:  0.6387,	Adjusted R-squared:  0.6356 
# F-statistic: 206.8 on 1 and 117 DF,  p-value: < 2.2e-16
# #ccl : locis avec HiSD sont mauvais, on doit les enlever

#bonne stat !!
#! faire deterQual2() avec ça !
#!mean : bizarre car Meanint plus lowConfscore mais pas pc1~library => ptete PC2 ?



pheatmap(logpvals.raw[,1:39],cluster_rows = F,cluster_cols = F,labels_col= paste("PC",1:39),display_numbers = T)




