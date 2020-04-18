#config et library
options(stringsAsFactors=F)
set.seed(12345)
library(data.table)
library(stringr)
library(limma)
library(plyr)
source("scripts/utils.R")
#output dir 
script_name <- "estim_pval_cutoff_with_permut"
outputDir <- file.path("analyses","withoutIUGR",script_name)
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

models<-list()
model<-14
names(batch)
varToModel<-c("Group_Sex",'batch',"Mat.Age","latino","Group_Complexity_Fac")

#samples = CTRL ou LGA
samples_F_F<-samples[rowSums(is.na(batch[samples,varToModel]))==0&batch[samples,"Group_name"]%in%c("C","L")] 
length(samples_F_F) # 70
table(batch[samples_F_F,"Group_name"])
# C  L 
# 34 36 

mat.age<-as.numeric(batch[samples_F_F,"Mat.Age"])

batches<-factor(batch[samples_F_F,"batch"])
group_sex<-factor(batch[samples_F_F,"Group_Sex"])
group_sex<-revalue(group_sex,c("1"="FC","2"="MC","3"="FI","4"="MI","5"="FL","6"="ML"))
group_complexity_fac<-factor(batch[samples_F_F,"Group_Complexity_Fac"])
latino<-factor(batch[samples_F_F,"latino"])
#complexity<-as.numeric(batch[samples_F_F,"Library_Complexity"])

#model 14 :
formule<- ~0 + group_sex  + batches  + latino + mat.age + group_complexity_fac
#on veut enregistrer dans un df mes pvalues des 800k locis pour chaque permut 
resP<-data.frame(row.names = locisF)
# et pour notre distrib observé, en 1er :

design<-model.matrix(formule)
colnames(design)<-make.names(colnames(design))
fit <- lmFit(data_F[,samples_F_F], design)  
cont.matrix <- makeContrasts(C.L = "(group_sexFC+group_sexMC)-(group_sexFL+group_sexML)",
                             MC.ML="group_sexMC-group_sexML",
                             FC.FL="group_sexFC-group_sexFL",
                             ML.FL="group_sexML-group_sexFL",
                             MC.FC="group_sexMC-group_sexFC",
                             F.M="(group_sexFC+group_sexFL)-(group_sexMC+group_sexML)",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
res<-topTable(fit2,coef="FC.FL",n=Inf)
nrow(res)
head(res)
resP[rownames(res),"pval.obs.FC.FL"]<-res$P.Value

res<-topTable(fit2,coef="MC.ML",n=Inf)
resP[rownames(res),"pval.obs.MC.ML"]<-res$P.Value


res<-topTable(fit2,coef="ML.FL",n=Inf)
resP[rownames(res),"pval.obs.ML.FL"]<-res$P.Value
head(resP)
#for 100 permut : on melange groupes_sex a lint de chaque sex, et on fait les memes compas de groupe
subBatch<-batch[batch$Group_name%in%c("C","L"),]
samples_F_M<-samples_F_F[batch[samples_F_F,"Gender"]=="M"]
samples_F_Fe<-samples_F_F[batch[samples_F_F,"Gender"]=="F"]
batch_M<-batch[samples_F_M,]
batch_F<-batch[samples_F_Fe,]
dfPerm<-data.frame(row.names = samples_F_F,obs=batch[samples_F_F,"Group_Sex"])
for(i in 1:100){
  print(paste("perm",i))
  subBatch[samples_F_M,"Group_Sex"]<-sample(batch_M$Group_Sex)
  subBatch[samples_F_Fe,"Group_Sex"]<-sample(batch_F$Group_Sex)
  
  dfPerm[,paste0("perm.",i)]<-subBatch[samples_F_F,"Group_Sex"]
  print(dfPerm[,c(1,i+1)])
  
  group_sex<-factor(subBatch[samples_F_F,"Group_Sex"])
  group_sex<-revalue(group_sex,c("1"="FC","2"="MC","5"="FL","6"="ML"))
  formule<- ~0 + group_sex  + batches  + latino + mat.age + group_complexity_fac
  design<-model.matrix(formule)
  colnames(design)<-make.names(colnames(design))
  fit <- lmFit(data_F[,samples_F_F], design)  
  cont.matrix <- makeContrasts(C.L = "(group_sexFC+group_sexMC)-(group_sexFL+group_sexML)",
                               MC.ML="group_sexMC-group_sexML",
                               FC.FL="group_sexFC-group_sexFL",
                               levels=design)
  
  
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  
  
  for (compa in c("FC.FL","MC.ML")){
    res<-topTable(fit2,coef=compa,n=Inf)
    print(head(res))
    resP[rownames(res),paste0("pval.",i,compa)]<-res$P.Value
    
  }
  
}
write.csv2(resP,paste0(output,"_100_permuts_within_sex.csv"))

resP<-resP[,!(str_detect(names(resP),"ML.FL"))]
head(resP)

#permut between sex

#obs : 
resP1<-data.frame(row.names = locisF)
group_sex<-factor(batch[samples_F_F,"Group_Sex"])
group_sex<-revalue(group_sex,c("1"="FC","2"="MC","5"="FL","6"="ML"))
formule<- ~0 + group_sex  + batches  + latino + mat.age + group_complexity_fac

# et pour notre distrib observé, en 1er :

design<-model.matrix(formule)
colnames(design)<-make.names(colnames(design))
fit <- lmFit(data_F[,samples_F_F], design)  
cont.matrix <- makeContrasts(C.L = "(group_sexFC+group_sexMC)-(group_sexFL+group_sexML)",
                             MC.ML="group_sexMC-group_sexML",
                             FC.FL="group_sexFC-group_sexFL",
                             ML.FL="group_sexML-group_sexFL",
                             MC.FC="group_sexMC-group_sexFC",
                             F.M="(group_sexFC+group_sexFL)-(group_sexMC+group_sexML)",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
res<-topTable(fit2,coef="MC.FC",n=Inf)
nrow(res)
head(res)
resP1[rownames(res),"pval.obs.MC.FC"]<-res$P.Value

res<-topTable(fit2,coef="ML.FL",n=Inf)
resP1[rownames(res),"pval.obs.ML.FL"]<-res$P.Value
head(resP1)

subBatch<-batch[batch$Group_name%in%c("C","L"),]
samples_F_C<-samples_F_F[batch[samples_F_F,"Group_name"]=="C"]
samples_F_L<-samples_F_F[batch[samples_F_F,"Group_name"]=="L"]
batch_C<-batch[samples_F_C,]
batch_L<-batch[samples_F_L,]
dfPerm1<-data.frame(row.names = samples_F_F,obs=batch[samples_F_F,"Group_Sex"])

for(i in 1:100){
  print(paste("perm",i))
  subBatch[samples_F_C,"Group_Sex"]<-sample(batch_C$Group_Sex)
  subBatch[samples_F_L,"Group_Sex"]<-sample(batch_L$Group_Sex)
  
  dfPerm1[,paste0("perm.",i)]<-subBatch[samples_F_F,"Group_Sex"]
  print(dfPerm1[,c(1,i+1)])
  
  group_sex<-factor(subBatch[samples_F_F,"Group_Sex"])
  group_sex<-revalue(group_sex,c("1"="FC","2"="MC","5"="FL","6"="ML"))
  formule<- ~0 + group_sex  + batches  + latino + mat.age + group_complexity_fac
  design<-model.matrix(formule)
  colnames(design)<-make.names(colnames(design))
  fit <- lmFit(data_F[,samples_F_F], design)  
  cont.matrix <- makeContrasts(C.L = "(group_sexFC+group_sexMC)-(group_sexFL+group_sexML)",
                               ML.FL="group_sexML-group_sexFL",
                               MC.FC="group_sexMC-group_sexFC",
                               levels=design)
  
  
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  
  
  for (compa in c("MC.FC","ML.FL")){
    res<-topTable(fit2,coef=compa,n=Inf)
    print(head(res))
    resP1[rownames(res),paste0("pval.",i,compa)]<-res$P.Value
    
  }
  
}
write.csv2(resP1,paste0(output,"_100_permuts_between_sex.csv"))

head(resP1[,1:10])

#distrib pval min foreach permut
resP.FC.FL<-resP[,str_detect(names(resP),"FC.FL")]
head(resP.FC.FL)
plot(density(apply(log10(as.matrix(resP.FC.FL[,-1])),2,min)))
abline(v=min(log10(resP.FC.FL[,1])))

resP.MC.ML<-resP[,str_detect(names(resP),"MC.ML")]
head(resP.MC.ML)
plot(density(apply(log10(as.matrix(resP.MC.ML[,-1])),2,min)))
abline(v=min(log10(resP.MC.ML[,1])))

resP.ML.FL<-resP[,str_detect(names(resP),"ML.FL")]
head(resP.ML.FL)
plot(density(apply(log10(as.matrix(resP.ML.FL[,-1])),2,min)))
abline(v=min(log10(resP.ML.FL[,1])))

#distrib pval min foreach locis
resP.FC.FL<-resP[,str_detect(names(resP),"FC.FL")]
head(resP.FC.FL)
plot(density(apply(log10(as.matrix(resP.FC.FL[,-1])),1,min)))
abline(v=log10(0.001))
lines(density(log10(as.matrix(resP.FC.FL[,1]))),col=2)

resP.MC.ML<-resP[,str_detect(names(resP),"MC.ML")]
head(resP.MC.ML)
plot(density(apply(log10(as.matrix(resP.MC.ML[,-1])),1,min)))
abline(v=log10(0.001))

resP.ML.FL<-resP[,str_detect(names(resP),"ML.FL")]
head(resP.ML.FL)
plot(density(apply(log10(as.matrix(resP.ML.FL[,-1])),1,min)))
abline(v=log10(0.001))


resP.FC.FL<-resP[,str_detect(names(resP),"FC.FL")]
head(resP.FC.FL)
plot(density(apply(as.matrix(resP.FC.FL[,-1]),2,mean)))
abline(v=mean(resP.FC.FL[,1]))

resP.MC.ML<-resP[,str_detect(names(resP),"MC.ML")]
head(resP.MC.ML)
plot(density(apply(as.matrix(resP.MC.ML[,-1]),2,mean)))
abline(v=mean(resP.MC.ML[,1]))

resP.ML.FL<-resP[,str_detect(names(resP),"ML.FL")]
head(resP.ML.FL)
plot(density(apply(as.matrix(resP.ML.FL[,-1]),2,mean)))
abline(v=mean(resP.ML.FL[,1]))
