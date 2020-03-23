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

#nas :
mat<-as.matrix(data_all[,samples])
sum(is.na(mat))/length(mat) #1%
sum(rowSums(is.na(mat))>0)/nrow(mat) #7.9% des locis
# data_all$pctNA<-rowSums(is.na(mat))/ncol(mat)
# head(data_all)
# fwrite(data_all,"../../ref/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt")

batch=read.csv2("../../ref/20-02-17_batch_119samples_IUGRLGACTRL_noNA_withHpaC.csv",header=T,row.names = 1)


#fct use : 
deterSeuilQC<-function(dataCpG,metrique,qTestes,qualMetriques=c(1,2),test="quantile",lowerThan=T){
  if (1%in% qualMetriques){
    quals<-rep(0,length(qTestes))
  }
  if (2%in% qualMetriques){
    quals2<-rep(0,length(qTestes))
    quals3<-rep(0,length(qTestes))
    quals4<-rep(0,length(qTestes))
  }
  
  print(paste("determine si",metrique,"peut permettre d'exclure des locis de faible confiance"))
  for (i in 1:length(qTestes)){
    q<-qTestes[i]
    print(paste("seuil",q))
    if(test=="quantile"){
      q<-quantile(dataCpG[,metrique],q,na.rm=T)
      
    }
    if(lowerThan){
      locisLo<-na.exclude(rownames(dataCpG)[dataCpG[,metrique]<=q])
      print(paste("analyses qualité des",length(locisLo),"locis avec",metrique,"<",q))
    }else{
      locisLo<-na.exclude(rownames(dataCpG)[dataCpG[,metrique]>=q])
      print(paste("analyses qualité des",length(locisLo),"locis avec",metrique,">",q))
    }
    
    
    
    if(length(locisLo)>100000){
      locisLo<-sample(locisLo,100000)
    }
    dataCpGLo<-as.matrix(dataCpG[locisLo,samples])
    
    if (1%in% qualMetriques){
      quals[i]<-deterQual(dataCpGLo,maxMethyl = 5)
    }
    
    
    if (2%in% qualMetriques){
      res<-deterQual2(dataCpGLo,batch)
      quals2[i]<-res$p
      quals3[i]<-res$pctPC
      quals4[i]<-res$r2
    }
    
  }
  
  print(plot(qTestes,quals,main=paste('pctLocis Avec Vrais Zeros en fonction',metrique)))
  print(plot(qTestes,quals2,main=paste('pval PC~library',metrique),log="y"))
  print(plot(qTestes,quals3,main=paste('pctPC Library',metrique)))
  print(plot(qTestes,quals4,main=paste('r2 PC~library',metrique)))
  return(list(qual1=quals,qual2.p=quals2,quals2.pctPC=quals3,qual2.r=quals4 ))
}

#je veux enlever les locis de mauvaise qual selon un critére bien définit
#ccl des 2 premiers tests : methyl entre 0 et 20, 
#bon locis = pctdescore entre 0 et 20 élévé ou pct0 faible (<0.7) (ou sd > 25 (locis ++ unmethylé), mais pas utilisé pour le moment)
#qualitéData=pct0*pctEntre0et20

#seuil optimale = maximiser ses 3 metriques avec msp1c, confscore, nbVraisZeros, pct0 

#comment : chercher seuil bas (d'exclusion) et haut (d'inclusion)
#1)seuil bas : pour quantile de 0.1 a 0.9, tant que locis < quantile x n'augmente pas pctBonlocis augmenter

# trouver seuil  a l'aide des 2 metriques de quals
#a) fct deterQual
source('scripts/deterDataQual.R')

#b) qualité data_all
mat<-as.matrix(data_all[,samples])
deterQual(mat) #0.4153642
deterQual2(mat,batch) #PC 1  ( 16.9 % de la variance a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3801368616723

#c) deter seuil en fct quantile 
#critere d'inclusion locis : on est sur que c'est un bon locis si :
# - pct0 faible : a qpartir de quel pct0 on  des bon/mauvais locis ?
pctTestes<-1:19/20
qualsPct0<-deterSeuilQC(data_all,"pct0",pctTestes,test = "brut",lowerThan = F) 
#locis avec 35% de 0 sont influencés par librarycomplexity => crtère d'inclusion des locis dans possiblement a enlever
#on veut enlelver les locis faux zeros parmi les locis pct0>0.35
locis0<-rownames(data_all)[data_all$pct0>0.35]
data_0<-data_all[locis0,]
dim(data_0) #1446462     131


#mm chose pour le pctNA
max(data_all$pctNA)
pctTestes<-0:7/50
qualsNAs<-deterSeuilQC(data_all,"pctNA",pctTestes, test= "brut",lowerThan = F) #exclusion : pctNA>0.02?
pctTestes<-0:5/200
qualsNAs<-deterSeuilQC(data_all,"pctNA",pctTestes, test= "brut") #exclusion : pctNA!=0


locis0<-rownames(data_all)[data_all$pct0>0.35]
data_0<-data_all[locis0,]
dim(data_0) #1446462     131

#tester, msp1c, conf

# =  confidencescorehi
qTestes<-1:19/20
qualsConfScore<-deterSeuilQC(data_all,"confidenceScore",qTestes,lowerThan = F) 
#inclusion : > 0.65

# =  confidencescoreNormhi
qTestes<-1:19/20
qualsConfScoreNorm<-deterSeuilQC(data_all,"confidenceScoreNorm",qTestes,lowerThan = F) #pas plus inclusif que confidenceScore



# mean 
qTestes<-1:19/20
qualsMeans<-deterSeuilQC(data_all,"mean",qTestes,lowerThan = F) 
#inclusion : >0.75



# sd
qTestes<-1:19/20
qualsMeans<-deterSeuilQC(data_all,"sd",qTestes,lowerThan = F) 
#inclusoin >0.75




locisArmNA<-rownames(data_all)[data_all$pctNA!=0]
length(locisArmNA) #135k


locisRmvable<-na.omit(rownames(data_all)[(data_all$pct0>0.35&
                                   data_all$confidenceScore<quantile(data_all$confidenceScore,0.65,na.rm=T)&
                                   
                                   data_all$mean<quantile(data_all$mean,0.75,na.rm=T)&
                                   data_all$sd<quantile(data_all$sd,0.75,na.rm=T))])


length(locisRmvable) # 992749      
head(locisRmvable)

data_rmvable<-data_all[locisRmvable,]
#enrichissement en conf et mean et sd
hist(data_rmvable$RankConfidenceScore,breaks=100)
hist(data_all[!(rownames(data_all)%in%locisRmvable),"RankConfidenceScore"],breaks=100)

locisKeep<-rownames(data_all)[!(rownames(data_all)%in%union(locisRmvable,locisArmNA))]
length(locisKeep)
#msp1c
qTestes<-1:9/10
qualsMsp1c<-deterSeuilQC(data_rmvable,"msp1c",qTestes,qualMetriques = 1)  
##entre 0.04 et 0.2
qTestes<-2:10/50
qualsMsp1c2<-deterSeuilQC(data_rmvable,"msp1c",qTestes,qualMetriques = 1) 
#exclusion 0.10


#confScore
qTestes<-1:9/10
qualsConf<-deterSeuilQC(data_rmvable,"confidenceScore",qTestes,qualMetriques = 1) #inclusion >0.7
#mais entre 0.1 et 0.2 ya un enrichissemnet ++ en vrai zeros
qTestes<-1:10/50
qualsConf<-deterSeuilQC(data_rmvable,"confidenceScore",qTestes,qualMetriques = 1) #exclusion : 0.08

#mean
qTestes<-1:9/10
qualsMean<-deterSeuilQC(data_rmvable,"mean",qTestes,qualMetriques = 1)  #inclusion : <0.1?

#0.1 inclusion  ? 
qTestes<-1:15/50
qualsMean<-deterSeuilQC(data_rmvable,"mean",qTestes,qualMetriques = 1) 
#0.04 pour etre exact



#sd
qTestes<-1:9/10
qualsConf<-deterSeuilQC(data_rmvable,"sd",qTestes,qualMetriques = 1) #inclusion <0.2 ?

qTestes<-1:10/50
qualsConf<-deterSeuilQC(data_rmvable,"sd",qTestes,qualMetriques = 1) #inclusion <0.04


#nbMethylNonzeros
qTestes<-1:9/10
qualsConf<-deterSeuilQC(data_rmvable,"nbMethylNonZeros",qTestes,qualMetriques = 1) #exclusion nbMethyl = 0



#ccl : locis exclu = nbMethylNonZeros = 0, msp1c<q0.14, et confscore <q0.12 (mais inutile car msp1c suffit), 
locisArm<-na.omit(locisRmvable[data_rmvable$nbMethylNonZeros==0|data_rmvable$msp1c<quantile(data_rmvable$msp1c,0.10)|data_rmvable$pctNA>0])
sum(is.na(locisArm))
head(locisArm)
length(locisArm) #361279



#locislomean bad qual locis ?
plot(density(data_all$mean))
lines(density(na.omit(data_all[locisRmvable,"mean"])),col=2)
lines(density(na.omit(data_all[locisArm,"mean"])),col=3) #plutot leauvais niveau mso1c
plot(density(data_all$sd))
lines(density(na.omit(data_all[locisRmvable,"sd"])),col=2)
lines(density(na.omit(data_all[locisArm,"sd"])),col=3)

hist(data_rmvable$RankConfidenceScore,breaks=100)
hist(data_rmvable[locisArm,"RankConfidenceScore"],breaks=100)


#check quals locisarm et loocis filtre
locis<-rownames(data_all)[!(rownames(data_all)%in%locisArm)]

deterQual(as.matrix(data_all[locisArm,samples]),maxMethyl = 5)#0.177
deterQual(as.matrix(data_all[locis,samples]),maxMethyl = 5) #0.45 au lieu de 0.415

deterQual2(as.matrix(data_all[locisArm,samples]),batch) # PC 1  ( 8.1 % de la variance a R2 avec Library_Complexity = 0.33 et pval = 10^ -11.4460774167366
deterQual2(as.matrix(data_all[locis,samples]),batch)
#PC 1  ( 21.1 % de la variance a R2 avec Library_Complexity = 0.74 et pval = 10^ -35.0833868846837
#au lieu de : PC 1  ( 16.9 % de la variance a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3801368616723
data_F<-data_all[locis,]

#les locis rmvable mais pas exclus
locisRmvable<-locisRmvable[!(locisRmvable%in%locisArm)]
length(locisRmvable) #818540

#enleve si ne passe pas critere inclusion de sd et mean < q0.04:
locisArm<-locisRmvable[data_all[locisRmvable,"sd"]>quantile(data_rmvable$sd,0.04, na.rm=T)&
                         data_all[locisRmvable,"mean"]>quantile(data_rmvable$mean,0.04, na.rm=T)&
                         data_all[locisRmvable,"msp1c"]<quantile(data_rmvable$msp1c,0.7, na.rm=T)&
                         data_all[locisRmvable,"confidenceScore"]<quantile(data_rmvable$confidenceScore,0.7, na.rm=T)&
                         data_all[locisRmvable,"nbMethylNonZeros"]<4]
length(locisArm) #477197
locis<-locisRmvable[!(locisRmvable%in%locisArm)]
length(locis)#192k
hist(data_rmvable$RankConfidenceScore,breaks=100)
hist(data_rmvable[locisArm,"RankConfidenceScore"],breaks=100) 


deterQual(as.matrix(data_all[locis,samples]),maxMethyl = 5)#0.37
deterQual(as.matrix(data_all[locisArm,samples]),maxMethyl = 5)#0.17


deterQual2(as.matrix(data_all[locis,samples]),batch) #PC 1  ( 6.4 % de la variance a R2 avec Library_Complexity = 0.35 et pval = 10^ -12.0958497599772
deterQual2(as.matrix(data_all[locisArm,samples]),batch) #"PC 1  ( 8.1 % de la variance a R2 avec Library_Complexity = 0.33 et pval = 10^ -11.446077416736

#ccl finale : 
locis<-setdiff(rownames(data_F),locisArm)
length(locis)#120k
head(locis)
sum(is.na(locis))
write.table(locis,"locisPassantQC.txt",sep = ",")
data_F<-data_all[locis,]

head(data_F)
dim(data_F)
summary(data_F)
sum(is.na(data_F[,samples]))

sum(rowSums(is.na(data_F[,samples])>0))/nrow(data_F) #44% ont des NAs.. wtf

