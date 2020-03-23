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
#tester, msp1c, conf

# =  confidencescorehi
qTestes<-1:19/20
qualsConfScore<-deterSeuilQC(data_all,"confidenceScore",qTestes,lowerThan = F) 
#inclusion : > 0.65

# =  confidencescoreNormhi
qTestes<-1:19/20
qualsConfScore<-deterSeuilQC(data_all,"confidenceScoreNorm",qTestes,lowerThan = F) #pas plus inclusif que confidenceScore



# nb methyl Non zeros 
qTestes<-1:19/20
qualsConfScore<-deterSeuilQC(data_all,"nbMethylNonZeros",qTestes,lowerThan = F) 
#inclusion : >0.65


# mean 
qTestes<-1:19/20
qualsMeans<-deterSeuilQC(data_all,"mean",qTestes,lowerThan = F) 
#inclusion : >0.75



# sd
qTestes<-1:19/20
qualsMeans<-deterSeuilQC(data_all,"sd",qTestes,lowerThan = F) 
#inclusoin >0.75


locisRmvable<-rownames(data_all)[data_all$pct0>0.35&
                                   data_all$confidenceScore<quantile(data_all$confidenceScore,0.65,na.rm=T)&
                                   data_all$nbMethylNonZeros<4&
                                   data_all$mean<quantile(data_all$mean,0.75,na.rm=T)&
                                   data_all$mean<quantile(data_all$sd,0.75,na.rm=T)]

length(locisRmvable) #701095        


data_rmvable<-data_all[locisRmvable,]
#enrichissement en conf et mean
plot(density(data_all$msp1c))
lines(density(data_rmvable$msp1c),col=2)
#msp1c
qTestes<-1:9/10
qualsMsp1c<-deterSeuilQC(data_rmvable,"msp1c",qTestes) 
##entre 0.04 et 0.2
qTestes<-2:10/50
qualsMsp1c2<-deterSeuilQC(data_rmvable,"msp1c",qTestes) 
#<0.12 => exclusion

#confScore
qTestes<-1:9/10
qualsConf<-deterSeuilQC(data_rmvable,"confidenceScore",qTestes,qualMetriques = 1) #confidence score nul car pval ~library de  plus en plus fort
#mais entre 0.1 et 0.2 ya un enrichissemnet ++ en vrai zeros
qTestes<-1:10/50
qualsConf<-deterSeuilQC(data_rmvable,"confidenceScore",qTestes,qualMetriques = 1) #seuil : 0.1

#mean
qTestes<-1:9/10
qualsMean<-deterSeuilQC(data_rmvable,"mean",qTestes,qualMetriques = 1)
#0.2 exclusion ? 
qTestes<-1:15/50
qualsMean<-deterSeuilQC(data_0,"mean",qTestes) 
#0.12 pour etre exact


# et <0.8 exclusion ?


#sd
qualsMean<-deterSeuilQC(data_0,"sd",qTestes)

#quals1 diminue jusqu'a 0.2 puis augmente ! 
#quals2 pval diminue jusya 0.1 (gros pic de 0.06 à 0.1)puis augmente, r2 le + fort entre 0.08 et 0.1, pval 
#donc peut etre enlver les locis entre 0.06 et 0.1 ?

plot(density(data_all$mean))
abline(v=0.06)
abline(v=0.1,col=2)
locisARm<-locis0[data_all[locis0,"mean"]>quantile(data_all$mean,0.06)&data_all[locis0,"mean"]<quantile(data_all$mean,0.1)]
length(locisARm)
#locislomean bad qual locis ?
plot(density(log10(data_all$msp1c)))
lines(density(na.omit(log10(data_all[locis0,"msp1c"]))),col=2)
lines(density(na.omit(log10(data_all[locisARm,"msp1c"]))),col=3) #plutot leauvais niveau mso1c

plot(density(na.omit(log10(data_all$confidenceScore))))
lines(density(na.omit(log10(data_all[locis0,"confidenceScore"]))),col=2)
lines(density(na.omit(log10(data_all[locisARm,"confidenceScore"]))),col=3)#bad locis

#check quals locisarm et loocis filtre
locis<-rownames(data_all)[!(rownames(data_all)%in%locisARm)]
deterQual(as.matrix(data_all[locisARm,samples]),maxMethyl = 5)#0.332
deterQual(as.matrix(data_all[locis,samples]),maxMethyl = 5) #0.42 au lieu de 0.415
deterQual2(as.matrix(data_all[locisARm,samples]),batch) # pas signif, donc la pval dépend des moins bon


deterQual2(as.matrix(data_all[locis,samples]),batch)
#PC 1  16.9 % de la variance a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3754771648479

#moins stringent : entre 0.02 et 0.1, on gagne en qual ? quantile(data_all$mean,0.02) = 0.09
locisARm<-locis0[data_all[locis0,"mean"]>quantile(data_all$mean,0.02)&data_all[locis0,"mean"]<quantile(data_all$mean,0.1)]

#check quals locisarm et loocis filtre
locis<-rownames(data_all)[!(rownames(data_all)%in%locisARm)]
deterQual(as.matrix(data_all[locisARm,samples]),maxMethyl = 5) #0.329
deterQual(as.matrix(data_all[locis,samples]),maxMethyl = 5) #0.425 au lieu de 0.42 au lieu de 0.415


deterQual2(as.matrix(data_all[locis,samples]),batch)#"PC 1  ( 16.9 % de la variance a R2 avec Library_Complexity = 0.76 et pval = 10^ -37.3736045681155"

#pas de diminution importante, de la pval/r2, on utilise pas mean comme filtre

#sd
qTestes<-1:19/20
qualsSD<-deterSeuilQC(data_all,"sd",qTestes)
#quals diminue jusqu'a 0.4, plateau jusqu'a 0.7 , et redescend ensuite 
#quals2 diminue formatement a 0.7 assi => enlver >0.7 ?
locisARm<-locis0[data_all[locis0,"sd"]>quantile(data_all$mean,0.7)]

#check quals locisarm et loocis filtre
locis<-rownames(data_all)[!(rownames(data_all)%in%locisARm)]
deterQual(as.matrix(data_all[locisARm,samples]),maxMethyl = 5) #0.40
deterQual(as.matrix(data_all[locis,samples]),maxMethyl = 5) #0.43 au lieu de 0.415 donc o gagne

deterQual2(as.matrix(data_all[locis,samples]),batch) #PC 1 21.1 % de la variance a R2 avec Library_Complexity = 0.74 et pval = 10^ -35.0833868846837

#donc on diminue en pval par rapport a avant c'est mieuc

plot(density(log10(data_all$msp1c)))

lines(density(na.omit(log10(data_all[locis0,"msp1c"]))),col=2)
lines(density(na.omit(log10(data_all[locisARm,"msp1c"]))),col=3)
lines(density(na.omit(log10(data_all[locis,"msp1c"]))),col=4)


plot(density(na.omit(log10(data_all$confidenceScore))))
lines(density(na.omit(log10(data_all[locis0,"confidenceScore"]))),col=2)
lines(density(na.omit(log10(data_all[locisARm,"confidenceScore"]))),col=3)
lines(density(na.omit(log10(data_all[locis,"confidenceScore"]))),col=4)  #bcp mieux !

#ccl sd : on filtre au dessus de 0.7



