
#get eQTL info : links SNP with RNA expr

#config et library
options(stringsAsFactors=F)
set.seed(12345)
#library(limma)
library(data.table)
library(stringr)
#source("scripts/utils.R")

#output dir 
script_name <- "eQTL_test"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
date<-Sys.Date()
output <- file.path(outputDir,date)

#I found 3 interesting eQTL db : 
#1) GTeX, https://gtexportal.org/home/, 
#2) eQTL in lymphoblastoid cell lines, https://www.hsph.harvard.edu/liming-liang/software/eqtl/
#3) dutch biobanks with Conditional eQTL https://eqtl.onderzoek.io/index.php?page=info

#1) GTeX : 
#get data in : https://gtexportal.org/home/datasets 
# info des fichiers sur : https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/README_eQTL_v8.txt
#les fichiers qui m'interessent : 
  #File extensions: *.signif_variant_gene_pairs.txt.gz, because contain Significant variant-gene pairs
wb_eQTL<-read.table("../../ref/eQTL/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt",
                                  header = T,
                                  sep = "\t")
head(wb_eQTL)
#stat :
nrow(wb_eQTL) #2 414 653 variant-gene pairs
length(unique(wb_eQTL$gene_id)) #dans 12360 gene differents 
plot(density(table(wb_eQTL$gene_id)))
summary(as.vector(table(wb_eQTL$gene_id)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    22.0    87.0   195.4   230.0  8571.0 
#en moy : 195 var par gene

#need convert ensg. in symbols
#with bitr
library(clusterProfiler)
library(org.Hs.eg.db)
#need remove the ".x" in order to be understood by the translator
library(stringr)
wb_eQTL$gene_idTrunked<-str_remove(wb_eQTL$gene_id,"\\..*")

genes.df<-bitr(unique(wb_eQTL$gene_idTrunked),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
#18.62% of input gene IDs are fail to map
head(genes.df)
head(unique(wb_eQTL$gene_idTrunked))
#add symbol
for(symbol in genes.df$SYMBOL){
  ensembl<-genes.df$ENSEMBL[genes.df$SYMBOL==symbol]
  wb_eQTL[which(wb_eQTL$gene_idTrunked%in%ensembl),"gene"]<-symbol
}
head(wb_eQTL,20)
write.csv2(wb_eQTL,"../../ref/Whole_Blood.eQTL.csv")

#integr in annot : make region eQTL : +/-x pb du locus signif
#1) deter x : regarder distribution distances entre locis consec
wb_eQTL$distDuSuivant<-as.vector(sapply(1:(nrow(wb_eQTL)),function(i){
  if(wb_eQTL[i+1,"gene_id"]==wb_eQTL[i,"gene_id"]&i!=nrow(wb_eQTL)){
    return(wb_eQTL[i+1,"tss_distance"]-wb_eQTL[i,"tss_distance"])
  }else{
    return(NA)
  }
  
}))

plot(density(na.omit(wb_eQTL$distDuSuivant[wb_eQTL$distDuSuivant<2000]))) #pique à 30 pb
abline(v=)
#zoom
plot(density(na.omit(wb_eQTL$distDuSuivant[wb_eQTL$distDuSuivant<600]))) #pas d'autre pique signif

summary(wb_eQTL$distDuSuivant)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       0     113     394    2679    1153 1852417   12360 

#grace à ça on estime que, les locis ld'une meme regions regulatrice sont à 30 pb les uns des autres, 
#Avec un incertitude de 50% on definit notre fenetre a +/-45pb pour définir une région regulant le gene
# donc 90pb autour du locis signif et une region régulant le gène
wb_eQTL$region<-paste(wb_eQTL$tss_distance-45,wb_eQTL$tss_distance+45,sep = ":")
head(wb_eQTL)
write.csv2(wb_eQTL,"../../ref/Whole_Blood.eQTL.csv")



##integration features 4/6 et eQTL 
library(data.table)
library(stringr)
wb_eQTL<-fread("../../ref/Whole_Blood.eQTL.csv",select = 2:17)
wb_eQTL
#1) need chr pos absolu
#get translator hg38>hg19
translator<-fread("../../ref/eQTL/GTEx_Analysis_v8_eQTL/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt",
                  select = c(1,8))
translator
wb_eQTL<-merge.data.table(wb_eQTL,translator,by="variant_id")

wb_eQTL[,chr:=str_extract(variant_id_b37,"^[0-9XY]{1,2}") ]
wb_eQTL<-wb_eQTL[!is.na(chr),]
wb_eQTL[,chr:=paste0("chr",chr) ]



wb_eQTL[,pos:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2)) ]
wb_eQTL

fwrite(wb_eQTL,"../../ref/Whole_Blood.eQTL.csv",sep = ";")

##In which featureRegion are this SNP ?

features<-fread("../Reference/feature_all_042820.bed",skip = 1,select = c(1:4),col.names = c("chr","start","end","type"))

features[,longueur:=end-start]
features[,type:=as.factor(type)]
library(stringr)
features[,chrReg:=paste(str_extract(chr,"[0-9XY]{1,2}$"),start,end,sep = ":")]
fwrite(features,"../../ref/CD34_all_chromatin_feature.csv",sep = ";")
#length region par type
library(ggplot2)
p<-ggplot(data = features,aes(x=type))
p+geom_boxplot(aes(y=longueur))+scale_y_continuous(trans = 'log10')  
p+geom_bar()
# les regions "0" majoriatire (~300k reg contre ~150k region pour les autre), correspondant a l'heterochromatine,
#sont de large région, (75% entre 3k et 5k), les autres sont à 2k de longueur

#loading
wb_eQTL<-fread("../../ref/Whole_Blood.eQTL.csv")
wb_eQTL
features<-fread("../../ref/CD34_all_chromatin_feature.csv")
features

#pour 1 region type,tout les eQtl dans celle ci sont annoter
reg<-features[3,]
wb_eQTL[chr==reg$chr&pos>reg$start&pos<reg$end,"type"]<-reg$type
i<-1
while(nrow(wb_eQTL[!is.na(type),])==0){
  print(i)
  reg<-features[i,]
  wb_eQTL[chr==reg$chr&pos>reg$start&pos<reg$end,"type"]<-reg$type
  i<-i+1
}
wb_eQTL[!is.na(type),] #ça marche, maintenant il faut subsetter sinon ce sera trop long

#test annot eQTL avec features
for (chrom in "chr21"){
  print(chrom)
  n<-ceiling(table(wb_eQTL$chr)[chrom]/5000)
  print(paste("partionner",n,"fois"))
  q2<-0
  for(i in 1:3){
    f<-round(i/n,2)
    q1<-quantile(wb_eQTL[chr==chrom,]$pos,f)
    sub_eQTL<-wb_eQTL[chr==chrom&pos<=q1&pos>q2,]
    
    if(nrow(sub_eQTL)>0){
      print(paste(chrom,"part",i,", annotation de",nrow(sub_eQTL),"locis.."))
      sub_features<-features[chr==chrom&end<=q1&start>q2,]
      for(j in 1:nrow(sub_features)){
        reg<-sub_features[j,]
        sub_eQTL[pos>=reg$start&pos<=reg$end,"type"]<-reg$type
        sub_eQTL[pos>=reg$start&pos<=reg$end,"chrRegion"]<-reg$chrReg
      }
      q2<-q1
      #fwrite(sub_eQTL,file = paste0(output,"annotFeat",chrom,"part",i,".txt"))
      print("..annotes et enregistrer")
    }
    
  }
}

sub_eQTL

for (chrom in unique(features$chr)){
  print(chrom)
  n<-ceiling(table(wb_eQTL$chr)[chrom]/5000)
  print(paste("partionner",n,"fois"))
  q2<-0
  for(i in 1:n){
    f<-round(i/n,2)
    q1<-quantile(wb_eQTL[chr==chrom,]$pos,f)
    sub_eQTL<-wb_eQTL[chr==chrom&pos<=q1&pos>q2,]
    
    if(nrow(sub_eQTL)>0){
      print(paste(chrom,"part",i,", annotation de",nrow(sub_eQTL),"locis.."))
      sub_features<-features[chr==chrom&end<=q1&start>q2,]
      for(j in 1:nrow(sub_features)){
        reg<-sub_features[j,]
        sub_eQTL[pos>=reg$start&pos<=reg$end,"type"]<-reg$type
        sub_eQTL[pos>=reg$start&pos<=reg$end,"chrRegion"]<-reg$chrReg
      }
      q2<-q1
      fwrite(sub_eQTL,file = paste0(output,"annotFeat",chrom,"part",i,".txt"))
      print("..annotes et enregistrer")
    }
    
  }
}
#check si bien annoté :
df<-fread("analyses/eQTL_test/2020-04-29annotFeatchr1part37.txt")
df[is.na(type),] #oui, pas de na

#test merge with wb_eqtl
#need put col type in the same type (factor)
df[,type:=as.factor(type)] 
test<-rbind(df,sub_eQTL)

#recup et merge res
sub_eQTL[,type:=as.character(type)]

wb_eQTL2<-sub_eQTL

for(chr in paste0("chr",c(1:22,"X"))){
  print(chr)
  for(file in list.files("analyses/eQTL_test/",pattern = paste0(chr,"part"))){
    print(strsplit(file,chr)[[1]][2])
    sub_eQTL<-fread(paste0("analyses/eQTL_test/",file))
    wb_eQTL2<-rbind(wb_eQTL2,sub_eQTL)
    
  }
}

#checker si bien 
nrow(wb_eQTL2) #2398562
wb_eQTL2<-wb_eQTL2[order(chr,pos)]
fwrite(wb_eQTL2,"../../ref/Whole_Blood.eQTL.csv",sep = ";")

#nb de count SNP par type
library(ggplot2)
ggplot(data = wb_eQTL2,aes(x=type))+geom_bar() 

#remove na et save : 
wb_eQTL2<-wb_eQTL2[!is.na(type)]
wb_eQTL3<-wb_eQTL2[!is.na(type)][,1:(ncol(wb_eQTL2)-1)]

fwrite(wb_eQTL3,"../../ref/Whole_Blood.eQTL.csv",sep = ";")
wb_eQTL<-wb_eQTL3
rm(wb_eQTL2,wb_eQTL3)
#SNP signif enrichit en 4 et 6 ?
wb_eQTL_count<-wb_eQTL[,n.snp:=.N,by=type][!duplicated(n.snp),c('type','n.snp')]
wb_eQTL_count<-wb_eQTL_count[!is.na(type)][order(type)]
wb_eQTL_count$n.region<-table(features$type)
features[,len.tot:=sum(longueur),by=type]
wb_eQTL_count$length.tot<-features[!duplicated(len.tot)][order(type)][,"len.tot"]
wb_eQTL_count[,cpm.snp:=n.snp/length.tot*1000000]
ggplot(data = wb_eQTL_count,aes(x=type,y=cpm.snp))+geom_col()+ggtitle("sig. blood eQTL / CD34+ Chromatin feature (in cpm)") #count par million SNP par type
fwrite(wb_eQTL_count,"summary_eQTL_Features.csv",sep = ";")

#next : eQTR : region regulé par df
eQTR<-features[,.(chr,start,end,type,longueur,chrReg)]

#1) add region lvl : nEQTL dans la region
wb_eQTL<-fread("../../ref/Whole_Blood.eQTL.csv",sep = ";",dec = ",")
wb_eQTL<-wb_eQTL[chrRegion!=""]
wb_eQTL[,nQTL.reg:=.N,chrRegion]

#add genelvl : 
#1) col gene : n.snp.reg.gene | pval_beta | 
wb_eQTL[,nQTL.gene:=.N,gene_id] #bonus mais pas forcément utile

wb_eQTL[,nQTL.reg.gene:=.N,by=.(chrRegion,gene_id)]
wb_eQTL

#2)SNP :    
#moy(pval_nom-threshold)
wb_eQTL[,avg.pval.thr:=mean(pval_nominal-pval_nominal_threshold),by=.(chrRegion,gene_id)]

#| mean(distDuTSS)/sd
wb_eQTL[,avg.distTSS:=mean(tss_distance),by=.(chrRegion,gene_id)]
wb_eQTL[,sd.distTSS:=sd(tss_distance),by=.(chrRegion,gene_id)]
#| distr.distDuTSS 
wb_eQTL[,distr.distTSS:=paste(summary(tss_distance)[-4],collapse="/"),by=.(chrRegion,gene_id)]

#| mean(posInReg)/sd 
wb_eQTL[,avg.pos:=mean(pos),by=.(chrRegion,gene_id)]
wb_eQTL[,sd.pos:=sd(pos),by=.(chrRegion,gene_id)]
#|distr(posInReg)
wb_eQTL[,distr.pos:=paste(summary(pos)[-4],collapse="/"),by=.(chrRegion,gene_id)]
#| distr(posInRegAllQTL)
wb_eQTL[,distr.pos.All:=paste(summary(pos)[-4],collapse="/"),chrRegion]
wb_eQTL
#| moy(distsEntreEux)/sd
dist_all<-function(vec){
  vec_dist<-sapply(1:(length(vec)-1), function(i){
    
    return(abs(vec[i]-vec[i+1:(length(vec)-i)]))
    
  })
  return(unlist(vec_dist))
}
dist_all(tail(wb_eQTL)$pos)

wb_eQTL[,avg.distBtwThem:=mean(dist_all(pos)),by=.(chrRegion,gene_id)]
wb_eQTL[,sd.distBtwThem:=sd(dist_all(pos)),by=.(chrRegion,gene_id)]
#| moy(distsAllQTL de la reg)/sd
wb_eQTL[,avg.distBtwAll:=mean(dist_all(pos)),by=.(chrRegion)]
wb_eQTL[,sd.distBtwAll:=sd(dist_all(pos)),by=.(chrRegion)]

#clear
eQTR<-unique(wb_eQTL,by = c("chrRegion","gene_id"))
eQTR<-eQTR[,.(chrRegion,type,nQTL.reg,gene_id,
              gene,pval_beta,nQTL.reg.gene,
              avg.pval.thr,avg.distTSS,sd.distTSS,distr.distTSS,
              avg.distBtwThem,sd.distBtwThem,
              avg.distBtwAll,sd.distBtwAll,
              avg.pos,sd.pos,
              distr.pos,
              distr.pos.All
              )]
eQTR[,chrReg:=chrRegion]
eQTR<-eQTR[,-c("chrRegion")]
features<-features[,.(chr,start,end,chrReg,longueur)]
eQTR<-merge(features,eQTR,by=c("chrReg"))
eQTR
#number of gene by region 
eQTR[,nGenes:=.N,by=.(chrReg)]
fwrite(eQTR,"../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv",sep = ";")

eQTR<-fread("../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv",sep = ";")

hist(eQTR$nGenes,breaks = 50) #majority of one gene by region,, but 100k region avec 2 genes, 75k avec 3genes...
hist(eQTR[type%in%c(4,6)]$nGenes,breaks = 50) 
#split region : si 2 clusters de SNPs distincts, en coupant au milieu de maxPosSNPGene1 et minPosSNPGene2 
eQTR1<-eQTR[nGenes==1]
eQTR2<-eQTR[nGenes>1]
eQTR2
#3 choses : 1) 1 SNP peut lier 2 genes (ex : R3HCC1L, et LOXL4)
#2)  1 SNP contenu dans un cluster de SNP liant un gene A peut lier un gene B
#3) les genes Symbols sont mal annoter
#nb de genes sans symbol
nrow(unique(eQTR,by = "gene_id")[gene==""]) #2300

ensemblSansSymbol<-unique(eQTR,by = "gene_id")[gene==""]$gene_id #2300
head(ensemblSansSymbol)
length(ensemblSansSymbol) #2300
library(biomartr)

library(stringr)
biomartr::getMarts()
dtset<-biomartr::getDatasets("ENSEMBL_MART_ENSEMBL")[]
dtset[str_detect(dtset$dataset,"hsapiens"),]
biomartr::getAttributes(mart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
result_BM <- biomartr::biomart( genes      = ensemblSansSymbol, # genes that we wanted info
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = "hgnc_symbol", # attributes were selected with biomartr::getAttributes()
                                filters    = "ensembl_gene_id_version")
#put in gene colonne
result_BM<-result_BM[result_BM$hgnc_symbol!="",]
nrow(result_BM) #☺ajout 712 symbol
for(i in 1:nrow(result_BM)){
  ensembl<-result_BM$ensembl_gene_id_version[i]
  symbol<-result_BM$hgnc_symbol[i]
  eQTR[gene_id==ensembl,"gene"]<-symbol
}
eQTR
fwrite(eQTR,"../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv",sep = ";")
eQTR<-fread("../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv")
eQTR
sapply(eQTR, class)
eQTR<-eQTR[,type:=as.factor(type)] #put type in factor instead of numeric
eQTR[,start.eQTR:=sapply(distr.pos,function(x){
  return(as.numeric(strsplit(x,"/")[[1]][1]))
})]
eQTR
eQTR[,end.eQTR:=sapply(distr.pos,function(x){
  return(as.numeric(strsplit(x,"/")[[1]][5]))
})]
eQTR[,length.eQTR:=end.eQTR-start.eQTR]
eQTR


#exploration des régions
eQTR1<-eQTR[nGenes==1]
eQTR2<-eQTR[nQTL.reg.gene>1]
eQTR2
plot(density(eQTR2$length.eQTR)),log="x") #pique à 
abline(v=quantile(eQTR2$length.eQTR,0.25))
median(eQTR2$length.eQTR) #1024pb => region +/-500pb d'un SNP

eQTR[length.eQTR<1000,start.eQTR2:=start.eQTR-(500-floor(length.eQTR/2))]
eQTR[length.eQTR<1000,end.eQTR2:=end.eQTR+(500-ceiling(length.eQTR/2))]

eQTR[,length.eQTR2:=end.eQTR2-start.eQTR2]
eQTR
fwrite(eQTR,"../../ref/CD34_chromatin_feature_region_annotated_with_BloodeQTL.csv",sep = ";")
#splitte région : 
#dans featuresReg, si eQTR <



#région 
#A Remove Car useless :
# #clivage region possible pour cluster de SNP A et B (median A < median B, et nA<nB) si :
# #1) maxPos A < minPos B
# #2) nA*-log10(pval) > 1/5e(nB*-log10(pval))
# eQTR2<-eQTR2[order(chr,start,avg.pos)]
# #1) maxPos A < minPos B
# #a) tester avec une region
# 
# #b) need a fonction ~= function dist_all()

#add CpG 
#ajout col CpG in this eQTR_df : |CpGid|pos|chr|typePrev(doit etre==type) | genePred |TSSdist | EnsemblAnnot
# annot<-fread("../../ref/annotation_CpG_HELP_ALL_070420.csv",dec = ",")
# annot[,locisID:=V1]
# annot<-annot[,-c(1)]
# ord<-names(annot)[c(length(annot),1:(length(annot)-1))]
# annot<-annot[,..ord]
# fwrite(annot, "../../ref/annotation_CpG_HELP_ALL_070420.csv",sep = ";")
annot<-fread("../../ref/annotation_CpG_HELP_ALL_070420.csv")
annoth<-head(annot,100)
#1) need to add chrineRegion for merge(eQTR,annot, by=c("chrReg"))
#test
for(chrom in unique(annoth$chr)){
  sub_eQTR<-eQTR[chr==chrom]
  annoth[,chrReg:=findReg(start,eQTR)sub_eQTR$chrReg[which(sub_eQTR$start<=start&sub_eQTR$end>start)],by="locisID"]
  
}

findReg<-function(poss,regions.dt){
  
  return()
}
annoth 



#check si bien annoté :
df<-fread("analyses/eQTL_test/2020-05-04annotRegchr12part10.txt")
df[is.na(chrRegion),] #oui, pas de na
sum(df$chrRegion=="")
df[chrReg!=""]
#work
annot2<-df
#collapse sub_annot
for(chr in paste0("chr",c(1:22,"X"))){
  print(chr)
  for(file in list.files("analyses/eQTL_test/",pattern = paste0(chr,"part"))){
    print(strsplit(file,chr)[[1]][2])
    sub_annot<-fread(paste0("analyses/eQTL_test/",file))
    annot2<-rbind(annot2,sub_annot)
    
  }
}
annot2
annotsToKeep<-c("LocisID","start","gene",)
eQTR.CpG<-merge(eQTR,annot2[,..annotsToKeep], by=c("chrReg"))
#