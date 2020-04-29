
#get eQTL info : links SNP with RNA expr

#config et library
options(stringsAsFactors=F)
set.seed(12345)
#library(limma)
source("scripts/utils.R")

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
#length region par type
library(ggplot2)
p<-ggplot(data = features,aes(x=type))
p+geom_boxplot(aes(y=longueur))+scale_y_continuous(trans = 'log10')  
p+geom_bar()
# les regions "0" majoriatire (~300k reg contre ~150k region pour les autre), correspondant a l'heterochromatine,
#sont de large région, (75% entre 3k et 5k), les autres sont à 2k de longueur


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
      sub_features<-features[chr==chrom&start<=q1&end>q2,]
      for(j in 1:nrow(sub_features)){
        reg<-sub_features[j,]
        sub_eQTL[pos>=reg$start&pos<=reg$end,"type"]<-reg$type
        
      }
      q2<-q1
      fwrite(sub_eQTL,file = paste0(output,"annotEnsembl",chrom,"part",i,".txt"))
      print("..annotes et enregistrer")
    }
    
  }
}

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
      sub_features<-features[chr==chrom&start<=q1&end>q2,]
      for(j in 1:nrow(sub_features)){
        reg<-sub_features[j,]
        sub_eQTL[pos>=reg$start&pos<=reg$end,"type"]<-reg$type
        
      }
      q2<-q1
      fwrite(sub_eQTL,file = paste0(output,"annotEnsembl",chrom,"part",i,".txt"))
      print("..annotes et enregistrer")
    }
    
  }
}
#check si bien annoté :
df<-fread("analyses/eQTL_test/2020-04-29annotEnsemblchr1part37.txt")
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
nrow(wb_eQTL2)
wb_eQTL2<-wb_eQTL2[order(chr,pos)]
fwrite(wb_eQTL2,"../../ref/Whole_Blood.eQTL.csv",sep = ";")
#SNP signif enrichit en 4 et 6 ?
ggplot(data = wb_eQTL2,aes(x=type))+geom_bar() #nb de count SNP par type

wb_eQTL2_count<-wb_eQTL2[,n.snp:=.N,by=type][!duplicated(n.snp),c('type','n.snp')]
wb_eQTL2_count<-wb_eQTL2_count[!is.na(type)][order(type)]
wb_eQTL2_count$n.region<-table(features$type)
features[,len.tot:=sum(longueur),by=type]
wb_eQTL2_count$length.tot<-features[!duplicated(len.tot)][order(type)][,"len.tot"]
wb_eQTL2_count[,cpm.snp:=n.snp/length.tot*1000000]
ggplot(data = wb_eQTL2_count,aes(x=type,y=cpm.snp))+geom_col()+ggtitle("eQTL significatif / features (en cpm)") #count par million SNP par type
fwrite(wb_eQTL2_count,"summary_eQTL_Features.csv",sep = ";")
