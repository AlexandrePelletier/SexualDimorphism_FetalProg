
#get eQTL info : links SNP with RNA expr

#config et library
options(stringsAsFactors=F)
set.seed(12345)
library(limma)
source("scripts/scRNA-seq integration.R")
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

#! faire df region reg : region stricte :si chevauchement. + Large si region sep de 90 pb

write.csv2(wb_eQTL,"../../ref/Whole_Blood.eQTL.csv")
