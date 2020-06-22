#2020-06-19 
#PATHWAYS ANALYSIS
#Based on GeneScore, i want perform pathway and gwas analysis to determine master gene regulator for future biological validation.

#let's first focus on Female Genes : 
library(data.table)
library(stringr)
library(clusterProfiler)
library(pathview)
library(ggplot2)
resF<-fread("analyses/withoutIUGR/2020-06-03_CF.LF_CpG_Gene_scorev3.2.csv")
source("scripts/utils.R")

#lets defined threshold ~ best compromise between keeping genes with hi Pvalue and genes with weak.
resF[,noSigGene:=all(pval>10e-3),by="gene"]
resF[,SigGene:=any(pval<10e-5),by="gene"]
plot(density(resF$GeneScore))
lines(density(unique(resF[noSigGene==TRUE],by="gene")$GeneScore),col=2)

lines(density(unique(resF[SigGene==TRUE],by="gene")$GeneScore),col=3)
#5there is pva >10^-3 but with hi g
abline(v=200)
gF<-unique(resF[GeneScore>200]$gene)
length(gF) #1123

kF <- enrichKEGG(gene         = bitr(gF,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05)
dkF<-as.data.frame(kF)
nrow(dkF) #13
dkF
dkF<-data.table(dkF)
dkF[,GeneScore.avg:=mean(unique(resF,by="gene")$GeneScore[unique(resF,by="gene")$gene %in% tr(geneID,tradEntrezInSymbol = T)],na.rm=T),.(ID)]
dkF
dotplot(kF,x=dkF$GeneScore,showCategory=13)
emapplot(kF)
enrichplot::upsetplot(kF,n=13)
#3 things :  
#1) MAPK have 13 genes exclusif, 
#2) 6 genes are shared between HPI, plurisigpathway, cush syndr, gastric cancer, hippo sig path, and basal cell carcinoma
#3) sig pathway and MODY share 2 exclu genes.
#what are this genes and plot pathview with genescore ?
dkF[,genes:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "/"),by="ID"]
#1)
genes1<-setdiff(tr(dkF[ID=="hsa04010"]$genes),unique(unlist(lapply(dkF[ID!="hsa04010"]$genes,tr))))
genes1
# [1] "MAPKAPK3" "DAXX"     "RPS6KA2"  "DUSP4"    "TRAF2"    "MAP3K8"   "IGF2"     "MAP3K11" 
# [9] "CACNB3"   "MAP2K5"   "TAOK2"    "NFATC3"   "NFATC1"   "MKNK2" 
#first make geneList of genescore
geneList<-unique(resF,by="gene")$GeneScore
names(geneList)<-unique(resF,by="gene")$gene
geneList<-sort(geneList,decreasing = T)
genes.df<-bitr(names(geneList),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneList[genes.df$SYMBOL]
geneList.Entrez<-genes.df$GeneScore
names(geneList.Entrez)<-genes.df$ENTREZID
pathview(gene.data  = geneList.Entrez,
         pathway.id = "hsa04010", 
         species    = "hsa",
         limit      = list(gene=max(abs(geneList.Entrez)), cpd=1),
         kegg.dir ="analyses/withoutIUGR/pathwaysMap")


#2020-06-22 ; trouver sex spe signif pathways /genes 
#I) by path enrichment of sex spé sig epigen affected genes 
#a) first, deter genescore cutoff x / deter EpigenAffGene (EAG) : 
# 1000 permuts of all samples CTRL / LGA => LIMMA => pval and FC => df GeneScore => save pct genescore permuts >= geneScore obs
#make a permutFunction
batch<-fread("../../ref/batch_CD34_library_date_090420.csv",dec = ",")

batch[,sample:=V1][,Group_Sex:=paste(Group_name,Gender,sep="_")]

cols<-names(batch)[c(ncol(batch),3:(ncol(batch)-1))]
batch<-batch[,..cols]
batch<-batch[,-c("newComplexityScore","newComplexityScore2","newComplexityScore3","newComplexityScore4","SeqDepthHi")]
fwrite(batch,"../../ref/cleaned_batch_CD34_library_date_220620.csv",sep = ";")

batch<-batch[Group_name%in%c("C","L")]
batch

varNumToModel<-c("Mat.Age")
varFacToModel<-c("Group_Sex",'batch',"latino","Group_Complexity_Fac")
varToModel<-c(varNumToModel,varFacToModel)

sample_F<-batch$sample[!(apply(is.na(batch[,..varToModel]),1,any))]
length(sample_F) #70
batch<-batch[sample%in%sample_F,c("sample",..varToModel),]
batch[,Mat.Age:=as.numeric(Mat.Age)]

batch[,Group_Sex:=as.factor(Group_Sex)]
batch[,batch:=as.factor(batch)]
batch[,latino:=as.factor(latino)]
batch[,Group_Complexity_Fac:=as.factor(Group_Complexity_Fac)]
str(batch)
batch

methyl_data<-fread("../../ref/2020-05-25_methyl_data_before_limma.csv")

methyl_df<-data.frame(methyl_data)
rownames(methyl_df)<-methyl_df$locisID
head(methyl_df)
formule<- ~0 + Group_Sex  + batch  + latino + Mat.Age + Group_Complexity_Fac
compas<-c("Group_SexC_M-Group_SexL_M","Group_SexC_F-Group_SexL_F")
library(stringr)
abbrev_compas<-str_remove_all(compas,"Group_Sex")
cpg.regs_ref<-fread("../../ref/2020-06-03_All_CpG-Gene_links.csv")

MethAnalysisExec<-function(methyl_df,batch_F,formule,compas,abbrev_compas,cpg.regs_ref){
  library(limma)
  design<-model.matrix(formule,data = batch_F)
 
  fit <- lmFit(methyl_df[,batch_F$sample], design)  
  
  cont.matrix <- makeContrasts(contrasts = compas,
                               levels=design)
  colnames(cont.matrix)<-abbrev_compas
  
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2) 
  res_list<-lapply(abbrev_compas,function(compa){
    res<-topTable(fit2,coef=compa,n =Inf)
    res<-CalcGeneScore(res,cpg.regs_ref,sumToGene=T)
    res<-res[order(gene)]
    genescore<-res$GeneScore
    names(genescore)<-res$gene
    return(genescore)
  })
  
  names(res_list)<-abbrev_compas
   return(res_list)
  
}

CalcGeneScore<-function(res,cpg.regs_ref,sumToGene=F){
  res<-data.table(locisID=as.numeric(rownames(res)),
                  meth.change=res$logFC,
                  pval=res$P.Value)[order(locisID)]
  
  res<-merge(res,cpg.regs_ref,by="locisID")
  res[,CpGScore:=((-log10(pval)/4)*meth.change)*RegWeight*LinksWeight]
  res[,nCpGWeight:=(1/sum(1/(abs(CpGScore)+1)))^(1/3),by="gene"]
  
  res[,GeneScore:=sum(CpGScore)*nCpGWeight,by="gene"]
  
  #bonus annot
  res[,nCpG.Gene:=.N,by=.(gene)]
  res[,nCpGSig.Gene:=sum(pval<10^-3),by=.(gene)]
  
  if(sumToGene){
    return(unique(res[order(-GeneScore,pval)],by="gene"))
  }else{
    return(res[order(-GeneScore,pval)])
  }
}

#test CalcGeneScore
design<-model.matrix(formule,data = batch)

fit <- lmFit(methyl_df[,batch$sample], design)  

cont.matrix <- makeContrasts(contrasts = compas,
                             levels=design)
colnames(cont.matrix)<-abbrev_compas

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2) 

res<-topTable(fit2,coef=abbrev_compas[2],n =Inf)
res<-CalcGeneScore(res,cpg.regs_ref,sumToGene=T)

res
plot(res$nCpG.Gene,res$nCpGWeight)
summary(res$nCpGWeight)
plot(res$nCpG.Gene,res$GeneScore)
plot(res[nCpG.Gene<20]$nCpG.Gene,res[nCpG.Gene<20]$GeneScore)
res[nCpGWeight>2]
res[gene=="NUP85"&pval<0.01]
res[locisID==1824925]


#test limma exec 
res_list<-MethAnalysisExec(methyl_df, batch,formule, compas,abbrev_compas,cpg.regs_ref )

res_list



MethAnalysisPermut<-function(batch_F,col_to_permut,n,sample_col="sample",cpg.regs_ref=NULL,formule=NULL,seed=123){
  set.seed(seed)
  library(data.table)
  
  
  if(is.null(cpg.regs_ref)){
    cpg.regs_ref_path<-"../../ref/2020-06-03_All_CpG-Gene_links.csv"
    print(paste("recup cpg.gene regul annotation in ",cpg.regs_ref_path))
    cpg.regs_ref<-fread("../../ref/2020-06-03_All_CpG-Gene_links.csv")
  }
  
  if(is.null(formule)){
    
    
    print("limma model selon cette formule :")
    print(formule)
  }
  
  print("Limma execution of observed batch..")
  res<-LimmaExec(methyl_data,batch_F,formule)
  res<-CalcGeneScore(res,cpg.regs_ref)
  resGenes<-unique(res[,.(gene,GeneScore)][order(gene,pval)],by="gene")
  resGenes_all<-data.table(resGenes)
  print(paste("PERMUTATION of",col_to_permut,n,"times..."))
  
  
  batch1<-data.table(batch)
  for(i in 1:n){
    print(paste("permutation",i))
    
    batch1[,Group_Sex:=sample(Group_Sex)]
    
    print(" ==> limma execution")
    res<-LimmaExec(methyl_data,batch_F,formule)
    print(" ==> Gene Score calculation")
    
    
  }
  print("END PERMUTATION")
  
  print("save  pct genescore permuts >= geneScore obs")
  
  
  return(batch[,])
  
}

#=> gene score cutoff = Accept Gene in EAF if appeared < x%  => genesF and genesM spe
#b) deter sex spé sig genes   : 
#=> 100 permutation of all samples CTRL / LGA => LIMMA => pval and FC => gene score => genescore cutoff of x (for each save nb of genes > x)
#=> setdiff(genesF, genesM) and inversely => save genes spéF and M => all genes obs dans <5 permut => sex spé sig genes

#b) pathway enrichment of this sex spé genes

#II) by path enrichment of of sig affected genes
