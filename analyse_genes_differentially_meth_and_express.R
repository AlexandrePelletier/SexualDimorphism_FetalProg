
#2020-08-12 explore GeneSig
#visual res of LGAvsCTRL

library(data.table)
res_meth1<-fread("../../../Alexandre_SC/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")

# show enrichment Genscore in DEG
#a) boxplot
library(ggplot2)
library(patchwork)
res_meth1[,DEG:=padj<0.05]
p1<-ggplot(res_meth1[!is.na(DEG)])+geom_boxplot(aes(x=DEG,y=GeneScore))
p2<-ggplot(res_meth1[!is.na(DEG)])+geom_boxplot(aes(x=DEG,y=meth.change))
p3<-ggplot(res_meth1[!is.na(DEG)])+geom_boxplot(aes(x=DEG,y=pval))+scale_y_log10()
p1+p2+p3
#b) plot FC ~ GeneScore, DEG vs all
ggplot(res_meth1[!is.na(DEG)])+geom_point(aes(x=log2FoldChange,y=GeneScore))+facet_wrap(vars(padj<0.2))
d1<-ggplot(res_meth1[padj<0.2])+geom_point(aes(x=log2FoldChange,y=GeneScore))
d2<-ggplot(res_meth1[padj<0.2])+geom_point(aes(x=log2FoldChange,y=meth.change))

d1+d2
library(ggrepel)
ggplot(res_meth1[padj<0.2],aes(log2FoldChange,GeneScore))+
  geom_point()+
  geom_label_repel(aes(label = ifelse(abs(GeneScore)>60&abs(log2FoldChange)>0.8,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+theme_classic()

ggplot(res_meth1[padj<0.2&pval<0.01],aes(log2FoldChange,meth.change))+
  geom_point()+
  geom_label_repel(aes(label = ifelse(abs(meth.change)>20&abs(log2FoldChange)>0.8,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+theme_classic()

#un gene a un methchange neg mais genescore pos : 
res_meth1[meth.change<0&padj<0.2&GeneScore>0] #KLF4


#analyse GENE BY GENE
resF<-fread("analyses/withoutIUGR/2020-07-30_resCF_LF.csv")
resF[gene=="KLF4"&pval<0.01]
resM<-fread("analyses/withoutIUGR/2020-07-30_resCM_LM.csv")
resM[gene=="KLF4"&tss_dist>0&tss_dist<1000][order(tss_dist)] #hypermeth sex spé

resF[,compa:="CLF"]
resM[,compa:="CLM"]
res<-merge(resF,resM,all = T)
res<-merge(res_meth1[,.(gene,baseMean,log2FoldChange,padj)],res,by="gene",all.x = T)

#♦++ hypermet et downreg
res[GeneScore>50&log2FoldChange<(-0.6)&padj<0.2&pval<0.01&meth.change>0,hyperMetDownReg:=(.N>2),by="gene"]
res[hyperMetDownReg==T]
res[hyperMetDownReg==T,changeSexSpe:=length(unique(compa))==1,by="gene"]
res[hyperMetDownReg==T&changeSexSpe==T]
res[hyperMetDownReg==T,gene_of_interest:=T]
#♦hypermet et upreg
res[GeneScore>30&log2FoldChange>0.6&padj<0.2&pval<0.05&meth.change>0]
res[GeneScore>30&log2FoldChange>0.6&padj<0.2&pval<0.05&meth.change>0,hyperMetUpReg:=(.N>2),by="gene"]
res[hyperMetUpReg==T]#LGALS
re
#just 1 cpg pval<0.01 for lgals1, interessant qd meme ?
source("scripts/utils/methyl_utils.R")
plotMeth(res[compa=="CLF"&gene=="LGALS1"&pval<0.05]$locisID,plot = "jitter")

res[hyperMetUpReg==T,changeSexSpe:=length(unique(compa))==1,by="gene"]
res[hyperMetUpReg==T&changeSexSpe==T]
res[hyperMetUpReg==T,gene_of_interest:=T]


#♦hypomet et upreg
res[GeneScore<(-50)&log2FoldChange>(0.6)&padj<0.2&pval<0.01,hypoMetUpReg:=(.N>2),by="gene"]
res[hypoMetUpReg==T]


#♦hypomet et downreg
res[GeneScore<(-50)&log2FoldChange<(-0.6)&padj<0.2&pval<0.01,hypoMetDownReg:=(.N>2),by="gene"]
res[hypoMetDownReg==T]

#♦diffmet et diffreg
res[abs(GeneScore)>50&(sign(meth.change)!=sign(GeneScore))&abs(meth.change)>20&abs(log2FoldChange)>0.6&padj<0.2&pval<0.05]
res[abs(GeneScore)>50&(sign(meth.change)!=sign(GeneScore))&abs(meth.change)>20&abs(log2FoldChange)>0.6&padj<0.2&pval<0.05,hyper_and_hypoMet_DEG:=T]
res[hyper_and_hypoMet_DEG==T,changeSexSpe:=length(unique(compa))==1,by="gene"]
res[hyper_and_hypoMet_DEG==T,gene_of_interest:=T]
res[gene_of_interest==T]$gene
fwrite(res[gene_of_interest==T],"analyses/withoutIUGR/genes_of_interest/cpg_genes_of_interest.csv",sep=";")

for(gen in unique(res[gene_of_interest==T]$gene)){
  print(gen)
  plotMeth(res[compa=="CLF"&gene_of_interest==T&gene==gen]$locisID)
  
  ggsave(paste0("analyses/withoutIUGR/genes_of_interest/",gen,"_boxplot.png"))
  plotMeth(res[compa=="CLF"&gene_of_interest==T&gene==gen]$locisID,plot = "jitter")
  ggsave(paste0("analyses/withoutIUGR/genes_of_interest/",gen,"_jitter.png"))
  
}



