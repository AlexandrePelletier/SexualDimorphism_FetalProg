
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
res[abs(GeneScore)>50&(sign(meth.change)!=sign(GeneScore))&abs(meth.change)>20&abs(log2FoldChange)>0.6&padj<0.2&pval<0.05,hyper_and_hypoMet:=T]
res[hyper_and_hypoMet==T,changeSexSpe:=length(unique(compa))==1,by="gene"]
res[hyper_and_hypoMet==T,gene_of_interest:=T]
res[gene_of_interest==T]$gene

fwrite(res[gene_of_interest==T],"analyses/withoutIUGR/genes_of_interest/cpg_genes_of_interest.csv",sep=";")

res[gene_of_interest==T,nCpGCandidat:=.N,by="gene"]
cols<-c("changeSexSpe","hypoMetDownReg","hyperMetUpReg","hypoMetUpReg","hyperMetDownReg","hyper_and_hypoMet")

res[gene_of_interest==T,regulation:=paste(..cols[which(sapply(.SD,function(x)return(any(x==TRUE))))],collapse = "/"),by="gene",.SDcols=cols]


#genes in which KEGG pathways ? stemness, longevity, nutrient sensing pathway ?
source("../../../Alexandre_SC/scripts/utils/FindInfosGenes.R")
resg<-unique(res[gene_of_interest==T,
                 .(gene,baseMean,log2FoldChange,padj,compa,GeneScore,pval1000perm,nCpGCandidat,regulation)],
             by=c("gene","compa"))[order(-GeneScore)]

resg[,fonction:=sapply(gene,function(x)return(find_fonction(x,pathFromRoot = T)$Fonction))]
resg

resg[,pathway:=sapply(gene,function(x)return(find_KEGGPathway(x,pathFromRoot = T,searchInLocal = F)$Pathway))]
resg

fwrite(resg,
       "analyses/withoutIUGR/genes_of_interest/genes_of_interest.csv",sep=";")

#plot meth and expr in cells

for(gen in unique(res[gene_of_interest==T]$gene)){
  print(gen)
  plotMeth(res[compa=="CLF"&gene_of_interest==T&gene==gen]$locisID)
  
  ggsave(paste0("analyses/withoutIUGR/genes_of_interest/",gen,"_boxplot.png"))
  plotMeth(res[compa=="CLF"&gene_of_interest==T&gene==gen]$locisID,plot = "jitter")
  ggsave(paste0("analyses/withoutIUGR/genes_of_interest/",gen,"_jitter.png"))
  
}

library(Seurat)
cbps<-readRDS("../../../Alexandre_SC/analyses/2020-07-08_all_cbps.rds")
head(cbps@meta.data)
Idents(cbps)<-"seurat_clusters"
new.cluster.ids <- c("LMPP",
                     "HSC-SOCS3",
                     "MPP-RP",
                     "HSC-Er",
                     "HSC-EGR4",
                     "HSC-AVP",
                     "MEP",
                     "Ly P",
                     "HSC-CD164",
                     "GMP",
                     "EMP",
                     "proB",
                     "MPP-DC",
                     "ErP",
                     "HSC-HSPA5",  
                     "B cell",
                     "pDC-Cycle",
                     "div_cells",
                     "Ba/Eo/Mas",
                     "LT-HSC",
                     "T cell",
                     "Ba",
                     "Mo",
                     "NK")
names(new.cluster.ids)<-levels(cbps)
cbps<-RenameIdents(cbps,new.cluster.ids)
cbps[['diff_lvl']]<-cbps[['cell_type']]
cbps[['cell_type']]<-Idents(cbps)
DimPlot(cbps,label = T)
DefaultAssay(cbps)<-"RNA"
for(gen in unique(res[gene_of_interest==T]$gene)){
  print(gen)
  u<-FeaturePlot(cbps,gen,label = T,label.size = 2,repel = T,max.cutoff = "q75")
  v<-VlnPlot(cbps,gen,log = T)
  u/v
  ggsave(paste0("analyses/withoutIUGR/genes_of_interest/",gen,"_expr_plots.png"))
  
}


#genes longevity
longevityPaths<-c("AMPK signaling pathway","mTOR signaling pathway","Longevity regulating pathway",
                  "PI3K-Akt signaling pathway","Wnt signaling pathway","FoxO signaling pathway","Insulin signaling pathway",
                  "Ras signaling pathway","Rap1 signaling pathway","Regulation of lipolysis in adipocytes")
longevityPathsID<-unique(data.table(Reduce(rbind,list(data.table(as.data.frame(resF_KEGG.GSEA))[,.(ID,Description)],
                                                      data.table(as.data.frame(resFM_KEGG))[,.(ID,Description)])))[Description%in%longevityPaths]$ID)


