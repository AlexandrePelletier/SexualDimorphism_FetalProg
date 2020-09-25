
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


# 2020-09-16 Make list of CpG
library(data.table)
library(stringr)
#I) CpG with DEG interesting


#DEG
resE<-list(
  "cl"=fread("../singlecell/analyses/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv"),
  "cflf"=fread("../singlecell/analyses/04-DEG_in_LGA/2020-09-07_pseudo_bulk_DEseq2_LgaFVsCtrlF_CBP1and3andCbp558andMaleandIugr_samples_excluded_regr_batch_all_genes.csv"),
  "cl_hsc"=fread("../singlecell/analyses/04-DEG_in_LGA/2020-09-07_pseudo_bulk_HSC_MPP_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv"),
  "cflf_hsc"=fread("../singlecell/analyses/04-DEG_in_LGA/2020-09-07_pseudo_bulk_HSC_MPP_DEseq2_LgaFVsCtrlF_CBP1and3andCbp558andMaleandIugr_samples_excluded_regr_batch_all_genes.csv")
)

lapply(names(resE),function(comp)resE[[comp]][,compa:=strsplit(..comp,"_")[[1]][1]])
lapply(names(resE),function(comp)resE[[comp]][,subpop:=ifelse(str_detect(..comp,"_"),"hsc/mpp","all_progen")])
resE<-Reduce(function(x, y) merge(x[,.(gene,baseMean,log2FoldChange,pvalue,padj,compa,subpop)],
                                  y[,.(gene,baseMean,log2FoldChange,pvalue,padj,compa,subpop)], all = T),resE)



#meth
resM<-list(
  "cl"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_C.L.csv"),
  "cfcm"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CF.CM.csv"),
  "cflm"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CF.LM.csv"),
  "cflf"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CF.LF.csv"),
  "cmlf"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CM.LF.csv"),
  "cmlm"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CM.LM.csv"),
  "lmlf"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_LM.LF.csv")
)
lapply(names(resM),function(comp)resM[[comp]][,compa:=..comp])
resM<-Reduce(function(x, y) merge(x, y, all = T),resM)


resM[,cpg_of_interest:=abs(meth.change>20)&pval<0.005&(abs(tss_dist)<10000|in_eQTR==TRUE)&type!=0]



distFromNext<-function(x){
  x_sorted<-sort(x)
  dists<-rep(NA,length(x_sorted))
  for(i in 1:(length(x_sorted)-1)){
    dists[i]<-abs(x_sorted[i]-x_sorted[i+1])
  }
  return(dists[match(x,x_sorted)])
}
distFromPrev<-function(x){
  x_sorted<-sort(x,decreasing =T )
  dists<-rep(NA,length(x_sorted))
  for(i in 1:(length(x_sorted)-1)){
    dists[i]<-abs(x_sorted[i]-x_sorted[i+1])
  }
  return(dists[match(x,x_sorted)])
}



resM[,dist_from_next_cpg:=distFromNext(tss_dist),by=c('gene',"compa")]
resM[,dist_from_prev_cpg:=distFromPrev(tss_dist),by=c('gene',"compa")]
resM[,cpg_close_to_int:=c(FALSE,c(dist_from_next_cpg<300&cpg_of_interest==TRUE)[-.N])|
       c(c(dist_from_next_cpg<300&cpg_of_interest==TRUE,FALSE)[-1]),by=c('gene',"compa")]
resM[cpg_of_interest==T|cpg_close_to_int==TRUE]

resM[,cpg_close_to_int_int:=cpg_close_to_int&pval<0.05&abs(meth.change)>20]

res_of_int<-merge(resM[cpg_of_interest==T|cpg_close_to_int_int==T],resE,by=c("gene","compa"),all.x = T)

res_of_int<-res_of_int[,.(locisID,chr,pos,gene,compa,meth.change,pval,cpg_of_interest,cpg_close_to_int_int,tss_dist,in_eQTR,type,feature_type_name,
                          GeneScore,gene,compa,subpop,baseMean,log2FoldChange,pvalue,padj,
                          dist_from_next_cpg,dist_from_prev_cpg)]
res_of_int[order(-GeneScore,pval)]


fwrite(res_of_int[order(gene,-GeneScore,pval)],"analyses/model14_without_iugr/2020-09-17_cpgs_of_interest_and_cpg_close_pval0.005_meth.change20_tss_dist10k_or_ineQTR_type123456.csv",sep=";")

#2020-09-25 make cpg list of int of genes of int (deg, module)
library(data.table)

degs<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_padj0.05_DEG.csv")
degs[abs(log2FoldChange)>0.5]
meth<-fread("analyses/model14_without_iugr/2020-09-24_all_res_with_perm.csv")
methg<-unique(meth[order(gene,compa,pval)],by=c("gene","compa"))
demet<-merge(degs[,.(gene,baseMean,log2FoldChange,pvalue,padj)],meth,by="gene")
genes_DEG<-unique(demet[abs(log2FoldChange)>0.5&abs(GeneScore)>60]$gene)

hormon<-fread("analyses/genes_of_interest/2020-09-24_genes_hub_modul_hormon.csv")
genes_hormon<-unique(hormon[abs(GeneScore)>60]$gene)

longev<-fread("analyses/genes_of_interest/2020-09-24_genes_hub_modul_longevity.csv")

genes_longev<-unique(longev[abs(GeneScore)>60]$gene)

stem<-fread("analyses/genes_of_interest/2020-09-25_genes_hub_modul_stem.csv")
genes_stem<-unique(stem[abs(GeneScore)>60]$gene)

nut<-fread("analyses/genes_of_interest/2020-09-25_genes_hub_modul_nutrient.csv")

genes_nutrient<-unique(nut[abs(GeneScore)>60]$gene)

prolif<-fread("analyses/genes_of_interest/2020-09-25_genes_hub_modul_prolif.csv")
genes_prolif<-unique(prolif[abs(GeneScore)>60]$gene)
genes_dint<-unique(c(genes_DEG,genes_hormon,genes_longev,genes_nutrient,genes_prolif,genes_stem))
length(genes_dint)
methg[,DEG:=gene%in%genes_DEG]
methg[,hormon:=gene%in%genes_hormon]
methg[,longev:=gene%in%genes_longev]
methg[,nutrient:=gene%in%genes_nutrient]
methg[,prolif:=gene%in%genes_prolif]
methg[,stem:=gene%in%genes_stem]
fwrite(methg[gene%in%genes_dint],"analyses/genes_of_interest/2020-09-25_all_genes_of_interest.csv",sep=";")

meth[,cpg_of_interest:=abs(meth.change>20)&pval<0.005&(abs(tss_dist)<10000|in_eQTR==TRUE)&type!=0]



distFromNext<-function(x){
  x_sorted<-sort(x)
  dists<-rep(NA,length(x_sorted))
  for(i in 1:(length(x_sorted)-1)){
    dists[i]<-abs(x_sorted[i]-x_sorted[i+1])
  }
  return(dists[match(x,x_sorted)])
}
distFromPrev<-function(x){
  x_sorted<-sort(x,decreasing =T )
  dists<-rep(NA,length(x_sorted))
  for(i in 1:(length(x_sorted)-1)){
    dists[i]<-abs(x_sorted[i]-x_sorted[i+1])
  }
  return(dists[match(x,x_sorted)])
}



meth[,dist_from_next_cpg:=distFromNext(tss_dist),by=c('gene',"compa")]
meth[,dist_from_prev_cpg:=distFromPrev(tss_dist),by=c('gene',"compa")]
meth[,cpg_close_to_int:=c(FALSE,c(dist_from_next_cpg<300&cpg_of_interest==TRUE)[-.N])|
       c(c(dist_from_next_cpg<300&cpg_of_interest==TRUE,FALSE)[-1]),by=c('gene',"compa")]
meth[cpg_of_interest==T|cpg_close_to_int==TRUE]

meth[,cpg_close_to_int_int:=cpg_close_to_int&pval<0.05&abs(meth.change)>20]

res_of_int<-merge(meth[cpg_of_interest==T|cpg_close_to_int_int==T],methg[gene%in%genes_dint])


fwrite(res_of_int[,.(locisID,chr,pos,pval,meth.change,cpg_of_interest,cpg_close_to_int,dist_from_next_cpg,gene,GeneScore,compa,DEG,stem,prolif,hormon,longev,nutrient)][order(gene,-GeneScore,pval)],
       "analyses/genes_of_interest/2020-09-25_all_cpgs_of_interest_and_cpg_close_in_genes_of_interest_pval0.005_meth.change20_tss_dist10k_or_ineQTR_type123456.csv",sep=";")




