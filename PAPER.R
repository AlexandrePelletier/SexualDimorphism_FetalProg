library(stringr)
library(data.table)
library(clusterProfiler)
library(Seurat)
library(DESeq2)
library(limma)
library(ggplot2)
library(ggrepel)
#library(biomaRt)
library(org.Hs.eg.db)
source("scripts/utils/methyl_utils.R")

set.seed(1234)
options(stringsAsFactors = T)

#I) CTRLF vs CTRLM
dir.create("analyses/paper/ctrlM_vs_ctrlF",recursive = T)
# 1) GENESCORE
res<-fread("analyses/model14_without_iugr/2020-09-24_all_res_with_perm.csv")
# #add entrez id
# listEnsembl()
# ensembl <- useEnsembl(biomart = "ensembl")
# searchDatasets(mart = ensembl, pattern = "hsapiens")
# ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
# searchAttributes(mart = ensembl, pattern = "hgnc") #hgnc_symbol
# searchAttributes*(mart = ensembl, pattern = "entrez") #entrezgene_id
# genes<-getBM(attributes = c('hgnc_symbol', 'entrezgene_id'),
#       filters = 'hgnc_symbol',
#       values = unique(res$gene), 
#       mart = ensembl)
# genes<-data.table(genes)[,gene:=hgnc_symbol][,entrez:=entrezgene_id]
# genes[,is.true:=entrez==min(entrez),by="gene"]
# genes<-genes[is.true==T]
# res<-merge(res[,-"entrez"],genes[,.(gene,entrez)],all.x=T,by="gene")
# fwrite(res,"analyses/model14_without_iugr/2020-09-24_all_res_with_perm.csv")

resgC<-unique(res[compa=="cfcm"],by="gene")

# 2) PATHWAY
#   2a) Over-repres test
ggplot(resgC)+geom_density(aes(fill=pval1000perm<0.005,y=GeneScore))+geom_hline(yintercept = 50)
resgC[GeneScore>50&pval1000perm<0.01] #78 genes

keggC<-enrichKEGG(resgC[GeneScore>50&pval1000perm<0.01]$entrez,organism = "hsa",pvalueCutoff = 0.2 )
resK_C<-data.table(keggC@result)[order(pvalue)]
fwrite(resK_C,"analyses/paper/ctrlM_vs_ctrlF/2020-10-05_over_repr_kegg_ctrlFvsCtrlM_genescore50_pval0.01.csv",sep=";")

GO_C<-enrichGO(resgC[GeneScore>50&pval1000perm<0.01]$entrez,OrgDb = org.Hs.eg.db ,pvalueCutoff = 0.2 ,universe = as.character(resgC$entrez))
resGO_C<-data.table(GO_C@result)[order(pvalue)]
fwrite(resGO_C,"analyses/paper/ctrlM_vs_ctrlF/2020-10-05_over_repr_GO_ctrlFvsCtrlM_genescore50_pval0.01.csv",sep=";")

#   2b) GSEA (or similar), pval / perm
#transform en rank pos neg 
resgC[GeneScore>0,GeneScore_rnk:=rank(GeneScore)]
resgC[GeneScore<0,GeneScore_rnk:=-rank(abs(GeneScore))]
resgC[GeneScore==0,GeneScore_rnk:=0,by="compa"]
genelist<-resgC$GeneScore_rnk
names(genelist)<-resgC$entrez
kegg_gsea_C<-gseKEGG(sort(genelist,decreasing = T),organism = "hsa",minGSSize    = 15,maxGSSize=300 )
resK_gsea_C<-data.table(kegg_gsea_C@result)[order(pvalue)]
fwrite(resK_gsea_C,"analyses/paper/ctrlM_vs_ctrlF/2020-10-05_gsea_kegg_ctrlFvsCtrlM_GeneScoreRankPosNeg.csv",sep=";")

GO_gsea_C<-gseGO(sort(genelist,decreasing = T),OrgDb = org.Hs.eg.db ,minGSSize    = 15,maxGSSize=300)
resGO_gsea_C<-data.table(GO_gsea_C@result)[order(pvalue)]
fwrite(resGO_gsea_C,"analyses/paper/ctrlM_vs_ctrlF/2020-10-05_gsea_GO_ctrlFvsCtrlM_GeneScoreRankPosNeg.csv",sep=";")


#   2c) metalongevity genes 
long_ids<-c("hsa04211","hsa04310","hsa04910","hsa04010","hsa04152","hsa04068","hsa04150")
long_genes<-lapply(long_ids,getGenesKEGGPathw)
lg<-ul(long_genes)
length(lg)
saveRDS(lg,"analyses/paper/longevity_genes.rds")

stem_paths<-c("Wnt","TGF-beta","JAK","pluripotency","MAPK")
stem_ids<-sapply(stem_paths,function(path_word)unique(resK_gsea_C[str_detect(Description,path_word)]$ID))
stem_ids$JAK<-"hsa04630"
stem_genes<-lapply(stem_ids,getGenesKEGGPathw)

stem_hub<-lapply(stem_genes, function(genes){
  return(genes[sapply(genes, function(gene)sum(sapply(stem_genes, function(path)gene%in%path))>=length(stem_ids)/4)])
})
stem_hub<-as.vector(unique(unlist(stem_hub)))
length(stem_hub)#151
ggplot(resgC)+geom_boxplot(aes(x=compa,y=GeneScore,fill=gene%in%stem_hub))

shapiro.test(resgC[gene%in%stem_hub]$GeneScore) #not normal
wilcox.test(resgC[gene%in%stem_hub]$GeneScore,resgC[!(gene%in%stem_hub)]$GeneScore)#p=1.519e-14

nutrient_paths<-c("Longevity","mTOR","PI3K","AMPK sig","Insulin sig","FoxO")
nutrient_ids<-sapply(nutrient_paths,function(path_word)unique(resK_gsea_C[str_detect(Description,path_word)]$ID))
nutrient_ids#no pi3k
nutrient_ids$PI3K<-"hsa04151"
nutrient_genes<-lapply(nutrient_ids,getGenesKEGGPathw)

nutrient_hub<-lapply(nutrient_genes, function(genes){
  return(genes[sapply(genes, function(gene)sum(sapply(nutrient_genes, function(path)gene%in%path))>=length(nutrient_ids)/4)])
})
nutrient_hub<-as.vector(unique(unlist(nutrient_hub)))
length(nutrient_hub)#146
ggplot(resgC)+geom_boxplot(aes(x=compa,y=GeneScore,fill=gene%in%nutrient_hub))

shapiro.test(resgC[gene%in%nutrient_hub]$GeneScore) #not normal
wilcox.test(resgC[gene%in%nutrient_hub]$GeneScore,resgC[!(gene%in%nutrient_hub)]$GeneScore)#p=2.96e-06

prolif_paths<-c("Breast cance","Gastric","Chronic mye","Acute myeloid","Endometrial ca","Colorectal can","Hepatocellular","Melanoma")
prolif_ids<-sapply(prolif_paths,function(path_word)unique(resK_gsea_C[str_detect(Description,path_word)]$ID))
prolif_ids#ok
prolif_genes<-lapply(prolif_ids,getGenesKEGGPathw)

prolif_hub<-lapply(prolif_genes, function(genes){
  return(genes[sapply(genes, function(gene)sum(sapply(prolif_genes, function(path)gene%in%path))>=length(prolif_ids)/4)])
})
prolif_hub<-as.vector(unique(unlist(prolif_hub)))
length(prolif_hub)#158
ggplot(resgC)+geom_boxplot(aes(x=compa,y=GeneScore,fill=gene%in%prolif_hub))

shapiro.test(resgC[gene%in%prolif_hub]$GeneScore) #not normal
wilcox.test(resgC[gene%in%prolif_hub]$GeneScore,resgC[!(gene%in%prolif_hub)]$GeneScore)#p=2.96e-06

#hormon  :
hormon_paths<-c("Estrogen signaling pathway","GnRH secretion","Morphine addiction","Adrenergic signaling in cardiomyocytes","Oxytocin signaling pathway",
                 "Prolactin signaling pathway","Aldosterone synthesis and secretion","Cholinergic synapse","Parathyroid hormone synthesis, secretion and action",
                 "Endocrine resistance","Growth hormone synthesis, secretion and action","Insulin secretion","Cortisol synthesis and secretion")

hormon_ids<-sapply(hormon_paths,function(path_word)unique(resK_gsea_C[str_detect(Description,path_word)]$ID))
hormon_ids#no
hormon_ids<-c("hsa04935","hsa04915","hsa04911","hsa04927","hsa04918","hsa04928","hsa04925")
hormon_genes<-lapply(hormon_ids,getGenesKEGGPathw)

hormon_hub<-lapply(hormon_genes, function(genes){
  return(genes[sapply(genes, function(gene)sum(sapply(hormon_genes, function(path)gene%in%path))>=length(hormon_ids)/4)])
})
hormon_hub<-as.vector(unique(unlist(hormon_hub)))
length(hormon_hub)#111
ggplot(resgC)+geom_boxplot(aes(x=compa,y=GeneScore,fill=gene%in%hormon_hub))
shapiro.test(resgC[gene%in%hormon_hub]$GeneScore) #not normal
wilcox.test(resgC[gene%in%hormon_hub]$GeneScore,resgC[!(gene%in%hormon_hub)]$GeneScore)#p=2.96e-06


ul<-function(x)as.vector(unique(unlist(x)))
modules<-list(prolif=ul(prolif_genes),stem=ul(stem_genes),nutrient=ul(nutrient_genes),hormon=ul(hormon_genes))
saveRDS(modules,"analyses/paper/modules_stem_nutrient_prolif.rds")

hub<-list(prolif=prolif_hub,stem=stem_hub,nutrient=nutrient_hub,hormon=hormon_hub)
saveRDS(hub,"analyses/paper/hub_stem_nutrient_prolif.rds")

# 3) Validation
#   a.PseudoBulk 
#   plot DEG~GeneScore et genes in Metalongevity
cbps<-readRDS("../singlecell/analyses/02-hematopo_datasets_integration/all_cbps/all_cbps.rds")
#with wilcoxon
#C vs LGA

Idents(cbps)<-"diff_state"
DimPlot(cbps)

#PSEUDO BULK DEG analysis
# Prep the bulk_datas
##[start prep data]

unique(cbps@meta.data$orig.ident)
cbps#56343 features across 19283 samples 

#filter for ambigous sample, and CBP1 sample, and get ctrl sample
counts <- subset(cbps,ambigous==FALSE&orig.ident!="cd34_hto1_0C1I1L"&group=="ctrl"&diff_state=="progen")@assays$RNA@counts 
dim(counts) #33600  8869
metadata <- subset(cbps,ambigous==FALSE&orig.ident!="cd34_hto1_0C1I1L"&group=="ctrl"&diff_state=="progen")@meta.data

table(unique(data.table(metadata),by="sample")$group_sex) #3F, 4M

## Remove lowly expressed genes which have less than 10 cells with any counts
counts <- counts[rowSums(counts > 0) >= 10, ]

dim(counts) #15563  8869




# Create a data frame with the sample IDs, cluster IDs and condition

mtsc<-unique(data.table(metadata),by=c("sample"))
mtsc <- mtsc[,.(sample,sex,orig.ident)] 


# Assign the rownames of the metadata to be the sample IDs
sample_metadata<-data.frame(mtsc)
rownames(sample_metadata) <- sample_metadata$sample
head(sample_metadata)

# Aggregate across cluster-sample groups
library(Matrix)
library(Matrix.utils)
sample_counts <- t(aggregate.Matrix(t(as.matrix(counts)), 
                       groupings = metadata$sample, fun = "sum") )

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
sample_counts<-sample_counts[,order(colnames(sample_counts))]
sample_metadata<-sample_metadata[order(rownames(sample_metadata)),]
all(rownames(sample_metadata) == colnames(sample_counts))        
#[end prep data]

#I) CTRL vs LGA  : regression group, batch et sex
dds1 <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = sample_metadata, 
                               design = ~ sex+orig.ident)

#Visualisation des datas
#[start visual]
# Transform counts for data visualization
rld <- rlog(dds1, blind=TRUE)

# Plot PCAby deseq2
DESeq2::plotPCA(rld, intgroup = "orig.ident")
DESeq2::plotPCA(rld, intgroup = "sex")
# Plot others PCs
pcs<-prcomp(t(assay(rld)))
pcs<-pcs$x
pcs<-data.frame(pcs)
pcs$sample<-as.factor(rownames(pcs))
pcs$group<-sample_metadata[rownames(pcs),"group"]
pcs$batch<-sample_metadata[rownames(pcs),"orig.ident"]
pcs$sex<-sample_metadata[rownames(pcs),"sex"]
pcs$group_sex<-sample_metadata[rownames(pcs),"group_sex"]
library(ggplot2)
ggplot(pcs)+geom_point(aes(x=PC2,y=PC3,col=sex))

ggplot(pcs)+geom_point(aes(x=PC3,y=PC4,col=sex))
ggplot(pcs)+geom_point(aes(x=PC4,y=PC5,col=sex)) #separation des group !
ggplot(pcs)+geom_point(aes(x=PC4,y=PC5,col=sex))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_cor <- cor(assay(rld))

# Plot heatmap
library(pheatmap)
pheatmap(rld_cor, annotation = sample_metadata[, c("sex","orig.ident"), drop=F])
#[end visual]


# Run DESeq2 differential expression analysis
dds1 <- DESeq(dds1)

# Plot dispersion estimates
plotDispEsts(dds1)

# Output results of Wald test for contrast 
resultsNames(dds1)
# resultsNames(dds)
res1 <- results(dds1, 
                contrast = c("sex","M","F"),
                alpha = 0.05)

res1 <- lfcShrink(dds1, 
                  contrast = c("group","M","F"),
                  res=res1)
?lfcShrink
?results

# Turn the results object into a tibble for use with tidyverse functions
library(tibble)
res_tbl1 <- res1 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl1
res_tbl1<-data.table(res_tbl1)

res_meth1<-merge(res_tbl1,resgC[,.(gene,GeneScore,pval1000perm,locisID,meth.change,pval,nCpG.Gene,nCpGSig.Gene)],by="gene",all.x=T)

res_meth1[padj<0.05]

res_meth1[padj<0.05&GeneScore>50]
fwrite(res_meth1[order(padj)],
       "analyses/paper/ctrlM_vs_ctrlF/2020-10-06_pseudo_bulk_DEseq2_ctrlMVsCtrlF_regr_on_batch_all_genes.csv",
       sep=";")


res_meth1[is.na(padj),padj:=1]
ggplot(res_meth1)+geom_point(aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.05))+ scale_color_manual(values = c("grey","red")) + theme_minimal() +theme(legend.position = "none")
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_volcano_plot_DEG.png")

ggplot(res_meth1[padj<1])+geom_point(aes(x=log2FoldChange,y=GeneScore,col=padj<0.05))+
  scale_x_continuous(limits = c(-2.5,2.5))+scale_color_manual(values = c("grey","red")) +
  theme_minimal() +theme(legend.position = "none")
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_plot_DEG_vs_DMG.png")


#♣with label
ggplot(res_meth1[padj<1],aes(log2FoldChange,-log10(pvalue),col=padj<0.05))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  geom_label_repel(aes(label = ifelse(padj<0.01&abs(log2FoldChange)>2,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  theme_minimal() +
  theme(legend.position = "none")
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_volcano_plot_DEG_anno_padj0.01_FC2.png")


ggplot(res_meth1[padj<1],aes(log2FoldChange,-log10(pvalue),col=padj<0.05))+
  geom_point()+
  scale_x_continuous(limits = c(-5,5))+scale_color_manual(values = c("grey","red"))+
  geom_label_repel(aes(label = ifelse(padj<0.05&abs(log2FoldChange)>1.5,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  theme_minimal() +
  theme(legend.position = "none")
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_volcano_plot_DEG_zoom_anno_padj0.05_FC1.5.png")


# plot FC ~ GeneScore, DEG
library(ggrepel)
ggplot(res_meth1[padj<1&!is.na(GeneScore)],aes(log2FoldChange,GeneScore,col=padj<0.05))+
  geom_point()+
  scale_x_continuous(limits = c(-2.5,2.5))+scale_color_manual(values = c("grey","red"))+
  geom_label_repel(aes(label = ifelse(padj<0.05&abs(GeneScore)>30&abs(log2FoldChange)>0.6,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+theme_classic()
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_plot_DEG_vs_DMG_anno.png")
res_meth1[padj<0.05&abs(GeneScore)>60&abs(log2FoldChange)>0.5]

#plot module
#volcano
p1<-ggplot(res_meth1[padj<1],aes(log2FoldChange,-log10(pvalue),col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$nutrient,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("nutrient")+
  theme_minimal() +
  theme(legend.position = "none")

p2<-ggplot(res_meth1[padj<1],aes(log2FoldChange,-log10(pvalue),col=padj<0.1))+
    geom_point()+
    scale_color_manual(values = c("grey","red"))+
    geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$stem,gene,"")),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50')+
    ggtitle("stem")+
  theme_minimal() +
    theme(legend.position = "none")
p3<-ggplot(res_meth1[padj<1],aes(log2FoldChange,-log10(pvalue),col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$prolif,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("prolif")+
theme_minimal() +
  theme(legend.position = "none")
p4<-ggplot(res_meth1[padj<1],aes(log2FoldChange,-log10(pvalue),col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$hormon,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("hormon")+
theme_minimal() +
  theme(legend.position = "none")
  
(p1+p2)/(p3+p4)
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_2020-10-06_volcano_plots_DEG_zoom_anno_modules.png.png")

#degvsDMG
p1<-ggplot(res_meth1[padj<1],aes(log2FoldChange,GeneScore,col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  scale_x_continuous(limits = c(-5,5))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$nutrient,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("nutrient")+
  theme_minimal() +
  theme(legend.position = "none")

p2<-ggplot(res_meth1[padj<1],aes(log2FoldChange,GeneScore,col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  scale_x_continuous(limits = c(-5,5))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$stem,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("stem")+
  theme_minimal() +
  theme(legend.position = "none")
p3<-ggplot(res_meth1[padj<1],aes(log2FoldChange,GeneScore,col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  scale_x_continuous(limits = c(-5,5))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$prolif,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("prolif")+
  theme_minimal() +
  theme(legend.position = "none")
p4<-ggplot(res_meth1[padj<1],aes(log2FoldChange,GeneScore,col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  scale_x_continuous(limits = c(-5,5))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%hub$hormon,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("hormon")+
  theme_minimal() +
  theme(legend.position = "none")

(p1+p2)/(p3+p4)
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_2020-10-06_plots_DEG_vs_DMG_anno_modules.png")



growthFactor_genes<-tr("9021/8651/8835/8660/3726/1387/2778/4214/775/111/5578/6774/776/5606/208/9564/5331/2767/112/5747/6755/3710/6416/3717/6752/5600/90993/1385/25759/3667/115/53358/5601/3265/2776/3486/2475/5295/6654/5566/2771/2885/10488/113/5579/6776/6777/109/6655/5595/5335/3708/5336/5332/8503/1398/5594/5293/2932/10000/2770/1388/84699/5567/2773/207/5604/5296/2033/399694/9586/107/5582/196883/1432/114",tradEntrezInSymbol = T)
ggplot(res_meth1[padj<1],aes(log2FoldChange,-log10(pvalue),col=padj<0.1))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+
  geom_label_repel(aes(label = ifelse(padj<0.1&gene%in%growthFactor_genes,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  ggtitle("growth_hormon")+
  theme_minimal() +
  theme(legend.position = "none")


#   b.diff male fem en in hsc

cbps_ctrl<-subset(cbps,group=="ctrl")
mtd<-data.table(cbps_ctrl@meta.data)
mtd[,n_cells.sample:=.N,by="sample"]
mtd[,pct.ct:=.N/n_cells.sample,by=c("cell_type","sample")]

mtd[,pct.lin:=.N/n_cells.sample,by=c("new.lineage","sample")]
unique(mtd[,.(sample,sex,cell_type,pct.ct)])
ggplot(unique(mtd[,.(sample,sex,cell_type,pct.ct)]))+
  geom_boxplot(aes(x=sex,y=pct.ct,col=sex))+
  facet_wrap("cell_type",scales = "free_y")+
  scale_y_continuous(labels = scales::percent)
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_distr_cell_type_by_sex.png")
ggplot(unique(mtd[,.(sample,sex,new.lineage,pct.lin)]))+
  geom_boxplot(aes(x=sex,y=pct.lin,col=sex))+
  facet_wrap("new.lineage",scales = "free_y")+
  scale_y_continuous(labels = scales::percent)
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_distr_lineage_by_sex.png")


#II) CTRL vs LGA

# 1) GENESCORE
# 2) PATHWAY
#   2a) Over-repres test
#   2b) GSEA (or similar)
#   2c) longevity effect increase in LGA
#     i) pval pathway increase
#     ii) methylScore increase 


# 
# 3) PseudoBulk
#   plot DEG~GeneScore et genes in Metalongevity

#III) Validation
# 1) Public data
#   1a) epigen clock
#   1b) cd34 m vs f or sc CTRL LGA
# 2) reponse aux stress nutri
