library(stringr)
library(data.table)
library(clusterProfiler)
library(Seurat)
library(DESeq2)
library(limma)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(org.Hs.eg.db)
source("scripts/utils/methyl_utils.R")

set.seed(1234)
options(stringsAsFactors = T)

#I) CTRLF vs CTRLM
dir.create("analyses/paper/ctrlM_vs_ctrlF",recursive = T)
# 1) GENESCORE
res<-fread("analyses/model14_without_iugr/2020-09-24_all_res_with_perm.csv")
#add entrez id
listEnsembl()
ensembl <- useEnsembl(biomart = "ensembl")
searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
searchAttributes(mart = ensembl, pattern = "hgnc") #hgnc_symbol
searchAttributes(mart = ensembl, pattern = "entrez") #entrezgene_id
genes<-getBM(attributes = c('hgnc_symbol', 'entrezgene_id'),
      filters = 'hgnc_symbol',
      values = unique(res$gene), 
      mart = ensembl)
genes<-data.table(genes)[,gene:=hgnc_symbol][,entrez:=entrezgene_id]
genes[,is.true:=entrez==min(entrez),by="gene"]
genes<-genes[is.true==T]
res<-merge(res[,-"entrez"],genes[,.(gene,entrez)],all.x=T,by="gene")
fwrite(res,"analyses/model14_without_iugr/2020-09-24_all_res_with_perm.csv")
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

stem_paths<-c("Wnt","TGF-beta","JAK-STAT","pluripotency","MAPK")
stem_ids<-sapply(stem_paths,function(path_word)unique(resK_gsea_C[str_detect(Description,path_word)]$ID))
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
prolif_ids#no pi3k
prolif_genes<-lapply(prolif_ids,getGenesKEGGPathw)

prolif_hub<-lapply(prolif_genes, function(genes){
  return(genes[sapply(genes, function(gene)sum(sapply(prolif_genes, function(path)gene%in%path))>=length(prolif_ids)/4)])
})
prolif_hub<-as.vector(unique(unlist(prolif_hub)))
length(prolif_hub)#158
ggplot(resgC)+geom_boxplot(aes(x=compa,y=GeneScore,fill=gene%in%prolif_hub))

shapiro.test(resgC[gene%in%prolif_hub]$GeneScore) #not normal
wilcox.test(resgC[gene%in%prolif_hub]$GeneScore,resgC[!(gene%in%prolif_hub)]$GeneScore)#p=2.96e-06


# 3) PseudoBulk 
#   plot DEG~GeneScore et genes in Metalongevity
cbps<-readRDS("../singlecell/analyses/02-hematopo_datasets_integration/all_cbps/all_cbps.rds")
#with wilcoxon
#C vs LGA
Idents(cbps)<-"cell_type"
levels(cbps)
# [1] "1-LMPP"           "2-HSC-AVP"        "3-HSC-SOCS3"      "4-MPP-EMP"       
# [5] "5-HSC-EGR1"       "6-HSC-AREG"       "7-LyP"            "8-EMP"           
# [9] "9-proT"           "10-HSC-CD164"     "11-GMP"           "12-ProB-2"       
# [13] "13-MEP"           "14-div_cells"     "15-B cell"        "16-pDC-Cycle"    
# [17] "17-Ba/Eo/Mas"     "18-LT-HSC"        "19-HSC/MPP-HSPA5" "20-T cell"       
# [21] "21-HSC/MPP-IRF1"  "22-proB-SPIB"     "23-Mas"           "24-Mo"           
# [25] "25-NK" 


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
mtsc<-data.frame(mtsc)
rownames(mtsc) <- mtsc$sample
head(mtsc)

# Aggregate across cluster-sample groups
library(Matrix)
library(Matrix.utils)
pb <- t(aggregate.Matrix(t(as.matrix(counts)), 
                       groupings = metadata, fun = "sum") )
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
cluster_counts
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
cluster_counts<-cluster_counts[,order(colnames(cluster_counts))]
cluster_metadata<-cluster_metadata[order(rownames(cluster_metadata)),]
all(rownames(cluster_metadata) == colnames(cluster_counts))        
#[end prep data]

#I) CTRL vs LGA  : regression group, batch et sex
dds1 <- DESeqDataSetFromMatrix(cluster_counts, 
                               colData = cluster_metadata, 
                               design = ~ group+sex+orig.ident)

#Visualisation des datas
#[start visual]
# Transform counts for data visualization
rld <- rlog(dds1, blind=TRUE)

# Plot PCAby deseq2
DESeq2::plotPCA(rld, intgroup = "orig.ident")
DESeq2::plotPCA(rld, intgroup = "group")
# Plot others PCs
pcs<-prcomp(t(assay(rld)))
pcs<-pcs$x
pcs<-data.frame(pcs)
pcs$sample<-as.factor(rownames(pcs))
pcs$group<-cluster_metadata[rownames(pcs),"group"]
pcs$batch<-cluster_metadata[rownames(pcs),"orig.ident"]
pcs$sex<-cluster_metadata[rownames(pcs),"sex"]
pcs$group_sex<-cluster_metadata[rownames(pcs),"group_sex"]
library(ggplot2)
ggplot(pcs)+geom_point(aes(x=PC2,y=PC3,col=group))

ggplot(pcs)+geom_point(aes(x=PC3,y=PC4,col=group))
ggplot(pcs)+geom_point(aes(x=PC4,y=PC5,col=group)) #separation des group !
ggplot(pcs)+geom_point(aes(x=PC4,y=PC5,col=batch))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_cor <- cor(assay(rld))

# Plot heatmap
library(pheatmap)
pheatmap(rld_cor, annotation = cluster_metadata[, c("group","sex","orig.ident"), drop=F])
#[end visual]


# Run DESeq2 differential expression analysis
dds1 <- DESeq(dds1)

# Plot dispersion estimates
plotDispEsts(dds1)

# Output results of Wald test for contrast 
resultsNames(dds1)
# resultsNames(dds)
res1 <- results(dds1, 
                contrast = c("group","lga","ctrl"),
                alpha = 0.05)

res1 <- lfcShrink(dds1, 
                  contrast = c("group","lga","ctrl"),
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
mLF<-unique(fread("../Sexual_Dimorphism/JANUARY2020/Alexandre_Methyl/analyses/withoutIUGR/2020-07-30_resCF_LF.csv")[order(-GeneScore,pval)],by="gene")
res_meth1<-merge(res_tbl1,mLF[,.(gene,GeneScore,pval1000perm,locisID,meth.change,pval,nCpG.Gene,nCpGSig.Gene)],by="gene",all.x=T)
res_meth1[padj<0.05]
res_meth1[padj<0.05&GeneScore>50]
fwrite(res_meth1,
       "analyses/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv",
       sep=";")

padj_cutoff <- 0.05
# Subset the significant results
sig_res <- dplyr::filter(res_meth1, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res
fwrite(sig_res,"analyses/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_padj0.05_DEG.csv",sep=";")

ggplot(res_meth1)+geom_point(aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.05))+ scale_color_manual(values = c("grey","red")) + theme_minimal() +theme(legend.position = "none")

ggplot(res_meth1)+geom_point(aes(x=log2FoldChange,y=GeneScore,col=padj<0.05))+ scale_color_manual(values = c("grey","red")) + theme_minimal() +theme(legend.position = "none")

res_meth1[padj<0.05&abs(GeneScore)>60&abs(log2FoldChange)>0.5]




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
