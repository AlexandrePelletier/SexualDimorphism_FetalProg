
#2020-07-08 explore GeneSig
#visual res of LGAvsCTRL


library(data.table)
res_meth1<-fread("../../../Alexandre_SC/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")

#1)see diff d'expr

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)

library(tidyr)
gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- inner_join(ei[, c("sample", "group" )], gathered_top20_sig, by = c("sample" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = group), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


#2) show enrichment Genscore in DEG
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

ggplot(res_meth1[padj<0.2],aes(log2FoldChange,GeneScore))+
  geom_point()+
  geom_label_repel(aes(label = ifelse(GeneScore>60&abs(log2FoldChange)>0.8,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+theme_classic()

#un gene a un methchange neg mais genescore pos : 
res_meth1[meth.change<0&padj<0.2&GeneScore>0] #KLF4

#KLF4
resF<-fread("analyses/withoutIUGR/2020-07-03_resCF_LF.csv")
resF[gene=="KLF4"&pval<0.01]

resM<-fread("analyses/withoutIUGR/2020-07-03_resCM_LM.csv")
resM[gene=="KLF4"&tss_dist>0&tss_dist<1000][order(tss_dist)] #hypermeth sex spé



#integr DEG & single cell 
scDEG<-fread("../../../Alexandre_SC/analyses/04-DEG_in_LGA/htos_cbp123_LF.CMFLM_DESeq2_cells_filtered.csv")


cols<-c("p_val","avg_logFC","pct.1","pct.2","p_val_adj")
scDEG[,(paste0(cols,"_sc")):=.SD,.SD=cols]

res_meth2<-merge(res_meth1,scDEG[,.(gene,p_val_sc,avg_logFC_sc,pct.1_sc,pct.2_sc,p_val_adj_sc)])
ggplot(res_meth2)+geom_point(aes(-log10(padj),-log10(p_val_adj_sc)))
ggplot(res_meth2)+geom_point(aes(log2FoldChange,avg_logFC_sc))
ggplot(res_meth2[DEG==TRUE])+geom_point(aes(padj,p_val_adj_sc))
ggplot(res_meth2[DEG==TRUE])+geom_point(aes(log2FoldChange,avg_logFC_sc))


#3) pathway analysis on regul DEG epigen regulated
lep_genes<-c("PLA2G4A",	"MAPK1",	"CRP",	"EDN1",	"EGR1",	"FOS",	 "GRB2",	"HIF1A",	"HRAS",	"JAK2",	"LEP",	"LEPR",
"POMC",	"MAPK3",	"MAP2K1",	"PTPN1",	"PTPN11",	"RAF1",	"SOS1",	"SOS2",	"STAT3",	"TIMP1",	"TRH", "VEGFA",
"PLA2G4C", 	"SOCS3",	"PIAS3",	"CYCS")

res_meth1[,lep_gene:=gene%in%lep_genes]
res_meth1[lep_gene==TRUE,]

res_meth1[lep_gene==TRUE&pvalue<0.05]


ggplot(res_meth1[!is.na(DEG)],aes(x=lep_gene,y=GeneScore))+
  geom_boxplot()+
  geom_point(aes(col=ifelse(lep_gene==TRUE,"red",element_blank())))+
  geom_label_repel(aes(label = ifelse(lep_gene==TRUE&GeneScore>60,gene,"")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  theme_classic()

ggplot(res_meth1[lep_gene==TRUE],aes(log2FoldChange,GeneScore))+
  geom_point()
median(res_meth1$GeneScore,na.rm = T)

#DISCOVER PATHWAY / BIOLOGICAL PROCESS AFFECTED IN FEMALE :
#based on meth F only, what are the pathways/BP/msigdb/gwas ?
resF<-fread("analyses/withoutIUGR/2020-07-03_resCF_LF.csv")
resM<-fread("analyses/withoutIUGR/2020-07-03_resCM_LM.csv")

library(patchwork)
p1<-ggplot(unique(resF,by="gene"))+geom_point(aes(-log10(pval1000perm),GeneScore))+ggtitle("Female")
p2<-ggplot(unique(resM,by="gene"))+geom_point(aes(-log10(pval1000perm),GeneScore))+ggtitle("Male")
p1+p2
plot(density(unique(resF,by="gene")$GeneScore))
genesF<-unique(resF[GeneScore>100&pval1000perm<0.01]$gene)
length(genesF)#1700
genesM<-unique(resM[GeneScore>100&pval1000perm<0.01]$gene)
length(genesM)#149

library(clusterProfiler)
library(org.Hs.eg.db)
resF_KEGG <- enrichKEGG(gene         = bitr(genesF,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05)
nrow(data.frame(resF_KEGG)) #55
dotplot(resF_KEGG,showCategory=55)
emapplot(resF_KEGG,showCategory=55)

resM_KEGG <- enrichKEGG(gene         = bitr(genesM,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)
nrow(data.frame(resM_KEGG)) #16
dotplot(resF_KEGG,showCategory=16)
emapplot(resF_KEGG,showCategory=16)

enrichplot::upsetplot(resF_KEGG)
enrichplot::upsetplot(resM_KEGG)

#then, on all this pathways, how many are female spe ?
pathKEGG_Fo<-setdiff(data.frame(resF_KEGG)$Description,data.frame(resM_KEGG)$Description)
pathKEGG_Fo
pathKEGG_Fo[]

#on this pathway could we see gene expression dysregulation in LGAF ?

#=> pathways pertubed in LGAF, what are the hypothesis for this pertubance, and what are the possible biological consequences?

#how validate this hypothesis ?

