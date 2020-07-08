
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

#3) pathway analysis on regul DEG epigen regulated

lep_genes<-c("PLA2G4A",	"MAPK1",	"CRP",	"EDN1",	"EGR1",	"FOS",	 "GRB2",	"HIF1A",	"HRAS",	"JAK2",	"LEP",	"LEPR",
"POMC",	"MAPK3",	"MAP2K1",	"PTPN1",	"PTPN11",	"RAF1",	"SOS1",	"SOS2",	"STAT3",	"TIMP1",	"TRH", "VEGFA",
"PLA2G4C", 	"SOCS3",	"PIAS3",	"CYCS")

res_meth1[,lep_gene:=gene%in%lep_genes]
res_meth1[lep_gene==TRUE,]
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


