#What is the epigen impact of stress exposure to overnutrition?
library(data.table)
library(ggplot2)
source("scripts/utils/methyl_utils.R")
library(limma)
#DATA EXPLO
methyl<-fread("datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
methyl

meta<-fread("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")
meta
samples<-colnames(methyl)sample%in%meta[Group_name%in%c("C","L")]$sample
methyl<-methyl[,.SD,.SDcols=c("locisID",samples)]
dim(methyl) #70 samples
meta<-meta[sample%in%colnames(methyl)]
nrow(meta) #70
table(meta$Group_Sex)
# C_F C_M L_F L_M 
# 16  18  20  16 

#QC
#see "main"

#PCA
pca<-prcomp(t(data.frame(methyl,row.names = 1)))
pct.pcs<-pctPC(pca)
barplot(pct.pcs)
pct.pcs>1 #○to pc47
#which covar is important ? correl



var_num<-c("SeqDepth","Library_Complexity","PI","HC..cm.",
           "Length..cm.","Weight..g.","GA..wk.",
           "Weight.at.term..lbs.","Mat.Age")
var_fac<-c("sequencing","DNA.extraction","date","batch",
           "Gender","Labor","Smoking",
           "Etoh","Drugs","GDM","Preterm","Ethnicity","latino","Group_name","Group_Sex")

plotCovarPCs(pca,1:47,meta,var_num,var_fac)

pcs<-pca$x
pcs<-data.table(data.frame(pcs),keep.rownames = "sample")
pcs

ggplot(meta)+geom_point(aes(x=PC3,y=PC42,col=GDM))

correl(meta$SeqDepth,meta$Group) #no
correl(meta$SeqDepth,meta$Sex) #no
correl(meta$SeqDepth,meta$Group_Sex)#no

correl(meta$Library_Complexity,meta$Group_name) #no
correl(meta$Library_Complexity,meta$Gender) #no
correl(meta$Library_Complexity,meta$Group_Sex)#no

correl(as.factor(meta$batch),as.factor(meta$Group_name))#no
correl(as.factor(meta$batch),as.factor(meta$Gender))#no
correl(as.factor(meta$batch),as.factor(meta$Group_Sex)) #no

correl(as.factor(meta$Drugs),as.factor(meta$Group_name))#no
correl(as.factor(meta$Drugs),as.factor(meta$Gender))#no
correl(as.factor(meta$Drugs),as.factor(meta$Group_Sex)) #no

correl(as.factor(meta$GDM),as.factor(meta$Group_name))#no
correl(as.factor(meta$GDM),as.factor(meta$Gender))#no
correl(as.factor(meta$GDM),as.factor(meta$Group_Sex)) #no

correl(as.factor(meta$latino),as.factor(meta$Group_name))#no
correl(as.factor(meta$latino),as.factor(meta$Gender))#no
correl(as.factor(meta$latino),as.factor(meta$Group_Sex)) #no

correl(as.factor(meta$Preterm),as.factor(meta$Group_name))#no
correl(as.factor(meta$Preterm),as.factor(meta$Gender))#no
correl(as.factor(meta$Preterm),as.factor(meta$Group_Sex)) #no

correl(as.factor(meta$sequencing),as.factor(meta$Group_name))#no
correl(as.factor(meta$sequencing),as.factor(meta$Gender))#no
correl(as.factor(meta$sequencing),as.factor(meta$Group_Sex)) #no

correl(meta$Group_Complexity_Fac,as.factor(meta$Group_Sex)) #no

meta[,(var_fac):=lapply(.SD, as.factor),.SDcols=var_fac]

meta[,Group_Complexity_Fac,as.factor(Group_Complexity_Fac)]

summary(lm(PC1~Library_Complexity:Group_name+Group_name+SeqDepth+latino,data = meta)) #pour PC1

#PC2, unknown variation
summary(lm(PC3~Group_name:Gender+GDM,data = meta)) #PC3, group:sex spe
meta[,Group_Complexity_Fac:=factor(Group_Complexity_Fac,levels = c("1", "2", "3","4"))]
head(meta$Group_Complexity_Fac)

ggplot(meta)+geom_boxplot(aes(x=latino,y=pct0ApresF))

correl(meta$Mat.Age,meta$pct0ApresF) #sig donc on inclut aussi


#MODEL => DMC between Ctrl and LGA, and within and between sex
varToModel<-c("Group_Sex","Mat.Age","latino","Group_Complexity_Fac")

metaF<-meta[rowSums(is.na(meta[,.SD,.SDcols=varToModel]))==0]
formule<- ~0 + Group_Sex  + latino + Mat.Age + Group_Complexity_Fac +batch

design<-model.matrix(formule,data = metaF)
fit <- lmFit(data.frame(methyl,row.names = "locisID"), design)

cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                             F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                             CF.CM="Group_SexC_F-Group_SexC_M",
                             CF.LF="Group_SexC_F-Group_SexL_F",
                             CF.LM="Group_SexC_F-Group_SexL_M",
                             CM.LM="Group_SexC_M-Group_SexL_M",
                             CM.LF="Group_SexC_M-Group_SexL_F",
                             LM.LF="Group_SexL_M-Group_SexL_F",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
results <- decideTests(fit2)
sum(abs(results)) #8
colSums(abs(results))
locisSig<-rownames(fit2$p.value)[apply(fit2$p.value<0.001,1,any)]
length(locisSig) #11561
compas<-colnames(fit2$p.value)

#epigenomics Sex dim : 
cpgs_genes<-fread("ref/2020-06-29_All_CpG-Gene_links.csv")
res_list<-list()
res<-topTable(fit2,coef = compa,n = Inf)
for(compa in compas){
  print(compa)
  res<-topTable(fit2,coef = compa,n = Inf)
  fwrite(res,paste0("analyses/model14_without_iugr/2020-09-09_res_limma_",compa,".csv"),sep=";")
  res<-CalcGeneScore(res,cpgs_genes)
  fwrite(res,paste0("analyses/model14_without_iugr/2020-09-09_res_with_genescore_",compa,".csv"),sep=";")
  res_list[[compa]]<-res
  print(sum(res$pval<0.001))
  
  
}


r<-topTable(fit2,coef = "C.L",n = Inf)
ggplot(r)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>20))+scale_color_manual(values = c("grey","red"))+
  ggtitle(compa,
          subtitle = paste(sum(r$P.Value<0.001),"DMC with pval<0.001,",sum(r$P.Value<0.001&r$logFC>20),"hyperM with logFC>20 and",sum(r$P.Value<0.001&r$logFC<(-20)),"hypoM with logFC<-20"))

for(compa in compas){
  r<-topTable(fit2,coef = compa,n = Inf)
  p<-ggplot(r)+
    geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>20))+
    scale_color_manual(values = c("grey","red"))+
    ggtitle(compa,
            subtitle = paste(sum(r$P.Value<0.001),
                             "DMC with pval<0.001,",
                             sum(r$P.Value<0.001&r$logFC>20),
                             "hyperM with logFC>20 and",
                             sum(r$P.Value<0.001&r$logFC<(-20)),
                             "hypoM with logFC<-20"))
  
  ggsave(paste0("analyses/model14_without_iugr/2020-09-09_volcano_plot_",compa,".png"),p)
  
}

lapply(names(res_list),function(comp)res_list[[comp]][,compa:=..comp])
res<-Reduce(function(x, y) merge(x, y, all = T),res_list)

#2020-09-11 venn diagr
library(data.table)
library(VennDiagram)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
source("scripts/utils/FonctionsUtiles.R")
#all CF
res<-list(
"cfcm"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CF.CM.csv"),
"cflm"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CF.LM.csv"),
"cflf"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CF.LF.csv"),
"cmlf"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CM.LF.csv"),
"cmlm"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_CM.LM.csv"),
"lmlf"=fread("analyses/model14_without_iugr/2020-09-09_res_with_genescore_LM.LF.csv")
)
lapply(names(res),function(comp)res[[comp]][,compa:=..comp])
res<-Reduce(function(x, y) merge(x, y, all = T),res)
#cpg with hi FC and pval
venn.diagram(list("CF.CM"=unique(res[compa=="cfcm"&pval<0.001&abs(meth.change)>20]$locisID),
                  "CF.LM"=unique(res[compa=="cflm"&pval<0.001&abs(meth.change)>20]$locisID),
                  "CF.LF"=unique(res[compa=="cflf"&pval<0.001&abs(meth.change)>20]$locisID)),
             filename = "analyses/model14_without_iugr/venn/venn_locis_CFvsOthers.tiff"
)


#genes with hi FC and pval
genes.list<-list("CF.CM"=unique(res[compa=="cfcm"&pval<0.001&abs(meth.change)>20]$gene),
                 "CF.LM"=unique(res[compa=="cflm"&pval<0.001&abs(meth.change)>20]$gene),
                 "CF.LF"=unique(res[compa=="cflf"&pval<0.001&abs(meth.change)>20]$gene))
venn.diagram(genes.list,
             filename = "analyses/model14_without_iugr/venn/venn_genes_CFvsOthers.tiff"
)   

for(comparison in names(genes.list) ){
  print(comparison)
  genes.list[[comparison]]<-bitr(genes.list[[comparison]],
                                 fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db)$ENTREZID #put in ENTREZID because KEGG Pathways is in ENTREZID in this package
}
res.compa<-compareCluster(genes.list,
                          fun = "enrichKEGG",
                          organism="hsa",
                          pvalueCutoff=0.05) 
emapplot(res.compa,showCategory = nrow(as.data.frame(res.compa)))

#ggenes with hi Genescore

p<-ggplot(unique(res,by=c("gene","compa")))+
  geom_jitter(aes(compa,GeneScore,col=abs(GeneScore)>60))+
  scale_color_manual(values = c("grey","red"))

p+geom_boxplot(aes(compa,GeneScore),alpha=0.5,outlier.shape = NA)
  
genes.list<-list("CF.LM"=unique(res[compa=="cflm"&abs(GeneScore)>60]$gene),
                 "CF.LF"=unique(res[compa=="cflf"&abs(GeneScore)>60]$gene),
                 "CF.CM"=unique(res[compa=="cfcm"&abs(GeneScore)>60]$gene))
venn.diagram(genes.list,
             filename = "analyses/model14_without_iugr/venn/venn_genescore_CFvsOthers.tiff") 
for(comparison in names(genes.list) ){
  print(comparison)
  genes.list[[comparison]]<-bitr(genes.list[[comparison]],
                                 fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db)$ENTREZID #put in ENTREZID because KEGG Pathways is in ENTREZID in this package
}
res.compa<-compareCluster(genes.list,
                          fun = "enrichKEGG",
                          organism="hsa",
                          pvalueCutoff=0.05) 
emapplot(res.compa,showCategory = nrow(as.data.frame(res.compa)),pie="count")

#pval/FC vs genescore
# in cfcm
genes.list<-list("genescoreCF.CM"=unique(res[compa=="cfcm"&abs(GeneScore)>60]$gene),
                 "genesCF.CM"=unique(res[compa=="cfcm"&pval<0.001&abs(meth.change)>20]$gene))
venn.diagram(genes.list,
             filename = "analyses/model14_without_iugr/venn/venn_genescore_vs_genes_cfcm.tiff" )


genes.list
for(comparison in names(genes.list) ){
  print(comparison)
  genes.list[[comparison]]<-bitr(genes.list[[comparison]],
                                          fromType = "SYMBOL",
                                          toType = "ENTREZID",
                                          OrgDb = org.Hs.eg.db)$ENTREZID #put in ENTREZID because KEGG Pathways is in ENTREZID in this package
}
genes.list

res.compa<-compareCluster(genes.list,
                          fun = "enrichKEGG",
                          organism="hsa",
                          pvalueCutoff=0.05) 
emapplot(res.compa,showCategory = nrow(as.data.frame(res.compa)))

resK<-enrichKEGG(setdiff(genes.list$genescoreCF.CM,genes.list$genesCF.CM),organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(resK)#no

# in cflf
genes.list<-list("genescoreCF.LF"=unique(res[compa=="cflf"&abs(GeneScore)>60]$gene),
                 "genesCF.LF"=unique(res[compa=="cflf"&pval<0.001&abs(meth.change)>20]$gene))
venn.diagram(genes.list,
             filename = "analyses/model14_without_iugr/venn/venn_genescore_vs_genes_cflf.tiff" )


for(comparison in names(genes.list) ){
  print(comparison)
  genes.list[[comparison]]<-bitr(genes.list[[comparison]],
                                 fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db)$ENTREZID #put in ENTREZID because KEGG Pathways is in ENTREZID in this package
}
res.compa<-compareCluster(genes.list,
                          fun = "enrichKEGG",
                          organism="hsa",
                          pvalueCutoff=0.05) 
emapplot(res.compa,showCategory = nrow(as.data.frame(res.compa)),pie="count")

resK<-enrichKEGG(setdiff(genes.list$genescoreCF.LF,genes.list$genesCF.LF),organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(resK)#no

# in cflf=m
genes.list<-list("genescoreCF.LM"=unique(res[compa=="cflm"&abs(GeneScore)>60]$gene),
                 "genesCF.LM"=unique(res[compa=="cflm"&pval<0.001&abs(meth.change)>20]$gene))
venn.diagram(genes.list,
             filename = "analyses/model14_without_iugr/venn/venn_genescore_vs_genes_cflm.tiff" )


for(comparison in names(genes.list) ){
  print(comparison)
  genes.list[[comparison]]<-bitr(genes.list[[comparison]],
                                 fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db)$ENTREZID #put in ENTREZID because KEGG Pathways is in ENTREZID in this package
}
res.compa<-compareCluster(genes.list,
                          fun = "enrichKEGG",
                          organism="hsa",
                          pvalueCutoff=0.05) 
emapplot(res.compa,showCategory = nrow(as.data.frame(res.compa)),pie="count")

resK<-enrichKEGG(setdiff(genes.list$genescoreCF.LF,genes.list$genesCF.LF),organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(resK)#no


#Common Stress response
#=> Intersection with intersection(CM.LM CM.LF) & Setdiff with LM.LF

p<-ggplot(unique(res,by=c("gene","compa")))+
  geom_jitter(aes(compa,GeneScore,col=abs(GeneScore)>60))+
  scale_color_manual(values = c("grey","red"))

p+geom_boxplot(aes(compa,GeneScore),alpha=0.5,outlier.shape = NA)

genes.list<-list("CM.LM"=unique(res[compa=="cmlm"&abs(GeneScore)>60]$gene),
                 "CM.LF"=unique(res[compa=="cmlf"&abs(GeneScore)>60]$gene),
                 "LM.LF"=unique(res[compa=="lmlf"&abs(GeneScore)>60]$gene))
venn.diagram(genes.list,
             filename = "analyses/model14_without_iugr/venn/venn_genescore_cmvslm_cmvslf_lgamvslgaF.tiff") 

genes<-intersect(unique(res[compa=="cmlm"&abs(GeneScore)>60]$gene),unique(res[compa=="cmlf"&abs(GeneScore)>60]$gene))
# Setdiff with LM.LF
genes<-setdiff(genes,unique(res[compa=="lmlf"&abs(GeneScore)>60]$gene))

resK<-enrichKEGG(bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(resK)#Hippo signaling pathway
tr(as.data.frame(resK)$geneID,tradEntrezInSymbol = T) #"PRKCZ"  "PPP1CA" "CCND1"  "CCND2"  "WNT1"   "SMAD3"  "AXIN1"  "ACTG1" "CSNK1D"
#intersection with cf
genes<-intersect(genes,setdiff(intersect(unique(res[compa=="cflm"&abs(GeneScore)>60]$gene),
                                  unique(res[compa=="cflf"&abs(GeneScore)>60]$gene)),
                        unique(res[compa=="cfcm"&abs(GeneScore)>60]$gene)))

length(genes) #173
resK<-enrichKEGG(bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(resK)#Hippo signaling pathway
tr(as.data.frame(resK)$geneID,tradEntrezInSymbol = T)



#GSEA [a simplifier] ? => functional impact different between female and male
library(clusterProfiler)
library(org.Hs.eg.db)
resg<-unique(res[order(gene,compa,pval)],by=c("gene","compa"))
genes<-bitr(resg$gene,"SYMBOL","ENTREZID",org.Hs.eg.db)
genes<-data.table(gene=genes$SYMBOL,entrez_id=genes$ENTREZID)
resg<-merge(resg,genes,by="gene")

#diff methylated pathway
geneList<-abs(resg[compa==compas[1]]$GeneScore)
names(geneList)<-resg[compa==compas[1]]$entrez_id
geneList<-sort(geneList,decreasing = T)
resK<- gseKEGG(geneList     = rank(geneList),
               organism     = 'hsa', 
               minGSSize    = 50,
               pvalueCutoff = 0.05,
               verbose = FALSE)
resK<-data.table(as.data.frame(resK))
resK[,compa:=compas[1]]

for(comp in compas[2:length(compas)]){
  geneList<-abs(resg[compa==comp]$GeneScore)
  names(geneList)<-resg[compa==comp]$entrez_id
  geneList<-sort(geneList,decreasing = T)
  resk<- gseKEGG(geneList     = rank(geneList),
                           organism     = 'hsa', 
                           minGSSize    = 50,
                           pvalueCutoff = 0.05,
                           verbose = FALSE)
  resk<-data.table(as.data.frame(resk))
  resk[,compa:=..comp]
  resK<-merge(resK,resk,all=T)
}

head(resK[order(p.adjust)],50)
resK[,score_type:="DMG"]

#hyper methylated pathway
geneList<-resg[compa==compas[1]]$GeneScore
names(geneList)<-resg[compa==compas[1]]$entrez_id
geneList<-sort(geneList,decreasing = T)
resK2<- gseKEGG(geneList     = rank(geneList),
               organism     = 'hsa', 
               minGSSize    = 50,
               pvalueCutoff = 0.05,
               verbose = FALSE)
resK2<-data.table(as.data.frame(resK2))
resK2[,compa:=compas[1]]

for(comp in compas[2:length(compas)]){
  geneList<-abs(resg[compa==comp]$GeneScore)
  names(geneList)<-resg[compa==comp]$entrez_id
  geneList<-sort(geneList,decreasing = T)
  resk<- gseKEGG(geneList     = rank(geneList),
                 organism     = 'hsa', 
                 minGSSize    = 50,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)
  resk<-data.table(as.data.frame(resk))
  resk[,compa:=..comp]
  resK2<-merge(resK2,resk,all=T)
}

head(resK2[order(p.adjust)],50)
resK2[,score_type:="hyperM"]
resK<-merge(resK,resK2,all=T,by=colnames(resK))

#hippo methylated pathway
geneList<-(-resg[compa==compas[1]]$GeneScore)
names(geneList)<-resg[compa==compas[1]]$entrez_id
geneList<-sort(geneList,decreasing = T)
resK2<- gseKEGG(geneList     = rank(geneList),
                organism     = 'hsa', 
                minGSSize    = 50,
                pvalueCutoff = 0.05,
                verbose = FALSE)
resK2<-data.table(as.data.frame(resK2))
resK2[,compa:=compas[1]]

for(comp in compas[2:length(compas)]){
  geneList<-abs(resg[compa==comp]$GeneScore)
  names(geneList)<-resg[compa==comp]$entrez_id
  geneList<-sort(geneList,decreasing = T)
  resk<- gseKEGG(geneList     = rank(geneList),
                 organism     = 'hsa', 
                 minGSSize    = 50,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)
  resk<-data.table(as.data.frame(resk))
  resk[,compa:=..comp]
  resK2<-merge(resK2,resk,all=T)
}

head(resK2[order(p.adjust)],50)
resK2[,score_type:="hypoM"]
resK<-merge(resK,resK2,all=T,by=colnames(resK))

plot(density(resg$GeneScore))
abline(v=40)
abline(v=-40)

#OverRepr
{
  candidat_genes.list<-lapply(compas, function())
  
  resFM_KEGG <- compareCluster(gene         =candidat_genes.list ,
                               fun = "enrichKEGG",
                               organism="hsa",
                               pvalueCutoff=0.1) 
}

#0) Pathway sig ? (permut)


#1) is there pathway spe ?
{
  unique(resK$compa)
  #for LGA F : pathway in all compa with LF, but not in others
  resK[(compa=="CF.LF"&score_type%in%c("hyperM","DMG"))|
         (compa=="CM.LF"&score_type%in%c("hyperM","DMG"))|
         (compa=="LM.LF"&score_type%in%c("hyperM","DMG"))|
         (compa=="C.L"&score_type%in%c("hyperM","DMG")),
       in.all.LGAF:=.N==8,"Description"]
  unique(resK[in.all.LGAF==T]$Description)
  resK[,exclu.LGAF:=all(in.all.LGAF==TRUE),"Description"]
  resK[exclu.LGAF==T]
  
  #for LGA M
  resK[(compa=="CF.LM"&score_type%in%c("hyperM","DMG"))|
         (compa=="CM.LM"&score_type%in%c("hyperM","DMG"))|
         (compa=="LM.LF"&score_type%in%c("hypoM","DMG"))|
         (compa=="C.L"&score_type%in%c("hyperM","DMG")),
       in.all.LGAM:=.N==8,"Description"]
  unique(resK[in.all.LGAM==T]$Description)
  
  #exclu femelle : 
  setdiff(unique(resK[in.all.LGAF==T]$Description),unique(resK[in.all.LGAM==T]$Description))
  #exclu male : 
  setdiff(unique(resK[in.all.LGAM==T]$Description),unique(resK[in.all.LGAF==T]$Description))
  
  resK[,exclu.LGAM:=all(in.all.LGAM),"Description"]
  resK[exclu.LGAM==T]
  
  #for Male
  resK[(compa=="CF.LM"&score_type%in%c("hyperM","DMG"))|
         (compa=="CM.LF"&score_type%in%c("hyperM","DMG"))|
         (compa=="CF.CM"&score_type%in%c("hyperM","DMG"))
       (compa=="LM.LF"&score_type%in%c("hypoM","DMG")),
       in.all.LGAM:=.N==6,"Description"]
  resK[in.all.LGAM==T]
  resK[,exclu.LGAM:=all(in.all.LGAM),"Description"]
  resK[exclu.LGAM==T]
  
  #pathway sex 
  resK[compa=="CF.CM"&score_type=="DMG"]$Description
  
}


#2) Is there LGAF ++ impacted via un pathway ?
{
  #a) see pval pathway
  path<-c("FoxO signaling pathway","PI3K-Akt signaling pathway",
          "Longevity regulating pathway","Signaling pathways regulating pluripotency of stem cells",
          "Insulin signaling pathway","Alzheimer disease","Growth hormone synthesis, secretion and action",
          "Cushing syndrome","Thyroid hormone signaling pathway")
  unique(resK$compa)
  #pathway lgaF
  
  #pathway lgaM
  
  #b) keep genes ++ diff ctrl vs lga; lgaf vs lgaM (~delta GeneScore)
  
}

#3) Ouverture : Why some ctrl male have a "lga" profile ?


#4) COnfirm by Regulation of expression of genes of the pathways/BP ? => Key Genes to validate
#genes pathways identif are DEG ?


