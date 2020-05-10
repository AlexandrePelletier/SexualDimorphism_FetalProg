#determiner comment utiliser GSEA pour nos données de methylation

#config et library
options(stringsAsFactors=F)
set.seed(12345)
library(limma)
source("scripts/scRNA-seq integration.R")
source("scripts/utils.R")
#output dir 
script_name <- "GSEA_test"
outputDir <- file.path("analyses",script_name)
dir.create(path = outputDir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
date<-Sys.Date()
output <- file.path(outputDir,date)

#selon gsea-msigDB, pour que GSEA donne des résultats relavant pour des datas autres que RNA-seq, il est préférable de 
#1) ranker les genes "in order of most (largest value) to least (smallest value) "of interest" 
#2) si metric pas relevant biologiquement, vaux mieux juste mettre le rank en metrique

#> donc, nous chercons a faire un score permettant de trier les CpG selon leur chance d'etre rellement impactant sur l'expression du gène relié.
#pour cela, un cpg est le plus susceptible d'influencer un gène si :
#1) il a un FC elevé : 
#2) la difference avec les CTRL est très significative 
#3) il est sur un prom ou un enh(4, 6)
#4) il est sur un zone définit comme régulatrice dans ENSEMBL(ou UCSC ?)
#5) il est a proximité d'un motifTF
#6) son changement de niveau de méhylation est  cross corrélé avec d'autres CpG 


#=> ça nous donne un score par CpG (CpGScore) susceptible d'influencé tel ou tel gène.
#Aves ces scores, nous créons ensuite un score par gène :
#1) en faisant la somme des CpGscores lié au gène* normalisé par le nb totale de CpG lié à ce gène
#2) pour valoriser les gènes réellement exprimé dans les HSPC** (i.e dont la methylation est reelement susceptible d'impacter son expression),
#bonus1) il est dans une DMR
#bonus2) si le genes expr est tres différents entre les pop ?

#ce score est ensuite amoindrit si le gène n'est pas du tout retrouvé exprimés dans nos datas RNA-seq de référence**

#*si améliore par la suite le lien CpG-gene, en liant les CpG a plusieurs gènes par esemble (pas seulement le plus proche), on pourra facielement donné
#**  datas RNA-seq de HSPC publique, et nos datas (sc)RNAseq de  CTRL et LGA.

#le FC, Pval, si sur prom ou sur enh,  la distance du TSS, 
#with clusterProf, 

fit2<-readRDS("analyses/main/fit2_model13.rds")
colnames(fit2$coefficients)
res<- topTable(fit2,coef="FC.FL",n =Inf)
locis<-rownames(res)

annot<-read.csv2("../../ref/annotation_CpG_HELP_ALL_070420.csv",row.names = 1)
head(annot)
nrow(annot)
#annot res :
res<-data.frame(row.names = locis,
                FC=res$logFC,
                pval=res$P.Value,
                pval.adj=res$adj.P.Val,
                annot[locis,c("chr","start","posAvant","gene","ENTREZID","type","feature_type_name")])
#on enleve les locis non rattaché a des genes
sum(is.na(res$gene)) #5971 > 340026 why ?
res<-res[!is.na(res$gene),]
nrow(res)

#now, we can create our CpGscore : 
#we want FC=pvalue > 4/6 > regulatory element > TF > wieght of correlation > RNAseq reference
#so principal weight is FC[0:1] + pval[0:1] (max=2) (the DMC weigtht)
#We multiply this weight by the "regulatory weight" [0-1]
# +1 if (4,6), 0.5 if 5, 0.1 otherwise
# +0.7 if(CTCF,enh,prom,)+ 0.3 if prom flanking, openchrine, #+0.5 if(TFbind site), or 0 otherwise
#rank this and norm the rank to obtain a score between 0 and 1
 
##DMC weigtht
#(FC[0:1] : toavoid overfitting, we use the rank of abs(FC) as metriks normalized
res$FCScore<-rank(abs(res$FC))/nrow(res)
head(res)
tail(res)
#pval[0:1] : the same (with rank), to have the same weight
res$PvalScore<-rank(1/res$pval)/nrow(res)
res$DMCWeight<-res$PvalScore+res$FCScore

##regulatory weight : je veux un poids entre 0.5 et 2, 0.5 etant aucune annotation regulatrice,
#et 2 etant toutes les annotations regulatrice (4/6, CTCF/prom/enh/, TF motif) 
#=> on part d'un score de 0.5 et +0.5 a chaque annot

# +1 if (4,6), +0.75 if 5, +0.5 otherwise
res$TypeScore<-sapply(res$type,function(x){
  if(x%in%c(4,6)){
    score<-1
  }else if(x==5){
    score<-0.75
  }else{
    score<-0.5
  }
  return(score)
})

# +0.5 if(CTCF,enh,prom,)+ 0.25 if prom flanking, openchrine, #+0.5 if(TFbind site), or 0 otherwise
res$EnsRegScore<-sapply(res$feature_type_name,function(x){
  vecX<-tr(x)
  if(any(c("CTCF Binding Site","Promoter","Enhancer")%in%vecX)){
    score<-0.5
  }else if(any(c("Open chromatin","Promoter Flanking Region")%in%vecX)){
    score<-0.25
  }else{
    score<-0
  }
  if("TF binding site"%in%vecX){
    score<-score+0.5
  }
  return(score)
})


# add the 2 score to make the regulWeight:
res$RegWeight<-res$TypeScore+res$EnsRegScore
head(res)
tail(res)

#ConfWeight :
# confWeight [0-1] : TSS, eQTL links, (cross correl ?)
#rank normalisé de sum de : 
#tss +1 si <1000pb, descend progressivement sinon  : 
res$TSSScore<-sapply(abs(res$posAvant),function(x){
  if(x<1000){
    score<-1
  }else {
    score<-1/(1+log10(x-1000))
  }
  return(score)
})
head(res)
tail(res)

#eQTL link : si locis +/-xpb rattaché à expression genes => +1
#with gTex : 
wb_eQTL<-fread("../../ref/Whole_Blood.eQTL.csv")
wb_eQTL
head(wb_eQTL)
genes_eQTL<-unique(wb_eQTL$gene)

##Remarque : refine the links CpG Gene with eQTL : 


# +1 si dans regionReg, descend progressivement sinon  :
#dabord checker si locis bien lié au gene : 
#trad en pos absolu : 
head(wb_eQTL) # dans le variant_id, il y a la position sous la forme "chr1_64764"
#c'est la pos de quel référence ? 


res$eQTLScore<-sapply(rownames(res),function(locus){
  gene<-res[locus,"gene"]
  if(gene%in%genes_eQTL){
    pos<-res[locus,"posAvant"]
    regions<-strsplit(wb_eQTL$region[which(wb_eQTL$gene==gene)],":")
    #!reflechir comment checker efficacement si dans les régions
    if(any(sapply(regions,function(region){
      return(pos>=as.numeric(region[1])&pos<=as.numeric(region[2]))
      
    }))){
      score<-1
    }else{
      distMin<-min(sapply(regions,function(region){
        start<-as.numeric(region[1])
        end<-as.numeric(region[2])
        if(pos<start){
          
          dist<-start-pos
        }else{
          dist<-pos-end
        }
        return(dist)
        
      }))
      score<-1/(1+log10(distMin))
    }
      
  }else{
    score<-0
  }
  
  
  return(score)
})
head(res,20)
annot[rownames(res),"eQTLScore"]<-res$eQTLScore
write.csv2(annot,"../../ref/annotation_CpG_HELP_ALL_070420.csv")

#cross correl : si cpg bouge pareillement a d'autres cpg a travers ech sans etre corrélé a group_complecity => poids 0>1



##CpG score: DMCweight*RegWeight*ConfWeight
res$CpGScore<-res$DMCWeight*res$RegWeight*res$ConfWeight

#now we can calculate de GenesScores by sum(CpGScore in a gene)/nCpGtot in the genes and add a weight of this new score depending of the gene expr: 

#first, det the res by CpG before to put res by genes :
head(res)
resGenes<-resLocisToGenes(res,withNCpGTot = T)
head(resGenes) #get also the number CpG detected in at least 1 sample HPA2c library(nCpG), compared the number of CpG tot detected by Msp1c ref
for(gene in rownames(resGenes)){
  resGenes[gene,"GeneScore"]<-sum(res$CpGScore[which(res$gene==gene)])
}
head(resGenes)
#GSEA without weitghed for non expr genes :
library(clusterProfiler)
KeggGS <- read.gmt("../../ref/c2.cp.kegg.v7.1.symbols.gmt") 
geneList<-makeGeneList(resGenes,score = "GeneScore",returnRank = T)
head(geneList,50)

kk<-GSEA(geneList,TERM2GENE=KeggGS,nPerm = 1000,pvalueCutoff = 0.1)
head(kk,20) #done
emapplot(kk)

#normalized with the number of CpG in data_F by genes
for(gene in rownames(resGenes)){
  resGenes[gene,"GeneScore2"]<-sum(res$CpGScore[which(res$gene==gene)])/resGenes[gene,"nCpG"]
}
head(resGenes)
#GSEA without weitghed for non expr genes :

geneList<-makeGeneList(resGenes,score = "GeneScore2",returnRank = T)
head(geneList,50)

kk2<-GSEA(geneList,TERM2GENE=KeggGS,nPerm = 1000,pvalueCutoff = 0.1)
head(kk2,20) #done
emapplot(kk2)


#normalized with the number of CpG detected by msp1 by genes
for(gene in rownames(resGenes)){
  resGenes[gene,"GeneScore3"]<-sum(res$CpGScore[which(res$gene==gene)])/resGenes[gene,"nCpGtot"]
}
head(resGenes)
#GSEA without weitghed for non expr genes :

geneList<-makeGeneList(resGenes,score = "GeneScore3",returnRank = T)
head(geneList,50)

kk3<-GSEA(geneList,TERM2GENE=KeggGS,nPerm = 1000,pvalueCutoff = 0.1)
head(kk3,20) #done
emapplot(kk3)

saveRDS(list(resLocis=res,resGenes=resGenes),paste0(output,"resCLF.rds"))



#need better norm than just divide by number of CpG tot, that too influencing
#moyenne pondéré /2: pour eviter inflaté genes avec beaucoup de locis, on veut ponderé le score des CpG :
#un gene est considéré comme tres impacté si il a 3 CpG important d'activer : (par exmple 1 sur enh, prom, genebody)
#donc, on notre moyenne pondéré favorisera surout les 3 topCpG score par gene : 
  #ScoreGene=CpG1+CpG2+CpG3 + 0.5*CpG4 + 0.2*CpG5 ... / (sommeCoeffTot)

#mais nb de CpG important dépend de cb CpG ont été rattaché au gene à la base : 
#donc plutot que 3 disons qu'il y a x CpG important par gene, pour trouver x, nous devons prendre en compte : 
#1) nb de CpG par gene : si peu de cpg rattaché, il ya moins de chance d'y avoir des cpg importants quand qd il yen a bcp
#2) taile du genes : plus il est grand plus on s'attend a qu'il est plus de cpg pour le reguler

#1) x=3 pour gene avec nCpG dans lespace interquartile , 2/4 si <q25/>q75, 1/5 si <q5/>q95
#2) +1/-1 si taille gene <q25/>q75

#ebayes/methylGSA ?


#we can add RNAseq expr to improve the score: 
CTRL_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInCTRL.csv",row.names = 1)
LGA_expr_by_pop<-read.csv2("analyses/test_scRNAseq_integration/geneExprInLGA.csv",row.names = 1)
listExprMat<-list(CTRL=CTRL_expr_by_pop,LGA=LGA_expr_by_pop)
#add genes Expr in bulk scrnaseq and in subpop
resGenes<-find.scGeneExpr(resGenes,listExprMat)
head(resGenes)#keep only cols nCpG and Bulk_expr in ctrl and LGA
#add HSPC pub expr in resGENEs
resGenes$HSPC1<-annot[match(rownames(resGenes),annot$gene),"HSPC1"]
resGenes$HSPC2<-annot[match(rownames(resGenes),annot$gene),"HSPC2"]
head(resGenes)




                
                


