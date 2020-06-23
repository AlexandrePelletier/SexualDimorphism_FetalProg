#2020-05-30
#config et library
options(stringsAsFactors=F)
set.seed(12345)
#library(limma)
library(data.table)

wb_eQTL<-fread("../../ref/eQTL/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt")
wb_eQTL
#stat :
nrow(wb_eQTL) #2 414 653 variant-gene pairs
length(unique(wb_eQTL$gene_id)) #dans 12360 gene differents 
plot(density(table(wb_eQTL$gene_id)))
summary(as.vector(table(wb_eQTL$gene_id)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    22.0    87.0   195.4   230.0  8571.0 
#en moy : 195 var par gene

#1) need convert ensg. in symbols
#with bitr
#need remove the ".x" in order to be understood by the translator
library(stringr)
library(biomartr)
length(unique(wb_eQTL$gene_id))
biomartr::getMarts()
dtset<-biomartr::getDatasets("ENSEMBL_MART_ENSEMBL")[]
dtset[str_detect(dtset$dataset,"hsapiens"),]
biomartr::getAttributes(mart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")


result_BM <- biomartr::biomart( genes      = sapply(unique(wb_eQTL$gene_id),function(x)strsplit(x,"\\.")[[1]][1]), # genes that we wanted info
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = "hgnc_symbol", # attributes were selected with biomartr::getAttributes()
                                filters    = "ensembl_gene_id")
#put in gene colonne
result_BM<-result_BM[result_BM$hgnc_symbol!="",]
result_BM<-data.table(result_BM)
result_BM<-result_BM[,gene_id:=ensembl_gene_id][,-"ensembl_gene_id"]
result_BM<-result_BM[,gene:=hgnc_symbol][,-"hgnc_symbol"]
result_BM
nrow(result_BM) #ajout 10730 symbol

wb_eQTL[,gene_id_version:=gene_id]

wb_eQTL[,gene_id:=sapply(gene_id,function(x)strsplit(x,"\\.")[[1]][1])]
wb_eQTL
wb_eQTL<-merge(wb_eQTL,result_BM,all.x=T,by="gene_id")
wb_eQTL[!is.na(gene)]


#2) need chr pos absolu
#get translator hg38>hg19
translator<-fread("../../ref/eQTL/GTEx_Analysis_v8_eQTL/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt",
                  select = c(1,8))
translator
wb_eQTL<-merge(wb_eQTL,translator,all.x=T,by="variant_id")

wb_eQTL[,chr:=str_extract(variant_id_b37,"^[0-9XY]{1,2}") ]
wb_eQTL<-wb_eQTL[!is.na(chr),]
wb_eQTL[,chr:=paste0("chr",chr) ]

wb_eQTL[,pos:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2)) ]
wb_eQTL

#3) def eQTR

#regarder de cb les eQTL sont espacés en moyenne :
distSuivant<-function(positions){
  return(c(sapply(1:(length(positions)-1), function(i)abs(positions[i]-positions[i+1])),NA))
}

wb_eQTL[,distDuSuivant:=distSuivant(tss_distance),by="gene_id"]

wb_eQTL

plot(density(na.omit(wb_eQTL$distDuSuivant[wb_eQTL$distDuSuivant<2000]))) #pique à 30 pb

#zoom
plot(density(na.omit(wb_eQTL$distDuSuivant[wb_eQTL$distDuSuivant<200&wb_eQTL$distDuSuivant>0]))) #pas d'autre pique signif

summary(wb_eQTL$distDuSuivant)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       0     113     394    2679    1153 1852417   12360 


plot(density(log10(unique(wb_eQTL$min_pval_nominal))))
summary(unique(wb_eQTL$min_pval_nominal))  #++0 =>  
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 0.000e+00 3.000e-10 6.107e-06 2.023e-06 3.395e-04 

head(wb_eQTL,100)

# def eQTR : a partir de most signif eQTL localement (a +/-5kb) : 
#def region de eQTL a <1500 pb les uns des autres sans chute de pvalue < 10-3, length reg min =1000, max reg =5000.

#i)def min local  : pval dans 1er quartile et  a plus de 5000pb d'un eqtl plus signif
wb_eQTL[,minLocal.cand:=pval_nominal<quantile(pval_nominal,0.25),by=c("gene_id")]

is.minLocal<-function(dists,pvals){
  isMins<-sapply(1:length(dists), function(i){
    locisA5kb<-which(dists>(dists[i]-5000)&dists<(dists[i]+5000))
    
    if(all(pvals[locisA5kb]>=pvals[i])){
      return(T)
    }else{
      return(F)
    }
  })
  
  return(isMins)
}
wb_eQTL[minLocal.cand==T,minLocal:=is.minLocal(tss_distance,pval_nominal),by=c("gene_id")]
head(wb_eQTL,100)
#ii)
findeQTR<-function(minLocal,pvals,dists,seed=750,log10.pval.dt.thr=3,length.max=5000,length.min=1000){
  reg<-0
  length.seed.reg<-1
  eQTRs<-rep(NA,length(pvals))
  iMinsLocals<-which(minLocal==T)
  for(i in iMinsLocals){
    #pour chaque minlocal, si pas deja affecté a une reg : 
    if(!(i%in%which(!is.na(eQTRs)))){
      # => une nouvelle reg, avec un nouveau pval.thr de seed : 
      reg<-reg+1
      
      pvals.thr<-10^(log10(pvals[i])+log10.pval.dt.thr)
      # si eQTL a une dist <seed =>  inclus dans reg
      idansSeedsReg<-which(dists>(dists[i]-seed)&dists<(dists[i]+seed))
      idansReg<-idansSeedsReg
      eQTRs[idansReg]<-reg
      #pour chaque nouvelle reg si length.seed.reg < length.max, regarde si eQTL a une dist <seed, si oui, :
      length.seed.reg<- max(dists[idansSeedsReg])-min(dists[idansSeedsReg])
      
      if(length.seed.reg>0){
        while(length.seed.reg<length.max){
          
          start.seed.reg<-min(dists[idansSeedsReg])
          end.seed.reg<-max(dists[idansSeedsReg])
          candIdansSeedsReg<-setdiff(which(dists>(start.seed.reg-seed)&dists<(end.seed.reg+seed)),idansSeedsReg)
          if(length(candIdansSeedsReg)>0){
            
            #si length.seed.reg <length.min => inclus dans seed.reg
            if(length.seed.reg<length.min){
              idansSeedsReg<-c(idansSeedsReg,candIdansSeedsReg)
              eQTRs[idansSeedsReg]<-reg
              start.seed.reg<-min(dists[idansSeedsReg])
              end.seed.reg<-max(dists[idansSeedsReg])
              length.seed.reg<-end.seed.reg-start.seed.reg
            }
            #si length.seed.reg > length.min => inclus mais stopseed si diff pval min <log10.pval.dt.thr
            else{
              #ajout eQTL dans reg mais pas oblogatoirement dans seeed : 
              idansReg<-c(idansSeedsReg,candIdansSeedsReg)
              
              #inclus dans seed si assez signif :
              sigIdansSeedsReg<-candIdansSeedsReg[which(pvals[candIdansSeedsReg]<pvals.thr)]
              #sil y a des signifs on continue le seed :
              if(length(sigIdansSeedsReg)>0){
                idansSeedsReg<-c(idansSeedsReg,sigIdansSeedsReg)
                eQTRs[idansSeedsReg]<-reg
                start.seed.reg<-min(dists[idansSeedsReg])
                end.seed.reg<-max(dists[idansSeedsReg])
                length.seed.reg<-end.seed.reg-start.seed.reg
              }else{
                #sinon, on inclus les eQTL dans reg et on stop seed
                eQTRs[idansReg]<-reg
                length.seed.reg<-length.max
                
              }
              
            }
            
            
          }else{
            length.seed.reg<-length.max
          }
        }
        #sinon, on inclus les eQTL dans reg et on stop seed
        eQTRs[idansReg]<-reg
        
        
      }
      
    }
    
    
  }
  
  
  
  
  return(eQTRs)
}

#test pour 1 
wb_eQTL[gene=="WASH6P",eQTR:=findeQTR(minLocal,pval_nominal,tss_distance,seed = 1500,log10.pval.dt.thr = 3,length.max=5000,length.min = 1000)]


#pour all
wb_eQTL[,eQTR:=findeQTR(minLocal,pval_nominal,tss_distance,seed = 1500,log10.pval.dt.thr = 3,length.max=5000,length.min = 1000),by="gene_id"]
wb_eQTL[!is.na(eQTR),is.start.eQTR:= tss_distance==min(tss_distance),by=.(gene_id,eQTR)]
wb_eQTL[!is.na(eQTR),is.end.eQTR:= tss_distance==max(tss_distance),by=.(gene_id,eQTR)]

head(wb_eQTL,100) #OK

plot(density(log10(unique(wb_eQTL$RegScore))))


#2) save pos most signif eQTL => already done


#3) add eQTRweight based on nb of CpG and pval
wb_eQTL[,RegScore:=mean(-log10(pval_nominal)),by=.(gene_id,eQTR)]

max(wb_eQTL$RegScore)#221
fwrite(wb_eQTL,"../../ref/2020-05-31_Whole_Blood.eQTL_with_eQTR.csv",sep=";")

#4) make eQTR :

wb_eQTL<-wb_eQTL[order(chr,gene_id,pos)]
head(wb_eQTL,100)
wb_eQTR<-wb_eQTL[!is.na(eQTR)]

wb_eQTR[,start.eQTR:=min(pos),by=.(gene_id,eQTR)]
wb_eQTR[,end.eQTR:=max(pos),by=.(gene_id,eQTR)]
wb_eQTR[,n.eQTL:=.N,by=.(gene_id,eQTR)]
wb_eQTR[,length.eQTR:=end.eQTR-start.eQTR]
summary(wb_eQTR$length.eQTR)
plot(density((wb_eQTR$length.eQTR)))
wb_eQTR<-wb_eQTR[minLocal==T]
wb_eQTR
wb_eQTRl<-wb_eQTR[,.(chr,pos,eQTR,start.eQTR,end.eQTR,n.eQTL,RegScore,gene,tss_distance,pval_nominal)][!is.na(gene)&gene!=""]
wb_eQTRl
#add +/-500pb pour eQTR with only 1 cpg
wb_eQTRl[n.eQTL==1,start.eQTR:=start.eQTR-500]
wb_eQTRl[n.eQTL==1,end.eQTR:=end.eQTR+500]
fwrite(wb_eQTRl,"../../ref/2020-05-31_Whole_Blood.eQTR_light.csv",sep=";")


#5) add locis
cpgs<-fread("../../ref/annotation_CpG_HELP_ALL_070420.csv")
cpgs<-cpgs[,.(locisID,chr,pos)]
cpgs
#remove empty annot
cpgs<-cpgs[chr!=""]
cpgs<-cpgs[pos!=""]
cpgs<-cpgs[!is.na(pos)]
cpgs
eQTRll<-wb_eQTRl[,.(chr,start.eQTR,end.eQTR,gene)]
rm(translator,wb_eQTL,eQTR,wb_eQTRl)
output<-"analyses/eQTL_integration/2020-05-31"
options(future.globals.maxSize = 4000 * 1024^2)
for(chrom in unique(cpgs$chr)){
  print(chrom)
  cpgsl<-cpgs[chr==chrom]
  start1<-0
  for(p in 1:100){
    f<-p/100
    print(paste(chrom,p,"%"))
    stop<-quantile(cpgsl$pos,f)
    cpgsll<-cpgsl[pos>start1&pos<=stop]
    print(paste(nrow(cpgsll),"locis"))
    eQTRlll<-eQTRll[chr==chrom&(start.eQTR%between%c(start1,stop)|end.eQTR%between%c(start1,stop))]
    print(paste(nrow(eQTRlll),"regions"))
    
    if(nrow(eQTRlll)>0){
      #initialise sRegs
      cpg<-cpgsll[1,]
      sRegs<-eQTRlll[(start.eQTR<cpg$pos)&(end.eQTR>cpg$pos)]
      sRegs[,locisID:=cpg$locisID]
      sRegs[,pos:=cpg$pos]
      
      for(i in 2:nrow(cpgsll)){
        cpg<-cpgsll[i,]
        sReg<-eQTRlll[(start.eQTR<cpg$pos)&(end.eQTR>cpg$pos)]
        sReg[,locisID:=cpg$locisID]
        sReg[,pos:=cpg$pos]
        sRegs<-rbind(sRegs,sReg)
      }
    }
    
    fwrite(sRegs,file=paste0(output,"temp",chrom,p))
    start1<-stop
  }
  
  
}
sRegs<-fread("analyses/eQTL_integration/2020-05-31tempchr101")

for(file in list.files("analyses/eQTL_integration/",pattern = "temp")){
  print(strsplit(file,"chr")[[1]][2])
  sReg<-fread(paste0("analyses/eQTL_integration/",file))
  sRegs<-rbind(sRegs,sReg)
  
}


sRegs
sRegs[duplicated(sRegs)]
sRegs[locisID==2190609]
sRegs<-unique(sRegs)
sRegs[locisID==2190609]

# add  eQTR info : 
#to avoid conflict, "pos" column of eQTR > "pos.eQTL"
wb_eQTRl<-fread("../../ref/2020-05-31_Whole_Blood.eQTR_light.csv")


wb_eQTRl<-wb_eQTRl[,pos.eQTL:=pos][,-"pos"]
#same for tss_distance
wb_eQTRl<-wb_eQTRl[,tss_dist.eQTL:=tss_distance][,-"tss_distance"]
wb_eQTRl<-unique(wb_eQTRl)
sRegs
CpG.regs<-merge(sRegs,wb_eQTRl,by=c('chr','start.eQTR','end.eQTR','gene'),all.x = T,allow.cartesian = T)

CpG.regs[locisID==2190609] #ok
CpG.regs

#add distTSS CpG
#1) dist from eQTL
CpG.regs[,eQTL_dist:=pos-pos.eQTL]
summary(CpG.regs$eQTL_dist)
plot(density(CpG.regs$eQTL_dist))
CpG.regs[eQTL_dist>50000] #il y a des cpg tres loins de l'eqtl, 
#on enleve CpGlinks si a > 5kb du most signif
CpG.regs<-CpG.regs[abs(eQTL_dist)<5000] 
CpG.regs #492262 cpg-gene links
plot(density(CpG.regs$eQTL_dist))
#2) dist from tss
CpG.regs[,tss_dist:=tss_dist.eQTL+eQTL_dist]
CpG.regs
summary(CpG.regs$tss_dist)
plot(density(CpG.regs$tss_dist))

#refine th dataframe before save
CpG.regs<-CpG.regs[,avg.mlog10.pv.eQTLs:=RegScore][,-"RegScore"]
plot(density(CpG.regs$avg.mlog10.pv.eQTLs))

summary(CpG.regs$avg.mlog10.pv.eQTLs)

#save this new CpG-Gene links
fwrite(CpG.regs,"../../ref/2020-06-01_CpG_Gene_links_based_on_whole_blood_eQTL.csv",sep=";")

#2020-06-01
library(data.table)
#4) add metaeQTL : on veut avoir une ref de asso eqtl-gene stable across tissues.
meta_eQTL<-fread("../../../../../Downloads/GTEx_Analysis_v8.metasoft.txt")
meta_eQTL 
# choose  RE2 as model stat : New random effects model optimized to detect associations under heterogeneity. (Han and Eskin, AJHG 2011)
#STAT2_RE2	RE2 statistic heterogeneity part
plot(density(meta_eQTL$STAT2_RE2))
# Q  stat is an Heterogeneity estimates : for each snp-gene pair, association is it heterogen between tissue of same sample ?
cor(meta_eQTL$Q,meta_eQTL$STAT2_RE2,use = "complete.obs") #0.985788
plot(density(na.omit(meta_eQTL$Q)))
summary(meta_eQTL$Q)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     0.0    57.6    81.8   123.5   130.0  9078.0  399262 

plot(density(na.omit(meta_eQTL[PVALUE_Q<0.01]$Q)))
summary(meta_eQTL[PVALUE_Q<0.01]$Q)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 6.635   82.215  108.403  161.000  168.524 9078.040 

#also Isquare to estimate heterogeneity
cor(meta_eQTL$I_SQUARE,meta_eQTL$STAT2_RE2,use = "complete.obs") #0.44
plot(density(na.omit(meta_eQTL$I_SQUARE)))
meta_eQTL[I_SQUARE>90] #Q, stat2_RE2 high too if ++ study, but not if peu de study => bonne metrique pour filitrer


#FILTRER POUR HETEROGENEITE TROP GRANDE
meta_eQTL[I_SQUARE>50]#6M / 13M
meta_eQTL<-meta_eQTL[I_SQUARE<50]

#select gighly conserved signif snp-gene pairs across tiisues with pvalue_RE2 
plot(density(log10(meta_eQTL$PVALUE_RE2))) 
#take only high significant asso : 
summary(-log10(meta_eQTL$PVALUE_RE2)) #3rd quarter will be on cutoff : ~10^-30
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000019  4.174217 11.871462       Inf 28.743995       Inf 

meta_eQTL<-meta_eQTL[PVALUE_RE2<10^-30,]
meta_eQTL #1599039 snp-gene pair
meta_eQTL[order(-PVALUE_RE2)]

#2020-06-19
library(data.table)
# add metaeQTL : on veut avoir une ref de asso eqtl-gene stable across tissues.
meta_eQTL<-fread("../../../../../Downloads/GTEx_Analysis_v8.metasoft.txt")
meta_eQTL 
#keep only no heterogen eQTL and signif
#heterogeneity thanks to Qvalue :
meta_eQTL_H<-meta_eQTL[PVALUE_Q>0.1] #2.4M
#with number of study>40
meta_eQTL_HS<-meta_eQTL_H[`#STUDY`>40] #1.8M
meta_eQTL_HS
#PVALUE_FE <-30
plot(density(-log10(meta_eQTL_HS$PVALUE_FE)))
abline(v=30)
summary(-log10(meta_eQTL_HS$PVALUE_FE))
meta_eQTL_HSP<-meta_eQTL_HS[PVALUE_FE<10e-30] #500k


#recover col gene, chr, pos (in hg37) 
#1) in hg38 :
meta_eQTL_HSP[,chr:=sapply(RSID,function(x){
  return(strsplit(x,"_")[[1]][1])
})]

meta_eQTL_HSP[,pos:=sapply(RSID,function(x){
  return(as.numeric(strsplit(x,"_")[[1]][2]))
})]
meta_eQTL_HSP

meta_eQTL_HSP[,gene_id:=sapply(RSID,function(x){
  return(strsplit(x,",")[[1]][2])
})]
meta_eQTL_HSP
length(unique(meta_eQTL_HSP$gene_id)) #6884 gene

fwrite(meta_eQTL_HSP,"../../ref/eQTL/2020-06-19_meta_eQTL_filtered_and_annotated.csv",sep=";")

meta_eQTL_HSP<-fread("../../ref/eQTL/2020-06-19_meta_eQTL_filtered_and_annotated.csv")

#trans in symbol : #to optimize with hsapiens_gene_ensembl_version when biomartr server will work
ref<-fread("../../../Alexandre_SC/ref/ENSEMBL_ID_to_SYMBOL.csv")
ref
meta_eQTL_HSP$gene<-ref$hgnc_symbol[match(sapply(meta_eQTL_HSP$gene_id,function(x)strsplit(x,"\\.")[[1]][1]),ref$ensembl_gene_id)]
meta_eQTL_HSP[is.na(gene)] #200k/500k

length(unique(meta_eQTL_HSP[is.na(gene)]$gene_id))#1170



#stat :
nrow(meta_eQTL_HSP) #596030 variant-gene pairs
length(unique(meta_eQTL_HSP$gene_id)) #dans 6884 gene differents 
plot(density(table(meta_eQTL_HSP$gene_id)))
summary(as.vector(table(meta_eQTL_HSP$gene_id)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    4.00   17.00   86.58   66.00 7656.00 
#en moy : 86 var par gene
meta_eQTL_HSP[!is.na(gene)] #395496
length(unique(meta_eQTL_HSP[!is.na(gene)]$gene))#4843


#PUT IN hg19 format
#get translator hg38>hg19
library(stringr)
translator<-fread("../../ref/eQTL/GTEx_Analysis_v8_eQTL/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt",
                  select = c(1,8))
translator
meta_eQTL_HSP[,variant_id:=sapply(RSID,function(x){
  return(strsplit(x,",")[[1]][1])
})]


meta_eQTL_HSP<-merge(meta_eQTL_HSP,translator,all.x=T,by="variant_id")

meta_eQTL_HSP[,chr:=str_extract(variant_id_b37,"^[0-9XY]{1,2}") ]
meta_eQTL_HSP<-meta_eQTL_HSP[!is.na(chr),]
meta_eQTL_HSP[,chr:=paste0("chr",chr) ]

meta_eQTL_HSP[,pos:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2)) ]
meta_eQTL_HSP
#save
fwrite(meta_eQTL_HSP,"../../ref/eQTL/2020-06-19_meta_eQTL_filtered_and_annotated.csv",sep=";")



#2020-06-23
library(data.table)
meta_eQTL_HSP<-fread("../../ref/eQTL/2020-06-19_meta_eQTL_filtered_and_annotated.csv")
#get only pvalue
meta_eQTL_HSP<-meta_eQTL_HSP[,.(PVALUE_FE,chr,pos,gene_id,gene)]
unique(meta_eQTL_HSP,by="gene_id") #◘880
sum(unique(meta_eQTL_HSP,by="gene_id")$gene!="") #4840
meta_eQTL_HSP[gene==""] # a lot
#so new translate 
library(biomartr)
biomartr::getAttributes(mart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
result_BM <- biomartr::biomart( genes      = unique(meta_eQTL_HSP$gene_id), # genes that we wanted info
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = "hgnc_symbol", # attributes were selected with biomartr::getAttributes()
                                filters    = "ensembl_gene_id_version")

nrow(result_BM) #2707, pas mieux
meta_eQTL_HSPG<-meta_eQTL_HSP[gene!=""]
meta_eQTL_HSPG<-meta_eQTL_HSPG[,-"gene_id"]
#3) def eQTR

#regarder de cb les eQTL sont espacés en moyenne :
distSuivant<-function(positions){
  return(c(sapply(1:(length(positions)-1), function(i)abs(positions[i]-positions[i+1])),NA))
}

meta_eQTL_HSPG[,distDuSuivant:=distSuivant(pos),by="gene"]



plot(density(na.omit(meta_eQTL_HSPG$distDuSuivant[meta_eQTL_HSPG$distDuSuivant<2000]))) #pique à 30 pb

#zoom
plot(density(na.omit(meta_eQTL_HSPG$distDuSuivant[meta_eQTL_HSPG$distDuSuivant<200&meta_eQTL_HSPG$distDuSuivant>0]))) #pas d'autre pique signif

summary(meta_eQTL_HSPG$distDuSuivant)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       0     113     394    2679    1153 1852417   12360 
meta_eQTL_HSPG[distDuSuivant==0]

plot(density(log10(unique(meta_eQTL_HSPG$PVALUE_FE))))
summary(unique(meta_eQTL_HSPG$PVALUE_FE))   
Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.000e+00 0.000e+00 0.000e+00 1.836e-31 3.000e-36 9.992e-30

head(meta_eQTL_HSPG,100)

# def eQTR : a partir de most signif eQTL localement (a +/-5kb) : 
#def region de eQTL a <1500 pb les uns des autres sans chute de pvalue < 10-3, length reg min =1000, max reg =5000.

#i)def min local  : pval dans 1er quartile et  a plus de 5000pb d'un eqtl plus signif
meta_eQTL_HSPG[,minLocal.cand:=PVALUE_FE<=quantile(PVALUE_FE,0.25),by=c("gene")]

is.minLocal<-function(dists,pvals){
  isMins<-sapply(1:length(dists), function(i){
    locisA5kb<-which(dists>(dists[i]-5000)&dists<(dists[i]+5000))
    
    if(all(pvals[locisA5kb]>=pvals[i])){
      return(T)
    }else{
      return(F)
    }
  })
  
  return(isMins)
}
meta_eQTL_HSPG[minLocal.cand==T,minLocal:=is.minLocal(pos,PVALUE_FE),by=c("gene")]
head(meta_eQTL_HSPG,100)
#ii)
findeQTR<-function(minLocal,pvals,dists,seed=750,log10.pval.dt.thr=3,length.max=5000,length.min=1000){
  reg<-0
  length.seed.reg<-1
  eQTRs<-rep(NA,length(pvals))
  iMinsLocals<-which(minLocal==T)
  for(i in iMinsLocals){
    #pour chaque minlocal, si pas deja affecté a une reg : 
    if(!(i%in%which(!is.na(eQTRs)))){
      # => une nouvelle reg, avec un nouveau pval.thr de seed : 
      reg<-reg+1
      
      pvals.thr<-10^(log10(pvals[i])+log10.pval.dt.thr)
      # si eQTL a une dist <seed =>  inclus dans reg
      idansSeedsReg<-which(dists>(dists[i]-seed)&dists<(dists[i]+seed))
      idansReg<-idansSeedsReg
      eQTRs[idansReg]<-reg
      #pour chaque nouvelle reg si length.seed.reg < length.max, regarde si eQTL a une dist <seed, si oui, :
      length.seed.reg<- max(dists[idansSeedsReg])-min(dists[idansSeedsReg])
      
      if(length.seed.reg>0){
        while(length.seed.reg<length.max){
          
          start.seed.reg<-min(dists[idansSeedsReg])
          end.seed.reg<-max(dists[idansSeedsReg])
          candIdansSeedsReg<-setdiff(which(dists>(start.seed.reg-seed)&dists<(end.seed.reg+seed)),idansSeedsReg)
          if(length(candIdansSeedsReg)>0){
            
            #si length.seed.reg <length.min => inclus dans seed.reg
            if(length.seed.reg<length.min){
              idansSeedsReg<-c(idansSeedsReg,candIdansSeedsReg)
              eQTRs[idansSeedsReg]<-reg
              start.seed.reg<-min(dists[idansSeedsReg])
              end.seed.reg<-max(dists[idansSeedsReg])
              length.seed.reg<-end.seed.reg-start.seed.reg
            }
            #si length.seed.reg > length.min => inclus mais stopseed si diff pval min <log10.pval.dt.thr
            else{
              #ajout eQTL dans reg mais pas oblogatoirement dans seeed : 
              idansReg<-c(idansSeedsReg,candIdansSeedsReg)
              
              #inclus dans seed si assez signif :
              sigIdansSeedsReg<-candIdansSeedsReg[which(pvals[candIdansSeedsReg]<pvals.thr)]
              #sil y a des signifs on continue le seed :
              if(length(sigIdansSeedsReg)>0){
                idansSeedsReg<-c(idansSeedsReg,sigIdansSeedsReg)
                eQTRs[idansSeedsReg]<-reg
                start.seed.reg<-min(dists[idansSeedsReg])
                end.seed.reg<-max(dists[idansSeedsReg])
                length.seed.reg<-end.seed.reg-start.seed.reg
              }else{
                #sinon, on inclus les eQTL dans reg et on stop seed
                eQTRs[idansReg]<-reg
                length.seed.reg<-length.max
                
              }
              
            }
            
            
          }else{
            length.seed.reg<-length.max
          }
        }
        #sinon, on inclus les eQTL dans reg et on stop seed
        eQTRs[idansReg]<-reg
        
        
      }
      
    }
    
    
  }
  
  
  
  
  return(eQTRs)
}

#test pour 1 
meta_eQTL_HSPG[gene=="ABCC2"]
meta_eQTL_HSPG[gene=="ABCC2",eQTR:=findeQTR(minLocal,PVALUE_FE,pos,seed = 1500,log10.pval.dt.thr = 3,length.max=5000,length.min = 1000)]
head(meta_eQTL_HSPG[gene=="ABCC2"],100)

#pour all
meta_eQTL_HSPG[,eQTR:=findeQTR(minLocal,PVALUE_FE,pos,seed = 1500,log10.pval.dt.thr = 3,length.max=5000,length.min = 1000),by="gene"]

meta_eQTL_HSPG[!is.na(eQTR),is.start.eQTR:= pos==min(pos),by=.(gene,eQTR)]
meta_eQTL_HSPG[!is.na(eQTR),is.end.eQTR:=pos==max(pos),by=.(gene,eQTR)]
head(meta_eQTL_HSPG,100) #OK

fwrite(meta_eQTL_HSPG,"../../ref/2020-06-23_meta.eQTL_with_eQTR.csv",sep=";")

#3) add eQTRweight based on nb of CpG and pval
meta_eQTL_HSPG[,RegScore:=mean(-log10(PVALUE_FE)),by=.(gene,eQTR)]
# [stop here : eviter inf dans regScore ]

max(meta_eQTL_HSPG$RegScore)#Inf

plot(density(log10(unique(meta_eQTL_HSPG$RegScore))))
fwrite(meta_eQTL_HSPG,"../../ref/2020-05-31_Whole_Blood.eQTL_with_eQTR.csv",sep=";")


#4) make eQTR :


meta_eQTL_HSPG<-meta_eQTL_HSPG[order(chr,gene_id,pos)]
head(meta_eQTL_HSPG,100)
meta_eQTR<-meta_eQTL_HSPG[!is.na(eQTR)]

meta_eQTR[,start.eQTR:=min(pos),by=.(gene_id,eQTR)]
meta_eQTR[,end.eQTR:=max(pos),by=.(gene_id,eQTR)]
meta_eQTR[,n.eQTL:=.N,by=.(gene_id,eQTR)]
meta_eQTR[,length.eQTR:=end.eQTR-start.eQTR]
summary(meta_eQTR$length.eQTR)
plot(density((meta_eQTR$length.eQTR)))
meta_eQTR<-meta_eQTR[minLocal==T]
meta_eQTR
meta_eQTRl<-meta_eQTR[,.(chr,pos,eQTR,start.eQTR,end.eQTR,n.eQTL,RegScore,gene,tss_distance,pval_nominal)][!is.na(gene)&gene!=""]
meta_eQTRl
#add +/-500pb pour eQTR with only 1 cpg
meta_eQTRl[n.eQTL==1,start.eQTR:=start.eQTR-500]
meta_eQTRl[n.eQTL==1,end.eQTR:=end.eQTR+500]
fwrite(meta_eQTRl,"../../ref/2020-05-31_Whole_Blood.eQTR_light.csv",sep=";")


#5) add locis
cpgs<-fread("../../ref/annotation_CpG_HELP_ALL_070420.csv")
cpgs<-cpgs[,.(locisID,chr,pos)]
cpgs
#remove empty annot
cpgs<-cpgs[chr!=""]
cpgs<-cpgs[pos!=""]
cpgs<-cpgs[!is.na(pos)]
cpgs
eQTRll<-meta_eQTRl[,.(chr,start.eQTR,end.eQTR,gene)]
rm(translator,meta_eQTL_HSPG,eQTR,meta_eQTRl)
output<-"analyses/eQTL_integration/2020-05-31"
options(future.globals.maxSize = 4000 * 1024^2)
for(chrom in unique(cpgs$chr)){
  print(chrom)
  cpgsl<-cpgs[chr==chrom]
  start1<-0
  for(p in 1:100){
    f<-p/100
    print(paste(chrom,p,"%"))
    stop<-quantile(cpgsl$pos,f)
    cpgsll<-cpgsl[pos>start1&pos<=stop]
    print(paste(nrow(cpgsll),"locis"))
    eQTRlll<-eQTRll[chr==chrom&(start.eQTR%between%c(start1,stop)|end.eQTR%between%c(start1,stop))]
    print(paste(nrow(eQTRlll),"regions"))
    
    if(nrow(eQTRlll)>0){
      #initialise sRegs
      cpg<-cpgsll[1,]
      sRegs<-eQTRlll[(start.eQTR<cpg$pos)&(end.eQTR>cpg$pos)]
      sRegs[,locisID:=cpg$locisID]
      sRegs[,pos:=cpg$pos]
      
      for(i in 2:nrow(cpgsll)){
        cpg<-cpgsll[i,]
        sReg<-eQTRlll[(start.eQTR<cpg$pos)&(end.eQTR>cpg$pos)]
        sReg[,locisID:=cpg$locisID]
        sReg[,pos:=cpg$pos]
        sRegs<-rbind(sRegs,sReg)
      }
    }
    
    fwrite(sRegs,file=paste0(output,"temp",chrom,p))
    start1<-stop
  }
  
  
}
sRegs<-sRegs[1,]

for(file in list.files("analyses/eQTL_integration/",pattern = "temp")){
  print(strsplit(file,"chr")[[1]][2])
  sReg<-fread(paste0("analyses/eQTL_integration/",file))
  sRegs<-rbind(sRegs,sReg)
  
}


sRegs
sRegs[duplicated(sRegs)]
sRegs[locisID==2190609]
sRegs<-unique(sRegs)
sRegs[locisID==2190609]

# add  eQTR info : 
#to avoid conflict, "pos" column of eQTR > "pos.eQTL"
meta_eQTRl<-fread("../../ref/2020-05-31_Whole_Blood.eQTR_light.csv")


meta_eQTRl<-meta_eQTRl[,pos.eQTL:=pos][,-"pos"]
#same for tss_distance
meta_eQTRl<-meta_eQTRl[,tss_dist.eQTL:=tss_distance][,-"tss_distance"]
meta_eQTRl<-unique(meta_eQTRl)
sRegs
CpG.regs<-merge(sRegs,meta_eQTRl,by=c('chr','start.eQTR','end.eQTR','gene'),all.x = T,allow.cartesian = T)

CpG.regs[locisID==2190609] #ok
CpG.regs

#add distTSS CpG
#1) dist from eQTL
CpG.regs[,eQTL_dist:=pos-pos.eQTL]
summary(CpG.regs$eQTL_dist)
plot(density(CpG.regs$eQTL_dist))
CpG.regs[eQTL_dist>50000] #il y a des cpg tres loins de l'eqtl, 
#on enleve CpGlinks si a > 5kb du most signif
CpG.regs<-CpG.regs[abs(eQTL_dist)<5000] 
CpG.regs #492262 cpg-gene links
plot(density(CpG.regs$eQTL_dist))
#2) dist from tss
CpG.regs[,tss_dist:=tss_dist.eQTL+eQTL_dist]
CpG.regs
summary(CpG.regs$tss_dist)
plot(density(CpG.regs$tss_dist))
#eQTL-CpG_Score = -log10(pval)* distScore [1 if +/-500, 0.9 if +/- 1000, 0.8 +/-1500, 0.7 +/-2k,0.6+/-2500]


CpG.regs[,distScore:=sapply(abs(eQTL_dist),function(x){
  if(x<500){
    distScore<-1 
  }else {
    distScore<-1-0.4*x/2500
  }
  
  return(distScore)
})]

plot(density(CpG.regs$eQTL.CpGLinks)) 
#GeneLinks-CpG_Score = log10((RegScore*distScore))
CpG.regs[,Gene.CpGLinks:=log10(RegScore*distScore)] #integrate with RegScore : mean of pval of eQTL in the region
plot(density(CpG.regs$Gene.CpGLinks)) 

#save this new CpG-Gene links
fwrite(CpG.regs,"../../ref/2020-06-01_CpG_Gene_links_based_on_whole_blood_eQTL.csv",sep=";")
