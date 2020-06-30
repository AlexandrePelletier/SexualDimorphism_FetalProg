#2020-06-30
#deter Sex-Spe Methhylated Genes by Permutation

#deter genescore cutoff x / deter EpigenAffGene (EAG) : 
# 1000 permuts of all samples CTRL / LGA => LIMMA => pval and FC => df GeneScore => save pct genescore permuts >= geneScore obs
#=> gene score cutoff = Accept Gene in EAF if appeared < x%  => genesF and genesM spe
#make a permutFunction

deterEpigenAffGene<-function(res_to_compas,var_to_permut,permut=1000,pvalThr=0.01,){
  library(data.table)
  source("scripts/calculeGeneScore.R")
  batch<-fread("../../ref/cleaned_batch_CD34_library_date_220620.csv")
  batch<-batch[Group_name%in%c("C","L")]
  varNumToModel<-c("Mat.Age")
  varFacToModel<-c("Group_Sex",'batch',"latino","Group_Complexity_Fac")
  varToModel<-c(varNumToModel,varFacToModel)
  sample_F<-batch$sample[!(apply(is.na(batch[,..varToModel]),1,any))]
  
  print(paste(length(sample_F),"samples to permut"))
  
  batch<-batch[sample%in%sample_F,c("sample",..varToModel),]
  
  batch[,(varNumToModel):=lapply(.SD,as.numeric),.SDcols=varNumToModel]
  batch[,(varFacToModel):=lapply(.SD,as.factor),.SDcols=varFacToModel]
  methyl_df<-fread("../../ref/2020-05-25_methyl_data_before_limma.csv")
  methyl_df<-data.frame(methyl_df)
  rownames(methyl_df)<-methyl_df$locisID
  
  formule<- ~0 + Group_Sex  + batch  + latino + Mat.Age + Group_Complexity_Fac
  compas<-unlist(lapply(strsplit(names(res_to_compas),"-"),function(x)paste(paste0("Group_Sex",x),collapse = "-")))
  cpg.regs_ref<-fread("../../ref/2020-06-29_All_CpG-Gene_links.csv")
  res_to_compas<-lapply(res_to_compas,
                        CalcGeneScore,cpg.regs_ref,sumToGene=T)
  
  
  
  
  return(res_to_compas)
}




MethAnalysisExec<-function(methyl_df,batch_F,formule,compas,abbrev_compas,cpg.regs_ref,analysis="obs"){
  library(limma)
  design<-model.matrix(formule,data = batch_F)
  
  fit <- lmFit(methyl_df[,batch_F$sample], design)  
  
  cont.matrix <- makeContrasts(contrasts = compas,
                               levels=design)
  colnames(cont.matrix)<-abbrev_compas
  
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2) 
  res_list<-lapply(abbrev_compas,function(compa){
    res<-topTable(fit2,coef=compa,n =Inf)
    res<-CalcGeneScore(res,cpg.regs_ref,sumToGene=T)
    res<-res[order(gene)][,.(gene,GeneScore)]
    col<-paste(compa,analysis,sep = ".")
    return(res[,(col):=as.integer(GeneScore)][,-"GeneScore"])
  })
  return(Reduce(merge,res_list))
  
}

#run methylAnalysis for obs :
res<-MethAnalysisExec(methyl_df, batch,formule, compas,abbrev_compas,cpg.regs_ref )

res_all<-res
res_all
batch1<-data.table(batch)
set.seed(123)
j<-1
for(i in 1:1000){
  print(paste(i,"/",1000))
  
  batch1[,Group_Sex:=sample(Group_Sex)]
  
  res_all<-merge(res_all,MethAnalysisExec(methyl_df, batch1,formule, compas,abbrev_compas,cpg.regs_ref,analysis = i))
  if(ncol(res)>200){
    print(paste("save",j))
    fwrite(res_all,paste0("2020-06-23_temp",j))
    res_all<-res
    j<-j+1
  }
  
}
fwrite(res_all,paste0("2020-06-23_temp",j))
res_all
print("save  pct genescore permuts >= geneScore obs")


MethAnalysisPermut<-function(batch_F,col_to_permut,n,sample_col="sample",cpg.regs_ref=NULL,formule=NULL,seed=123){
  set.seed(seed)
  library(data.table)
  
  
  if(is.null(cpg.regs_ref)){
    cpg.regs_ref_path<-"../../ref/2020-06-03_All_CpG-Gene_links.csv"
    print(paste("recup cpg.gene regul annotation in ",cpg.regs_ref_path))
    cpg.regs_ref<-fread("../../ref/2020-06-03_All_CpG-Gene_links.csv")
  }
  
  if(is.null(formule)){
    formule<- ~0 + Group_Sex  + batch  + latino + Mat.Age + Group_Complexity_Fac
    print("limma model selon cette formule :")
    print(formule)
  }
  
  print("Limma execution of observed batch..")
  res<-LimmaExec(methyl_data,batch_F,formule)
  res<-CalcGeneScore(res,cpg.regs_ref)
  resGenes<-unique(res[,.(gene,GeneScore)][order(gene,pval)],by="gene")
  resGenes_all<-data.table(resGenes)
  print(paste("PERMUTATION of",col_to_permut,n,"times..."))
  
  
  batch1<-data.table(batch)
  for(i in 1:n){
    print(paste("permutation",i))
    
    batch1[,Group_Sex:=sample(Group_Sex)]
    
    print(" ==> limma execution")
    res<-LimmaExec(methyl_data,batch_F,formule)
    print(" ==> Gene Score calculation")
    
    
  }
  print("END PERMUTATION")
  
  print("save  pct genescore permuts >= geneScore obs")
  
  
  return(batch[,])
  
}


