#2020-06-30
#deter Sex-Spe Methhylated Genes by Permutation

#deter genescore cutoff x / deter EpigenAffGene (EAG) : 
# 1000 permuts of all samples CTRL / LGA => LIMMA => pval and FC => df GeneScore => save pct genescore permuts >= geneScore obs
#=> gene score cutoff = Accept Gene in EAF if appeared < x%  => genesF and genesM spe
#make a permutFunction

permutGeneScore<-function(res_to_permut,methyl_df,cpg.regs_ref,batch,var_to_permut,n_perm=1000,seed=1234,
                             varNumToModel=c("Mat.Age"),varFacToModel=c("Group_Sex",'batch',"latino","Group_Complexity_Fac"),
                             formule= ~0 + Group_Sex  + batch  + latino + Mat.Age + Group_Complexity_Fac,
                             verbose=FALSE
                             ){
  library(data.table)
  library(stringr)
  source("scripts/calculeGeneScore.R")
  
  if(!is.list(res_to_permut)){
    print("need results in a  list named by comparison")
  }
  
  varToModel<-c(varNumToModel,varFacToModel)
  sample_F<-batch$sample[!(apply(is.na(batch[,..varToModel]),1,any))]
  
  print(paste(length(sample_F),"samples to permut"))
  
  batch<-batch[sample%in%sample_F,c("sample",..varToModel),]
  
  batch[,(varNumToModel):=lapply(.SD,as.numeric),.SDcols=varNumToModel]
  batch[,(varFacToModel):=lapply(.SD,as.factor),.SDcols=varFacToModel]
  
  compas_df<-data.frame(compa=unlist(lapply(strsplit(names(res_to_permut),"-"),function(x)paste(paste0(var_to_permut,x),collapse = "-"))),
                        abbrev_compa=names(res_to_permut))
  print(compas_df)
  

  res_all<-data.table(gene=sort(unique(res_to_permut[[1]]$gene)))
 
  
  batch_sim<-copy(batch)
  print(paste("set seed to",seed))
  set.seed(seed)
  
  print(paste("permutation de",var_to_permut,n_perm,"fois..."))
  for(i in 1:n_perm){
    print(paste(i,"/",n_perm))
    
    batch_sim[,Group_Sex:=sample(Group_Sex)]
    res_list<-RunMethAnalysis(methyl_df = methyl_df,
                              batch_F = batch_sim,
                              formule = formule,
                              compas_df =compas_df,
                              cpg.regs_ref = cpg.regs_ref,
                              sumToGene = T,
                              verbose=verbose)
    
  
      res_list<-lapply(names(res_list),function(compa){
        res<-res_list[[compa]]
        res<-res[order(gene)][,.(gene,GeneScore)]
        col<-paste0(str_sub(compa,str_length(compa)),i)
        return(res[,(col):=as.integer(GeneScore)][,-"GeneScore"])
      })
        
      res_list<-Reduce(merge,res_list)
      res_all<-merge(res_all,res_list,by="gene")
      
  }
  
  res_to_permut<-lapply(1:length(res_to_permut),function(i){

      res<-unique(res_to_permut[[i]],by="gene")
      compa<-names(res_to_permut)[i]
      cols<-paste0(str_sub(compa,str_length(compa)),1:n_perm)
      col<-paste0("pval",n_perm,"perm")
      res<-res[order(gene)][,(col):=rowSums(..res_all[order(gene)][,.SD,.SD=cols]>GeneScore)/(..n_perm)]
      cols<-c("gene",col)
      return(merge(res_to_permut[[i]],res[,.SD,.SD=cols],by="gene"))
  })
  names(res_to_permut)<-compas_df$abbrev_compa
  return(res_to_permut)
}

RunMethAnalysis<-function(methyl_df,batch_F,formule,compas_df,cpg.regs_ref,sumToGene=FALSE,verbose=TRUE){
  library(limma)
  design<-model.matrix(formule,data = batch_F)
  
  fit <- lmFit(methyl_df[,batch_F$sample], design)  
  
  cont.matrix <- makeContrasts(contrasts = compas_df$compa,
                               levels=design)
  colnames(cont.matrix)<-compas_df$abbrev_compa
  
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2) 
  if(nrow(compas_df)>1){
    
    res_list<-lapply(compas_df$abbrev_compa,function(compa){
      res<-topTable(fit2,coef=compa,n =Inf)
      res<-CalcGeneScore(res,cpg.regs_ref,sumToGene=sumToGene,verbose = verbose)
      return(res[order(-GeneScore)])
    })
    names(res_list)<-compas_df$abbrev_compa
    return(res_list)
    
  }else{
    res<-topTable(fit2,coef=compas_df$abbrev_compa,n =Inf)
    res<-CalcGeneScore(res,cpg.regs_ref,sumToGene=sumToGene,verbose=verbose)
    return(res[order(-GeneScore)])
  }
  
  
  
}




