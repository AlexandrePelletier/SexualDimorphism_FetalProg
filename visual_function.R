plotPCA<-function(pca,PCx,PCy,colorBatch,batch=NULL,showSampleIDs=F){
  library(stringr)
  if(is.null(batch)){
    batch<-read.csv2("../../ref/20-02-17_batch_119samples_IUGRLGACTRL_noNA_withHpaC.csv",header=T,row.names = 1)
  }
  pc<-pca$x
  namesFac<-levels(as.factor(batch[,colorBatch]))
  if(colorBatch=="Group"){
    namesFac<-c("CTRL","IUGR","LGA")
  }
  if(colorBatch=="Sex"){
    namesFac<-c("M","F")
  }
  if(colorBatch=="Group_Sex"){
    namesFac<-c("female_CTRL","male_CTRL","female_IUGR","male_IUGR","female_LGA","male_LGA")
  }
  
  pct.varX<-pca$sdev[PCx]^2/sum(pca$sdev^2)*100
  pct.varY<-pca$sdev[PCy]^2/sum(pca$sdev^2)*100
  plot(pc[,PCx],
       pc[,PCy],
       col=batch[rownames(pc),colorBatch]+1,
       pch=19,
       xlab = paste0("PC",PCx," (",round(pct.varX,1),"% of variance)"),
       ylab = paste0("PC",PCy," (",round(pct.varY,1),"% of variance)"))
  
  if(showSampleIDs){
    textxy(pc[,PCx],pc[,PCy],rownames(pc),offset=-.7)
  }
  colors<-as.numeric(levels(as.factor(na.omit(batch[,colorBatch]+1 ))))
  legend("topleft", legend = namesFac ,col =colors,pch=19, title=str_sub(colorBatch,1,7))
  
}

plotPCVarExplain<-function(pca,rngPCs,lineSeuilPct=1,returnPCSeuils=T){
  source("scripts/utils.R")
  pct.varPCs<-pctPC(pca,rngPCs)
  barplot(pct.varPCs[rngPCs],names.arg = rngPCs,ylab = "Percent of variance explained")
  abline(h=lineSeuilPct)
  return(as.numeric(names(pct.varPCs)[pct.varPCs>lineSeuilPct]))
  
}

plotCovarPCs<-function(pca,rngPCs,batch,var_num,var_fac,exclude=NULL,seuilP=0.1){
  source("scripts/utils.R")
  pc<-data.frame(pca$x)
  batch_num=batch[rownames(pc),var_num]
  batch_fac=batch[rownames(pc),var_fac]
  batch_num=t(batch_num)
  batch_fac=t(batch_fac)
  split_num=split(batch_num,rownames(batch_num))
  split_fac=split(batch_fac,rownames(batch_fac))
  pv_num=lapply(split_num,function(x){
    FAC1.p<-rep(0,length(rngPCs))
    for (i in rngPCs){
      FAC1<-as.numeric(x)
      FAC1<-lm(pc[,i]~FAC1)
      FAC1.p[i]<-anova(FAC1)$Pr[1]
    }
    return(FAC1.p)})
  
  pv_fac=lapply(split_fac,function(x){
    FAC1.p<-rep(0,length(rngPCs))
    for (i in rngPCs){
      FAC1<-as.factor(x)
      FAC1<-lm(pc[,i]~FAC1)
      FAC1.p[i]<- anova(FAC1)$Pr[1]
    }
    return(FAC1.p)})
  
  
  pvals.num<-do.call(rbind,pv_num)
  pvals.fac<-do.call(rbind,pv_fac)
  final_pv=rbind(pvals.num,pvals.fac)
  pv2=data.matrix(final_pv)
  pv2[which(pv2>seuilP)]<-1 ####here I basicaly put them to 1 if less than 0.05
  logpvals.raw<--log10(pv2)
  
  pct.varPCs<-pctPC(pca,rngPCs)
  vars<-rownames(logpvals.raw)[!(rownames(logpvals.raw)%in%exclude)]
  pheatmap(logpvals.raw[vars,rngPCs],cluster_rows = F,cluster_cols = F,
           labels_col= paste0("PC",rngPCs,"(",round(pct.varPCs[as.character(rngPCs)],0),"%)"),
           display_numbers = T,
           color = colorRampPalette(c("white", "red"))(13), breaks = c(40,20,10:1, 0.5,0.1))
  
  return(pv2)
  
  
}
