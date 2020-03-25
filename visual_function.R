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
