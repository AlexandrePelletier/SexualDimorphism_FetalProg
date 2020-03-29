
#plot multi plot de densite
plot.multi.dens <- function(s)
 {
 junk.x = NULL
 junk.y = NULL
 for(i in 1:length(s))
 {
 junk.x = c(junk.x, density(s[[i]])$x)
 junk.y = c(junk.y, density(s[[i]])$y)
 }
 xr <- range(junk.x)
 yr <- range(junk.y)
 plot(density(s[[1]]), xlim = xr, ylim = yr, main = "",)
 for(i in 1:length(s))
 {
 lines(density(s[[i]]), xlim = xr, ylim = yr, col = i, lwd=1)
 }
 }


plusgrandque<-function(x,n){
  return(x>n)
}


#enregistrer dans : "C:/Users/alexh/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/FonctionsUtiles.R"

is.empty<-function(x){
  if(is.null(x)){
    return(T)
  }else if(is.character(x)){
    
    if(length(x)==0){
      return(T)
    }else if(all(is.na(x))|all(x=="NA")){
      
      return(T)} else if(all(x=="")){
        return(T)
      }else{
        return(F)
      }
  }else if(is.vector(x)){
    
    if(length(x)==0){
      return(T)
    }else if(all(is.na(x))){
      
      return(T)} else if(all(x=="")){
        return(T)
      }else{
        return(F)
      }
  }else if(is.numeric(x)){
    if(x==0){
      return(T)
    }else{
      return(F)
    }
  }else if(is.data.frame(x)){
    if(all(dim(x)==0)){
      return(T)
    }else{
      return(F)
    }
  }else{
    return(F)
  }
  
}


maj<-function(x){
  library(stringr)
  return(str_to_upper(x))
  
}

zeroIfNull<-function(x){
  if(is.null(x)){
    return(0)
  }else{
    return (x)
  }
}


DemanderDeContinuer<-function(msg="continuer ?"){
  # renvoit True si user say o
  if (readline(prompt=paste(msg,"(o/n)")) == "o"){
    return (T)
  }else {
    return (F)
  }
}

choisirParmiListe<-function(vecteur,indices=F,AfficheNum=T){
  
  if(AfficheNum){
    lignes<-paste(paste(1:length(vecteur),vecteur,sep = ") "),collapse = "\n" )
  }else{
    lignes<-paste(vecteur,collapse = "\n" )
  }
  choix<-readline(prompt=cat(paste("choisir les elements ( 1-3 == 1,2,3 ; 'a'=tous)",lignes,sep = "\n"),sep="\n"))
  library(stringr)
  if(choix=='a'){
    choix<-1:length(vecteur)
  }
  else if (str_detect(choix,"-")){
    choix<-as.numeric(strsplit(choix,'-')[[1]])
    choix<-choix[1]:choix[2]
  }else {
    choix<-as.numeric(strsplit(choix,',')[[1]])
  }
  
  
  if(indices){
    return(choix)
  }else{
    return(vecteur[choix])
  }
  
  
}

Afficher_vecteur_lignes<-function(lignes){
  if(is.data.frame(lignes)){
    temp<-c()
    for (variable in colnames(lignes)){
      temp<-c(temp,paste(variable, paste(lignes[,variable],sep = ", "),sep = " : "))
    }
    
    lignes<-temp
  }
  cat(paste(lignes,collapse = "\n"))
}

tradGenesID<-function(genes,from="ENSEMBL",to="SYMBOL",return="vecteur"){
  if(from=="ENSEMBL" & to == "SYMBOL"){
    genesID<-read.csv2("ENSEMBL_ID_to_SYMBOL.csv",row.names = 1)

  }
  if(from=="SYMBOL" & to == "ENSEMBL"){
    genesID<-read.csv2("ENSEMBL_ID_to_SYMBOL.csv",row.names = 2)
    
    
  }
  
  if(return == "vecteur"){
    return(genesID[genes,to])
  }else{
    return(data.frame(rownames=rownames(genesID)[genes],genesID[genes,]))
  }
    
    
}
