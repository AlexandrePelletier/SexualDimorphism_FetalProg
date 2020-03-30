
source("scripts/FonctionsUtiles.R")

#fonctions ajouter/suppr gene dans tableaux de marqueurs voir abbrev_pop.csv pour les noms exactes
suppr_marqueurs<-function(genes, pops="all",importance=2){
  if (importance==2){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursRestreints.csv",colClasses = "character",row.names = 1)
  }else if (importance==1){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTresRestreints.csv",colClasses = "character",row.names = 1)
  } else {
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTous.csv",colClasses = "character",row.names = 1)
  }
  
  for (gene in maj(genes)){
    if (pops=="all"){
      marqueurs<-marqueurs[rownames(marqueurs)!=gene,]
      
    }else{
      i<-1
      for (pop in GetPopDe(gene)){
        if (!(pop %in% pops )){
          marqueurs[gene,paste("pop",i,sep="")]<-pop
          i<-i+1
        }
        
      }
      while(i<ncol(marqueurs)-1){
        marqueurs[gene,paste("pop",i,sep="")]<-""
        i<-i+1
      }
      
      
    }
      
    }
  if (importance==2){
    write.csv2(marqueurs, file = "D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursRestreints.csv",na = "")
  }else if (importance==1){
    write.csv2(marqueurs, file = "D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTresRestreints.csv",na = "")
  } else {
    write.csv2(marqueurs, file = "D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueurs.csv",na = "")
  }
  
  return(marqueurs)
}


add_new_marqueurs<-function(genes, new_pops, ref=NULL, importance=2, retourner=T){
  
  if (importance==2){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursRestreints.csv",colClasses = "character",row.names = 1)
  }else if (importance==1){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTresRestreints.csv",colClasses = "character",row.names = 1)
  } else {
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTous.csv",colClasses = "character",row.names = 1)
  }
  
  pops_abbrev<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/abbrev_pop.csv",colClasses="character")
  
  for (gene in maj(genes)){
    pops<-c()
    if (gene %in% rownames(marqueurs)){
      
      #on recupere la liste des pops
      lastColPop<-tail(names(marqueurs),2)[1]
      numLastColPop<-tail(strsplit(lastColPop,"")[[1]],1)
      nb_pop_marked<-sum(!sapply(marqueurs[gene,paste("pop",1:numLastColPop,sep="")],is.empty))
      for (i in 1:nb_pop_marked){
        pops <- c(pops,as.character(marqueurs[gene,paste("pop",i,sep="")]))
        
      }
      print(paste("gene :",gene,"deja connu comme marqueur de :", paste(pops, collapse = ", ")))
      
    }else{
      nb_pop_marked<-0
      
    }
    
    
      
      for (new_pop in new_pops){
        if (!(new_pop %in% pops)){
          if(new_pop %in% pops_abbrev$abbrev){
            nb_pop_marked<-nb_pop_marked+1
            marqueurs[gene,paste("pop",nb_pop_marked,sep="")]<-new_pop
          }else {
            if (readline(prompt=paste("'",new_pop,"'","n'est pas repertori?, voulez vous l'ajouter ? (o/n)"  )) =="o"){
              
              pops_abbrev<-addAbbrevPop(new_pop,readline(prompt="ecrire la description : " ))
              nb_pop_marked<-nb_pop_marked+1
              marqueurs[gene,paste("pop",nb_pop_marked,sep="")]<-new_pop
            }
          }
          
        }
      }
      
    if(!is.null(ref)){
      if(is.empty(marqueurs[gene,"ref"])){
        marqueurs[gene,"ref"]<-ref
      }else{
        if(!(ref %in% strsplit(marqueurs[gene,"ref"],";")[[1]]) ){
          marqueurs[gene,"ref"]<-paste(marqueurs[gene,"ref"], ref,sep=";")
        }
        
      }
      
    }
    
  } 
  if (importance==2){
    write.csv2(marqueurs, file = "D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursRestreints.csv",na = "")
  }else if (importance==1){
    write.csv2(marqueurs, file = "D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTresRestreints.csv",na = "")
  } else {
    write.csv2(marqueurs, file = "D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueurs.csv",na = "")
  }
  
  
  
  if (retourner == T){
    return(marqueurs) 
  }
  
}

addAbbrevPop<-function(newpop, description){
  pops_abbrev<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/abbrev_pop.csv",colClasses="character")
  newligne<-nrow(pops_abbrev)+1
  
  pops_abbrev[newligne,"abbrev"]<-newpop
  pops_abbrev[newligne,"population"]<-description
  write.csv2(pops_abbrev,file = "D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/abbrev_pop.csv",row.names = F)
  return(pops_abbrev)
  
}

GetMarqueursDe<-function(pop, strict=TRUE, importance=2){
  if (importance==2){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursRestreints.csv",colClasses = "character",row.names = 1)
  }else if (importance==1){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTresRestreints.csv",colClasses = "character",row.names = 1)
  } else {
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTous.csv",colClasses = "character",row.names = 1)
  }
  
  if (strict){
    genes<-c()
    for (gene in rownames(marqueurs)){
      if (any(marqueurs[gene,]==pop)){
        
        genes<-c(genes, gene)
        
        return(genes)

      }

    }
    
  }else {
    library(stringr)
    genes<-list()
    for (gene in rownames(marqueurs)){
      for (popnum in colnames(marqueurs)){
        
        if (str_detect(marqueurs[gene,popnum],pop)){
          popu<-marqueurs[gene,popnum]
          
          if (is.null(genes[[popu]])){
            
            
            genes[[popu]]<-gene
          }else{
            
            genes[[popu]]<-c(genes[[popu]],gene)
            
          }
          
          
        }
      }
    }
    return(genes)
    
  }
  
}

GetAllPopDe<-function(marqueurs,strict=F){
  liste<-list()
  allmarqueurs<-c()
  if(!strict){
    dfMarqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursRestreints.csv",colClasses = "character",row.names = 1)
    library(stringr)
    for(marqueur in maj(marqueurs)){
      marqueursTrouves<-rownames(dfMarqueurs)[str_detect(rownames(dfMarqueurs),marqueur)]
      allmarqueurs<-c(allmarqueurs,marqueursTrouves)
    }
    
    marqueurs<-allmarqueurs
  }
  
  for (marqueur in maj(marqueurs)){
    liste[[marqueur]]<-GetPopDe(marqueur)
  }
  return(liste)
}

GetPopDe<-function(marqueur,importance=2){
  library(stringr)
  marqueur<-maj(marqueur)
  if (importance==2){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursRestreints.csv",colClasses = "character",row.names = 1)
  }else if (importance==1){
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTresRestreints.csv",colClasses = "character",row.names = 1)
  } else {
    marqueurs<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/marqueursTous.csv",colClasses = "character",row.names = 1)
  }
  
  
  if (!all(is.na(marqueurs[marqueur,]))){
    pops<-c()
    for (pop in marqueurs[marqueur,]){
      
      if (str_length(pop)>1){
        pops<-c(pops,pop)
        
      }
      
      
    }
    return(pops)
    
  }else{
    return(NULL)
  }
}


AfficherAbbrevPops<-function(pop="all",retourner=TRUE){
  dfAbbrev<-read.csv2("D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/abbrev_pop.csv",colClasses="character")
  
  if (pop=='all'){
    
    print(dfAbbrev)
    if (retourner){
      return(dfAbbrev)
    }
    
    
  }else {
    library(stringr)
    subsetdf<-data.frame()
    for (description in dfAbbrev$population){
      
      if(str_detect(description,pop)){
        ligne<-dfAbbrev[which(dfAbbrev$population==description),]
        print(ligne)
        subsetdf<-merge(subsetdf,ligne,all=T)
      }
    }
    
  }
}
