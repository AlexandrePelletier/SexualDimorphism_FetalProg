## !!! encore ? optimiser :
## - genes infos : 
#
##        - description entrez et protein atlas blood (de type https://www.proteinatlas.org/ENSG00000132465-JCHAIN/blood),
##http://www.proteinatlas.org/api/search_download.php?search=JCHAIN&format=json&columns=gs,rnabcs,rnabcd,rnabls,rnabld&compress=no
##        -lien single cell db si on trouve
##        - ecrire dependencies : KEGG rest, reticulate..
#         _ ajouter try pour eviter que ça plante a chaque fois que le serveur repond pas
#-temps de calcul : faire le min de requete biomartr.


#         

source("scripts/FonctionsUtiles.R")
library(stringr)
options(stringsAsFactors = F)
library(biomartr)


annotation_auto_genes<-function(genes,annotations="all",df=NULL,save=T,indep=F,onlydejasave=F, reset=F, resetFull=F, comment=NULL){
  ### fonction qui annote* les gènes et les enregistre dans le fichier genes_infos.csv pour qu'il puisse être reutilisé ensuite
  
    genes<-unique(maj(genes))
    if (annotations != "all"){
      #! a faire
    }
  
  ## si un commentaire a ajouter sur les genes demander :
  
  
  #on recupere les genes deja annoter dans le fichier genes_infos.csv :
    
    if (is.null(df)&!indep){
      print("recuperation des donnes deja recuperees...")
      genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1,colClasses = "character")
    }else if (!is.null(df)){
      genes_infos<-df
    }
  
  if (indep ){
    
    if(!is.null(df)){
      genes_a_annoter<-genes[!(genes %in% rownames(genes_infos))]
    }else{
      genes_a_annoter<-genes
      genes_infos<-data.frame(row.names = genes_a_annoter)
    }
    
    
  }else if (!reset & !resetFull){
    # si on veut juste annoter les nouveaux genes :
    genes_a_annoter<-genes[!(genes %in% rownames(genes_infos))]
    
    
  }else if (reset) {
    #reannoter tout les genes demandes:
    genes_a_annoter<-genes
    
  }else{
    #reanntoter tout les genes du tableau :
    genes_a_annoter<-union(genes,rownames(genes_infos))
    genes_infos<-data.frame(row.names = genes_a_annoter)
  }
  
  if(!is.null(comment)){
    genes_infos<-add_comment(genes,comment,df = genes_infos,saveOnFile = F)
    
  }
  print("reussit !")
  if (length(genes_a_annoter)>0 & !onlydejasave){
    print("recherche des descriptions des nouveaux genes...")
    try(genes_infos<-find_description(genes_a_annoter,df=genes_infos))
    print("reussit!")
    print("recherche de la fonction des genes...")
    try(genes_infos<-find_fonction(genes_a_annoter,df=genes_infos,save=T) )
    print("reussit!")
    print("recherche des process biologiques associes...")
    try(genes_infos<-find_GoBioProcess(genes_a_annoter,df=genes_infos,save=T))
    print("reussit!")
    print("recherche des pathways associes...")
    try(genes_infos<-find_KEGGPathway(genes_a_annoter,df=genes_infos,save=T))
    print("reussit!")

    print("recherche si c'est un facteur de transcription...")
    try(genes_infos<-findIfTF(genes_a_annoter,genes_infos))
    print("reussit!")
    
    print("recherche si ils sont lies a la differentiation du tissu...")
    try(genes_infos<-findCbIn(genes = genes_a_annoter,listeKeywords = bioKeywords,genes_infos = genes_infos))
    print("reussit!")
  }
    
    ##on actualise la col "MarqueurPops"
    print("recherche si ils sont des marqueurs de sous population deja connu...")
    genes_infos<-findMarqueursPops(rownames(genes_infos),df=genes_infos)
    print("reussit!")
    if(save){
      print("sauvegarde des informations...")
      write.csv2(genes_infos,file="ref/genes_infos.csv",na = "")
      print("reussit !")
    }
    
    
  return(genes_infos[genes,])
}


find_description<-function(genes,df=NULL,searchInLocal=T){
  if(is.null(df)){
    
    df<-data.frame(row.names = genes)
    
  }
  if(searchInLocal){
    genes_desc<-read.csv2("ref/AllGeneSymbolsAndDescriptions.csv",row.names = 1)
    genes_desc<-data.frame(row.names =genes,Description= genes_desc[genes,"Approved.name"])
    
  }else{
    print("a coder")
  }

  
  
  if (is.null(df)){
    return(genes_desc)
  }
  else{
    df[genes,"Description"]<-genes_desc$Description
    return (df)
  }
}


find_fonction<-function(genes,df=NULL,save=F,searchInLocal=T){
  if(is.null(df)){
    
    df<-data.frame(row.names = genes)
    
  }
  if(searchInLocal){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1,na.strings = "")
    genesDejaAnnot<-genes[genes%in%rownames(genes_infos)]
    genesDejaAnnot<-genesDejaAnnot[(!is.na(genesInfos[genesDejaAnnot,"Fonction"]))&(genesInfos[genesDejaAnnot,"Fonction"]!="")]
    if(length(genesDejaAnnot>0)){
      print(paste("genes deja annote :", paste(genesDejaAnnot,collapse = ", ")))
      df[genesDejaAnnot,"Fonction"]<-genes_infos[genesDejaAnnot,"Fonction"]
      genesAAnnot<-genes[!(genes%in%genesDejaAnnot)]
    }else{
      genesAAnnot<-genes
    }
    
  }else{
    genesAAnnot<-genes
  }
  
  
  #I) trouver fonction description de UNIPROT
  #1) trouver uniprotID des genes
  genesUniprots<-uniprotID(genesAAnnot)
  
  #2) chercher la description de la fonction par requete url grace une fonction python
  library(reticulate)
  #try(use_python("C:/ProgramData/Anaconda3/",required = T))
  
  source_python("scripts/GetUniprotFunction.py")
  if (length(genesUniprots$uniprotswissprot)>0){
    
    for (uniprotID in genesUniprots$uniprotswissprot){
      
      
      fonction<-getFonctionUniProt(uniprotID)
      
      genesUniprots[ which(genesUniprots$uniprotswissprot==uniprotID),"Fonction"]<-fonction
      
    }
  }
  
  if (is.null(df)){
    return(genesUniprots)
  }
  else{
    #df[genesUniprots$hgnc_symbol,"uniprotID"]<-genesUniprots$uniprotswissprot
    df[genesUniprots$hgnc_symbol,"Fonction"]<-genesUniprots$Fonction
    
    if(save){
      print(paste("sauvegarde de ",length(genesAAnnot),"nouveau genes annotes pour leur fonction dans ref/genes_infos.csv"))
      genes_infos[genesAAnnot,"uniprotID"]<-genesUniprots[genesAAnnot,"uniprotswissprot"] 
      genes_infos[genesAAnnot,"Fonction"]<-df[genesAAnnot,"Fonction"]
      write.csv2(genes_infos,file="D:/Profils/apelletier/Dropbox/LAB_Delahaye/Alexandre_SC/analyses/ref/genes_infos.csv",na = "")
    }
    return (df)
  }
}


findIfTF<-function(genes,genes_infos=NULL,save=F,large=F){
  
  
  if(is.null(genes_infos)){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1)
    
  }
  
  for (gene in genes){
    if(!is.empty(genes_infos[gene,"Fonction"])){
      
      fonction<-str_to_lower(genes_infos[gene,"Fonction"])
      if(!is.empty(genes_infos[gene,"Description"])){
        fonction<-paste(fonction,str_to_lower(genes_infos[gene,"Description"]))
        
      }
      
      if(large){
        if(str_detect(fonction,"^.{0,70}transcription.{0,2} regulator|^.{0,15}transcriptional.{0,70} activator|^.{0,15}transcriptional.{0,70} (repressor|regulator)")){
          genes_infos[gene,"TF"]<-T
        }else{
          genes_infos[gene,"TF"]<-F
        }
      }else{
        if(str_detect(fonction,"^.{0,70}transcription.{0,2} regulator")){
          genes_infos[gene,"TF"]<-T
        }else{
          genes_infos[gene,"TF"]<-F
        }
        if(str_detect(fonction,"^.{0,15}transcriptional.{0,70} activator")){
          genes_infos[gene,"TranscriptionalActivator"]<-T
        }else{
          genes_infos[gene,"TranscriptionalActivator"]<-F
        }
        if(str_detect(fonction,"^.{0,15}transcriptional.{0,70} (repressor|regulator)")){
          genes_infos[gene,"TranscriptionalRepressor"]<-T
        }else{
          genes_infos[gene,"TranscriptionalRepressor"]<-F
        }
        
      }
      
      
      
      
    }
    
  }
  
  
  if (save){
    write.csv2(genes_infos,file="ref/genes_infos.csv",na = "")
  }
  
  return(genes_infos)
}





findCbIn<-function(genes,listeKeywords=NULL,genes_infos=NULL,save=F){
  #Ly NK Bcell Tcell DC Gran Neu Mas Eo Bas Mac Ery Mk MyeloP LyP HematoP Immunity Metabo Stem CycleCell
  if(is.null(listeKeywords)){
    source("ref/bioKeywords.R")
    listeKeywords<-bioKeywords
    
  }
  
  if(is.null(genes_infos)){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1)
    
  }
  
  for (gene in genes){
    
    infos<-c()
    if(any(names(genes_infos)=="Fonction")){
      infos<-c(infos,strsplit(str_to_lower(genes_infos[gene,"Fonction"]),"\\.")[[1]])
    }
    
    
    if(any(names(genes_infos)=="Description")){
      if(!is.empty(genes_infos[gene,"Description"])){
        infos<-c(infos,str_to_lower(genes_infos[gene,"Description"]),str_to_lower(genes_infos[gene,"Description"]))
      }
      
    }
    
    if(any(names(genes_infos)=="BioProcess")){
      infos<-c(infos,strsplit(str_to_lower(genes_infos[gene,"BioProcess"]),";")[[1]])
    }
    if(any(names(genes_infos)=="Pathway")){
      infos<-c(infos,strsplit(str_to_lower(genes_infos[gene,"Pathway"]),";")[[1]])
    }

    
    
    for(concept in names(listeKeywords)){
      
      genes_infos[gene,concept]<-sum(str_detect(infos, listeKeywords[[concept]]))
      
    }
    
  }
  if (save){
    write.csv2(genes_infos,file="ref/genes_infos.csv",na = "")
  }
  
  return(genes_infos)
}



findCbStemRelated<-function(genes,genes_infos=NULL,save=F){ #! a finir
  StemKeywords<-"stem|..."
  
  options(stringsAsFactors = F)
  library(stringr)
  if(is.null(genes_infos)){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1)
    
  }
  
  for (gene in genes){
    fonction<-strsplit(str_to_lower(genes_infos[gene,"Fonction"]),".")[[1]]
    bioprocess<-strsplit(str_to_lower(genes_infos[gene,"BioProcess"]),";")[[1]]
    pathway<-strsplit(str_to_lower(genes_infos[gene,"Pathway"]),";")[[1]]
    infos<-c(fonction,bioprocess,pathway)
    
    genes_infos[gene,"HematopoRelated"]<-sum(str_detect(infos, stemKeywords))
  }
  if (save){
    write.csv2(genes_infos,file="ref/genes_infos.csv",na = "")
  }
  
  return(genes_infos)
}

uniprotID<-function(genes){
  
  result_BM <- biomartr::biomart( genes      = genes, # genes that we wanted info
                                  mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                  dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                  attributes = "uniprotswissprot", # attributes were selected with biomartr::getAttributes()
                                  filters    = "hgnc_symbol",) # specify what ID type was stored in the fasta file retrieved with biomartr::getGenome()
  
  result_BM<-result_BM[which(result_BM$uniprotswissprot!=""),]
  return(result_BM[!duplicated(result_BM$hgnc_symbol),])
  
}


find_GoBioProcess<-function(genes,df=NULL,save=F,searchInLocal=T){
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  if(is.null(df)){
    
    df<-data.frame(row.names = genes)
    
  }
  if(searchInLocal){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1,na.strings = "")
    genesDejaAnnot<-genes[genes%in%rownames(genes_infos)[!is.na(genes_infos$BioProcess)]]
    print(paste("genes deja annote :", paste(genesDejaAnnot,collapse = ", ")))
    df[genesDejaAnnot,"BioProcess"]<-genes_infos[genesDejaAnnot,"BioProcess"]
    genesAAnnot<-genes[!(genes%in%genesDejaAnnot)]
  }else{
    genesAAnnot<-genes
  }
  
  for (gene in genesAAnnot){
    if (gene %in% keys(org.Hs.eg.db,keytype = "SYMBOL")){
      
      gene_infos<- AnnotationDbi::select(org.Hs.eg.db, gene ,columns =  c("ONTOLOGY","GO"), keytype = 'SYMBOL')
      
      GO_IDs<-subset(gene_infos, ONTOLOGY=="BP" &EVIDENCE != "IEA", select=GO)$GO
      GO_infos<- AnnotationDbi::select(GO.db,GO_IDs,c("TERM"),"GOID")$TERM
      
      df[gene,"BioProcess"]<-paste(unique(GO_infos),collapse = ";")
      
    }
  }
  if(save&length(genesAAnnot)>0){
    print(paste("sauvegarde de ",length(genesAAnnot),"nouveau genes annotes pour les BioProcess dans ref/genes_infos.csv"))
    genes_infos[genesAAnnot,"BioProcess"]<-df[genesAAnnot,"BioProcess"]
    write.csv2(genes_infos,file="ref/genes_infos.csv",na = "")
  }
  return(df)
  
}

cleanBioProcess<-function(genes,genes_infos=NULL,nbMotsSignifMax=2,nbTermesSimilaireMax=1,save=F){
  #! ajouter 1ere etape = enlever termes non pertinent : faire fct pour calculer score de pertinence des termes 
  #dependant de mots choisis (relatif a hematopoKeyWord, StemKeyword,MetaboKeyword,precision(negative regulation
  #of mieux que regulation of)),et de la raret? des mots (par rapport a l'ensemble des termes DO)
  options(stringsAsFactors = F)
  library(stringr)
  if(is.null(genes_infos)){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1)
    
  }
  nbgeneSansChgt<-0
  for (gene in genes){
    
    AllTerms<-unique(strsplit(genes_infos[gene,"BioProcess"],";")[[1]])
    
    keep<-rep(T,length(AllTerms))
    listeTermsDecompo<-strsplit(str_to_lower(AllTerms)," |-")
    termesAExclure<-c()
    i<-0
    for (motsDuTerme in listeTermsDecompo){
      i<-i+1
      terme<-AllTerms[i]
      if(!(terme%in%termesAExclure)){
        
        vecNbMotsCommuns<-sapply(listeTermsDecompo,nbDeMotsCommun,motsDuTerme,exceptions=c("regulation","protein","activation","pathway","process","response"))
        
        if (sum(vecNbMotsCommuns>nbMotsSignifMax)>nbTermesSimilaireMax){
          termesSimilaires<-AllTerms[vecNbMotsCommuns>nbMotsSignifMax]
          
          if(length(termesSimilaires)>3){
            termesGardes<-c(termesSimilaires[c(which.min(str_length(termesSimilaires)),which.max(str_length(termesSimilaires)))])
          }else{
            termesGardes<-termesSimilaires[which.max(str_length(termesSimilaires))]
          }
          
          keep<-keep&((vecNbMotsCommuns<=nbMotsSignifMax)|AllTerms%in%termesGardes)
          
          termesAExclure<-AllTerms[(!keep)|AllTerms%in%termesGardes]
          
        }
      }
      
    }
    termesAvant<-length(AllTerms)
    AllTerms<-AllTerms[keep]
    termesApres<-length(AllTerms)
    genes_infos[gene,"BioProcess"]<-paste(AllTerms,collapse = ";")
    if(termesApres<termesAvant){
      
    }else{
      nbgeneSansChgt<-nbgeneSansChgt+1
    }
    
  }
  print(paste(nbgeneSansChgt,"genes n'ont pas ?t? affect?"))
  
  return(genes_infos)
  
  
}

nbDeMotsCommun<-function(vecDeMots,vecDeMots2,exceptions=NULL,plusGrandQue=3){
  library(stringr)
  vecDeMots<-vecDeMots[str_length(vecDeMots)>plusGrandQue]
  vecDeMots2<-vecDeMots2[str_length(vecDeMots2)>plusGrandQue]
  
  if(!is.null(exceptions)){
    for(exception in exceptions){
      vecDeMots<-vecDeMots[vecDeMots!=exception]
      vecDeMots2<-vecDeMots2[vecDeMots2!=exception]
    }
    
  }
  return(sum(vecDeMots%in%vecDeMots2))
}



find_KEGGPathway<-function(genes,df=NULL,save=F,searchInLocal=T){
  
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(KEGGREST)
  if(is.null(df)){
    
    df<-data.frame(row.names = genes)
    
  }
  if(searchInLocal){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1,na.strings = "")
    genesDejaAnnot<-genes[genes%in%rownames(genes_infos)[!is.na(genes_infos$Pathway)]]
    print(paste("genes deja annote :", paste(genesDejaAnnot,collapse = ", ")))
    df[genesDejaAnnot,"Pathway"]<-genes_infos[genesDejaAnnot,"Pathway"]
    genesAAnnot<-genes[!(genes%in%genesDejaAnnot)]
  }else{
    genesAAnnot<-genes
  }
  
  for (gene in genesAAnnot){

    df[gene,"Pathway"]<- ""
    KEGGs_infos<-c()
    KEGG_IDs<-data.frame()
    if (gene %in% keys(org.Hs.eg.db,keytype = "SYMBOL")){
      KEGG_IDs<-AnnotationDbi::select(org.Hs.eg.db, gene,columns =  c("PATH"), keytype = 'SYMBOL') 
      if(!is.null(KEGG_IDs$PATH)){
        KEGG_IDs<-KEGG_IDs$PATH
      }else{
        KEGG_IDs<-NULL
      }
      
    }else {
      KEGG_IDs<-NULL
    }
    if(!is.null(KEGG_IDs)){
      
      for (KEGGID in KEGG_IDs){
          
          try(KEGG_info<-keggGet(paste("path:hsa",KEGGID,sep="")),silent=T)
          try(KEGGs_infos<-c(KEGGs_infos,strsplit(KEGG_info[[1]]$NAME,"-")[[1]][1]),silent = T)
        }
        df[gene,"Pathway"]<- paste(KEGGs_infos, collapse =";") 
      }
    
  }
  if(save){
    print(paste("sauvegarde de ",length(genesAAnnot),"nouveau genes annotes pour les pathways dans ref/genes_infos.csv"))
    genes_infos[genesAAnnot,"Pathway"]<-df[genesAAnnot,"Pathway"]
    write.csv2(genes_infos,file="ref/genes_infos.csv",na = "")
  }
  return(df)
}

  


add_comment<-function(genes, comment, df=NULL, saveOnFile=T,retourner=T){
  options(stringsAsFactors = F)
  if(is.null(df)){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1)
    
  }else {
    genes_infos<-df
  }
  
  for (gene in genes){
    if(is.empty(genes_infos[gene,"commentaire"])){
        genes_infos[gene,"commentaire"]<-comment
      }else{
        if(!(comment %in% strsplit(genes_infos[gene,"commentaire"],", ")[[1]])){
          genes_infos[gene,"commentaire"]<-paste(genes_infos[gene,"commentaire"],comment,sep=", ")
        }
      }
  }
  if (saveOnFile){
    write.csv2(genes_infos, file = "ref/genes_infos.csv",na = "")
  }
  if(retourner){
    return(genes_infos)
  }
}


suppr_comment<-function(genes, comment, df=NULL, saveOnFile=T,retourner=T){
  
  
  if(is.null(df)){
    genes_infos<-read.csv2("ref/genes_infos.csv",row.names = 1, na.strings = NULL,colClasses = "character")
    
  }else {
    genes_infos<-df
  }
  
  if(genes=="all"){
    genes<-rownames(genes_infos)
  }
  
  for (gene in genes){
    if(!is.empty(genes_infos[gene,"commentaire"])){
      
      if((comment %in% strsplit(genes_infos[gene,"commentaire"],", ")[[1]])){
        comments<-strsplit(genes_infos[gene,"commentaire"],", ")[[1]]
        commentairesGardes<-c()
        for (commenta in comments){
          if (!str_detect(comment, commenta)){
            commentairesGardes<-c(commentairesGardes,comment)
          }
        }
        genes_infos[gene,"commentaire"]<-paste(commentairesGardes,collapse = ", ")
      }
      
    }
    
  }
  if (saveOnFile){
    write.csv2(genes_infos, file = "ref/genes_infos.csv",na = "")
  }
  if(retourner){
    return(genes_infos)
  }
}


findMarqueursPops<-function(genes,df=NULL){
  source("scripts/ManagePopMarqueurs.R")
  if(is.null(df)){
    df<-data.frame()
  }
  
  for (gene in genes){
    pops<-GetPopDe(gene)
    if(!is.null(pops)){
      df[gene,"MarqueurPops"]<-paste(pops,collapse = ";")
    }else{
      df[gene,"MarqueurPops"]<-""
    }
      
    }
  
  return(df)
}



  
getAnnot<-function(gene,annotationType,df=NULL){
  
  if(is.null(df)){
    df<-read.csv2("ref/genes_infos.csv",row.names = 1,colClasses = "character")
  }
  
  if (!is.empty(df[gene,annotationType])){
    
    return(strsplit(df[gene,annotationType],";")[[1]])
    }else {
    return(NULL)
  }
}

FindAnnotCommun<-function(gene1,gene2,annotationType,df=NULL){

annotsGene1<-getAnnot(gene=gene1,annotationType=annotationType,df=df )
annotsGene2<-getAnnot(gene=gene2,annotationType=annotationType,df=df)
return(intersect(annotsGene1,annotsGene2))
}

TrouveGenesAmis<-function(gene,genes,genes_infos, max=2,annotationType="all"){
  #calcule a cb de g?ne, le g?ne est reli?, d'annotation bio partag? avec suffisament de marqueur du cluster
  #(pour commencer, disons etre r?li? a un g?ne <=> 1annot partag? pour 3 annots tot,
  # et peut marquer jusqu'a 2 point : 2/nbgenesTots)
  amis<-c()
  if (annotationType=="all"){
    annotations<-colnames(genes_infos)[colnames(genes_infos)%in%c("BioProcess","Pathway","MarqueurPops")]
  }else{
    annotations<-annotation
  }

  for (annot in annotations){
    for (gene_ami in genes){
      annot_commun<-FindAnnotCommun(gene1 = gene,gene2 = gene_ami,annotationType = annotation,df=genes_infos)
      if (length(annot_commun)/length(getAnnot(gene,annotation))){
        amis<-c(amis,ami)
      }
    }



  }
  return(amis)
}

ancienTrouveGenesAmis<-function(gene,genes,genes_infos,annotationType="all"){ 
   
  
  genesAmi=c()
  
  if (annotationType=="all"){
    annotations<-colnames(genes_infos)[colnames(genes_infos)%in%c("BioProcess","Pathway","MarqueurPops")]
  }else{
    annotations<-annotationType
  }
  for (annotation in annotations){
    if (annotation=="BioProcess"){
      
      annotsMin<-0.6+length(getAnnot(gene,"BioProcess"))/5
      
    }else if(annotation=="Pathway"){
      annotsMin<-0.6+length(getAnnot(gene,"Pathway"))/5
    }else{
      annotsMin<-0.6+length(getAnnot(gene,"MarqueurPops"))/5
    }
    
    for (gene_ami in genes){
      
      annot_commun<-FindAnnotCommun(gene1 = gene,gene2 = gene_ami,annotationType = annotation,df=genes_infos)
      
      
      if( length(annot_commun) >= annotsMin){
        genesAmi<-union(genesAmi, gene_ami)
      }
    }
    
  }
  
  return(genesAmi)
}

TrouveGenesLinks<-function(gene,genes,genes_infos,annotationType="all"){
  #retourne les genes relies au gene par au moins 1 2 ou 3 annotation bio partag? avec suffisament de marqueur du cluster
  #(pour commencer, disons etre r?li? a un g?ne <=> 1annot partag? pour 3 annots tot,
  # et peut marquer jusqu'a 2 point : 2/nbgenesTots)
  amis<-c()
  if (annotationType=="all"){
    annotations<-colnames(genes_infos)[colnames(genes_infos)%in%c("BioProcess","Pathway","MarqueurPops")]
  }else{
    annotations<-annotationType
  }
  
  for (annot in annotations){
    for (gene_ami in setdiff(genes,gene)){
      #trouve annot min a voir pour etre considerer comme lié, depend du nb d'annot tot du gene le moins annot
      #si<10=> 1, si <20 2, sinon 3
      nAnnotsTotGeneMoinsAnnote<-min(c(length(getAnnot(gene,annot,df = genes_infos)),length(getAnnot(gene_ami,annot,df = genes_infos))))
      #print(nAnnotsTotGeneMoinsAnnote)
      if(nAnnotsTotGeneMoinsAnnote<7){
        annotMin<-1
      }else if(nAnnotsTotGeneMoinsAnnote<15){
        annotMin<-2
      }else if (nAnnotsTotGeneMoinsAnnote<25){
        annotMin<-3
      }else{
        annotMin<-4
      }
      
      # si gene a beaucoup trop d'annot, l'handicapé :
      nAnnotsTotGenePlusAnnote<-max(c(length(getAnnot(gene,annot,df = genes_infos)),length(getAnnot(gene_ami,annot,df=genes_infos))))
      if(nAnnotsTotGenePlusAnnote>25){
        annotMin<-annotMin+1
      }

      else if(nAnnotsTotGeneMoinsAnnote>50){
        annotMin<-annotMin+2
      }

      annot_commun<-FindAnnotCommun(gene1 = gene,gene2 = gene_ami,annotationType = annot,df=genes_infos)

      
      
      if (length(annot_commun)>=annotMin){
        print(paste("linked phenotype for ",gene,"and",gene_ami,"with :",paste(annot_commun,collapse = ", ")))
        amis<-c(amis,gene_ami)
      }
    }
  }
    
    
  
  return(amis)
}

getDfFrequenceAnnots<-function(genes, genes_infos, annotationType="all"){
  
  # fonction qui renvoit un df des annotations trie par nb qu'il apparait parmi les genes
  #marqueurs du cluster avec un importance d'au moins 3
  
  if (annotationType=="all"){
    annotations<-c("BioProcess","Pathway","MarqueurPops")
  }else{
    annotations<-annotationType
  }
  annot_df<-data.frame()
  
  for (gene in genes){
    
    annots<-getAnnot(gene, annotationType,df = genes_infos)
    for(typeDannotation in annotations){
      for (annot in annots){
        if (!(annot %in% rownames(annot_df))){
          annot_df[annot,"nb"]<-1
          annot_df[annot,"genes"]<-gene
          
        } else {
          annot_df[annot,"nb"]<-annot_df[annot,"nb"]+1
          annot_df[annot,"genes"]<-paste(annot_df[annot,"genes"],gene,sep=";")
        }
      }
    }
      
    }
    
  return(annot_df[order(annot_df[,"nb"],decreasing = T),])
}  
