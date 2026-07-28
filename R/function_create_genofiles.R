#' Create files of genotype and indiv/marker stats from genotype .Rdata files
#'
#' @param dir.name.save directory with .Rdata files with genotype
#' @param new.dir new directory for saving files
#' @param ploidy ploidy of the population
#' @param corres_ATCG dataframe with correspondance between Allele_A/Allele_B to ATCG
#'
#' @importFrom rlang .data
#' @returns names of the indivCR and markerCR files
#' @export
#' 
#' @examples
#' \dontrun{
#' ploidy = 3
#' dir.name.save = system.file("extdata/output_genotyping/1", package = "GenoTriplo")
#' new.dir = "1"
#' create_genofiles(dir.name.save=dir.name.save,new.dir=new.dir,ploidy=ploidy)
#' }
#' 
create_genofiles = function(dir.name.save,new.dir,ploidy,corres_ATCG=NULL){
  cat("Creating files...\n")
  t1=Sys.time()
  if (!is.null(corres_ATCG)){
    if (is.null(corres_ATCG$probeset_id) | is.null(corres_ATCG$Allele_A) | is.null(corres_ATCG$Allele_B)){stop("Le fichier de correspondance doit au moins avoir probeset_id, Allele_A et Allele_B comme colonne.")}
  }
  # Fonction pour changement des 012 en AB et ATCG
  fx=function(x){
    x['geno']=gsub(pattern = "A",replacement = x['Allele_A'],x = x['geno']) # on remplace les A d'abord car A est aussi une base nucleic (eviter de remplacer quelque chose qui a deja ete remplace)
    x['geno']=gsub(pattern = "B",replacement = x['Allele_B'],x = x['geno'])
    return(x['geno'])
  }
  fx2=function(x){
    gsub(pattern = "/",replacement = " ",x = gsub(pattern = 'NA',replacement = '0',x = x))
  }
  # Gestion des fichiers de sortie
  # .map file
  write.table(file = paste0(dir.name.save,"/",new.dir,".map"),x = "#MarkerID",quote = F,row.names = F,col.names = F)
  # for indivCR file
  n_no_na=0
  n_tot=0
  # indivCR and markerCR
  file.indivCR = paste0(dir.name.save,"/",new.dir,"_indivCR.csv")
  file.markerCR = paste0(dir.name.save,"/",new.dir,"_markerCR.csv")
  # All genotype files
  listF=list.files(path=paste0(dir.name.save,"/all_geno"))
  res_genotyping=NULL # initialisation pour eviter 'no visible binding'
  for (fi in listF){
    load(paste0(dir.name.save,"/all_geno/",fi)) # res_genotyping
    # .map file
    write.table(x = data.frame(MarkerID=res_genotyping[[3]]$MarkerName),file = paste0(dir.name.save,"/",new.dir,".map"),quote = F,row.names = F,col.names = F,append = T)
    
    if (is.null(corres_ATCG)){
      write.table(x = cbind(add_categories(X=res_genotyping[[3]],ploidy=ploidy),data.frame(locat.file=rep(x = paste0(dir.name.save,"/all_geno/",fi),nrow(res_genotyping[[3]])))),file = file.markerCR,
                  sep = ";",quote = F,row.names = F,col.names = (fi==listF[1]),append=(fi!=listF[1]))
    } else {
      
      df.ajout = cbind(add_categories(X=res_genotyping[[3]],ploidy=ploidy),data.frame(locat.file=rep(x = paste0(dir.name.save,"/all_geno/",fi),nrow(res_genotyping[[3]])))) %>%
        left_join(corres_ATCG %>% select(.data$probeset_id,.data$Allele_A,.data$Allele_B),by=c("MarkerName"="probeset_id")) # ajoute Allele_A et Allele_B
      write.table(x = df.ajout,file = file.markerCR,
                  sep = ";",quote = F,row.names = F,col.names = (fi==listF[1]),append=(fi!=listF[1]))
    }
    n_no_na = n_no_na + apply(X = res_genotyping[[1]],MARGIN = 2,FUN = function(x) sum(!is.na(x)))
    n_tot = n_tot + nrow(res_genotyping[[1]])
  }
  # Compte le nombre d'individu total pour ecrire ensuite ligne par ligne les fichiers ped et genoAPIS
  n_indiv = ncol(res_genotyping[[1]])
  # indivCR file
  write.table(x = data.frame(SampleName=names(n_no_na),CR=round(n_no_na/n_tot,5)),file = file.indivCR,sep = ";",quote = F,row.names = F,col.names = T)
  # .ped file => En ligne les individus et leur genotype (les markers sont dans le meme ordre que le .map)
  write.table(file = paste0(dir.name.save,"/",new.dir,".ped"),x = c("#SampleName Genotypes"),quote = F,row.names = F,col.names = F)
  # .geno file => comme .ped mais en 0 a ploidy
  write.table(file = paste0(dir.name.save,"/",new.dir,".geno"),x = c("#SampleName Genotypes"),quote = F,row.names = F,col.names = F)
  
  # _genoATCG ou _genoAB file => Comme .ped mais ATCG ou AB avec / en separteur. Il y a aussi colnames (marker) et rownames(indiv)
  # _genoAPIS file => markerCR.csv as res_marker + genoATCG ou genoAB as data_APIS
  if (!is.null(corres_ATCG)){extension.name.geno = "_genoATCG.csv"}else{extension.name.geno = "_genoAB.csv"}
  list.marker=c()
  for (k in seq_len(n_indiv)){
    geno_indiv_k=c()
    for (fi in listF){
      load(paste0(dir.name.save,"/all_geno/",fi)) # res_genotyping
      geno_indiv_k = c(geno_indiv_k,res_genotyping[[1]][,k])
      if (k==1){
        list.marker=c(list.marker,rownames(res_genotyping[[1]]))
      }
    }
    samplenames = colnames(res_genotyping[[1]])
    # .geno file
    write.table(x=matrix(data = c(samplenames[k],geno_indiv_k),nrow = 1),file=paste0(dir.name.save,"/",new.dir,".geno"),append=TRUE,col.names = FALSE,quote = FALSE,row.names = F,sep = " ")
    
    # Corres ATCG et changement de -1 0 1 2 en AB
    # Changement en AB
    for (p in 1:(ploidy-1)){
      geno_indiv_k[geno_indiv_k==p]=paste(paste(rep("B",p),collapse="/"),paste(rep("A",ploidy-p),collapse="/"),sep = "/")
    }
    geno_indiv_k[geno_indiv_k==0]=paste(rep("A",ploidy),collapse="/")
    geno_indiv_k[geno_indiv_k==ploidy]=paste(rep("B",ploidy),collapse="/")
    
    if (!is.null(corres_ATCG)){ # si on a les correspndance entre A/B et ATCG
      corres_ATCG = corres_ATCG[match(list.marker,corres_ATCG$probeset_id),] %>%
        mutate(geno=geno_indiv_k)
      
      geno_indiv_k=as.vector(apply(X = corres_ATCG,MARGIN = 1,FUN = fx))
    }
    geno_indiv_k[geno_indiv_k==-1]=paste(rep("NA",ploidy),collapse="/")
    
    geno.apis=geno_indiv_k
    geno.ped=as.vector(sapply(X = geno_indiv_k,FUN = fx2))
    
    # .ped file
    write.table(x=matrix(data = c(samplenames[k],geno.ped),nrow = 1),file=paste0(dir.name.save,"/",new.dir,".ped"),append=TRUE,col.names = FALSE,quote = FALSE,row.names = F,sep = "  ")
    # genoATCG ou genoAB file
    if (k==1){
      write.table(x=matrix(data=c("SampleName",list.marker),nrow = 1),file=paste0(dir.name.save,"/",new.dir,extension.name.geno),col.names = FALSE,quote = FALSE,row.names = F,sep = ";")
    }
    write.table(x=matrix(data = c(samplenames[k],geno.apis),nrow = 1),file=paste0(dir.name.save,"/",new.dir,extension.name.geno),col.names = FALSE,quote = FALSE,row.names = F,sep = ";",append = T)
  }
  delay2 = Sys.time()-t1
  cat(paste0("End creating output files : ",round(delay2,3)," ",units(delay2),"\n"))
  # ------remove duplicat here-------
  # db_ind = c()
  # db_nam = c()
  # count_na = function(x){
  #   round(1-(length(which(x==-1))/length(x)),3)
  # }
  # for (pat in c("_2$","\\.2$","_BIS$","\\.BIS$")){ # add more si necessary
  #   indice = regexpr(pattern = pat,text = df_cr_ind$SampleName,fixed = FALSE) # fixed : so that . is not a regular expression replacing all character
  #   db_ind = c(db_ind,df_cr_ind$SampleName[which(indice!=-1)])
  #   db_nam = c(db_nam,df_cr_ind$SampleName[which(df_cr_ind$SampleName %in% substr(x = df_cr_ind$SampleName[which(indice!=-1)],
  #                                                                                 start = 1,
  #                                                                                 stop = nchar(df_cr_ind$SampleName[which(indice!=-1)])-nchar(pat,)))])
  # }
  # if (length(db_ind)>1){
  #   for (k in 1:length(db_ind)){
  #     cr1=df_cr_ind$CR[df_cr_ind$SampleName==db_ind[k]]
  #     cr2=df_cr_ind$CR[df_cr_ind$SampleName==db_nam[k]]
  #     if (cr1>cr2){
  #       data_APIS = as.matrix(data_APIS[-which(rownames(data_APIS)==db_nam[k]),]) # suppress
  #       rownames(data_APIS)[which(rownames(data_APIS)==db_ind[k])] = db_nam[k] # rename because there is a _2 or else
  #     } else {
  #       data_APIS = as.matrix(data_APIS[-which(rownames(data_APIS)==db_ind[k]),]) # just suppress, the good name is still here
  #     }
  #   }
  # }
  # save(data_APIS,res_marker,file = paste0("./output_genotyping/",save_n,"/",save_n,"_genoAPIS.Rdata"))
  # tab = res_geno[[3]][,c("toKeep","nClus")] # nest plus un table
  # df_cr_marker = data.frame(MarkerName = rownames(res_geno[[3]]),CR=res_geno[[3]][,"CR"])
  # list_return = list(df_cr_ind,df_cr_marker,tab)
  
  return(list(file.indivCR,file.markerCR))
}