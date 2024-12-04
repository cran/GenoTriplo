#' Launch genotyping phase in parallel
#'
#' Function that launch the genotyping phase from the dataset with SampleName, Contrast and SigStren for each markers and the result of the 'Run_clustering' function.
#'
#' @param data_clustering dataframe result from create dataset phase
#' @param res_clust object from clustering phase
#' @param ploidy ploidy of offspring
#' @param SeuilNoCall threshold of the probability of belonging to a cluster
#' @param SeuilNbSD threshold for the distance between an individuals and his cluster (x=Contrast)
#' @param SeuilSD threshold for the standard deviation of a cluster (SeuilSD*(1+0.5*abs(mean_contrast_cluster)))
#' @param n_core number of cores used for parallelization
#' @param corres_ATCG dataframe with the correspondence between A/B of AXAS and A/T/C/G (three columns : probeset_id, Allele_A, Allele_B)
#' @param pop Yes or No : are individuals from a same population
#' @param cr_marker call rate threshold
#' @param fld_marker FLD threshold
#' @param hetso_marker HetSO threshold
#' @param save_n name of the saving file. If '' no auto save and return value is changed
#' @param batch batch number in case of parallelization else ignore
#' @param ALL TRUE/FALSE whether the dataset has been cut or not (from the shiny app)
#' @param path_log path for log file when run by the shiny app
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import foreach
#' @import dplyr
#'
#' @return if save_n != '' : 3 objects list : dataframe with call rate by individuals, dataframe with call rate and other metrics of markers and another dataframe -- Automatically save results. Else : return list with genotype
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(GenoTriplo_to_clust)
#' data(GenoTriplo_to_geno)
#' res = Run_Genotyping(data_clustering=GenoTriplo_to_clust,
#'                      res_clust=GenoTriplo_to_geno,
#'                      ploidy=3)
#' }
#'

Run_Genotyping = function(data_clustering,res_clust,ploidy,SeuilNoCall=0.85,SeuilNbSD=2.8,SeuilSD=0.28,n_core=1,corres_ATCG=NULL,pop="Yes",cr_marker=0.97,fld_marker=3.4,hetso_marker=-0.3,save_n='',batch="",ALL=TRUE,path_log=''){
  # n_core = parallel::detectCores() - 2
  if (!is.numeric(ploidy)){
    stop("ploidy must be numeric.")
  } else if (ploidy>3 | ploidy<2){
    stop("ploidy must be 2 or 3 !")
  }
  if (!is.numeric(n_core)){
    stop("n_core must be numeric.")
  }else if (parallel::detectCores()<n_core){
    n_core=parallel::detectCores()
    warning("The number of core asked is to high : will be set to maximum.")
  }
  if (save_n!=''){
    if (!dir.exists("./output_genotyping")){dir.create("./output_genotyping")}
    if (!dir.exists(paste0("./output_genotyping/",save_n))){
      dir.create(paste0("./output_genotyping/",save_n),recursive = TRUE)
    }
  }
  
  MarkerId = names(res_clust) # stockage des noms des differents marqueurs
  nb_of_marker = length(MarkerId) # stockage du nombre total de marker
  if (nb_of_marker<500 & pop!="Yes" & save_n!='' & path_log!=''){
    write(x = "Warning : The number of marker might not be enough.",file=path_log,append=T)
  } else if (nb_of_marker<500 & pop!="Yes" & save_n==''){
    warning("The number of marker might not be enough.")
  }
  batchsize = ifelse(nb_of_marker<n_core,nb_of_marker,nb_of_marker%/%n_core)
  i_max = nb_of_marker %/% batchsize # quit a ce que le dernier batch soit un peu plus grand (de quelques markers)
  if (batchsize>2000){
    batchsize = batchsize%/%2
    i_max = 2*i_max
    while (batchsize>2000){
      batchsize = batchsize%/%2
      i_max = 2*i_max
    }
  } else if (batchsize<500 & i_max>1 & pop != "Yes"){
    batchsize = batchsize*2
    i_max = i_max%/%2
    while(batchsize<500 & i_max>1){
      batchsize = batchsize*2
      i_max = ifelse(i_max%%2==0,i_max/2,i_max%/%2+1)
    }
  }

  if (save_n!='' & path_log!=''){
    write(x = paste0("Batchsize : ",batchsize," ; Nb Batch : ",i_max),file = path_log,append=T)
  }
  SampleName = unique(data_clustering$SampleName)
  t0 = Sys.time()
  # Parametrage des clusters
  clust_name = parallel::makeCluster(min(n_core,i_max))
  doParallel::registerDoParallel(clust_name)
  iter=NULL
  # Boucle foreach qui permet de paralleliser les taches
  res_paral_geno = foreach(iter=1:i_max) %dopar% {
    if (pop=="Yes"){ # Lance genotypage avec la fonction adequate
      if (iter != i_max){ # petite particularite si iter == imax => le dernier batch ne fait pas forcement la taille maximum de batchsize
        GenoAssign_pop_same(resClustering = res_clust[(1+(iter-1)*batchsize):(iter*batchsize)],
                            SampleName = SampleName,
                            SeuilNoCall = SeuilNoCall,SeuilNbSD = SeuilNbSD,SeuilSD = SeuilSD,
                            NbClustMax = ploidy+1,
                            Dataset=data_clustering[data_clustering$MarkerName %in% names(res_clust)[(1+(iter-1)*batchsize):(iter*batchsize)],],
                            cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
      } else {
        # Lance le genotypage avec le dernier batch qui va en realiste jusquau dernier marker (au cas ou le batch ne soit pas complet)
        GenoAssign_pop_same(resClustering = res_clust[(1+(iter-1)*batchsize):nb_of_marker],
                            SampleName = SampleName,
                            SeuilNoCall = SeuilNoCall,SeuilNbSD = SeuilNbSD,SeuilSD = SeuilSD,
                            NbClustMax = ploidy+1,
                            Dataset=data_clustering[data_clustering$MarkerName %in% names(res_clust)[(1+(iter-1)*batchsize):nb_of_marker],],
                            cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
      }
    } else if (pop!="Yes"){ # Lance genotypage avec la fonction adequate
      if (iter != i_max){
        GenoAssign_pop_dif(resClustering = res_clust[(1+(iter-1)*batchsize):(iter*batchsize)],
                           SampleName = SampleName,
                           SeuilNoCall = SeuilNoCall,SeuilNbSD = SeuilNbSD,SeuilSD = SeuilSD,
                           NbClustMax = ploidy+1,
                           Dataset=data_clustering[data_clustering$MarkerName %in% names(res_clust)[(1+(iter-1)*batchsize):(iter*batchsize)],],
                           cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
      } else {
        # Lance le genotypage avec le dernier batch qui va en realiste jusquau dernier marker (au cas ou le batch ne soit pas complet)
        GenoAssign_pop_dif(resClustering = res_clust[(1+(iter-1)*batchsize):nb_of_marker],
                           SampleName = SampleName,
                           SeuilNoCall = SeuilNoCall,SeuilNbSD = SeuilNbSD,SeuilSD = SeuilSD,
                           NbClustMax = ploidy+1,
                           Dataset=data_clustering[data_clustering$MarkerName %in% names(res_clust)[(1+(iter-1)*batchsize):nb_of_marker],],
                           cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
      }
    }
  }
  # Arrete des clusters (indispensable)
  parallel::stopCluster(clust_name)
  doParallel::stopImplicitCluster()
  t1=Sys.time()
  delay = t1-t0

  res_geno=res_paral_geno[[1]] # stockage de la premiere iteration de parallelisation
  if (i_max>1){ # si il y a eu plus dun lancement (parallelisation)
    for (k in 2:i_max){ # stockage du reste des resultats de la maniere dont il convient pour chaque list/df/vecteur
      res_geno[[1]]=rbind(res_geno[[1]],res_paral_geno[[k]][[1]])
      res_geno[[2]]=c(res_geno[[2]],res_paral_geno[[k]][[2]])
      res_geno[[3]]=rbind(res_geno[[3]],res_paral_geno[[k]][[3]])
    }
  }
  res_geno[[4]]=res_geno[[1]] # creation dun 4e element a la liste de resultat
  # Changement des genotype 0,1,2,... en A/A A/B ou A/A/A, A/A/B,... suivant la ploidy
  for (k in 1:(ploidy-1)){
    res_geno[[4]][res_geno[[4]]==k]=paste(paste(rep("B",k),collapse="/"),paste(rep("A",ploidy-k),collapse="/"),sep = "/")
  }
  res_geno[[4]][res_geno[[4]]==0]=paste(rep("A",ploidy),collapse="/")
  res_geno[[4]][res_geno[[4]]==ploidy]=paste(rep("B",ploidy),collapse="/")

  fx=function(vect){
    L=length(vect)
    A=vect[L-1]
    B=vect[L]
    vect=gsub(pattern = "A",replacement = A,x = vect) # on remplace les A d'abord car A est aussi une base nucleic (eviter de remplacer quelque chose qui a deja ete remplace)
    vect=gsub(pattern = "B",replacement = B,x = vect)
    return(vect)
  }

  res_geno[[4]]$probeset_id=rownames(res_geno[[4]])
  if (!is.null(corres_ATCG)){ # si on a les correspndance entre A/B et ATCG
    dta=left_join(x = res_geno[[4]],y = corres_ATCG,by = "probeset_id") %>% # left_join sur le probeset_id pour avoir les genotypes dun marker et la correspondence allele A et B a cote
      select(all_of(c(colnames(res_geno[[4]]),"Allele_A","Allele_B"))) %>% # on ne garde que les colonnes qui nousinteressent
      select(-c("probeset_id")) %>% # on retire probeset_id
      apply(FUN = fx,MARGIN = 1) # on change les A/A par T/T (exemple) via la fonction fx (un peu plus haut)
    dta = dta[-which(rownames(dta) %in% c('Allele_A','Allele_B')),] # on retire les alleles de reference
    colnames(dta)=res_geno[[4]]$probeset_id # nomme les colonnes avec le nom des marqueurs
    dta[dta==-1]=paste(rep("NA",ploidy),collapse="/")
    res_geno[[4]] = dta
    data_APIS = dta
  } else {
    res_geno[[4]][res_geno[[4]]==-1]=paste(rep("NA",ploidy),collapse="/")
    res_geno[[4]] = res_geno[[4]] %>%
      select(-c("probeset_id")) %>%
      t()
    data_APIS = res_geno[[4]]
  }
  # Sauvegarde de tous les resultats sous differents nom suivant le resultat (APIS/geno reel/CR/all)
  res_marker = add_categories(X=res_geno[[3]],ploidy=ploidy) # %>% filter(.,toKeep)
  res_geno[[3]]=res_marker

  count_na = function(x){
    round(1-(length(which(x==-1))/length(x)),3)
  }
  df_cr_ind = data.frame(SampleName=names(res_geno[[1]]),CR=apply(X = res_geno[[1]],MARGIN = 2,FUN = count_na))
  if (ALL){ # if dataset has not been cut -> remove duplicat here
    db_ind = c()
    db_nam = c()
    for (pat in c("_2$","\\.2$","_BIS$","\\.BIS$")){ # add more si necessary
      indice = regexpr(pattern = pat,text = df_cr_ind$SampleName,fixed = FALSE) # fixed : so that . is not a regular expression replacing all character
      db_ind = c(db_ind,df_cr_ind$SampleName[which(indice!=-1)])
      db_nam = c(db_nam,df_cr_ind$SampleName[which(df_cr_ind$SampleName %in% substr(x = df_cr_ind$SampleName[which(indice!=-1)],
                                                                                    start = 1,
                                                                                    stop = nchar(df_cr_ind$SampleName[which(indice!=-1)])-nchar(pat,)))])
    }
    if (length(db_ind)>1){
      for (k in 1:length(db_ind)){
        cr1=df_cr_ind$CR[df_cr_ind$SampleName==db_ind[k]]
        cr2=df_cr_ind$CR[df_cr_ind$SampleName==db_nam[k]]
        if (cr1>cr2){
          data_APIS = as.matrix(data_APIS[-which(rownames(data_APIS)==db_nam[k]),]) # suppress
          rownames(data_APIS)[which(rownames(data_APIS)==db_ind[k])] = db_nam[k] # rename because there is a _2 or else
        } else {
          data_APIS = as.matrix(data_APIS[-which(rownames(data_APIS)==db_ind[k]),]) # just suppress, the good name is still here
        }
      }
    }
    res_geno[[4]]=data_APIS
  }
  if (ALL & save_n!=''){
    write.table(x = res_marker,file = paste0("./output_genotyping/",save_n,"/",save_n,"_markerCR.csv"),sep = ";",row.names = FALSE)
    save(data_APIS,res_marker,file = paste0("./output_genotyping/",save_n,"/",save_n,"_genoAPIS.Rdata"))
    if (!is.null(corres_ATCG)){
      write.table(x = res_geno[[4]],file = paste0("./output_genotyping/",save_n,"/",save_n,"_genoATCG.csv"),sep = ";",row.names = TRUE,col.names = NA)
    } else {
      write.table(x = res_geno[[4]],file = paste0("./output_genotyping/",save_n,"/",save_n,"_genoAB.csv"),sep = ";",row.names = TRUE,col.names = NA)
    }
    save(res_geno,data_clustering,delay,file=paste0("./output_genotyping/",save_n,"/",save_n,"_genotyped.Rdata"))
  } else if (!ALL & save_n!='') {
    write.table(x = res_marker,file = paste0("./output_genotyping/",save_n,"/",save_n,"_",batch,"_markerCR.csv"),sep = ";",row.names = FALSE)
    save(data_APIS,res_marker,file = paste0("./output_genotyping/",save_n,"/",save_n,"_",batch,"_genoAPIS.Rdata"))
    if (!is.null(corres_ATCG)){
      write.table(x = res_geno[[4]],file = paste0("./output_genotyping/",save_n,"/",save_n,"_",batch,"_genoATCG.csv"),sep = ";",row.names = TRUE,col.names = NA)
    } else {
      write.table(x = res_geno[[4]],file = paste0("./output_genotyping/",save_n,"/",save_n,"_",batch,"_genoAB.csv"),sep = ";",row.names = TRUE,col.names = NA)
    }
    save(res_geno,data_clustering,delay,file=paste0("./output_genotyping/",save_n,"/",save_n,"_",batch,"_genotyped.Rdata"))
  } else if (save_n==''){
    return(res_geno)
  }
  if (save_n!=''){
    tab = res_geno[[3]][,c("toKeep","nClus")] # nest plus un table
    df_cr_marker = data.frame(MarkerName = rownames(res_geno[[3]]),CR=res_geno[[3]][,"CR"])
    list_return = list(df_cr_ind,df_cr_marker,tab)
    return(list_return)
  }
}

