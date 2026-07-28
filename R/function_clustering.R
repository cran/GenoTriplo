#' Clustering function from file
#'
#' Clustering function to run clustering with no parallelization process and autosave
#' Used by Clustering_from_dir. Shouldn't be used else.
#'
#' @param dirname directory of the file to clusterize.
#' @param filename filename of the data to clusterize
#' @param ploidy ploidy of individuals
#' @param n_iter number of iterations to perform for clustering
#' @param Dmin minimal distance between two clusters
#' @param new.dir new directory name for saving
#'
#' @import dplyr
#' @import Rmixmod
#' @importFrom stringr str_extract_all
#' @importFrom rlang .data
#'
#' @export
#' 
#' @return Autosave of the results
#' 
#' @examples
#' \dontrun{
#' dir.create = system.file("extdata/output_create/1", package = "GenoTriplo")
#' file.create = system.file("extdata/output_create/1","1_1_to_clust.Rdata", package = "GenoTriplo")
#' ploidy=3
#' Clustering_from_file(dirname=dir.create,filename=file.create,ploidy=ploidy)
#' }
#' 

Clustering_from_file = function(dirname,filename,ploidy,n_iter=5,Dmin=0.28,new.dir=NULL){
  if (is.null(new.dir)){
    new.dir=length(list.dirs("./output_clustering",recursive = F))+1
  }
  if (length(grep(pattern = "output_clustering",x = new.dir))==0){
    dir.name.save = paste0("./output_clustering/",new.dir)
  } else if (length(grep(pattern = "^./",x = new.dir,fixed = F))==0){
    dir.name.save=paste0("./",new.dir)
  } else {dir.name.save=new.dir}
  
  if (!dir.exists(dir.name.save)){
    dir.create(path = dir.name.save,recursive = TRUE)
    cat(paste0("Directory created : ",dir.name.save,"\n"))
  }
  path.data_clust = paste0(dirname,"/",filename)
  tmp = load(path.data_clust) # data_clustering
  # Permet d'eviter la note 'no visible bunding'
  data_clustering = eval(parse(text = tmp))
  new.filename=sub(pattern = "to_clust",replacement = "to_geno",x = filename)
  res_clust = Clustering(data_clustering=data_clustering,nb_clust_possible=ploidy+1,n_iter=n_iter,Dmin=Dmin)
  if (file.exists(paste0(dir.name.save,"/",new.filename))){cat(paste0("Attention, le fichier ",paste0(dir.name.save,"/",new.filename),"vient d'etre remplace!\n"))}
  save(res_clust,data_clustering,file = paste0(dir.name.save,"/",new.filename))
}

#' Clustering function
#'
#' Clustering function to run clustering with no parallelization process nor auto save
#'
#' @param data_clustering dataset with Contrast and SigStren for each individuals (as SampleName) and each markers (as MarkerName)
#' @param nb_clust_possible number of cluster possible (ploidy+1)
#' @param n_iter number of iterations to perform for clustering
#' @param Dmin minimal distance between two clusters
#'
#' @import dplyr
#' @import Rmixmod
#' @importFrom rlang .data
#'
#' @return list of results of clustering
#' @export
#'
#' @examples
#' \dontrun{
#' file.create = system.file("extdata/output_create/1","1_1_to_clust.Rdata", package = "GenoTriplo")
#' load(file.create)
#' ploidy=3
#' Clustering(data_clustering=data_clustering,nb_clust_possible=ploidy+1)
#' }
#' 
#'
Clustering=function(data_clustering,nb_clust_possible,n_iter=5,Dmin=0.28){
  if (is.null(data_clustering$SampleName) | is.null(data_clustering$MarkerName) | is.null(data_clustering$Contrast) | is.null(data_clustering$SigStren)){
    stop("One of SampleName, MarkerName, Contrast or SigStren is missing from dataset.")
  }
  SampleName = unique(data_clustering$SampleName)
  listM=unique(data_clustering$MarkerName) # vecteur avec le nom des marqueurs
  n_tot=length(listM) # taille du vecteur
  
  res_clust=vector("list", n_tot) # creation de la list de sortie de la fonction (le resultat du clustering)
  names(res_clust)=listM # chaque element de la liste prend le nom dun marker
  for (k in (1:n_tot)){ # boucle sur lensemble du nombre des marqueurs (1<=k<=n_tot)
    m=listM[k] # le nom du marqueur actuel
    dta=data_clustering[data_clustering$MarkerName==m & !is.infinite(data_clustering$Contrast) & !is.na(data_clustering$Contrast),]
    if (nrow(dta)==0){
      df = data.frame(partition=rep(NA,length(SampleName)),proba=rep(0,length(SampleName)))
      means = NA
      res_clust[[k]] = list(df=df,means=means) # on stock le resultat
    } else {
      sn = unique(dta$SampleName)
      dta = dta[,c("Contrast","SigStren")] # on ne garde que les deux colonnes qui nous interessent
      clust_opt = min(2*nb_clust_possible,nrow(dta)) # on copie le nombre de cluster possible dans un nouvel objet et on sassure quil y a assez dinvididu
      mixmodICL=MixmodBoucle(dta=dta,nb_clust_opt = clust_opt,iter=n_iter) # lancement de la fonction de clustering
      while (clust_opt>1 & mixmodICL@error){  # on refait tourner des modeles si il y a eu que des erreurs mais on diminue le nombre d'iterations (n-2) et le nombre de cluster possible
        clust_opt=clust_opt-1 # diminution du nombre de cluster possible max
        mixmodICL=MixmodBoucle(dta=dta,nb_clust_opt = clust_opt,iter=max(1,n_iter-2)) # on relance
      }
      if (mixmodICL@error){
        res_clust[[k]]='Error'
      } else { # Si algorithme de clustering a reussi
        Opt=OptimizeCluster(mixmodICL = mixmodICL,nbClustMax = nb_clust_possible,distmin = Dmin) # on lance l'optimisation
        partition = Opt[[1]]
        proba = apply(Opt[[2]],MARGIN=1,FUN=max)
        df = data.frame(partition=partition,proba=proba,SampleName=sn)
        df = left_join(data.frame(SampleName=SampleName),df,by=c("SampleName"="SampleName")) %>% select(partition,proba)
        df$proba[is.na(df$partition)]=0
        means = Opt[[3]]
        res_clust[[k]] = list(df=df,means=means) # on stock le resultat
      }
    }
  }
  return(res_clust)
}

#' Loop of clustering
#'
#' @param dta dataset with Contrast and SigStren only
#' @param nb_clust_opt (ploidy +1) * 2 by default
#' @param iter number of iterations to perform the clustering algorithm
#'
#' @import Rmixmod
#'
#' @return the best clustering model among the iter tried
#'
#' @keywords internal
#' @noRd

MixmodBoucle = function(dta,nb_clust_opt,iter=5){
  mod = mixmodCluster(data=dta,nbCluster=nb_clust_opt,criterion="ICL") # lancement du clustering avec 2*nb_clust_possible
  # on lance avec bcp de cluster possible pour maximiser les chances de trouver les genotypes avec peu d'individus
  if (iter>1){ # si on a demande plus dune iteration
    for (i in 1:(iter-1)){  # 'iter' models sont fittes au total
      mod2 = mixmodCluster(data=dta,nbCluster=nb_clust_opt,criterion="ICL")
      if (!mod2@error & !mod@error){ # si on a pas eu derreur
        if (mod2@bestResult@likelihood>mod@bestResult@likelihood){ # et que likelihood meilleur avec le nouveau model
          mod=mod2 # on le stock (sinon on garde l'ancien et on recommence pour faire le bon nombre d'iteration)
        }
      } else if (!mod2@error & mod@error){ # si on avait une erreur avant, on garde le nouveau si il na pas derreur (on ne peut pas comparer les likelihood vu uqil ny en avait pas avant)
        mod=mod2
      }
    }
  }
  return(mod) # on retourne le model qui maximise la likelihood (vraissemblance)
}

#' Optimize the clustering algorithm
#'
#' @param mixmodICL the best model
#' @param nbClustMax number of clusters maximum (ploidy+1)
#' @param distmin minimal distance between two cluster
#'
#' @import Rmixmod
#' @importFrom stats weighted.mean
#' @return the optimized model
#'
#' @keywords internal
#' @noRd

OptimizeCluster = function(mixmodICL,nbClustMax,distmin=0.28){
  TropProche=TRUE # on part du postulat qu'il y a des clusters trop proche (puisquon en demande plus que le nombre de genotype maximum)
  nbClust=mixmodICL@bestResult@nbCluster # on stock le nombre de cluster du modele
  # Verification que chaque cluster contient au moins un individu
  means = mixmodICL@bestResult@parameters@mean[,1]
  proba = mixmodICL@bestResult@proba
  partition = mixmodICL@bestResult@partition
  gp_non_vide = unique(partition)
  if (nbClust != length(gp_non_vide)){
    gp_vide = c(1:nbClust)[which(!1:nbClust %in% gp_non_vide)] # deja dans ordre croissant
    i = 0 # indice_a_rajouter au fur et a mesure de la boucle for (si on supprime gp 2, le 5e gp devient en fait le 4e)
    for (gp_suppr in gp_vide){
      nbClust = nbClust-1
      
      means = means[-(gp_suppr-i)]
      proba=matrix(proba[,-(gp_suppr-i)],nrow=nrow(proba),ncol=nbClust)
      
      tmp=partition
      tmp[tmp>(gp_suppr-i)]=tmp[tmp>(gp_suppr-i)]-1
      partition=as.integer(tmp)
      
      i = i+1
    }
  }
  while(TropProche & nbClust>1){
    M=means
    ord = order(M)
    sor = sort(M)
    dis = c()
    mdis=c()
    coef=c()
    for (k in 1:(length(M)-1)){
      d = sor[k+1]-sor[k]
      dis=c(dis,d)
      mdis=c(mdis,(sor[k+1]+sor[k])/2)
    }
    coef=1+abs(mdis) # ifelse(test = mdis<0,-mdis+1,mdis+1) # ce coefficient permet de tenir compte du fait quavec le contrast, le milieu est plus resserre et a linvese les bords sont plus etires
    threshold = coef*distmin
    
    isProxy = dis-threshold
    if (min(isProxy)<0 | nbClust>nbClustMax){ # si on a deux groupes trop proche ou quon a tjr trop de groupe : on fusionne
      n=which(isProxy==min(isProxy))
      a=ord[n]
      b=ord[n+1]
      nbClust=nbClust-1
      if (a<b){
        clustproche.min=a
        clustproche.max=b
      } else {
        clustproche.min=b
        clustproche.max=a
      }
      # Change la valeur de la moyenne du groupe fusionne (ContrastCCS) : pondere par la taille des groupes qui fusionnes
      means[clustproche.min] = weighted.mean(x = means[c(clustproche.min,clustproche.max)],
                                             w = c(length(partition[partition==clustproche.min]),
                                                   length(partition[partition==clustproche.max])))
      # Retire la ligne de l'ancien groupe dans les mean
      means = means[-clustproche.max]
      # Change la repartition des individus (on regroupe les deux plus proches & on diminue de 1 le numero de groupe des gp superieur)
      tmp=partition
      tmp[tmp==clustproche.max]=clustproche.min
      tmp[tmp>clustproche.max]=tmp[tmp>clustproche.max]-1
      partition=as.integer(tmp)
      # Met a jour des proba dappartenance a un groupe
      proba[,clustproche.min] = proba[,clustproche.min]+proba[,clustproche.max]
      # Supprime la colonne correspondant au groupe supprime
      proba=matrix(proba[,-clustproche.max],nrow=nrow(proba),ncol=nbClust)
    } else {
      TropProche=FALSE
    }
  }
  return(list(partition,proba,means))
}
