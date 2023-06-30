#' Clustering function
#'
#' Clustering function to run clustering with no parallelization process nor auto save
#'
#' @param dataset dataset with Contrast and SigStren for each individuals and each markers
#' @param nb_clust_possible number of cluster possible (ploidy+1)
#' @param n_iter number of iterations to perform for clustering
#' @param Dmin minimal distance between two clusters
#'
#' @import dplyr
#' @import Rmixmod
#' @importFrom rlang .data
#'
#' @return list of results of clustering
#'
#' @export
#'
#' @examples
#' data(GenoTriplo_to_clust)
#' ploidy=3
#' res = Clustering(dataset=GenoTriplo_to_clust,
#'                  nb_clust_possible=ploidy+1,n_iter=5)
#'

Clustering = function(dataset,nb_clust_possible,n_iter=5,Dmin=0.28){

  listM=unique(dataset$MarkerName) # vecteur avec le nom des marqueurs
  n_tot=length(listM) # taille du vecteur

  res=list() # creation de la list de sortie de la fonction (le resultat du clustering)
  res[[n_tot]]=NA # initialisation d'une liste avec le bon nombre d'element (pas utile je crois en fait, mais ?a permet de gagner du temps en general)
  names(res)=listM # chaque element de la liste prend le nom dun marker
  for (k in (1:n_tot)){ # boucle sur lensemble du nombre des marqueurs (1<=k<=n_tot)
    m=listM[k] # le nom du marqueur actuel
    dta=dataset %>%
      filter(.data$MarkerName==m) %>% # on ne garde que les lignes du bon marqueur
      select(c("Contrast","SigStren")) # on ne garde que les deux colonnes qui nous interessent
    clust_opt = 2*nb_clust_possible # on copie le nombre de cluster possible dans un nouvel objet
    mixmodICL=MixmodBoucle(dta=dta,nb_clust_opt = clust_opt,iter=n_iter) # lancement de la fonction de clustering
    while (clust_opt>1 & mixmodICL@error){  # on refait tourner des modeles si il y a eu que des erreurs mais on diminue le nombre d'iterations (n-2) et le nombre de cluster possible
      clust_opt=clust_opt-1 # diminution du nombre de cluster possible max
      mixmodICL=MixmodBoucle(dta=dta,nb_clust_opt = clust_opt,iter=max(1,n_iter-2)) # on relance
    }
    if (mixmodICL@error){
      res[[k]]='Error'
    } else { # Si algorithme de clustering a reussi
      Opt=OptimizeCluster(mixmodICL = mixmodICL,nbClustMax = nb_clust_possible,distmin = Dmin) # on lance l'optimisation
      res[[k]] = Opt@bestResult # on stock le resultat
    }
  }
  return(res)
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

# fonction OptimizeCluster : permet de rassembler les clusters qui sont trop proches pour etre deux genotypes differents
# mixmodICL : le meilleur resultat du clustering
# nbClustMax : ploidy+1, le nombre de cluster max attendu
# distmin : la distance minimale entre 2 clusters (unite presque arbitraire)
#' Optimize the clustering algorithm
#'
#' @param mixmodICL the best model
#' @param nbClustMax number of clusters maximum (ploidy+1)
#' @param distmin minimal distance between two cluster
#'
#' @import Rmixmod
#'
#' @return the optimized model
#'
#' @keywords internal
#' @noRd

OptimizeCluster = function(mixmodICL,nbClustMax,distmin=0.28){
  TropProche=TRUE # on part du postulat qu'il y a des clusters trop proche (puisquon en demande plus que le nombre de genotype maximum)
  nbClust=mixmodICL@bestResult@nbCluster # on stock le nombre de cluster du modele
  # Verification que chaque cluster contient au moins un individu
  gp_non_vide = unique(mixmodICL@bestResult@partition)
  if (nbClust != length(gp_non_vide)){
    gp_vide = c(1:nbClust)[which(!1:nbClust %in% gp_non_vide)] # deja dans ordre croissant
    i = 0 # indice_a_rajouter au fur et a mesure de la boucle for (si on supprime gp 2, le 5e gp devient en fait le 4e)
    for (gp_suppr in gp_vide){
      nbClust = nbClust-1

      mixmodICL@bestResult@parameters@mean = matrix(mixmodICL@bestResult@parameters@mean[-(gp_suppr-i),],nrow=nbClust,ncol=2)
      mixmodICL@bestResult@parameters@variance=mixmodICL@bestResult@parameters@variance[-(gp_suppr-i)]
      mixmodICL@bestResult@parameters@proportions=mixmodICL@bestResult@parameters@proportions[-(gp_suppr-i)]
      mixmodICL@bestResult@proba=matrix(mixmodICL@bestResult@proba[,-(gp_suppr-i)],nrow=nrow(mixmodICL@bestResult@proba),ncol=nbClust)

      tmp=mixmodICL@bestResult@partition
      tmp[tmp>(gp_suppr-i)]=tmp[tmp>(gp_suppr-i)]-1
      mixmodICL@bestResult@partition=as.integer(tmp)

      mixmodICL@results[[1]]@nbCluster=mixmodICL@bestResult@nbCluster
      mixmodICL@results[[1]]@parameters=mixmodICL@bestResult@parameters
      mixmodICL@results[[1]]@proba=mixmodICL@bestResult@proba
      mixmodICL@results[[1]]@partition=as.integer(mixmodICL@bestResult@partition)

      i = i+1
    }

  }
  while(TropProche & nbClust>1){
    M=mixmodICL@bestResult@parameters@mean[,1]
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
    coef=ifelse(test = mdis<0,-mdis+1,mdis+1) # ce coefficient permet de tenir compte du fait quavec le contrast, le milieu est plus resserre et a linvese les bords sont plus etires
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
      # Change nb de cluster
      mixmodICL@bestResult@nbCluster=nbClust
      mixmodICL@nbCluster=nbClust
      # Change la valeur de la moyenne du groupe fusionne (ContrastCCS) : pondere par la taille des groupes qui fusionnes
      mixmodICL@bestResult@parameters@mean[clustproche.min,1] = mean(c(rep(mixmodICL@bestResult@parameters@mean[clustproche.min,1],length(mixmodICL@bestResult@partition[mixmodICL@bestResult@partition==clustproche.min])),
                                                                       rep(mixmodICL@bestResult@parameters@mean[clustproche.max,1],length(mixmodICL@bestResult@partition[mixmodICL@bestResult@partition==clustproche.max]))))
      # Change la valeur de la moyenne du gp fusionne (SigStren) : pondere par la taille des groupes qui fusionnes
      mixmodICL@bestResult@parameters@mean[clustproche.min,2] = mean(c(rep(mixmodICL@bestResult@parameters@mean[clustproche.min,2],length(mixmodICL@bestResult@partition[mixmodICL@bestResult@partition==clustproche.min])),
                                                                       rep(mixmodICL@bestResult@parameters@mean[clustproche.max,2],length(mixmodICL@bestResult@partition[mixmodICL@bestResult@partition==clustproche.max]))))
      # Retire la ligne de l'ancien groupe dans les mean
      mixmodICL@bestResult@parameters@mean = matrix(mixmodICL@bestResult@parameters@mean[-clustproche.max,],nrow=nbClust,ncol=2)
      # Change les valeurs des ecart-types (addition des sigmas pour mieux correspondre)
      mixmodICL@bestResult@parameters@variance[[clustproche.min]]=mixmodICL@bestResult@parameters@variance[[clustproche.min]]+
        mixmodICL@bestResult@parameters@variance[[clustproche.max]]
      # Enleve les valeurs decart-types du groupe supprime
      mixmodICL@bestResult@parameters@variance=mixmodICL@bestResult@parameters@variance[-clustproche.max]
      # Change la repartition des individus (on regroupe les deux plus proches & on diminue de 1 le numero de groupe des gp superieur)
      tmp=mixmodICL@bestResult@partition
      tmp[tmp==clustproche.max]=clustproche.min
      tmp[tmp>clustproche.max]=tmp[tmp>clustproche.max]-1
      mixmodICL@bestResult@partition=as.integer(tmp)
      # Met a jour les proportions dans les groupe
      mixmodICL@bestResult@parameters@proportions[clustproche.min] = mixmodICL@bestResult@parameters@proportions[clustproche.min]+
        mixmodICL@bestResult@parameters@proportions[clustproche.max]
      # Supprime le groupe qui nexiste plus
      mixmodICL@bestResult@parameters@proportions=mixmodICL@bestResult@parameters@proportions[-clustproche.max]
      # Met a jour des proba dappartenance a un groupe
      mixmodICL@bestResult@proba[,clustproche.min] = mixmodICL@bestResult@proba[,clustproche.min]+mixmodICL@bestResult@proba[,clustproche.max]
      # Supprime la colonne correspondant au groupe supprime
      mixmodICL@bestResult@proba=matrix(mixmodICL@bestResult@proba[,-clustproche.max],nrow=nrow(mixmodICL@bestResult@proba),ncol=nbClust)

      # Fait pareil avec result
      mixmodICL@results[[1]]@nbCluster=mixmodICL@bestResult@nbCluster
      mixmodICL@results[[1]]@parameters=mixmodICL@bestResult@parameters
      mixmodICL@results[[1]]@proba=mixmodICL@bestResult@proba
      mixmodICL@results[[1]]@partition=as.integer(mixmodICL@bestResult@partition)


    } else {
      TropProche=FALSE
    }
  }
  return(mixmodICL)
}
