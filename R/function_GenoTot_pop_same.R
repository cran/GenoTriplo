#' Genotype clusters for individuals from a same population
#'
#' @param resClustering object from clustering phase
#' @param Dataset dataset with Contrast and SigStren for each marker and individuals
#' @param SampleName list of samples name
#' @param NbClustMax number of cluster maximum (ploidy+1)
#' @param SeuilNoCall threshold of the probability of belonging to a cluster
#' @param SeuilSD threshold for the standard deviation of a cluster
#' @param SeuilNbSD thresold for the distance between an individuals and his cluster (x=Contrast)
#' @param cr_marker call rate threshold
#' @param fld_marker FLD threshold
#' @param hetso_marker HetSO threshold
#' @param tryEqHW bool. If TRUE, will delete cluster that shouldn't happen under H-W equilibrium (heterozygote with less individuals than both homozygote in triploid for example). Additionnaly, the cluster should contain less than 5% of the individuals. Else it stays.
#'
#' @import Rmixmod
#'
#' @return list with result of genotyping, vector of markers with maximum clusters and dataframe with metrics for each markers
#'
#' @keywords internal
#' @noRd

GenoAssign_pop_same = function(resClustering,Dataset,SampleName,NbClustMax,SeuilNoCall,SeuilSD,SeuilNbSD,cr_marker=NULL,fld_marker=NULL,hetso_marker=NULL,verbose=FALSE,
                               tryEqHW=TRUE){
  # Voir function_Geno_Tot_pop_dif.R pour les parties equivalentes commentees entre les deux fonctions
  if (is.null(Dataset$MarkerName) | is.null(Dataset$SampleName) | is.null(Dataset$Contrast) | is.null(Dataset$SigStren)){
    stop("One of SampleName, MarkerName, Contrast or SigStren is missing from Dataset.")
  }
  if (is.null(tryEqHW)){
    tryEqHW=FALSE
  } else if (is.na(tryEqHW)){
    tryEqHW=FALSE
  } else if (!is.logical(tryEqHW)){
    cat("tryEqHW set to FALSE.\n")
    tryEqHW=FALSE
  }
  MarkerName = names(resClustering) # Stockage des noms des marker
  list_max_clust=c()
  resGenoAssign = as.data.frame(matrix(nrow = length(MarkerName),ncol=length(SampleName))) # creation df de resultats (avant autant de ligne que de marqueur et de colonne que d'individu)
  colnames(resGenoAssign)=SampleName  # les colonnes (variables) prennent le noms des individus
  rownames(resGenoAssign)=MarkerName  # les lignes prennet le nom des marqueurs
  df_classif=data.frame(toKeep=rep(NA,length(MarkerName)),CR=NA,FLD=NA,HetSO=NA,HomRO=NA,nClus=NA,MAF=NA,p.HW=NA,Message=NA)
  rownames(df_classif)=MarkerName
  for (mn in MarkerName){
    if (verbose){cat(paste0("Marker : ",mn,"\n"))}
    if (length(resClustering[[mn]])==1){ # signifie quil ne contient que 'Error'
      df_classif[mn,]=c(FALSE,NA,NA,NA,NA,NA,NA,NA,"No clustering convergence.")
      resGenoAssign[mn,]=-1
    } else {
      tmp_partition = resClustering[[mn]]$df$partition
      tmp_means = resClustering[[mn]]$means
      # On verifie quil ny a pas dindividu trop eloignes sinon on les supprime
      nb_sd_indiv = NbSD_gp(genotype = tmp_partition,
                            data = Dataset[Dataset$MarkerName==mn,],
                            SampleName=SampleName)
      if (length(which(nb_sd_indiv>SeuilNbSD))>0){
        tmp_partition[nb_sd_indiv>SeuilNbSD] = NA
      }
      # On verifie malgre la suppression des individus trop eloignes quun gp n'est pas trop etale
      verif = verif_sd(genotype = tmp_partition,
                       data = Dataset[Dataset$MarkerName==mn,],
                       SeuilSD=SeuilSD)
      tmp_partition=verif[[1]]
      
      tmp_partition[resClustering[[mn]]$df$proba<SeuilNoCall] = NA # normalement n'en retire pas trop sauf si gros pb
      
      # Remise a niveau des partitions pour eviter qu'il y ait un trou (si un cluster a ete supprime completement)
      tmp_means = tmp_means[sort(unique(tmp_partition))] # sort() retire les NA
      tmp_partition = match(tmp_partition,sort(unique(tmp_partition)))
      
      # On range les clusters dans l'ordre croissant
      tmp_partition = match(tmp_means,sort(tmp_means))[tmp_partition]
      tmp_means = sort(tmp_means)
      
      # Check qu'il n'y a pas de cluster parasite (difficile a controler mais dans l'idee, les proportions doivent suivre des proportions en panmixie)
      # i.e. Equilibre de Hardy-Weinberg : il ne devrait pas y avoir un cluster intermediaire avec moins d'individu que les clusters qui l'entoure
      # Un rare cas serait une combinaison heterozygote deletere
      tmp_tab = table(tmp_partition)
      tmp_tryEqHW = tryEqHW
      while(tmp_tryEqHW){
        tmp_tryEqHW=FALSE
        lessHW=NULL
        if (length(tmp_tab)>2){
          for (k in seq_len(length(tmp_tab)-2)){
            if (tmp_tab[k+1]<tmp_tab[k] & tmp_tab[k+1]<tmp_tab[k+2] & tmp_tab[k+1]<sum(tmp_tab)*0.05){
              lessHW=k+1
            }
          }
        }
        if (!is.null(lessHW)){
          # On met NA si probleme puis on rearrange
          tmp_partition[tmp_partition==lessHW] = NA
          tmp_means = tmp_means[sort(unique(tmp_partition))] # sort() retire les NA
          tmp_partition = match(tmp_partition,sort(unique(tmp_partition)))
          tmp_tab = table(tmp_partition)
          tmp_tryEqHW = TRUE
        }
      }
      
      if (! is.character(resClustering[[mn]])){
        n_clust_restant = length(unique(tmp_partition[!is.na(tmp_partition)]))
        if (n_clust_restant>0){
          # On verifie le nombre de cluster restant apres avoir enlever les potentiels individus trop eloigne et les groupes trop etendus
          # val_clust=as.numeric(as.character(unique(tmp_partition[which(!is.na(tmp_partition))])))
          
          i=which.max(abs(tmp_means))
          inverse=FALSE
          if (n_clust_restant!=NbClustMax & n_clust_restant!=1){ # On verifie qu'il n'y a pas de probleme de proportion ()
            if (i==n_clust_restant){ # l'homozygote le plus extreme est celui avec un contrast positif
              i_min = which.min(tmp_means)
              if (tmp_tab[i_min]>2*tmp_tab[i]){ # ca veut dire que l'homozygote est probablement decale par rapport a 0
                inverse=TRUE
              }
            } else { # l'homozygote le plus extreme est celui avec un contrast negatif
              i_min = which.max(tmp_means)
              if (tmp_tab[i_min]>2*tmp_tab[i]){
                inverse=TRUE
              }
            }
          }
          if (inverse==TRUE){
            i = ifelse(i==1,n_clust_restant,1) # on change le sens de remplissage
          }
          if (i==1){ # on remplit de gauche a droite
            resGenoAssign[mn,] = NbClustMax - tmp_partition
          } else { # on remplit de gauche a droite
            resGenoAssign[mn,] = n_clust_restant - tmp_partition
          }
          
          resGenoAssign[mn,is.na(tmp_partition)]=-1 # tous les NA mis pendant l'analyse en -1 (no call)
        }
      }
      # On test si ce marker sera garde ou non
      ToKeep = keepMarkerPloidy(which = NbClustMax,
                                marker = mn,
                                genotypePop = resGenoAssign[mn,],
                                data = Dataset[Dataset$MarkerName==mn,],
                                cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
      df_classif[mn,]=ToKeep
      if (!is.na(df_classif[mn,"Message"])){
        if (df_classif[mn,"Message"]=='CR' & verif[[2]]=="SD_gp"){
          df_classif[mn,"Message"]="CR-SDgp"
        }
      }
      if (n_clust_restant==NbClustMax){
        list_max_clust=c(list_max_clust,mn)
      }
    }
  }
  df_classif = cbind(data.frame(MarkerName=rownames(df_classif)),df_classif)
  resGenoAssign[is.na(resGenoAssign)]=-1
  return(list(resGenoAssign,list_max_clust,df_classif))
}
