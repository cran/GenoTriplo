#' Genotype clusters for individuals from a same population
#'
#' @param resClustering object from clustering phase
#' @param SampleName list of samples name
#' @param NbClustMax number of cluster maximum (ploidy+1)
#' @param SeuilNoCall threshold of the probability of belonging to a cluster
#' @param SeuilSD threshold for the standard deviation of a cluster
#' @param SeuilNbSD thresold for the distance between an individuals and his cluster (x=Contrast)
#' @param Dataset dataset with Contrast and SigStren for each marker and individuals
#' @param cr_marker call rate threshold
#' @param fld_marker FLD threshold
#' @param hetso_marker HetSO threshold
#'
#' @import Rmixmod
#'
#' @return list with result of genotyping, vector of markers with maximum clusters and dataframe with metrics for each markers
#'
#' @keywords internal
#' @noRd

GenoAssign_pop_same = function(resClustering,SampleName,NbClustMax,SeuilNoCall,SeuilSD,SeuilNbSD,Dataset,cr_marker=NULL,fld_marker=NULL,hetso_marker=NULL,verbose=FALSE){
  # Voir function_Geno_Tot_pop_dif.R pour les parties equivalentes commentees entre les deux fonctions
  MarkerName = names(resClustering) # Stockage des noms des marker
  list_max_clust=c()
  resGenoAssign = as.data.frame(matrix(nrow = length(MarkerName),ncol=length(SampleName))) # creation df de resultats (avant autant de ligne que de marqueur et de colonne que d'individu)
  names(resGenoAssign)=SampleName  # les colonnes (variables) prennent le noms des individus
  rownames(resGenoAssign)=MarkerName  # les lignes prennet le nom des marqueurs
  df_classif=data.frame(toKeep=rep(NA,length(MarkerName)),CR=NA,FLD=NA,HetSO=NA,HomRO=NA,nClus=NA,MAF=NA,Message=NA)
  rownames(df_classif)=MarkerName
  if (NbClustMax==3){
    if (is.null(Dataset$MarkerName)){
      stop("Your dataset must contain MarkerName as variable.")
    }
    for (k in MarkerName){ # on trouve deja ceux avec 3 groupes : assignation automatique + moyenne pour les suivants
      # On verifie quil ny a pas dindividu trop eloignes sinon on les supprime
      nb_sd_indiv = NbSD_gp(genotype = resClustering[[k]]@partition,
                            data = Dataset[Dataset$MarkerName==k,],
                            SampleName=SampleName)
      if (length(which(nb_sd_indiv>SeuilNbSD))>0){
        resClustering[[k]]@partition=as.factor(resClustering[[k]]@partition)
        resClustering[[k]]@partition[nb_sd_indiv>SeuilNbSD] = NA
      }
      # On verifie malgre la suppression des individus trop eloignes quun gp n'est pas trop etale
      verif = verif_sd(genotype = resClustering[[k]]@partition,
                       data = Dataset[Dataset$MarkerName==k,],
                       SeuilSD=SeuilSD)
      resClustering[[k]]@partition=as.factor(verif[[1]])
      if (verbose){print(paste0("Marker : ",k))}
      if (! is.character(resClustering[[k]])){
        # On verifie le nombre de cluster restant apres avoir enlever les potentiels individus trop eloigne et les groupes trop etendus
        val_clust=as.numeric(as.character(unique(resClustering[[k]]@partition[which(!is.na(resClustering[[k]]@partition))])))
        n_clust_restant = length(val_clust)
        if (n_clust_restant==3){ # si nombre max (diploide) : assignation des genotype dans lordre
          ordre = order(resClustering[[k]]@parameters@mean[,1])
          tmp = as.factor(resClustering[[k]]@partition)
          tmp2=rep(NA,length(tmp))
          for (l in 1:3){
            tmp2[tmp==ordre[l]]=l
          }
          resGenoAssign[k,]=-as.numeric(as.character(tmp2))+4-1 # y=ax+b avec y=1=>x=3;y=2=>x=2;y=1=>x=1 puis -1 car valeur = 0,1,2 et non 1,2,3
        } else if (n_clust_restant==2){ # si 2 : on regarde le plus extreme qui devient homoz et lautre heteroz
          m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
          m2=resClustering[[k]]@parameters@mean[val_clust[2],1]
          if (abs(m1)>abs(m2)){
            resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=1
            if (m1>m2){
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=0
            } else {
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=2
            }
          } else {
            resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=1
            if (m1>m2){
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=2
            } else {
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=0
            }
          }
        } else if (n_clust_restant==1){ # si 1 : homoz associe a la valeur moyenne du cluster (fct de si positif ou negatif)
          m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
          if (m1>0){
            resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=0
          } else {
            resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=2
          }
        }
        tmp3 = apply(resClustering[[k]]@proba,MARGIN=1,FUN=max)
        resGenoAssign[k,tmp3<SeuilNoCall]=-1
        resGenoAssign[k,is.na(resClustering[[k]]@partition)]=-1

        # On test si ce marker sera garde ou non
        ToKeep = keepMarkerdiplo(marker = k,
                                 genotypePop = resGenoAssign[k,],
                                 data = Dataset[Dataset$MarkerName==k,],
                                 cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
        df_classif[k,]=ToKeep
        if (!is.na(df_classif[k,"Message"])){
          if (df_classif[k,"Message"]=='CR' & verif[[2]]=="SD_gp"){
            df_classif[k,"Message"]="CR-SDgp"
          }
        }
        if (n_clust_restant==3){
          list_max_clust=c(list_max_clust,k)
        }
      }
    }
  } else if (NbClustMax==4){ # si triploide
    if (is.null(Dataset$MarkerName)){
      stop("Your dataset must contain MarkerName as variable.")
    }
    for (k in MarkerName){ # on trouve deja ceux avec 3 groupes : assignation automatique + moyenne pour les suivants
      # On verifie quil ny a pas dindividu trop eloignes sinon on les supprime
      nb_sd_indiv = NbSD_gp(genotype = resClustering[[k]]@partition,
                            data = Dataset[Dataset$MarkerName==k,],
                            SampleName=SampleName)
      if (length(which(nb_sd_indiv>SeuilNbSD))>0){
        resClustering[[k]]@partition=as.factor(resClustering[[k]]@partition)
        resClustering[[k]]@partition[nb_sd_indiv>SeuilNbSD] = NA
      }
      # On verifie malgre la suppression des individus trop eloignes quun gp n'est pas trop etale
      verif = verif_sd(genotype = resClustering[[k]]@partition,
                       data = Dataset[Dataset$MarkerName==k,],
                       SeuilSD=SeuilSD)
      resClustering[[k]]@partition=as.factor(verif[[1]])

      if (verbose){print(paste0("Marker : ",k))}
      if (! is.character(resClustering[[k]])){
        val_clust = as.numeric(as.character(unique(resClustering[[k]]@partition[which(!is.na(resClustering[[k]]@partition))])))
        n_clust_restant = length(val_clust)
        if (n_clust_restant==4){ # si 4 clusters : dans lordre
          ordre = order(resClustering[[k]]@parameters@mean[,1])
          tmp = as.factor(resClustering[[k]]@partition)
          tmp2=rep(NA,length(tmp))
          for (l in 1:4){
            tmp2[tmp==ordre[l]]=l
          }
          resGenoAssign[k,]=-as.numeric(as.character(tmp2))+5-1
        } else if (n_clust_restant==3){ # si 3 clusters : clusters avec la valuer la plus extreme assigne homoz associe puis dans lordre (les deux autre heteroz)
          m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
          m2=resClustering[[k]]@parameters@mean[val_clust[2],1]
          m3=resClustering[[k]]@parameters@mean[val_clust[3],1]
          prop1=sum(resClustering[[k]]@partition==val_clust[1],na.rm = TRUE)
          prop2=sum(resClustering[[k]]@partition==val_clust[2],na.rm = TRUE)
          prop3=sum(resClustering[[k]]@partition==val_clust[3],na.rm = TRUE)
          m=c(m1,m2,m3)
          prop=c(prop1,prop2,prop3)
          i=which(abs(m)==max(abs(m)))

          # i_min=which(abs(m)==min(abs(m))) # marche pas

          inverse=FALSE

          if (m[i]>0){

            i_min = which(m == min(m))
            if (prop[i_min]>2*prop[i]){
              inverse=TRUE
            }

            ordre = order(resClustering[[k]]@parameters@mean[val_clust,1])
            tmp = as.factor(resClustering[[k]]@partition)
            tmp2=rep(NA,length(tmp))
            for (l in 1:3){
              tmp2[tmp==val_clust[ordre[l]]]=l
            }
            if (inverse){
              resGenoAssign[k,]=-as.numeric(as.character(tmp2))+4
            } else {
              resGenoAssign[k,]=-as.numeric(as.character(tmp2))+3 # car 1 doit etre 2 ; 2e doit etre 1 et 3e doit etre 0 (ordre croissant)
            }
          } else{

            i_min = which(m == max(m))
            if (prop[i_min]>2*prop[i]){
              inverse=TRUE
            }

            ordre = order(resClustering[[k]]@parameters@mean[val_clust,1])
            tmp = as.factor(resClustering[[k]]@partition)
            tmp2=rep(NA,length(tmp))
            for (l in 1:3){
              tmp2[tmp==val_clust[ordre[l]]]=l
            }
            if (inverse){
              resGenoAssign[k,]=-as.numeric(as.character(tmp2))+3
            } else {
              resGenoAssign[k,]=-as.numeric(as.character(tmp2))+4
            }
          }
          # if (m[i]>0){
          #   ordre = order(resClustering[[k]]@parameters@mean[val_clust,1])
          #   tmp = as.factor(resClustering[[k]]@partition)
          #   tmp2=rep(NA,length(tmp))
          #   for (l in 1:3){
          #     tmp2[tmp==val_clust[ordre[l]]]=l
          #   }
          #   resGenoAssign[k,]=-as.numeric(as.character(tmp2))+3 # car 1 doit etre 1 ; 2e doit etre 2 et 3e doit etre 3 (ordre croissant)
          # } else{
          #   ordre = order(resClustering[[k]]@parameters@mean[val_clust,1])
          #   tmp = as.factor(resClustering[[k]]@partition)
          #   tmp2=rep(NA,length(tmp))
          #   for (l in 1:3){
          #     tmp2[tmp==val_clust[ordre[l]]]=l
          #   }
          #   resGenoAssign[k,]=-as.numeric(as.character(tmp2))+4
          # }
        } else if (n_clust_restant==2){ # si 2 clusters : valeur la plus extreme homoz, lautre heteroz : le plus proche
          m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
          m2=resClustering[[k]]@parameters@mean[val_clust[2],1]
          if (abs(m1)>abs(m2)){
            if (m1>m2){
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=0
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=1
            } else {
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=3
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=2
            }
          } else {
            if (m1>m2){
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=3
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=2
            } else {
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=0
              resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=1
            }
          }
        } else if (n_clust_restant==1){ # si que 1 : homoz associe (pos ou neg)
          m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
          if (m1>0){
            resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=0
          } else {
            resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=3
          }
        }

        tmp3 = apply(resClustering[[k]]@proba,MARGIN=1,FUN=max)
        resGenoAssign[k,tmp3<SeuilNoCall]=-1
        resGenoAssign[k,is.na(resClustering[[k]]@partition)]=-1
        # On test si ce marker sera garde ou non
        ToKeep = keepMarkertriplo(marker = k,
                                  genotypePop = resGenoAssign[k,],
                                  data = Dataset[Dataset$MarkerName==k,],
                                  cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
        df_classif[k,]=ToKeep
        if (!is.na(df_classif[k,"Message"])){
          if (df_classif[k,"Message"]=='CR' & verif[[2]]=="SD_gp"){
            df_classif[k,"Message"]="CR-SDgp"
          }
        }
        if (n_clust_restant==4){
          list_max_clust=c(list_max_clust,k)
        }
      }
    }
  }
  df_classif = cbind(data.frame(MarkerName=rownames(df_classif)),df_classif)
  return(list(resGenoAssign,list_max_clust,df_classif))
}
