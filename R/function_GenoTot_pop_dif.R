#' Genotype clusters for individuals from different populations
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

GenoAssign_pop_dif = function(resClustering,SampleName,NbClustMax,SeuilNoCall,SeuilSD,SeuilNbSD,Dataset,cr_marker=NULL,fld_marker=NULL,hetso_marker=NULL,verbose=FALSE){
  MarkerName = names(resClustering) # Stockage des noms des marker
  resGenoAssign = as.data.frame(matrix(nrow = length(MarkerName),ncol=length(SampleName))) # creation df de resultats (avant autant de ligne que de marqueur et de colonne que d'individu)
  names(resGenoAssign)=SampleName # les colonnes (variables) prennent le noms des individus
  rownames(resGenoAssign)=MarkerName # les lignes prennet le nom des marqueurs
  list_max_clust=c() # liste stock les numeros des marker avec un nombre de genotype maximum
  list_max_clust_false=c() # pareil mais ceux qui ont perdu un cluster apres verification (pour eviter quils passent 2 fois dans les boucles)
  df_classif=data.frame(toKeep=rep(NA,length(MarkerName)),CR=NA,FLD=NA,HetSO=NA,HomRO=NA,nClus=NA,MAF=NA,Message=NA) # creation df de classification des marqueurs (pour trier selon des criteres)
  rownames(df_classif)=MarkerName # les lignes prennent le nom des marqueurs
  if (NbClustMax==3){ # on regarde le nombre de genotype maximum : ici diploide
    if (is.null(Dataset$MarkerName)){
      stop("Your dataset must contain MarkerName as variable.")
    }
    Mean_max_clust=data.frame(Gp0=NA,Gp1=NA,Gp2=NA,P0=NA,P1=NA,P2=NA) # creation dun df qui contiendra les moyenne et proportion des genotype pour chaque marqueur ou il y a le nombre maximum de genotype
    for (k in MarkerName){ # on trouve deja ceux avec 3 groupes : assignation automatique + moyenne pour les suivants
      # On verifie quil ny a pas dindividu trop eloignes sinon on les mets en NoCall
      nb_sd_indiv = NbSD_gp(genotype = resClustering[[k]]@partition,
                            data = Dataset[Dataset$MarkerName==k,],
                            SampleName=SampleName)
      if (length(which(nb_sd_indiv>SeuilNbSD))>0){
        resClustering[[k]]@partition=as.factor(resClustering[[k]]@partition)
        resClustering[[k]]@partition[nb_sd_indiv>SeuilNbSD] = NA # individus trop eloignes mis en NA
      }
      # On verifie malgre la suppression des individus trop eloignes quun gp n'est pas trop etale
      verif = verif_sd(genotype = resClustering[[k]]@partition,
                       data = Dataset[Dataset$MarkerName==k,],
                       SeuilSD=SeuilSD)
      resClustering[[k]]@partition=as.factor(verif[[1]])

      if (verbose){print(paste0("Marker : ",k))}
      if (! is.character(resClustering[[k]])){
        val_clust = unique(resClustering[[k]]@partition[which(!is.na(resClustering[[k]]@partition))]) # stockage des numeros de cluster restant
        n_clust_restant = length(val_clust) # calcul du nombre de cluster restant
        if (n_clust_restant==3){ # on regarde si il y a 3 clusters (le nombre maximum en diploide)
          # Les genotypes sont ordonnes pour que le bon numero leur soient assigne
          # Verification du SeuilNoCall => donne valeur -1 (qui sera transforme en NA par la suite)
          ordre = order(resClustering[[k]]@parameters@mean[,1])
          tmp = as.factor(resClustering[[k]]@partition)
          tmp2=rep(NA,length(tmp))
          for (l in 1:3){
            tmp2[tmp==ordre[l]]=l
          }
          resGenoAssign[k,]=-as.numeric(as.character(tmp2))+4-1 # y=ax+b avec y=1=>x=3;y=2=>x=2;y=1=>x=1 puis -1 car valeur = 0,1,2 et non 1,2,3
          tmp3 = apply(resClustering[[k]]@proba,MARGIN=1,FUN=max)
          resGenoAssign[k,tmp3<SeuilNoCall]=-1
          resGenoAssign[k,is.na(resClustering[[k]]@partition)]=-1

          # On test si ce marker sera garde ou non
          ToKeep = keepMarkerdiplo(marker = k,
                                genotypePop = resGenoAssign[k,],
                                data = Dataset[Dataset$MarkerName==k,],
                                cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
          df_classif[k,]=ToKeep # Stockage du resultat
          if (!is.na(df_classif[k,"Message"])){
            if (df_classif[k,"Message"]=='CR' & verif[[2]]=="SD_gp"){
              df_classif[k,"Message"]="CR-SDgp" # changement dun des messages de rejet suivant une condition (si trop de NA et du a un gp avec trop grand sd)
            }
          }
          if (ToKeep$toKeep){ # si ce marker est bon
            list_max_clust=c(list_max_clust,k) # on stock le fait que ce marker soit bon
            Mean_max_clust[length(list_max_clust),]=c(resClustering[[k]]@parameters@mean[ordre], # on garde ses valeurs de centre de cluster comme reference
                                                      resClustering[[k]]@parameters@proportions[ordre])
          } else { # si marker pas bon
            list_max_clust_false=c(list_max_clust_false,k) # on stock le fait quil ne soit pas bon (mais on stock quand meme pour pas quil repasse dans la prochaine boucle)
          }


        }
      }
    }

    # Calcul des moyennes de reference des differents genotypes possibles (avec prise en compte de leur proportion au sein de leur proche marker)
    Mu_gp0 = sum(Mean_max_clust$Gp0*Mean_max_clust$P0)/sum(Mean_max_clust$P0)
    Mu_gp1 = sum(Mean_max_clust$Gp1*Mean_max_clust$P1)/sum(Mean_max_clust$P1)
    Mu_gp2 = sum(Mean_max_clust$Gp2*Mean_max_clust$P2)/sum(Mean_max_clust$P2)
    Mu_tot=c(Mu_gp0,Mu_gp1,Mu_gp2)
    if (is.na(Mu_tot[1])){
      stop("Changer les parametres de restrictions : pas de marker polyhigh")
    }
    # Nouvel boucle pour genotyper les marker avec moins que le nombre maximum de genotype (cest pourquoi on enleve les marker presents dans les deux list_max_clust(_false))
    for (k in MarkerName[!MarkerName %in% c(list_max_clust,list_max_clust_false)]){ # mtn quon a les moyennes, on peut assigner le reste
      if (verbose){print(paste0("Marker : ",k))}
      if (! is.character(resClustering[[k]])){ # Verification que le clustering a bien fonctionne pour ce marker
        # Comme plus haut, on verifie le nombre de cluster restant apres avoir enlever les potentiels individus trop eloigne et les groupes trop etendus
        val_clust = unique(resClustering[[k]]@partition[which(!is.na(resClustering[[k]]@partition))])
        n_clust_restant = length(val_clust)
        if (n_clust_restant==1){ # si il reste plus quun groupe
          m1=resClustering[[k]]@parameters@mean[val_clust[1],1] # stockage de son centre
          D1=abs(Mu_tot-m1) # stockage des distance aux moyennes de reference
          g1=which(D1==min(D1)) # groupe (genotype) de reference le plus proche
          resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=-g1+4-1 # assignation de ce groupe
        } else if (n_clust_restant==2){
          # reproduction du chemin pour 2 groupe mais cette fois ci, si les 2 groupes pointes vers le meme groupe de reference,
          # lassignation se fait suivant le groupe de reference designe (si cest l'omoz negatif, le cluster l eplus negatif devient homoz et lautre heteroz,...)
          m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
          m2=resClustering[[k]]@parameters@mean[val_clust[2],1]
          D1=abs(Mu_tot-m1)
          D2=abs(Mu_tot-m2)
          g1=which(D1==min(D1))
          g2=which(D2==min(D2))
          if (g1==g2){ # si meme groupe plus proche -> NoCall pour le gp le + eloigne
            # resGenoAssign[k,]=-1
            if (g1==1){
              if (m1<m2){#D1 est donc le vrai g1
                g2=g2+1
              } else {
                g1=g1+1
              }
            } else if (g1==3){
              if (m1>m2){#D1 est donc le vrai g4
                g2=g2-1
              } else {
                g1=g1-1
              }
            } else if (g1==2){ # si heteroz
              if (D1[g1]<D2[g2]){ # si clust 1 est plus proche de 2 que clust 2 ne lest
                if (m2<Mu_tot[2]){ # si plutot negatif : prend homoz negatif
                  g2=g2-1
                } else {
                  g2=g2+1
                }
              } else { # si clust 2 est plus proche
                if (m1<Mu_tot[2]){ # si plutot negatif : prend homoz negatif
                  g1=g1-1
                } else {
                  g1=g1+1
                }
              }
            }
          }
          resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=-g1+4-1 # y=ax+b avec y=1=>x=3;y=2=>x=2;y=1=>x=1 puis -1 car valeur = 0,1,2 et non 1,2,3
          resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=-g2+4-1
        }
        tmp = apply(resClustering[[k]]@proba,MARGIN=1,FUN=max)
        resGenoAssign[k,tmp<SeuilNoCall]=-1
        resGenoAssign[k,is.na(resClustering[[k]]@partition)]=-1
      } else {
        resGenoAssign[k,]=-1
      }
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
    }
  } else if (NbClustMax==4){ # on refait pareil mais pour les triploides
    if (is.null(Dataset$MarkerName)){
      stop("Your dataset must contain MarkerName as variable.")
    }
    Mean_max_clust=data.frame(Gp0=NA,Gp1=NA,Gp2=NA,Gp3=NA,P0=NA,P1=NA,P2=NA,P3=NA)
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
      # print(k)
      val_clust = unique(resClustering[[k]]@partition[which(!is.na(resClustering[[k]]@partition))])
      n_clust_restant = length(val_clust)
      if (n_clust_restant==4){
        ordre = order(resClustering[[k]]@parameters@mean[,1])
        tmp = as.factor(resClustering[[k]]@partition)
        tmp2=rep(NA,length(tmp))
        for (l in 1:4){
          tmp2[tmp==ordre[l]]=l
        }
        resGenoAssign[k,]=-as.numeric(as.character(tmp2))+5-1
        tmp3 = apply(resClustering[[k]]@proba,MARGIN=1,FUN=max)
        resGenoAssign[k,tmp3<SeuilNoCall]=-1
        resGenoAssign[k,is.na(resClustering[[k]]@partition)]=-1
        if (!is.null(Dataset)){
          # On test si on garde ou non
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
          if (ToKeep$toKeep){
            list_max_clust=c(list_max_clust,k)
            Mean_max_clust[length(list_max_clust),]=c(resClustering[[k]]@parameters@mean[ordre],
                                                      resClustering[[k]]@parameters@proportions[ordre])
          } else {
            list_max_clust_false=c(list_max_clust_false,k)
          }

        } else {
          list_max_clust=c(list_max_clust,k)
          Mean_max_clust[length(list_max_clust),]=c(resClustering[[k]]@parameters@mean[ordre],
                                                    resClustering[[k]]@parameters@proportions[ordre]) # on save les mean des cluster pour faire une moyenne generale
        }
      }
    }

    Mu_gp0 = sum(Mean_max_clust$Gp0*Mean_max_clust$P0)/sum(Mean_max_clust$P0)
    Mu_gp1 = sum(Mean_max_clust$Gp1*Mean_max_clust$P1)/sum(Mean_max_clust$P1)
    Mu_gp2 = sum(Mean_max_clust$Gp2*Mean_max_clust$P2)/sum(Mean_max_clust$P2)
    Mu_gp3 = sum(Mean_max_clust$Gp3*Mean_max_clust$P3)/sum(Mean_max_clust$P3)
    Mu_tot=c(Mu_gp0,Mu_gp1,Mu_gp2,Mu_gp3)
    if (is.na(Mu_tot[1])){
      stop("Changer les parametres de restrictions : pas de marker polyhigh")
    }
    for (k in MarkerName[!MarkerName %in% c(list_max_clust,list_max_clust_false)]){ # mtn quon a les moyennes, on peut assigner le reste
      # print(k)
      val_clust = unique(resClustering[[k]]@partition[which(!is.na(resClustering[[k]]@partition))])
      n_clust_restant = length(val_clust)
      if (n_clust_restant==1){
        m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
        D1=abs(Mu_tot-m1)
        g1=which(D1==min(D1))
        resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=-g1+5-1
      } else if (n_clust_restant==2){
        m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
        m2=resClustering[[k]]@parameters@mean[val_clust[2],1]
        D1=abs(Mu_tot-m1)
        D2=abs(Mu_tot-m2)
        g1=which(D1==min(D1))
        g2=which(D2==min(D2))
        if (g1==g2){ # si meme groupe plus proche -> NoCall pour le gp le + eloigne
          if (g1==1){
            if (m1<m2){#D1 est donc le vrai g1
              g2=g2+1
            } else {
              g1=g1+1
            }
          } else if (g1==4){
            if (m1>m2){#D1 est donc le vrai g4
              g2=g2-1
            } else {
              g1=g1-1
            }
          } else if (g1==2){
            if (D1[g1]<D2[g2]){ # si clust 1 est plus proche de 2 que clust 2 ne lest
              if (m2<Mu_tot[2]){ # si plutot negatif : prend homoz negatif
                g2=g2-1
              } else {
                g2=g2+1
              }
            } else { # si clust 2 est plus proche
              if (m1<Mu_tot[2]){ # si plutot negatif : prend homoz negatif
                g1=g1-1
              } else {
                g1=g1+1
              }
            }
          } else if (g1==3){
            if (D1[g1]<D2[g2]){ # si clust 1 est plus proche de 2 que clust 2 ne lest
              if (m2<Mu_tot[3]){ # si plutot negatif : prend homoz negatif
                g2=g2-1
              } else {
                g2=g2+1
              }
            } else { # si clust 2 est plus proche
              if (m1<Mu_tot[3]){ # si plutot negatif : prend homoz negatif
                g1=g1-1
              } else {
                g1=g1+1
              }
            }
          }
        }
        resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=-g1+5-1 # y=ax+b avec y=1=>x=3;y=2=>x=2;y=1=>x=1 puis -1 car valeur = 0,1,2 et non 1,2,3
        resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=-g2+5-1
      } else if (n_clust_restant==3){
        m1=resClustering[[k]]@parameters@mean[val_clust[1],1]
        m2=resClustering[[k]]@parameters@mean[val_clust[2],1]
        m3=resClustering[[k]]@parameters@mean[val_clust[3],1]
        D1=abs(Mu_tot-m1)
        D2=abs(Mu_tot-m2)
        D3=abs(Mu_tot-m3)
        g1=which(D1==min(D1))
        g2=which(D2==min(D2))
        g3=which(D3==min(D3))
        if (g1==g2 & g1==g3){ # on met tout NA
          g1=5
          g2=5
          g3=5
        } else if (g1==g2){
          if (g1==1){
            if (m1<m2){#D1 est donc le vrai g1
              g2=g2+1
              if (g2==g3){
                g3=g3+1
              }
            } else {
              g1=g1+1
              if (g1==g3){
                g3=g3+1
              }
            }
          } else if (g1==4){
            if (m1>m2){#D1 est donc le vrai g4
              g2=g2-1
              if (g2==g3){
                g3=g3-1
              }
            } else {
              g1=g1-1
              if (g1==g3){
                g3=g3-1
              }
            }
          } else if (g1==2){
            if (g3==1){
              if (m1<m2){ #g1 reste 2 et g2 passe 3
                g2=g2+1
              } else {
                g1=g1+1
              }
            } else if (g3==3){
              if (m1<m2){ #g1 est donc peut etre trop neg
                g1=g1-1
              } else {
                g2=g2-1
              }
            } else if (g3==4){
              if (m1<m2){
                g2=g2+1
              } else {
                g1=g1+1
              }
            }
          } else if (g1==3){
            if (g3==1){
              if (m1<m2){ #g1 reste 2 et g2 passe 3
                g1=g1-1
              } else {
                g2=g2-1
              }
            } else if (g3==2){
              if (m1<m2){ #g1 est donc peut etre trop neg
                g2=g2+1
              } else {
                g1=g1+1
              }
            } else if (g3==4){
              if (m1<m2){
                g1=g1-1
              } else {
                g2=g2-1
              }
            }
          }
        } else if (g2==g3){
          if (g2==1){
            if (m3<m2){#D1 est donc le vrai g1
              g2=g2+1
              if (g2==g1){
                g1=g1+1
              }
            } else {
              g3=g3+1
              if (g3==g1){
                g1=g1+1
              }
            }
          } else if (g2==4){
            if (m3>m2){#D1 est donc le vrai g4
              g2=g2-1
              if (g2==g1){
                g1=g1-1
              }
            } else {
              g3=g3-1
              if (g3==g1){
                g1=g1-1
              }
            }
          } else if (g2==2){
            if (g1==1){
              if (m2<m3){ #g1 reste 2 et g2 passe 3
                g3=g3+1
              } else {
                g2=g2+1
              }
            } else if (g1==3){
              if (m2<m3){ #g1 est donc peut etre trop neg
                g2=g2-1
              } else {
                g3=g3-1
              }
            } else if (g1==4){
              if (m2<m3){
                g3=g3+1
              } else {
                g2=g2+1
              }
            }
          } else if (g2==3){
            if (g1==1){
              if (m2<m3){ #g1 reste 2 et g2 passe 3
                g2=g2-1
              } else {
                g3=g3-1
              }
            } else if (g1==2){
              if (m2<m3){ #g1 est donc peut etre trop neg
                g3=g3+1
              } else {
                g2=g2+1
              }
            } else if (g1==4){
              if (m2<m3){
                g2=g2-1
              } else {
                g3=g3-1
              }
            }
          }
        } else if (g1==g3){
          if (g1==1){
            if (m1<m3){#D1 est donc le vrai g1
              g3=g3+1
              if (g3==g2){
                g2=g2+1
              }
            } else {
              g1=g1+1
              if (g1==g2){
                g2=g2+1
              }
            }
          } else if (g1==4){
            if (m1>m3){#D1 est donc le vrai g4
              g3=g3-1
              if (g3==g2){
                g2=g2-1
              }
            } else {
              g1=g1-1
              if (g1==g2){
                g2=g2-1
              }
            }
          } else if (g1==2){
            if (g2==1){
              if (m1<m3){ #g1 reste 2 et g2 passe 3
                g3=g3+1
              } else {
                g1=g1+1
              }
            } else if (g2==3){
              if (m1<m3){ #g1 est donc peut etre trop neg
                g1=g1-1
              } else {
                g3=g3-1
              }
            } else if (g2==4){
              if (m1<m3){
                g3=g3+1
              } else {
                g1=g1+1
              }
            }
          } else if (g1==3){
            if (g2==1){
              if (m1<m3){ #g1 reste 2 et g2 passe 3
                g1=g1-1
              } else {
                g3=g3-1
              }
            } else if (g2==2){
              if (m1<m3){ #g1 est donc peut etre trop neg
                g3=g3+1
              } else {
                g1=g1+1
              }
            } else if (g1==4){
              if (m1<m2){
                g1=g1-1
              } else {
                g2=g2-1
              }
            }
          }
        }
        resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[1])]=-g1+5-1 # y=ax+b avec y=1=>x=3;y=2=>x=2;y=1=>x=1 puis -1 car valeur = 0,1,2 et non 1,2,3
        resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[2])]=-g2+5-1
        resGenoAssign[k,which(resClustering[[k]]@partition==val_clust[3])]=-g3+5-1
      }
      tmp = apply(resClustering[[k]]@proba,MARGIN=1,FUN=max)
      resGenoAssign[k,tmp<SeuilNoCall]=-1
      resGenoAssign[k,is.na(resClustering[[k]]@partition)]=-1

      if (!is.null(Dataset)){
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
      }
    }
  }
  df_classif = cbind(data.frame(MarkerName=rownames(df_classif)),df_classif)
  return(list(resGenoAssign,list_max_clust,df_classif,Mu_tot))
}
