#' Calculate standard deviation of clusters
#'
#' @param genotype genotype of the individuals (0,1,2,3) for a given marker
#' @param data dataframe with Contrast & SampleName
#' @param SeuilSD Thresold for standars deviation of a cluster/group
#' @importFrom stats sd
#'
#' @return list with genotype and whether a group has been deleted
#'
#' @keywords internal
#' @noRd
verif_sd=function(genotype,data,SeuilSD=0.28){
  res=genotype # on stock le genotype dans un vecteur qui portera les modifications
  message="" # rien a signaler au depart
  for (i in unique(as.numeric(genotype))){ # pour chacun des clusters (genotypes)
    if (!is.na(i)){
      indiv = data$SampleName[which(genotype==i)] # stock les indivs du genotype i
      if (length(indiv)>1){
        sd_marqueur=sd(data$Contrast[data$SampleName %in% indiv]) # calcul du sd du groupe
        m=mean(data$Contrast[data$SampleName %in% indiv]) # calcul de la moyenne du groupe
        if (sd_marqueur>SeuilSD*(1+0.5*abs(m))){ # comparaison au seuil (multiplication par la deriveCCS de la moyenne pour lisser un peu la courbe sans que les valeurs soient entre -infini:+infini)
          res[which(genotype==i)] = NA # si passe pas le sueil, le groupe entier est misNA
          message="SD_gp" # indique quil y a eu un changement
        }
      }
    }
  }
  return(list(res,message))
}
