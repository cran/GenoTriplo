#' Calcul distance between individuals and its cluster center
#'
#' @param genotype the cluster of each individuals of the population for a given marker
#' @param data dataset with Contrast for each individuals
#' @param SampleName the names of each individuals
#'
#' @importFrom stats sd
#'
#' @return the cluster of each individuals of the population for a given marker modified if necessary
#'
#' @keywords internal
#' @noRd

NbSD_gp = function(genotype,data,SampleName){
  res=rep(NA,length(SampleName)) # creation dun vecteur de resultat de la meme taille que le vecteur de noms (une distance par individu)
  for (l in unique(as.numeric(genotype))){ # on regarde genotype par genotype
    indiv = SampleName[which(genotype==l)] # on stock le nom des individus ayant le genotype l
    sd_marqueur=sd(data$Contrast[data$SampleName %in% indiv]) # on calcule lecart type de ce genotype
    m_marqueur=mean(data$Contrast[data$SampleName %in% indiv])# et sa moyenne

    dist_to_mean_gp_l = abs(data$Contrast[data$SampleName %in% indiv]-m_marqueur) # enfin on calcule la distance entre chacun des individus et la moyenne du culster
    nb_sd_gp_l = dist_to_mean_gp_l/sd_marqueur # on divise par lecart type pour avoir une distance en nombre d'ecart type
    if (is.na(sd_marqueur)){ # si NA (notamment si il y a que 1 individu)
      nb_sd_gp_l = 0 # on met a 0 (si que 1 indiv, il est pas loin de son centre...)
    }
    res[which(genotype==l)] = nb_sd_gp_l # on stock les nb de sd du genotype l dans le vecteur de resultat
  }
  return(res)
}
