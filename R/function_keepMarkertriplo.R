#' Discriminate markers (triploid version)
#'
#' @param marker marker name
#' @param genotypePop genotype of the population for teh marker
#' @param data dataframe with Contrast and SigStren for each individuals for the given marker
#' @param cr_marker threshold call rate for the marker
#' @param fld_marker thresholf FLD for the marker
#' @param hetso_marker threshold HetSO for the marker
#'
#' @importFrom stats var
#'
#' @return dataframe with different metrics for the given marker
#'
#' @keywords internal
#' @noRd

keepMarkertriplo = function(marker,genotypePop,data,cr_marker=NULL,fld_marker=NULL,hetso_marker=NULL){
  if (is.null(data$SigStren)){
    stop("Need SigStren")
  }
  if (is.null(data$Contrast)){
    stop("Need Contrast")
  }
  alleles = unique(as.numeric(genotypePop)) # les differents genotypes identifies pour le marker
  nClus = length(which(alleles != -1)) # le nombre de genotypes identifies
  CR = length(which(genotypePop!=-1))/length(genotypePop) # le callrate pour ce marker
  FLD=NA # initialisation de variable
  HetSO=NA # initialisation de variable
  HomRO=NA # initialisation de variable
  MAF = NA # initialisation de variable

  # Initialisation des differents seuils (voir manuel AXAS)
  SeuilCR=ifelse(test = is.null(cr_marker),yes = 0.97,no=cr_marker)
  SeuilFLD=ifelse(test = is.null(fld_marker),yes = 3.4,no=fld_marker)
  SeuilHetSO=ifelse(test = is.null(hetso_marker),yes = -0.3,no=hetso_marker)
  SeuilHomRO=c(0.6,0.3,0.3,-0.9)

  # Calcul des moyennes en contrast et sigstren pour chaque genotype
  if (0 %in% alleles){
    p0=data$Contrast[data$SampleName %in% colnames(genotypePop[which(genotypePop==0)])]
    m0=mean(p0)
    a0=mean(data$SigStren[data$SampleName %in% colnames(genotypePop[which(genotypePop==0)])])
  }
  if (1 %in% alleles){
    p1=data$Contrast[data$SampleName %in% colnames(genotypePop[which(genotypePop==1)])]
    m1=mean(p1)
    a1=mean(data$SigStren[data$SampleName %in% colnames(genotypePop[which(genotypePop==1)])])
  }
  if (2 %in% alleles){
    p2=data$Contrast[data$SampleName %in% colnames(genotypePop[which(genotypePop==2)])]
    m2=mean(p2)
    a2=mean(data$SigStren[data$SampleName %in% colnames(genotypePop[which(genotypePop==2)])])
  }
  if (3 %in% alleles){
    p3=data$Contrast[data$SampleName %in% colnames(genotypePop[which(genotypePop==3)])]
    m3=mean(p3)
    a3=mean(data$SigStren[data$SampleName %in% colnames(genotypePop[which(genotypePop==3)])])
  }

  # Suivant les genotypes, differentes metrics peuvent etre calcules
  if (nClus==1){
    if (0 %in% alleles){
      HomRO=m0
    } else if (1 %in% alleles){
      HomRO=m1
    } else if (2 %in% alleles){
      HomRO=-m2
    } else if (3 %in% alleles){
      HomRO=-m3
    }
  }

  if (nClus==2){
    if (all(c(0,1) %in% alleles)){
      poolVar = ((length(p0)-1)*var(p0)+(length(p1)-1)*var(p1))/(length(c(p0,p1))-2)
      FLD = abs(m1-m0)/sqrt(poolVar)

      HetSO = a1-a0
      HomRO = m0
    } else if (all(c(0,2) %in% alleles)){
      poolVar = ((length(p0)-1)*var(p0)+(length(p2)-1)*var(p2))/(length(c(p0,p2))-2)
      FLD = abs(m2-m0)/sqrt(poolVar)

      HetSO = a2-a0
      HomRO = m0
    } else if (all(c(0,3) %in% alleles)){
      poolVar = ((length(p0)-1)*var(p0)+(length(p3)-1)*var(p3))/(length(c(p0,p3))-2)
      FLD = abs(m3-m0)/sqrt(poolVar)

      # HetSO = NA # pas de HetSO parce que pas d'heteroz
      if (m0>0 & m3<0){
        HomRO = min(m0,abs(m3))
      } else if (m0>0 & m3>0){
        HomRO = -m3
      } else if (m0<0 & m3<0){
        HomRO = m0
      }
    } else if (all(c(1,2) %in% alleles)){
      poolVar = ((length(p1)-1)*var(p1)+(length(p2)-1)*var(p2))/(length(c(p1,p2))-2)
      FLD = abs(m2-m1)/sqrt(poolVar)

      # HetSO = NA # pas de HetSO parce que pas d'homoz pour comparer
      if (m1>0 & m2<0){
        HomRO = min(m1,abs(m2))
      } else if (m1>0 & m2>0){
        HomRO = -m2
      } else if (m1<0 & m2<0){
        HomRO = m1
      }
    } else if (all(c(1,3) %in% alleles)){
      poolVar = ((length(p3)-1)*var(p3)+(length(p1)-1)*var(p1))/(length(c(p3,p1))-2)
      FLD = abs(m1-m3)/sqrt(poolVar)

      HetSO = a1-a3
      HomRO = -m3
    } else if (all(c(2,3) %in% alleles)){
      poolVar = ((length(p3)-1)*var(p3)+(length(p2)-1)*var(p2))/(length(c(p3,p2))-2)
      FLD = abs(m2-m3)/sqrt(poolVar)

      HetSO = a2-a3
      HomRO = -m3
    }
  }
  if (nClus==3){
    if (all(c(0,1,2) %in% alleles)){
      poolVar = ((length(p0)-1)*var(p0)+(length(p1)-1)*var(p1)+(length(p2)-1)*var(p2))/(length(c(p0,p1,p2))-3)
      if (abs(m1-m0)<abs(m1-m2)){ # calcul FLD avec 0 et 1
        FLD = abs(m1-m0)/sqrt(poolVar)
      } else { # calcul FLD avec 1 et 2
        FLD = abs(m1-m2)/sqrt(poolVar)
      }

      HetSO = min(a1-a0,a2-a0)
      HomRO = m0
    } else if (all(c(1,2,3) %in% alleles)){
      poolVar = ((length(p1)-1)*var(p1)+(length(p2)-1)*var(p2)+(length(p3)-1)*var(p3))/(length(c(p1,p2,p3))-3)
      if (abs(m2-m1)<abs(m2-m3)){ # calcul FLD avec 0 et 1
        FLD = abs(m2-m1)/sqrt(poolVar)
      } else { # calcul FLD avec 1 et 2
        FLD = abs(m2-m3)/sqrt(poolVar)
      }

      HetSO = min(a1-a3,a2-a3)
      HomRO = -m3
    } else if (all(c(0,1,3) %in% alleles)){
      poolVar = ((length(p0)-1)*var(p0)+(length(p1)-1)*var(p1)+(length(p3)-1)*var(p3))/(length(c(p0,p1,p3))-3)
      if (abs(m1-m0)<abs(m1-m3)){ # calcul FLD avec 0 et 1
        FLD = abs(m1-m0)/sqrt(poolVar)
      } else { # calcul FLD avec 1 et 2
        FLD = abs(m1-m3)/sqrt(poolVar)
      }

      HetSO = a1-a3-(a0-a3)*(m1-m3)/(m0-m3) # pas sur de la formule. Voir ce que ?a donne.

      if (m0>0 & m3<0){
        HomRO = min(m0,abs(m3))
      } else if (m0>0 & m3>0){
        HomRO = -m3
      } else if (m0<0 & m3<0){
        HomRO = m0
      }
    } else if (all(c(0,2,3) %in% alleles)){
      poolVar = ((length(p0)-1)*var(p0)+(length(p2)-1)*var(p2)+(length(p3)-1)*var(p3))/(length(c(p0,p2,p3))-3)
      if (abs(m2-m0)<abs(m2-m3)){ # calcul FLD avec 0 et 1
        FLD = abs(m2-m0)/sqrt(poolVar)
      } else { # calcul FLD avec 1 et 2
        FLD = abs(m2-m3)/sqrt(poolVar)
      }

      HetSO = a2-a3-(a0-a3)*(m2-m3)/(m0-m3) # pas sur de la formule. Voir ce que ?a donne.

      if (m0>0 & m3<0){
        HomRO = min(m0,abs(m3))
      } else if (m0>0 & m3>0){
        HomRO = -m3
      } else if (m0<0 & m3<0){
        HomRO = m0
      }
    }
  }
  if (nClus==4){
    poolVar = ((length(p0)-1)*var(p0)+(length(p1)-1)*var(p1)+(length(p2)-1)*var(p2)+(length(p3)-1)*var(p3))/(length(c(p0,p1,p2,p3))-4)
    if (abs(m1-m0)<abs(m2-m1) & abs(m1-m0)<abs(m3-m2)){ # calcul FLD avec 0 et 1
      FLD = abs(m1-m0)/sqrt(poolVar)
    } else if (abs(m2-m1)<abs(m1-m0) & abs(m2-m1)<abs(m3-m2)){ # calcul FLD avec 1 et 2
      FLD = abs(m2-m1)/sqrt(poolVar)
    } else if (abs(m3-m2)<abs(m1-m0) & abs(m3-m2)<abs(m2-m1)){
      FLD = abs(m3-m2)/sqrt(poolVar)
    }

    if (a1<a2){
      HetSO = a1-a3-(a0-a3)*(m1-m3)/(m0-m3) # pas sur de la formule. Voir ce que ?a donne.
    } else {
      HetSO = a2-a3-(a0-a3)*(m2-m3)/(m0-m3) # pas sur de la formule. Voir ce que ?a donne.
    }


    if (m0>0 & m3<0){
      HomRO = min(m0,abs(m3))
    } else if (m0>0 & m3>0){
      HomRO = -m3
    } else if (m0<0 & m3<0){
      HomRO = m0
    }
  }
  # Calcul MAF
  f1 = length(which(genotypePop==0))/length(genotypePop)
  f2 = length(which(genotypePop==1))/length(genotypePop)
  f3 = length(which(genotypePop==2))/length(genotypePop)
  f4 = length(which(genotypePop==3))/length(genotypePop)
  MAF = min((3*f1+2*f2+f3)/3,(3*f4+2*f3+f2)/3)
  # Stockage des metrics dans un df
  res=data.frame(toKeep=TRUE,CR=round(CR,2),FLD=round(FLD,2),
                 HetSO=round(HetSO,2),HomRO=round(HomRO,2),
                 nClus=nClus,MAF=round(MAF,2),Message=NA)

  # Test contre les seuils : si un seuil ne passe pas, toKeep prend la valeur FALSE et on rajoute au message le parametre qui na pas passe le seuil
  if (!is.na(CR)){
    if (CR < SeuilCR){
      res$toKeep=FALSE
      if (is.na(res$Message)){
        res$Message="CR"
      } else {
        res$Message = paste(res$Message,"CR")
      }
    }
  }
  if (!is.na(FLD)){
    if (FLD < SeuilFLD){
      res$toKeep=FALSE
      if (is.na(res$Message)){
        res$Message="FLD"
      } else {
        res$Message = paste(res$Message,"FLD")
      }
    }
  }
  if (!is.na(HetSO)){
    if (HetSO < SeuilHetSO){
      res$toKeep=FALSE
      if (is.na(res$Message)){
        res$Message="HetSO"
      } else {
        res$Message = paste(res$Message,"HetSO")
      }
    }
  }
  if (!is.na(HomRO)){
    SeuilHomRO[nClus]
    if (HomRO < SeuilHomRO[nClus]){
      res$toKeep=FALSE
      if (is.na(res$Message)){
        res$Message="HomRO"
      } else {
        res$Message = paste(res$Message,"HomRO")
      }
    }
  }
  return(res)
}
