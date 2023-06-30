#' Discriminate markers (diploid version)
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

keepMarkerdiplo = function(marker,genotypePop,data,cr_marker=NULL,fld_marker=NULL,hetso_marker=NULL){
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
  SeuilHomRO=c(0.6,0.3,-0.9)

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

  # Suivant les genotypes, differentes metrics peuvent etre calcules
  if (1 %in% alleles & nClus>1){
    # calcul FLD et HetSO
    if (nClus==3){
      poolVar = ((length(p0)-1)*var(p0)+(length(p1)-1)*var(p1)+(length(p2)-1)*var(p2))/(length(c(p0,p1,p2))-3)
      if (abs(m1-m0)<abs(m1-m2)){ # calcul FLD avec 0 et 1
        FLD = abs(m1-m0)/sqrt(poolVar)
      } else { # calcul FLD avec 1 et 2
        FLD = abs(m1-m2)/sqrt(poolVar)
      }
      HetSO = a1-a2-(a0-a2)*(m1-m2)/(m0-m2)
    } else { # donc nClus==2
      if (0 %in% alleles){
        poolVar = ((length(p0)-1)*var(p0)+(length(p1)-1)*var(p1))/(length(c(p0,p1))-2)
        FLD = abs(m1-m0)/sqrt(poolVar)

        HetSO = a1-a0
      } else { # cest que cest 2
        poolVar = ((length(p1)-1)*var(p1)+(length(p2)-1)*var(p2))/(length(c(p1,p2))-2)
        FLD = abs(m1-m2)/sqrt(poolVar)

        HetSO = a1-a2
      }
    }
  }
  # HomRO
  if (0 %in% alleles & 2 %in% alleles){
    if (m0>0 & m2<0){
      HomRO = min(m0,abs(m2))
    } else if (m0>0 & m2>0){
      HomRO = -m2
    } else if (m0<0 & m2<0){
      HomRO = m0
    }
  } else if (0 %in% alleles){
    HomRO = m0
  } else if (2 %in% alleles){
    HomRO = -m2
  }
  # Z-scores
  # if (0 %in% alleles){
  #   x0=data$Contrast[data$SampleName %in% colnames(genotypePop[which(genotypePop==0)])]
  #   y0=data$SigStren[data$SampleName %in% colnames(genotypePop[which(genotypePop==0)])]
  #
  #   zscoreX0 = abs((x0-mean(x0))/sd(x0)) # max ? je ne sais pas quelle valeur prendre
  # }

  # Calcul MAF
  f1 = length(which(genotypePop==0))/length(genotypePop)
  f2 = length(which(genotypePop==1))/length(genotypePop)
  f3 = length(which(genotypePop==2))/length(genotypePop)
  MAF = min((2*f1+f2)/2,(2*f3+f2)/2)
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
