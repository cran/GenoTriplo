#' Add Categories
#'
#' @param X dataframe with metrics of markers
#' @param ploidy ploidy level of the population
#'
#' @return the same dataframe with a suppleme,tary column indicating the categorie of the marker
#'
#' @keywords internal
#' @noRd

add_categories = function(X,ploidy){
  cate=c()
  for (k in 1:nrow(X)){
    if (is.na(X[k,'Message'])){
      if (X[k,'nClus']==ploidy+1){
        cate=c(cate,"PolyHighResolution")
      } else if (X[k,'nClus']==1){
        cate=c(cate,"MonoHighResolution")
      } else {
        cate=c(cate,"NoMinorHomozygote")
      }
    } else if (regexpr(pattern = "CR",text = X[k,"Message"])!=-1){
      cate=c(cate,"CallRateBelowThreshold")
    } else if (regexpr(pattern = "HetSO",text = X[k,"Message"])!=-1){
      cate=c(cate,"OffTargetVariant")
    } else if (!is.na(X[k,"Message"])){
      cate=c(cate,"Other")
    }
  }
  X$Categorie = cate
  return(X)
}
