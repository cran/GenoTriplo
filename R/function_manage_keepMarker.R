#' Manage the discrimination of markers (apply keepMarkerXXX depending which)
#'
#' @param which NbClustMax (ploidy +1)
#' @param marker marker name
#' @param genotypePop genotype of the population for teh marker
#' @param data dataframe with Contrast and SigStren for each individuals for the given marker
#' @param cr_marker threshold call rate for the marker
#' @param fld_marker thresholf FLD for the marker
#' @param hetso_marker threshold HetSO for the marker
#'
#' @return dataframe with different metrics for the given marker
#'
#' @keywords internal
#' @noRd
#' 
keepMarkerPloidy = function(which,
                            marker,
                            genotypePop,
                            data,
                            cr_marker=NULL,fld_marker=NULL,hetso_marker=NULL){
  if (which==3){
    res=keepMarkerdiplo(marker = marker,
                        genotypePop = genotypePop,
                        data = data,
                        cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
  } else if (which==4){
    res=keepMarkertriplo(marker = marker,
                         genotypePop = genotypePop,
                         data = data,
                         cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
  } else if (which==5){
    res=keepMarkertetra(marker = marker,
                        genotypePop = genotypePop,
                        data = data,
                        cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
  } else {
    stop(paste0("ploidy==",which-1," not supported !"))
  }
  return(res)
}
