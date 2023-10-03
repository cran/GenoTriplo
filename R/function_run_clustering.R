#' Launch parallel clustering
#'
#' Launch the clustering phase in parallel from the dataset with SampleName, Contrast and SigStren for each markers.
#'
#' @param data_clustering dataframe result from create dataset phase
#' @param ploidy ploidy of offspring
#' @param save_n name of the saving file
#' @param n_iter number of iterations of clustering
#' @param D_min threshold distance between two clusters
#' @param n_core number of cores used for parallelization
#' @param path_log path for log file when run by the shiny app
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom rlang .data
#' @import foreach
#'
#' @return the result of clustering or automatically save a list of objects if a saving name has been provided
#'
#' @export
#'
#' @examples
#'
#' data(GenoTriplo_to_clust)
#' res = Run_Clustering(data_clustering=GenoTriplo_to_clust,
#'                      ploidy=3,n_iter=5,n_core=1)
#' # or if you want to automatically save the result
#' # This will automatically create a folder and save the result in it
#' # Run_Clustering(data_clustering=GenoTriplo_to_clust,
#' #                ploidy=3,n_iter=5,n_core=1,save_n='exemple')
#'
#'

Run_Clustering = function(data_clustering,ploidy,save_n='',n_iter=5,D_min=0.28,n_core=1,path_log=''){

  if (!is.numeric(ploidy)){
    stop("ploidy must be a numeric value !")
  } else if (ploidy>3 | ploidy<2){
    stop("ploidy must be 2 or 3 !")
  }
  if (!is.numeric(n_core)){
    stop("n_core must be numeric.")
  }else if (parallel::detectCores()<n_core){
    n_core=parallel::detectCores()
    warning("The number of core asked is to high : will be set to maximum.")
  }

  nb_indiv=length(unique(data_clustering$SampleName))
  if (nb_indiv<50 & save_n!='' & path_log!=''){
    write(x = "Warning : The number of individuals might not be enough.",file=path_log,append=T)
  } else if (nb_indiv<50 & save_n==''){
    warning("The number of individuals might not be enough.")
  }

  MarkerId = unique(data_clustering$MarkerName) # stockage des noms des differents marqueurs
  nb_of_marker = length(MarkerId) # stockage du nombre total de marker
  batchsize = nb_of_marker%/%n_core
  i_max = nb_of_marker %/% batchsize # quit a ce que le dernier batch soit un peu plus grand (de quelques marker)
  if (batchsize>2000){
    batchsize = batchsize%/%2
    i_max = 2*i_max
    while (batchsize>2000){
      batchsize = batchsize%/%2
      i_max = 2*i_max
    }
  }

  if (save_n!='' & path_log!=''){
    write(x = paste0("Batchsize : ",batchsize," ; Nb Batch : ",i_max),file = path_log,append=T)
  }
  t0 = Sys.time()
  # Parametrage des clusters
  clust_name = parallel::makeCluster(min(n_core,i_max))
  doParallel::registerDoParallel(clust_name)
  # Boucle foreach qui permet de paralleliser les taches
  i=NULL
  res_clust = foreach(i=1:i_max,.combine = 'c') %dopar% {
    gc()
    # Load packages and functions
    if (i!=i_max){ # petite particularite si i == imax => le dernier batch ne fait pas forcement la taille maximum de batchsize
      # Lance le clustering
      clust=Clustering(dataset = data_clustering[data_clustering$MarkerName %in% MarkerId[(1+(i-1)*batchsize):(i*batchsize)],],
                       nb_clust_possible = ploidy+1,
                       n_iter = n_iter,Dmin = D_min)
    } else {
      # Lance le clustering avec le dernier batch qui va en realiste jusquau dernier marker (au cas ou le batch ne soit pas complet)
      clust=Clustering(dataset = data_clustering[data_clustering$MarkerName %in% MarkerId[(1+(i-1)*batchsize):nb_of_marker],],
                       nb_clust_possible = ploidy+1,
                       n_iter = n_iter,Dmin = D_min)
    }
    clust
  }
  # Arrete des clusters (indispensable)
  parallel::stopCluster(clust_name)
  doParallel::stopImplicitCluster()
  t1=Sys.time()
  delay = t1-t0
  if (!dir.exists("./output_clustering") & save_n!=''){
    dir.create("./output_clustering")
  }
  if (save_n!=''){
    save(res_clust,data_clustering,delay,file=paste0("./output_clustering/",save_n,"_to_geno.Rdata"))
  } else {
    return(res_clust)
  }
}

