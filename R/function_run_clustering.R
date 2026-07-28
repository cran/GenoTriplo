#' Launch parallel clustering from directory
#'
#' Launch the clustering phase in parallel from the dataset with SampleName, Contrast and SigStren for each markers (MarkerName).
#'
#' @param dir_from_create directory name with dataset files saved in './output_create'
#' @param new.dir directory name where output is saved as in './output_clustering'
#' @param ploidy ploidy of individuals
#' @param n_iter number of iterations of clustering
#' @param Dmin threshold distance between two clusters
#' @param n_core number of cores used for parallelization
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
#' \dontrun{
#' dir.create = system.file("extdata/output_create/1", package = "GenoTriplo")
#' Clustering_parallele_from_dir(dir_from_create=dir.create,ploidy=3)
#' }
#' 

Clustering_parallele_from_dir = function(dir_from_create,new.dir=NULL,ploidy,n_iter=5,Dmin=0.28,n_core=1){
  if (is.null(new.dir)){
    new.dir=length(list.dirs("./output_clustering",recursive = F))+1
  }
  if (length(grep(pattern = "output_clustering",x = new.dir))==0){
    dir.name.save = paste0("./output_clustering/",new.dir)
  } else if (length(grep(pattern = "^./",x = new.dir,fixed = F))==0) {
    dir.name.save=paste0("./",new.dir)
  } else {dir.name.save=new.dir}
  if (!dir.exists(dir.name.save)){
    dir.create(dir.name.save,recursive = T)
    cat(paste0("Directory created : ",dir.name.save,"\n"))
  }
  if (!is.numeric(ploidy)){
    stop("ploidy must be a numeric value !")
  }
  if (!is.numeric(n_core)){
    stop("n_core must be numeric.")
  }else if (parallel::detectCores()<n_core){
    n_core=parallel::detectCores()
    warning("The number of core asked is to high : will be set to maximum.")
  }
  t0 = Sys.time()
  # Parametrage des clusters
  cat("Launching clustering...\n")
  cat(paste0("Core used : ",n_core,"\n"))
  clust_name = parallel::makeCluster(n_core)
  doParallel::registerDoParallel(clust_name)
  # Boucle foreach qui permet de paralleliser les taches
  datafiles = list.files(dir_from_create)
  f=NULL
  res_clust = foreach(f=datafiles,.combine = 'c') %dopar% {
    gc()
    # Lance le clustering
    Clustering_from_file(dirname=dir_from_create,filename=f,ploidy=ploidy,
                         n_iter=n_iter,Dmin=Dmin,
                         new.dir=dir.name.save)
  }
  # Arrete des clusters (indispensable)
  parallel::stopCluster(clust_name)
  doParallel::stopImplicitCluster()
  t1=Sys.time()
  delay = t1-t0
  cat(paste0("Total clustering time : ",round(delay,3),units(delay),"\n"))
}

#' Clustering function from directory
#'
#' Clustering function to run clustering with no parallelization process and autosave
#'
#' @param dir_from_create directory of the file to clusterize.
#' @param ploidy ploidy of individuals
#' @param n_iter number of iterations to perform for clustering
#' @param Dmin minimal distance between two clusters
#' @param new.dir new directory name for saving in './output_clustering'
#'
#' @import Rmixmod
#' @importFrom rlang .data
#'
#' @export
#'
#' @return Autosave of the results in new.dir
#' 
#' @examples
#' \dontrun{
#' dir.create = system.file("extdata/output_create/1", package = "GenoTriplo")
#' ploidy = 3
#' Clustering_from_dir(dir_from_create=dir.create,ploidy=ploidy)
#' }

Clustering_from_dir = function(dir_from_create,ploidy,n_iter=5,Dmin=0.28,new.dir=NULL){
  if (!dir.exists("./output_clustering")){
    dir.create("./output_clustering")
    cat(paste0("Directory created : ","./output_clustering","\n"))
  }
  if (is.null(new.dir)){
    new.dir=length(list.dirs("./output_clustering",recursive = F))+1
  }
  if (length(grep(pattern = "output_clustering",x = new.dir))==0){
    new.dir.tot = paste0("./output_clustering/",new.dir)
  } else if (length(grep(pattern = "^./",x = new.dir,fixed = F))==0) {
    new.dir.tot=paste0("./",new.dir)
  } else {new.dir.tot=new.dir}
  if (!dir.exists(new.dir.tot)){
    dir.create(path = new.dir.tot,recursive = TRUE)
    cat(paste0("Directory created : ",new.dir.tot,"\n"))
  }
  all.files=list.files(dir_from_create)
  for (f in all.files){
    cat(paste0("Clustering from file : ",dir_from_create,"/",f," ...\n"))
    Clustering_from_file(dirname=dir_from_create,
                         filename=f,
                         ploidy=ploidy,n_iter=n_iter,Dmin=Dmin,
                         new.dir=new.dir.tot)
  }
}
