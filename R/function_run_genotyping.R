#' Launch genotyping phase in parallel from directory
#'
#' Function that launch the genotyping phase from the dataset with SampleName, Contrast and SigStren for each markers and the result of the 'Run_clustering' function.
#'
#' @param dir_from_clustering directory with files from clustering phase
#' @param new.dir name of the new directory for saving results in './output_genotyping'
#' @param ploidy ploidy of individuals
#' @param SeuilNoCall threshold of the probability of belonging to a cluster
#' @param SeuilNbSD threshold for the distance between an individuals and his cluster (x=Contrast)
#' @param SeuilSD threshold for the standard deviation of a cluster (SeuilSD*(1+0.5*abs(mean_contrast_cluster)))
#' @param n_core number of cores used for parallelization
#' @param corres_ATCG dataframe with the correspondence between A/B of AXAS and A/T/C/G (three columns : probeset_id, Allele_A, Allele_B)
#' @param same.pop Boolean : are individuals from a same population
#' @param cr_marker call rate threshold
#' @param fld_marker FLD threshold
#' @param hetso_marker HetSO threshold
#' @param tryEqHW Only if same.pop==TRUE else ignore. If TRUE, will delete cluster that shouldn't happen under H-W equilibrium (heterozygote with less individuals than both homozygote in triploid for example). Additionnaly, the cluster should contain less than 5% of the individuals. Else it stays.

#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import foreach
#' @import dplyr
#'
#' @return Void. Save the results.
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' ploidy=3
#' dir_from_clustering=system.file("extdata/output_clustering/1", package = "GenoTriplo")
#' Genotyping_parallele_from_dir(dir_from_clustering=dir_from_clustering,ploidy=ploidy)
#' }
#' 
Genotyping_parallele_from_dir = function(dir_from_clustering,new.dir=NULL,
                                         ploidy,
                                         SeuilNoCall=0.85,SeuilNbSD=2.8,SeuilSD=0.28,
                                         n_core=1,
                                         corres_ATCG=NULL,
                                         same.pop=TRUE,cr_marker=0.97,fld_marker=3.4,
                                         hetso_marker=-0.3,tryEqHW=TRUE){
  # n_core = parallel::detectCores() - 2
  if (is.null(new.dir)){
    new.dir=length(list.dirs("./output_genotyping",recursive = F))+1
  }
  if (length(grep(pattern = "output_genotyping",x = new.dir))==0){
    dir.name.save = paste0("./output_genotyping/",new.dir)
  } else if (length(grep(pattern = "^./",x = new.dir,fixed = F))==0) {
    dir.name.save=paste0("./",new.dir)
  } else {dir.name.save=new.dir}
  if (!dir.exists(dir.name.save)){
    dir.create(dir.name.save,recursive = T)
    cat(paste0("Directory created : ",new.dir,"\n"))
  }
  if (!is.numeric(ploidy)){
    stop("ploidy must be numeric.")
  } else if (ploidy !=4 & ploidy!=3 & ploidy!=2){
    stop("ploidy must be 2 or 3 or 4!")
  }
  if (!is.numeric(n_core)){
    stop("n_core must be numeric.")
  }else if (parallel::detectCores()<n_core){
    n_core=parallel::detectCores()
    warning("The number of core asked is to high : will be set to maximum.")
  }
  t0 = Sys.time()
  # Parametrage des clusters
  cat("Launching genotyping\n")
  cat(paste0("Core used : ",n_core,"\n"))
  clust_name = parallel::makeCluster(n_core)
  doParallel::registerDoParallel(clust_name)
  datafiles = list.files(dir_from_clustering)
  f=NULL
  # Boucle foreach qui permet de paralleliser les taches
  foreach(f=datafiles) %dopar% {
    Genotyping_from_file(dirname=dir_from_clustering,filename=f,new.dir=dir.name.save,
                         SeuilNoCall = SeuilNoCall,SeuilNbSD = SeuilNbSD,SeuilSD = SeuilSD,
                         ploidy = ploidy,
                         cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker,
                         same.pop=same.pop,tryEqHW = tryEqHW)
  }
  # Arrete des clusters (indispensable)
  parallel::stopCluster(clust_name)
  doParallel::stopImplicitCluster()
  t1=Sys.time()
  delay = t1-t0
  cat(paste0("Total genotyping time : ",round(delay,3)," ",units(delay),"\n"))
  out = create_genofiles(dir.name.save=dir.name.save,new.dir=new.dir,ploidy=ploidy,corres_ATCG=corres_ATCG)
  return(out)
}

#' Run genotyping from file
#'
#' @param dirname name of the directory where the file is 
#' @param filename name of the file with the data for genotyping
#' @param new.dir new directory name for saving
#' @param ploidy ploidy of individuals
#' @param SeuilNoCall threshold of the probability of belonging to a cluster
#' @param SeuilNbSD threshold for the distance between an individuals and his cluster (x=Contrast)
#' @param SeuilSD threshold for the standard deviation of a cluster (`SeuilSD*(1+0.5*abs(mean_contrast_cluster))`)
#' @param cr_marker call rate threshold
#' @param fld_marker FLD threshold
#' @param hetso_marker HetSO threshold
#' @param same.pop Boolean : are individuals from a same population
#' @param tryEqHW Only if same.pop==TRUE else ignore. If TRUE, will delete cluster that shouldn't happen under H-W equilibrium (heterozygote with less individuals than both homozygote in triploid for example). Additionnaly, the cluster should contain less than 5% of the individuals. Else it stays.
#'
#' @returns Void. Save the results
#' @export
#' 
#' @examples
#' \dontrun{
#' ploidy=3
#' dir_from_clustering=system.file("extdata/output_clustering/1", package = "GenoTriplo")
#' file_from_clustering=system.file(
#'                       "extdata/output_clustering/1",
#'                       "1_1_to_geno.Rdata", package = "GenoTriplo")
#' Genotyping_from_file(dirname=dir_from_clustering,ploidy=ploidy,filename=file_from_clustering)
#' }
Genotyping_from_file = function(dirname,filename,new.dir=NULL,
                                ploidy,
                                SeuilNoCall=0.85,SeuilNbSD=2.8,SeuilSD=0.28,
                                cr_marker=0.97,fld_marker=3.4,hetso_marker=-0.3,
                                same.pop=TRUE,tryEqHW=TRUE){
  if (is.null(new.dir)){
    dir.nb = length(list.dirs("./output_genotyping",recursive = F))+1
    new.dir = paste0("./output_genotyping/",dir.nb)
  }
  if (length(grep(pattern = "output_genotyping",x = new.dir))==0){
    new.dir = paste0("./output_genotyping/",dir.nb)
  } else if (length(grep(pattern = "^./",x = new.dir,fixed = F))==0) {
    new.dir=paste0("./",new.dir)
  } else {new.dir=new.dir}
  if (!dir.exists(new.dir)){
    dir.create(new.dir,recursive = T)
    cat(paste0("Directory created : ",new.dir,"\n"))
  }
  
  path.data_clust = paste0(dirname,"/",filename)
  tmp=load(path.data_clust) # res_clust data_clustering
  # Permet d'eviter la note 'no visible bunding'
  names(tmp)=tmp
  res_clust=eval(parse(text = tmp["res_clust"]))
  data_clustering=eval(parse(text = tmp["data_clustering"]))
  new.filename=sub(pattern = "to_geno",replacement = "genotyped",x = filename)
  SampleName=unique(data_clustering$SampleName)
  if (same.pop==TRUE){ # Lance genotypage avec la fonction adequate
    res_genotyping=GenoAssign_pop_same(resClustering = res_clust,
                                       Dataset=data_clustering,
                                       SampleName = SampleName,
                                       SeuilNoCall = SeuilNoCall,SeuilNbSD = SeuilNbSD,SeuilSD = SeuilSD,
                                       NbClustMax = ploidy+1,
                                       cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker,
                                       tryEqHW=tryEqHW)
  } else if (same.pop!=TRUE){ # Lance genotypage avec la fonction adequate
    res_genotyping=GenoAssign_pop_dif(resClustering = res_clust,
                                      Dataset=data_clustering,
                                      SampleName = SampleName,
                                      SeuilNoCall = SeuilNoCall,SeuilNbSD = SeuilNbSD,SeuilSD = SeuilSD,
                                      NbClustMax = ploidy+1,
                                      cr_marker=cr_marker,fld_marker=fld_marker,hetso_marker=hetso_marker)
  }
  if (!dir.exists(paste0(new.dir,"/all_geno"))){ # cree un dossier temporaire pour les resultats de chaque coeur
    dir.create(paste0(new.dir,"/all_geno"),recursive = TRUE)
    cat(paste0("Directory created for genotype: ",new.dir,"/all_geno","\n"))
  }
  if (file.exists(paste0(new.dir,"/all_geno/",new.filename))){cat(paste0("Warning, the file ",paste0(new.dir,"/all_geno/",new.filename),"has been replaced!\n"))}
  save(path.data_clust,res_genotyping,file = paste0(new.dir,"/all_geno/",new.filename))
}
