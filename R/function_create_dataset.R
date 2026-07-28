#' Create dataset from filename in appropriate format
#'
#' Create SigStren and Contrast variables from luminescence values of probeset A and B of each markers and return a dataframe to be used for clustering or save the result if a saving name is given
#'
#' @param filename filename from AXAS : `AxiomGT1.summary.txt`
#' @param dir.name saving name (the number of the chunk and _to_clust is automatically added)
#' @param indiv_name vector with name of individuals to keep
#' @param marker_name vector with name of marker to keep
#' @param n_marker_chunk nombre de lignes par chunk (nombre pair car chaque marker est par pair de ligne)
#' @param simplify.names Samplenames might be like so : IndivXX_a1_a1.CEL (if from Axiom) -> allows to select only IndivXX as SampleName
#' 
#' @import dplyr
#' @importFrom data.table fread
#' @return save file in `./output_create/dir.name` as `k_to_clust.Rdata` with k the chunk number.
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' filepath = system.file("extdata", "AxiomGT1.summary_examples.txt", package = "GenoTriplo")
#' Create_Dataset_from_file(filename=filepath)
#' }
Create_Dataset_from_file = function(filename="./AxiomGT1.summary.txt",dir.name=NULL,indiv_name=NULL,marker_name=NULL,n_marker_chunk=4000,simplify.names=TRUE){
  if (is.null(dir.name)){
    dir.name=length(list.dirs("./output_create",recursive = F))+1
  }
  if(length(grep(pattern = "output_create",x = dir.name))==0){
    dir.save = paste0("./output_create/",dir.name)
  } else if (length(grep(pattern = "^./",x = dir.name,fixed = F))==0){
    dir.save=paste0("./",dir.name)
  } else {dir.save=dir.name}
  
  if (!dir.exists(dir.save)){
    dir.create(dir.save,recursive = TRUE)
    cat(paste0("Directory created : ",dir.save,"\n"))
  }
  # Determination du nombre de ligne de description
  con = file(filename, open = "r")
  coln_skip = 0
  
  repeat {
    line = readLines(con, n = 1)
    if (length(line) == 0 || !startsWith(line, "#")) break
    coln_skip = coln_skip + 1
  }
  
  close(con)
  
  # Lecture du header des colonnes
  if (simplify.names){
    coln = as.vector(sapply(
      X   = colnames(fread(file = filename, skip = coln_skip, nrows = 1)),
      FUN = function(x) strsplit(x, "_")[[1]][1]
    ))
  } else {
    coln = as.vector(sapply(
      X   = colnames(fread(file = filename, skip = coln_skip, nrows = 1)),
      FUN = function(x) x
    ))
  }
  # Gestion des doublons -> on donne un nom different si ca arrive
  len.coln = length(coln)
  if (len.coln!=length(unique(coln))){
    l.dub = c()
    for (m in 1:(len.coln-1)){
      if (length(which(coln[m+1:len.coln]==coln[m]))>0){
        l.dub=c(l.dub,coln[m])
      }
    }
    tab = table(l.dub)
    for (id.name in names(tab)){
      for (id.name.num in 1:tab[id.name]){
        coln[match(id.name,coln)]=paste0(id.name,"_",id.name.num+1)
      }
    }
  }
  
  # Rownames du fichier
  first_col = fread(
    file   = filename,
    skip   = coln_skip,
    select = 1,
    header = FALSE
  )$V1
  
  # Verifie les lignes allant par paire
  has_suffix  = grepl("-[AB]$", first_col)
  base_names  = sub("-[AB]$", "", first_col)
  
  bases_A     = base_names[has_suffix & grepl("-A$", first_col)]
  bases_B     = base_names[has_suffix & grepl("-B$", first_col)]
  if (is.null(marker_name)){
    valid_bases = intersect(bases_A, bases_B)
  } else {
    valid_bases = intersect(intersect(bases_A, bases_B), marker_name)
  }
  if (length(valid_bases)==0){
    stop("valid_bases is empty. Maybe the list supplied is not the good one.")
  }
  useful_idx  = which(base_names %in% valid_bases & has_suffix)+coln_skip
  if (length(useful_idx)==0){
    stop("useful_idx is empty.")
  }
  header_skip = min(useful_idx)-1   # nombre de lignes à sauter avant les données
  
  # Ouvrir une connexion persistante
  con = file(filename, open = "r")
  
  # Sauter les lignes d'en-tete une seule fois
  readLines(con, n = header_skip)
  if (n_marker_chunk%%2!=0){stop("n_marker_chunk provided is not even !")}
  # Boucle de lecture par chunks + creation du data_clustering
  k      = 1
  tot_n=header_skip
  repeat {
    t0    = Sys.time()
    lines = readLines(con, n = n_marker_chunk)
    tmp_name = paste0(dir.name,"_",k)
    # Fin de fichier
    if (length(lines) == 0) break
    chunk = fread(text = paste(lines, collapse = "\n"), header = FALSE)
    chunk_rows_id = (tot_n+1):(tot_n+nrow(chunk))
    chunk = chunk[chunk_rows_id %in% useful_idx,]
    colnames(chunk)=coln
    if (!is.null(indiv_name)){
      indiv_name=c(coln[1],indiv_name)
      chunk = chunk %>% select(all_of(indiv_name))
    }
    if (nrow(chunk)>1 & nrow(chunk)%%2==0 & ncol(chunk)>1){
      data_clustering = Create_Dataset(data = chunk)
      if (file.exists(paste0(dir.save,"/",tmp_name,"_to_clust.Rdata"))){cat(paste0("Attention, le fichier ",paste0(dir.save,"/",tmp_name,"_to_clust.Rdata")," vient d'etre remplace!\n"))}
      save(data_clustering,file=paste0(dir.save,"/",tmp_name,"_to_clust.Rdata"))
    } else if (nrow(chunk)%%2!=0){
      stop("Nombre de ligne impair : des lignes indesirables doivent se trouver au milieu du fichier. Cela decale le cadre de lecture.\n
           Le probleme peut aussi venir de la liste d'individu ou de marker.")
    } else {
      # Nothing
    }
    t        = Sys.time()
    cat(sprintf("Chunk %d : %.2f sec -- ok\n", k, t-t0))
    k = k + 1
    
    # Si chunk incomplet = dernier chunk
    if (length(lines) < n_marker_chunk) break
    tot_n=tot_n+n_marker_chunk
  }
  close(con)
}

#' Create dataset in appropriate format
#'
#' Create SigStren and Contrast variables from luminescence values of probeset A and B of each markers and return a dataframe to be used for clustering or save the result if a saving name is given
#' Shouldn't be used outside `Create_Dataset_from_file`
#'
#' @param data dataframe with probeset_id or equivalent as first variable (markername finishing by -A or -B depending on the probeset) and individuals as variable with luminescence values for each probeset (dataset created by bash code by shiny app)
#' @import dplyr
#' @importFrom rlang .data
#'
#' @keywords internal
#' @noRd
#' @return number of individuals and markers (automatically save the dataset)
Create_Dataset = function(data){
  # if (is.null(data$probeset_id)){stop("Must have 'probeset_id' as first variable !")}
  # Séparer A et B
  idx_A = grepl("-A$", data[[1]])
  idx_B = grepl("-B$", data[[1]])
  
  data_A = data[idx_A, ]
  data_B = data[idx_B, ]
  
  # Nettoyer les noms
  probeset_id = sub("-A$", "", data_A[[1]])
  
  # Extraire matrices numeriques
  mat_A = as.matrix(data_A[, -1])
  mat_B = as.matrix(data_B[, -1])
  
  # Calculs vectorises
  SigStren = round((log2(mat_A) + log2(mat_B)) / 2, 3)
  Contrast = round(log2(mat_A / mat_B), 3)
  
  # Gerer Inf
  SigStren[is.infinite(SigStren)] = NA
  Contrast[is.infinite(Contrast)] = NA
  
  # Remise en format long propre (optionnel)
  data_clustering = data.frame(
    SampleName = rep(colnames(mat_A), each = nrow(mat_A)),
    MarkerName = rep(probeset_id, times = ncol(mat_A)),
    SigStren = as.vector(SigStren),
    Contrast = as.vector(Contrast)
  )
  
  return(data_clustering)
}
