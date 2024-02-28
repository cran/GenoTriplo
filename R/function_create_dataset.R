#' Create dataset in appropriate format
#'
#' Create SigStren and Contrast variables from luminescence values of probeset A and B of each markers and return a dataframe to be used for clustering or save the result if a saving name is given
#'
#' @param data dataframe with probeset_id as first variable (markername finishing by -A or -B depending on the probeset) and individuals as variable with luminescence values for each probeset (dataset created by bash code by shiny app)
#' @param save_name saving name
#' @import dplyr
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom utils memory.limit
#' @importFrom rlang .data
#'
#' @return number of individuals and markers (automatically save the dataset)
#'
#' @export

Create_Dataset = function(data,save_name=NULL){
  if (is.null(data$probeset_id)){stop("Must have 'probeset_id' as first variable !")}
  dta = tidyr::pivot_longer(data = data,cols = 2:dim(data)[2]) # on passe les SampleName en tant que variable unique (tableau moins large mais plus long : plus que 3 colonnes)
  rm(list='data')
  gc()
  dta$AorB = substr(dta$probeset_id,nchar(dta$probeset_id),nchar(dta$probeset_id)) # creation dune variable AorB suivant que le Marker soit A ou B
  dta$probeset_id = substr(dta$probeset_id,1,(nchar(dta$probeset_id)-2)) # retire le suffixe A/B suite a la creation de la variable AorB (-2 pour le '-A' ou '-B' des probesets)
  data_clustering = dta %>% group_by(.data$probeset_id,.data$name) %>%
    summarise(A=.data$value[.data$AorB=='A'],B=.data$value[.data$AorB=='B']) %>%
    mutate(SigStren=((log2(.data$A)+log2(.data$B))/2), # creation variable SigStren
           Contrast=log2(.data$A/.data$B)) %>%
    rename(SampleName=.data$name,MarkerName=.data$probeset_id) %>%
    select(.data$SampleName,.data$MarkerName,.data$SigStren,.data$Contrast)

  # dta$Id_unique=paste(dta$name,dta$probeset_id,sep = "_sep_") # creation dun identifiant pour chaque ligne en vue de la prochaine fonction
  # dta=dta[,c(3:5)]
  # gc()
  # if (getRversion()<'4.2.0'){
  #   memory.limit(size = 40000) # gestion de lespace memoire automatique a partir de la version 4.2.0
  # }
  # dta_wide= dta %>% tidyr::pivot_wider(  # inverse de pivot_longer : but est de creer une colonne A avec les valeur de A et B avec les valeurs de B
  #   id_cols=c(.data$Id_unique), # identifiant unique MesureIndiv et je laisse aussi les autres
  #   names_from=.data$AorB,  # nouvelles colonnes nommees a partir de cette variable
  #   values_from=.data$value) # valeur prise de cette variable
  # rm(list='dta')
  # gc()
  #
  # dta_final=as.data.frame(dta_wide)  # repasse en df plutot que tibble
  # rm(list='dta_wide')
  # gc()
  # pos_sep = as.numeric(regexpr(pattern = "_sep_",text = dta_final$Id_unique)) # on stock la position du separateur pour couper plus tard
  # data_clustering = dta_final %>%
  #   mutate(SigStren=((log2(.data$A)+log2(.data$B))/2), # creation variable SigStren
  #          Contrast=log2(.data$A/.data$B), # creation variable Contrast
  #          SampleName=substr(.data$Id_unique,1,pos_sep-1), # creation variable SampleName (on coupe Id_unique)
  #          MarkerName=substr(.data$Id_unique,pos_sep+5,nchar(.data$Id_unique)) # creation variable MarkerName (on coupe Id_unique)
  #          # ContrastCCS=asinh(4*((A-B)/(A+B)))/asinh(4) # creation de la variable ContrastCCS (utile lors du calcul decart type dun genotype)
  #   ) %>%
  #   select(-c(.data$Id_unique,.data$A,.data$B)) # retire les variable Id_unique, A et B
  if (!dir.exists("./output_create")){
    dir.create("./output_create")
  }
  if (!is.null(save_name)){
    save(data_clustering,file=paste0("./output_create/",save_name,"_to_clust.Rdata"))
    n_ind=length(unique(data_clustering$SampleName))
    n_marker=length(unique(data_clustering$MarkerName))
    return(c(n_ind,n_marker))
  } else {
    return(data_clustering)
  }
}

