#' Shiny App for genotyping
#'
#' Launch a shiny interface to use GenoTriplo. Really easy to use and user friendly, this will help you gain time !
#'
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#' @import processx
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @import shinyBS
#' @import shinyFiles
#' @import shinythemes
#' @import htmltools
#' @import bslib
#' @importFrom utils read.delim read.table
#' @importFrom DT datatable dataTableOutput renderDataTable
#' @importFrom utils head
#' @importFrom rlang .data
#' @return void : most results are automatically saved
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' launch_GenoShiny()
#' }
#' 

launch_GenoShiny = function(){
  addResourcePath("www",system.file("www",package = "GenoTriplo"))
  # Define UI for application that draws a histogram
  ui <- fluidPage(
    theme = shinytheme("cerulean"),
    # Application title
    fluidRow(column(12,titlePanel(title = div("Genotyping with R",style="font-family:Arial;font-weight:bold;margin-bottom:-1em",
                                              img(src="/www/sysaaf.png",height=70,width=70),
                                              img(src="/www/inrae.png",height=50,width=120),
                                              img(src="/www/feamp2.png",height=70,width=110))))),
    fluidRow(column(12,titlePanel(title = div("By J.Roche",style="color:#000000;font-size:10px;height:20px;")))),
    fluidRow(
      column(6,
             h2("Choose your step",style="font-size:20px;font-weight:bold"),
             navset_pill_list(widths = c(3,9),id = 'nav',selected = 1,
                              nav_panel("1. Create dataset",value = 1,
                                        h3("AxiomGT1.summary.txt file from AXAS",style="font-size:16px;font-weight:bold"),
                                        shinyFilesButton(
                                          id="axas_txt_file",
                                          label="Choose file",
                                          title="Choose file",
                                          multiple=FALSE
                                        ),
                                        fileInput(inputId = "marker_file",
                                                  label = div("File with the list of the marker to keep (optional)",
                                                              bsButton(inputId = "q3",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                  accept = ".txt"),
                                        fileInput(inputId = "indiv_file",
                                                  label = div("File with the list of the individuals to keep (optional)",
                                                              bsButton(inputId = "q2",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                  accept = ".txt"),
                                        textInput(inputId = "save_name1",
                                                  label = div("Name of the directory & file to save for clustering",
                                                              bsButton(inputId = "q5",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                        actionButton(inputId = "launch_create",label = "Create dataset")
                              ),
                              nav_panel("2. Clustering",value = 2,
                                        conditionalPanel(condition = "input.showProClust == input.hideProClust",
                                                         actionButton(inputId = "showProClust",label = "Add more control")),
                                        conditionalPanel(condition = "input.showProClust != input.hideProClust",
                                                         actionButton(inputId = "hideProClust",label = "Remove added control"),
                                                         sliderInput(inputId = "n_iter",
                                                                     label = div("Number of iterations of clustering for each marker",
                                                                                 bsButton(inputId = "q6",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = 1,max = 10,value = 5,step = 1),
                                                         sliderInput(inputId = "Dmin",
                                                                     label = div("Minimal distance between 2 clusters (adjusted contrast scale)",
                                                                                 bsButton(inputId = "q7",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = 0,max = 0.5,
                                                                     value = 0.28,step = 0.01)),
                                        h3("Directory in './output_create' with files for clustering",style="font-size:16px;font-weight:bold"),
                                        shinyDirButton(id='dir_file_clust',label = "Choose directory",title = "Choose directory"),
                                        selectInput(inputId = "Ploidy1",
                                                    label = div("Choose the number of ploidy",
                                                                bsButton(inputId = "q9",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                    choices = c(2,3,4),
                                                    multiple = FALSE,
                                                    selected = 3),
                                        sliderInput(inputId = "n_core1",
                                                    label = div("Number of cores for parallelization",
                                                                bsButton(inputId = "q10",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                    min = 1,max = parallel::detectCores() - 1,value = parallel::detectCores() - 2,step = 1),
                                        textInput(inputId = "save_name2",
                                                  label = div("Name of the directory & file to save for genotyping",
                                                              bsButton(inputId = "q12",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                        actionButton(inputId = "launch_clust",label = "Run clustering")
                              ),
                              nav_panel("3. Genotyping",value = 3,
                                        conditionalPanel(condition = "input.showProGeno == input.hideProGeno",
                                                         actionButton(inputId = "showProGeno",label = "Add more control")),
                                        conditionalPanel(condition = "input.showProGeno != input.hideProGeno",
                                                         actionButton(inputId = "hideProGeno",label = "Remove added control"),
                                                         selectInput(inputId = "tryEqHW",
                                                                     label = div("Delete non-possible cluster size",
                                                                                 bsButton(inputId = "q190",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     choices = c("Yes","No"),selected = "Yes",multiple = FALSE),
                                                         sliderInput(inputId = "SeuilNoCall",
                                                                     label = div("No Call threshold (individual)",
                                                                                 bsButton(inputId = "q13",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = 0,max = 1,value = 0.85,step = 0.01),
                                                         sliderInput(inputId = "SeuilNbSD",
                                                                     label = div("Distance from cluster s center",
                                                                                 bsButton(inputId = "q14",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = 0,max = 4,value = 2.8,step = 0.1),
                                                         sliderInput(inputId = "SeuilSD",
                                                                     label = div("SD threshold for a cluster",
                                                                                 bsButton(inputId = "q15",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = 0,max = 0.5,value = 0.28,step = 0.01),
                                                         sliderInput(inputId = "fld_marker",
                                                                     label = div("FLD threshold",
                                                                                 bsButton(inputId = "q27",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = 0,max = 7,value = 3.4,step = 0.05),
                                                         sliderInput(inputId = "hetso_marker",
                                                                     label = div("HetSO threshold",
                                                                                 bsButton(inputId = "q28",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = -1.5,max = 1,value = -0.3,step = 0.01),
                                                         sliderInput(inputId = "cr_marker",
                                                                     label = div("Call Rate threshold (marker)",
                                                                                 bsButton(inputId = "q26",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     min = 0,max = 1,value = 0.97,step = 0.01)
                                                         
                                        ),
                                        h3("Directory in './output_clustering' with files for genotyping",style="font-size:16px;font-weight:bold"),
                                        shinyDirButton(id='dir_file_geno',label = "Choose directory",title = "Choose directory"),
                                        fileInput(inputId = "corres_ATCG",
                                                  label = div("Correspondance between Allele and ATCG (optional)",
                                                              bsButton(inputId = "q17",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                  accept = c(".csv",".txt")),
                                        conditionalPanel(condition = "output.corres_load",
                                                         uiOutput(outputId = "probeset"),
                                                         uiOutput(outputId = "colA"),
                                                         uiOutput(outputId = "colB")),
                                        selectInput(inputId = "Ploidy2",
                                                    label = div("Choose the number of ploidy",
                                                                bsButton(inputId = "q18",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                    choices = c(2,3,4),
                                                    multiple = FALSE,
                                                    selected = 3),
                                        selectInput(inputId = "pop",
                                                    label = div("Are individuals from a same population ?",
                                                                bsButton(inputId = "q19",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                    choices = c("Yes","No"),selected = "Yes",multiple = FALSE),
                                        sliderInput(inputId = "n_core2",
                                                    label = div("Number of cores for parallelization",
                                                                bsButton(inputId = "q20",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                    min = 1,max = parallel::detectCores() - 1,value = parallel::detectCores() - 2,step = 1),
                                        textInput(inputId = "save_name3",
                                                  label = div("Name of the result file to save (no extension required) (click the '?' for help)",
                                                              bsButton(inputId = "q22",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                        actionButton(inputId = "launch_geno",label = "Run genotyping")
                              ),
                              nav_panel("4.Visualisation",value=4,
                                        fileInput(inputId = "file_visu",
                                                  label = div("File finishing by _markerCR.csv containing markers to visualize",
                                                              bsButton(inputId = "q23",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                  accept = ".csv"),
                                        selectizeInput(
                                          inputId = "marker_name",
                                          label = "Choose (or type) the marker to visualize",
                                          choices = NULL,
                                          multiple = FALSE,
                                          options = list(
                                            maxOptions = 50,
                                            placeholder = "Type marker name..."
                                          )
                                        ),
                                        conditionalPanel(condition = "input.marker_name !== undefined && input.marker_name !== ''",#"input.marker_name != null",
                                                         selectInput(inputId = "categories",
                                                                     label = div("Select the category of marker",
                                                                                 bsButton(inputId = "q240",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                     choices = "All",selected = "All",multiple = FALSE),
                                                         #c("All","PolyHighResolution","NoMinorHomozygote","MonoHighResolution","OffTargetVariant","CallRateBelowThreshold","Other")
                                                         fileInput(inputId = "select_indiv",
                                                                   label = div("Select file of individuals to highligth (optionnal)",
                                                                               bsButton(inputId = "q24",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                   accept = ".txt"),
                                                         uiOutput(outputId = "which_indiv"),
                                                         fileInput(inputId = "select_plate",
                                                                   label = div("Select file with individuals and their plate (optionnal)",
                                                                               bsButton(inputId = "q25",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                   accept = ".txt"),
                                                         uiOutput(outputId = "which_plate"),
                                                         textInput(inputId = "save_name4",
                                                                   label = div("Name of the plot to save for clustering",
                                                                               bsButton(inputId = "q29",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                                         actionButton(inputId = "SavePlot",label = "Save plot"))
                                        
                              )
             ),
      ),
      column(5,
             h5(div("Warnings : Make sure your app is launched on the wanted directory. If not, quit the app and use setwd() !",
                    style="color:red;font-weight:bold;margin-bottom:3em")),
             conditionalPanel(condition = "input.nav==1",
                              h5("You can find './output_create/log_create_time.log' for more info."),
                              textOutput(outputId = "end1")),
             conditionalPanel(condition = "input.nav==2",
                              h5("You can find './output_create/log_clust_time.log' for more info."),
                              textOutput(outputId = "end2")),
             conditionalPanel(condition = "input.nav==3",
                              h5("You can find './output_genotyping/log_geno_time.log' for more info."),
                              dataTableOutput(outputId = "head_corres"),
                              textOutput(outputId = "end3"),
                              plotOutput(outputId = "plot_CR"),
                              dataTableOutput(outputId = "Table_geno")),
             conditionalPanel(condition = "input.nav==4 & input.marker_name !== undefined",
                              h5(div("Warnings : The file can take some time to be effectivelly uploaded !",
                                     style="color:red;font-weight:bold;margin-bottom:3em")),
                              plotOutput(outputId = "plot_visu",
                                         # width = "100%",
                                         # height = 700,
                                         click = clickOpts("plot_click"),
                                         brush = brushOpts("plot_brush")),
                              h4("Informations on the selected marker"),
                              dataTableOutput(outputId = "Tab4"),
                              h4("Informations on the selected individual (click)"),
                              verbatimTextOutput("click_info"),
                              h4("Informations on the selected individuals (brush)"),
                              verbatimTextOutput("brush_info"),
                              uiOutput(outputId = "Select_Geno"),
                              actionButton(inputId = "Change_Geno",label = "Change Genotype"),
                              br(),
                              h4("/!\\ Before recreating files : takes time if large dataset !",style="color:#B80000"),
                              actionButton("recreate_file", "Recreate files of genotype from .Rdata"),
             )
      ),
      column(1)
    ),
    bsTooltip(id = "q2",title = "The .txt file must be an unique column as follow (no header) :<br/>Indiv1<br/>Indiv2<br/>Indiv3<br/>..."),
    bsTooltip(id = "q3",title = "The .txt file must be an unique column as follow (no header) :<br/>marker1<br/>marker2<br/>marker3<br/>..."),
    bsTooltip(id = "q4",title = "Type the number of expected markers<br/>Useless if you have a list of marker<br/>Type 0 for all"),
    bsTooltip(id = "q5",title = "Type the name you want for saving ;<br/>_to_clust.Rdata will automatically be added"),
    bsTooltip(id = "q6",title = "Number of iteration of clustering the algorithm have to perform for a single marker before choosing the best ;<br/>5 recommended"),
    bsTooltip(id = "q7",title = "The clustering algorithm asks for more than the espected number of cluster and merges those that do not respect a minimal distance set by this parameter ;<br/>0.28 recommended for diploid and triploid"),
    bsTooltip(id = "q9",title = "The number of ploidy (2 if diploid, 3 if triploid)"),
    bsTooltip(id = "q10",title = "The number of core the script can use<br/>Max is number of core of the computer -1 so you can t block your computer<br/>Max-1 recommended (default)"),
    bsTooltip(id = "q12",title = "Type the name you want for saving ;<br/>_to_geno.Rdata will automatically be added"),
    bsTooltip(id = "q13",title = "Threshold for the probability of an individual to belong to a cluster (else NoCall)<br/>0.85 recommended"),
    bsTooltip(id = "q14",title = "Distance (in number of standard deviation) to his cluster s center allowed for an individual (else NoCall)<br/>2.8 recommended"),
    bsTooltip(id = "q15",title = "Maximal standard deviation allowed for a cluster (else NoCall for the entire cluster) after removing individuals through the previous parameter<br/>0.15 recommended"),
    bsTooltip(id = "q17",title = "A .csv file with at least 3 columns named probeset_id, Allele_A and Allele_B containing respectively probeset names as in the dataset, ATCG correspondance to allele A of AXAS and ATCG correspondance to allele B of AXAS"),
    bsTooltip(id = "q18",title = "The number of ploidy (2 if diploid, 3 if triploid)"),
    bsTooltip(id = "q19",title = "Are all individuals in the dataset from a same population or not ?<br/>To know if it is expected to have markers with individuals only homozygous but some AA (from a population) and others BB (from another population)"),
    bsTooltip(id = "q190",title = "Do not change if you don't know what it means"),
    bsTooltip(id = "q20",title = "The number of core the script can use<br/>Max is number of core of the computer -1 so you can t block your computer<br/>Max-1 recommended (default)"),
    bsTooltip(id = "q22",title = "Base of saving name for multiple dataset<br/> Will be added :<br/>_genoAPIS.Rdata for the dataset ready for APIS assignation<br/>_indivCR.csv and _markerCR.csv respectively for Call Rate of individuals and marker<br/>_genoATCG.csv for real genotype of individuals (if a link file from AB to ATCG was given)<br/>_genotyped.Rdata for the R output with all data",trigger = "click"),
    bsTooltip(id = "q23",title = "The file ending by markerCR.csv created during the 3. Genotyping phase"),
    bsTooltip(id = "q24",title = "The .txt file must be an unique column as follow (no header) :<br/>Indiv1<br/>Indiv2<br/>Indiv3<br/>..."),
    bsTooltip(id = "q240",title = "Choose the category of marker you want to visualize. Select All if the category does not matter."),
    bsTooltip(id = "q25",title = "The .txt file (Sample_QCFilteredCR.txt from AXAS) must have its first two columns with header as follow (WITH header):<br/>Indiv_Name Plate_Number<br/>Indiv1 Plate1<br/>Indiv2 Plate2<br/>Indiv3 Plate3<br/>... ...<br/>It is not a problem if individual names have their cel attached (as Indiv1_K8.CEL)"),
    bsTooltip(id = "q26",title = "Call rate needed for a marker to be kept as good marker (below : bad ; above : good)"),
    bsTooltip(id = "q27",title = "FLD : Fisher s Linear Discriminant<br/>Measurement of the quality of a marker by taking into account the distance between the two closest genotypes and their standard deviation<br/>If under : bad ; if higher : good"),
    bsTooltip(id = "q28",title = "HetSO : Heterozygous Strength Offset<br/>Measurement of the position of the heterozygous cluster compared to homozygous clusters (y axis) as the heterozygous cluster is expected to be above homozygous clusters.<br/>If under : bad ; if higher : good"),
    bsTooltip(id = "q29",title = "Type the name you want for saving ;<br/>The marker name will automatically be added at the end of the name file;<br/>The plot is saved in the plot folder (which is automatically be created)")
  )
  
  # Define server logic
  server <- function(input, output,session) {
    
    # Info Fin Genotyping
    output$Table_geno = renderDataTable({
      if (length(fin_geno$cr_marker)>0){
        df=as.data.frame.matrix(table(fin_geno$cr_marker[,c("toKeep","nClus")],useNA = "ifany"))
        names(df)=paste0("nClus : ",names(df))
        rownames(df)[which(rownames(df)=="TRUE")]="Good marker"
        rownames(df)[which(rownames(df)=="FALSE")]="Bad marker"
        datatable(df, rownames = TRUE,options = list(dom = 't'))
      }
    })
    
    output$plot_CR = renderPlot({
      if (length(fin_geno$cr_indiv>0)){
        indiv = data.frame(CR=fin_geno$cr_indiv[,"CR"],Type="Indiv")
        marker = data.frame(CR=fin_geno$cr_marker[,"CR"],Type="Marker")
        # plot.scale = c("Indiv"="#3C70FC","Marker"="#DC3030")
        p1 = ggplot(data=indiv,aes(x=.data$CR))+
          geom_histogram(position="identity",binwidth = 0.005,col="#3C70FC",fill="#3C70FC")+
          theme_bw()+
          labs(x="Call Rate",y="Count",title = "Individuals")
        p2 = ggplot(data=marker,aes(x=.data$CR))+
          geom_histogram(position="identity",binwidth = 0.01,col="#DC3030",fill="#DC3030")+
          theme_bw()+
          labs(x="Call Rate",y="Count",title = "Marker")
        fin_geno$plot = plot_grid(p1,p2,ncol=2)
        fin_geno$plot
      }
    })
    
    ##### Info Visualisation #####
    data = reactiveValues(raw=data.frame(),
                          geno_all=data.frame(),
                          raw_all=data.frame(),
                          raw_marker=data.frame(),
                          geno=data.frame(),
                          valid=data.frame(),
                          marker0=data.frame(),
                          marker=c(),
                          plot=ggplot(),
                          data_plot=data.frame(),
                          plate=data.frame(),
                          indiv=c(),
                          SampleName=c())
    
    # To Save new genotype
    observeEvent(input$recreate_file,{
      if (length(data$valid)!=0 && !is.null(data$valid$locat.file)){
        dir.name.save.tmp=strsplit(x = data$valid$locat.file[1],split = "/")[[1]]
        dir.name.save=paste(dir.name.save.tmp[1:(length(dir.name.save.tmp)-2)],collapse = "/") # -2 car /all_geno est deja dans la fonction et on ne veut pas le nom d'un fichier
        new.dir=sub(pattern = "_markerCR.csv",replacement = "",x = input$file_visu$name)
        if (! is.null(data$valid$Allele_A) && ! is.null(data$valid$Allele_B) && ! is.null(data$valid$probeset_id)){
          create_genofiles(dir.name.save=dir.name.save,new.dir=new.dir,ploidy=max(data$valid$nClus,na.rm = T)-1,
                           corres_ATCG=data.frame(probeset_id=data$valid$probeset_id,Allele_A=data$valid$Allele_A,Allele_B=data$valid$Allele_B))
        } else {
          create_genofiles(dir.name.save=dir.name.save,new.dir=new.dir,ploidy=max(data$valid$nClus,na.rm = T)-1,corres_ATCG=NULL)
        }
      }
    })
    # round2=function(x){round(x,2)}
    observeEvent(input$file_visu,{
      if (!is.null(input$file_visu$datapath)){
        data$valid = read.table(file = input$file_visu$datapath,header = T,sep = ";")
        data$marker0 = data$valid %>%
          arrange(desc(.data$toKeep),desc(.data$MAF))
        # Categories remise a 0 pour eviter les problemes + changement des choix
        updateSelectInput(
          session,
          choices = c("All",unique(data$valid$Categorie)),
          inputId = "categories",
          selected = "All"
        )
      }
    })
    
    observeEvent(c(input$categories,data$marker0),{
      if (length(data$marker0)>0){
        if (input$categories=='All'){
          data$marker = data$marker0 %>% pull('MarkerName')
          # rownames(.data)
        } else {
          data$marker = data$marker0 %>%
            filter(.data$Categorie==input$categories) %>% pull('MarkerName')
          # rownames(.data)
        }
        updateSelectizeInput(
          session,
          "marker_name",
          choices = data$marker,
          server = TRUE
        )
      }
    })
    
    output$which_plate = renderUI({
      if (!is.null(input$select_plate$datapath)){
        indiv_plate = read.delim(file = input$select_plate$datapath,header=TRUE,sep = "\t",numerals = "no.loss")
        if (length(indiv_plate)==1){
          indiv_plate = read.delim(file = input$select_plate$datapath,header=TRUE,sep = ";",numerals = "no.loss")
        }
        names(indiv_plate)[1:2]=c("Indiv","Plate")
        indice1=regexpr(pattern="_[A-Z][0-9][0-9].CEL$",indiv_plate$Indiv)
        indice2=regexpr(pattern="_[A-Z][0-9].CEL$",indiv_plate$Indiv)
        indiv_plate$Indiv[indice1!=-1]=substr(x = indiv_plate$Indiv[indice1!=-1],start = 1,stop = (indice1[indice1!=-1]-1))
        indiv_plate$Indiv[indice2!=-1]=substr(x = indiv_plate$Indiv[indice2!=-1],start = 1,stop = (indice2[indice2!=-1]-1))
        data$plate = format(indiv_plate[,1:2],scientific=FALSE)
        selectInput(inputId = "which_plate_select",
                    label = div("Choose (or type) the plate to visualize",
                                tipify(el = bsButton(inputId = "2",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "You can choose the plate you want to visualize<br/>Select All if you want to see all individuals no matter the plate")),
                    choices = c("All",unique(format(indiv_plate$Plate,scientific=FALSE))),selected = "All")
      }
    })
    
    output$which_indiv = renderUI({
      if (!is.null(input$select_indiv$datapath)){
        data$indiv = read.delim(file = input$select_indiv$datapath,header=FALSE,sep = "\t")
        selectInput(inputId = "which_indiv_select",
                    label = div("Choose the individuals to visualize",
                                tipify(el = bsButton(inputId = "3",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "You can choose either to see the selected individuals written on the .txt file (Selected) or all individuals (All)")),
                    choices = c("All","Selected"),selected = "Selected")
      }
    })
    
    # Display possible genotype to change
    output$Select_Geno = renderUI({
      if (length(data$marker)>0){ # verify that a datasset has been loaded
        if (max(data$valid$nClus,na.rm=T)==4){
          radioButtons(inputId = "Select_Geno3",
                       label = div("Choose the genotype to change",
                                   tipify(el = bsButton(inputId = "30",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                          title = "Choose the genotype you want for the selected individuals on the graph to update them via the button below.")),
                       choices = c("3","2","1","0","-1"),selected = "3",inline = TRUE)
        } else if (max(data$valid$nClus,na.rm=T)==3){
          radioButtons(inputId = "Select_Geno2",
                       label = div("Choose the genotype to change",
                                   tipify(el = bsButton(inputId = "31",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                          title = "Choose the genotype you want for the selected individuals on the graph to update them via the button below.")),
                       choices = c("2","1","0","-1"),selected = "2",inline = TRUE)
        }
      }
    })
    
    # Update the genotype of the selected individuals
    observeEvent(input$Change_Geno,{
      if (max(data$valid$nClus,na.rm=T)==4){
        data$geno[,which(colnames(data$geno) %in% data$SampleName)] = as.numeric(input$Select_Geno3) #rownames(data$geno) == input$marker_name
        data$valid[data$valid$MarkerName==input$marker_name,c('toKeep','CR','FLD','HetSO','HomRO','nClus','MAF','Message')] = keepMarkertriplo(marker = input$marker_name,
                                                                                                                                               genotypePop = data$geno[which(rownames(data$geno) == input$marker_name),],
                                                                                                                                               data = data$raw_marker)
      } else if (max(data$valid$nClus,na.rm=T)==3){
        data$geno[,which(colnames(data$geno) %in% data$SampleName)] = as.numeric(input$Select_Geno2) #rownames(data$geno) == input$marker_name
        data$valid[data$valid$MarkerName==input$marker_name,c('toKeep','CR','FLD','HetSO','HomRO','nClus','MAF','Message')] = keepMarkerdiplo(marker = input$marker_name,
                                                                                                                                              genotypePop = data$geno[which(rownames(data$geno) == input$marker_name),],
                                                                                                                                              data = data$raw_marker)
      }
      data$geno_all[rownames(data$geno_all)==input$marker_name,] = data$geno
      load(data$valid$locat.file[data$valid$MarkerName==input$marker_name])
      res_genotyping[[1]] = data$geno_all
      res_genotyping[[3]][res_genotyping[[3]]$MarkerName == input$marker_name,]=data$valid[data$valid$MarkerName==input$marker_name,c('MarkerName','toKeep','CR','FLD','HetSO','HomRO','nClus','MAF','Message')]
      save(path.data_clust,res_genotyping,file = data$valid$locat.file[data$valid$MarkerName==input$marker_name])
    })
    
    ##### Graph #####
    output$plot_visu = renderPlot({
      if (!is.null(input$file_visu$datapath) & ! is.null(input$marker_name)){
        if (input$marker_name %in% data$marker){
          if (!input$marker_name %in% data$raw_marker$MarkerName){
            path.data_clust=NULL # avoid no visible binding
            data_clustering=NULL # avoid no visible binding
            load(data$valid$locat.file[data$valid$MarkerName==input$marker_name]) # res_genotyping & path.data_clust
            load(path.data_clust) # data_clustering
            data$geno_all=res_genotyping[[1]]
            data$raw_all=data_clustering
          }
          data$raw_marker = data$raw_all[data$raw_all$MarkerName==input$marker_name,]
          data$geno=data$geno_all[rownames(data$geno_all)==input$marker_name,]
          do_plot=TRUE
          alpha_manual = c("TRUE"=1,"FALSE"=0.2,"Selected"=1)
          col_manual = c("1"="black","2"="black")
          pch_manual = c("1"=21,"2"=23)
          xmin = min(data$raw_marker$Contrast,na.rm=T)
          xmax = max(data$raw_marker$Contrast,na.rm=T)
          ymin = min(data$raw_marker$SigStren,na.rm=T)
          ymax = max(data$raw_marker$SigStren,na.rm=T)
          if (! is.null(input$which_plate_select)){
            if (! input$which_plate_select=="All"){
              indiv_in_plate = data$plate$Indiv[data$plate$Plate==input$which_plate_select]
              ind = which(names(data$geno) %in% indiv_in_plate)
              if (length(ind)>0){
                data$data_plot=data.frame(Contrast=data$raw_marker$Contrast[data$raw_marker$SampleName %in% indiv_in_plate],
                                          SigStren=data$raw_marker$SigStren[data$raw_marker$SampleName %in% indiv_in_plate],
                                          fill=as.factor(as.numeric(data$geno[input$marker_name,ind])),
                                          alpha=factor(data$valid[input$marker_name,"toKeep"],levels=c("TRUE","FALSE","Selected")),
                                          col_pch=factor(x = 1,levels = c(1,2)),
                                          SampleName=data$raw_marker$SampleName[data$raw_marker$SampleName %in% indiv_in_plate])
              } else {
                do_plot=FALSE
              }
            } else {
              data$data_plot=data.frame(Contrast=data$raw_marker$Contrast,
                                        SigStren=data$raw_marker$SigStren,
                                        fill=as.factor(as.numeric(data$geno[input$marker_name,])),
                                        alpha=factor(data$valid[input$marker_name,"toKeep"],levels=c("TRUE","FALSE","Selected")),
                                        col_pch=factor(x = 1,levels = c(1,2)),
                                        SampleName=data$raw_marker$SampleName)
            }
          } else {
            data$data_plot=data.frame(Contrast=data$raw_marker$Contrast,
                                      SigStren=data$raw_marker$SigStren,
                                      fill=as.factor(as.numeric(data$geno[input$marker_name,])),
                                      alpha=factor(data$valid[input$marker_name,"toKeep"],levels=c("TRUE","FALSE","Selected")),
                                      col_pch=factor(x = 1,levels = c(1,2)),
                                      SampleName=data$raw_marker$SampleName)
          }
          if (!is.null(input$which_indiv_select)){
            if (input$which_indiv_select=='Selected'){
              selected_indiv = read.delim(file = input$select_indiv$datapath,header = FALSE,sep="\t")
              names(selected_indiv)="Indiv"
              data$data_plot$col_pch[data$data_plot$SampleName %in% selected_indiv$Indiv]=2
              data$data_plot$alpha[data$data_plot$SampleName %in% selected_indiv$Indiv]="Selected"
              alpha_manual = c("TRUE"=0.1,"FALSE"=0.1,"Selected"=1)
            }
          }
          if (max(data$valid$nClus,na.rm=T)==4){
            fill_manual=c("-1"="grey50","0"="red","1"="orange","2"="#15D409","3"="blue")
          } else if (max(data$valid$nClus,na.rm=T)==3) {
            fill_manual=c("-1"="grey50","0"="red","1"="#BCE30C","2"="blue")
          } else if (max(data$valid$nClus,na.rm=T)==5){
            fill_manual=c("-1"="grey50","0"="red","1"="orange","2"="purple","3"="#15D409","4"="blue")
          } else {
            message("Fonctionne avec diplo et triplo ; pas au dessus pour l'instant !")
          }
          if (do_plot){
            data$plot=ggplot(data=data$data_plot,aes(x=.data$Contrast,y=.data$SigStren,fill=.data$fill,shape=.data$col_pch,alpha=.data$alpha,col=.data$col_pch))+
              geom_point(pch=21,cex=3)+
              scale_fill_manual(values = fill_manual)+
              scale_alpha_manual(values = alpha_manual)+
              scale_color_manual(values = col_manual)+
              scale_shape_manual(values = pch_manual)+
              guides(alpha="none",col="none",shape="none")+
              theme_bw()+theme(title = element_text(size=22,face = "bold"),
                               axis.title = element_text(size=20),
                               legend.text = element_text(size=15),
                               legend.title = element_text(size=17),axis.text = element_text(size=15))+
              labs(x='Contrast AXAS',y='Signal Strength',fill='Genotype',title = input$marker_name)+
              xlim(c(xmin,xmax))+ylim(c(ymin,ymax))
            data$plot
          } else {
            data$plot=ggplot(data=data.frame(x=0,y=0),aes(x=.data$x,y=.data$y))+
              geom_text(x=0,y=0,label="No Individuals on that plate\nAre you sure it is the good .txt file ?")+
              xlim(-1,1)+ylim(-1,1)+
              theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    panel.background = element_blank(),
                    axis.text = element_blank())
            data$plot
          }
        }
      }
    })
    output$click_info <- renderPrint({
      if (!is.null(input$file_visu$datapath) & ! is.null(input$marker_name)){
        if (input$marker_name %in% data$marker){
          tmp = nearPoints(df = data$data_plot[,c('SampleName','Contrast','SigStren')],xvar = "Contrast",yvar = "SigStren",
                           coordinfo = input$plot_click, addDist = TRUE,threshold = 10,maxpoints = 1)
          print(tmp)
          data$SampleName = tmp$SampleName
        }
      }
    })
    
    output$brush_info <- renderPrint({
      if (!is.null(input$file_visu$datapath) & ! is.null(input$marker_name)){
        if (input$marker_name %in% data$marker){
          tmp = brushedPoints(df = data$data_plot[,c('SampleName','Contrast','SigStren')],xvar = "Contrast",yvar = "SigStren",
                              brush = input$plot_brush)
          print(tmp)
          data$SampleName = tmp$SampleName
        }
      }
    })
    
    observeEvent(input$SavePlot,{
      if (!dir.exists("./plot")){
        dir.create("./plot")
      }
      ggsave(filename = paste0("./plot/",input$save_name4,"_",input$marker_name,".png"),plot = data$plot,width = 9,height = 9)
    })
    output$Tab4 <- renderDataTable({
      if (!is.null(input$marker_name)){
        if (input$marker_name %in% data$marker){
          nb_indiv_clus = data.frame(X1=NA) # initialisation
          id_clus = unique(as.numeric(data$geno[input$marker_name,]))[which(unique(as.numeric(data$geno[input$marker_name,])) != -1)]
          for (k in sort(id_clus,decreasing = TRUE)){
            eval(parse(text = paste0("nb_indiv_clus$geno",k,"=length(which(data$geno[input$marker_name,]==k))")))
          }
          nb_indiv_clus=nb_indiv_clus[,-1]
          dta_display = data$valid %>% select(-c("FLD","HetSO","HomRO"))
          datatable(cbind(dta_display[dta_display$MarkerName==input$marker_name,],nb_indiv_clus,data.frame(Rank=which(data$marker0$MarkerName==input$marker_name))), rownames = TRUE,options = list(dom = 't'))
        }
      }
    })
    ##### End process #####
    end_process = reactiveValues(x="",
                                 y="",
                                 z="",
                                 datapath_create="",
                                 datapath_clust="")
    output$end1 = renderText({
      end_process$x
    })
    output$end2 = renderText({
      end_process$y
    })
    output$end3 = renderText({
      end_process$z
    })
    
    ##### Create dataset section #####
    roots = c(home=getwd())
    shinyFileChoose(input, "axas_txt_file", roots = roots, filetypes = c("", "txt", "csv"))
    
    observeEvent(input$axas_txt_file, {
      fileinfo = parseFilePaths(roots, input$axas_txt_file)
      end_process$x = paste0("File selected to create raw data : ",fileinfo$datapath)
      end_process$datapath_create=fileinfo$datapath
    })
    
    observeEvent(input$launch_create,{
      if (input$save_name1 == ""){
        end_process$x = paste0("No name -- Set to ",length(list.dirs("./output_create",recursive = F))," in ./output_create")
      }
      if (! end_process$datapath_create==""){
        if (!dir.exists("./output_create")){dir.create("./output_create")}
        
        # has_cmd = TRUE
        
        # os = Sys.info()["sysname"] # Operating system -> to check after if it is Windows or else
        t0=Sys.time()
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        # date_time=gsub(pattern = " CEST",replacement = "_",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 15)
        path_log = paste0("./output_create/log_create_",date_time,".log")
        write(x = "-----Create dataset-----",file = path_log)
        
        if (!is.null(input$marker_file)){
          marker = read.table(file = input$marker_file$datapath,header = FALSE)[,1]
        } else {marker=NULL}
        
        if (!is.null(input$indiv_file)){
          indiv = read.table(file = input$indiv_file$datapath,header = FALSE)[,1]
        } else {indiv=NULL}
        
        write(x = paste0("Start formating... Time : ",Sys.time()),file = path_log,append = TRUE)
        message(paste0("Start formating... Time : ",Sys.time()))
        Create_Dataset_from_file(filename=end_process$datapath_create,indiv_name=indiv,marker_name=marker,dir.name = input$save_name1)
        message(paste0("End process... Time : ",Sys.time()))
        write(x=paste0("End process... Time : ",Sys.time()),file = path_log,append=T)
        t_create = Sys.time()-t0
        end_process$x = paste0("Dataset saved ! Runing time : ",round(t_create,2)," ",units(t_create))
        gc()
      }
    })
    ##### Clustering section #####
    shinyDirChoose(input, "dir_file_clust", roots = roots) # roots set juste before create dataset part
    
    observeEvent(input$dir_file_clust, {
      end_process$datapath_clust = parseDirPath(roots, input$dir_file_clust)
      end_process$y = paste0("Folder selected for clustering : ",end_process$datapath_clust)
    })
    
    observeEvent(input$launch_clust,{
      if (end_process$datapath_clust!="" & input$save_name2 != ""){
        t0=Sys.time()
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 15)
        if (!dir.exists("./output_clustering")){dir.create("./output_clustering")}
        path_log = paste0("./output_clustering/log_clust_",date_time,".log")
        write(x = "-----Clustering-----",file = path_log)
        message(paste0("Start clustering... Time : ",Sys.time()))
        write(x = paste0("Start clustering... Time : ",Sys.time()),file = path_log,append = T)
        write(x = "--Parameters--",file = path_log,append = T)
        param = c(input$Ploidy1,input$save_name2,input$n_iter,input$Dmin,input$n_core1)
        param_names = c("Ploidy","Dir save name in './output_clustering'","N iter","D min","N core")
        for (k in 1:length(param)){
          write(x=paste0(param_names[k]," : ",param[k]),file = path_log,append=T)
        }
        Clustering_parallele_from_dir(dir_from_create=end_process$datapath_clust,
                                      new.dir=input$save_name2,
                                      ploidy=as.numeric(input$Ploidy1),
                                      n_iter=as.numeric(input$n_iter),
                                      Dmin=as.numeric(input$Dmin),
                                      n_core=as.numeric(input$n_core1))
        
        message(paste0("End clustering... Time : ",Sys.time()))
        write(x = paste0("End clustering... Time : ",Sys.time()),file = path_log,append = T)
        t_clust = Sys.time()-t0
        end_process$y = paste0("Clustering saved ! Runing time : ",round(t_clust,2)," ",units(t_clust))
      }
    })
    
    ##### Genotyping section #####
    fin_geno = reactiveValues(cr_indiv=data.frame(),cr_marker=data.frame(),plot=ggplot(),
                              corresLoad=FALSE,corres_ATCGtmp=data.frame(),corres_ATCG=data.frame(),atcg_choices=c(),
                              datapath_geno="")
    shinyDirChoose(input, "dir_file_geno", roots = roots) # roots set juste before create dataset part
    
    observeEvent(input$dir_file_geno, {
      fin_geno$datapath_geno = parseDirPath(roots, input$dir_file_geno)
      end_process$z = paste0("Folder selected for genotyping : ",fin_geno$datapath_geno)
    })
    
    observeEvent(input$corres_ATCG,{
      if (!is.null(input$corres_ATCG$datapath)){
        fin_geno$corres_ATCGtmp = read.table(file = input$corres_ATCG$datapath,header=T,sep=";",check.names = FALSE)
        if (length(fin_geno$corres_ATCGtmp[1,])==1){
          fin_geno$corres_ATCGtmp = read.table(file = input$corres_ATCG$datapath,header=T,sep="\t",check.names = FALSE)
        }
        fin_geno$atcg_choices = colnames(fin_geno$corres_ATCGtmp)
        if (length(fin_geno$atcg_choices)<3){
          stop("The correspondance file must have at least 3 columns ! If it has, try to change the separator by ';' or '\t'")
        }
        if (all(c("Allele_A","Allele_B","probeset_id") %in% fin_geno$atcg_choices)){
          fin_geno$corres_ATCG = fin_geno$corres_ATCGtmp
        } else {
          fin_geno$corresLoad = TRUE
        }
      }
    })
    
    output$head_corres=renderDataTable({
      if (fin_geno$corresLoad){
        datatable(head(fin_geno$corres_ATCGtmp),rownames = FALSE,options = list(dom = 't'),
                  caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Head correspondance table'))
      }
    })
    output$corres_load <- reactive({
      fin_geno$corresLoad
    })
    outputOptions(output, 'corres_load', suspendWhenHidden=FALSE)
    
    observeEvent(c(input$colA,input$colB,input$probeset),{
      if ((input$colA != input$colB) & (input$colA != input$probeset) & (input$colB != input$probeset)){
        fin_geno$corres_ATCG = fin_geno$corres_ATCGtmp %>% rename(Allele_A=input$colA,Allele_B=input$colB,probeset_id=input$probeset)
      }
    })
    
    output$colA = renderUI({
      if (fin_geno$corresLoad){
        selectInput(inputId = "colA",
                    label = div("Choose the variable corresponding to Allele_A",
                                tipify(el = bsButton(inputId = "4",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "You dont have at least one of these column names : Allele_A, Allele_B, probeset_id. Please select the corresponding column name.")),
                    choices = fin_geno$atcg_choices,multiple = FALSE)
      }
    })
    output$colB = renderUI({
      if (fin_geno$corresLoad){
        selectInput(inputId = "colB",
                    label = div("Choose the variable corresponding to Allele_B",
                                tipify(el = bsButton(inputId = "5",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "You dont have at least one of these column names : Allele_A, Allele_B, probeset_id. Please select the corresponding column name.")),
                    choices = fin_geno$atcg_choices,multiple = FALSE)
      }
    })
    output$probeset = renderUI({
      if (fin_geno$corresLoad){
        selectInput(inputId = "probeset",
                    label = div("Choose the variable corresponding to probeset_id",
                                tipify(el = bsButton(inputId = "6",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "You dont have at least one of these column names : Allele_A, Allele_B, probeset_id. Please select the corresponding column name.")),
                    choices = fin_geno$atcg_choices,multiple = FALSE)
      }
    })
    
    observeEvent(input$launch_geno,{
      fin_geno$cr_indiv = data.frame()
      fin_geno$cr_marker = data.frame()
      if (fin_geno$datapath_geno!="" & input$save_name3 != ""){
        if (!dir.exists("./output_genotyping")){dir.create("./output_genotyping")}
        dirpath = paste0("./output_genotyping/",input$save_name3,"/")
        if (!dir.exists(dirpath)){
          dir.create(dirpath,recursive = TRUE)
        }
        t0=Sys.time()
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 15)
        path_log = paste0("./output_genotyping/log_geno_",date_time,".log")
        write(x = "-----Genotyping-----",file = path_log)
        if (! is.null(input$corres_ATCG$datapath)){
          write(x = "Correspondance A-B to ATCG provided",file = path_log,append = T)
        } else {
          write(x = "Correspondance A-B to ATCG NOT provided",file = path_log,append = T)
          fin_geno$corres_ATCG = NULL
        }
        message(paste0("Start genotyping... Time : ",Sys.time()))
        write(x = paste0("Start genotyping... Time : ",Sys.time()),file = path_log,append = T)
        write(x = "--Parameters--",file = path_log,append = T)
        param = c(input$Ploidy2,input$save_name3,input$SeuilNoCall,input$SeuilNbSD,input$SeuilSD,
                  input$n_core2,input$pop,input$cr_marker,input$fld_marker,input$hetso_marker)
        param_names = c("Ploidy","Save name","Threshold No Call","Threshold Nb SD","Threshold SD","N core","Same pop","CR marker","FLD marker","HetSO marker")
        for (k in 1:length(param)){
          write(x=paste0(param_names[k]," : ",param[k]),file = path_log,append=T)
        }
        get.file=Genotyping_parallele_from_dir(dir_from_clustering = fin_geno$datapath_geno,new.dir = input$save_name3,
                                               ploidy = as.numeric(input$Ploidy2),
                                               SeuilNoCall = as.numeric(input$SeuilNoCall),
                                               SeuilNbSD = as.numeric(input$SeuilNbSD),
                                               SeuilSD = as.numeric(input$SeuilSD),
                                               n_core=as.numeric(input$n_core2),
                                               corres_ATCG = fin_geno$corres_ATCG,
                                               same.pop = ifelse(input$pop=='Yes',TRUE,FALSE),
                                               cr_marker = as.numeric(input$cr_marker),
                                               fld_marker = as.numeric(input$fld_marker),
                                               hetso_marker = as.numeric(input$hetso_marker),
                                               tryEqHW = ifelse(input$tryEqHW=='Yes',TRUE,FALSE))
        
        fin_geno$cr_indiv = read.table(file = get.file[[1]],header = T,sep = ";")
        fin_geno$cr_marker = read.table(file = get.file[[2]],header = T,sep = ";")
        t_geno=Sys.time()-t0
        fin_geno$corres_ATCGtmp = data.frame()
        end_process$z = paste0("Genotype saved ! Runing time : ",round(t_geno,2)," ",units(t_geno))
      }
    })
  }
  
  # Set Max size used in Shiny (heavy files for example)
  options(shiny.maxRequestSize=10000*1024^2) # 10000 pour 10Go => augmenter si besoin
  shinyApp(ui = ui, server = server)
}

# For image in packages :
# https://shiny.rstudio.com/reference/shiny/1.0.2/addresourcepath
