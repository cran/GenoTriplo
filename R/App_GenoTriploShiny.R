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
#' @import shinythemes
#' @import htmltools
#' @importFrom utils read.delim read.table
#' @importFrom DT datatable dataTableOutput renderDataTable
#' @importFrom utils head
#' @importFrom rlang .data
#' @return void : most results are automatically saved
#'
#' @export

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
             navlistPanel(widths = c(3,9),id = 'nav',selected = 1,
                          tabPanel("1. Create dataset",value = 1,
                                   fileInput(inputId = "txt_file",
                                             label = div("AxiomGT1.summary.txt file from AXAS",
                                                         bsButton(inputId = "q1",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = ".txt"),
                                   fileInput(inputId = "marker_file",
                                             label = div("File with the list of the marker to keep (optional)",
                                                         bsButton(inputId = "q3",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = ".txt"),
                                   fileInput(inputId = "indiv_file",
                                             label = div("File with the list of the individuals to keep (optional)",
                                                         bsButton(inputId = "q2",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = ".txt"),
                                   textInput(inputId = "save_name1",
                                             label = div("Name of the file to save for clustering (no extension required)",
                                                         bsButton(inputId = "q5",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                   actionButton(inputId = "launch_create",label = "Create dataset")
                          ),
                          tabPanel("2. Clustering",value = 2,
                                   conditionalPanel(condition = "input.showProClust == input.hideProClust",
                                                    actionButton(inputId = "showProClust",label = "Add more control")),
                                   conditionalPanel(condition = "input.showProClust != input.hideProClust",
                                                    actionButton(inputId = "hideProClust",label = "Remove added control"),
                                                    sliderInput(inputId = "n_iter",
                                                                label = div("Number of iterations of clustering for each marker",
                                                                            bsButton(inputId = "q6",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                min = 1,max = 10,value = 5,step = 1),
                                                    sliderInput(inputId = "D_min",
                                                                label = div("Minimal distance between 2 clusters (adjusted contrast scale)",
                                                                            bsButton(inputId = "q7",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                min = 0,max = 0.5,
                                                                value = 0.28,step = 0.01)),
                                   fileInput(inputId = "data_file_clust",
                                             label = div("Data saved for clustering",
                                                         bsButton(inputId = "q8",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = c(".rda",".Rdata")),
                                   selectInput(inputId = "Ploidy1",
                                               label = div("Choose the number of ploidy",
                                                           bsButton(inputId = "q9",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                               choices = c(2,3),
                                               multiple = FALSE,
                                               selected = 3),
                                   sliderInput(inputId = "n_core1",
                                               label = div("Number of cores for parallelization",
                                                           bsButton(inputId = "q10",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                               min = 1,max = parallel::detectCores() - 1,value = parallel::detectCores() - 2,step = 1),
                                   textInput(inputId = "save_name2",
                                             label = div("Name of the file to save for genotyping (no extension required)",
                                                         bsButton(inputId = "q12",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                   actionButton(inputId = "launch_clust",label = "Run clustering")
                          ),
                          tabPanel("3. Genotyping",value = 3,
                                   conditionalPanel(condition = "input.showProGeno == input.hideProGeno",
                                                    actionButton(inputId = "showProGeno",label = "Add more control")),
                                   conditionalPanel(condition = "input.showProGeno != input.hideProGeno",
                                                    actionButton(inputId = "hideProGeno",label = "Remove added control"),
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
                                   fileInput(inputId = "clust_file_geno",
                                             label = div("Data saved for genotyping",
                                                         bsButton(inputId = "q16",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = ".Rdata"),
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
                                               choices = c(2,3),
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
                          tabPanel("4.Visualisation",value=4,
                                   fileInput(inputId = "file_visu",
                                             label = div("Data genotyped to visualize (can take some time, especially if it is a large dataset)",
                                                         bsButton(inputId = "q23",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = ".Rdata"),
                                   uiOutput(outputId = "marker_select"),
                                   conditionalPanel(condition = "input.marker_name !== undefined",
                                                    selectInput(inputId = "categories",
                                                                label = div("Select the category of marker",
                                                                            bsButton(inputId = "q240",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c("All","PolyHighResolution","NoMinorHomozygote","MonoHighResolution","OffTargetVariant","CallRateBelowThreshold","Other"),selected = "All",multiple = FALSE),
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
             h5(div("Warnings : Make sure yout app is launched on the wanted directory. If not, quit the app and use setwd() !",
                    style="color:red;font-weight:bold;margin-bottom:3em")),
             h5(div("Warnings : Wait for 'Upload complete' when uploading a file before running (won't work)",
                    style="color:red;font-weight:bold;margin-bottom:3em")),
             conditionalPanel(condition = "input.nav==1",
                              h5("You can find './output_create/log_create.log' for more info."),
                              textOutput(outputId = "end1")),
             conditionalPanel(condition = "input.nav==2",
                              h5("You can find './output_clustering/log_clust.log' for more info."),
                              textOutput(outputId = "end2")),
             conditionalPanel(condition = "input.nav==3",
                              h5("You can find './output_genotyping/saving_name/log_geno.log' for more info."),
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
                                         click = "plot_click",
                                         brush = "plot_brush"),
                              h4("Informations on the selected marker"),
                              dataTableOutput(outputId = "Tab4"),
                              h4("Informations on the selected individual (click)"),
                              verbatimTextOutput("click_info"),
                              h4("Informations on the selected individuals (brush)"),
                              verbatimTextOutput("brush_info"))
      ),
      column(1)
    ),
    bsTooltip(id = "q1",title = "Exit file of AXAS after running from .CEL datasets"),
    bsTooltip(id = "q2",title = "The .txt file must be an unique column as follow (no header) :<br/>Indiv1<br/>Indiv2<br/>Indiv3<br/>..."),
    bsTooltip(id = "q3",title = "The .txt file must be an unique column as follow (no header) :<br/>marker1<br/>marker2<br/>marker3<br/>..."),
    bsTooltip(id = "q4",title = "Type the number of expected markers<br/>Useless if you have a list of marker<br/>Type 0 for all"),
    bsTooltip(id = "q5",title = "Type the name you want for saving ;<br/>_to_clust.Rdata will automatically be added"),
    bsTooltip(id = "q6",title = "Number of iteration of clustering the algorithm have to perform for a single marker before choosing the best ;<br/>5 recommended"),
    bsTooltip(id = "q7",title = "The clustering algorithm asks for more than the espected number of cluster and merges those that do not respect a minimal distance set by this parameter ;<br/>0.22 recommended for diploid and triploid"),
    bsTooltip(id = "q8",title = "The file ending by _to_clust.Rdata created during the 1. Create dataset phase. Can also be a .Rdata/.rda file with one dataframe named data_clustering with four variables: SigStren, Contrast, SampleName and MarkerName."),
    bsTooltip(id = "q9",title = "The number of ploidy (2 if diploid, 3 if triploid)"),
    bsTooltip(id = "q10",title = "The number of core the script can use<br/>Max is number of core of the computer -1 so you can t block your computer<br/>Max-1 recommended (default)"),
    bsTooltip(id = "q12",title = "Type the name you want for saving ;<br/>_to_geno.Rdata will automatically be added"),
    bsTooltip(id = "q13",title = "Threshold for the probability of an individual to belong to a cluster (else NoCall)<br/>0.85 recommended"),
    bsTooltip(id = "q14",title = "Distance (in number of standard deviation) to his cluster s center allowed for an individual (else NoCall)<br/>2.8 recommended"),
    bsTooltip(id = "q15",title = "Maximal standard deviation allowed for a cluster (else NoCall for the entire cluster) after removing individuals through the previous parameter<br/>0.15 recommended"),
    bsTooltip(id = "q16",title = "The file ending by _to_geno.Rdata created during the 2. Clustering phase"),
    bsTooltip(id = "q17",title = "A .csv file with at least 3 columns named probeset_id, Allele_A and Allele_B containing respectively probeset names as in the dataset, ATCG correspondance to allele A of AXAS and ATCG correspondance to allele B of AXAS"),
    bsTooltip(id = "q18",title = "The number of ploidy (2 if diploid, 3 if triploid)"),
    bsTooltip(id = "q19",title = "Are all individuals in the dataset from a same population or not ?<br/>To know if it is expected to have markers with individuals only homozygous but some AA (from a population) and others BB (from another population)"),
    bsTooltip(id = "q20",title = "The number of core the script can use<br/>Max is number of core of the computer -1 so you can t block your computer<br/>Max-1 recommended (default)"),
    bsTooltip(id = "q22",title = "Base of saving name for multiple dataset<br/> Will be added :<br/>_genoAPIS.Rdata for the dataset ready for APIS assignation<br/>_indivCR.csv and _markerCR.csv respectively for Call Rate of individuals and marker<br/>_genoATCG.csv for real genotype of individuals (if a link file from AB to ATCG was given)<br/>_genotyped.Rdata for the R output with all data",trigger = "click"),
    bsTooltip(id = "q23",title = "The file ending by _genotyped.Rdata created during the 3. Genotyping phase"),
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
      if (length(fin_geno$bon_mauvais)>0){
        df=as.data.frame.matrix(table(fin_geno$bon_mauvais))
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
                          raw_marker=data.frame(),
                          geno=data.frame(),
                          valid=data.frame(),
                          marker0=data.frame(),
                          marker=c(),
                          plot=ggplot(),
                          plate=data.frame(),
                          indiv=c())

    # round2=function(x){round(x,2)}
    observeEvent(input$file_visu,{
      if (!is.null(input$file_visu$datapath)){
        load(file = input$file_visu$datapath)
        data$raw=data_clustering
        data$geno=res_geno[[1]]
        data$valid=res_geno[[3]] # mutate_if(res_geno[[3]],is.numeric,round2) # jai change dans la fonction de genotypage pour eviter de devoir le faire ici
        data$marker0 = data$valid %>%
          arrange(desc(.data$toKeep),desc(.data$MAF))
      }
    })

    observeEvent(c(input$categories,data$marker0),{
      if (length(data$marker0)>0){
        if (input$categories=='All'){
          data$marker = data$marker0 %>%
            rownames(.data)
        } else {
          data$marker = data$marker0 %>%
            filter(.data$Categorie==input$categories) %>%
            rownames(.data)
        }
      }
    })

    output$marker_select = renderUI({
      if (length(data$marker)>0){
        selectizeInput(inputId = "marker_name",
                       label = div("Choose (or type) the marker to visualize",
                                   tipify(el = bsButton(inputId = "1",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                          title = "You can choose the marker you want to visualize or type its name easily thank to autocompletion<br/>The marker are ordered by decreasing MAF with good marker first then bad marker")),
                       choices = data$marker,multiple = FALSE,options=list(maxOptions=500)) # length(data$marker)
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

    ##### Graph #####
    output$plot_visu = renderPlot({
      if (!is.null(input$file_visu$datapath) & ! is.null(input$marker_name)){
        if (input$marker_name %in% data$marker){
          data$raw_marker = data$raw[data$raw$MarkerName==input$marker_name,]
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
                df_plot=data.frame(Contrast=data$raw_marker$Contrast[data$raw_marker$SampleName %in% indiv_in_plate],
                                   SigStren=data$raw_marker$SigStren[data$raw_marker$SampleName %in% indiv_in_plate],
                                   fill=as.factor(as.numeric(data$geno[input$marker_name,ind])),
                                   alpha=factor(data$valid[input$marker_name,"toKeep"],levels=c("TRUE","FALSE","Selected")),
                                   col_pch=factor(x = 1,levels = c(1,2)),
                                   Name=data$raw_marker$SampleName[data$raw_marker$SampleName %in% indiv_in_plate])
              } else {
                do_plot=FALSE
              }
            } else {
              df_plot=data.frame(Contrast=data$raw_marker$Contrast,
                                 SigStren=data$raw_marker$SigStren,
                                 fill=as.factor(as.numeric(data$geno[input$marker_name,])),
                                 alpha=factor(data$valid[input$marker_name,"toKeep"],levels=c("TRUE","FALSE","Selected")),
                                 col_pch=factor(x = 1,levels = c(1,2)),
                                 Name=data$raw_marker$SampleName)
            }
          } else {
            df_plot=data.frame(Contrast=data$raw_marker$Contrast,
                               SigStren=data$raw_marker$SigStren,
                               fill=as.factor(as.numeric(data$geno[input$marker_name,])),
                               alpha=factor(data$valid[input$marker_name,"toKeep"],levels=c("TRUE","FALSE","Selected")),
                               col_pch=factor(x = 1,levels = c(1,2)),
                               Name=data$raw_marker$SampleName)
          }
          if (!is.null(input$which_indiv_select)){
            if (input$which_indiv_select=='Selected'){
              selected_indiv = read.delim(file = input$select_indiv$datapath,header = FALSE,sep="\t")
              names(selected_indiv)="Indiv"
              df_plot$col_pch[df_plot$Name %in% selected_indiv$Indiv]=2
              df_plot$alpha[df_plot$Name %in% selected_indiv$Indiv]="Selected"
              alpha_manual = c("TRUE"=0.1,"FALSE"=0.1,"Selected"=1)
            }
          }
          if (max(data$valid$nClus,na.rm=T)==4){
            fill_manual=c("-1"="grey50","0"="red","1"="orange","2"="#15D409","3"="blue")
          } else if (max(data$valid$nClus,na.rm=T)==3) {
            fill_manual=c("-1"="grey50","0"="red","1"="#BCE30C","2"="blue")
          } else {
            message("Fonctionne avec diplo et triplo ; pas au dessus pour l'instant !")
          }
          if (do_plot){
            data$plot=ggplot(data=df_plot,aes(x=.data$Contrast,y=.data$SigStren,fill=.data$fill,shape=.data$col_pch,alpha=.data$alpha,col=.data$col_pch))+
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
              labs(Contrast='Contrast AXAS',SigStren='Signal Strength',fill='Genotype',title = input$marker_name)+
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
          nearPoints(df = data$raw_marker[,c('SampleName','Contrast','SigStren')],coordinfo = input$plot_click, addDist = TRUE,threshold = 10,maxpoints = 10,)
        }
      }
    })

    output$brush_info <- renderPrint({
      if (!is.null(input$file_visu$datapath) & ! is.null(input$marker_name)){
        if (input$marker_name %in% data$marker){
          brushedPoints(df = data$raw_marker[,c('SampleName','Contrast','SigStren')],brush = input$plot_brush)
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
          datatable(cbind(dta_display[input$marker_name,],nb_indiv_clus,data.frame(Rank=which(data$marker==input$marker_name))), rownames = TRUE,options = list(dom = 't'))
        }
      }
    })
    ##### End process #####
    end_process = reactiveValues(x="",
                                 y="",
                                 z="")
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
    observeEvent(input$launch_create,{
      if (input$save_name1 == ""){
        end_process$x = "Saving name required !"
      }
      if (! is.null(input$txt_file) & input$save_name1 != ""){
        if (!dir.exists("./output_create")){dir.create("./output_create")}

        has_cmd = TRUE

        os = Sys.info()["sysname"] # Operating system -> to check after if it is Windows or else
        t0=Sys.time()
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        # date_time=gsub(pattern = " CEST",replacement = "_",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 16)
        path_log = paste0("./output_create/",date_time,"_log_create.log")
        write(x = "-----Create dataset-----",file = path_log)
        path=input$txt_file$datapath
        # path=gsub("\\","/",path,fixed=TRUE) # inutile
        write(x = paste0("Removing firsts lines from : ",input$txt_file$name," -- Time : ",Sys.time()),file = path_log,append = TRUE)
        sh_path = system.file('sh', package='GenoTriplo')
        if (os=="Windows"){
          run(command = "cmd.exe",c("/c","call","Extraction.sh",path),wd=sh_path)
        } else {
          run(paste0("sh ./Extraction.sh ",path),wd=sh_path)
        }

        if (!file.exists("./tmp.txt")){
          has_cmd = FALSE
          write(x = "Running from R as bash command couldn't run (see what to do to run bash command on Windows if wanted)",file = path_log,append = TRUE)
          raw_data = read.table(file = path,header = TRUE,comment.char = "#",check.names = FALSE)
          colnames(raw_data) = gsub(pattern = "_[A-Z][0-9].CEL",replacement = "",x = colnames(raw_data))
          colnames(raw_data) = gsub(pattern = "_[A-Z][0-9][0-9].CEL",replacement = "",x = colnames(raw_data))
        }

        if (!is.null(input$marker_file)){
          write(x = paste0("Selecting marker from : ",input$marker_file$name),file = path_log,append = TRUE)
          path=input$marker_file$datapath
          if (has_cmd){
            if (os=="Windows"){
              run(command = "cmd.exe",c("/c","call","Select_marker.sh",path),wd=sh_path)
            } else {
              run(paste0("sh ./Select_marker.sh ",path),wd=sh_path)
            }
          } else {
            marker = read.table(file = path,header = FALSE)
            colnames(marker)='MarkerName'
            MarkerName = gsub(pattern = "-A$",replacement = "",x = raw_data[,1])
            MarkerName = gsub(pattern = "-B$",replacement = "",x = MarkerName)
            to_keep = which(MarkerName %in% marker$MarkerName)
            raw_data = raw_data[to_keep,]
          }
        }

        if (!is.null(input$indiv_file)){
          write(x = paste0("Selecting individuals from : ",input$indiv_file$name),file = path_log,append = TRUE)
          path=input$indiv_file$datapath
          if (has_cmd){
            if (os=="Windows"){
              run(command = "cmd.exe",c("/c","call","Select_indiv.sh",path),wd=sh_path)
            } else {
              run(paste0("sh ./Select_indiv.sh ",path),wd=sh_path)
            }
          } else {
            indiv = read.table(file = path,header = FALSE)
            colnames(indiv)='SampleName'
            indiv$SampleName = gsub(pattern = "_[A-Z][0-9].CEL",replacement = "",x = indiv$SampleName)
            indiv$SampleName = gsub(pattern = "_[A-Z][0-9][0-9].CEL",replacement = "",x = indiv$SampleName)
            to_keep = which(colnames(raw_data) %in% indiv$SampleName)
            raw_data = raw_data[,c(1,to_keep)]
          }
        }

        if (has_cmd){
          write(x = paste0("Loading dataset... Time : ",Sys.time()),file = path_log,append = TRUE)
          message(paste0("Loading dataset... Time : ",Sys.time()))
          if (getRversion()<'4.2.0'){
            memory.limit(size = 40000) # gestion de lespace memoire automatique a partir de la version 4.2.0
          }
          dta0=read.delim(file = paste0(sh_path,"/tmp.txt"),header = TRUE,sep="\t")
          file.remove(paste0(sh_path,"/tmp.txt"))
        } else {
          dta0=raw_data
          rm(raw_data)
          gc()
        }

        write(x = paste0("Start formating... Time : ",Sys.time()),file = path_log,append = TRUE)
        message(paste0("Start formating... Time : ",Sys.time()))
        out=Create_Dataset(data = dta0,save_name = input$save_name1)
        write(x = paste0("Number of individuals : ",out[1]),file = path_log,append = TRUE)
        write(x = paste0("Number of markers : ",out[2]),file = path_log,append = TRUE)
        message(paste0("End process... Time : ",Sys.time()))
        write(x=paste0("End process... Time : ",Sys.time()),file = path_log,append=T)
        t_create = Sys.time()-t0
        end_process$x = paste0("Dataset saved ! Runing time : ",round(t_create,2)," ",units(t_create))
        rm(list='dta0')
        gc()
      }
    })
    ##### Clustering section #####
    observeEvent(input$launch_clust,{
      if (! is.null(input$data_file_clust) & input$save_name2 != ""){
        if (!dir.exists("./output_clustering")){dir.create("./output_clustering")}
        t0=Sys.time()
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        # date_time=gsub(pattern = " CEST",replacement = "_",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 16)
        path_log = paste0("./output_clustering/",date_time,"_log_clust.log")
        write(x = "-----Clustering-----",file = path_log)
        write(x = paste0("Load dataset : ",input$data_file_clust$name," -- ",Sys.time()),file = path_log,append = T)
        load(file = input$data_file_clust$datapath)
        l_marker = unique(data_clustering$MarkerName)
        nb_marker = length(l_marker)
        taille_batch_marker = 1000*as.numeric(input$n_core1)
        if (nb_marker>(1.5*taille_batch_marker)){
          write(x = paste0("The dataset is separated in smaller dataset to optimize running time !"),file = path_log,append = T)
          write(x = paste0("Cut dataset : ",input$data_file_clust$name," - Time : ",Sys.time()),file = path_log,append = T)
          message(paste0("Cut dataset... Time : ",Sys.time()))
          if (!dir.exists("./tmp")){
            dir.create("./tmp")
          }
          i_coupe = ifelse(nb_marker %% taille_batch_marker==0,nb_marker/taille_batch_marker,nb_marker %/% taille_batch_marker +1) # on va creer des sous jeu de donnes avec 10.000 marker
          for (k in 1:i_coupe){
            if (k!=i_coupe){
              dta=data_clustering[data_clustering$MarkerName %in% l_marker[(1+(k-1)*taille_batch_marker):(k*taille_batch_marker)],]
            } else {
              dta=data_clustering[data_clustering$MarkerName %in% l_marker[(1+(k-1)*taille_batch_marker):nb_marker],]
            }
            save(dta,file = paste0("./tmp/dta_",k,".Rdata"))
          }
          rm(list=c("data_clustering","dta"))
          gc()
          message(paste0("Start clustering... Time : ",Sys.time()))
          write(x = paste0("Start clustering... Time : ",Sys.time()),file = path_log,append = T)
          write(x = "--Parameters--",file = path_log,append = T)
          param = c(input$Ploidy1,input$save_name2,input$n_iter,input$D_min,input$n_core1)
          param_names = c("Ploidy","Save name","N iter","D min","N core")
          for (k in 1:length(param)){
            write(x=paste0(param_names[k]," : ",param[k]),file = path_log,append=T)
          }
          for (k in 1:i_coupe){
            load(paste0("./tmp/dta_",k,".Rdata"))
            Run_Clustering(data_clustering = dta,
                           ploidy = as.numeric(input$Ploidy1),
                           save_n = paste0(input$save_name2,"_",k),
                           n_iter=as.numeric(input$n_iter),
                           D_min=as.numeric(input$D_min),
                           n_core=as.numeric(input$n_core1),
                           path_log=path_log)
            rm(list="dta")
            gc()
          }
          unlink("./tmp",recursive = TRUE)
          message(paste0("Aggregating results... Time : ",Sys.time()))
          write(x = paste0("Aggregating results... Time : ",Sys.time()),file = path_log,append = T)
          load(paste0("./output_clustering/",input$save_name2,"_",1,"_to_geno.Rdata"))
          file.remove(paste0("./output_clustering/",input$save_name2,"_",1,"_to_geno.Rdata"))
          dta=data_clustering
          res=res_clust
          tim=delay
          for (k in 2:i_coupe){
            load(paste0("./output_clustering/",input$save_name2,"_",k,"_to_geno.Rdata"))
            file.remove(paste0("./output_clustering/",input$save_name2,"_",k,"_to_geno.Rdata"))
            gc()
            dta=rbind(dta,data_clustering)
            res=c(res,res_clust)
            tim=tim+delay
          }
          data_clustering=dta
          res_clust=res
          delay=tim
          save(res_clust,data_clustering,delay,file=paste0("./output_clustering/",input$save_name2,"_to_geno.Rdata"))
        } else {
          message(paste0("Start clustering... Time : ",Sys.time()))
          write(x = paste0("Start clustering... Time : ",Sys.time()),file = path_log,append = T)
          write(x = "--Parameters--",file = path_log,append = T)
          param = c(input$Ploidy1,input$save_name2,input$n_iter,input$D_min,input$n_core1)
          param_names = c("Ploidy","Save name","N iter","D min","N core")
          for (k in 1:length(param)){
            write(x=paste0(param_names[k]," : ",param[k]),file = path_log,append=T)
          }
          Run_Clustering(data_clustering = data_clustering,
                         ploidy = as.numeric(input$Ploidy1),
                         save_n = input$save_name2,
                         n_iter=as.numeric(input$n_iter),
                         D_min=as.numeric(input$D_min),
                         n_core=as.numeric(input$n_core1),
                         path_log=path_log)
        }
        message(paste0("End clustering... Time : ",Sys.time()))
        write(x = paste0("End clustering... Time : ",Sys.time()),file = path_log,append = T)
        t_clust = Sys.time()-t0
        end_process$y = paste0("Clustering saved ! Runing time : ",round(t_clust,2)," ",units(t_clust))
        # in_process$y = data.frame()
      }
    })
    ##### Genotyping section #####
    fin_geno = reactiveValues(cr_indiv=data.frame(),cr_marker=data.frame(),bon_mauvais=data.frame(),plot=ggplot(),
                              corresLoad=FALSE,corres_ATCGtmp=data.frame(),corres_ATCG=data.frame(),atcg_choices=c())

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
      fin_geno$bon_mauvais = data.frame()
      if (! is.null(input$clust_file_geno) & input$save_name3 != ""){
        if (!dir.exists("./output_genotyping")){dir.create("./output_genotyping")}
        if (!dir.exists(paste0("./output_genotyping/",input$save_name3))){
          dir.create(paste0("./output_genotyping/",input$save_name3),recursive = TRUE)
        }
        dirpath = paste0("./output_genotyping/",input$save_name3,"/")
        t0=Sys.time()
        message(paste0("Load dataset... Time : ",Sys.time()))
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        # date_time=gsub(pattern = " CEST",replacement = "_",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 16)
        path_log = paste0(dirpath,date_time,"_log_geno.log")
        write(x = "-----Genotyping-----",file = path_log)
        write(x = paste0("Load dataset : ",input$clust_file_geno$name," -- ",Sys.time()),file = path_log,append = T)
        load(file = input$clust_file_geno$datapath)
        if (! is.null(input$corres_ATCG$datapath)){
          write(x = "Correspondance A-B to ATCG provided",file = path_log,append = T)
        } else {
          write(x = "Correspondance A-B to ATCG NOT provided",file = path_log,append = T)
          fin_geno$corres_ATCG = NULL
        }
        l_marker = unique(data_clustering$MarkerName)
        nb_marker = length(l_marker)
        taille_batch_marker = 1000*as.numeric(input$n_core2)
        if (nb_marker>(1.5*taille_batch_marker)){
          message(paste0("Cut dataset... Time : ",Sys.time()))
          write(x = paste0("The dataset is separated in smaller dataset to optimize running time !"),file = path_log,append = T)
          write(x = paste0("Cut dataset - Time : ",Sys.time()),file = path_log,append = T)
          if (!dir.exists("./tmp")){
            dir.create("./tmp")
          }
          i_coupe = ifelse(nb_marker %% taille_batch_marker==0,nb_marker/taille_batch_marker,nb_marker %/% taille_batch_marker +1) # on va creer des sous jeu de donnes avec 10.000 marker
          for (k in 1:i_coupe){
            if (k!=i_coupe){
              dta_c=data_clustering[data_clustering$MarkerName %in% l_marker[(1+(k-1)*taille_batch_marker):(k*taille_batch_marker)],]
              res_c=res_clust[l_marker[(1+(k-1)*taille_batch_marker):(k*taille_batch_marker)]]
            } else {
              dta_c=data_clustering[data_clustering$MarkerName %in% l_marker[(1+(k-1)*taille_batch_marker):nb_marker],]
              res_c=res_clust[l_marker[(1+(k-1)*taille_batch_marker):nb_marker]]
            }
            save(dta_c,res_c,file = paste0("./tmp/dta_",k,".Rdata"))
          }
          save(data_clustering,file="./tmp/all.Rdata")
          rm(list=c("data_clustering","dta_c","res_clust","res_c"))
          gc()
          message(paste0("Start genotyping... Time : ",Sys.time()))
          write(x = paste0("Start genotyping... Time : ",Sys.time()),file = path_log,append = T)
          write(x = "--Parameters--",file = path_log,append = T)
          param = c(input$Ploidy2,input$save_name3,input$SeuilNoCall,input$SeuilNbSD,input$SeuilSD,
                    input$n_core2,input$pop,input$cr_marker,input$fld_marker,input$hetso_marker)
          param_names = c("Ploidy","Save name","Threshold No Call","Threshold Nb SD","Threshold SD","N core","Same pop","CR marker","FLD marker","HetSO marker")
          for (k in 1:length(param)){
            write(x=paste0(param_names[k]," : ",param[k]),file = path_log,append=T)
          }
          for (k in 1:i_coupe){
            load(paste0("./tmp/dta_",k,".Rdata"))
            file.remove(paste0("./tmp/dta_",k,".Rdata"))
            list_return_geno=Run_Genotyping(data_clustering = dta_c,
                                            res_clust = res_c,
                                            ploidy = as.numeric(input$Ploidy2),
                                            SeuilNoCall = as.numeric(input$SeuilNoCall),
                                            SeuilNbSD = as.numeric(input$SeuilNbSD),
                                            SeuilSD = as.numeric(input$SeuilSD),
                                            n_core=as.numeric(input$n_core2),
                                            corres_ATCG = fin_geno$corres_ATCG,
                                            pop = input$pop,
                                            cr_marker = as.numeric(input$cr_marker),
                                            fld_marker = as.numeric(input$fld_marker),
                                            hetso_marker = as.numeric(input$hetso_marker),
                                            save_n = input$save_name3,
                                            batch=k,
                                            ALL=FALSE,
                                            path_log=path_log)
            # fin_geno$cr_indiv = rbind(fin_geno$cr_indiv,list_return_geno[[1]])
            fin_geno$cr_marker = rbind(fin_geno$cr_marker,list_return_geno[[2]])
            fin_geno$bon_mauvais = rbind(fin_geno$bon_mauvais,list_return_geno[[3]])
            rm(list=c("dta_c","res_c"))
            gc()
          }
          message(paste0("Aggregating results... Time : ",Sys.time()))
          write(x = paste0("Aggregating results... Time : ",Sys.time()),file = path_log,append = T)
          load(paste0(dirpath,input$save_name3,"_",1,"_genoAPIS.Rdata"))
          file.remove(paste0(dirpath,input$save_name3,"_",1,"_genoAPIS.Rdata"))
          apis=data_APIS
          marker_cr = read.table(file = paste0(dirpath,input$save_name3,"_",1,"_markerCR.csv"),header = TRUE,sep = ";",check.names = FALSE)
          file.remove(paste0(dirpath,input$save_name3,"_",1,"_markerCR.csv"))
          if (!is.null(fin_geno$corres_ATCG)){
            atcg = read.table(file = paste0(dirpath,input$save_name3,"_",1,"_genoATCG.csv"),header = TRUE,sep = ";",check.names = FALSE)
            file.remove(paste0(dirpath,input$save_name3,"_",1,"_genoATCG.csv"))
          } else {
            ab = read.table(file = paste0(dirpath,input$save_name3,"_",1,"_genoAB.csv"),header = TRUE,sep = ";",check.names = FALSE)
            file.remove(paste0(dirpath,input$save_name3,"_",1,"_genoAB.csv"))
          }
          for (k in 2:i_coupe){
            load(paste0(dirpath,input$save_name3,"_",k,"_genoAPIS.Rdata"))
            file.remove(paste0(dirpath,input$save_name3,"_",k,"_genoAPIS.Rdata"))
            apis = cbind(apis,data_APIS)
            marker_cr2 = read.table(file = paste0(dirpath,input$save_name3,"_",k,"_markerCR.csv"),header = TRUE,sep = ";",check.names = FALSE)
            file.remove(paste0(dirpath,input$save_name3,"_",k,"_markerCR.csv"))
            marker_cr = rbind(marker_cr,marker_cr2)
            if (!is.null(fin_geno$corres_ATCG)){
              atcg2 = read.table(file = paste0(dirpath,input$save_name3,"_",k,"_genoATCG.csv"),header = TRUE,sep = ";",check.names = FALSE)
              file.remove(paste0(dirpath,input$save_name3,"_",k,"_genoATCG.csv"))
              atcg = cbind(atcg,atcg2[,-1]) #remove SampleName column
            } else {
              ab2 = read.table(file = paste0(dirpath,input$save_name3,"_",k,"_genoAB.csv"),header = TRUE,sep = ";",check.names = FALSE)
              file.remove(paste0(dirpath,input$save_name3,"_",k,"_genoAB.csv"))
              ab = cbind(ab,ab2[,-1]) #remove SampleName column
            }
            gc()
          }
          write.table(x = marker_cr,file = paste0(dirpath,input$save_name3,"_markerCR.csv"),sep = ";",row.names = FALSE)
          data_APIS = apis
          if (!is.null(fin_geno$corres_ATCG)){
            rownames(atcg)=atcg[,1]
            SN = atcg[,1]
            atcg=atcg[,-1]
            write.table(x = atcg,file = paste0(dirpath,input$save_name3,"_genoAB.csv"),sep = ";",row.names = TRUE,col.names = NA)
            # .ped and .map output
            markers = colnames(atcg)
            dataped = atcg %>%
              sapply(FUN = gsub,pattern='/',replacement=' ') %>%
              as.data.frame(.data) %>%
              sapply(FUN = gsub,pattern='NA',replacement='0') %>%
              as.data.frame(.data)
            rownames(dataped) = SN
            write(x = "#SampleName Genotypes",file = paste0(dirpath,input$save_name3,".ped"))
            write.table(x=dataped,file=paste0(dirpath,input$save_name3,".ped"),append=TRUE,col.names = FALSE,quote = FALSE,sep = "  ")
            write(x = "#MarkerID",file=paste0(dirpath,input$save_name3,".map"))
            write(x=markers,file = paste0(dirpath,input$save_name3,".map"),append=TRUE)
          } else {
            rownames(ab)=ab[,1]
            SN = ab[,1]
            ab=ab[,-1]
            write.table(x = ab,file = paste0(dirpath,input$save_name3,"_genoAB.csv"),sep = ";",row.names = TRUE,col.names = NA)
            # .ped and .map output
            markers = colnames(ab)
            dataped = ab %>%
              sapply(FUN = gsub,pattern='/',replacement=' ') %>%
              as.data.frame(.data) %>%
              sapply(FUN = gsub,pattern='NA',replacement='0') %>%
              as.data.frame(.data)
            rownames(dataped) = SN
            write(x = "#SampleName Genotypes",file = paste0(dirpath,input$save_name3,".ped"))
            write.table(x=dataped,file=paste0(dirpath,input$save_name3,".ped"),append=TRUE,col.names = FALSE,quote = FALSE,sep = "  ")
            write(x = "#MarkerID",file=paste0(dirpath,input$save_name3,".map"))
            write(x=markers,file = paste0(dirpath,input$save_name3,".map"),append=TRUE)
          }
          load(paste0(dirpath,input$save_name3,"_",1,"_genotyped.Rdata"))
          file.remove(paste0(dirpath,input$save_name3,"_",1,"_genotyped.Rdata"))
          res_g = res_geno
          tim = delay
          for (k in 2:i_coupe){
            load(paste0(dirpath,input$save_name3,"_",k,"_genotyped.Rdata"))
            file.remove(paste0(dirpath,input$save_name3,"_",k,"_genotyped.Rdata"))
            res_g[[1]] = rbind(res_g[[1]],res_geno[[1]])
            res_g[[2]] = c(res_g[[2]],res_geno[[2]])
            res_g[[3]] = rbind(res_g[[3]],res_geno[[3]])
            res_g[[4]] = cbind(res_g[[4]],res_geno[[4]])
            tim = tim + delay
            gc()
          }
          res_geno = res_g
          res_marker = res_geno[[3]] %>% filter(.data$toKeep)
          count_na = function(x){
            round(1-(length(which(x==-1))/length(x)),3)
          }
          fin_geno$cr_indiv = data.frame(SampleName=names(res_geno[[1]]),CR=apply(X = res_geno[[1]],MARGIN = 2,FUN = count_na))

          # Gestion doublon pour APIS
          db_ind = c()
          db_nam = c()
          for (pat in c("_2",".2","_BIS",".BIS")){ # add more si necessary
            indice = regexpr(pattern = pat,text = fin_geno$cr_indiv$SampleName,fixed = TRUE) # fixed : so that . is not a regular expression replacing all character
            db_ind = c(db_ind,fin_geno$cr_indiv$SampleName[which(indice!=-1)])
            db_nam = c(db_nam,fin_geno$cr_indiv$SampleName[which(fin_geno$cr_indiv$SampleName %in% substr(x = fin_geno$cr_indiv$SampleName[which(indice!=-1)],
                                                                                                          start = 1,
                                                                                                          stop = nchar(fin_geno$cr_indiv$SampleName[which(indice!=-1)])-nchar(pat,)))])
          }
          if (length(db_ind)>1){
            for (k in 1:length(db_ind)){
              cr1=fin_geno$cr_indiv$CR[fin_geno$cr_indiv$SampleName==db_ind[k]]
              cr2=fin_geno$cr_indiv$CR[fin_geno$cr_indiv$SampleName==db_nam[k]]
              if (cr1>cr2){
                data_APIS = data_APIS[-which(rownames(data_APIS)==db_nam[k]),] # suppress
                rownames(data_APIS)[which(rownames(data_APIS)==db_ind[k])] = db_nam[k] # rename because there is a _2 or else
              } else {
                data_APIS = data_APIS[-which(rownames(data_APIS)==db_ind[k]),] # just suppress, the good name is still here
              }
            }
          }

          save(data_APIS,res_marker,file = paste0(dirpath,input$save_name3,"_genoAPIS.Rdata"))
          load("./tmp/all.Rdata") # load data_clustering
          delay=tim
          save(res_geno,data_clustering,delay,file=paste0(dirpath,input$save_name3,"_genotyped.Rdata"))
          unlink("./tmp",recursive = TRUE)
        } else {
          message(paste0("Start genotyping... Time : ",Sys.time()))
          write(x = paste0("Start genotyping... Time : ",Sys.time()),file = path_log,append = T)
          write(x = "--Parameters--",file = path_log,append = T)
          param = c(input$Ploidy2,input$save_name3,input$SeuilNoCall,input$SeuilNbSD,input$SeuilSD,
                    input$n_core2,input$pop,input$cr_marker,input$fld_marker,input$hetso_marker)
          param_names = c("Ploidy","Save name","Seuil No Call","Seuil Nb SD","Seuil SD","N core","Same pop","CR marker","FLD marker","HetSO marker")
          for (k in 1:length(param)){
            write(x=paste0(param_names[k]," : ",param[k]),file = path_log,append=T)
          }
          list_return_geno=Run_Genotyping(data_clustering = data_clustering,
                                          res_clust = res_clust,
                                          ploidy = as.numeric(input$Ploidy2),
                                          SeuilNoCall = as.numeric(input$SeuilNoCall),
                                          SeuilNbSD = as.numeric(input$SeuilNbSD),
                                          SeuilSD = as.numeric(input$SeuilSD),
                                          n_core=as.numeric(input$n_core2),
                                          corres_ATCG = fin_geno$corres_ATCG,
                                          pop = input$pop,
                                          cr_marker = as.numeric(input$cr_marker),
                                          fld_marker = as.numeric(input$fld_marker),
                                          hetso_marker = as.numeric(input$hetso_marker),
                                          save_n = input$save_name3,
                                          batch="",
                                          ALL=TRUE,
                                          path_log=path_log)
          fin_geno$cr_indiv = list_return_geno[[1]]
          fin_geno$cr_marker = list_return_geno[[2]]
          fin_geno$bon_mauvais = list_return_geno[[3]]

          # .ped and .map output
          if (!is.null(fin_geno$corres_ATCG)){
            atcg_ab = read.table(file = paste0(dirpath,input$save_name3,"_genoATCG.csv"),sep=";",header=TRUE,check.names = FALSE)
          } else {
            atcg_ab = read.table(file = paste0(dirpath,input$save_name3,"_genoAB.csv"),sep=";",header=TRUE,check.names = FALSE)
          }
          SN = atcg_ab[,1]
          atcg_ab=atcg_ab[,-1]
          # .ped and .map output
          markers = colnames(atcg_ab)
          dataped = atcg_ab %>%
            sapply(FUN = gsub,pattern='/',replacement=' ') %>%
            as.data.frame(.data) %>%
            sapply(FUN = gsub,pattern='NA',replacement='0') %>%
            as.data.frame(.data)
          rownames(dataped) = SN
          write(x = "#SampleName Genotypes",file = paste0(dirpath,input$save_name3,".ped"))
          write.table(x=dataped,file=paste0(dirpath,input$save_name3,".ped"),append=TRUE,col.names = FALSE,quote = FALSE,sep = "  ")
          write(x = "#MarkerID",file=paste0(dirpath,input$save_name3,".map"))
          write(x=markers,file = paste0(dirpath,input$save_name3,".map"),append=TRUE)
        }
        write.table(x = fin_geno$cr_indiv,file = paste0(dirpath,input$save_name3,"_indivCR.csv"),sep = ";",row.names = FALSE)
        message(paste0("End genotyping... Time : ",Sys.time()))
        write(x = paste0("End genotyping... Time : ",Sys.time()),file = path_log,append = T)
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
