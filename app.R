

# Example Data Files -----------------------------------------------------------

Example_Matrix_File <- "Example_Data/Immgen_Microarray_Combined_SubsetNo3_v3_20231229.txt"
Example_Meta_File <- "Example_Data/Immgen_Microarray_Combined_SubsetNo5_meta.txt"

# Housekeeping gene libraries --------------------------------------------------

Eisenberg_hkg_File <- "default_housekeeping/eisenberg.txt"
Lin500_hkg_File  <- "default_housekeeping/lin500.txt"
HSIAO_hkg_File  <- "default_housekeeping/hsiao.txt"
Eisenberg_hkg_mm_File <- "default_housekeeping/eisenberg_mouse.txt"
Lin500_hkg_mm_File  <- "default_housekeeping/lin500_mouse.txt"
HSIAO_hkg_mm_File  <- "default_housekeeping/hsiao_mouse.txt"

# Gene set data ----------------------------------------------------------------
GeneSet_RData_File <- "GeneSet_Data/GeneSet_List_HS_v6.RData"
GeneSetTableIn_File <- "GeneSet_Data/GeneSet_CatTable_v6.txt"

# Gene Conversion Table --------------------------------------------------------
MM_HS_Conv_File <- "Mouse_To_Human_Gene_Conversion/hs_mm_homo_updated_v3_20231213.txt"

# Load Libraries ---------------------------------------------------------------
options(shiny.maxRequestSize = 10000*1024^2)

# Cran
library(shiny)
library(shinyjqui)
library(limma)
library(bapred)
library(DT)
library(shinycssloaders)
library(metafolio)
library(corrplot)
library(readr)
library(dplyr)
library(tidyverse)
library(zip)
library(glue)
library(factoextra)
library(FactoMineR)
library(magick)
library(ggfortify)
library(ggplotify)
library(draw)
library(clValid)
library(shinyWidgets)
library(tidyr)
library(umap)
library(plotly)
library(statVisual)
library(svglite)
# BioConductor
library(sva)
library(preprocessCore)#
library(scater)#
library(affyPLM)
library(Harman)
library(RUVSeq)
library(GSVA)#
library(ComplexHeatmap)#
library(pvca)#
library(edgeR)#
# Github
library(immunedeconv)


# Function to remove genes that are lowly expressed
ExprFilter2 <- function(vec,criteria,proportion) {
  Samp2meet <- length(vec)*proportion
  meet <- sum(vec>criteria)
  if (meet > Samp2meet) {
    return(TRUE)
  } else { return(FALSE) }
}

cv <- function(x){
  (sd(x)/mean(x))*100
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

MouseToHuman <- function(data,conv,gene_col = TRUE) {
  require(dplyr)
  if (!gene_col) {
    data <- cbind(Genes = rownames(data),data)
  }
  data <- merge(conv,data,by.x = "Mouse", by.y = colnames(data)[1], all.y = T)
  data[which(is.na(data[,"Human"])),"Human"] <- toupper(data[which(is.na(data[,"Human"])),"Mouse"])
  data <- data[,-1]
  colnames(data)[1] <- "Genes"
  if (TRUE %in% duplicated(data[, 1])) {
    data_dup <- data %>% dplyr::group_by(Genes) %>% dplyr::filter(n() > 1) %>% as.data.frame()
    data_nodup <- data %>% dplyr::group_by(Genes) %>% dplyr::filter(n() == 1) %>% as.data.frame()
    data_dup_summ <- data_dup %>%
      dplyr::group_by(Genes) %>%
      dplyr::summarise_all(max) %>%
      as.data.frame()
    data <- rbind(data_nodup,data_dup_summ)
  }
  rownames(data) <- data[,1]
  return(data)
}





# Read in back end data --------------------------------------------------------

## Housekeeping genes ----------------------------------------------------------
Eisenberg_hkg <- read.delim(Eisenberg_hkg_File, header = F, sep = '\t')
Lin500_hkg <- read.delim(Lin500_hkg_File, header = F, sep = '\t')
HSIAO_hkg <- read.delim(HSIAO_hkg_File, header = F, sep = '\t')
Eisenberg_hkg_mm <- read.delim(Eisenberg_hkg_mm_File, header = F, sep = '\t')
Lin500_hkg_mm <- read.delim(Lin500_hkg_mm_File, header = F, sep = '\t')
HSIAO_hkg_mm <- read.delim(HSIAO_hkg_mm_File, header = F, sep = '\t')

## Gene set data ---------------------------------------------------------------
geneset <- loadRData(GeneSet_RData_File)
geneset_df <- as.data.frame(read_delim(GeneSetTableIn_File, delim = '\t', col_names = T))
### Remove ER and TCGA genesets ------------------------------------------------
geneset_df <- geneset_df[which(!geneset_df[,1] %in% c("TCGA","ER Stress")),]

# Mouse to human conversion ----------------------------------------------------
MM_HS_Conv <- as.data.frame(read_delim(MM_HS_Conv_File,delim = '\t', col_names = F))
colnames(MM_HS_Conv) <- c("Human","Mouse")

# Just moved this from near the save observations in the server as temp_directory was needed earlier in the code
temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
dir.create(temp_directory)

umap.mod <- umap.defaults

# UI ---------------------------------------------------------------------------
ui <-
  shiny::navbarPage("{ Batch-Flex }",
                    ## Data Input Tab ------------------------------------------
                    shiny::tabPanel("Step 1 - Data Input",
                                    shiny::sidebarLayout(
                                      ### Sidebar ------------------------------
                                      shiny::sidebarPanel(
                                        h3("Matrix Parameters"),
                                        # File to be corrected
                                        shiny::fluidRow(
                                          shiny::column(8,
                                                        shiny::fileInput("uncorrected_matrix_input",
                                                                         "Matrix for Batch Correction"
                                                                         )
                                                        ),
                                          ## Deliminator for input matrix
                                          shiny::column(2, style = 'padding-left:2px;padding-right:2px',
                                                        shiny::selectInput("matrix_delim", "Deliminator",
                                                                           c("Tab" = '\t',"Comma" = ',', "Space" = ' ',"Semicolon" = ';', "Colon" = ':')
                                                                           )
                                                        ),
                                          shiny::column(2, style = "margin-top:15px",
                                                        shiny::radioButtons("HumanOrMouse",NULL,c("Human","Mouse"))
                                                        )
                                        ),
                                        # Need to log?
                                        shiny::fluidRow(
                                          shiny::column(4, style = "margin-top:-15px;padding-right:2px",
                                                        shiny::checkboxInput("RawCountCheck", "Raw counts input matrix"),
                                                        conditionalPanel(condition = "input.RawCountCheck",
                                                                         div(selectInput("RawCountNorm","Normalize Raw Counts By:",
                                                                                     c("No Normalization" = "none","TMM" = "TMM","Upper Quartile" = "upperquartile")),
                                                                             style = "margin-top:-10px")
                                                                         )
                                                        ),
                                          shiny::column(2, style = "margin-top:-15px",
                                                        shiny::checkboxInput("Log_Choice","Log2+1", value = T)
                                                        ),
                                          shiny::column(2, style = "margin-top:-25px",
                                                        numericInput("SeedSet","Set Seed:", value = 101, min = 1, step = 1)
                                                        ),
                                          shiny::column(4,
                                                        actionButton("UseExpData","Load Example Data")
                                                        )
                                          ),
                                        div(hr(),style = "margin-top:-20px;margin-bottom:-10px"),
                                        h3("Meta Data Parameters"),
                                        # Batch information in colnames?
                                        shiny::selectInput(
                                          "batch_info",
                                          "Is Batch Information Found in Column Names?",
                                          c("Yes", "No"),
                                          selected = "No"
                                        ),
                                        # Name derived or user provided meta file
                                        # If batch information is found in names
                                        shiny::uiOutput("user_batch_info_file"),
                                        shiny::uiOutput("batch_delim"),
                                        shiny::uiOutput("batch_names"),
                                        # if batch information is provided by user
                                        fluidRow(
                                          column(12, style = "margin-top:-20px",
                                                 shiny::uiOutput("rendInit_Batch_Select")
                                                 )
                                        )
                                      ),
                                      ### Main ---------------------------------
                                      shiny::mainPanel(
                                        shiny::uiOutput("rendMatrixHeader"),
                                        DT::dataTableOutput("uncorrected_matrix_output_input"),
                                        shiny::uiOutput("rendMetaHeader"),
                                        DT::dataTableOutput("meta_file")
                                      )
                                    )
                    ),
                    ## Batch Correction Tab ------------------------------------
                    shiny::tabPanel("Step 2 - Batch Correction",
                                    # Side Panel
                                    shiny::sidebarLayout(
                                      ### Sidebar ------------------------------
                                      shiny::sidebarPanel(
                                        width = 3,
                                        shiny::tabsetPanel(
                                          shiny::tabPanel("Batch Criteria",
                                                          p(),
                                                          #### Batch Crit Input----------------------------
                                                          shiny::fluidRow(
                                                            shiny::column(1,
                                                                          h2("1)")
                                                                          ),
                                                            shiny::column(8,
                                                                          uiOutput("rendbatch_correction_method"),
                                                                          ),
                                                            shiny::column(3, style = "margin-top:10px",
                                                                          uiOutput("rendQuantNorm")
                                                                          )
                                                            ),
                                                          ##### Limma ----------------------------
                                                          shiny::conditionalPanel(condition = "input.batch_correction_method == 'Limma'",
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("2)")
                                                                                                  ),
                                                                                    shiny::column(10,
                                                                                                  shiny::column(6,
                                                                                                                shiny::uiOutput("batch1_selection_limma")
                                                                                                                ),
                                                                                                  shiny::column(6,
                                                                                                                shiny::uiOutput("batch2_selection_limma")
                                                                                                                )
                                                                                                  )
                                                                                    ),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("3)")
                                                                                                  ),
                                                                                    shiny::column(11,
                                                                                                  shiny::uiOutput("covariate_choices_limma")
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          ##### ComBat ----------------------------
                                                          shiny::conditionalPanel(condition = "input.batch_correction_method == 'ComBat'",
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("2)")
                                                                                                  ),
                                                                                    shiny::column(7,
                                                                                                  shiny::uiOutput("batch_selection_ComBat")
                                                                                                  ),
                                                                                    shiny::column(4, style = "margin-top:20px",
                                                                                                  shiny::checkboxInput("combat_parametric", "Parametric?", value = TRUE)
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          ##### Mean Centering ----------------------------
                                                          shiny::conditionalPanel(condition = "input.batch_correction_method == 'Mean Centering'",
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("2)")
                                                                                                  ),
                                                                                    shiny::column(11,
                                                                                                  shiny::uiOutput("batch_selection_mean_centering")
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          ##### ComBatseq ----------------------------
                                                          shiny::conditionalPanel(condition = "input.batch_correction_method == 'ComBatseq'",
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("2)")
                                                                                                  ),
                                                                                    shiny::column(11,
                                                                                                  shiny::uiOutput("batch1_selection_ComBatseq")
                                                                                                  )
                                                                                    ),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("3)")
                                                                                                  ),
                                                                                    shiny::column(11,
                                                                                                  shiny::uiOutput("covariate_choices_ComBatseq")
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          ##### Harman -------------------------------
                                                          shiny::conditionalPanel(condition = "input.batch_correction_method == 'Harman'",
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("2)")
                                                                                                  ),
                                                                                    shiny::column(10,
                                                                                                  shiny::column(6,
                                                                                                                shiny::uiOutput("batch_selection_harman")
                                                                                                                ),
                                                                                                  shiny::column(6,
                                                                                                                shiny::uiOutput("treatment_selection_harman")
                                                                                                                )
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          ##### RUVg -------------------------------
                                                          conditionalPanel(condition = "input.batch_correction_method == 'RUVg'",
                                                                           shiny::fluidRow(
                                                                             shiny::column(1,
                                                                                           h2("2)")
                                                                                           ),
                                                                             shiny::column(11,
                                                                                           selectInput(
                                                                                             "RUVg_housekeeping_selection",
                                                                                             "Select Housekeeping Genes for Control",
                                                                                             c("Eisenberg", "Lin500", "HSIAO", "UserInput")),
                                                                                           conditionalPanel(condition = "input.RUVg_housekeeping_selection == 'UserInput'",
                                                                                                            fileInput("RUVg_user_control_genes", "Please Provide List of Control Genes"))
                                                                                           )
                                                                             ),
                                                                           fluidRow(
                                                                             column(4, style = 'padding-right:0px;',
                                                                                    numericInput("RUVg_estimate_factors","Factors to Estimate",
                                                                                                 value = 2,min = 1,max = 10,step = 1)
                                                                                    ),
                                                                             column(4,
                                                                                    numericInput("RUVg_drop_factors","Factors to Drop",
                                                                                                 value = 0,min = 0,max = 9,step = 1, width = "95%")
                                                                                    ),
                                                                             column(4, style = 'padding-left:0px',
                                                                                    selectInput("RUVg_tolerance","RUVg Tolerance",
                                                                                      c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10),selected = 1e-8)
                                                                                    )
                                                                             ),
                                                                           fluidRow(
                                                                             column(5,
                                                                                    checkboxInput("RUVg_mean_centered","Mean Centered",value = FALSE)
                                                                                    ),
                                                                             column(4,
                                                                                    checkboxInput("RUVg_rounded","Rounded",value = FALSE)
                                                                                    )
                                                                             )
                                                                           ),
                                                          ##### SVA --------------------------------
                                                          shiny::conditionalPanel(condition = "input.batch_correction_method == 'SVA'",
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(1,
                                                                                                  h2("2)")
                                                                                                  ),
                                                                                    shiny::column(10,
                                                                                                  shiny::column(6,
                                                                                                                shiny::selectInput("svaMethod_bc","Method",c("be","leek"))
                                                                                                                ),
                                                                                                  shiny::column(6,
                                                                                                                shiny::uiOutput("SVA_variable_of_interest_bc")
                                                                                                                )
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          div(hr(),style = "margin-top:-20px;margin-bottom:-10px"),
                                                          #### Plot Inputs ------------------------------------
                                                          ##### PCA Plot ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca'",
                                                                                  shiny::h3("PCA Parameters"),
                                                                                  shiny::radioButtons("PCA_type",NULL,c("Meta Annotation","Cluster Annotation"),inline = T),
                                                                                  shiny::conditionalPanel(condition = "input.PCA_type == 'Meta Annotation'",
                                                                                                          shiny::uiOutput("rendbatch_choices_PCA"),
                                                                                                          shiny::uiOutput("rendPCAhover")
                                                                                  ),
                                                                                  shiny::conditionalPanel(condition = "input.PCA_type == 'Cluster Annotation'",
                                                                                                          shiny::numericInput("cluster_number","Number of Clusters",
                                                                                                                              value = 2,step = 1,min = 1, width = "50%")
                                                                                                          )
                                                                                  ),
                                                          ##### PCA MC ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca_mc'",
                                                                                  shiny::h3("Multiple Component PCA Parameters"),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  shiny::uiOutput("PCA_mc_color_choice")
                                                                                                  ),
                                                                                    shiny::column(6,
                                                                                                  shiny::numericInput("PCA_mc_slider","Number of PCs",
                                                                                                                      value = 5,step = 1,min = 1)
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          ##### PCA Details ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca_dt'",
                                                                                  shiny::h3("PCA Detail Parameters"),
                                                                                  shiny::uiOutput("PCA_factors_choices")
                                                                                  ),
                                                          ##### UMAP ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'umap'",
                                                                                  shiny::h3("UMAP Parameters"),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(4,
                                                                                                  selectInput("UMAPmetricSelec","Metric:",
                                                                                                              choices = c("euclidean","manhattan","cosine","pearson","pearson2","hamming","correlation"))
                                                                                                  ),
                                                                                    shiny::column(4,
                                                                                                  numericInput("UMAPminDist","Min Distance:", value = 0.1,step = 0.05, min = 0)
                                                                                                  ),
                                                                                    shiny::column(4,
                                                                                                  numericInput("UMAPnnb","N Neighbors:", value = 15,step = 1, min = 5, max = 50)
                                                                                                  )
                                                                                    ),
                                                                                  selectInput("UMAPFeatureCategory","Select Feature Type:",
                                                                                              c("Meta Features","Matrix Features","PCA Projections","Immune Deconvolution Features","Gene Set Pathways")),
                                                                                  uiOutput("rendUMAPImmuneDeconvMethods"),
                                                                                  shiny::conditionalPanel(condition = "input.UMAPFeatureCategory != 'Gene Set Pathways'",
                                                                                                          fluidRow(
                                                                                                            column(9, style = 'padding-right:2px;',
                                                                                                                   selectizeInput("UMAPFeatSelection", label = "Select Feature:", choices = NULL,
                                                                                                                                  multiple = F, selected = 1,width = "100%")
                                                                                                                   ),
                                                                                                            column(3, style = 'padding-left:2px;',
                                                                                                                   selectInput("UMAPlogOpt","Log:", choices = c("No Log","Log2","Log2+1","Log10","Log10+1"))
                                                                                                                   )
                                                                                                            )
                                                                                                          ),
                                                                                  shiny::conditionalPanel(condition = "input.UMAPFeatureCategory == 'Gene Set Pathways'",
                                                                                                          selectInput("UMAPGeneSetCat","Gene Set Database:",choices = unique(geneset_df[,1])),
                                                                                                          div(DT::dataTableOutput("UMAPGeneSetTableUI"), style = "font-size:10px")
                                                                                                          )
                                                                                  ),
                                                          ##### Cluster ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'cluster_main'",
                                                                                  shiny::h3("Cluster Parameters"),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  numericInput("cluster_n_MV_features","Top Most Variable Features:",
                                                                                                               value = 1000)
                                                                                                  ),
                                                                                    shiny::column(6,
                                                                                                  selectInput("VarianceMeasure", "Select Variance Measure:",
                                                                                                              choices = c("MAD","CV","VAR"))
                                                                                                  )
                                                                                  )
                                                          ),
                                                          ##### Heatmap ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'cluster_main' & input.Cluster_main_pan == 'heatmap'",
                                                                                  div(shiny::h3("Heatmap Parameters"), style = "margin-top:-20px"),
                                                                                  uiOutput("rendHeatmapAnnoSel")
                                                          ),
                                                          ##### RLE ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'rle'",
                                                                                  shiny::h3("RLE Parameters"),
                                                                                  shiny::uiOutput("batch_choices_RLE")
                                                                                  ),
                                                          ##### Exp Var ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'exp_var'",
                                                                                  shiny::h3("Explanatory Variables Parameters"),
                                                                                  shiny::conditionalPanel(condition = "input.ExpVar_Plots == 'pvca'",
                                                                                                          shiny::fluidRow(
                                                                                                            shiny::column(6,
                                                                                                                          numericInput("pvcacluster_n_MV_features","Top Most Variable Features:",
                                                                                                                                       value = 20000)
                                                                                                                          ),
                                                                                                            shiny::column(6,
                                                                                                                          selectInput("pvcaVarianceMeasure", "Select Variance Measure:",
                                                                                                                                      choices = c("MAD","CV","VAR"))
                                                                                                                          )
                                                                                                            )
                                                                                                          ),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(9,
                                                                                                  shiny::uiOutput("variable_choices_EV")
                                                                                                  ),
                                                                                    shiny::column(3, style = "padding-left:2px",
                                                                                                  uiOutput("rendExpPlotLog"),
                                                                                                  uiOutput("rendpvcaPct")
                                                                                                  )
                                                                                    )
                                                                                  ),
                                                          ##### SVA ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'sva'",
                                                                                  shiny::h3("SVA Parameters"),
                                                                                  shiny::uiOutput("uncorrected_SVA_variable_of_interest"),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  shiny::selectInput("svaMethod","Method",c("leek","be"))
                                                                                                  ),
                                                                                    shiny::column(6,
                                                                                                  shiny::numericInput("SVAvarNum","Number of Variables",
                                                                                                                      value = 300,min = 0, step =1)
                                                                                                  )
                                                                                    ),
                                                                                  shiny::h4("Download Surrogate Variables"),
                                                                                  shiny::fluidRow(
                                                                                    column(6,
                                                                                           actionButton("save_SVA_surrogate_variables", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(6,
                                                                                           shiny::downloadButton("dnldsave_SVA_surrogate_variables","Dowload Meta with SVA")
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### Box Plot ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'box'",
                                                                                  style = "overflow-y:scroll;overflow-x:hidden;max-height:70vh;position:relative;",
                                                                                  shiny::h3("Box Plot Parameters"),
                                                                                  uiOutput("rendBPsampSubset"),
                                                                                  uiOutput("rendBPsampCriteria"),
                                                                                  uiOutput("rendBPgroupCriteria"),
                                                                                  uiOutput("rendBPgroupSelection"),
                                                                                  ## Limits the size of the select input box for group selection
                                                                                  tags$head(
                                                                                    tags$style(HTML(
                                                                                      '.selectize-input {
                                                                                      max-height: 102px;
                                                                                      overflow-y: auto;
                                                                                      }'
                                                                                      )
                                                                                      )
                                                                                    ),
                                                                                  selectInput("BPFeatureCategory","Select Feature Type:",
                                                                                              c("Matrix Features","Meta Features","PCA Projections","Immune Deconvolution Features","Gene Set Pathways")),
                                                                                  uiOutput("rendImmuneDeconvMethods"),
                                                                                  shiny::conditionalPanel(condition = "input.BPFeatureCategory != 'Gene Set Pathways'",
                                                                                                          fluidRow(
                                                                                                            column(9, style = 'padding-right:2px;',
                                                                                                                   selectizeInput("BPFeatSelection", label = "Select Feature:", choices = NULL,
                                                                                                                                  multiple = F, selected = 1,width = "100%")
                                                                                                                   ),
                                                                                                            column(3, style = 'padding-left:2px;',
                                                                                                                   selectInput("BPlogOpt","Log:", choices = c("No Log","Log2","Log2+1","Log10","Log10+1"))
                                                                                                                   )
                                                                                                            )
                                                                                                          ),
                                                                                  shiny::conditionalPanel(condition = "input.BPFeatureCategory == 'Gene Set Pathways'",
                                                                                                          selectInput("BPGeneSetCat","Gene Set Database:",choices = c(unique(geneset_df[,1]),"User Upload")),
                                                                                                          shiny::conditionalPanel(condition = "input.BPGeneSetCat == 'User Upload'",
                                                                                                                                  fileInput("UserGeneset","Upload Geneset",
                                                                                                                                            accept = c(".txt",".tsv",".csv",".zip",".RData",".gmt"))
                                                                                                                                  ),
                                                                                                          div(DT::dataTableOutput("GeneSetTableUI"), style = "font-size:10px")
                                                                                                          ),
                                                                                  fluidRow(
                                                                                    column(5, style = 'padding-right:2px;',
                                                                                           selectInput("BPplotstatComp","Stat Test Method:",
                                                                                                       choices = c("none","wilcox.test","t.test","kruskal.test","anova"))
                                                                                           ),
                                                                                    column(4, style = 'padding-left:4px;padding-right:4px;',
                                                                                           checkboxInput("BPplotsampledots","Include Dot Annotation", value = F)
                                                                                           ),
                                                                                    column(3, style = 'padding-left:4px;',
                                                                                           numericInput("BPplotDotSize","Dot Size:", value = 2, step = 0.25)
                                                                                           )
                                                                                    ),
                                                                                  fluidRow(
                                                                                    column(5, style = 'padding-right:2px;',
                                                                                           selectInput("BPplotXaxOrder","X-Axis Group Order",
                                                                                                       choices = c("Ascending","Descending","Not Specificed"))
                                                                                           ),
                                                                                    column(4, style = 'padding-left:4px;padding-right:4px;',
                                                                                           checkboxInput("BPremoveSingles","Remove groups 1 or less samples", value = T)
                                                                                           ),
                                                                                    column(3, style = 'padding-left:4px;',
                                                                                           checkboxInput("BPflipBP","Flip Axis", value = F),
                                                                                           radioButtons("BPorViolin",NULL,choices = c("Box Plot","Violin Plot"), selected = "Box Plot")
                                                                                           )
                                                                                    )
                                                          )
                                                          ),
                                          #### Figure Params ------------------------------------
                                          shiny::tabPanel("Figure Settings",
                                                          ##### PCA ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca'",
                                                                                  p(),
                                                                                  h4("PCA Figure Parameters"),
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           numericInput("pcaFontSize","Font Size:",
                                                                                                        value = 12, step = 1)
                                                                                    ),
                                                                                    column(6,
                                                                                           numericInput("PCAdotSize","Dot Size",value = 5)
                                                                                           )
                                                                                  ),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("pcaHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("pcaWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("pcaUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### PCA MC ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca_mc'",
                                                                                  p(),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("pca_mcHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("pca_mcWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("pca_mcUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### PCA Details ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca_dt'",
                                                                                  p(),
                                                                                  h3("PCA Details Figure Parameters"),
                                                                                  h4("Font Sizes"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("pcaDetAxisTkSize","Axis Tick:",
                                                                                                        value = 14, step = 1)
                                                                                           ),
                                                                                    column(4,
                                                                                           numericInput("pcaDetAxisTtSize","Axis Title:",
                                                                                                        value = 16, step = 1)
                                                                                           ),
                                                                                    column(4,
                                                                                           numericInput("pcaDetLabelSize","Bar Label:",
                                                                                                        value = 4.5, step = 0.5)
                                                                                           )
                                                                                    ),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("pca_dtHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("pca_dtWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("pca_dtUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### UMAP ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'umap'",
                                                                                  p(),
                                                                                  h3("UMAP Figure Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("umapAxisTkSize","Axis Tick Font:",
                                                                                                        value = 12, step = 1)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("umapAxisTtSize","Axis Title Font:",
                                                                                                        value = 14, step = 1)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("umapdotSize","Dot Size",value = 1)
                                                                                    )
                                                                                  ),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("umapHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("umapWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("umapUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### Cluster ------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'cluster_main' & input.Cluster_main_pan == 'cluster'",
                                                                                  p(),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("cluster_mainHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("cluster_mainWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("cluster_mainUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### Heatmap ---------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'cluster_main' & input.Cluster_main_pan == 'heatmap'",
                                                                                  p(),
                                                                                  h4("Heatmap Figure Parameters"),
                                                                                  checkboxGroupInput("HeatRowColNames",NULL,c("Turn on Row Names","Turn on Column Names"),inline = T),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           numericInput("heatmapHeight","Height (in)",value = 8)
                                                                                    ),
                                                                                    column(6,
                                                                                           numericInput("heatmapWidth","Width (in)",value = 10)
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### RLE ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'rle'",
                                                                                  p(),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("rleHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("rleWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("rleUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### Exp Var ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'exp_var'",
                                                                                  p(),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("exp_varHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("exp_varWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("exp_varUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### PVCA ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'pvca'",
                                                                                  p(),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("pvcaHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("pvcaWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("pvcaUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                          ),
                                                          ##### SVA ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'sva'",
                                                                                  p(),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("svaHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("svaWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("svaUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  ),
                                                          ##### Box Plot ------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'box'",
                                                                                  p(),
                                                                                  h4("Boxplot Figure Parameters"),
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           selectInput("BPTheme","Theme:",
                                                                                                       choices = c("Minimal" = "theme_minimal","Grey" = "theme_grey","BW" = "theme_bw",
                                                                                                                   "Linedraw" = "theme_linedraw","Light" = "theme_light","Dark" = "theme_dark",
                                                                                                                   "Classic" = "theme_classic","Void" = "theme_void","Test" = "theme_test"))
                                                                                    ),
                                                                                    column(6,
                                                                                           selectInput("BPxAxisOrient","X-Axis Label Angle",
                                                                                                       choices = c(90,45,0))
                                                                                    )
                                                                                  ),
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           textInput("BPplotHeight","Plot Height:",value = "500px")
                                                                                    ),
                                                                                    column(6,
                                                                                           textInput("BPplotWidth","Plot Width:",value = "100%")
                                                                                    )
                                                                                  ),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("BPplot1XAxisSize","X-Axis Font Size:",
                                                                                                        value = 16, step = 1)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("BPplot1YAxisSize","Y-Axis Font Size:",
                                                                                                        value = 16, step = 1)
                                                                                    ),
                                                                                    column(4,
                                                                                           textInput("BPplot1YAxisLim","Y-Axis Limit:",
                                                                                                     placeholder = "min,max")
                                                                                    )
                                                                                  ),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("BPHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("BPWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("BPUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                                                  )
                                                          ),

                                          #### File Export ------------------------------------
                                          shiny::tabPanel("Zip File Export",
                                                          shiny::p(),
                                                          shiny::textInput("file_name", "Please Name the Zipped File", value = "BatchFlex"),
                                                          shiny::downloadButton(
                                                            outputId = "download_btn",
                                                            label = "Download All Saved Files as Zip",
                                                            icon = icon("file-download")
                                                            ),
                                                          shiny::p(),
                                                          pickerInput(
                                                            "select_save_files",
                                                            "Choose Files to Zip",
                                                            list.files(temp_directory),
                                                            options = list(`actions-box` = TRUE),
                                                            multiple = TRUE
                                                            )
                                                          )
                                          )
                                        ),
                                      ### Main ---------------------------------
                                        mainPanel(
                                          shiny::tabsetPanel(
                                            id = "uncorrected_panel",
                                            #### Matrix ---------------------------------
                                            shiny::tabPanel(
                                              "Matrix",
                                              shiny::p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected Matrix"),
                                                              DT::dataTableOutput("uncorrected_matrix_output"),
                                                              p(),
                                                              shiny::fluidRow(
                                                                column(4, style = 'padding-right:0px;',
                                                                       shiny::actionButton("save_uncorrected_matrix", "Add to Zip File Export")
                                                                       ),
                                                                column(4, style = 'padding-left:0px;',
                                                                       shiny::downloadButton("dnldsave_uncorrected_matrix","Dowload Single Table")
                                                                )
                                                              )
                                                ),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected Matrix"),
                                                              DT::dataTableOutput("corrected_matrix"),
                                                              p(),
                                                              shiny::uiOutput("rendsave_corrected_matrix")
                                                )
                                              ),
                                              value = "mat"
                                            ),
                                            #### PCA -------------------------------------
                                            shiny::tabPanel(
                                              "PCA",
                                              p(),
                                              tabsetPanel(
                                                id = "PCA_main_pan",
                                                ##### PCA Plot ---------------------------------
                                                shiny::tabPanel(
                                                  "PCA Plot",
                                                  p(),
                                                  shiny::fluidRow(
                                                    shiny::column(6,
                                                                  shiny::h3("Uncorrected PCA"),
                                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotlyOutput("uncorrected_PCA")), type = 6),
                                                                  p(),
                                                                  shiny::fluidRow(
                                                                    column(4, style = 'padding-right:0px;',
                                                                           shiny::actionButton("save_uncorrected_PCA_plot","Add to Zip File Export")
                                                                           ),
                                                                    column(4, style = 'padding-left:0px;',
                                                                           shiny::downloadButton("dnldsave_uncorrected_PCA_plot","Dowload Single Plot")
                                                                           )
                                                                  )
                                                    ),
                                                    shiny::column(6,
                                                                  shiny::h3("Batch Corrected PCA"),
                                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotlyOutput("corrected_PCA")), type = 6),
                                                                  p(),
                                                                  shiny::fluidRow(
                                                                    column(4, style = 'padding-right:0px;',
                                                                           shiny::actionButton("save_corrected_PCA_plot", "Add to Zip File Export")
                                                                           ),
                                                                    column(4, style = 'padding-left:0px;',
                                                                           downloadButton("dnldsave_corrected_PCA_plot","Dowload Single Plot")
                                                                           )
                                                                    )
                                                                  )
                                                    ),
                                                  value = "pca"
                                                ),
                                                ##### PCA MC ---------------------------------
                                                shiny::tabPanel(
                                                  "Multiple Components PCA",
                                                  shiny::p(),
                                                  shiny::fluidRow(
                                                    shiny::column(6,
                                                                  shiny::h3("Uncorrected Multiple Components PCA"),
                                                                  shinycssloaders::withSpinner( shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_PCA_multiple_components")), type = 6),
                                                                  p(),
                                                                  shiny::fluidRow(
                                                                    column(4, style = 'padding-right:0px;',
                                                                           shiny::actionButton("save_uncorrected_PCA_mc_plot","Add to Zip File Export")
                                                                    ),
                                                                    column(4, style = 'padding-left:0px;',
                                                                           downloadButton("dnldsave_uncorrected_PCA_mc_plot","Dowload Single Plot")
                                                                    )
                                                                  )

                                                    ),
                                                    shiny::column(6,
                                                                  shiny::h3("Batch Corrected Multiple Components PCA"),
                                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_PCA_multiple_components")), type = 6),
                                                                  p(),
                                                                  shiny::fluidRow(
                                                                    column(4, style = 'padding-right:0px;',
                                                                           shiny::actionButton("save_corrected_PCA_mc_plot", "Add to Zip File Export")
                                                                    ),
                                                                    column(4, style = 'padding-left:0px;',
                                                                           downloadButton("dnldsave_corrected_PCA_mc_plot","Dowload Single Plot")
                                                                    )
                                                                  )

                                                    )
                                                  ),
                                                  value = "pca_mc"),
                                                ##### PCA Details ---------------------------------
                                                shiny:: tabPanel(
                                                  "PCA Details",
                                                  shiny::p(),
                                                  shiny::fluidRow(
                                                    shiny::column(6,
                                                                  shiny::h3("Uncorrected PCA Details"),
                                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_scree_plot")), type = 6),
                                                                  shiny::actionButton("save_uncorrected_scree_plot", "Save Uncorrected Scree Plot"),
                                                                  shiny::actionButton("save_uncorrected_PCA_components", "Save Uncorrected PCA Components"),
                                                                  shiny::p(),
                                                                  DT::dataTableOutput("uncorrected_contribution_table"),
                                                                  shiny::actionButton("save_uncorrected_contribution_table", "Save Uncorrected Contribution Table"),
                                                                  shiny::hr(),
                                                                  DT::dataTableOutput("uncorrected_contribution_counts"),
                                                                  shiny::actionButton("save_uncorrected_contribution_counts", "Save Uncorrected Contribution Counts")
                                                    ),
                                                    shiny::column(6,
                                                                  shiny::h3("Batch Corrected PCA Details"),
                                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_scree_plot")), type = 6),
                                                                  shiny::actionButton("save_corrected_scree_plot", "Save Corrected Scree Plot"),
                                                                  shiny::actionButton("save_corrected_PCA_components", "Save Corrected PCA Components"),
                                                                  shiny::p(),
                                                                  DT::dataTableOutput("corrected_contribution_table"),
                                                                  shiny::actionButton("save_corrected_contribution_table", "Save Corrected Contribution Table"),
                                                                  shiny::hr(),
                                                                  DT::dataTableOutput("corrected_contribution_counts"),
                                                                  shiny::actionButton("save_corrected_contribution_counts", "Save Corrected Contribution Counts")
                                                    )
                                                  ),
                                                  value = "pca_dt"
                                                ),
                                                ##### UMAP ---------------------------------
                                                shiny::tabPanel(
                                                  "UMAP",
                                                  shiny::p(),
                                                  shiny::fluidRow(
                                                    shiny::column(6,
                                                                  shiny::h3("Uncorrected UMAP"),
                                                                  shinycssloaders::withSpinner( shinyjqui::jqui_resizable(plotly::plotlyOutput("uncorrected_UMAP")), type = 6),
                                                                  p(),
                                                                  shiny::fluidRow(
                                                                    column(4, style = 'padding-right:0px;',
                                                                           shiny::actionButton("save_uncorrected_UMAP","Add to Zip File Export")
                                                                    ),
                                                                    column(4, style = 'padding-left:0px;',
                                                                           downloadButton("dnldsave_uncorrected_UMAP","Dowload Single Plot")
                                                                    )
                                                                  )

                                                    ),
                                                    shiny::column(6,
                                                                  shiny::h3("Batch Corrected UMAP"),
                                                                  shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotly::plotlyOutput("corrected_UMAP")), type = 6),
                                                                  p(),
                                                                  shiny::fluidRow(
                                                                    column(4, style = 'padding-right:0px;',
                                                                           shiny::actionButton("save_corrected_UMAP", "Add to Zip File Export")
                                                                    ),
                                                                    column(4, style = 'padding-left:0px;',
                                                                           downloadButton("dnldsave_corrected_UMAP","Dowload Single Plot")
                                                                    )
                                                                  )

                                                    )
                                                  ),
                                                  value = "umap"),
                                              ),
                                              value = "pca_main"
                                            ),
                                            #### Cluster ---------------------------------
                                            ##### Cluster Plots ---------------------------------
                                            shiny::tabPanel(
                                              "Cluster Plots",
                                              p(),
                                              tabsetPanel(
                                                id = "Cluster_main_pan",
                                                tabPanel("Cluster Plots",
                                                         p(),
                                                         shiny::fluidRow(
                                                           shiny::column(6,
                                                                         shiny::h3("Uncorrected Cluster Plots"),
                                                                         shiny::h4("Elbow Plot"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_elbow_plot", height = "250px")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_uncorrected_elbow_plot", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_uncorrected_elbow_plot","Dowload Single Plot")
                                                                           )
                                                                         ),
                                                                         shiny::hr(),
                                                                         shiny::h4("Silhouette Plot"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_silhouette_plot", height = "250px")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_uncorrected_silhouette_plot", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_uncorrected_silhouette_plot","Dowload Single Plot")
                                                                           )
                                                                         ),
                                                                         shiny::hr(),
                                                                         shiny::h4("Dunn Index Plot"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_dunn_index_plot", height = "250px")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_uncorrected_dunn_index_plot", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_uncorrected_dunn_index_plot","Dowload Single Plot")
                                                                           )
                                                                         )

                                                                         ),
                                                           shiny::column(6,
                                                                         shiny::h3("Batch Corrected Cluster Plots"),
                                                                         shiny::h4("Elbow Plot"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_elbow_plot", height = "250px")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_corrected_elbow_plot", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_corrected_elbow_plot","Dowload Single Plot")
                                                                           )
                                                                         ),
                                                                         shiny::hr(),
                                                                         shiny::h4("Silhouette Plot"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_silhouette_plot", height = "250px")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_corrected_silhouette_plot", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_corrected_silhouette_plot","Dowload Single Plot")
                                                                           )
                                                                         ),
                                                                         shiny::hr(),
                                                                         shiny::h4("Dunn Index Plot"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_dunn_index_plot", height = "250px")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_corrected_dunn_index_plot", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_corrected_dunn_index_plot","Dowload Single Plot")
                                                                           )
                                                                         )

                                                                         )
                                                           ),
                                                         value = "cluster"),
                                                ##### Heatmaps ---------------------------------
                                                tabPanel("Heatmaps",
                                                         p(),
                                                         shiny::fluidRow(
                                                           shiny::column(6,
                                                                         shiny::h3("Uncorrected Heatmap"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_heatmap", height = "700px", width = "100%")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_uncorrected_heatmap", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_uncorrected_heatmap","Dowload Single Plot")
                                                                           )
                                                                         )

                                                                         ),
                                                           shiny::column(6,
                                                                         shiny::h3("Batch Corrected Heatmap"),
                                                                         shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_heatmap", height = "700px", width = "100%")), type = 6),
                                                                         p(),
                                                                         shiny::fluidRow(
                                                                           column(4, style = 'padding-right:0px;',
                                                                                  shiny::actionButton("save_corrected_heatmap", "Add to Zip File Export")
                                                                           ),
                                                                           column(4, style = 'padding-left:0px;',
                                                                                  downloadButton("dnldsave_corrected_heatmap","Dowload Single Plot")
                                                                           )
                                                                         )

                                                           )
                                                         ),
                                                         value = "heatmap"
                                                         )
                                              ),
                                              value = "cluster_main"
                                            ),
                                            #### RLE ---------------------------------
                                            shiny::tabPanel(
                                              "Relative Log Expression Plot",
                                              p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected RLE Plot"),
                                                              shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_RLE_plot")), type = 6),
                                                              p(),
                                                              shiny::fluidRow(
                                                                column(4, style = 'padding-right:0px;',
                                                                       shiny::actionButton("save_uncorrected_RLE_plot", "Save Uncorrected RLE Plot")
                                                                ),
                                                                column(4, style = 'padding-left:0px;',
                                                                       downloadButton("dnldsave_uncorrected_RLE_plot","Dowload Single Plot")
                                                                )
                                                              )
                                                              ),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected RLE Plot"),
                                                              shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_RLE_plot")), type = 6),
                                                              p(),
                                                              shiny::fluidRow(
                                                                column(4, style = 'padding-right:0px;',
                                                                       shiny::actionButton("save_corrected_RLE_plot", "Save Corrected RLE Plot")
                                                                ),
                                                                column(4, style = 'padding-left:0px;',
                                                                       downloadButton("dnldsave_corrected_RLE_plot","Dowload Single Plot")
                                                                )
                                                              )

                                                )
                                              ),
                                              value = "rle"
                                            ),
                                            #### Exp Var ---------------------------------
                                            shiny::tabPanel(
                                              "Explanatory Variables",
                                              p(),
                                              shiny::tabsetPanel(
                                                id = "ExpVar_Plots",
                                                ##### Exp Var ---------------------------------
                                                shiny::tabPanel("Explanatory Variables Plot",
                                                                shiny::fluidRow(
                                                                  shiny::column(6,
                                                                                shiny::h3("Uncorrected Explanatory Variables Plot"),
                                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_EV_plot")), type = 6),
                                                                                p(),
                                                                                shiny::fluidRow(
                                                                                  column(4, style = 'padding-right:0px;',
                                                                                         shiny::actionButton("save_uncorrected_EV_plot", "Add to Zip File Export")
                                                                                  ),
                                                                                  column(4, style = 'padding-left:0px;',
                                                                                         downloadButton("dnldsave_uncorrected_EV_plot","Dowload Single Plot")
                                                                                  )
                                                                                )),
                                                                  shiny::column(6,
                                                                                shiny::h3("Batch Corrected Explanatory Variables Plot"),
                                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_EV_plot")), type = 6),
                                                                                p(),
                                                                                shiny::fluidRow(
                                                                                  column(4, style = 'padding-right:0px;',
                                                                                         shiny::actionButton("save_corrected_EV_plot", "Add to Zip File Export")
                                                                                  ),
                                                                                  column(4, style = 'padding-left:0px;',
                                                                                         downloadButton("dnldsave_corrected_EV_plot","Dowload Single Plot")
                                                                                  )
                                                                                )

                                                                  )
                                                                ),
                                                                value = "exp_plot"
                                                                ),
                                                ##### PVCA ---------------------------------
                                                shiny::tabPanel("PVCA",
                                                                shiny::fluidRow(
                                                                  shiny::column(6,
                                                                                shiny::h3("Uncorrected PVCA Plot"),
                                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_PVCA_plot")), type = 6),
                                                                                p(),
                                                                                shiny::fluidRow(
                                                                                  column(4, style = 'padding-right:0px;',
                                                                                         shiny::actionButton("save_uncorrected_pvca_plot", "Add to Zip File Export")
                                                                                  ),
                                                                                  column(4, style = 'padding-left:0px;',
                                                                                         downloadButton("dnldsave_uncorrected_pvca_plot","Dowload Single Plot")
                                                                                  )
                                                                                )
                                                                                ),
                                                                  shiny::column(6,
                                                                                shiny::h3("Batch Corrected PVCA"),
                                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_PVCA_plot")), type = 6),
                                                                                p(),
                                                                                shiny::fluidRow(
                                                                                  column(4, style = 'padding-right:0px;',
                                                                                         shiny::actionButton("save_corrected_pvca_plot", "Add to Zip File Export")
                                                                                  ),
                                                                                  column(4, style = 'padding-left:0px;',
                                                                                         downloadButton("dnldsave_corrected_pvca_plot","Dowload Single Plot")
                                                                                  )
                                                                                )

                                                                  )
                                                                ),
                                                                value = "pvca")
                                              ),
                                              value = "exp_var"
                                              ),
                                            #### SVA ---------------------------------
                                            shiny::tabPanel(
                                              "SVA Analysis",
                                              p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected SVA Analysis"),
                                                              shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_SVA_probability_density")), type = 6),
                                                              p(),
                                                              shiny::fluidRow(
                                                                column(4, style = 'padding-right:0px;',
                                                                       shiny::actionButton("save_uncorrected_SVA_probability_density", "Add to Zip File Export")
                                                                ),
                                                                column(4, style = 'padding-left:0px;',
                                                                       downloadButton("dnldsave_uncorrected_SVA_probability_density","Dowload Single Plot")
                                                                )
                                                              ),
                                                              p(),
                                                              verbatimTextOutput("uncorrected_SVA_nsv_print"),
                                                              ),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected SVA Analysis"),
                                                              shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_SVA_probability_density")), type = 6),
                                                              p(),
                                                              shiny::fluidRow(
                                                                column(4, style = 'padding-right:0px;',
                                                                       shiny::actionButton("save_corrected_SVA_probability_density", "Add to Zip File Export")
                                                                ),
                                                                column(4, style = 'padding-left:0px;',
                                                                       downloadButton("dnldsave_corrected_SVA_probability_density","Dowload Single Plot")
                                                                )
                                                              ),
                                                              p(),
                                                              verbatimTextOutput("corrected_SVA_nsv_print"),
                                                              )
                                                ),
                                              value = "sva"
                                              ),
                                            #### Box Plot ---------------------------------
                                            shiny::tabPanel(
                                              "Box Plot",
                                              p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected Box Plot"),
                                                              shiny::uiOutput("renduncorrected_Box_plot"),
                                                              p(),
                                                              shiny::fluidRow(
                                                                column(4, style = 'padding-right:0px;',
                                                                       shiny::actionButton("save_uncorrected_Box_plot", "Add to Zip File Export")
                                                                ),
                                                                column(4, style = 'padding-left:0px;',
                                                                       downloadButton("dnldsave_uncorrected_Box_plot","Dowload Single Plot")
                                                                )
                                                              )),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected Box Plot"),
                                                              shiny::uiOutput("rendcorrected_Box_plot"),
                                                              p(),
                                                              shiny::fluidRow(
                                                                column(4, style = 'padding-right:0px;',
                                                                       shiny::actionButton("save_corrected_Box_plot", "Add to Zip File Export")
                                                                ),
                                                                column(4, style = 'padding-left:0px;',
                                                                       downloadButton("dnldsave_corrected_Box_plot","Dowload Single Plot")
                                                                )
                                                              )

                                                )
                                              ),
                                              value = "box"
                                            )
                                          )
                                          )
                                      )
                                    )
                    )


# Read and Manipulate Input Files

# Server -----------------------------------------------------------------------
server <- function(input, output, session) {


  output$rendbatch_correction_method <- renderUI({

    if (input$RawCountCheck) {
      if (input$RawCountNorm == "none") {
        shiny::selectInput("batch_correction_method",
                           "Method of Batch Correction",
                           c("Select","ComBatseq"))
      } else {
        shiny::selectInput("batch_correction_method",
                           "Method of Batch Correction",
                           c("Select", "Limma", "ComBat", "Mean Centering", "Harman", "RUVg","SVA"))
      }
    } else {
      shiny::selectInput("batch_correction_method",
                         "Method of Batch Correction",
                         c("Select", "Limma", "ComBat", "Mean Centering", "Harman", "RUVg","SVA"))
    }


  })

  output$rendQuantNorm <- renderUI({

    if (input$RawCountCheck) {
      if (input$RawCountNorm != "none") {
        checkboxInput("QuantNorm","Quantile Normalize",value = T)
      }
    } else {
      checkboxInput("QuantNorm","Quantile Normalize",value = T)
    }

  })

  ## START OF UNCORRECTED PLOT DATA --------------------------------------------
  ### Input File Processing ----------------------------------------------------

  # reactive value
  files_to_download <- shiny::reactiveValues()

  #### Matrix Processing -------------------------------------------------------
  # developing reactive variable for the input expression matrix to be used downstream

  uncorrected_matrix <- reactiveVal()
  user_provided_batch_info <- reactiveVal()
  boxplot_meta_file <- reactiveVal()
  uncorr_boxplot_meta_file <- reactiveVal()
  corr_boxplot_meta_file <- reactiveVal()

  shiny::observeEvent(input$UseExpData, {

    # Load example matrix
    uncorrected_matrix <- readr::read_delim(Example_Matrix_File,
                                            delim = '\t',
                                            name_repair = "minimal",
                                            na = c("", "NA", "N/A"))
    uncorrected_matrix <- uncorrected_matrix[, !duplicated(colnames(uncorrected_matrix))]
    uncorrected_matrix(uncorrected_matrix)

    # Load example batch data
    user_provided_batch_info <- as.data.frame(readr::read_delim(Example_Meta_File,
                                                                delim = '\t',
                                                                na = c("", "NA", "N/A")))
    user_provided_batch_info(user_provided_batch_info)
    updateRadioButtons(session,"HumanOrMouse",selected = "Mouse")
    updateSelectInput(session,"Init_Batch_Select",selected = "Study")

  })

  # Observe if matrix uploaded by user
  shiny::observe({

    req(input$uncorrected_matrix_input)
    uncorrected_matrix_read <- readr::read_delim(input$uncorrected_matrix_input$datapath,
                                            delim = input$matrix_delim,
                                            name_repair = "minimal",
                                            na = c("", "NA", "N/A"))
    uncorrected_matrix_read <- uncorrected_matrix_read[, !duplicated(colnames(uncorrected_matrix_read))]
    uncorrected_matrix(uncorrected_matrix_read)
    updateRadioButtons(session,"HumanOrMouse",selected = "Human")

  })

  # Data processing for input matrix of gene expression data
  uncorrected_matrix_filtered <- shiny::reactive({
    req(uncorrected_matrix())
    uncorrected_matrix_unfiltered <- uncorrected_matrix()
    colnames(uncorrected_matrix_unfiltered)[1] <- "Genes"
    if (TRUE %in% duplicated(uncorrected_matrix_unfiltered[, 1])) {
      data_dup <- uncorrected_matrix_unfiltered %>% dplyr::group_by(Genes) %>% dplyr::filter(n() > 1) %>% as.data.frame()
      data_nodup <- uncorrected_matrix_unfiltered %>% dplyr::group_by(Genes) %>% dplyr::filter(n() == 1) %>% as.data.frame()
      data_dup_summ <- data_dup %>%
        dplyr::group_by(Genes) %>%
        dplyr::summarise_all(max) %>%
        as.data.frame()
      uncorrected_matrix_unfiltered <- rbind(data_nodup,data_dup_summ)
    } else {
      uncorrected_matrix_nodupes <- uncorrected_matrix_unfiltered
    }
    uncorrected_matrix_nodupes$filter <- apply(
      uncorrected_matrix_nodupes[,-1],
      1,
      function(x) ExprFilter2(x, 1, 0.05)
    )
    uncorrected_matrix_filtered <- uncorrected_matrix_nodupes[which(uncorrected_matrix_nodupes$filter == TRUE),]
    uncorrected_matrix_filtered <- uncorrected_matrix_filtered[,-ncol(uncorrected_matrix_filtered)]
    uncorrected_matrix_filtered <- as.data.frame(uncorrected_matrix_filtered)
    uncorrected_matrix_filtered
  })

  # preprocessing steps
  uncorrected_numeric_matrix <- shiny::reactive({
    req(uncorrected_matrix_filtered())
    uncorrected_numeric_matrix <- uncorrected_matrix_filtered()
    rownames(uncorrected_numeric_matrix) <- uncorrected_numeric_matrix[,1]
    uncorrected_numeric_matrix <- uncorrected_numeric_matrix[,-1]

    if (input$Log_Choice){
      uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix + 1)
    }
    uncorrected_numeric_matrix <- cbind(uncorrected_numeric_matrix, uncorrected_matrix_filtered()[,1,drop=FALSE])
    uncorrected_numeric_matrix <- uncorrected_numeric_matrix %>% dplyr::relocate(Genes)
    uncorrected_numeric_matrix
  })
  # rendering input table for visualization
  output$uncorrected_matrix_output <- DT::renderDataTable({
    req(uncorrected_numeric_matrix())
    DT::datatable(uncorrected_numeric_matrix(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)

  })
  output$rendMatrixHeader <- renderUI({
    req(uncorrected_numeric_matrix())
    h3("Uncorrected Matrix Preview")
  })
  output$uncorrected_matrix_output_input <- DT::renderDataTable({
    req(uncorrected_numeric_matrix())
    DT::datatable(uncorrected_numeric_matrix(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })

  #### Batch Meta Processing ---------------------------------------------------
  output$user_batch_info_file <- shiny::renderUI({

    if (input$batch_info == "No") {
      shiny::fluidRow(
        shiny::column(9,
                      shiny::fileInput("user_provided_batch_info", "Please Provide Meta File")
        ),
        shiny::column(3,
                      shiny::selectInput("meta_delim", "Deliminator", c("Tab" = '\t',"Comma" = ',', "Space" = ' ',"Semicolon" = ';', "Colon" = ':'))
        )
      )
    }

  })

  # determining if meta table should be derived or inputted
  output$batch_delim <- shiny::renderUI({
    if(input$batch_info == "Yes"){
      shiny::textInput("batch_name_delim",
                       "Deliminator to Derive Batch Information from Column Names")
    }
  })

  # User to provide batch names if deliminator is entered
  output$batch_names <- shiny::renderUI({
    if(shiny::isTruthy(input$batch_name_delim)){
      if(input$batch_info == "Yes"){
        shiny::textInput("batch_names",
                         "Please Provide Batch Names Separated by Commas")
      }
    }
  })

  # creating variable with user provided delim and batch names
  meta_file_delim <- shiny::reactive({
    print(input$batch_name_delim)
  })
  meta_file_names <- shiny::reactive({
    print(unlist(strsplit(input$batch_names, ",")))
  })

  shiny::observe({

    req(input$user_provided_batch_info)
    # Reading user provided meta file
    user_provided_batch_info <- as.data.frame(readr::read_delim(input$user_provided_batch_info$datapath,
                                                                delim = input$meta_delim,
                                                                na = c("", "NA", "N/A")))
    user_provided_batch_info(user_provided_batch_info)

  })

  # Attaching and rendering user generated or user provided meta file
  meta_file_react <- shiny::reactive({
    if (input$batch_info == "Yes"){
      meta_names <- colnames(uncorrected_matrix())
      meta_names <- as.data.frame(meta_names[-1])
      meta_names <- meta_names %>% dplyr::rename( "Original Sample Name" = "meta_names[-1]")
      meta_names <- cbind(
        meta_names,
        tidyr::separate_wider_delim(
          meta_names,
          cols = 1,
          delim = meta_file_delim(),
          names = meta_file_names()
        )
      )
      meta_names
    } else if (input$batch_info == "No"){
      user_provided_batch_info()
    }
  })

  aligned_meta_file <- shiny::reactive({
    req(uncorrected_numeric_matrix())
    req(meta_file_react())
    unaligned_meta_file <- meta_file_react()
    # getting meta column name that contains sample names
    column_containing_ids <- colnames(unaligned_meta_file)[1]
    # filter to samples in matrix
    unaligned_meta_file_filtered <- unaligned_meta_file[which(unaligned_meta_file[,column_containing_ids] %in% colnames(uncorrected_numeric_matrix())),]
    # align sample names in meta file with matrix, make sure same order
    aligned_meta_file_filtered <- unaligned_meta_file_filtered[match(colnames(uncorrected_numeric_matrix()[,-1]), unaligned_meta_file_filtered[,column_containing_ids]),]
    aligned_meta_file_filtered <- as.data.frame(aligned_meta_file_filtered)
    aligned_meta_file_filtered <- aligned_meta_file_filtered %>% relocate(any_of(column_containing_ids))
    boxplot_meta_file(aligned_meta_file_filtered)
    uncorr_boxplot_meta_file(aligned_meta_file_filtered)
    corr_boxplot_meta_file(aligned_meta_file_filtered)
    aligned_meta_file_filtered
  })

  output$rendMetaHeader <- renderUI({
    req(uncorrected_matrix())
    req(user_provided_batch_info())
    h3("Meta Preview")
  })
  output$meta_file <- DT::renderDataTable({
    req(uncorrected_matrix())
    DT::datatable(aligned_meta_file(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)

  })

  batch_names_from_meta <- shiny::reactive({
    batch_names_from_meta <- colnames(aligned_meta_file())
    batch_names_from_meta
  })

  #### PCA ---------------------------------------------------------------------
  # Generating PCA plot
  output$rendbatch_choices_PCA <- shiny::renderUI({
    shiny::selectInput("batch_choices_PCA",
                       "Group Color",
                       choices = c("Select", unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       selected = batch1_choices())
  })

  uncorrected_batch_choices_PCA2 <- shiny::reactive({
    if (isTruthy(input$batch_choices_PCA)) {
      if (input$batch_choices_PCA != "Select"){
        uncorrected_batch_choices_PCA2 <- input$batch_choices_PCA
        uncorrected_batch_choices_PCA2
      }else {
        uncorrected_batch_choices_PCA2 <- NULL
        uncorrected_batch_choices_PCA2
      }
    }

  })

  output$rendPCAhover <- shiny::renderUI({
    batch_choice <- input$batch_choices_PCA
    shiny::selectInput("PCAhover",
                       "Information to Display on Hover:",
                       choices = unlist(strsplit(batch_names_from_meta(), ",")),
                       selected = c(batch_names_from_meta()[1],batch_choice),
                       multiple = T)
  })

  PCA_data <- shiny::reactive({
    PCA_data <- cbind((t(uncorrected_numeric_matrix()[,-1])), aligned_meta_file())
    PCA_data
  })

  uncorrected_PCA_react <- reactive({

    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "PCA Object")
      mat <- uncorrected_numeric_matrix()
      rownames(mat) <- mat[,1]
      if (input$PCA_type == "Cluster Annotation"){
        #cluster::pam(as.data.frame(t(uncorrected_numeric_matrix()[,-1])), input$cluster_number)
        pca <- cluster::pam(as.data.frame(t(mat[,-1])), input$cluster_number)
      } else if (input$PCA_type == "Meta Annotation"){
        #stats::prcomp(as.data.frame(t(uncorrected_numeric_matrix()[,-1])), scale. = TRUE)
        pca <- stats::prcomp(as.data.frame(t(mat[,-1])), scale. = TRUE)
      }
      incProgress(0.5, detail = "Complete!")
    })
    pca

  })

  uncorrected_PCA_plot_df <- reactive({

    meta <- aligned_meta_file()
    pca <- uncorrected_PCA_react()
    batch_choice <- input$batch_choices_PCA
    hover_choice <- input$PCAhover
    metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice))
    metaCols <- metaCols[which(metaCols %in% colnames(meta))]
    PC_x <- "PC1"
    PC_y <- "PC2"
    uncorrected_PCA_points <- as.data.frame(pca$x[,c(1,2)])
    uncorrected_PCA_plot_df <- merge(uncorrected_PCA_points,meta[,metaCols,drop = F],by.x = 0, by.y = colnames(meta)[1])
    colnames(uncorrected_PCA_plot_df)[1] <- colnames(meta)[1]
    PCX <- gsub("^PC","",PC_x)
    PCY <- gsub("^PC","",PC_y)
    uncorrected_PCA_plot_df[,PC_x] <- as.numeric(uncorrected_PCA_plot_df[,PC_x]) / (pca$sdev[as.numeric(PCX)] * sqrt(nrow(uncorrected_PCA_plot_df)))
    uncorrected_PCA_plot_df[,PC_y] <- as.numeric(uncorrected_PCA_plot_df[,PC_y]) / (pca$sdev[as.numeric(PCY)] * sqrt(nrow(uncorrected_PCA_plot_df)))
    uncorrected_PCA_plot_df$text <- NA
    for (r in seq(nrow(uncorrected_PCA_plot_df))) {
      text_vec <- c()
      for (c in metaCols) {
        text_vec <- c(text_vec,
                      paste0("\\<br\\>\\<b\\>", c, ":\\</b\\> ",uncorrected_PCA_plot_df[r,c]))
      }
      text = paste(text_vec,collapse = '')
      text = gsub("\\\\", '', text)
      text = paste0(text,"<extra></extra>") #removes outside of box text
      uncorrected_PCA_plot_df$text[r] = text
    }
    uncorrected_PCA_plot_df

  })

  uncorrected_PCA_plot_forDlnd_react <- shiny::reactive({
    if (input$PCA_type == "Cluster Annotation"){
      autoplot(
        uncorrected_PCA_react(),
        frame = TRUE,
        frame.type = 'norm')
    } else if (input$PCA_type == "Meta Annotation"){
      autoplot(
        uncorrected_PCA_react(),
        data = PCA_data(),
        color = uncorrected_batch_choices_PCA2()
      )
    }
  })

  uncorrected_PCA_plot_react <- shiny::reactive({

    meta <- aligned_meta_file()
    plot_df <- uncorrected_PCA_plot_df()
    batch_choice <- input$batch_choices_PCA
    hover_choice <- input$PCAhover
    metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice))
    PC_x <- "PC1"
    PC_y <- "PC2"
    colnames(plot_df)[which(colnames(plot_df) == PC_x)] <- "X"
    colnames(plot_df)[which(colnames(plot_df) == PC_y)] <- "Y"
    dotSize <- input$PCAdotSize
    fontSize <- input$pcaFontSize

    if (isTruthy(batch_choice)) {
      if (batch_choice != "Select") {
        colnames(plot_df)[which(colnames(plot_df) == batch_choice)] <- "ColorBatch"
        p4 <- plot_ly() %>%
          add_trace(
            data = plot_df,
            x = ~X,
            y = ~Y,
            type = "scatter",
            mode = "markers",
            color = ~ColorBatch,
            legendgroup=batch_choice,
            marker=list(size=dotSize,symbol = 'circle'),
            hovertemplate = ~text
          ) %>%
          layout(xaxis = list(zeroline = FALSE,title = PC_x),
                 yaxis = list(zeroline = FALSE,title = PC_y),
                 font = list(size = fontSize)) %>%
          config(
            toImageButtonOptions = list(
              format = "svg"
            )
          )
        p4
      } else {
        p4 <- plot_ly() %>%
          add_trace(
            data = plot_df,
            x = ~X,
            y = ~Y,
            type = "scatter",
            mode = "markers",
            marker = list(color = "black"),
            marker=list(size=dotSize,symbol = 'circle'),
            hovertemplate = ~text
          ) %>%
          layout(xaxis = list(zeroline = FALSE,title = PC_x),
                 yaxis = list(zeroline = FALSE,title = PC_y),
                 font = list(size = fontSize)) %>%
          config(
            toImageButtonOptions = list(
              format = "svg"
            )
          )
        p4
      }
    }


  })

  output$uncorrected_PCA <- renderPlotly({
    uncorrected_PCA_plot_react()
  })


  #### PCA MC ------------------------------------------------------------------
  # Generating multiple PC plot
  output$PCA_mc_color_choice <- shiny::renderUI({
    shiny::selectInput("PCA_mc_color_choice",
                       "Group color",
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  uncorrected_PCA_mc_color_choice2 <- shiny::reactive({
    if (isTruthy(input$PCA_mc_color_choice)) {
      if (input$PCA_mc_color_choice == "Select"){
        NULL
      }else {
        uncorrected_PCA_mc_color_choice2 <- input$PCA_mc_color_choice
      }
    }
  })
  uncorrected_PCA_multiple_components <- shiny::reactive({
    mat_PCA_mc <- uncorrected_numeric_matrix()[,-1]
    names(mat_PCA_mc) <- NULL
    uncorrected_PCA_mc_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_PCA_mc)),
      colData = aligned_meta_file(),
      rowData = uncorrected_numeric_matrix()[,1]
    )
    SummarizedExperiment::assay(uncorrected_PCA_mc_SCE, "logcounts") <- SingleCellExperiment::counts(uncorrected_PCA_mc_SCE)
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.25, detail = "Running multiple components PCA")
      uncorrected_PCA_mc_SCE <- scater::runPCA(uncorrected_PCA_mc_SCE, ncomponents = 50)
      incProgress(0.25, detail = "Plotting multiple components PCA")
      p <- scater::plotPCA(
        uncorrected_PCA_mc_SCE,
        ncomponents = input$PCA_mc_slider,
        colour_by = uncorrected_PCA_mc_color_choice2()
      )
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$uncorrected_PCA_multiple_components <- shiny::renderPlot({
    uncorrected_PCA_multiple_components()
  })

  #### PCA Details--------------------------------------------------------------
  # Uncorrected PCA details plots
  uncorrected_PCA_details <- shiny::reactive({
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running FactoMineR PCA")
      pc_obj <- FactoMineR::PCA(t(uncorrected_numeric_matrix()[,-1]), graph = F)
      incProgress(0.5, detail = "Complete!")
    })
    pc_obj
  })
  uncorrected_PCA_details2 <- shiny::reactive({
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running stats prcomp PCA")
      pc_obj <- stats::prcomp(t(uncorrected_numeric_matrix()[,-1]))
      incProgress(0.5, detail = "Complete!")
    })
    pc_obj
  })
  uncorrected_scree_plot_react <- shiny::reactive({
    uncorrected_scree_eig <- as.data.frame(get_eig(uncorrected_PCA_details()))
    uncorrected_scree_eig$Dimensions <- gsub("Dim\\.","",rownames(uncorrected_scree_eig))
    uncorrected_scree_eig$`Variance Percent` <- paste0(round(uncorrected_scree_eig$variance.percent,1),"%")
    uncorrected_scree_eig_top10 <- uncorrected_scree_eig[1:10,]
    AxisTickFont <- input$pcaDetAxisTkSize
    AxisTitleFont <- input$pcaDetAxisTtSize
    LabelFont <- input$pcaDetLabelSize

    p <- ggplot(data=uncorrected_scree_eig_top10, aes(x=reorder(Dimensions,-variance.percent), y=variance.percent)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_minimal() +
      geom_text(aes(label=`Variance Percent`), vjust=-0.3,size=LabelFont) +
      labs(x = "Dimensions",y = "Variance Percent") +
      theme(axis.text = element_text(size = AxisTickFont),
            axis.title = element_text(size = AxisTitleFont))
    p
  })
  output$uncorrected_scree_plot <- shiny::renderPlot({
    uncorrected_scree_plot_react()
  })
  output$PCA_factors_choices <- shiny::renderUI({
    shiny::selectInput("PCA_factors_choices",
                       "Select Factor for PCA Details",
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))))
  })
  uncorrected_PCA_factors_choices2 <- shiny::reactive({
    if (isTruthy(input$PCA_factors_choices)) {
      if (input$PCA_factors_choices == "Select"){
        NULL
      }else {
        uncorrected_PCA_factors_choices2 <- input$PCA_factors_choices
      }
    }
  })
  uncorrected_PCA_individuals <- shiny::reactive({
    uncorrected_PCA_individuals <- uncorrected_PCA_details()
    uncorrected_PCA_factors <- as.data.frame(uncorrected_PCA_individuals$ind$contrib)
    uncorrected_PCA_factors$factors <- as.vector(aligned_meta_file()[,input$PCA_factors_choices])
    uncorrected_PCA_factors_final_sum <- stats::aggregate(. ~ factors, uncorrected_PCA_factors, sum)
    uncorrected_PCA_factors_final_mean <- stats::aggregate(. ~ factors, uncorrected_PCA_factors, mean)
    uncorrected_PCA_factors_final_sdv <- stats::aggregate(. ~ factors, uncorrected_PCA_factors, sd)
    uncorrected_PCA_factors_final_longer_sum <- tidyr::pivot_longer(uncorrected_PCA_factors_final_sum, !factors, names_to = "PC_Components")
    uncorrected_PCA_factors_final_longer_mean <- tidyr::pivot_longer(uncorrected_PCA_factors_final_mean, !factors, names_to = "PC_Components")
    uncorrected_PCA_factors_final_longer_sdv <- tidyr::pivot_longer(uncorrected_PCA_factors_final_sdv, !factors, names_to = "PC_Components")
    colnames(uncorrected_PCA_factors_final_longer_sum)[3] <- "Sum"
    colnames(uncorrected_PCA_factors_final_longer_mean)[3] <- "Mean"
    colnames(uncorrected_PCA_factors_final_longer_sdv)[3] <- "SDV"
    uncorrected_PCA_factors_final_longer <- uncorrected_PCA_factors_final_longer_sum[,1:2]
    uncorrected_PCA_factors_final_longer <- cbind(
      uncorrected_PCA_factors_final_longer,
      uncorrected_PCA_factors_final_longer_sum[,3],
      uncorrected_PCA_factors_final_longer_mean[,3],
      uncorrected_PCA_factors_final_longer_sdv[,3]
    )
    colnames(uncorrected_PCA_factors_final_longer)[1] <- input$PCA_factors_choices
    uncorrected_PCA_factors_final_longer
  })
  output$uncorrected_contribution_table <- DT::renderDataTable({
    if (isTruthy(input$PCA_factors_choices)) {
      DT::datatable(uncorrected_PCA_individuals(),
                    options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                   pageLength = 10,
                                   scrollX = T),
                    rownames = F) %>%
        formatRound(columns = c(3:ncol(uncorrected_PCA_individuals())), digits = 4)
    }

  })
  uncorrected_contribution_counts <- shiny::reactive({
    uncorrected_PCA_individuals <- uncorrected_PCA_details()
    uncorrected_PCA_factors <- as.data.frame(uncorrected_PCA_individuals$ind$contrib)
    uncorrected_PCA_factors$factors <- as.vector(aligned_meta_file()[,input$PCA_factors_choices])
    uncorrected_PCA_factors_final_count <- uncorrected_PCA_factors %>% dplyr::group_by(factors) %>% dplyr::summarise(individuals = length(factors))
    colnames(uncorrected_PCA_factors_final_count)[1] <- input$PCA_factors_choices
    uncorrected_PCA_factors_final_count
  })
  output$uncorrected_contribution_counts <- DT::renderDataTable({
    if (isTruthy(input$PCA_factors_choices)) {
      DT::datatable(uncorrected_contribution_counts(),
                    options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                   pageLength = 10,
                                   scrollX = T),
                    rownames = F)
    }

  })

  #### UMAP --------------------------------------------------------------------

  output$rendUMAPImmuneDeconvMethods <- renderUI({

    if (input$UMAPFeatureCategory == "Immune Deconvolution Features") {
      if (input$HumanOrMouse == "Human") {
        selectInput("UMAPImmuneDeconvMethods","Immune Deconvolution Methods:",
                    c("mcp_counter","estimate","quantiseq","xcell","epic","abis"))
      } else {
        selectInput("UMAPImmuneDeconvMethods","Immune Deconvolution Methods:",
                    #c("mmcp_counter","seqimmucc","dcq","base"))
                    c("mmcp_counter","dcq","base"))
      }
    }

  })

  UMAP_ImmDeconv_uncorr_react <- reactive({

    req(input$UMAPImmuneDeconvMethods)
    deconvMethod <- input$UMAPImmuneDeconvMethods
    mat <- uncorrected_numeric_matrix()
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running Immune Deconvolution")
      if (input$HumanOrMouse == "Human") {
        deconv <- as.data.frame(deconvolute(mat, deconvMethod))
      } else {
        deconv <- as.data.frame(deconvolute_mouse(mat, deconvMethod))
      }
      incProgress(0.5, detail = "Complete!")
    })
    deconv


  })

  UMAP_PCA_Proj_uncorr_Samples <- reactive({

    if (input$uncorrected_panel == "pca_main") {
      if (input$PCA_main_pan == "umap") {
        withProgress(message = "Processing Uncorrected", value = 0, {
          incProgress(0.5, detail = "Running PCA")
          mat_uncorr <- uncorrected_numeric_matrix()
          rownames(mat_uncorr) <- mat_uncorr[,1]
          mat_uncorr <- mat_uncorr[,-1]
          mat_uncorr_t <- as.data.frame(t(mat_uncorr))
          pca_uncorr <- prcomp(t(mat_uncorr_t), scale = TRUE)
          rot_uncorr <- as.data.frame(t(pca_uncorr[["rotation"]]))
          rot_uncorr$`Principal Component` <- rownames(rot_uncorr)
          rot_uncorr <- rot_uncorr %>% relocate(`Principal Component`)
          incProgress(0.5, detail = "Complete!")
        })
        rot_uncorr
      }
    }
  })

  UMAP_Feature_Choices <- reactive({

    req(aligned_meta_file())
    Features <- NULL
    FeatCat <- input$UMAPFeatureCategory
    if (FeatCat == "Matrix Features") {
      Features <- uncorrected_numeric_matrix()[,1]
      Features
    } else if (FeatCat == "Meta Features") {
      Features <- colnames(aligned_meta_file())[-1]
      Features
    } else if (FeatCat == "PCA Projections") {
      Features <- UMAP_PCA_Proj_uncorr_Samples()[,1]
      Features
    } else if (FeatCat == "Immune Deconvolution Features") {
      Features <- UMAP_ImmDeconv_uncorr_react()[,1]
      Features
    }

  })

  shiny::observe({

    updateSelectizeInput(session = session, inputId = "UMAPFeatSelection",
                         choices = UMAP_Feature_Choices(),
                         server = T)

  })

  output$UMAPGeneSetTableUI <- DT::renderDataTable({

    GeneSetTable_sub <- GeneSetTableBack_react()
    DT::datatable(GeneSetTable_sub,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)

  })

  UMAP_uncorr_anno_df <- shiny::reactive({

    FeatCat <- input$UMAPFeatureCategory
    Feature <- input$UMAPFeatSelection
    meta <- aligned_meta_file()
    NameCol <- colnames(meta)[1]
    if (FeatCat == "Matrix Features") {
      mat <- uncorrected_numeric_matrix()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
      meta
    } else if (FeatCat == "PCA Projections") {
      mat <- UMAP_PCA_Proj_uncorr_Samples()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
      meta
    } else if (FeatCat == "Immune Deconvolution Features") {
      mat <- UMAP_ImmDeconv_uncorr_react()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
      meta
    } else if (FeatCat == "Gene Set Pathways") {
      gs_name <- GeneSetTableBack_react()[input$UMAPGeneSetTableUI_rows_selected,ncol(GeneSetTableBack_react())]
      gs <- geneset[gs_name]
      mat <- uncorrected_numeric_matrix()
      rownames(mat) <- mat[,1]
      mat <- mat[,-1]
      mat <- as.matrix(mat)
      withProgress(message = "Processing Uncorrected", value = 0, {
        incProgress(0.5, detail = "Running ssGSEA")
        ssGSEA_param <- GSVA::ssgseaParam(mat,gs)
        ssGSEA <- GSVA::gsva(ssGSEA_param)
        incProgress(0.5, detail = "Complete!")
      })
      ssGSEA <- as.data.frame(t(ssGSEA))
      ssGSEA[,NameCol] <- rownames(ssGSEA)
      meta <- merge(meta,ssGSEA)
      meta
    } else {
      meta <- meta
      meta
    }

  })

  uncorrected_umap_coord <- reactive({

    AnnoMeta <- UMAP_uncorr_anno_df()
    uncorrected_PCA_react_mat <- as.matrix(uncorrected_PCA_react()$x)
    umapNN <- input$UMAPnnb
    umapMinDist <- input$UMAPminDist
    umapMetric <- input$UMAPmetricSelec
    if (umapNN >= nrow(uncorrected_PCA_react_mat)) {
      umapNN <- nrow(uncorrected_PCA_react_mat) - 1
    }
    if (all(isTruthy(c(umapMetric,umapMinDist,umapNN)))) {

      otherMets <- c("euclidean","manhattan","cosine","pearson","pearson2")

      if (umapMetric == "hamming" || umapMetric == "correlation") {
        tdata_fit_df <- as.data.frame(uwot::umap(uncorrected_PCA_react_mat,metric = umapMetric, n_neighbors = umapNN, min_dist = umapMinDist))
      } else if (umapMetric %in% otherMets) {
        umap.mod$metric <- umapMetric
        umap.mod$min_dist <- umapMinDist
        umap.mod$n_neighbors <- umapNN
        withProgress(message = "Processing Uncorrected", value = 0, {
          incProgress(0.5, detail = "Running UMAP")
          tdata_fit <- umap::umap(uncorrected_PCA_react_mat,config = umap.mod)
          incProgress(0.5, detail = "Complete!")
        })
        tdata_fit_df <- as.data.frame(tdata_fit$layout)
      }
      colnames(tdata_fit_df) <- c("UMAP1","UMAP2")
      tdata_fit_df <- tdata_fit_df %>%
        mutate(ID=row_number())
      tdata_fit_df$SampleName <- rownames(tdata_fit_df)
      tdata_fit_df <- tdata_fit_df %>%
        relocate(SampleName)
      tdata_fit_df
    }

  })

  uncorrected_UMAP_react <- reactive({

    req(input$UMAPFeatSelection)
    plot_df <- uncorrected_umap_coord()
    rownames(plot_df) <- plot_df[,1]
    UMAPdotSize <- 2
    umap_annoCol <- input$UMAPFeatSelection
    meta <- UMAP_uncorr_anno_df()
    NameCol <- colnames(meta)[1]
    umapdotSize <- input$umapdotSize
    umapAxisTkSize <- input$umapAxisTkSize
    umapAxisTtSize <- input$umapAxisTtSize

    plot_df <- merge(plot_df,meta[,c(NameCol,umap_annoCol)], by.x = colnames(plot_df)[1], by.y = NameCol)

    if (is.null(umap_annoCol)) {
      k <- plot_df %>%
        ggplot(aes(UMAP1, UMAP2,
                   text = paste("</br><b>Sample Name:</b> ", SampleName,
                                sep = "")))
    } else {
      k <- plot_df %>%
        ggplot(aes(UMAP1, UMAP2, color=!!sym(umap_annoCol),
                   text = paste("</br><b>Sample Name:</b> ", SampleName,
                                "</br><b>",umap_annoCol,":</b> ", !!sym(umap_annoCol),
                                sep = "")))
    }

    k <- k + geom_point(shape = 19, size = umapdotSize) +
      theme_minimal()

    k <- k + theme(axis.text.x = element_text(size = umapAxisTkSize),
                   axis.title.x = element_text(size = umapAxisTtSize),
                   axis.text.y = element_text(size = umapAxisTkSize),
                   axis.title.y = element_text(size = umapAxisTtSize))

    if (is.numeric(plot_df[,umap_annoCol])) {
      k <- k + scale_colour_gradient(low = "#56B1F7",high = "#132B43")
    }
    k


  })
  output$uncorrected_UMAP <- plotly::renderPlotly({

    p <- uncorrected_UMAP_react()
    ply <- ggplotly(p, tooltip = "text")
    ply

  })

  #### Cluster -----------------------------------------------------------------
  # uncorrected cluster analysis
  cluster_mv_features_uncorr_matrix <- reactive({
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$cluster_n_MV_features
      var_type <- input$VarianceMeasure
      mat <- uncorrected_numeric_matrix()
      featColName <- colnames(mat)[1]
      rownames(mat) <- mat[,1]
      mat <- mat[,-1]
      mad <- NULL
      var <- NULL
      cv <- NULL
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(log2(mat + 1), 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(log2(mat + 1), 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(log2(mat + 1), 1, cv))
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = topN)
        out <- cbind(names(cv), cv[names(cv)], mat[names(cv),])
        colnames(out) <- c("Gene", "CV", colnames(mat))
        dataset <- mat[names(cv),]
      }
      dataset[,featColName] <- rownames(dataset)
      dataset <- dataset %>% dplyr::relocate(any_of(featColName))
      incProgress(0.5, detail = "Complete!")
    })
    dataset

  })
  uncorrected_elbow_analysis <- reactive({
    req(cluster_mv_features_uncorr_matrix())
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running Elbow Analysis")
      uncorrected_elbow_analysis <- fviz_nbclust(x = t(cluster_mv_features_uncorr_matrix()[,-1]), kmeans, method = "wss",verbose = T)
      incProgress(0.5, detail = "Complete!")
    })
    uncorrected_elbow_analysis
  })
  output$uncorrected_elbow_plot <- renderPlot({
    uncorrected_elbow_analysis()
  })
  uncorrected_silhouette_analysis <- reactive({
    req(cluster_mv_features_uncorr_matrix())
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running Silhouette Analysis")
      uncorrected_silhouette_analysis <- fviz_nbclust(x = t(cluster_mv_features_uncorr_matrix()[,-1]), kmeans, method = "silhouette",verbose = T)
      incProgress(0.5, detail = "Complete!")
    })
    uncorrected_silhouette_analysis
  })
  output$uncorrected_silhouette_plot <- renderPlot({
    uncorrected_silhouette_analysis()
  })
  uncorrected_dunn_index_analysis <- reactive({
    req(cluster_mv_features_uncorr_matrix())
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running Dunn Index Analysis")
      dunn_k <- c(2:10)
      dunnin <- c()
      for (i in dunn_k){
        dunnin[i] <- dunn(
          distance = dist(t(cluster_mv_features_uncorr_matrix()[,-1])),
          clusters = kmeans(t(cluster_mv_features_uncorr_matrix()[,-1]), i)$cluster
        )
      }
      uncorrected_dunn_index_analysis <- as.data.frame(cbind(dunn_k, dunnin[-1]))
      colnames(uncorrected_dunn_index_analysis) <- c("cluster_number", "dunn_index")
      p <- ggplot(data = uncorrected_dunn_index_analysis, mapping = aes(x = cluster_number, y = dunn_index))+
        geom_point(color = "dodgerblue1")+
        geom_line(color = "dodgerblue1")+
        geom_vline(
          xintercept = uncorrected_dunn_index_analysis$cluster_number[which(max(uncorrected_dunn_index_analysis$dunn_index) == uncorrected_dunn_index_analysis$dunn_index)],
          color = "dodgerblue1",
          linetype = 2
        ) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14))
      incProgress(0.5, detail = "Complete!")
    })
    p

  })
  output$uncorrected_dunn_index_plot <- renderPlot({
    uncorrected_dunn_index_analysis()
  })


  #### Heatmap -----------------------------------------------------------------
  output$rendHeatmapAnnoSel <- renderUI({

    shiny::selectInput("HeatmapAnnoSel",
                       "Column Annoation:",
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       selected = batch1_choices(),
                       multiple = T)

  })

  heat_colAnn <- reactive({

    if (isTruthy(input$HeatmapAnnoSel)) {
      meta <- aligned_meta_file()
      rownames(meta) <- meta[,1]
      meta_sub <- meta[,input$HeatmapAnnoSel, drop = F]
      colAnn <- ComplexHeatmap::HeatmapAnnotation(df = meta_sub,
                                                  name = input$HeatmapAnnoSel,
                                                  which = 'col'
      )
      colAnn
    } else {
      colAnn <- NULL
      colAnn
    }

  })


  # Generating heatmap of uncorrected data
  uncorrected_heatmap <- shiny::reactive({
    req(cluster_mv_features_uncorr_matrix())

    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Generating Heatmap")

      colAnn <- heat_colAnn()
      HeatRowNames <- ifelse("Turn on Row Names" %in% input$HeatRowColNames,TRUE,FALSE)
      HeatColNames <- ifelse("Turn on Column Names" %in% input$HeatRowColNames,TRUE,FALSE)

      uncorrected_matrix_heatmap <- as.matrix(cluster_mv_features_uncorr_matrix()[,-1])
      uncorrected_matrix_heatmap_cols <- colnames(uncorrected_matrix_heatmap)
      #uncorrected_matrix_heatmap <- as.matrix(uncorrected_numeric_matrix()[,-1])
      #uncorrected_matrix_heatmap <- as.matrix(log2(uncorrected_numeric_matrix()[,-1] + 1))
      uncorrected_matrix_scaled <- t(apply(uncorrected_matrix_heatmap, 1, scale))
      colnames(uncorrected_matrix_scaled) <- uncorrected_matrix_heatmap_cols
      p <- suppressMessages(ComplexHeatmap::Heatmap(uncorrected_matrix_scaled, top_annotation = colAnn,
                                                    show_row_names = HeatRowNames, show_column_names = HeatColNames,
                                                    heatmap_legend_param = list(title = "Expression")))
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$uncorrected_heatmap <- shiny::renderPlot({
    uncorrected_heatmap()
  })

  #### RLE ---------------------------------------------------------------------
  # Generating uncorrected RLE plot
  output$batch_choices_RLE <- shiny::renderUI({
    shiny::selectInput("batch_choices_RLE",
                       "Group Color by Batch",
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))),
                       selected = batch1_choices())
  })
  uncorrected_batch_choices_RLE2 <- shiny::reactive({
    if (isTruthy(input$batch_choices_RLE)) {
      if (input$batch_choices_RLE == "Select"){
        NULL
      }else {
        uncorrected_batch_choices_RLE2 <- input$batch_choices_RLE
      }
    }
  })
  RLE_Obj_Uncorr <- reactive({

    mat_RLE <- uncorrected_numeric_matrix()[,-1]
    names(mat_RLE) <- NULL
    uncorrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_RLE)),
      colData = aligned_meta_file()[order(aligned_meta_file()$Study),],
      rowData = uncorrected_numeric_matrix()[,1]
    )
    uncorrected_RLE_SCE

  })
  uncorrected_RLE <- shiny::reactive({
    rle_obj <- RLE_Obj_Uncorr()
    req(rle_obj)
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running Relative Log Expression")
      if (isTruthy(input$batch_choices_RLE)) {
        p <- scater::plotRLE(
          rle_obj,
          exprs_values = "counts",
          color_by = uncorrected_batch_choices_RLE2()
        )
      }
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$uncorrected_RLE_plot <- shiny::renderPlot({
    uncorrected_RLE()
  })

  #### Exp Var -----------------------------------------------------------------
  # Generating EV plot
  output$variable_choices_EV <- shiny::renderUI({
    shiny::selectInput("variable_choices_EV",
                       "Select Variables to Plot",
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       selected = head(c(unlist(strsplit(batch_names_from_meta()[-1], ","))), n = 3),
                       #c(unlist(strsplit(batch_names_from_meta(), ","))),
                       multiple = TRUE
    )
  })
  output$rendExpPlotLog <- renderUI({

    if (input$ExpVar_Plots == "exp_plot") {
      checkboxInput("ExpPlotLog","Log Scale X-Axis", value = T)
    }

  })

  output$rendpvcaPct <- renderUI({

    if (input$ExpVar_Plots == "pvca") {
      numericInput("pvcaPct","% Threshold", value = 0.8, min = 0, max = 1)
    }

  })
  uncorrected_EV_df <- shiny::reactive({
    req(input$variable_choices_EV)
    my_colors <- metafolio::gg_color_hue(length(input$variable_choices_EV))
    mat_EV <- uncorrected_numeric_matrix()[,-1]
    names(mat_EV) <- NULL
    uncorrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_EV)),
      colData = aligned_meta_file(),
      rowData = uncorrected_numeric_matrix()[,1]
    )
    SummarizedExperiment::assay(uncorrected_EV_SCE, "logcounts") <- SingleCellExperiment::counts(uncorrected_EV_SCE)
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.25, detail = "Running PCA")
      uncorrected_EV_SCE_PCA <- scater::runPCA(uncorrected_EV_SCE)
      incProgress(0.25, detail = "Calculating per-gene variance")
      exp_mat_uncorr <- getVarianceExplained(uncorrected_EV_SCE_PCA,
                                             exprs_values = "logcounts",
                                             variables = input$variable_choices_EV)
      exp_mat_uncorr_melt <- reshape2::melt(exp_mat_uncorr)
      incProgress(0.5, detail = "Complete!")
    })
    exp_mat_uncorr_melt
  })

  uncorrected_EV <- reactive({

    exp_mat_uncorr_melt <- uncorrected_EV_df()
    p <- ggplot(exp_mat_uncorr_melt, aes(x = value,color = Var2)) +
      geom_density(size = 1)
    if (input$ExpPlotLog) {
      p <- p + scale_x_log10(limit = c(0.0001,100),labels = ~ format(.x, scientific = FALSE), breaks = c(0.001,0.01,0.1,1,10,100)) +
        geom_vline(xintercept = 1, linetype="dashed")
    }
    p <- p +
      theme_classic() +
      labs(x = "% Variance Explained", y = "Density", color = "Variable")
    p

  })


  output$uncorrected_EV_plot <- shiny::renderPlot({
    uncorrected_EV()
  })

  #### PVCA --------------------------------------------------------------------

  pvca_mv_features_uncorr_matrix <- reactive({

    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$pvcacluster_n_MV_features
      var_type <- input$pvcaVarianceMeasure
      mat <- uncorrected_numeric_matrix()
      featColName <- colnames(mat)[1]
      rownames(mat) <- mat[,1]
      mat <- mat[,-1]

      mad <- NULL
      var <- NULL
      cv <- NULL
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(log2(mat + 1), 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(log2(mat + 1), 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(log2(mat + 1), 1, cv))
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = topN)
        out <- cbind(names(cv), cv[names(cv)], mat[names(cv),])
        colnames(out) <- c("Gene", "CV", colnames(mat))
        dataset <- mat[names(cv),]
      }
      incProgress(0.5, detail = "Complete!")
    })

    dataset

  })

  pvca_uncorr_react <- reactive({

    req(pvca_mv_features_uncorr_matrix())
    meta <- aligned_meta_file()
    mat <- as.matrix(pvca_mv_features_uncorr_matrix())

    if (isTruthy(input$variable_choices_EV) & isTruthy(input$pvcaPct)) {
      vars <- input$variable_choices_EV
      if (length(vars) > 2) {
        withProgress(message = "Processing Uncorrected", value = 0, {
          incProgress(0.5, detail = "Running PCVA")
          pvca_res <- statVisual::PVCA(
            clin_data = meta,                # clinical
            clin_subjid = colnames(meta)[1], # sample name column
            gene_data = mat,                 # expression data
            batch.factors = vars,
            pct_threshold = input$pvcaPct)            # batch columns
          incProgress(0.5, detail = "Complete!")
        })
        pvca_res
      }
    }


  })

  output$uncorrected_PVCA_plot <- renderPlot({

    p <- pvca_uncorr_react()
    p

  })

  #### SVA ---------------------------------------------------------------------
  output$uncorrected_SVA_variable_of_interest <- renderUI({
    VarChoices <- c(unlist(strsplit(batch_names_from_meta()[-1], ",")))
    if (isTruthy(SVA_variable_of_interest_bc2())) {
      VarSelected <- SVA_variable_of_interest_bc2()
    } else {
      VarSelected <- VarChoices[which(!VarChoices %in% c(input$batch1_choices,input$batch2_choices))]
    }
    selectInput(
      "uncorrected_SVA_variable_of_interest",
      "Select Variable of Interest",
      VarChoices,
      selected = VarSelected
    )
  })
  uncorrected_SVA_variable_of_interest2 <- reactive({
    if (isTruthy(input$uncorrected_SVA_variable_of_interest)) {
      if (input$uncorrected_SVA_variable_of_interest == "Select"){
        NULL
      }else {
        uncorrected_SVA_variable_of_interest2 <- input$uncorrected_SVA_variable_of_interest
      }
    }
  })
  uncorrected_SVA_nsv <- reactive({
    if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      withProgress(message = "Processing Uncorrected", value = 0, {
        incProgress(0.5, detail = "Number of surrogate variables")
        uncorrected_mod <- model.matrix(reformulate(uncorrected_SVA_variable_of_interest2()), data = aligned_meta_file())
        uncorrected_mod_null <- model.matrix(~1, data = aligned_meta_file())
        uncorrected_SVA_nsv <- sva::num.sv(uncorrected_numeric_matrix()[,-1], uncorrected_mod, method=svaMethod, vfilter = svaVarNum)
        print(uncorrected_SVA_nsv)
        incProgress(0.5, detail = "Complete!")
      })
      uncorrected_SVA_nsv
    }
  })
  uncorrected_SVA_object <- reactive({
    if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      req(uncorrected_SVA_nsv())
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      withProgress(message = "Processing Uncorrected", value = 0, {
        incProgress(0.5, detail = "Running SVA")
        uncorrected_mod <- model.matrix(reformulate(uncorrected_SVA_variable_of_interest2()), data = aligned_meta_file())
        uncorrected_mod_null <- model.matrix(~1, data = aligned_meta_file())
        uncorrected_SVA_matrix <- as.matrix(uncorrected_numeric_matrix()[,-1])
        uncorrected_SVA_object <- sva::sva(uncorrected_SVA_matrix, uncorrected_mod, uncorrected_mod_null, n.sv = uncorrected_SVA_nsv(),
                                           numSVmethod = svaMethod, vfilter = svaVarNum)
        incProgress(0.5, detail = "Complete!")
      })
      df <- as.data.frame(uncorrected_SVA_object$sv)
      colnames(df) <- paste0("SVA_Uncorrected_Surrogate_Vars_",seq(ncol(df)))
      df <- cbind(aligned_meta_file(),df)
      uncorr_boxplot_meta_file(df)
      uncorrected_SVA_object
    }
  })
  observe({

    main_meta <- aligned_meta_file()
    uncorr_meta <- uncorr_boxplot_meta_file()
    corr_meta <- corr_boxplot_meta_file()
    main_meta_new <- merge(main_meta,uncorr_meta)
    main_meta_new <- merge(main_meta_new,corr_meta)
    boxplot_meta_file(main_meta_new)

  })
  uncorrected_SVA_probability_df <- reactive({
    if (isTruthy(uncorrected_SVA_object())) {
      uncorrected_SVA_probability_df <- data.frame(
        Genes = 1:length(uncorrected_SVA_object()$pprob.gam),
        latent_variable = uncorrected_SVA_object()$pprob.gam,
        variable_of_intrest = uncorrected_SVA_object()$pprob.b
      )
      uncorrected_SVA_probability_df_longer <- pivot_longer(uncorrected_SVA_probability_df, !Genes, names_to = "variable_type")
      uncorrected_SVA_probability_df_longer <- dplyr::rename(uncorrected_SVA_probability_df_longer, "probability_association_of_each_gene" = "value")
      uncorrected_SVA_probability_df_longer
    }
  })

  uncorrected_SVA_probability_ggplot <- reactive({
    if (isTruthy(uncorrected_SVA_probability_df())) {
      ggplot(uncorrected_SVA_probability_df(), aes(x = probability_association_of_each_gene,
                                                   fill = variable_type)) +
        geom_density(alpha = 0.5)
    }
  })

  output$uncorrected_SVA_probability_density <- renderPlot({
    uncorrected_SVA_probability_ggplot()
  })
  output$uncorrected_SVA_nsv_print <- renderPrint({
    if (isTruthy(uncorrected_SVA_nsv())) {
      print(paste("The Number of Estimated Surrogate Variables are:", uncorrected_SVA_nsv()))
    }
  })
  #### Box Plot ----------------------------------------------------------------

  output$rendBPsampSubset <- renderUI({

    FeatChoices <- c("Select All Samples",colnames(aligned_meta_file())[-1])
    selectInput("BPsampSubset","Subset Samples By:",choices = FeatChoices, selected = FeatChoices[1])

  })

  output$rendBPsampCriteria <- renderUI({

    req(input$BPsampSubset)
    subSelect <- input$BPsampSubset
    if (subSelect != "Select All Samples") {
      sampCrit <- unique(aligned_meta_file()[,subSelect])
      selectInput("BPsampCrit","Sample Criteria:",choices = sampCrit)
    }

  })

  output$rendBPgroupCriteria <- renderUI({

    GroupChoices <- colnames(aligned_meta_file())[-1]
    if (isTruthy(input$BPsampSubset)) {
      GroupChoices <- GroupChoices[which(GroupChoices!=input$BPsampSubset)]
    }
    selectInput("BPgroupCriteria","Grouping Criteria:",choices = GroupChoices, selected = GroupChoices[1])

  })

  output$rendBPgroupSelection <- renderUI({

    req(input$BPgroupCriteria)
    req(aligned_meta_file())
    meta <- aligned_meta_file()
    groupCrit <- input$BPgroupCriteria
    if (input$BPremoveSingles == T) {
      tab <- table(meta[,groupCrit])
      meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
    }
    GroupSelec <- unique(meta[,groupCrit])
    selectInput("BPgroupSelection","Select Groups:",choices = GroupSelec, selected = GroupSelec, multiple = T)

  })

  PCA_Proj_uncorr_Samples <- reactive({
    if (input$uncorrected_panel == "box") {
      withProgress(message = "Processing Uncorrected", value = 0, {
        incProgress(0.5, detail = "Running PCA")
        mat_uncorr <- uncorrected_numeric_matrix()
        rownames(mat_uncorr) <- mat_uncorr[,1]
        mat_uncorr <- mat_uncorr[,-1]
        mat_uncorr_t <- as.data.frame(t(mat_uncorr))
        pca_uncorr <- prcomp(t(mat_uncorr_t), scale = TRUE)
        rot_uncorr <- as.data.frame(t(pca_uncorr[["rotation"]]))
        rot_uncorr$`Principal Component` <- rownames(rot_uncorr)
        rot_uncorr <- rot_uncorr %>% relocate(`Principal Component`)
        incProgress(0.5, detail = "Complete!")
      })
      rot_uncorr
    }
  })

  output$rendImmuneDeconvMethods <- renderUI({

    if (input$BPFeatureCategory == "Immune Deconvolution Features") {
      if (input$HumanOrMouse == "Human") {
        selectInput("ImmuneDeconvMethods","Immune Deconvolution Methods:",
                    c("mcp_counter","estimate","quantiseq","xcell","epic","abis"))
      } else {
        selectInput("ImmuneDeconvMethods","Immune Deconvolution Methods:",
                    #c("mmcp_counter","seqimmucc","dcq","base"))
                    c("mmcp_counter","dcq","base"))
      }
    }

  })

  ImmDeconv_uncorr_react <- reactive({

    req(input$ImmuneDeconvMethods)
    deconvMethod <- input$ImmuneDeconvMethods
    mat <- uncorrected_numeric_matrix()
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]
    withProgress(message = "Processing Uncorrected", value = 0, {
      incProgress(0.5, detail = "Running Immune Deconvolution")
      if (input$HumanOrMouse == "Human") {
        deconv <- as.data.frame(deconvolute(mat, deconvMethod))
      } else {
        deconv <- as.data.frame(deconvolute_mouse(mat, deconvMethod))
      }
      incProgress(0.5, detail = "Complete!")
    })
    deconv

  })

  Mat_for_ssGSEA_uncorr <- reactive({

    mat <- uncorrected_numeric_matrix()
    if (input$HumanOrMouse == "Mouse") {
      mat_conv <- MouseToHuman(mat,MM_HS_Conv)
      mat_conv
    } else {
      mat <- uncorrected_numeric_matrix()
      mat
    }

  })

  UserGeneset_react <- reactive({

    if (input$BPGeneSetCat == "User Upload") {
      gs.u <- input$UserGeneset
      ext <- tools::file_ext(gs.u$datapath)
      req(gs.u)
      validate(need(ext == c("txt","tsv","csv","zip","RData","gmt"), "Please upload .txt, .tsv, .csv .zip, .RData, or .gmt file"))
      if (ext == "gmt") {
        gmt <- clusterProfiler::read.gmt(gs.u$datapath)
        gmt
      } else if (ext == "RData") {
        gs_u <- loadRData(gs.u$datapath)
        gs_u
      } else if (ext == "csv") {
        gmt <- as.data.frame(read_delim(gs.u$datapath, delim = ','))
        gmt
      } else {
        gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t'))
        gmt
      }
    }

  })

  GeneSetTableBack_react <- reactive({

    if (input$BPGeneSetCat == "User Upload") {
      req(input$UserGeneset)
      ext <- tools::file_ext(input$UserGeneset$datapath)
      if (ext == "RData") {
        gs_u <- UserGeneset_react()
        uGS_table <- as.data.frame(names(gs_u))
        colnames(uGS_table)[1] <- "GeneSet"
        uGS_table
      } else {
        gmt <- UserGeneset_react()
        uGS_table <- as.data.frame(unique(gmt[,1]))
        colnames(uGS_table)[1] <- "GeneSet"
        uGS_table
      }
    } else {
      GeneSetTable_sub <- geneset_df[which(geneset_df[,1] == input$BPGeneSetCat),]
      GeneSetTable_sub <- GeneSetTable_sub[,-c(1,2)]
      colnames(GeneSetTable_sub) <- c("Gene Set Category","Gene Set")
      GeneSetTable_sub
    }

  })

  GeneSetSelected <- reactive({

    gs_name <- as.character(GeneSetTableBack_react()[input$GeneSetTableUI_rows_selected,ncol(GeneSetTableBack_react())])
    if (input$BPGeneSetCat == "User Upload") {
      ext <- tools::file_ext(input$UserGeneset$datapath)
      gs <- UserGeneset_react()
      if (ext == "RData") {
        gs_u <- UserGeneset_react()
        gs <- gs_u[gs_name]
        gs
      } else {
        gs_u <- UserGeneset_react()
        gs_u_gs <- gs_u[which(gs_u[,1] == gs_name),]
        gs <- list()
        for (i in unique(gs_u_gs[,1])){
          gs[[i]] <- gs_u_gs[gs_u_gs[,1] == i,][,2]
        }
        gs
      }
    } else {
      gs <- geneset[gs_name]
      gs
    }

  })


  output$GeneSetTableUI <- DT::renderDataTable({

    GeneSetTable_sub <- GeneSetTableBack_react()
    DT::datatable(GeneSetTable_sub,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)

  })

  BP_Feature_Choices <- reactive({

    Features <- NULL
    FeatCat <- input$BPFeatureCategory
    if (FeatCat == "Matrix Features") {
      Features <- uncorrected_numeric_matrix()[,1]
      Features
    } else if (FeatCat == "Meta Features") {
      #Features <- colnames(aligned_meta_file())[-1]
      Features <- colnames(boxplot_meta_file())[-1]
      Features
    } else if (FeatCat == "PCA Projections") {
      Features <- PCA_Proj_uncorr_Samples()[,1]
      Features
    } else if (FeatCat == "Immune Deconvolution Features") {
      Features <- ImmDeconv_uncorr_react()[,1]
      Features
    } else {
      Features <- NULL
      Features
    }

  })

  output$rendBPlogOpt <- renderUI({

    if (input$BPFeatureCategory != "Gene Set Pathways") {
      selectInput("BPlogOpt","Log:", choices = c("No Log","Log2","Log2+1","Log10","Log10+1"))
    }

  })

  shiny::observe({

    updateSelectizeInput(session = session, inputId = "BPFeatSelection",
                         choices = BP_Feature_Choices(),
                         server = T)

  })

  output$renduncorrected_Box_plot <- renderUI({
    plotHeight <- input$BPplotHeight
    plotWidth <- input$BPplotWidth
    shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_Box_plot",height = plotHeight, width = plotWidth))
  })

  CohortBPPlot_df_react <- reactive({

    req(input$BPsampSubset)
    #meta <- aligned_meta_file()
    meta <- boxplot_meta_file()
    sampSubset <- input$BPsampSubset
    sampCrit <- input$BPsampCrit
    groupCrit <- input$BPgroupCriteria
    FeatCat <- input$BPFeatureCategory
    Feature <- input$BPFeatSelection
    NameCol <- colnames(meta)[1]
    removeSingles <- input$BPremoveSingles
    BPlog <- input$BPlogOpt
    if (FeatCat == "Matrix Features") {
      mat <- uncorrected_numeric_matrix()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
    } else if (FeatCat == "PCA Projections") {
      mat <- PCA_Proj_uncorr_Samples()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
    } else if (FeatCat == "Immune Deconvolution Features") {
      mat <- ImmDeconv_uncorr_react()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
    } else if (FeatCat == "Gene Set Pathways") {
      gs_name <- as.character(GeneSetTableBack_react()[input$GeneSetTableUI_rows_selected,ncol(GeneSetTableBack_react())])
      if (isTruthy(gs_name)) {
        gs <- geneset[gs_name]
        Feature <- gs_name
        #mat <- uncorrected_numeric_matrix()
        mat <- Mat_for_ssGSEA_uncorr()
        rownames(mat) <- mat[,1]
        mat <- mat[,-1]
        mat <- as.matrix(mat)
        withProgress(message = "Processing Uncorrected", value = 0, {
          incProgress(0.5, detail = "Running ssGSEA")
          ssGSEA_param <- GSVA::ssgseaParam(mat,gs)
          ssGSEA <- GSVA::gsva(ssGSEA_param)
          incProgress(0.5, detail = "Complete!")
        })
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA[,NameCol] <- rownames(ssGSEA)
        meta <- merge(meta,ssGSEA)
      }
    }

    if (isTruthy(Feature)) {
      if (Feature %in% colnames(meta)) {
        if (sampSubset != "Select All Samples") {
          meta <- meta[which(meta[,sampSubset] == sampCrit),]
          meta <- meta %>% select(any_of(c(NameCol,sampSubset,groupCrit,Feature)))
        } else {
          meta <- meta %>% select(any_of(c(NameCol,groupCrit,Feature)))
        }
        if (removeSingles == T) {
          tab <- table(meta[,groupCrit])
          meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
        }
        meta[,Feature] <- as.numeric(meta[,Feature])
        meta <- meta[which(!is.na(meta[,Feature])),]
        meta <- meta[which(meta[,Feature]!=Inf & meta[,Feature]!=-Inf),]
        meta[,groupCrit] <- as.factor(meta[,groupCrit])
        if (BPlog == "Log2") {
          meta[,Feature] <- log2(meta[,Feature])
        } else if (BPlog == "Log2+1") {
          meta[,Feature] <- log2(meta[,Feature]+1)
        } else if (BPlog == "Log10") {
          meta[,Feature] <- log10(meta[,Feature])
        } else if (BPlog == "Log10+1") {
          meta[,Feature] <- log10(meta[,Feature]+1)
        }

        meta
      }
      }


  })

  CohortBPPlot_react <- reactive({

    plotdf_full <- CohortBPPlot_df_react()
    sampSubset <- input$BPsampSubset
    sampCrit <- input$BPsampCrit
    groupCrit <- input$BPgroupCriteria
    FeatCat <- input$BPFeatureCategory
    Feature <- input$BPFeatSelection
    NameCol <- colnames(plotdf_full)[1]
    BPplottheme <- input$BPTheme
    StatMethod <- input$BPplotstatComp
    dotChoice <- input$BPplotsampledots
    dotSize <- input$BPplotDotSize
    bpFlip <- input$BPflipBP
    BPorVI <- input$BPorViolin
    if (FeatCat == "Gene Set Pathways") {
      gs_name <- as.character(GeneSetTableBack_react()[input$GeneSetTableUI_rows_selected,ncol(GeneSetTableBack_react())])
      Feature <- gs_name
    }
    Xaxis_font <- input$BPplot1XAxisSize              # Axis font size
    Yaxis_font <- input$BPplot1YAxisSize              # Axis font size
    Yaxis_lim <- input$BPplot1YAxisLim
    hjust_orient <- 1                                # Initial hjust
    axis_orient <- as.numeric(input$BPxAxisOrient)  # X-axis label orientation
    if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
      hjust_orient <- 0.5
    }
    BPorder <- input$BPplotXaxOrder
    BPGroupSelect <- input$BPgroupSelection

    if (isTruthy(Feature)) {
      plotdf <- plotdf_full[,c(NameCol,groupCrit,Feature)]
      plotdf <- plotdf[which(plotdf[,groupCrit] %in% BPGroupSelect),]
      colnames(plotdf) <- c("SampleName","Group","Feature")

      if (BPorder == "Descending"){
        barp <- ggplot(data = plotdf, aes(x=reorder(Group,-Feature, FUN = median),y=Feature, fill=Group))
        plotdf_dots <- plotdf
        plotdf_dots$Group <- reorder(plotdf_dots$Group,-plotdf_dots$Feature, FUN = median)
        plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
      }
      if (BPorder == "Ascending"){
        barp <- ggplot(data = plotdf, aes(x=reorder(Group,Feature, FUN = median),y=Feature, fill=Group))
        plotdf_dots <- plotdf
        plotdf_dots$Group <- reorder(plotdf_dots$Group,plotdf_dots$Feature, FUN = median)
        plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
      }
      if (BPorder == "Not Specificed"){
        barp <- ggplot(data = plotdf, aes(x=Group,y=Feature, fill=Group))
        plotdf_dots <- plotdf
        plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
      }
      if (BPorVI == "Box Plot") {
        barp <- barp + geom_boxplot(width = 0.5, lwd = 1)
      }
      if (BPorVI == "Violin Plot") {
        barp <- barp + geom_violin() +
          stat_summary(fun=median, geom="crossbar", width=0.5, color="black")
      }
      if (isTruthy(Yaxis_lim)) {
        barp <- barp +
          ylim(paste0(as.numeric(strsplit(Yaxis_lim,",")[[1]][1]),as.numeric(strsplit(Yaxis_lim,",")[[1]][2])))
      }

      barp <- barp +
        get(BPplottheme)() +
        labs(#title = BPTitle_in,
          x = groupCrit, y = Feature,
          fill = groupCrit)
      if (StatMethod != "none") {
        barp <- barp + ggpubr::stat_compare_means(method = StatMethod)
      }
      if (dotChoice) {
        barp <- barp + geom_point(data = plotdf_dots, aes(x=xj), col="grey14", size=dotSize)
      }
      barp <- barp + theme(axis.text.x = element_text(size = Xaxis_font,angle = axis_orient, hjust = hjust_orient),
                           axis.title.x = element_text(size = Xaxis_font),
                           axis.text.y = element_text(size = Yaxis_font),
                           axis.title.y = element_text(size = Yaxis_font),
                           legend.position = "none")
      if (bpFlip) {
        barp <- barp + coord_flip()
      }
      barp
    }

  })

  output$uncorrected_Box_plot <- renderPlot({

    barp <- CohortBPPlot_react()
    barp

  })

  ### END OF UNCORRECTED PLOT DATA ---------------------------------------------


  ### START OF BATCH CORRECTION CRITERIA AND PLOT DATA -------------------------


  #### Batch Criteria Selection ------------------------------------------------
  # Initial batch selection
  output$rendInit_Batch_Select <- shiny::renderUI({
    req(aligned_meta_file())
    batch_choices <- strsplit(batch_names_from_meta()[-1], ",")
    shiny::selectInput("Init_Batch_Select",
                       "Select Batch Variable",
                       c("Select", batch_choices),
                       selected = 2)

  })
  InitialBatchSelected <- reactive({
    if (isTruthy(input$Init_Batch_Select)) {
      if (input$Init_Batch_Select != "Select") {
        Batch <- input$Init_Batch_Select
        Batch
      } else {
        Batch <- NULL
        Batch
      }
    } else {
      Batch <- NULL
      Batch
    }
  })
  # Selection of batch criteria for correction.
  output$batch1_selection_limma <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c("Select", unlist(strsplit(batch_names_from_meta()[-1], ",")))[2]
    }
    shiny::selectInput("batch1_choices_limma",
                       "Batch 1 Variable",
                       c("Select", unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       selected = BatchSelected)
  })
  output$batch2_selection_limma <- shiny::renderUI({
    shiny::selectInput("batch2_choices_limma",
                       "Batch 2 Variable",
                       c("Select", unlist(strsplit(batch_names_from_meta()[-1], ","))))
  })
  output$covariate_choices_limma <- shiny::renderUI({
    shiny::selectInput("covariate_choices_limma",
                       "Select Any Covariates for Correction",
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       multiple = TRUE
    )
  })
  output$batch1_selection_ComBatseq <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(unlist(strsplit(batch_names_from_meta()[-1], ",")))[1]
    }
    shiny::selectInput("batch1_choices_ComBatseq",
                       "Select Variable for Batch",
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       selected = BatchSelected)
  })
  output$covariate_choices_ComBatseq <- shiny::renderUI({
    shiny::selectInput("covariate_choices_ComBatseq",
                       "Select Biological Variables",
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       multiple = TRUE
    )
  })
  output$batch_selection_ComBat <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(unlist(strsplit(batch_names_from_meta()[-1], ",")))[1]
    }
    shiny::selectInput(
      "batch_selection_ComBat",
      "Select Batch",
      c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
      selected = BatchSelected
    )
  })
  output$batch_selection_mean_centering <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(unlist(strsplit(batch_names_from_meta()[-1], ",")))[1]
    }
    shiny::selectInput(
      "batch1_choices_mean_centering",
      "Select Batch for Mean Centering",
      c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
      selected = BatchSelected
    )
  })
  output$batch_selection_harman <- renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(unlist(strsplit(batch_names_from_meta()[-1], ",")))[1]
    }
    selectInput(
      "batch1_choices_harman",
      "Select Batch",
      c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
      selected = BatchSelected
    )
  })
  output$treatment_selection_harman <- renderUI({
    treatOpt <- c(unlist(strsplit(batch_names_from_meta()[-1], ",")))
    harmBatch <- input$batch1_choices_harman
    selectInput(
      "treatment_selection_harman",
      "Variable of Interest",
      c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
      treatOpt[which(treatOpt != harmBatch)][1]
    )
  })
  output$SVA_variable_of_interest_bc <- renderUI({
    VarChoices <- c(unlist(strsplit(batch_names_from_meta()[-1], ",")))
    VarSelect <- VarChoices[which(VarChoices != InitialBatchSelected())]
    selectInput(
      "SVA_variable_of_interest_bc",
      "Variable of Interest",
      VarChoices,
      selected = VarSelect[1]
    )
  })


  # Allowing select input to provide NULL input into the functions
  batch1_choices <- shiny::reactive({
    if (input$batch_correction_method == "Limma"){
      if (isTruthy(input$batch1_choices_limma)) {
        if (input$batch1_choices_limma == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_limma
        }
      }
    }else if (input$batch_correction_method == "ComBat"){
      if (isTruthy(input$batch_selection_ComBat)) {
        if (input$batch_selection_ComBat == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch_selection_ComBat
        }
      }
    }else if (input$batch_correction_method == "ComBatseq"){
      if (isTruthy(input$batch1_choices_ComBatseq)) {
        if (input$batch1_choices_ComBatseq == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_ComBatseq
        }
      }
    }else if (input$batch_correction_method == "Mean Centering"){
      if (isTruthy(input$batch1_choices_mean_centering)) {
        if (input$batch1_choices_mean_centering == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_mean_centering
        }
      }
    }else if (input$batch_correction_method == "Harman"){
      if (isTruthy(input$batch1_choices_harman)) {
        if (input$batch1_choices_harman == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_harman
        }
      }
    }
  })
  batch2_choices <- shiny::reactive({
    if (isTruthy(input$batch2_choices_limma)) {
      if (input$batch2_choices_limma == "Select"){
        batch2_choices <- NULL
      } else{
        batch2_choices <- input$batch2_choices_limma
      }
    }
  })
  harman_treatment_input <- reactive({
    if (isTruthy(input$treatment_selection_harman)) {
      if (input$treatment_selection_harman == "Select"){
        NULL
      }else {
        harman_treatment_input <- input$treatment_selection_harman
      }
    }
  })
  SVA_variable_of_interest_bc2 <- reactive({
    if (isTruthy(input$SVA_variable_of_interest_bc)) {
      if (input$SVA_variable_of_interest_bc == "Select"){
        NULL
      }else {
        SVA_variable_of_interest_bc2 <- input$SVA_variable_of_interest_bc
      }
    }
  })

  #### Model Matrix ------------------------------------------------
  # creation of model matrix from the meta data (used in combatseq and limma)
  model_matrix <- shiny::reactive({
    if (input$batch_correction_method == "Limma"){
      if (is.null(input$covariate_choices_limma)){
        model_matrix <- NULL
      } else {
        total_covariates <- paste0(input$covariate_choices_limma,collapse = "+")
        model_matrix <- stats::model.matrix(reformulate(total_covariates), data = as.data.frame(aligned_meta_file()))
      }
    } else if (input$batch_correction_method == "ComBatseq"){
      if (is.null(input$covariate_choices_ComBatseq)){
        model_matrix <- NULL
      } else {
        counter <- 0
        for (covariate in 1:length(input$covariate_choices_ComBatseq)){
          variable_object <- input$covariate_choices_ComBatseq[covariate]
          if (counter == 0){
            total_covariates <- paste(variable_object, sep = "")
            counter <- counter + 1
          } else {
            total_covariates <- paste(total_covariates,"+", variable_object, sep = "")
            counter <- counter + 1
          }
        }
        model_matrix <- stats::model.matrix(reformulate(total_covariates), data = as.data.frame(aligned_meta_file()))
      }
    }
  })

  #RUVg housekeeping genes
  RUVg_housekeeping <- reactive({
    if (input$HumanOrMouse == "Human") {
      if(input$RUVg_housekeeping_selection == "Eisenberg"){
        RUVg_housekeeping <- Eisenberg_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection == "Lin500"){
        RUVg_housekeeping <- Lin500_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection == "HSIAO"){
        RUVg_housekeeping <- HSIAO_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection == "UserInput"){
        RUVg_housekeeping_input <- read.delim(input$RUVg_user_control_genes$datapath, sep = '\t')
        RUVg_housekeeping <- RUVg_housekeeping_input[,1]
        RUVg_housekeeping
      }
    } else {
      if(input$RUVg_housekeeping_selection == "Eisenberg"){
        RUVg_housekeeping <- Eisenberg_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection == "Lin500"){
        RUVg_housekeeping <- Lin500_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection == "HSIAO"){
        RUVg_housekeeping <- HSIAO_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection == "UserInput"){
        RUVg_housekeeping_input <- read.delim(input$RUVg_user_control_genes$datapath, sep = '\t')
        RUVg_housekeeping <- RUVg_housekeeping_input[,1]
        RUVg_housekeeping
      }
    }

  })

  Log_Norm_Matrix <- reactive({

    set.seed(input$SeedSet)
    req(uncorrected_matrix_filtered())
    uncorrected_numeric_matrix <- uncorrected_matrix_filtered()
    rownames(uncorrected_numeric_matrix) <- uncorrected_numeric_matrix[,1]
    uncorrected_numeric_matrix <- uncorrected_numeric_matrix[,-1]
    if (input$RawCountCheck) {
      if (input$RawCountNorm == "none") {
        updateCheckboxInput(session,"Log_Choice", value = F)
      } else {
        mat <- uncorrected_numeric_matrix
        mat_dgeList <- DGEList(counts = as.matrix(mat))
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Normalizing Factors")
          mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
          incProgress(0.25, detail = "CPM")
          uncorrected_numeric_matrix <- edgeR::cpm(mat_dgeList_Norm)
          incProgress(0.75, detail = "Complete!")
        })
      }
    }
    if (input$Log_Choice){
      withProgress(message = "Processing", value = 0, {
        incProgress(0.5, detail = "Logging Matrix")
        uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix + 1)
        incProgress(0.5, detail = "Complete!")
      })
    }
    if (isTruthy(input$QuantNorm)) {
      if (input$QuantNorm){
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Quantile Normalization")
          uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
          uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
          incProgress(0.5, detail = "Complete!")
        })
      }
    }
    uncorrected_numeric_matrix

  })

  #### Batch Correction ------------------------------------------------
  # The actual batch correction
  batch_correction <- shiny::reactive({

    uncorrected_numeric_matrix <- Log_Norm_Matrix()

    if(input$batch_correction_method == "Limma"){
      if (isTruthy(input$batch1_choices_limma) & isTruthy(input$batch2_choices_limma)) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Running Limma batch correction")
          batch_correction <- limma::removeBatchEffect(
            uncorrected_numeric_matrix,
            batch = c(unlist(aligned_meta_file()[,batch1_choices()])),
            batch2 = c(unlist(aligned_meta_file()[,batch2_choices()])),
            covariates = model_matrix()
          )
          incProgress(0.5, detail = "Complete!")
        })
        batch_correction
      }
    }else if(input$batch_correction_method == "ComBat"){
      if (isTruthy(batch1_choices())) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Running ComBat batch correction")
          batch_combat <- c(unlist(aligned_meta_file()[,batch1_choices()]))
          incProgress(0.25, detail = "Running ComBat batch correction")
          modcombat <-  stats::model.matrix(~1, data = as.data.frame(aligned_meta_file()))
          incProgress(0.25, detail = "Running ComBat batch correction")
          batch_correction <- sva::ComBat(
            dat = uncorrected_numeric_matrix,
            batch = batch_combat,
            mod = modcombat,
            par.prior = input$combat_parametric
          )
          incProgress(0.25, detail = "Complete!")
        })
        batch_correction
      }
    }else if(input$batch_correction_method == "Mean Centering"){
      if (isTruthy(batch1_choices())) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Running Mean Centering Batch Correction")
          mean_centering_batch = c(unlist(aligned_meta_file()[,batch1_choices()]))
          mean_centering_data = t(uncorrected_numeric_matrix)
          incProgress(0.25, detail = "Running Mean Centering Batch Correction")
          mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
          mean_center_correction <- as.data.frame(t(mean_center$xadj))
          batch_correction <- mean_center_correction
          incProgress(0.5, detail = "Complete!")
        })
        batch_correction
      }
    }else if(input$batch_correction_method == "ComBatseq"){
      if (isTruthy(batch1_choices())) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Running ComBatseq Batch Correction")
          combatseq_corrected <- sva::ComBat_seq(
            as.matrix((2^uncorrected_numeric_matrix)),
            batch = c(unlist(aligned_meta_file()[,batch1_choices()])),
            covar_mod = model_matrix()
          )
          incProgress(0.25, detail = "Running ComBatseq Batch Correction")
          log2(combatseq_corrected + 1)
          batch_correction <- combatseq_corrected
          incProgress(0.5, detail = "Complete!")
        })
        batch_correction
      }
    }else if(input$batch_correction_method == "Harman"){
      req(harman_treatment_input())
      req(batch1_choices())
      withProgress(message = "Processing", value = 0, {
        incProgress(0.25, detail = "Running Harman Batch Correction")
        harman_correction_PCA <- Harman::harman(
          uncorrected_numeric_matrix,
          expt = aligned_meta_file()[,harman_treatment_input()],
          batch = aligned_meta_file()[,batch1_choices()],
          limit = 0.95,
          printInfo = T,
          randseed = input$SeedSet
        )
        incProgress(0.25, detail = "Reconstructing Corrected Data")
        harman_correction <- Harman::reconstructData(harman_correction_PCA)
        batch_correction <- harman_correction
        incProgress(0.5, detail = "Complete!")
      })
      batch_correction
    }else if(input$batch_correction_method == "RUVg"){
      if(input$RUVg_housekeeping_selection == "UserInput"){
        req(input$RUVg_user_control_genes)
      }
      if (any(RUVg_housekeeping() %in% rownames(uncorrected_numeric_matrix))) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Running RUVg Batch Correction")
          RUVg_housekeeping_genes <- RUVg_housekeeping()[which(RUVg_housekeeping() %in% rownames(uncorrected_numeric_matrix))]
          RUVg_matrix <- uncorrected_numeric_matrix
          RUVg_correction <- RUVSeq::RUVg(
            as.matrix(RUVg_matrix),
            cIdx = RUVg_housekeeping_genes,
            k = input$RUVg_estimate_factors,
            drop = input$RUVg_drop_factors,
            center = input$RUVg_mean_centered,
            round = input$RUVg_rounded,
            tolerance = input$RUVg_tolerance,
            isLog = T
          )
          RUVg_correction_matrix <- as.data.frame(RUVg_correction$normalizedCounts)
          batch_correction <- RUVg_correction_matrix
          incProgress(0.5, detail = "Complete!")
        })
        batch_correction
      }
    }else if(input$batch_correction_method == "SVA"){
      if (isTruthy(SVA_variable_of_interest_bc2())) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Number of surrogate variables")
          expression_data <-  as.matrix(uncorrected_numeric_matrix()[,-1])
          mod <-  model.matrix(reformulate(SVA_variable_of_interest_bc2()), data = aligned_meta_file())
          mod0 <- model.matrix(~1, data = aligned_meta_file())
          n.sv <- sva::num.sv(expression_data, mod, method = input$svaMethod_bc)
          incProgress(0.25, detail = "Running SVA")
          svobj <- sva::sva(expression_data, mod , mod0, n.sv = n.sv)
          incProgress(0.25, detail = "Running Frozen SVA for batch correction")
          fsvaobj <- sva::fsva(expression_data, mod, svobj, expression_data)
          batch_correction <- fsvaobj$db
          incProgress(0.25, detail = "Complete!")
        })
        batch_correction
      }
    }else if(input$batch_correction_method == "Select"){
      uncorrected_numeric_matrix()[,-1]
    }
  })

  # Test text for troubleshooting. Will delete in final product.
  output$test_print <- shiny::renderText({
    print(batch1_choices())
  })

  # Putting the data frame back together for downstream analysis
  corrected_numeric_matrix <- shiny::reactive({
    if (isTruthy(batch_correction())) {
      corrected_numeric_matrix <- as.data.frame(batch_correction())
      row.names(corrected_numeric_matrix) <- uncorrected_numeric_matrix()[,1]
      corrected_numeric_matrix
    }
  })
  corrected_numeric_matrix2 <- shiny::reactive({
    if (isTruthy(batch_correction())) {
      corrected_numeric_matrix2 <- as.data.frame(batch_correction())
      corrected_numeric_matrix2$Genes <- uncorrected_numeric_matrix()[,1]
      corrected_numeric_matrix2 <- relocate(corrected_numeric_matrix2, Genes)
    }
  })

  # Rendering the corrected matrix
  output$corrected_matrix <- DT::renderDataTable({
    if (isTruthy(input$batch_correction_method)) {
      if (input$batch_correction_method != "Select") {
        DT::datatable(corrected_numeric_matrix2(),
                      options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                     pageLength = 10,
                                     scrollX = T),
                      rownames = F)
      }
    }
  })
  output$rendsave_corrected_matrix <- renderUI({
    if (isTruthy(input$batch_correction_method)) {
      if (input$batch_correction_method != "Select") {
        shiny::fluidRow(
          column(4, style = 'padding-right:0px;',
                 shiny::actionButton("save_corrected_matrix", "Add to Zip File Export")
          ),
          column(4, style = 'padding-left:0px;',
                 shiny::downloadButton("dnldsave_corrected_matrix","Dowload Single Table")
          )
        )

      }
    }
  })

  #### PCA ------------------------------------------------
  # Generating the corrected PCA
  corrected_batch_choices_PCA2 <- shiny::reactive({
    if (!is.null(input$batch_choices_PCA)) {
      if (input$batch_choices_PCA != "Select"){
        corrected_batch_choices_PCA2 <- input$batch_choices_PCA
        corrected_batch_choices_PCA2
      }else {
        corrected_batch_choices_PCA2 <- NULL
        corrected_batch_choices_PCA2
      }
    }

  })

  PCA_data_corrected <- shiny::reactive({
    PCA_data_corrected <- cbind((t(batch_correction())), aligned_meta_file())
    PCA_data_corrected
  })

  corrected_PCA_react <- shiny::reactive({

    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "PCA Object")
      if (input$PCA_type == "Cluster Annotation"){
        pca <- cluster::pam(as.data.frame(t(batch_correction())),input$cluster_number)
      } else if (input$PCA_type == "Meta Annotation"){
        pca <- stats::prcomp(as.data.frame(t(batch_correction())), scale. = TRUE)
      }
      incProgress(0.5, detail = "Complete!")
    })
    pca

  })

  corrected_PCA_plot_df <- reactive({

    meta <- aligned_meta_file()
    pca <- corrected_PCA_react()
    batch_choice <- input$batch_choices_PCA
    hover_choice <- input$PCAhover
    metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice))
    metaCols <- metaCols[which(metaCols %in% colnames(meta))]
    PC_x <- "PC1"
    PC_y <- "PC2"
    corrected_PCA_points <- as.data.frame(pca$x[,c(1,2)])
    corrected_PCA_plot_df <- merge(corrected_PCA_points,meta[,metaCols,drop = F],by.x = 0, by.y = colnames(meta)[1])
    colnames(corrected_PCA_plot_df)[1] <- colnames(meta)[1]
    PCX <- gsub("^PC","",PC_x)
    PCY <- gsub("^PC","",PC_y)
    corrected_PCA_plot_df[,PC_x] <- as.numeric(corrected_PCA_plot_df[,PC_x]) / (pca$sdev[as.numeric(PCX)] * sqrt(nrow(corrected_PCA_plot_df)))
    corrected_PCA_plot_df[,PC_y] <- as.numeric(corrected_PCA_plot_df[,PC_y]) / (pca$sdev[as.numeric(PCY)] * sqrt(nrow(corrected_PCA_plot_df)))
    corrected_PCA_plot_df$text <- NA
    for (r in seq(nrow(corrected_PCA_plot_df))) {
      text_vec <- c()
      for (c in metaCols) {
        text_vec <- c(text_vec,
                      paste0("\\<br\\>\\<b\\>", c, ":\\</b\\> ",corrected_PCA_plot_df[r,c]))
      }
      text = paste(text_vec,collapse = '')
      text = gsub("\\\\", '', text)
      text = paste0(text,"<extra></extra>") #removes outside of box text
      corrected_PCA_plot_df$text[r] = text
    }
    corrected_PCA_plot_df

  })

  corrected_PCA_plot_forDlnd_react <- shiny::reactive({
    if (input$PCA_type == "Cluster Annotation"){
      ggplot2::autoplot(
        corrected_PCA_react(),
        frame = TRUE,
        frame.type = 'norm')
    } else if (input$PCA_type == "Meta Annotation"){
      ggplot2::autoplot(
        corrected_PCA_react(),
        data = PCA_data_corrected(),
        color = corrected_batch_choices_PCA2()
      )
    }
  })

  corrected_PCA_plot_react <- shiny::reactive({

    meta <- aligned_meta_file()
    plot_df <- corrected_PCA_plot_df()
    batch_choice <- input$batch_choices_PCA
    hover_choice <- input$PCAhover
    metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice))
    PC_x <- "PC1"
    PC_y <- "PC2"
    colnames(plot_df)[which(colnames(plot_df) == PC_x)] <- "X"
    colnames(plot_df)[which(colnames(plot_df) == PC_y)] <- "Y"
    dotSize <- input$PCAdotSize
    dotSize <- input$PCAdotSize
    fontSize <- input$pcaFontSize

    if (isTruthy(batch_choice)) {
      if (batch_choice != "Select") {
        colnames(plot_df)[which(colnames(plot_df) == batch_choice)] <- "ColorBatch"
        p4 <- plot_ly() %>%
          add_trace(
            data = plot_df,
            x = ~X,
            y = ~Y,
            type = "scatter",
            mode = "markers",
            color = ~ColorBatch,
            legendgroup=batch_choice,
            marker=list(size=dotSize,symbol = 'circle'),
            hovertemplate = ~text
          ) %>%
          layout(xaxis = list(zeroline = FALSE,title = PC_x),
                 yaxis = list(zeroline = FALSE,title = PC_y),
                 font = list(size = fontSize)) %>%
          config(
            toImageButtonOptions = list(
              format = "svg"
            )
          )
        p4
      } else {
        p4 <- plot_ly() %>%
          add_trace(
            data = plot_df,
            x = ~X,
            y = ~Y,
            type = "scatter",
            mode = "markers",
            marker = list(color = "black"),
            marker=list(size=dotSize,symbol = 'circle'),
            hovertemplate = ~text
          ) %>%
          layout(xaxis = list(zeroline = FALSE,title = PC_x),
                 yaxis = list(zeroline = FALSE,title = PC_y),
                 font = list(size = fontSize)) %>%
          config(
            toImageButtonOptions = list(
              format = "svg"
            )
          )
        p4
      }
    }

  })
  output$corrected_PCA <- renderPlotly({
    corrected_PCA_plot_react()
  })

  #### PCA MC ------------------------------------------------
  corrected_PCA_mc_color_choice2 <- shiny::reactive({
    if (isTruthy(input$PCA_mc_color_choice)) {
      if (input$PCA_mc_color_choice == "Select"){
        NULL
      }else {
        corrected_PCA_mc_color_choice2 <- input$PCA_mc_color_choice
      }
    }
  })
  corrected_PCA_multiple_components <- shiny::reactive({
    mat_PCA_mc <- batch_correction()
    names(mat_PCA_mc) <- NULL
    corrected_PCA_mc_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_PCA_mc)),
      colData = aligned_meta_file(),
      rowData = corrected_numeric_matrix2()[,1]
    )
    SummarizedExperiment::assay(corrected_PCA_mc_SCE, "logcounts") <- SingleCellExperiment::counts(corrected_PCA_mc_SCE)
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.25, detail = "Running multiple components PCA")
      corrected_PCA_mc_SCE <- scater::runPCA(corrected_PCA_mc_SCE, ncomponents = 50)
      incProgress(0.25, detail = "Plotting multiple components PCA")
      p <- scater::plotPCA(
        corrected_PCA_mc_SCE,
        ncomponents = input$PCA_mc_slider,
        colour_by = corrected_PCA_mc_color_choice2()
      )
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$corrected_PCA_multiple_components <- shiny::renderPlot({
    corrected_PCA_multiple_components()
  })

  #### PCA Details ------------------------------------------------
  #Corrected PCA Details plots
  corrected_PCA_details <- shiny::reactive({
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Running FactoMineR PCA")
      pc_obj <- FactoMineR::PCA(as.data.frame(t(batch_correction())), graph = F)
      incProgress(0.5, detail = "Complete!")
    })
    pc_obj
  })
  corrected_PCA_details2 <- shiny::reactive({
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Running stats prcomp PCA")
      pc_obj <- stats::prcomp(as.data.frame(t(batch_correction())))
      incProgress(0.5, detail = "Complete!")
    })
    pc_obj
  })
  corrected_scree_plot_react <- shiny::reactive({
    corrected_scree_eig <- as.data.frame(get_eig(corrected_PCA_details()))
    corrected_scree_eig$Dimensions <- gsub("Dim\\.","",rownames(corrected_scree_eig))
    corrected_scree_eig$`Variance Percent` <- paste0(round(corrected_scree_eig$variance.percent,1),"%")
    corrected_scree_eig_top10 <- corrected_scree_eig[1:10,]
    AxisTickFont <- input$pcaDetAxisTkSize
    AxisTitleFont <- input$pcaDetAxisTtSize
    LabelFont <- input$pcaDetLabelSize

    p <- ggplot(data=corrected_scree_eig_top10, aes(x=reorder(Dimensions,-variance.percent), y=variance.percent)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_minimal() +
      geom_text(aes(label=`Variance Percent`), vjust=-0.3,size=LabelFont) +
      labs(x = "Dimensions",y = "Variance Percent") +
      theme(axis.text = element_text(size = AxisTickFont),
            axis.title = element_text(size = AxisTitleFont))
    p
  })
  output$corrected_scree_plot <- shiny::renderPlot({
    corrected_scree_plot_react()
  })

  corrected_PCA_factors_choices2 <- shiny::reactive({
    if (isTruthy(input$PCA_factors_choices)) {
      if (input$PCA_factors_choices == "Select"){
        NULL
      }else {
        corrected_PCA_factors_choices2 <- input$PCA_factors_choices
      }
    }
  })
  corrected_PCA_individuals <- shiny::reactive({
    corrected_PCA_individuals <- corrected_PCA_details()
    corrected_PCA_factors <- as.data.frame(corrected_PCA_individuals$ind$contrib)
    corrected_PCA_factors$factors <- as.vector(aligned_meta_file()[,input$PCA_factors_choices])
    corrected_PCA_factors_final_sum <- stats::aggregate(. ~ factors, corrected_PCA_factors, sum)
    corrected_PCA_factors_final_mean <- stats::aggregate(. ~ factors, corrected_PCA_factors, mean)
    corrected_PCA_factors_final_sdv <- stats::aggregate(. ~ factors, corrected_PCA_factors, sd)
    corrected_PCA_factors_final_longer_sum <- tidyr::pivot_longer(corrected_PCA_factors_final_sum, !factors, names_to = "PC_Components")
    corrected_PCA_factors_final_longer_mean <- tidyr::pivot_longer(corrected_PCA_factors_final_mean, !factors, names_to = "PC_Components")
    corrected_PCA_factors_final_longer_sdv <- tidyr::pivot_longer(corrected_PCA_factors_final_sdv, !factors, names_to = "PC_Components")
    colnames(corrected_PCA_factors_final_longer_sum)[3] <- "Sum"
    colnames(corrected_PCA_factors_final_longer_mean)[3] <- "Mean"
    colnames(corrected_PCA_factors_final_longer_sdv)[3] <- "SDV"
    corrected_PCA_factors_final_longer <- corrected_PCA_factors_final_longer_sum[,1:2]
    corrected_PCA_factors_final_longer <- cbind(
      corrected_PCA_factors_final_longer,
      corrected_PCA_factors_final_longer_sum[,3],
      corrected_PCA_factors_final_longer_mean[,3],
      corrected_PCA_factors_final_longer_sdv[,3]
    )
    colnames(corrected_PCA_factors_final_longer)[1] <- input$PCA_factors_choices
    corrected_PCA_factors_final_longer
  })
  output$corrected_contribution_table <- DT::renderDataTable({
    if (isTruthy(input$PCA_factors_choices)) {
      DT::datatable(corrected_PCA_individuals(),
                    options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                   pageLength = 10,
                                   scrollX = T),
                    rownames = F) %>%
        formatRound(columns = c(3:ncol(corrected_PCA_individuals())), digits = 4)
    }

  })
  corrected_contribution_counts <- shiny::reactive({
    corrected_PCA_individuals <- corrected_PCA_details()
    corrected_PCA_factors <- as.data.frame(corrected_PCA_individuals$ind$contrib)
    corrected_PCA_factors$factors <- as.vector(aligned_meta_file()[,input$PCA_factors_choices])
    corrected_PCA_factors_final_count <- corrected_PCA_factors %>% dplyr::group_by(factors) %>% dplyr::summarise(individuals = length(factors))
    colnames(corrected_PCA_factors_final_count)[1] <- input$PCA_factors_choices
    corrected_PCA_factors_final_count
  })
  output$corrected_contribution_counts <- DT::renderDataTable({
    if (isTruthy(input$PCA_factors_choices)) {
      DT::datatable(corrected_contribution_counts(),
                    options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                   pageLength = 10,
                                   scrollX = T),
                    rownames = F)
    }

  })

  #### UMAP --------------------------------------------------------------------
  UMAP_ImmDeconv_corr_react <- reactive({

    req(input$UMAPImmuneDeconvMethods)
    deconvMethod <- input$UMAPImmuneDeconvMethods
    mat <- corrected_numeric_matrix2()
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]
    if (input$HumanOrMouse == "Human") {
      deconv <- as.data.frame(deconvolute(mat, deconvMethod))
      deconv
    } else {
      deconv <- as.data.frame(deconvolute_mouse(mat, deconvMethod))
      deconv
    }

  })

  UMAP_PCA_Proj_corr_Samples <- reactive({
    if (input$uncorrected_panel == "pca_main") {
      if (input$PCA_main_pan == "umap") {
        withProgress(message = "Processing Batch Corrected", value = 0, {
          incProgress(0.5, detail = "Running PCA")
          mat_uncorr <- corrected_numeric_matrix2()
          rownames(mat_uncorr) <- mat_uncorr[,1]
          mat_uncorr <- mat_uncorr[,-1]
          mat_uncorr_t <- as.data.frame(t(mat_uncorr))
          pca_uncorr <- prcomp(t(mat_uncorr_t), scale = TRUE)
          rot_uncorr <- as.data.frame(t(pca_uncorr[["rotation"]]))
          rot_uncorr$`Principal Component` <- rownames(rot_uncorr)
          rot_uncorr <- rot_uncorr %>% relocate(`Principal Component`)
          incProgress(0.5, detail = "Complete!")
        })
        rot_uncorr
      }
    }
  })

  UMAP_corr_anno_df <- reactive({

    FeatCat <- input$UMAPFeatureCategory
    Feature <- input$UMAPFeatSelection
    meta <- aligned_meta_file()
    NameCol <- colnames(meta)[1]
    if (FeatCat == "Matrix Features") {
      mat <- corrected_numeric_matrix2()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
      meta
    } else if (FeatCat == "PCA Projections") {
      mat <- UMAP_PCA_Proj_corr_Samples()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
      meta
    } else if (FeatCat == "Immune Deconvolution Features") {
      mat <- UMAP_ImmDeconv_corr_react()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
      meta
    } else if (FeatCat == "Gene Set Pathways") {
      gs_name <- GeneSetTableBack_react()[input$UMAPGeneSetTableUI_rows_selected,ncol(GeneSetTableBack_react())]
      gs <- geneset[gs_name]
      mat <- corrected_numeric_matrix2()
      rownames(mat) <- mat[,1]
      mat <- mat[,-1]
      mat <- as.matrix(mat)
      withProgress(message = "Processing Batch Corrected", value = 0, {
        incProgress(0.5, detail = "Running ssGSEA")
        ssGSEA_param <- GSVA::ssgseaParam(mat,gs)
        ssGSEA <- GSVA::gsva(ssGSEA_param)
        incProgress(0.5, detail = "Complete!")
      })
      ssGSEA <- as.data.frame(t(ssGSEA))
      ssGSEA[,NameCol] <- rownames(ssGSEA)
      meta <- merge(meta,ssGSEA)
      meta
    } else {
      meta <- meta
      meta
    }

  })

  corrected_umap_coord <- reactive({

    corrected_PCA_react_mat <- as.matrix(corrected_PCA_react()$x)
    umapNN <- input$UMAPnnb
    umapMinDist <- input$UMAPminDist
    umapMetric <- input$UMAPmetricSelec
    if (umapNN >= nrow(corrected_PCA_react_mat)) {
      umapNN <- nrow(corrected_PCA_react_mat) - 1
    }
    if (all(isTruthy(c(umapMetric,umapMinDist,umapNN)))) {

      otherMets <- c("euclidean","manhattan","cosine","pearson","pearson2")

      if (umapMetric == "hamming" || umapMetric == "correlation") {
        tdata_fit_df <- as.data.frame(uwot::umap(corrected_PCA_react_mat,metric = umapMetric, n_neighbors = umapNN, min_dist = umapMinDist))
      } else if (umapMetric %in% otherMets) {
        umap.mod$metric <- umapMetric
        umap.mod$min_dist <- umapMinDist
        umap.mod$n_neighbors <- umapNN
        withProgress(message = "Processing Batch Corrected", value = 0, {
          incProgress(0.5, detail = "Running UMAP")
          tdata_fit <- umap::umap(corrected_PCA_react_mat,config = umap.mod)
          incProgress(0.5, detail = "Complete!")
        })
        tdata_fit_df <- as.data.frame(tdata_fit$layout)
      }
      colnames(tdata_fit_df) <- c("UMAP1","UMAP2")
      tdata_fit_df <- tdata_fit_df %>%
        mutate(ID=row_number())
      tdata_fit_df$SampleName <- rownames(tdata_fit_df)
      tdata_fit_df <- tdata_fit_df %>%
        relocate(SampleName)
      tdata_fit_df
    }

  })

  corrected_UMAP_react <- reactive({

    req(input$UMAPFeatSelection)
    plot_df <- corrected_umap_coord()
    rownames(plot_df) <- plot_df[,1]
    UMAPdotSize <- 2
    umap_annoCol <- input$UMAPFeatSelection
    meta <- UMAP_corr_anno_df()
    NameCol <- colnames(meta)[1]
    umapdotSize <- input$umapdotSize
    umapAxisTkSize <- input$umapAxisTkSize
    umapAxisTtSize <- input$umapAxisTtSize

    plot_df <- merge(plot_df,meta[,c(NameCol,umap_annoCol)], by.x = colnames(plot_df)[1], by.y = NameCol)

    if (is.null(umap_annoCol)) {
      k <- plot_df %>%
        ggplot(aes(UMAP1, UMAP2,
                   text = paste("</br><b>Sample Name:</b> ", SampleName,
                                sep = "")))
    } else {
      k <- plot_df %>%
        ggplot(aes(UMAP1, UMAP2, color=!!sym(umap_annoCol),
                   text = paste("</br><b>Sample Name:</b> ", SampleName,
                                "</br><b>",umap_annoCol,":</b> ", !!sym(umap_annoCol),
                                sep = "")))
    }

    k <- k + geom_point(shape = 19, size = umapdotSize) +
      theme_minimal()


    k <- k + theme(axis.text.x = element_text(size = umapAxisTkSize),
                   axis.title.x = element_text(size = umapAxisTtSize),
                   axis.text.y = element_text(size = umapAxisTkSize),
                   axis.title.y = element_text(size = umapAxisTtSize))

    if (is.numeric(plot_df[,umap_annoCol])) {
      k <- k + scale_colour_gradient(low = "#56B1F7",high = "#132B43")
    }
    k

  })

  output$corrected_UMAP <- renderPlotly({

    p <- corrected_UMAP_react()
    ggplotly(p, tooltip = "text")

  })

  #### Cluster -----------------------------------------------------------------
  cluster_mv_features_corr_matrix <- reactive({
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$cluster_n_MV_features
      var_type <- input$VarianceMeasure
      mat <- corrected_numeric_matrix2()
      rownames(mat) <- mat[,1]
      mat <- mat[,-1]

      mad <- NULL
      var <- NULL
      cv <- NULL
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(log2(mat + 1), 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(log2(mat + 1), 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(log2(mat + 1), 1, cv))
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = topN)
        out <- cbind(names(cv), cv[names(cv)], mat[names(cv),])
        colnames(out) <- c("Gene", "CV", colnames(mat))
        dataset <- mat[names(cv),]
      }
      incProgress(0.5, detail = "Complete!")
    })

    dataset

  })
  corrected_elbow_analysis <- reactive({
    req(cluster_mv_features_corr_matrix())
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Running Elbow Analysis")
      corrected_elbow_analysis <- fviz_nbclust(x = t(cluster_mv_features_corr_matrix()), kmeans, method = "wss",verbose = T)
      incProgress(0.5, detail = "Complete!")
    })
    corrected_elbow_analysis
  })
  output$corrected_elbow_plot <- renderPlot({
    corrected_elbow_analysis()
  })
  corrected_silhouette_analysis <- reactive({
    req(cluster_mv_features_corr_matrix())
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Running Silhouette Analysis")
      corrected_silhouette_analysis <- fviz_nbclust(x = t(cluster_mv_features_corr_matrix()), kmeans, method = "silhouette",verbose = T)
      incProgress(0.5, detail = "Complete!")
    })
    corrected_silhouette_analysis
  })
  output$corrected_silhouette_plot <- renderPlot({
    corrected_silhouette_analysis()
  })
  corrected_dunn_index_analysis <- reactive({
    req(cluster_mv_features_corr_matrix())
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Running Dunn Index Analysis")
      dunn_k <- c(2:10)
      dunnin <- c()
      for (i in dunn_k){
        dunnin[i] <- dunn(
          distance = dist(t(cluster_mv_features_corr_matrix())),
          clusters = kmeans(t(cluster_mv_features_corr_matrix()), i)$cluster
        )
      }
      corrected_dunn_index_analysis <- as.data.frame(cbind(dunn_k, dunnin[-1]))
      colnames(corrected_dunn_index_analysis) <- c("cluster_number", "dunn_index")
      p <- ggplot(data = corrected_dunn_index_analysis, mapping = aes(x = cluster_number, y = dunn_index))+
        geom_point(color = "dodgerblue1")+
        geom_line(color = "dodgerblue1")+
        geom_vline(
          xintercept = corrected_dunn_index_analysis$cluster_number[which(max(corrected_dunn_index_analysis$dunn_index) == corrected_dunn_index_analysis$dunn_index)],
          color = "dodgerblue1",
          linetype = 2
        ) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14))
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$corrected_dunn_index_plot <- renderPlot({
    corrected_dunn_index_analysis()
  })

  #### Heatmap -----------------------------------------------------------------
  # Generating heatmap of uncorrected data
  corrected_heatmap <- shiny::reactive({
    req(cluster_mv_features_corr_matrix())
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Generating Heatmap")
      colAnn <- heat_colAnn()
      HeatRowNames <- ifelse("Turn on Row Names" %in% input$HeatRowColNames,TRUE,FALSE)
      HeatColNames <- ifelse("Turn on Column Names" %in% input$HeatRowColNames,TRUE,FALSE)

      corrected_matrix_heatmap <- as.matrix(cluster_mv_features_corr_matrix())
      corrected_matrix_heatmap_cols <- colnames(corrected_matrix_heatmap)
      corrected_matrix_scaled <- t(apply(corrected_matrix_heatmap, 1, scale))
      colnames(corrected_matrix_scaled) <- corrected_matrix_heatmap_cols
      p <- suppressMessages(ComplexHeatmap::Heatmap(corrected_matrix_scaled, top_annotation = colAnn,
                                                    show_row_names = HeatRowNames, show_column_names = HeatColNames,
                                                    heatmap_legend_param = list(title = "Expression")))
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$corrected_heatmap <- shiny::renderPlot({
    corrected_heatmap()
  })

  #### RLE ------------------------------------------------
  corrected_batch_choices_RLE2 <- shiny::reactive({
    if (isTruthy(input$batch_choices_RLE)) {
      if (input$batch_choices_RLE == "Select"){
        NULL
      }else {
        corrected_batch_choices_RLE2 <- input$batch_choices_RLE
      }
    }
  })
  RLE_Obj_Corr <- reactive({

    mat_RLE <- batch_correction()
    names(mat_RLE) <- NULL
    corrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_RLE)),
      colData = aligned_meta_file()[order(aligned_meta_file()$Study),],
      rowData = corrected_numeric_matrix2()[,1]
    )
    corrected_RLE_SCE

  })
  corrected_RLE <- shiny::reactive({
    rle_obj <- RLE_Obj_Corr()
    req(rle_obj)
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Running Relative Log Expression")
      if (isTruthy(input$batch_choices_RLE)) {
        p <- scater::plotRLE(
          rle_obj,
          exprs_values = "counts",
          color_by = corrected_batch_choices_RLE2()
        )
      }
      incProgress(0.5, detail = "Complete!")
    })
    p

  })
  output$corrected_RLE_plot <- shiny::renderPlot({
    corrected_RLE()
  })

  corrected_EV_df <- shiny::reactive({
    req(input$variable_choices_EV)
    mat_EV <- batch_correction()
    names(mat_EV) <- NULL
    my_colors <- metafolio::gg_color_hue(length(input$variable_choices_EV))
    corrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_EV)),
      colData = aligned_meta_file(),
      rowData = corrected_numeric_matrix2()[,1]
    )
    SummarizedExperiment::assay(corrected_EV_SCE, "logcounts") <- SingleCellExperiment::counts(corrected_EV_SCE)
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.25, detail = "Running PCA")
      corrected_EV_SCE_PCA <- scater::runPCA(corrected_EV_SCE)
      incProgress(0.25, detail = "Calculating per-gene variance")
      exp_mat_corr <- getVarianceExplained(corrected_EV_SCE_PCA,
                                           exprs_values = "logcounts",
                                           variables = input$variable_choices_EV)
      exp_mat_corr_melt <- reshape2::melt(exp_mat_corr)
      incProgress(0.5, detail = "Complete!")
    })
    exp_mat_corr_melt
  })

  corrected_EV <- reactive({

    exp_mat_corr_melt <- corrected_EV_df()
    p <- ggplot(exp_mat_corr_melt, aes(x = value,color = Var2)) +
      geom_density(size = 1)
    if (input$ExpPlotLog) {
      p <- p + scale_x_log10(limit = c(0.0001,100),labels = ~ format(.x, scientific = FALSE), breaks = c(0.001,0.01,0.1,1,10,100)) +
        geom_vline(xintercept = 1, linetype="dashed")
    }
    p <- p +
      theme_classic() +
      labs(x = "% Variance Explained", y = "Density", color = "Variable")
    p

  })


  output$corrected_EV_plot <- shiny::renderPlot({
    corrected_EV()
  })

  #### PVCA --------------------------------------------------------------------

  pvca_mv_features_corr_matrix <- reactive({
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$pvcacluster_n_MV_features
      var_type <- input$pvcaVarianceMeasure
      mat <- corrected_numeric_matrix2()
      rownames(mat) <- mat[,1]
      mat <- mat[,-1]

      mad <- NULL
      var <- NULL
      cv <- NULL
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(log2(mat + 1), 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(log2(mat + 1), 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(log2(mat + 1), 1, cv))
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = topN)
        out <- cbind(names(cv), cv[names(cv)], mat[names(cv),])
        colnames(out) <- c("Gene", "CV", colnames(mat))
        dataset <- mat[names(cv),]
      }
      incProgress(0.5, detail = "Complete!")
    })

    dataset

  })

  pvca_corr_react <- reactive({

    req(pvca_mv_features_corr_matrix())
    meta <- aligned_meta_file()
    mat <- as.matrix(pvca_mv_features_corr_matrix())
    if (isTruthy(input$variable_choices_EV) & isTruthy(input$pvcaPct)) {
      vars <- input$variable_choices_EV
      if (length(vars) > 2) {
        withProgress(message = "Processing Batch Corrected", value = 0, {
          incProgress(0.5, detail = "Running PCVA")
          pvca_res <- statVisual::PVCA(
            clin_data = meta,                # clinical
            clin_subjid = colnames(meta)[1], # sample name column
            gene_data = mat,                 # expression data
            batch.factors = vars,
            pct_threshold = input$pvcaPct)            # batch columns
          incProgress(0.5, detail = "Complete!")
        })

        pvca_res
      }
    }

  })

  output$corrected_PVCA_plot <- renderPlot({

    p <- pvca_corr_react()
    p

  })

  #### SVA ---------------------------------------------------------------------
  corrected_SVA_variable_of_interest2 <- reactive({
    if (isTruthy(input$uncorrected_SVA_variable_of_interest)) {
      if (input$uncorrected_SVA_variable_of_interest == "Select"){
        NULL
      }else {
        corrected_SVA_variable_of_interest2 <- input$uncorrected_SVA_variable_of_interest
      }
    }
  })
  corrected_SVA_nsv <- reactive({
    if (isTruthy(corrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      withProgress(message = "Processing Batch Corrected", value = 0, {
        incProgress(0.5, detail = "Number of surrogate variables")
        corrected_mod <- model.matrix(reformulate(corrected_SVA_variable_of_interest2()), data = aligned_meta_file())
        corrected_mod_null <- model.matrix(~1, data = aligned_meta_file())
        corrected_SVA_nsv <- sva::num.sv(batch_correction(), corrected_mod, method = svaMethod, vfilter = svaVarNum)
        incProgress(0.5, detail = "Complete!")
      })
      corrected_SVA_nsv
    }
  })
  corrected_SVA_object <- reactive({
    if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      req(corrected_SVA_nsv())
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      withProgress(message = "Processing Batch Corrected", value = 0, {
        incProgress(0.5, detail = "Running SVA")
        corrected_mod <- model.matrix(reformulate(corrected_SVA_variable_of_interest2()), data = aligned_meta_file())
        corrected_mod_null <- model.matrix(~1, data = aligned_meta_file())
        corrected_SVA_matrix <- as.matrix(batch_correction())
        corrected_SVA_object <- sva::sva(corrected_SVA_matrix, corrected_mod, corrected_mod_null, n.sv = corrected_SVA_nsv(),
                                         numSVmethod = svaMethod, vfilter = svaVarNum)
        incProgress(0.5, detail = "Complete!")
      })
      df <- as.data.frame(corrected_SVA_object$sv)
      colnames(df) <- paste0("SVA_",input$batch_correction_method,"_Corrected_Surrogate_Vars_",seq(ncol(df)))
      df <- cbind(aligned_meta_file(),df)
      print(head)
      corr_boxplot_meta_file(df)
      corrected_SVA_object
    }
  })
  corrected_SVA_probability_df <- reactive({
    if (isTruthy(uncorrected_SVA_object())) {
      corrected_SVA_probability_df <- data.frame(
        Genes = 1:length(corrected_SVA_object()$pprob.gam),
        latent_variable = corrected_SVA_object()$pprob.gam,
        variable_of_intrest = corrected_SVA_object()$pprob.b
      )
      corrected_SVA_probability_df_longer <- pivot_longer(corrected_SVA_probability_df, !Genes, names_to = "variable_type")
      corrected_SVA_probability_df_longer <- dplyr::rename(corrected_SVA_probability_df_longer, "probability_association_of_each_gene" = "value")
      corrected_SVA_probability_df_longer
    }
  })

  corrected_SVA_probability_ggplot <- reactive({
    if (isTruthy(uncorrected_SVA_probability_df())) {
      ggplot(corrected_SVA_probability_df(), aes(x = probability_association_of_each_gene,fill = variable_type)) +
        geom_density(alpha = 0.5)
    }
  })

  output$corrected_SVA_probability_density <- renderPlot({
    corrected_SVA_probability_ggplot()
  })

  output$corrected_SVA_nsv_print <- renderPrint({
    if (isTruthy(uncorrected_SVA_nsv())) {
      print(paste("The Number of Estimated Surrogate Variables are:", corrected_SVA_nsv()))
    }
  })

  #### Box Plot ----------------------------------------------------------------

  PCA_Proj_corr_Samples <- reactive({
    if (input$uncorrected_panel == "box") {
      withProgress(message = "Processing Batch Corrected", value = 0, {
        incProgress(0.5, detail = "Running PCA")
        mat_corr <- corrected_numeric_matrix2()
        rownames(mat_corr) <- mat_corr[,1]
        mat_corr <- mat_corr[,-1]
        mat_corr_t <- as.data.frame(t(mat_corr))
        pca_corr <- prcomp(t(mat_corr_t), scale = TRUE)
        rot_corr <- as.data.frame(t(pca_corr[["rotation"]]))
        rot_corr$`Principal Component` <- rownames(rot_corr)
        rot_corr <- rot_corr %>% relocate(`Principal Component`)
        incProgress(0.5, detail = "Complete!")
      })
      rot_corr
    }
  })

  ImmDeconv_corr_react <- reactive({

    req(input$ImmuneDeconvMethods)
    deconvMethod <- input$ImmuneDeconvMethods
    mat <- corrected_numeric_matrix2()
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]
    withProgress(message = "Processing Batch Corrected", value = 0, {
      incProgress(0.5, detail = "Running Immune Deconvolution")
      if (input$HumanOrMouse == "Human") {
        deconv <- as.data.frame(deconvolute(mat, deconvMethod))
      } else {
        deconv <- as.data.frame(deconvolute_mouse(mat, deconvMethod))
      }
      incProgress(0.5, detail = "Complete!")
    })
    deconv

  })

  Mat_for_ssGSEA_corr <- reactive({

    mat <- corrected_numeric_matrix2()
    if (input$HumanOrMouse == "Mouse") {
      mat_conv <- MouseToHuman(mat,MM_HS_Conv)
      mat_conv
    } else {
      mat <- corrected_numeric_matrix2()
      mat
    }

  })

  output$rendcorrected_Box_plot <- renderUI({
    plotHeight <- input$BPplotHeight
    plotWidth <- input$BPplotWidth
    shinyjqui::jqui_resizable(shiny::plotOutput("corrected_Box_plot",height = plotHeight, width = plotWidth))
  })

  CohortBPPlot_corr_df_react <- reactive({

    req(input$BPsampSubset)
    mat <- corrected_numeric_matrix2()
    #meta <- aligned_meta_file()
    meta <- boxplot_meta_file()
    sampSubset <- input$BPsampSubset
    sampCrit <- input$BPsampCrit
    groupCrit <- input$BPgroupCriteria
    FeatCat <- input$BPFeatureCategory
    Feature <- input$BPFeatSelection
    NameCol <- colnames(meta)[1]
    removeSingles <- input$BPremoveSingles
    BPlog <- input$BPlogOpt
    if (FeatCat == "Matrix Features") {
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
    } else if (FeatCat == "PCA Projections") {
      mat <- PCA_Proj_corr_Samples()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
    } else if (FeatCat == "Immune Deconvolution Features") {
      mat <- ImmDeconv_corr_react()
      featdf <- mat[which(mat[,1] == Feature),]
      rownames(featdf) <- featdf[,1]
      featdf <- featdf[,-1]
      featdf <- as.data.frame(t(featdf))
      featdf[,NameCol] <- rownames(featdf)
      meta <- merge(meta,featdf)
    } else if (FeatCat == "Gene Set Pathways") {
      gs_name <- as.character(GeneSetTableBack_react()[input$GeneSetTableUI_rows_selected,ncol(GeneSetTableBack_react())])
      if (isTruthy(gs_name)) {
        gs <- geneset[gs_name]
        Feature <- gs_name
        mat <- Mat_for_ssGSEA_corr()
        rownames(mat) <- mat[,1]
        mat <- mat[,-1]
        mat <- as.matrix(mat)
        withProgress(message = "Processing Batch Corrected", value = 0, {
          incProgress(0.5, detail = "Running ssGSEA")
          ssGSEA_param <- GSVA::ssgseaParam(mat,gs)
          ssGSEA <- GSVA::gsva(ssGSEA_param)
          incProgress(0.5, detail = "Complete!")
        })
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA[,NameCol] <- rownames(ssGSEA)
        meta <- merge(meta,ssGSEA)
      }
    }

    if (isTruthy(Feature)) {
      if (Feature %in% colnames(meta)) {
        if (sampSubset != "Select All Samples") {
          meta <- meta[which(meta[,sampSubset] == sampCrit),]
          meta <- meta %>% select(any_of(c(NameCol,sampSubset,groupCrit,Feature)))
        } else {
          meta <- meta %>% select(any_of(c(NameCol,groupCrit,Feature)))
        }
        if (removeSingles == T) {
          tab <- table(meta[,groupCrit])
          meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
        }
        meta[,Feature] <- as.numeric(meta[,Feature])
        meta <- meta[which(!is.na(meta[,Feature])),]
        meta <- meta[which(meta[,Feature]!=Inf & meta[,Feature]!=-Inf),]
        meta[,groupCrit] <- as.factor(meta[,groupCrit])
        if (BPlog == "Log2") {
          meta[,Feature] <- log2(meta[,Feature])
        } else if (BPlog == "Log2+1") {
          meta[,Feature] <- log2(meta[,Feature]+1)
        } else if (BPlog == "Log10") {
          meta[,Feature] <- log10(meta[,Feature])
        } else if (BPlog == "Log10+1") {
          meta[,Feature] <- log10(meta[,Feature]+1)
        }

        meta
      }
    }

  })

  CohortBPPlot_corr_react <- reactive({

    plotdf_full <- CohortBPPlot_corr_df_react()
    sampSubset <- input$BPsampSubset
    sampCrit <- input$BPsampCrit
    groupCrit <- input$BPgroupCriteria
    FeatCat <- input$BPFeatureCategory
    Feature <- input$BPFeatSelection
    NameCol <- colnames(plotdf_full)[1]
    BPplottheme <- input$BPTheme
    StatMethod <- input$BPplotstatComp
    dotChoice <- input$BPplotsampledots
    dotSize <- input$BPplotDotSize
    bpFlip <- input$BPflipBP
    BPorVI <- input$BPorViolin
    if (FeatCat == "Gene Set Pathways") {
      gs_name <- as.character(GeneSetTableBack_react()[input$GeneSetTableUI_rows_selected,ncol(GeneSetTableBack_react())])
      Feature <- gs_name
    }
    Xaxis_font <- input$BPplot1XAxisSize              # Axis font size
    Yaxis_font <- input$BPplot1YAxisSize              # Axis font size
    Yaxis_lim <- input$BPplot1YAxisLim
    hjust_orient <- 1                                # Initial hjust
    axis_orient <- as.numeric(input$BPxAxisOrient)  # X-axis label orientation
    if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
      hjust_orient <- 0.5
    }
    BPorder <- input$BPplotXaxOrder
    BPGroupSelect <- input$BPgroupSelection

    if (isTruthy(Feature)) {
      plotdf <- plotdf_full[,c(NameCol,groupCrit,Feature)]
      plotdf <- plotdf[which(plotdf[,groupCrit] %in% BPGroupSelect),]
      colnames(plotdf) <- c("SampleName","Group","Feature")

      if (BPorder == "Descending"){
        barp <- ggplot(data = plotdf, aes(x=reorder(Group,-Feature, FUN = median),y=Feature, fill=Group))
        plotdf_dots <- plotdf
        plotdf_dots$Group <- reorder(plotdf_dots$Group,-plotdf_dots$Feature, FUN = median)
        plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
      }
      if (BPorder == "Ascending"){
        barp <- ggplot(data = plotdf, aes(x=reorder(Group,Feature, FUN = median),y=Feature, fill=Group))
        plotdf_dots <- plotdf
        plotdf_dots$Group <- reorder(plotdf_dots$Group,plotdf_dots$Feature, FUN = median)
        plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
      }
      if (BPorder == "Not Specificed"){
        barp <- ggplot(data = plotdf, aes(x=Group,y=Feature, fill=Group))
        plotdf_dots <- plotdf
        plotdf_dots$xj <- jitter(as.numeric(factor(plotdf_dots$Group)))
      }
      if (BPorVI == "Box Plot") {
        barp <- barp + geom_boxplot(width = 0.5, lwd = 1)
      }
      if (BPorVI == "Violin Plot") {
        barp <- barp + geom_violin() +
          stat_summary(fun=median, geom="crossbar", width=0.5, color="black")
      }
      if (isTruthy(Yaxis_lim)) {
        barp <- barp +
          ylim(paste0(as.numeric(strsplit(Yaxis_lim,",")[[1]][1]),as.numeric(strsplit(Yaxis_lim,",")[[1]][2])))
      }
      barp <- barp +
        get(BPplottheme)() +
        labs(#title = BPTitle_in,
          x = groupCrit, y = Feature,
          fill = groupCrit)
      if (StatMethod != "none") {
        barp <- barp + ggpubr::stat_compare_means(method = StatMethod)
      }
      if (dotChoice) {
        barp <- barp + geom_point(data = plotdf_dots, aes(x=xj), col="grey14", size=dotSize)
      }
      barp <- barp + theme(axis.text.x = element_text(size = Xaxis_font, angle = axis_orient, hjust = hjust_orient),
                           axis.title.x = element_text(size = Xaxis_font),
                           axis.text.y = element_text(size = Yaxis_font),
                           axis.title.y = element_text(size = Yaxis_font),
                           legend.position = "none")
      if (bpFlip) {
        barp <- barp + coord_flip()
      }
      barp
    }

  })

  output$corrected_Box_plot <- renderPlot({

    barp <- CohortBPPlot_corr_react()
    barp

  })


  ### END OF THE OUTPUT GENERATING CODE ----------------------------------------

  ### START OF THE SAVING FILES TO ZIP CODE ------------------------------------
  #### Uncorrected -------------------------------------------------------------
  # Saving files to ZIP
  dir.create(temp_directory)
  shiny::observeEvent(input$save_uncorrected_matrix, {
    file_name <- paste("uncorrected_matrix", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(uncorrected_matrix(), file.path(temp_directory, file_name))
    print("Uncorrected_matrix Ready for Zip")
  })
  output$dnldsave_uncorrected_matrix <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Matrix","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- uncorrected_matrix()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorrected_PCA_plot, {
    file_name <- paste("uncorrected_PCA_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_PCA_plot_forDlnd_react(),
                    height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
    print("Uncorrected_PCA_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_PCA_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_PCA","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_PCA_plot_forDlnd_react()
      ggplot2::ggsave(file,p, height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_PCA_mc_plot, {
    file_name <- paste("uncorrected_PCA_mc_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_PCA_multiple_components(),
                    height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
    print("Uncorrected_PCA_mc_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_PCA_mc_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_PCA_MultipleComponents","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_PCA_multiple_components()
      ggplot2::ggsave(file,p, height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_scree_plot, {
    file_name <- paste("uncorrected_Scree_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_scree_plot_react(),
                    height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
    print("Uncorrected_Scree_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_scree_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Scree_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_scree_plot_react()
      ggplot2::ggsave(file,p, height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_PCA_components, {
    file_name <- paste("uncorrected_PC_components", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(as.data.frame(uncorrected_PCA_details2()$x), file.path(temp_directory, file_name))
    print("uncorrected_PCA_details Ready for Zip")
  })
  output$dnld <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_PC_Components","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- save_uncorrected_PCA_components()
      write.table(as.data.frame(df$x),file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorrected_contribution_table, {
    file_name <- paste("uncorrected_contribution_table", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(uncorrected_PCA_individuals(), file.path(temp_directory, file_name))
    print("uncorrected_contribution_table Ready for Zip")
  })
  output$dnldsave_uncorrected_contribution_table <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Contribution_Table","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- uncorrected_PCA_individuals()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorrected_contribution_counts, {
    file_name <- paste("uncorrected_contribution_counts", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(uncorrected_contribution_counts(), file.path(temp_directory, file_name))
    print("uncorrected_contribution_counts Ready for Zip")
  })
  output$dnldsave_uncorrected_contribution_counts <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Contribution_Counts","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- uncorrected_contribution_counts()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorrected_UMAP, {
    file_name <- paste("uncorrected_UMAP", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_UMAP_react(),
           height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
    print("uncorrected_UMAP Ready for Zip")
  })
  output$dnldsave_uncorrected_UMAP <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_UMAP","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_UMAP_react()
      ggplot2::ggsave(file,p, height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
    }
  )
  observeEvent(input$save_uncorrected_elbow_plot, {
    file_name <- paste("uncorrected_elbow_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_elbow_analysis(),
           height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    print("uncorrected_elbow_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_elbow_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Elbow_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_elbow_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_uncorrected_silhouette_plot, {
    file_name <- paste("uncorrected_silhouette_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_silhouette_analysis(),
                    height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    print("uncorrected_silhouette_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_silhouette_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Silhouette_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_silhouette_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_uncorrected_dunn_index_plot, {
    file_name <- paste("uncorrected_dunn_index_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_dunn_index_analysis(),
                    height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    print("uncorrected_dunn_index_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_dunn_index_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Dunn_Index_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_dunn_index_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_heatmap, {
    file_name <- paste(temp_directory, "/", "uncorrected_heatmap", Sys.Date(),".svg", sep = "")
    svg(filename = file_name, height = input$heatmapHeight, width = input$heatmapWidth)
    ComplexHeatmap::draw(uncorrected_heatmap())
    dev.off()
    print("Uncorrected_heatmap Ready for Zip")
  })

  output$dnldsave_uncorrected_heatmap <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Heatmap","_",Sys.Date(),".svg")
    },
    content = function(file) {
      svg(filename = file, height = input$heatmapHeight, width = input$heatmapWidth)
      ComplexHeatmap::draw(uncorrected_heatmap())
      dev.off()
    }
  )
  shiny::observeEvent(input$save_uncorrected_RLE_plot, {
    file_name <- paste("uncorrected_RLE_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_RLE(),
                    height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
    print("Uncorrected_RLE_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_RLE_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_RLE_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_RLE()
      ggplot2::ggsave(file,p, height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_EV_plot, {
    file_name <- paste("uncorrected_EV_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_EV(),
                    height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
    print("Uncorrected_EV_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_EV_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_EV_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_EV()
      ggplot2::ggsave(file,p, height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_pvca_plot, {
    file_name <- paste("Uncorrected_PVCA_Plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = pvca_uncorr_react(),
                    height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
    print("Uncorrected_PVCA_Plot Ready for Zip")
  })
  output$dnldsave_uncorrected_pvca_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_PVCA_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- pvca_uncorr_react()
      ggplot2::ggsave(file,p, height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
    }
  )
  observeEvent(input$save_uncorrected_SVA_probability_density, {
    file_name <- paste("uncorrected_SVA_probability_density", Sys.Date(),".svg", sep = "")
    ggsave(filename = file_name, path = temp_directory, plot = uncorrected_SVA_probability_ggplot(),
           height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
    print("Uncorrected_SVA_probability_density Ready for Zip")
  })
  output$dnldsave_uncorrected_SVA_probability_density <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_SVA_Probability_Density","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- uncorrected_SVA_probability_ggplot()
      ggplot2::ggsave(file,p, height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
    }
  )
  observeEvent(input$save_SVA_surrogate_variables, {
    file_name <- paste("Meta_withSV_", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(boxplot_meta_file(), file.path(temp_directory, file_name))
    print("surrogate_variables Ready for Zip")
  })
  output$dnldsave_SVA_surrogate_variables <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Meta_withSV","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- boxplot_meta_file()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  observeEvent(input$save_uncorrected_Box_plot, {
    file_name <- paste("uncorrected_Box_plot", Sys.Date(),".svg", sep = "")
    ggsave(filename = file_name, path = temp_directory, plot = CohortBPPlot_react(),
           height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
    print("uncorrected_Box_plot Ready for Zip")
  })
  output$dnldsave_uncorrected_Box_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Uncorrected_Box_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- CohortBPPlot_react()
      ggplot2::ggsave(file,p, height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
    }
  )

  #### Corrected ---------------------------------------------------------------
  shiny::observeEvent(input$save_corrected_matrix, {
    file_name <- paste(input$batch_correction_method, "_corrected_matrix", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(corrected_numeric_matrix2(), file.path(temp_directory, file_name))
    print("Corrected_matrix Ready for Zip")
  })
  output$dnldsave_corrected_matrix <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Matrix","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- corrected_numeric_matrix2()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corrected_PCA_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_PCA_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_PCA_plot_forDlnd_react(),
           height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
    print("Corrected_PCA_plot Ready for Zip")
  })
  output$dnldsave_corrected_PCA_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_PCA_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_PCA_plot_forDlnd_react()
      ggplot2::ggsave(file,p, height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_PCA_mc_plot, {
    file_name <- paste(input$batch_correction_method,"corrected_PCA_mc_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_PCA_multiple_components(),
                    height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
    print("Corrected_PCA_mc_plot Ready for Zip")
  })
  output$dnldsave_corrected_PCA_mc_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_PCA_Multiple_Components_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_PCA_multiple_components()
      ggplot2::ggsave(file,p, height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_scree_plot, {
    file_name <- paste(input$batch_correction_method,"corrected_Scree_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_scree_plot_react(),
                    height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
    print("Corrected_Scree_plot Ready for Zip")
  })
  output$dnldsave_corrected_scree_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Scree_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_scree_plot_react()
      ggplot2::ggsave(file,p, height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_PCA_components, {
    file_name <- paste(input$batch_correction_method, "_corrected_PC_components", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(as.data.frame(corrected_PCA_details2()$x), file.path(temp_directory, file_name))
    print("Corrected_PCA_details Ready for Zip")
  })
  output$dnldsave_corrected_PCA_components <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","_Corrected_PC_Components","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- corrected_PCA_details2()
      write.table(as.data.frame(df$x),file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corrected_contribution_table, {
    file_name <- paste(input$batch_correction_method,"corrected_contribution_table", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(corrected_PCA_individuals(), file.path(temp_directory, file_name))
    print("corrected_contribution_table Ready for Zip")
  })
  output$dnldsave_corrected_contribution_table <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Contribution_Table","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- corrected_PCA_individuals()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corrected_contribution_counts, {
    file_name <- paste(input$batch_correction_method,"corrected_contribution_counts", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(corrected_contribution_counts(), file.path(temp_directory, file_name))
    print("corrected_contribution_counts Ready for Zip")
  })
  output$dnldsave_corrected_contribution_counts <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Contribution_Counts","_",Sys.Date(),".tsv")
    },
    content = function(file) {
      df <- corrected_contribution_counts()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corrected_UMAP, {
    file_name <- paste(input$batch_correction_method,"corrected_UMAP", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_UMAP_react(),
           height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
    print("corrected_UMAP Ready for Zip")
  })
  output$dnldsave_corrected_UMAP <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_UMAP","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_UMAP_react()
      ggplot2::ggsave(file,p, height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
    }
  )
  observeEvent(input$save_corrected_elbow_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_elbow_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_elbow_analysis(),
                    height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    print("corrected_elbow_plot Ready for Zip")
  })
  output$dnldsave_corrected_elbow_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Elbow_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_elbow_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_corrected_silhouette_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_silhouette_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_silhouette_analysis(),
                    height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    print("corrected_silhouette_plot Ready for Zip")
  })
  output$dnldsave_corrected_silhouette_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Silhouette_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_silhouette_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_corrected_dunn_index_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_dunn_index_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_dunn_index_analysis(),
                    height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    print("corrected_dunn_index_plot Ready for Zip")
  })
  output$dnldsave_corrected_dunn_index_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Dunn_Index_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_dunn_index_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_heatmap, {
    file_name <- paste(temp_directory,"/",input$batch_correction_method, "_corrected_heatmap", Sys.Date(),".svg", sep = "")
    svg(filename = file_name, height = input$heatmapHeight, width = input$heatmapWidth)
    ComplexHeatmap::draw(corrected_heatmap())
    dev.off()
    print("Corrected_heatmap Ready for Zip")
  })
  output$dnldsave_corrected_heatmap <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Heatmap","_",Sys.Date(),".svg")
    },
    content = function(file) {
      svg(filename = file, height = input$heatmapHeight, width = input$heatmapWidth)
      ComplexHeatmap::draw(corrected_heatmap())
      dev.off()
    }
  )
  shiny::observeEvent(input$save_corrected_RLE_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_RLE_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_RLE(),
                    height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
    print("Corrected_RLE_plot Ready for Zip")
  })
  output$dnldsave_corrected_RLE_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_RLE_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_RLE()
      ggplot2::ggsave(file,p, height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_EV_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_EV_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_EV(),
                    height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
    print("Corrected_EV_plot Ready for Zip")
  })
  output$dnldsave_corrected_EV_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_EV_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_EV()
      ggplot2::ggsave(file,p, height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_pvca_plot, {
    file_name <- paste("Corrected_PVCA_Plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = pvca_corr_react(),
                    height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
    print("Corrected_PVCA_Plot Ready for Zip")
  })
  output$dnldsave_corrected_pvca_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Corrected_PVCA_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- pvca_corr_react()
      ggplot2::ggsave(file,p, height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
    }
  )
  observeEvent(input$save_corrected_SVA_probability_density, {
    file_name <- paste(input$batch_correction_method, "_corrected_SVA_probability_density", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_SVA_probability_ggplot(),
                    height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
    print("Corrected_SVA_probability_density Ready for Zip")
  })
  output$dnldsave_corrected_SVA_probability_density <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_SVA_Probability_Density","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- corrected_SVA_probability_ggplot()
      ggplot2::ggsave(file,p, height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
    }
  )
  observeEvent(input$save_corrected_Box_plot, {
    file_name <- paste(input$batch_correction_method,"_corrected_Box_plot", Sys.Date(),".svg", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = CohortBPPlot_corr_react(),
                    height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
    print("corrected_Box_plot Ready for Zip")
  })
  output$dnldsave_corrected_Box_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",input$batch_correction_method,"_","Corrected_Box_Plot","_",Sys.Date(),".svg")
    },
    content = function(file) {
      p <- CohortBPPlot_corr_react()
      ggplot2::ggsave(file,p, height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
    }
  )
  has.new.files <- function() {
    unique(list.files(temp_directory))
  }
  get.files <- function() {
    list.files(temp_directory)
  }
  my_files <- reactivePoll(10, session, checkFunc=has.new.files, valueFunc=get.files)
  observeEvent(my_files(),ignoreInit = T,ignoreNULL = T, {
    print(my_files())
    updatePickerInput(session, "select_save_files", choices = my_files())
  })
  selected_files_to_download <- reactive({
    selected_files_to_download <- c()
    for (i in input$select_save_files){
      file_to_download <- paste(temp_directory, "/", i, sep = "")
      selected_files_to_download <- append(selected_files_to_download, file_to_download)
    }
    selected_files_to_download
  })
  output$download_btn <- downloadHandler(
    filename = function(){
      paste(input$file_name, "_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file){
      zip::zipr(
        zipfile = file,
        files = selected_files_to_download(),
        root = temp_directory
      )
    },
    contentType = "application/zip"
  )
}

# Run the application
shinyApp(ui = ui, server = server)
