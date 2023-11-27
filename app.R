

# Load Libraries ---------------------------------------------------------------
options(shiny.maxRequestSize = 10000*1024^2)
#packages <- c("shiny","shinyjqui","limma", "bapred","DT","shinycssloaders",
#              "metafolio", "corrplot","readr","dplyr",
#              "cluster","zip", "glue", "factoextra", "FactoMineR")
#installed_packages <- packages %in% rownames(installed.packages())
#invisible(lapply(packages, library, character.only = TRUE))
#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], ask = F)
#}
#bioCpacks <- c("sva", "preprocessCore", "scater", "affyPLM")
#installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
#if (any(installed_packages_BIOC == FALSE)) {
#  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
#}
#invisible(lapply(bioCpacks, library, character.only = TRUE))
#
#library(ggfortify)

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
library(sva)
library(preprocessCore)
library(scater)
library(affyPLM)
library(ggfortify)
library(ggplotify)
library(draw)


# Function to remove genes that are lowly expressed
ExprFilter2 <- function(vec,criteria,proportion) {
  Samp2meet <- length(vec)*proportion
  meet <- sum(vec>criteria)
  if (meet > Samp2meet) {
    return(TRUE)
  } else { return(FALSE) }
}

# UI ---------------------------------------------------------------------------
ui <- 
  shiny::navbarPage("{ Batch-Flex }",
                    ## Data Input Tab ------------------------------------------
                    shiny::tabPanel("Step 1 - Data Input",
                                    shiny::sidebarLayout(
                                      ### Sidebar ------------------------------
                                      shiny::sidebarPanel(
                                        # File to be corrected
                                        shiny::fluidRow(
                                          shiny::column(9,
                                                        shiny::fileInput("uncorrected_matrix_input", 
                                                                           "Matrix for Batch Correction"
                                                                         )
                                                        ),
                                          ## Deliminator for input matrix
                                          shiny::column(3,shiny::selectInput("matrix_delim", 
                                                                             "Deliminator", 
                                                                             c("Tab" = '\t',"Comma" = ',', "Space" = ' ',"Semicolon" = ';', "Colon" = ':')
                                                                             )
                                                        )
                                          ),
                                        # Need to log?
                                        shiny::fluidRow(
                                          shiny::column(4, style = "margin-top:-20px",
                                                        shiny::checkboxInput("log_transform", "Log2 + 1 Transform?")
                                          ),
                                          shiny::column(4, style = "margin-top:-20px",
                                                        shiny::checkboxInput("quantile_normalization", "Quantile Normalize?")
                                          )
                                        ),
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
                                        shiny::uiOutput("batch_names")
                                        # if batch information is provided by user
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
                                                            shiny::column(11,
                                                                          shiny::selectInput("batch_correction_method",
                                                                                             "Select the Method of Batch Correction",
                                                                                             c("Select", "Limma", "ComBat", "Mean Centering", "ComBatseq"))
                                                                          )
                                                          ),
                                                          ##### Limma Sidebar----------------------------
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
                                                          ##### ComBat Sidebar----------------------------
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
                                                          ##### Mean Centering Sidebar----------------------------
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
                                                          ##### ComBatseq Sidebar----------------------------
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
                                                          hr(),
                                                          #### Plot Inputs ------------------------------------
                                                          ##### PCA Sidebar------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca'",
                                                                                  shiny::h3("PCA Parameters"),
                                                                                  shiny::radioButtons("PCA_type",NULL,c("PCA with clustering","PCA without clustering"),inline = T),
                                                                                  shiny::conditionalPanel(condition = "input.PCA_type == 'PCA with clustering'",
                                                                                                          shiny::numericInput("cluster_number","Number of Clusters",
                                                                                                                              value = 2,step = 1,min = 1, width = "50%")
                                                                                                          ),
                                                                                  shiny::conditionalPanel(condition = "input.PCA_type == 'PCA without clustering'",
                                                                                                          shiny::fluidRow(
                                                                                                            shiny::column(6,
                                                                                                                           shiny::uiOutput("batch_choices_PCA")
                                                                                                                          ),
                                                                                                            shiny::column(6,
                                                                                                                           shiny::uiOutput("biological_choices_PCA")
                                                                                                                          )
                                                                                                            )
                                                                                                          )
                                                                                  ),
                                                          ##### PCA MC Sidebar------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_mc'",
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
                                                          ##### PCA Details Sidebar------------------------------------
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_dt'",
                                                                                  shiny::h3("PCA Detail Parameters"),
                                                                                  shiny::uiOutput("PCA_factors_choices")
                                                                                  ),
                                                          ##### RLE Sidebar------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'rle'",
                                                                                  shiny::h3("RLE Parameters"),
                                                                                  shiny::uiOutput("batch_choices_RLE")
                                                                                  ),
                                                          ##### Exp Var Sidebar------------------------------------
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'exp_var'",
                                                                                  shiny::h3("Explanatory Variables Parameters"),
                                                                                  shiny::uiOutput("variable_choices_EV")
                                                                                  )
                                                          #shiny::uiOutput("batch_selection")
                                                          #textOutput("test_print")
                                                          ),
                                          #### File Export Sidebar------------------------------------
                                          shiny::tabPanel("File Export",
                                                          shiny::p(),
                                                          shiny::textInput("file_name", "Please Name the Zipped File"),
                                                          shiny::downloadButton(
                                                            outputId = "download_btn",
                                                            label = "Download All Saved Files as Zip",
                                                            icon = icon("file-download")
                                                            )
                                                          )
                                          )
                                        ),
                                      ### Main ---------------------------------
                                        mainPanel(
                                          
                                          shiny::tabsetPanel(
                                            id = "uncorrected_panel",
                                            #### Matrix Main ---------------------------------
                                            shiny::tabPanel(
                                              "Matrix", 
                                              shiny::p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected Matrix"),
                                                              DT::dataTableOutput("uncorrected_matrix_output"),
                                                              shiny::actionButton("save_uncorrected_matrix", "Save Uncorrected Matrix")
                                                ),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected Matrix"),
                                                              DT::dataTableOutput("corrected_matrix"),
                                                              shiny::uiOutput("rendsave_corrected_matrix")
                                                )
                                              ),
                                              value = "mat"
                                            ),
                                            #### PCA Main ---------------------------------
                                            shiny::tabPanel(
                                              "PCA Plot", 
                                              p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected PCA"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("uncorrected_PCA"), type = 6)),
                                                              shiny::actionButton("save_uncorrected_PCA_plot","Save Uncorrected PCA Plot")
                                                ),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected PCA"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("corrected_PCA"), type = 6)),
                                                              shiny::actionButton("save_corrected_PCA_plot", "Save Corrected PCA Plot")
                                                )
                                              ),
                                              value = "pca"
                                            ),
                                            #### PCA MC Main ---------------------------------
                                            shiny::tabPanel(
                                              "Multiple Components PCA",
                                              shiny::p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected Multiple Components PCA"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("uncorrected_PCA_multiple_components"), type = 6)),
                                                              shiny::actionButton("save_uncorrected_PCA_mc_plot","Save Uncorrected PCA Multiple Components Plot")
                                                ),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected Multiple Components PCA"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("corrected_PCA_multiple_components"), type = 6)),
                                                              shiny::actionButton("save_corrected_PCA_mc_plot", "Save Corrected PCA Multiple Components Plot")
                                                )
                                              ),
                                              value = "pca_mc"),
                                            #### PCA Details Main ---------------------------------
                                            shiny:: tabPanel(
                                              "PCA Details",
                                              shiny::p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected PCA Details"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("uncorrected_scree_plot"), type = 6)),
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
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("corrected_scree_plot"), type = 6)),
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
                                            #shiny::tabPanel(
                                            #  "Heatmap", 
                                            #  shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_heatmap")),
                                            #  shiny::actionButton("save_uncorrected_heatmap", "Save Uncorrected Heatmap")
                                            #  ),
                                            #### RLE Main ---------------------------------
                                            shiny::tabPanel(
                                              "RLE Plot",
                                              p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected RLE Plot"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("uncorrected_RLE_plot"), type = 6)),
                                                              shiny::actionButton("save_uncorrected_RLE_plot", "Save Uncorrected RLE Plot"),),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected RLE Plot"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("corrected_RLE_plot"), type = 6)),
                                                              shiny::actionButton("save_corrected_RLE_plot", "Save Corrected RLE Plot")
                                                )
                                              ),
                                              value = "rle"
                                            ),
                                            #### Exp Var Main ---------------------------------
                                            shiny::tabPanel(
                                              "Explanatory Variables Plot",
                                              p(),
                                              shiny::fluidRow(
                                                shiny::column(6,
                                                              shiny::h3("Uncorrected Explanatory Variables Plot"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("uncorrected_EV_plot"), type = 6)),
                                                              shiny::actionButton("save_uncorrected_Ev_plot", "Save Uncorrected EV Plot"),),
                                                shiny::column(6,
                                                              shiny::h3("Batch Corrected Explanatory Variables Plot"),
                                                              shinyjqui::jqui_resizable(shinycssloaders::withSpinner(shiny::plotOutput("corrected_EV_plot"), type = 6)),
                                                              shiny::actionButton("save_corrected_EV_plot", "Save Corrected EV Plot")
                                                )
                                              ),
                                              value = "exp_var"
                                            )
                                          )
                                          )
                                      )
                                    )
                    )


# Read and Manipulate Input Files

# Server -----------------------------------------------------------------------
server <- function(input, output) {
  
  ## START OF UNCORRECTED PLOT DATA --------------------------------------------
  ### Input File Processing ----------------------------------------------------
  # Converting letter input to symbol for deliminator
  #deliminator <- shiny::reactive({
  #  if(input$matrix_delim == "Comma") {
  #    deliminator <-  ","
  #  }else if (input$matrix_delim == c("Tab")) {
  #    deliminator <-  "\t"
  #  }else if (input$matrix_delim == c("Semicolon")) {
  #    deliminator <-  ";"
  #  }else if (input$matrix_delim == c("Colon")) {
  #    deliminator <-  ":"
  #  }
  #})
  
  # reactive value
  files_to_download <- shiny::reactiveValues()
  
  #### Matrix Processing -------------------------------------------------------
  # developing reactive variable for the input expression matrix to be used downstream
  uncorrected_matrix <- shiny::reactive({
    req(input$uncorrected_matrix_input)
    #uncorrected_matrix <- readr::read_delim(input$uncorrected_matrix_input$datapath,
    #                                        delim = deliminator(),
    #                                        name_repair = "minimal",
    #                                        na = c("", "NA", "N/A"))
    uncorrected_matrix <- as.data.frame(readr::read_delim(input$uncorrected_matrix_input$datapath,
                                            delim = input$matrix_delim,
                                            name_repair = "minimal",
                                            na = c("", "NA", "N/A")))
    colnames(uncorrected_matrix) <- gsub("[[:punct:]]","_",colnames(uncorrected_matrix))
    colnames(uncorrected_matrix) <- gsub(" ","_",colnames(uncorrected_matrix))
    uncorrected_matrix <- uncorrected_matrix[, !duplicated(colnames(uncorrected_matrix))]
    uncorrected_matrix
  })
  
  # Data processing for input matrix of gene expression data
  uncorrected_matrix_filtered <- shiny::reactive({
    req(uncorrected_matrix())
    uncorrected_matrix_unfiltered <- uncorrected_matrix()
    colnames(uncorrected_matrix_unfiltered)[1] <- "Genes"
    if (TRUE %in% duplicated(uncorrected_matrix_unfiltered[, 1])) {
      uncorrected_matrix_nodupes <- uncorrected_matrix_unfiltered %>%
        dplyr::group_by(Genes) %>%
        dplyr::summarise_all(max)
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
    if (input$log_transform == TRUE){
      uncorrected_numeric_matrix <- log2(uncorrected_matrix_filtered()[,-1] + 1)
    }else if(input$quantile_normalization == TRUE){
      uncorrected_numeric_matrix <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_matrix_filtered()[,-1]))
      colnames(uncorrected_numeric_matrix) <- colnames(uncorrected_matrix_filtered()[,-1])
      uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix)
    }else {
      uncorrected_numeric_matrix <- uncorrected_matrix_filtered()[,-1]
    }
    uncorrected_numeric_matrix <- cbind(uncorrected_numeric_matrix, uncorrected_matrix_filtered()[,1,drop=FALSE])
    uncorrected_numeric_matrix <- uncorrected_numeric_matrix %>% dplyr::relocate(Genes)
    uncorrected_numeric_matrix
  })
  # rendering input table for visualization
  output$uncorrected_matrix_output <- DT::renderDataTable({
    uncorrected_matrix_output <- uncorrected_numeric_matrix()[,-1]
    uncorrected_matrix_output$Genes <- uncorrected_numeric_matrix()[,1]
    uncorrected_matrix_output <- uncorrected_matrix_output %>% dplyr::relocate(Genes)
    uncorrected_matrix_output
    DT::datatable(uncorrected_numeric_matrix(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
    
  })
  output$rendMatrixHeader <- renderUI({
    req(input$uncorrected_matrix_input)
    h3("Matrix Preview")
  })
  output$uncorrected_matrix_output_input <- DT::renderDataTable({
    req(input$uncorrected_matrix_input)
    uncorrected_matrix_output <- uncorrected_numeric_matrix()[,-1]
    uncorrected_matrix_output$Genes <- uncorrected_numeric_matrix()[,1]
    uncorrected_matrix_output <- uncorrected_matrix_output %>% dplyr::relocate(Genes)
    uncorrected_matrix_output
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
                      shiny::selectInput("meta_delim", "Deliminator?", c("Tab" = '\t',"Comma" = ',', "Space" = ' ',"Semicolon" = ';', "Colon" = ':'))
        )
      )
    }
    
  })
  
  # determining if meta table should be derived or inputted
  output$batch_delim <- shiny::renderUI({
    if(input$batch_info == "Yes"){
      shiny::textInput("batch_name_delim", 
                       "Deliminator to Derive Batch Information from Column Names")
    } #else if(input$batch_info == "No"){
    #fileInput("user_provided_batch_info", 
    #          "Please Provide Meta File")
    #}
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
  
  # Deliminator if user to provide meta file directly
  #output$meta_delim <- shiny::renderUI({
  #  if(input$batch_info == "No"){
  #selectInput("meta_delim", 
  #            "Deliminator?", 
  #            c("Select", "Comma", "Tab", "Semicolon", "Colon"))
  #  }
  #})
  #deliminator2 <- shiny::reactive({
  #  if(input$meta_delim == "Comma") {
  #    deliminator <-  ","
  #  }else if (input$meta_delim == c("Tab")) {
  #    deliminator <-  "\t"
  #  }else if (input$meta_delim == c("Semicolon")) {
  #    deliminator <-  ";"
  #  }else if (input$meta_delim == c("Colon")) {
  #    deliminator <-  ":"
  #  }
  #})
  
  # Reading user provided meta file
  user_provided_batch_info <- shiny::reactive({
    req(input$user_provided_batch_info)
    #user_provided_batch_info <- readr::read_delim(input$user_provided_batch_info$datapath,
    #                                              delim = deliminator2(),
    #                                              na = c("", "NA", "N/A"))
    user_provided_batch_info <- as.data.frame(readr::read_delim(input$user_provided_batch_info$datapath,
                                                  delim = input$meta_delim,
                                                  na = c("", "NA", "N/A")))
    #colnames(user_provided_batch_info) <- gsub("[[:punct:]]","_",colnames(user_provided_batch_info))
    #colnames(user_provided_batch_info) <- gsub(" ","_",colnames(user_provided_batch_info))
    #user_provided_batch_info <- apply(
    #  user_provided_batch_info, 
    #  2, 
    #  function(x) gsub("(?<=[a-z]|[A-Z])-(?=[a-z]|[A-Z])", "_", x, perl = TRUE)))
    user_provided_batch_info
  })
  
  # Attaching and rendering user generated or user provided meta file
  meta_file_react <- shiny::reactive({
    if (input$batch_info == "Yes"){
      if (isTruthy(meta_file_delim())) {
        meta_names <- colnames(uncorrected_matrix())
        meta_names <- as.data.frame(meta_names[-1])
        if (length(meta_file_names()) == ncol(meta_names)+1) {
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
        }
      }
    } else if (input$batch_info == "No"){
      user_provided_batch_info()
    }
  })
  aligned_meta_file <- shiny::reactive({
    unaligned_meta_file <- meta_file_react()
    #print(head(colnames(unaligned_meta_file)))
    #print(head(unaligned_meta_file[,1]))
    # getting meta column name that contains sample names
    column_containing_ids <- colnames(unaligned_meta_file)[grepl(paste(gsub("[[:punct:]]","",colnames(uncorrected_numeric_matrix)),collapse = "|"), gsub("[[:punct:]]","",unaligned_meta_file))][1]
    # reformat sample names in meta file
    unaligned_meta_file[,column_containing_ids] <- gsub("[[:punct:]]","_",unaligned_meta_file[,column_containing_ids])
    unaligned_meta_file[,column_containing_ids] <- gsub(" ","_",unaligned_meta_file[,column_containing_ids])
    # filter to samples in matrix
    unaligned_meta_file_filtered <- unaligned_meta_file[which(unaligned_meta_file[,column_containing_ids] %in% colnames(uncorrected_numeric_matrix())),]
    # align sample names in meta file with matrix, make sure same order
    aligned_meta_file_filtered <- unaligned_meta_file_filtered[match(colnames(uncorrected_numeric_matrix()[,-1]), unaligned_meta_file_filtered[, column_containing_ids]),]
    # relocate sample names to first column
    aligned_meta_file_filtered <- aligned_meta_file_filtered %>% dplyr::relocate(any_of(column_containing_ids))
    #unaligned_meta_file$filter <- apply(unaligned_meta_file, 1, function(x) any(x %in% as.vector(colnames(uncorrected_numeric_matrix()[,-1]))))
    ##print(head(colnames(unaligned_meta_file)))
    ##print(head(unaligned_meta_file[,1]))
    #unaligned_meta_file_filtered <- unaligned_meta_file[which(unaligned_meta_file$filter == TRUE),]
    ##print(head(colnames(unaligned_meta_file_filtered)))
    ##print(head(unaligned_meta_file_filtered[,1]))
    #unaligned_meta_file_filtered <- unaligned_meta_file_filtered[,-ncol(unaligned_meta_file_filtered)]
    ##print(head(colnames(unaligned_meta_file_filtered)))
    ##print(head(unaligned_meta_file_filtered[,1]))
    #column_containing_ids <- colnames(unaligned_meta_file_filtered)[grepl(colnames(uncorrected_numeric_matrix())[2], unaligned_meta_file_filtered)]
    ##print(head(column_containing_ids))
    #aligned_meta_file_filtered <- unaligned_meta_file_filtered[match(colnames(uncorrected_numeric_matrix()[,-1]), unaligned_meta_file_filtered[, column_containing_ids[1]]),]
    ##print(head(colnames(aligned_meta_file_filtered)))
    ##print(head(aligned_meta_file_filtered[,1]))
    aligned_meta_file_filtered
  })
  
  output$rendMetaHeader <- renderUI({
    req(input$user_provided_batch_info)
    h3("Meta Preview")
  })
  output$meta_file <- DT::renderDataTable({
    
    if (input$batch_info == "Yes"){
      if (isTruthy(meta_file_delim())) {
        if (length(meta_file_names()) == ncol(aligned_meta_file())+1) {
          DT::datatable(aligned_meta_file(),
                        options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                       pageLength = 10,
                                       scrollX = T),
                        rownames = F)
        }
      }
    } else {
      DT::datatable(aligned_meta_file(),
                    options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                   pageLength = 10,
                                   scrollX = T),
                    rownames = F)
    }
    
  })
  
  #### PCA ---------------------------------------------------------------------
  # Generating PCA plot
  output$batch_choices_PCA <- shiny::renderUI({
    shiny::selectInput("batch_choices_PCA", 
                       "Group Color", 
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  uncorrected_batch_choices_PCA2 <- shiny::reactive({
    if (!is.null(input$batch_choices_PCA)) {
      if (input$batch_choices_PCA == "Select"){
        NULL
      }else {
        uncorrected_batch_choices_PCA2 <- input$batch_choices_PCA
      }
    }
    
  })
  output$biological_choices_PCA <- shiny::renderUI({
    shiny::selectInput("biological_choices_PCA", 
                       "Group Shape", 
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  uncorrected_biological_choices_PCA2 <- shiny::reactive({
    if (input$biological_choices_PCA == "Select"){
      NULL
    }else {
      uncorrected_biological_choices_PCA2 <- input$biological_choices_PCA
    }
  })
  PCA_data <- shiny::reactive({
    PCA_data <- cbind((t(uncorrected_numeric_matrix()[,-1])), aligned_meta_file())
    PCA_data
  })
  uncorrected_PCA <- shiny::reactive({
    if (input$PCA_type == "PCA with clustering"){
      autoplot(
        cluster::pam(as.data.frame(t(uncorrected_numeric_matrix()[,-1])), input$cluster_number),
        #pam(as.data.frame(log2(t(uncorrected_numeric_matrix()[,-1] + 1))), input$cluster_number), 
        frame = TRUE,
        frame.type = 'norm')
    } else if (input$PCA_type == "PCA without clustering"){
      autoplot(
        stats::prcomp(as.data.frame(t(uncorrected_numeric_matrix()[,-1])), scale. = TRUE), 
        #prcomp(as.data.frame(log2(t(uncorrected_numeric_matrix()[,-1] + 1))), scale. = TRUE), 
        data = PCA_data(), 
        color = uncorrected_batch_choices_PCA2(),
        shape = uncorrected_biological_choices_PCA2()
      ) +
        ggplot2::scale_shape_manual(values = seq(0, length(aligned_meta_file()[,uncorrected_biological_choices_PCA2()])))
    }
  })
  output$uncorrected_PCA <- shiny::renderPlot({
    uncorrected_PCA()
  })
  
  
  #### PCA MC ------------------------------------------------------------------
  # Generating multiple PC plot
  output$PCA_mc_color_choice <- shiny::renderUI({
    shiny::selectInput("PCA_mc_color_choice", 
                       "Group color", 
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  uncorrected_PCA_mc_color_choice2 <- shiny::reactive({
    if (input$PCA_mc_color_choice == "Select"){
      NULL
    }else {
      uncorrected_PCA_mc_color_choice2 <- input$PCA_mc_color_choice
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
    uncorrected_PCA_mc_SCE <- scater::runPCA(uncorrected_PCA_mc_SCE, ncomponents = 50)
    scater::plotPCA(
      uncorrected_PCA_mc_SCE, 
      ncomponents = input$PCA_mc_slider, 
      colour_by = uncorrected_PCA_mc_color_choice2()
    )
  })
  output$uncorrected_PCA_multiple_components <- shiny::renderPlot({
    uncorrected_PCA_multiple_components()
  })
  
  
  #### PCA Details--------------------------------------------------------------
  # Uncorrected PCA details plots
  uncorrected_PCA_details <- shiny::reactive({
    FactoMineR::PCA(t(uncorrected_numeric_matrix()[,-1]))
  })
  uncorrected_PCA_details2 <- shiny::reactive({
    stats::prcomp(t(uncorrected_numeric_matrix()[,-1]))
  })
  uncorrected_scree_plot <- shiny::reactive({
    factoextra::fviz_eig(uncorrected_PCA_details(), addlabels = TRUE)
  })
  output$uncorrected_scree_plot <- shiny::renderPlot({
    uncorrected_scree_plot()
  })
  output$PCA_factors_choices <- shiny::renderUI({
    shiny::selectInput("PCA_factors_choices", 
                       "Select Factor for PCA Details", 
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))))
                       #c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  uncorrected_PCA_factors_choices2 <- shiny::reactive({
    if (input$PCA_factors_choices == "Select"){
      NULL
    }else {
      uncorrected_PCA_factors_choices2 <- input$PCA_factors_choices
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
    DT::datatable(uncorrected_PCA_individuals(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F) %>%
      formatRound(columns = c(3:ncol(uncorrected_PCA_individuals())), digits = 4)
    
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
    DT::datatable(uncorrected_contribution_counts(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
    
  })
  
  ## Generating heatmap of uncorrected data
  #uncorrected_heatmap <- shiny::reactive({
  #  uncorrected_matrix_heatmap <- as.matrix(uncorrected_numeric_matrix()[,-1])
  #  #uncorrected_matrix_heatmap <- as.matrix(log2(uncorrected_numeric_matrix()[,-1] + 1))
  #  uncorrected_matrix_scaled <- t(apply(uncorrected_matrix_heatmap, 1, scale))
  #  ComplexHeatmap::Heatmap(uncorrected_matrix_scaled)
  #})
  #output$uncorrected_heatmap <- shiny::renderPlot({
  #  uncorrected_heatmap()
  #})
  output$test <- DT::renderDataTable(
    print(meta_file_react()[input$batch_choices])
  )
  batch_names_from_meta <- shiny::reactive({
    batch_names_from_meta <- colnames(aligned_meta_file())
  })
  
  #### RLE ---------------------------------------------------------------------
  # Generating uncorrected RLE plot
  output$batch_choices_RLE <- shiny::renderUI({
    shiny::selectInput("batch_choices_RLE", 
                       "Group Color by Batch", 
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  uncorrected_batch_choices_RLE2 <- shiny::reactive({
    if (input$batch_choices_RLE == "Select"){
      NULL
    }else {
      uncorrected_batch_choices_RLE2 <- input$batch_choices_RLE
    }
  })
  uncorrected_RLE <- shiny::reactive({
    mat_RLE <- uncorrected_numeric_matrix()[,-1]
    names(mat_RLE) <- NULL
    uncorrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_RLE)),
      colData = aligned_meta_file(),
      rowData = uncorrected_numeric_matrix()[,1]
    )
    scater::plotRLE(
      uncorrected_RLE_SCE,
      exprs_values = "counts",
      color_by = uncorrected_batch_choices_RLE2()
    )
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
                       selected = c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       #c(unlist(strsplit(batch_names_from_meta(), ","))),
                       multiple = TRUE
    )
  })
  uncorrected_EV <- shiny::reactive({
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
    uncorrected_EV_SCE_PCA <- scater::runPCA(uncorrected_EV_SCE)
    scater::plotExplanatoryVariables(
      uncorrected_EV_SCE_PCA, 
      exprs_values = "logcounts", 
      variables = input$variable_choices_EV
    ) +
      ggplot2::scale_color_manual(values = my_colors)
  })
  output$uncorrected_EV_plot <- shiny::renderPlot({
    uncorrected_EV()
  })
  
  ### END OF UNCORRECTED PLOT DATA ---------------------------------------------
  
  
  ### START OF BATCH CORRECTION CRITERIA AND PLOT DATA -------------------------
  
  
  #### Batch Crtieria Selection ------------------------------------------------
  # Selection of batch criteria for correction. 
  output$batch1_selection_limma <- shiny::renderUI({
    shiny::selectInput("batch1_choices_limma", 
                       "Batch 1 Variable", 
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       selected = c(unlist(strsplit(batch_names_from_meta()[-1], ",")[1])))
                       #c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  output$batch2_selection_limma <- shiny::renderUI({
    shiny::selectInput("batch2_choices_limma", 
                       "Batch 2 Variable",
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       selected = c(unlist(strsplit(batch_names_from_meta()[-1], ",")[2])))
                       #c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  output$covariate_choices_limma <- shiny::renderUI({
    shiny::selectInput("covariate_choices_limma", 
                       "Select Any Covariates for Correction", 
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       #c("Select", unlist(strsplit(batch_names_from_meta(), ","))),
                       multiple = TRUE
    )
  })
  output$batch1_selection_ComBatseq <- shiny::renderUI({
    shiny::selectInput("batch1_choices_ComBatseq", 
                       "Select Variable for Batch", 
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))))
                       #c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  output$covariate_choices_ComBatseq <- shiny::renderUI({
    shiny::selectInput("covariate_choices_ComBatseq", 
                       "Select Biological Variables", 
                       c(unlist(strsplit(batch_names_from_meta()[-1], ","))),
                       #c("Select", unlist(strsplit(batch_names_from_meta(), ","))),
                       multiple = TRUE
    )
  })
  output$batch_selection_ComBat <- shiny::renderUI({
    shiny::selectInput(
      "batch_selection_ComBat",
      "Select Batch",
      c(unlist(strsplit(batch_names_from_meta()[-1], ",")))
      #c("Select", unlist(strsplit(batch_names_from_meta(), ",")))
    )
  })
  output$batch_selection_mean_centering <- shiny::renderUI({
    shiny::selectInput(
      "batch1_choices_mean_centering",
      "Select Batch for Mean Centering",
      c(unlist(strsplit(batch_names_from_meta()[-1], ",")))
      #c("Select", unlist(strsplit(batch_names_from_meta(), ",")))
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
  batch_selection_ComBat_input <- shiny::reactive({
    if (isTruthy(input$batch_selection_ComBat)) {
      if (input$batch_selection_ComBat == "Select"){
        NULL
      }else {
        batch_selection_Combat_input <- input$batch_selection_ComBat
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
        counter <- 0
        for (covariate in 1:length(input$covariate_choices_limma)){
          variable_object <- input$covariate_choices_limma[covariate]
          if (counter == 0){
            total_covariates <- paste(variable_object)
            counter <- counter + 1
          } else {
            total_covariates <- paste(total_covariates,"+", variable_object, sep = "")
            counter <- counter + 1
          }
        }
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
  
  #### Batch Correction ------------------------------------------------
  # The actual batch correction
  batch_correction <- shiny::reactive({
    if(input$batch_correction_method == "Limma"){
      if (isTruthy(batch1_choices()) & isTruthy(batch2_choices())) {
        limma::removeBatchEffect(
          uncorrected_numeric_matrix()[,-1],
          #log2(uncorrected_numeric_matrix()[,-1] + 1), 
          batch = c(unlist(aligned_meta_file()[,batch1_choices()])),
          batch2 = c(unlist(aligned_meta_file()[,batch2_choices()])),
          covariates = model_matrix()
        )
      }
    }else if(input$batch_correction_method == "ComBat"){
      if (isTruthy(batch_selection_ComBat_input())) {
        batch_combat <- c(unlist(aligned_meta_file()[,batch_selection_ComBat_input()]))
        modcombat <-  stats::model.matrix(~1, data = as.data.frame(aligned_meta_file()))
        ComBat_correction <- sva::ComBat(
          dat = uncorrected_numeric_matrix()[,-1],
          #dat = log2(uncorrected_numeric_matrix()[,-1] + 1), 
          batch = batch_combat, 
          mod = modcombat,
          par.prior = input$combat_parametric
        )
        ComBat_correction
      }
    }else if(input$batch_correction_method == "Mean Centering"){
      if (isTruthy(batch1_choices())) {
        mean_centering_batch = c(unlist(aligned_meta_file()[,batch1_choices()]))
        mean_centering_data = t(uncorrected_numeric_matrix()[,-1])
        mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
        mean_center_correction <- as.data.frame(t(mean_center$xadj))
        mean_center_correction
      }
    }else if(input$batch_correction_method == "ComBatseq"){
      if (isTruthy(batch1_choices())) {
        combatseq_corrected <- sva::ComBat_seq(
          as.matrix((2^uncorrected_numeric_matrix()[,-1])),
          batch = c(unlist(aligned_meta_file()[,batch1_choices()])),
          covar_mod = model_matrix()
        )
        log2(combatseq_corrected + 1)
        combatseq_corrected
      }
    }else if(input$batch_correction_method == "Select"){
      uncorrected_numeric_matrix()[,-1]
    }
  })
  
  # Test text for troubleshooting. Will delete in final product.
  output$test_print <- shiny::renderText({
    print(batch1_choices())
    #input$batch1_choices_mean_centering
    #c(unlist(aligned_meta_file()[,batch1_choices()]))
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
    if (input$batch_correction_method != "Select") {
      DT::datatable(corrected_numeric_matrix2(),
                    options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                   pageLength = 10,
                                   scrollX = T),
                    rownames = F)
    }
    
  })
  output$rendsave_corrected_matrix <- renderUI({
    if (input$batch_correction_method != "Select") {
      shiny::actionButton("save_corrected_matrix", "Save Corrected Matrix")
    }
  })
  
  #### PCA ------------------------------------------------
  # Generating the corrected PCA
  output$batch_choices_PCA <- shiny::renderUI({ 
    shiny::selectInput("batch_choices_PCA", 
                       "Group Color", 
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  corrected_batch_choices_PCA2 <- shiny::reactive({
    if (!is.null(input$batch_choices_PCA)) {
      if (input$batch_choices_PCA == "Select"){
        NULL
      }else {
        corrected_batch_choices_PCA2 <- input$batch_choices_PCA
      }
    }
    
  })
  output$biological_choices_PCA <- shiny::renderUI({ 
    shiny::selectInput("biological_choices_PCA", 
                       "Group Shape", 
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  corrected_biological_choices_PCA2 <- reactive({
    if (input$biological_choices_PCA == "Select"){
      NULL
    }else {
      corrected_biological_choices_PCA2 <- input$biological_choices_PCA
    }
  })
  PCA_data_corrected <- shiny::reactive({
    PCA_data_corrected <- cbind((t(batch_correction())), aligned_meta_file())
    PCA_data_corrected
  })
  corrected_PCA_react <- shiny::reactive({
    if (input$PCA_type == "PCA with clustering"){
      ggplot2::autoplot(
        cluster::pam(as.data.frame(t(batch_correction())),input$cluster_number), 
        frame = TRUE, 
        frame.type = 'norm')
    } else if (input$PCA_type == "PCA without clustering"){
      ggplot2::autoplot(
        stats::prcomp(as.data.frame(t(batch_correction())), scale. = TRUE), 
        data = PCA_data_corrected(), 
        color = corrected_batch_choices_PCA2(),
        shape = corrected_biological_choices_PCA2()
      ) +
        ggplot2::scale_shape_manual(values = seq(0, length(aligned_meta_file()[,corrected_biological_choices_PCA2()])))
    }
  })
  output$corrected_PCA <- shiny::renderPlot({
    corrected_PCA_react()
  })
  
  #Corrected PCA multiple components plot
  #output$corrected_PCA_mc_color_choice <- shiny::renderUI({
  #  shiny::selectInput("corrected_PCA_mc_color_choice", 
  #                     "Group color", 
  #                     c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  #})
  #### PCA MC ------------------------------------------------
  corrected_PCA_mc_color_choice2 <- shiny::reactive({
    if (input$PCA_mc_color_choice == "Select"){
      NULL
    }else {
      corrected_PCA_mc_color_choice2 <- input$PCA_mc_color_choice
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
    corrected_PCA_mc_SCE <- scater::runPCA(corrected_PCA_mc_SCE, ncomponents = 50)
    scater::plotPCA(
      corrected_PCA_mc_SCE, 
      ncomponents = input$PCA_mc_slider, 
      colour_by = corrected_PCA_mc_color_choice2()
    )
  })
  output$corrected_PCA_multiple_components <- shiny::renderPlot({
    corrected_PCA_multiple_components()
  })
  
  #### PCA Details ------------------------------------------------
  #Corrected PCA Details plots
  corrected_PCA_details <- shiny::reactive({
    FactoMineR::PCA(as.data.frame(t(batch_correction())))
  })
  corrected_PCA_details2 <- shiny::reactive({
    stats::prcomp(as.data.frame(t(batch_correction())))
  })
  corrected_scree_plot <- shiny::reactive({
    factoextra::fviz_eig(corrected_PCA_details(), addlabels = TRUE)
  })
  output$corrected_scree_plot <- shiny::renderPlot({
    corrected_scree_plot()
  })
  #output$PCA_factors_choices <- shiny::renderUI({
  #  shiny::selectInput("PCA_factors_choices", 
  #                     "Select Factor for PCA Details", 
  #                     c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  #})
  corrected_PCA_factors_choices2 <- shiny::reactive({
    if (input$PCA_factors_choices == "Select"){
      NULL
    }else {
      corrected_PCA_factors_choices2 <- input$PCA_factors_choices
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
    DT::datatable(corrected_PCA_individuals(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F) %>%
      formatRound(columns = c(3:ncol(corrected_PCA_individuals())), digits = 4)
    
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
    DT::datatable(corrected_contribution_counts(),
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
    
  })
  
  ## Generating the corrected heatmap
  #corrected_heatmap <- shiny::reactive({
  #  ComplexHeatmap::Heatmap(as.matrix(corrected_numeric_matrix2()[,-1]))
  #})
  #output$corrected_heatmap <- shiny::renderPlot({
  #  corrected_heatmap()
  #})
  
  # Generating corrected RLE plot
  #output$corrected_batch_choices_RLE <- shiny::renderUI({
  #  shiny::selectInput("corrected_batch_choices_RLE", 
  #                     "Group Color by Batch", 
  #                     c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  #})
  #### RLE ------------------------------------------------
  corrected_batch_choices_RLE2 <- shiny::reactive({
    if (input$batch_choices_RLE == "Select"){
      NULL
    }else {
      corrected_batch_choices_RLE2 <- input$batch_choices_RLE
    }
  })
  corrected_RLE <- shiny::reactive({
    mat_RLE <- batch_correction()
    names(mat_RLE) <- NULL
    corrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_RLE)),
      colData = aligned_meta_file(),
      rowData = corrected_numeric_matrix2()[,1]
    )
    scater::plotRLE(
      corrected_RLE_SCE,
      exprs_values = "counts",
      color_by = corrected_batch_choices_RLE2()
    )
  })
  output$corrected_RLE_plot <- shiny::renderPlot({
    corrected_RLE()
  })
  
  # Generating corrected EV plot
  #output$corrected_variable_choices_EV <- shiny::renderUI({
  #  shiny::selectInput("corrected_variable_choices_EV", 
  #                     "Select Variables to Plot", 
  #                     c(unlist(strsplit(batch_names_from_meta(), ","))),
  #                     multiple = TRUE
  #  )
  #})
  #### Exp Var ------------------------------------------------
  corrected_EV <- shiny::reactive({
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
    corrected_EV_SCE_PCA <- scater::runPCA(corrected_EV_SCE)
    scater::plotExplanatoryVariables(
      corrected_EV_SCE_PCA, 
      exprs_values = "logcounts", 
      variables = input$variable_choices_EV
    ) + 
      ggplot2::scale_color_manual(values = my_colors)
  })
  output$corrected_EV_plot <- shiny::renderPlot({
    corrected_EV()
  })
  ### END OF THE OUTPUT GENERATING CODE ----------------------------------------
  
  ### START OF THE SAVING FILES TO ZIP CODE ------------------------------------
  
  # Saving files to ZIP
  temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
  dir.create(temp_directory)
  shiny::observeEvent(input$save_uncorrected_matrix, {
    file_name <- paste("uncorrected_matrix", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(uncorrected_matrix(), file.path(temp_directory, file_name))
    print("Uncorrected_matrix Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_meta_file, {
    file_name <- paste("uncorrected_meta_file", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(aligned_meta_file(), file.path(temp_directory, file_name))
    print("Uncorrected_meta_file Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_PCA_plot, {
    file_name <- paste("uncorrected_PCA_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_PCA())
    print("Uncorrected_PCA_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_PCA_mc_plot, {
    file_name <- paste("uncorrected_PCA_mc_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_PCA_multiple_components())
    print("Uncorrected_PCA_mc_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_scree_plot, {
    file_name <- paste("uncorrected_Scree_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_scree_plot())
    print("Uncorrected_Scree_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_PCA_components, {
    file_name <- paste("uncorrected_PC_components", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(as.data.frame(uncorrected_PCA_details2()$x), file.path(temp_directory, file_name))
    print("uncorrected_PCA_details Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_contribution_table, {
    file_name <- paste("uncorrected_contribution_table", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(uncorrected_PCA_individuals(), file.path(temp_directory, file_name))
    print("uncorrected_contribution_table Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_contribution_counts, {
    file_name <- paste("uncorrected_contribution_counts", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(uncorrected_contribution_counts(), file.path(temp_directory, file_name))
    print("uncorrected_contribution_counts Ready for Zip")
  })
  #shiny::observeEvent(input$save_uncorrected_heatmap, {
  #  file_name <- paste(temp_directory, "/", "uncorrected_heatmap", Sys.Date(),".png", sep = "")
  #  png(filename = file_name)
  #  ComplexHeatmap::draw(uncorrected_heatmap())
  #  dev.off()
  #  print("Uncorrected_heatmap Ready for Zip")
  #})
  shiny::observeEvent(input$save_uncorrected_RLE_plot, {
    file_name <- paste("uncorrected_RLE_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_RLE())
    print("Uncorrected_RLE_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_uncorrected_EPC_plot, {
    file_name <- paste("uncorrected_EV_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_EV())
    print("Uncorrected_EV_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_matrix, {
    file_name <- paste(input$batch_correction_method, "_corrected_matrix", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(corrected_numeric_matrix2(), file.path(temp_directory, file_name))
    print("Corrected_matrix Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_PCA_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_PCA_plot", Sys.Date(),".png", sep = "")
    ggsave(filename = file_name, path = temp_directory, plot = corrected_PCA())
    print("Corrected_PCA_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_PCA_mc_plot, {
    file_name <- paste("corrected_PCA_mc_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_PCA_multiple_components())
    print("Corrected_PCA_mc_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_scree_plot, {
    file_name <- paste("corrected_Scree_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_scree_plot())
    print("Corrected_Scree_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_PCA_components, {
    file_name <- paste(input$batch_correction_method, "_corrected_PC_components", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(as.data.frame(corrected_PCA_details2()$x), file.path(temp_directory, file_name))
    print("Corrected_PCA_details Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_contribution_table, {
    file_name <- paste("corrected_contribution_table", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(corrected_PCA_individuals(), file.path(temp_directory, file_name))
    print("corrected_contribution_table Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_contribution_counts, {
    file_name <- paste("corrected_contribution_counts", Sys.Date(),".tsv", sep = "")
    readr::write_tsv(corrected_contribution_counts(), file.path(temp_directory, file_name))
    print("corrected_contribution_counts Ready for Zip")
  })
  #shiny::observeEvent(input$save_corrected_heatmap, {
  #  file_name <- paste(temp_directory,"/",input$batch_correction_method,  "_corrected_heatmap", Sys.Date(),".png", sep = "")
  #  png(filename = file_name)
  #  ComplexHeatmap::draw(corrected_heatmap())
  #  dev.off()
  #  print("Corrected_heatmap Ready for Zip")
  #})
  shiny::observeEvent(input$save_corrected_RLE_plot, {
    file_name <- paste(input$batch_correction_method, "_corrected_RLE_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_RLE())
    print("Corrected_RLE_plot Ready for Zip")
  })
  shiny::observeEvent(input$save_corrected_EV_plot, {
    file_name <- paste(input$batch_correction_method, "corrected_EV_plot", Sys.Date(),".png", sep = "")
    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_EV())
    print("Corrected_EV_plot Ready for Zip")
  })
  output$download_btn <- shiny::downloadHandler(
    filename = function(){
      paste(input$file_name, "_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file){
      zip::zip(
        zipfile = file,
        files = dir(temp_directory),
        root = temp_directory
      )
    },
    contentType = "application/zip"
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
