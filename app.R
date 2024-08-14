
# Back-End Input Files ---------------------------------------------------------

# Example Data Files
Example_Matrix_File <- "Example_Data/Example_Immgen_Microarray_Expression.txt"
Example_Meta_File <- "Example_Data/Example_Immgen_Microarray_Meta.txt"
Example_MatrixCounts_File <- "Example_Data/Example_Thapsigargin_treated_RNASeqCounts.txt"
Example_MetaCounts_File <- "Example_Data/Example_Thapsigargin_treated_Meta.txt"

# Homepage Files
#homepage_txt <- "Batchflex_Homepage/Intro_placeholder.txt"
homepage_filepath <- "Batchflex_Homepage/text_files/"
homepage_tutorial_list <- list.files(homepage_filepath)

# Housekeeping gene libraries
Eisenberg_hkg_File <- "default_housekeeping/eisenberg.txt"
Lin500_hkg_File  <- "default_housekeeping/lin500.txt"
HSIAO_hkg_File  <- "default_housekeeping/hsiao.txt"
Eisenberg_hkg_mm_File <- "default_housekeeping/eisenberg_mouse.txt"
Lin500_hkg_mm_File  <- "default_housekeeping/lin500_mouse.txt"
HSIAO_hkg_mm_File  <- "default_housekeeping/hsiao_mouse.txt"

# Gene set data
GeneSet_RData_File <- "GeneSet_Data/GeneSet_List_HS_v6.RData"
GeneSetTableIn_File <- "GeneSet_Data/GeneSet_CatTable_v6.txt"

# Gene Conversion Table
MM_HS_Conv_File <- "Mouse_To_Human_Gene_Conversion/hs_mm_homo_updated_v3_20231213.txt"


# Read in back end data --------------------------------------------------------

# Homepage Files
#homepage_text <- readtext::readtext(homepage_txt)
homepage_tutorial_text_list <- list()
for (file in homepage_tutorial_list){
  filename <- gsub("\\..*", "", file)
  homepage_tutorial_text_list[[filename]] <- readtext::readtext(paste0(homepage_filepath,file))
}

## Housekeeping genes
Eisenberg_hkg <- read.delim(Eisenberg_hkg_File, header = F, sep = '\t')
Lin500_hkg <- read.delim(Lin500_hkg_File, header = F, sep = '\t')
HSIAO_hkg <- read.delim(HSIAO_hkg_File, header = F, sep = '\t')
Eisenberg_hkg_mm <- read.delim(Eisenberg_hkg_mm_File, header = F, sep = '\t')
Lin500_hkg_mm <- read.delim(Lin500_hkg_mm_File, header = F, sep = '\t')
HSIAO_hkg_mm <- read.delim(HSIAO_hkg_mm_File, header = F, sep = '\t')

## Gene set data
geneset <- loadRData(GeneSet_RData_File)
geneset_df <- as.data.frame(read_delim(GeneSetTableIn_File, delim = '\t', col_names = T))
### Remove ER and TCGA genesets
geneset_df <- geneset_df[which(!geneset_df[,1] %in% c("TCGA","ER Stress")),]

# Mouse to human conversion
MM_HS_Conv <- as.data.frame(read_delim(MM_HS_Conv_File,delim = '\t', col_names = F))
colnames(MM_HS_Conv) <- c("Human","Mouse")

# Just moved this from near the save observations in the server as temp_directory was needed earlier in the code
temp_directory <- file.path(tempdir(), UUIDgenerate())
#temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
dir.create(temp_directory)

umap.mod <- umap.defaults


# UI ---------------------------------------------------------------------------
ui <-
  shiny::navbarPage("{ Batch-Flex }",
                    id = "navbar_id",
                    theme = bslib::bs_theme(bootswatch = "yeti"),
                    ## Homepage Tab --------------------------------------------
                    shiny::tabPanel("Homepage",
                                    bslib::page_fillable(
                                      bslib::layout_column_wrap(
                                        height = 900,
                                        bslib::layout_columns(
                                          bslib::card(bslib::card_header("Introduction"),
                                                      bslib::card_body(shiny::htmlOutput("homepage_text"))
                                          ),
                                          bslib::card(bslib::card_header("BatchFLEX Workflow"),
                                                      bslib::card_body(shiny::imageOutput("homepage"))
                                          ),
                                          col_widths = c(4,8)
                                        ),
                                        bslib::layout_column_wrap(
                                          bslib::card(bslib::card_header(shiny::HTML("BatchFLEX<br/>Shiny App")),
                                                      bslib::card_body(shiny::tags$a(shiny::imageOutput("shiny_app", width = 180, height = 200), href="https://shawlab-moffitt.shinyapps.io/batchflex_v2/", target = "_blank")),
                                                      shiny::actionButton("get_started", "Get Started", icon = tags$i(class = "fa-solid fa-chart-simple")),
                                                      shiny::actionButton("need_tutorial", "Need Tutorial?", icon = tags$i(class = "fa-solid fa-lines-leaning"))
                                          ),
                                          bslib::card(bslib::card_header(shiny::HTML("BatchFLEX<br/>Shiny Github")),
                                                      bslib::card_body(shiny::tags$a(shiny::imageOutput("shiny_github", width = 180, height = 200), href="https://github.com/shawlab-moffitt/BATCH-FLEX-ShinyApp", target = "_blank")),
                                                      shiny::actionButton("shiny_github_button", "Link to GitHub", icon = tags$i(class = "fa-brands fa-github"),
                                                                          onclick ="window.open('https://github.com/shawlab-moffitt/BATCH-FLEX-ShinyApp', '_blank')")
                                          ),
                                          bslib::card(bslib::card_header(shiny::HTML("BatchFLEX<br/>Function Github")),
                                                      bslib::card_body(shiny::tags$a(shiny::imageOutput("function_github", width = 180, height = 200), href="https://github.com/shawlab-moffitt/BATCHFLEX", target = "_blank")),
                                                      shiny::actionButton("function_github_button", "Link to GitHub", icon = tags$i(class = "fa-brands fa-github"),
                                                                          onclick ="window.open('https://github.com/shawlab-moffitt/BATCHFLEX', '_blank')")
                                          ),
                                          #bslib::card(bslib::card_header(shiny::HTML("Journal<br/>Link")),
                                          #            #tags$script("Shiny.addCustomMessageHandler('txt', function (txt) {navigator.clipboard.writeText(txt);});"),
                                          #            bslib::card_body(shiny::tags$a(shiny::imageOutput("journal", width = 180, height = 200), href="https://academic.oup.com/bioinformatics", target = "_blank")),
                                          #            shiny::actionButton("citation", "Cite", icon = tags$i(class = "fa-solid fa-arrow-up-right-from-square"),
                                          #                                onclick ="window.open('https://academic.oup.com/bioinformatics', '_blank')")
                                          #),
                                          bslib::card(bslib::card_header(shiny::HTML("MergeQC<br/>Shiny App")),
                                                      bslib::card_body(shiny::tags$a(shiny::imageOutput("shiny_mergeqc", width = 180, height = 200), href="https://shawlab-moffitt.shinyapps.io/mergeqc/", target = "_blank")),
                                                      shiny::actionButton("shiny_mergeqc_button", "Link to App", icon = tags$i(class = "fa-solid fa-arrow-up-right-from-square"),
                                                                          onclick ="window.open('https://shawlab-moffitt.shinyapps.io/mergeqc/', '_blank')")
                                          ),
                                          bslib::card(bslib::card_header(shiny::HTML("MergeQC<br/> Github")),
                                                      bslib::card_body(shiny::tags$a(shiny::imageOutput("mergeqc_github", width = 180, height = 200), href="https://github.com/shawlab-moffitt/mergeQC", target = "_blank")),
                                                      shiny::actionButton("mergeqc_github_button", "Link to GitHub", icon = tags$i(class = "fa-brands fa-github"),
                                                                          onclick ="window.open('https://github.com/shawlab-moffitt/mergeQC', '_blank')")
                                          )
                                        )
                                      )
                                    )
                    ),
                    ## Data Input Tab ------------------------------------------
                    shiny::tabPanel("Step 1 - Data Input",
                                    shiny::sidebarLayout(
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
                                          shiny::column(5, style = "margin-top:-15px",
                                                        checkboxInput("RawCountInput","Input is RNAseq Count", value = F),
                                                        #conditionalPanel(condition = "input.RawCountInput == true",
                                                        #                 div(selectInput("RawCountNorm","Normalize RNAseq count Counts By:",
                                                        #                             c("No Normalization" = "none","TMM" = "TMM","Upper Quartile" = "upperquartile"), selected = "TMM"),
                                                        #                     style = "margin-top:-15px")
                                                        #                 )
                                          ),
                                          shiny::column(2, style = "margin-top:-15px",
                                                        conditionalPanel(condition = "input.RawCountInput == false",
                                                                         shiny::checkboxInput("Log_Choice","Log2+1", value = T)
                                                        )
                                          ),
                                          shiny::column(2, style = "margin-top:-15px",
                                                        numericInput("SeedSet","Set Seed:", value = 101, min = 1, step = 1)
                                          )
                                          #shiny::column(3, style = "margin-top:15px",
                                          #              actionButton("UseExpData","Load Example")
                                          #)
                                        ),
                                        shiny::fluidRow(
                                          column(5, style = "padding-right:0px",
                                                 actionButton("UseExpData","Load Micro Array Example")
                                                 ),
                                          column(5, style = "padding-right:0px;padding-left:0px",
                                                 actionButton("UseRawExpData","Load RNASeq Counts Example")
                                                 ),
                                          column(2, style = "padding-left:0px",
                                                 tags$a(href="https://github.com/shawlab-moffitt/BATCH-FLEX-ShinyApp/tree/main/Example_Data", "Download example", target='_blank'))
                                        ),
                                        div(hr(),style = "margin-bottom:-10px"),
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
                                        shiny::uiOutput("rendbatch_names"),
                                        # if batch information is provided by user
                                        fluidRow(
                                          column(12, style = "margin-top:-15px",
                                                 shiny::uiOutput("rendInit_Batch_Select")
                                          )
                                        )
                                      ),
                                      shiny::mainPanel(
                                        uiOutput("rendDataLodedHelpText"),
                                        shiny::uiOutput("rendMatrixHeader"),
                                        uiOutput("rendExprHead"),
                                        DT::dataTableOutput("uncorrected_matrix_output_input"),
                                        shiny::uiOutput("rendMetaHeader"),
                                        uiOutput("rendClinHead"),
                                        DT::dataTableOutput("meta_file")
                                      )
                                    ),
                                    value = "step_1"
                    ),
                    ## Batch Correction Tab ------------------------------------
                    shiny::tabPanel("Step 2 - Batch Correction",
                                    # Side Panel
                                    shiny::sidebarLayout(
                                      ### Sidebar ------------------------------
                                      shiny::sidebarPanel(
                                        width = 3,
                                        shiny::tabsetPanel(
                                          shiny::tabPanel("Data Input",
                                                          p(),
                                                          #### Batch Crit Input ----------------------------
                                                          #conditionalPanel(condition = "input.uncorrected_panel == 'mat' & input.RawMatTabs == '1'",
                                                          #                 selectInput("RNAcountsMethod","RNASeq Counts Correction Method:",
                                                          #                             choices = c("ComBatseq","RUVg")),
                                                          #                 conditionalPanel(condition = "input.RNAcountsMethod == 'ComBatseq'",
                                                          #                                  uiOutput("rendbatch1_choices_ComBatseq_raw"),
                                                          #                                  uiOutput("rendcovariate_choices_ComBatseq_raw")
                                                          #                                  ),
                                                          #                 conditionalPanel(condition = "input.RNAcountsMethod == 'RUVg'",
                                                          #                                  selectInput(
                                                          #                                    "RUVg_housekeeping_selection1",
                                                          #                                    "Select Housekeeping Genes for Control",
                                                          #                                    c("Eisenberg", "Lin500", "HSIAO", "UserInput")),
                                                          #                                  conditionalPanel(condition = "input.RUVg_housekeeping_selection1 == 'UserInput'",
                                                          #                                                   fileInput("RUVg_user_control_genes1", "Please Provide List of Control Genes")),
                                                          #                                  fluidRow(
                                                          #                                    column(6, style = 'padding-right:0px;',
                                                          #                                           numericInput("RUVg_estimate_factors1","Factors to Estimate",
                                                          #                                                        value = 2,min = 1,max = 10,step = 1)
                                                          #                                    ),
                                                          #                                    column(6,
                                                          #                                           numericInput("RUVg_drop_factors1","Factors to Drop",
                                                          #                                                        value = 0,min = 0,max = 9,step = 1)
                                                          #                                    )
                                                          #                                  ),
                                                          #                                  fluidRow(
                                                          #                                    column(6,
                                                          #                                           selectInput("RUVg_tolerance1","RUVg Tolerance",
                                                          #                                                       c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10),selected = 1e-8)
                                                          #                                    ),
                                                          #                                    column(6,
                                                          #                                           checkboxInput("RUVg_mean_centered1","Mean Centered",value = TRUE)
                                                          #                                    )
                                                          #                                  )
                                                          #                                  )
                                                          #),
                                                          #conditionalPanel(condition = "input.RawMatTabs != '1'",
                                                                           #conditionalPanel(condition = "input.uncorrected_panel == 'mat' & input.RawMatTabs != '1'",
                                                                           bslib::accordion(
                                                                             bslib::accordion_panel(id = "CorrMethod",
                                                                                                    "Correction Method",
                                                                                                    h3("Method 1:"),
                                                                                                    shiny::fluidRow(
                                                                                                      shiny::column(8,
                                                                                                                    uiOutput("rendbatch_correction_method")
                                                                                                      ),
                                                                                                      shiny::column(4, style = "margin-top:10px",
                                                                                                                    uiOutput("rendQuantNorm")
                                                                                                      )
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method == 'Limma'",
                                                                                                                            shiny::fluidRow(
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::uiOutput("batch1_selection_limma")
                                                                                                                              ),
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::uiOutput("batch2_selection_limma")
                                                                                                                              )
                                                                                                                            ),
                                                                                                                            shiny::uiOutput("covariate_choices_limma")
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method == 'ComBat'",
                                                                                                                            shiny::fluidRow(
                                                                                                                              shiny::column(8,
                                                                                                                                            shiny::uiOutput("batch_selection_ComBat")
                                                                                                                              ),
                                                                                                                              shiny::column(4, style = "margin-top:20px",
                                                                                                                                            shiny::checkboxInput("combat_parametric", "Parametric?", value = TRUE)
                                                                                                                              )
                                                                                                                            )
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method == 'Mean Centering'",
                                                                                                                            shiny::uiOutput("batch_selection_mean_centering")
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method == 'ComBatseq'",
                                                                                                                            shiny::uiOutput("batch1_selection_ComBatseq"),
                                                                                                                            shiny::uiOutput("covariate_choices_ComBatseq")
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method == 'Harman'",
                                                                                                                            shiny::uiOutput("batch_selection_harman"),
                                                                                                                            shiny::fluidRow(
                                                                                                                              shiny::column(6,
                                                                                                                                            numericInput("HarmanCL1","Confidence Limit:", value = 0.95, min = 0, max = 1, step = 0.05)
                                                                                                                              ),
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::uiOutput("treatment_selection_harman")
                                                                                                                              )
                                                                                                                            )
                                                                                                    ),
                                                                                                    conditionalPanel(condition = "input.batch_correction_method == 'RUVg'",
                                                                                                                     selectInput(
                                                                                                                       "RUVg_housekeeping_selection",
                                                                                                                       "Select Housekeeping Genes for Control",
                                                                                                                       c("Eisenberg", "Lin500", "HSIAO", "UserInput")),
                                                                                                                     conditionalPanel(condition = "input.RUVg_housekeeping_selection == 'UserInput'",
                                                                                                                                      fileInput("RUVg_user_control_genes", "Please Provide List of Control Genes")),
                                                                                                                     fluidRow(
                                                                                                                       column(6, style = 'padding-right:0px;',
                                                                                                                              numericInput("RUVg_estimate_factors","Factors to Estimate",
                                                                                                                                           value = 2,min = 1,max = 10,step = 1)
                                                                                                                       ),
                                                                                                                       column(6,
                                                                                                                              numericInput("RUVg_drop_factors","Factors to Drop",
                                                                                                                                           value = 0,min = 0,max = 9,step = 1, width = "95%")
                                                                                                                       )
                                                                                                                       #column(4, style = 'padding-left:0px',
                                                                                                                       #       selectInput("RUVg_tolerance","RUVg Tolerance",
                                                                                                                       #                   c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10),selected = 1e-8)
                                                                                                                       #)
                                                                                                                     ),
                                                                                                                     fluidRow(
                                                                                                                       column(6,
                                                                                                                              selectInput("RUVg_tolerance","RUVg Tolerance",
                                                                                                                                          c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10),selected = 1e-8)
                                                                                                                       ),
                                                                                                                       column(6,
                                                                                                                              checkboxInput("RUVg_mean_centered","Mean Centered",value = TRUE)
                                                                                                                       )
                                                                                                                     )
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method == 'SVA'",
                                                                                                                            shiny::fluidRow(
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::selectInput("svaMethod_bc","Method",c("be","leek"))
                                                                                                                              ),
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::uiOutput("SVA_variable_of_interest_bc")
                                                                                                                              )
                                                                                                                            )
                                                                                                    ),
                                                                                                    div(hr(),style = "margin-top:-20px;margin-bottom:-10px"),
                                                                                                    h3("Method 2:"),
                                                                                                    shiny::fluidRow(
                                                                                                      shiny::column(8,
                                                                                                                    uiOutput("rendbatch_correction_method2")
                                                                                                      ),
                                                                                                      shiny::column(4, style = "margin-top:10px",
                                                                                                                    uiOutput("rendQuantNorm2")
                                                                                                      )
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method2 == 'Limma'",
                                                                                                                            shiny::fluidRow(
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::uiOutput("batch1_selection_limma2")
                                                                                                                              ),
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::uiOutput("batch2_selection_limma2")
                                                                                                                              )
                                                                                                                            ),
                                                                                                                            shiny::uiOutput("covariate_choices_limma2")
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method2 == 'ComBat'",
                                                                                                                            shiny::fluidRow(
                                                                                                                              shiny::column(8,
                                                                                                                                            shiny::uiOutput("batch_selection_ComBat2")
                                                                                                                              ),
                                                                                                                              shiny::column(4, style = "margin-top:20px",
                                                                                                                                            shiny::checkboxInput("combat_parametric2", "Parametric?", value = TRUE)
                                                                                                                              )
                                                                                                                            )
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method2 == 'Mean Centering'",
                                                                                                                            shiny::uiOutput("batch_selection_mean_centering2")
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method2 == 'ComBatseq'",
                                                                                                                            shiny::uiOutput("batch1_selection_ComBatseq2"),
                                                                                                                            shiny::uiOutput("covariate_choices_ComBatseq2")
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method2 == 'Harman'",
                                                                                                                            shiny::uiOutput("batch_selection_harman2"),
                                                                                                                            shiny::column(6,
                                                                                                                                          numericInput("HarmanCL2","Confidence Limit:", value = 0.95, min = 0, max = 1, step = 0.05)
                                                                                                                            ),
                                                                                                                            shiny::column(6,
                                                                                                                                          shiny::uiOutput("treatment_selection_harman2")
                                                                                                                            )
                                                                                                    ),
                                                                                                    conditionalPanel(condition = "input.batch_correction_method2 == 'RUVg'",
                                                                                                                     selectInput(
                                                                                                                       "RUVg_housekeeping_selection2",
                                                                                                                       "Select Housekeeping Genes for Control",
                                                                                                                       c("Eisenberg", "Lin500", "HSIAO", "UserInput")),
                                                                                                                     conditionalPanel(condition = "input.RUVg_housekeeping_selection == 'UserInput'",
                                                                                                                                      fileInput("RUVg_user_control_genes", "Please Provide List of Control Genes")),
                                                                                                                     fluidRow(
                                                                                                                       column(6, style = 'padding-right:0px;',
                                                                                                                              numericInput("RUVg_estimate_factors2","Factors to Estimate",
                                                                                                                                           value = 2,min = 1,max = 10,step = 1)
                                                                                                                       ),
                                                                                                                       column(6,
                                                                                                                              numericInput("RUVg_drop_factors2","Factors to Drop",
                                                                                                                                           value = 0,min = 0,max = 9,step = 1, width = "95%")
                                                                                                                       )
                                                                                                                     ),
                                                                                                                     fluidRow(
                                                                                                                       column(6,
                                                                                                                              selectInput("RUVg_tolerance2","RUVg Tolerance",
                                                                                                                                          c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10),selected = 1e-8)
                                                                                                                       ),
                                                                                                                       column(6,
                                                                                                                              checkboxInput("RUVg_mean_centered2","Mean Centered",value = TRUE)
                                                                                                                       )
                                                                                                                     )
                                                                                                    ),
                                                                                                    shiny::conditionalPanel(condition = "input.batch_correction_method2 == 'SVA'",
                                                                                                                            shiny::fluidRow(
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::selectInput("svaMethod_bc2","Method",c("be","leek"))
                                                                                                                              ),
                                                                                                                              shiny::column(6,
                                                                                                                                            shiny::uiOutput("SVA_variable_of_interest_bc2")
                                                                                                                              )
                                                                                                                            )
                                                                                                    ),
                                                                                                    conditionalPanel(condition = "input.RawCountInput == true",
                                                                                                                     div(hr(),style = "margin-top:-10px;margin-bottom:-10px"),
                                                                                                                     h4("RNAseq Count Normalization"),
                                                                                                                     fluidRow(
                                                                                                                       column(9,
                                                                                                                              selectInput("RawCountNorm","Normalize RNAseq Counts By:",
                                                                                                                                          c("TMM" = "TMM","Upper Quartile" = "upperquartile"),
                                                                                                                                          selected = "TMM")
                                                                                                                       ),
                                                                                                                       column(3,
                                                                                                                              div(shiny::checkboxInput("Log_Choice_Raw","Log2+1", value = T),
                                                                                                                                  style = "margin-top:30px")
                                                                                                                       )
                                                                                                                     ))
                                                                             )
                                                                           ),
                                                          #),
                                                          #### Plot Inputs ------------------------------------
                                                          p(),
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca'",
                                                                                  shiny::h3("PCA Parameters"),
                                                                                  shiny::radioButtons("PCA_type",NULL,c("Meta Annotation","Cluster Annotation"),inline = T),
                                                                                  shiny::conditionalPanel(condition = "input.PCA_type == 'Meta Annotation'",
                                                                                                          shiny::uiOutput("rendbatch_choices_PCA"),
                                                                                                          shiny::uiOutput("rendPCAhover")
                                                                                  ),
                                                                                  shiny::conditionalPanel(condition = "input.PCA_type == 'Cluster Annotation'",
                                                                                                          shiny::fluidRow(
                                                                                                            shiny::column(6,
                                                                                                                          shiny::numericInput("cluster_number","Number of Clusters",
                                                                                                                                              value = 2,step = 1,min = 1)
                                                                                                            ),
                                                                                                            shiny::column(6, style = "margin-top:15px",
                                                                                                                          shiny::checkboxInput("FrameClusters","Frame Clusters",value = T)
                                                                                                            )
                                                                                                          )
                                                                                  )
                                                          ),
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
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca_dt'",
                                                                                  shiny::h3("PCA Detail Parameters"),
                                                                                  shiny::uiOutput("PCA_factors_choices")
                                                          ),
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
                                                                                                          div(DT::dataTableOutput("GeneSetTableUIUMAP"), style = "font-size:10px")
                                                                                  )
                                                          ),
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'cluster_main'",
                                                                                  shiny::h3("Cluster Parameters"),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  numericInput("cluster_n_MV_features","Top Most Variable Features:",
                                                                                                               value = 2000)
                                                                                    ),
                                                                                    shiny::column(6,
                                                                                                  selectInput("VarianceMeasure", "Select Variance Measure:",
                                                                                                              choices = c("MAD","CV","VAR"))
                                                                                    )
                                                                                  )
                                                          ),
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'cluster_main' & input.Cluster_main_pan == 'heatmap'",
                                                                                  div(shiny::h3("Heatmap Parameters"), style = "margin-top:-20px"),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  selectInput("ClusterMethodHeat","Cluster Method:",
                                                                                                              choices = c("ward.D2","ward.D", "complete","single","average","mcquitty","median","centroid"))
                                                                                    ),
                                                                                    shiny::column(6,
                                                                                                  uiOutput("rendHeatmapAnnoSel")
                                                                                    )
                                                                                  )
                                                          ),
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'cluster_main' & input.Cluster_main_pan == 'clinfo'",
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  selectInput("ClusterMethod","Cluster Method:",
                                                                                                              choices = c("ward.D2","ward.D", "complete","single","average","mcquitty","median","centroid")),
                                                                                                  shiny::uiOutput("rendBarPFillCol")
                                                                                    ),
                                                                                    shiny::column(6,
                                                                                                  numericInput("NumOfCluster","Number of Clusters:",value = 3,min = 1,step = 1),
                                                                                                  div(shiny::checkboxInput("barPfill","Fill Bars as Percentage",value = F),style = "margin-top:40px")
                                                                                    )
                                                                                  ),
                                                          ),
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'rle'",
                                                                                  shiny::h3("RLE Parameters"),
                                                                                  shiny::uiOutput("batch_choices_RLE")
                                                          ),
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
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'sva'",
                                                                                  shiny::h3("SVA Parameters"),
                                                                                  shiny::uiOutput("uncorrected_SVA_variable_of_interest"),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  shiny::selectInput("svaMethod","Method",c("leek","be"))
                                                                                    ),
                                                                                    shiny::column(6,
                                                                                                  shiny::numericInput("SVAvarNum","Number of Variables",
                                                                                                                      value = 300,min = 0, step = 1)
                                                                                    )
                                                                                  ),
                                                                                  shiny::fluidRow(
                                                                                    shiny::column(6,
                                                                                                  shiny::numericInput("NSV1","Update Estimated Num of Surrogate Variables Matrix A",
                                                                                                                      value = 2,min = 2, step = 1)
                                                                                    ),
                                                                                    shiny::column(6,
                                                                                                  shiny::numericInput("NSV2","Update Estimated Num of Surrogate Variables Matrix B",
                                                                                                                      value = 2,min = 2, step = 1)
                                                                                    )
                                                                                  ),
                                                                                  selectizeInput("SVAcolor","Color Plot By:",choices = NULL),
                                                                                  uiOutput("rendSVAhover"),
                                                                                  actionButton("CalcSurVar","Calculate Number of Estimated Surrogate Variables"),
                                                                                  uiOutput("rendSurVarEstimates")
                                                                                  #shiny::h4("Download Surrogate Variables"),
                                                                                  #shiny::fluidRow(
                                                                                  #  column(6,
                                                                                  #         actionButton("save_SVA_surrogate_variables", "Add to Zip File Export")
                                                                                  #  ),
                                                                                  #  column(6,
                                                                                  #         shiny::downloadButton("dnldsave_SVA_surrogate_variables","Dowload Meta with SVA")
                                                                                  #  )
                                                                                  #)
                                                          ),
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'box'",
                                                                                  #style = "overflow-y:scroll;overflow-x:hidden;max-height:70vh;position:relative;",
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
                                                          shiny::conditionalPanel(condition = "input.uncorrected_panel == 'pca_main' & input.PCA_main_pan == 'pca'",
                                                                                  p(),
                                                                                  h4("PCA Figure Parameters"),
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           numericInput("pcaFontSize","Font Size:",
                                                                                                        value = 14, step = 1)
                                                                                    ),
                                                                                    column(6,
                                                                                           numericInput("PCAdotSize","Dot Size",value = 7)
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
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'cluster_main' & input.Cluster_main_pan == 'clinfo'",
                                                                                  p(),
                                                                                  h3("Bar Plot Figure Parameters"),
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           numericInput("barPAxisTkSize","Axis Tick Font:",
                                                                                                        value = 12, step = 1)
                                                                                    ),
                                                                                    column(6,
                                                                                           numericInput("barPAxisTtSize","Axis Title Font:",
                                                                                                        value = 14, step = 1)
                                                                                    )
                                                                                  ),
                                                                                  h4("Figure Download Parameters"),
                                                                                  fluidRow(
                                                                                    column(4,
                                                                                           numericInput("barPHeight","Height",value = 8)
                                                                                    ),
                                                                                    column(4,
                                                                                           numericInput("barPWidth","Width",value = 10)
                                                                                    ),
                                                                                    column(4,
                                                                                           selectInput("barPUnits","Units",choices = c("in","cm","mm","px"))
                                                                                    )
                                                                                  )
                                                          ),
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
                                                          shiny::conditionalPanel("input.uncorrected_panel == 'sva'",
                                                                                  p(),
                                                                                  h4("SVA Scatter Plot Figure Parameters"),
                                                                                  fluidRow(
                                                                                    column(6,
                                                                                           numericInput("SVAfontSize","Font Size:",
                                                                                                        value = 14, step = 1)
                                                                                    ),
                                                                                    column(6,
                                                                                           numericInput("SVAdotSize","Dot Size",value = 7)
                                                                                    )
                                                                                  ),
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
                                                                                                       choices = c(0,45,90))
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
                                          )
                                        )
                                      ),
                                      ### Main ---------------------------------
                                      mainPanel(
                                        shiny::tabsetPanel(
                                          id = "uncorrected_panel",
                                          #### Matrix ---------------------------------
                                          shiny::tabPanel(id = "MatPanel",
                                                          "Matrix",
                                                          shiny::p(),
                                                          shiny::uiOutput("rendMatrixTabs"),
                                                          value = "mat"
                                          ),
                                          #### PCA -------------------------------------
                                          shiny::tabPanel(
                                            "PCA",
                                            p(),
                                            tabsetPanel(
                                              id = "PCA_main_pan",
                                              shiny::tabPanel(
                                                "PCA Plot",
                                                p(),
                                                shiny::fluidRow(
                                                  shiny::column(6,
                                                                p(),
                                                                shiny::uiOutput("PCA1title"),
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
                                                                p(),
                                                                shiny::uiOutput("PCA2title"),
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
                                              shiny::tabPanel(
                                                "Multiple Components PCA",
                                                shiny::p(),
                                                shiny::fluidRow(
                                                  shiny::column(6,
                                                                p(),
                                                                shiny::uiOutput("PCAMC1title"),
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
                                                                p(),
                                                                shiny::uiOutput("PCAMC2title"),
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
                                              shiny:: tabPanel(
                                                "PCA Details",
                                                shiny::p(),
                                                shiny::fluidRow(
                                                  shiny::column(6,
                                                                p(),
                                                                shiny::uiOutput("PCADet1title"),
                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_scree_plot")), type = 6),
                                                                shiny::actionButton("save_uncorrected_scree_plot", "Add Plot to Zip File Export"),
                                                                shiny::actionButton("save_uncorrected_PCA_components", "Add PCA Components to Zip File Export"),
                                                                shiny::p(),
                                                                DT::dataTableOutput("uncorrected_contribution_table"),
                                                                shiny::actionButton("save_uncorrected_contribution_table", "Add table to Zip File Export"),
                                                                shiny::hr(),
                                                                DT::dataTableOutput("uncorrected_contribution_counts"),
                                                                shiny::actionButton("save_uncorrected_contribution_counts", "Add table to Zip File Export")
                                                  ),
                                                  shiny::column(6,
                                                                p(),
                                                                shiny::uiOutput("PCADet2title"),
                                                                shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_scree_plot")), type = 6),
                                                                shiny::actionButton("save_corrected_scree_plot", "Add Plot to Zip File Export"),
                                                                shiny::actionButton("save_corrected_PCA_components", "Add PCA Components to Zip File Export"),
                                                                shiny::p(),
                                                                DT::dataTableOutput("corrected_contribution_table"),
                                                                shiny::actionButton("save_corrected_contribution_table", "Add table to Zip File Export"),
                                                                shiny::hr(),
                                                                DT::dataTableOutput("corrected_contribution_counts"),
                                                                shiny::actionButton("save_corrected_contribution_counts", "Add table to Zip File Export")
                                                  )
                                                ),
                                                value = "pca_dt"
                                              ),
                                              shiny::tabPanel(
                                                "UMAP",
                                                shiny::p(),
                                                shiny::fluidRow(
                                                  shiny::column(6,
                                                                p(),
                                                                shiny::uiOutput("UMAP1title"),
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
                                                                p(),
                                                                shiny::uiOutput("UMAP2title"),
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
                                          shiny::tabPanel(
                                            "Cluster Plots",
                                            p(),
                                            tabsetPanel(
                                              id = "Cluster_main_pan",
                                              tabPanel("Cluster Plots",
                                                       p(),
                                                       shiny::fluidRow(
                                                         shiny::column(6,
                                                                       p(),
                                                                       shiny::uiOutput("Cluster1title"),
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
                                                                       p(),
                                                                       shiny::uiOutput("Cluster2title"),
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
                                              tabPanel("Heatmaps",
                                                       p(),
                                                       shiny::fluidRow(
                                                         shiny::column(6,
                                                                       p(),
                                                                       shiny::uiOutput("Heatmap1title"),
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
                                                                       p(),
                                                                       shiny::uiOutput("Heatmap2title"),
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
                                              ),
                                              tabPanel("Diversity Evaluation",
                                                       p(),
                                                       shiny::fluidRow(
                                                         shiny::column(6,
                                                                       p(),
                                                                       shiny::uiOutput("DivEval1title"),
                                                                       shinyjqui::jqui_resizable(shiny::plotOutput("uncorr_bar_plot", height = "300px")),
                                                                       shiny::fluidRow(
                                                                         column(4, style = 'padding-right:0px;',
                                                                                shiny::actionButton("save_uncorrected_barplot", "Add to Zip File Export")
                                                                         ),
                                                                         column(4, style = 'padding-left:0px;',
                                                                                downloadButton("dnldsave_uncorrected_barplot","Dowload Single Plot")
                                                                         )
                                                                       ),
                                                                       shiny::uiOutput("DivEvalHEAvg1title"),
                                                                       p(),
                                                                       uiOutput("rendErrorMessageUncorr"),
                                                                       tabsetPanel(
                                                                         tabPanel("Heterogeneity",
                                                                                  DT::dataTableOutput("uncorr_avg_het"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_uncorr_avg_het_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_uncorr_avg_het_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         ),
                                                                         tabPanel("Eveness",
                                                                                  DT::dataTableOutput("uncorr_avg_evn"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_uncorr_avg_evn_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_uncorr_avg_evn_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         )
                                                                       ),
                                                                       shiny::uiOutput("DivEvalHE1title"),
                                                                       p(),
                                                                       tabsetPanel(
                                                                         tabPanel("Heterogeneity",
                                                                                  div(DT::dataTableOutput("uncorr_het"), style = "overflow-X: scroll"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_uncorr_het_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_uncorr_het_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         ),
                                                                         tabPanel("Eveness",
                                                                                  div(DT::dataTableOutput("uncorr_evn"), style = "overflow-X: scroll"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_uncorr_evn_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_uncorr_evn_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         )
                                                                       )
                                                         ),
                                                         shiny::column(6,
                                                                       p(),
                                                                       shiny::uiOutput("DivEval2title"),
                                                                       shinyjqui::jqui_resizable(shiny::plotOutput("corr_bar_plot", height = "300px")),
                                                                       shiny::fluidRow(
                                                                         column(4, style = 'padding-right:0px;',
                                                                                shiny::actionButton("save_corrected_barplot", "Add to Zip File Export")
                                                                         ),
                                                                         column(4, style = 'padding-left:0px;',
                                                                                downloadButton("dnldsave_corrected_barplot","Dowload Single Plot")
                                                                         )
                                                                       ),
                                                                       shiny::uiOutput("DivEvalHEAvg2title"),
                                                                       p(),
                                                                       uiOutput("rendErrorMessageCorr"),
                                                                       tabsetPanel(
                                                                         tabPanel("Heterogeneity",
                                                                                  DT::dataTableOutput("corr_avg_het"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_corr_avg_het_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_corr_avg_het_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         ),
                                                                         tabPanel("Eveness",
                                                                                  DT::dataTableOutput("corr_avg_evn"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_corr_avg_evn_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_corr_avg_evn_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         )
                                                                       ),
                                                                       shiny::uiOutput("DivEvalHE2title"),
                                                                       p(),
                                                                       tabsetPanel(
                                                                         tabPanel("Heterogeneity",
                                                                                  div(DT::dataTableOutput("corr_het"), style = "overflow-X: scroll"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_corr_het_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_corr_het_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         ),
                                                                         tabPanel("Eveness",
                                                                                  div(DT::dataTableOutput("corr_evn"), style = "overflow-X: scroll"),
                                                                                  shiny::fluidRow(
                                                                                    column(4, style = 'padding-right:0px;',
                                                                                           shiny::actionButton("save_corr_evn_df", "Add to Zip File Export")
                                                                                    ),
                                                                                    column(4, style = 'padding-left:0px;',
                                                                                           downloadButton("dnldsave_corr_evn_df","Dowload Table")
                                                                                    )
                                                                                  )
                                                                         )
                                                                       )
                                                         )
                                                       ),
                                                       h4("Meta Information with Cluster Data"),
                                                       DT::dataTableOutput("ClusterInfoTab"),
                                                       p(),
                                                       shiny::fluidRow(
                                                         column(2, style = 'padding-right:0px;',
                                                                shiny::actionButton("save_ClusterInfoTab", "Add to Zip File Export")
                                                         ),
                                                         column(2, style = 'padding-left:0px;',
                                                                downloadButton("dnldsave_ClusterInfoTab","Dowload Table")
                                                         )
                                                       ),
                                                       value = "clinfo",
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
                                                            p(),
                                                            shiny::uiOutput("RLE1title"),
                                                            shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_RLE_plot")), type = 6),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_uncorrected_RLE_plot", "Add to Zip File Export")
                                                              ),
                                                              column(4, style = 'padding-left:0px;',
                                                                     downloadButton("dnldsave_uncorrected_RLE_plot","Dowload Single Plot")
                                                              )
                                                            )
                                              ),
                                              shiny::column(6,
                                                            p(),
                                                            shiny::uiOutput("RLE2title"),
                                                            shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_RLE_plot")), type = 6),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_corrected_RLE_plot", "Add to Zip File Export")
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
                                              shiny::tabPanel("Explanatory Variables Plot",
                                                              shiny::fluidRow(
                                                                shiny::column(6,
                                                                              p(),
                                                                              shiny::uiOutput("ExpVar1title"),
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
                                                                              p(),
                                                                              shiny::uiOutput("ExpVar2title"),
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
                                              shiny::tabPanel("PVCA",
                                                              shiny::fluidRow(
                                                                shiny::column(6,
                                                                              p(),
                                                                              shiny::uiOutput("PVCA1title"),
                                                                              uiOutput("rendErrorMessagePVCAUncorr"),
                                                                              shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_PVCA_plot")),
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
                                                                              p(),
                                                                              shiny::uiOutput("PVCA2title"),
                                                                              uiOutput("rendErrorMessagePVCACorr"),
                                                                              shinyjqui::jqui_resizable(shiny::plotOutput("corrected_PVCA_plot")),
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
                                                            p(),
                                                            shiny::uiOutput("SVA1title"),
                                                            shiny::fluidRow(
                                                              shiny::column(6,
                                                                            shiny::selectInput("SVAxAxis1","X-Axis Surrogate Variable",
                                                                                               choices = c("SVA_1","SVA_2"), width = "70%")
                                                              ),
                                                              shiny::column(6,
                                                                            shiny::selectInput("SVAyAxis1","X-Axis Surrogate Variable",
                                                                                               choices = c("SVA_1","SVA_2"), selected = "SVA_2", width = "70%")
                                                              )
                                                            ),
                                                            uiOutput("rendSVAtextError1"),
                                                            #shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_SVA_probability_density")), type = 6),
                                                            shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotly::plotlyOutput("uncorrected_SVA_scatter")), type = 6),
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
                                                            DT::dataTableOutput("uncorrected_SVA_scatter_df"),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_uncorrected_SVA_scatter_df", "Add to Zip File Export")
                                                              ),
                                                              column(4, style = 'padding-left:0px;',
                                                                     downloadButton("dnldsave_uncorrected_SVA_scatter_df","Dowload Table")
                                                              )
                                                            )
                                                            #p(),
                                                            #verbatimTextOutput("uncorrected_SVA_nsv_print"),
                                              ),
                                              shiny::column(6,
                                                            p(),
                                                            shiny::uiOutput("SVA2title"),
                                                            shiny::fluidRow(
                                                              shiny::column(6,
                                                                            shiny::selectInput("SVAxAxis2","X-Axis Surrogate Variable",
                                                                                               choices = c("SVA_1","SVA_2"), width = "70%")
                                                              ),
                                                              shiny::column(6,
                                                                            shiny::selectInput("SVAyAxis2","X-Axis Surrogate Variable",
                                                                                               choices = c("SVA_1","SVA_2"), selected = "SVA_2", width = "70%")
                                                              )
                                                            ),
                                                            uiOutput("rendSVAtextError2"),
                                                            #shinycssloaders::withSpinner(shinyjqui::jqui_resizable(shiny::plotOutput("corrected_SVA_probability_density")), type = 6),
                                                            shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotly::plotlyOutput("corrected_SVA_scatter")), type = 6),
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
                                                            DT::dataTableOutput("corrected_SVA_scatter_df"),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_corrected_SVA_scatter_df", "Add to Zip File Export")
                                                              ),
                                                              column(4, style = 'padding-left:0px;',
                                                                     downloadButton("dnldsave_corrected_SVA_scatter_df","Dowload Table")
                                                              )
                                                            )
                                                            #p(),
                                                            #verbatimTextOutput("corrected_SVA_nsv_print"),
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
                                                            p(),
                                                            shiny::uiOutput("BP1title"),
                                                            shiny::uiOutput("renduncorrected_Box_plot"),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_uncorrected_Box_plot", "Add to Zip File Export")
                                                              ),
                                                              column(4, style = 'padding-left:0px;',
                                                                     downloadButton("dnldsave_uncorrected_Box_plot","Dowload Single Plot")
                                                              )
                                                            ),
                                                            p(),
                                                            DT::dataTableOutput("uncorrected_Box_plot_df"),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_uncorrected_Box_plot_df", "Add to Zip File Export")
                                                              ),
                                                              column(4, style = 'padding-left:0px;',
                                                                     downloadButton("dnldsave_uncorrected_Box_plot_df","Dowload Table")
                                                              )
                                                            )
                                              ),
                                              shiny::column(6,
                                                            p(),
                                                            shiny::uiOutput("BP2title"),
                                                            shiny::uiOutput("rendcorrected_Box_plot"),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_corrected_Box_plot", "Add to Zip File Export")
                                                              ),
                                                              column(4, style = 'padding-left:0px;',
                                                                     downloadButton("dnldsave_corrected_Box_plot","Dowload Single Plot")
                                                              )
                                                            ),
                                                            p(),
                                                            DT::dataTableOutput("corrected_Box_plot_df"),
                                                            p(),
                                                            shiny::fluidRow(
                                                              column(4, style = 'padding-right:0px;',
                                                                     shiny::actionButton("save_corrected_Box_plot_df", "Add to Zip File Export")
                                                              ),
                                                              column(4, style = 'padding-left:0px;',
                                                                     downloadButton("dnldsave_corrected_Box_plot_df","Dowload Table")
                                                              )
                                                            )
                                                            
                                              )
                                            ),
                                            value = "box"
                                          )
                                        )
                                      )
                                    )
                    ),
                    ## Zip Export Tab----------------------------------------------
                    shiny::tabPanel("Step 3 - Zip File Export",
                                    shiny::div(
                                      id = "panel",
                                      style = "width: 600px; max-width: 100%; margin: 0 auto; padding-top: 100px;",
                                      shiny::wellPanel(
                                        shiny::p(),
                                        shiny::textInput("file_name", "Please Name the Zipped File", value = "BatchFlex", width = "100%"),
                                        shiny::p(),
                                        pickerInput(
                                          "select_save_files",
                                          "Choose Files to Zip",
                                          list.files(temp_directory),
                                          options = list(`actions-box` = TRUE),
                                          multiple = TRUE, width = "100%"
                                        ),
                                        p(),
                                        shiny::downloadButton(
                                          outputId = "download_btn",
                                          label = "Download All Selected Files as Zip",
                                          icon = icon("file-download"), width = "100%"
                                        ),
                                        shiny::downloadButton(
                                          outputId = "download_logs",
                                          label = "Download workflow log transcript",
                                          icon = icon("file-download"), width = "100%"
                                        )
                                      )
                                    )
                    ),
                    ## Tutorial Tab --------------------------------------------
                    shiny::tabPanel("Tutorial", 
                                    bslib::page_fillable(
                                      bslib::layout_sidebar(
                                        sidebar = bslib::accordion(
                                          bslib::accordion_panel("Data Input",
                                                                 shiny::fluidRow(
                                                                   shiny::column(12,
                                                                                 actionButton("load_microarray_example", "Load Microarray Example")),
                                                                   shiny::column(12,
                                                                                 actionButton("load_rnacounts_example", "Load RNACount Example")),
                                                                   shiny::column(12, 
                                                                                 actionButton("input_data", "Input Data")),
                                                                   shiny::column(12,
                                                                                 actionButton("merge_data", "Merge Data")),
                                                                   shiny::column(12,
                                                                                 actionButton("log_data", "Log Data")),
                                                                   shiny::column(12,
                                                                                 actionButton("select_species", "Select Species")),
                                                                   shiny::column(12, 
                                                                                 actionButton("select_batch", "Select Batch Variable"))
                                                                 )
                                          ),
                                          bslib::accordion_panel("Set Up Comparison",
                                                                 shiny::fluidRow(
                                                                   shiny::column(12,
                                                                                 actionButton("evaluate_batch", "Set Up Batch Eval")),
                                                                   shiny::column(12,
                                                                                 actionButton("evaluate_correction", "Set Up Correction Eval"))
                                                                 )
                                          ),
                                          bslib::accordion_panel("Batch Correction",
                                                                 shiny::fluidRow(
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_limma", "Limma")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_combat", "ComBat")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_combatseq", "ComBatSeq")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_qnorm", "Quantile Normalization")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_meancentering", "Mean Centering")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_harman", "Harman")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_ruvg_lognorm", "RUVg LogNorm")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_ruvg_counts", "RUVg Counts")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_correction_sva", "SVA"))
                                                                 )
                                          ),
                                          bslib::accordion_panel("Batch Evaluation",
                                                                 shiny::fluidRow(
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_pca", "PCA")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_mcpca", "MC PCA")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_pcadetails", "PCA Details")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_umap", "UMAP")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_cluster", "Cluster Plots")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_heatmap", "Heatmap")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_diversity", "Diversity Evaluation")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_rle", "Relative Log Expression")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_ev", "Explanatory Variables")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_pvca", "PVCA")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_sva", "SVA")),
                                                                   shiny::column(12,
                                                                                 actionButton("batch_evaluation_box", "Box Plot")),
                                                                 )
                                          ),
                                          bslib::accordion_panel("Figure Adjustment",
                                                                 shiny::fluidRow(
                                                                   shiny::column(12,
                                                                                 actionButton("figure_settings", "Figure Settings"))
                                                                 )
                                          ),
                                          bslib::accordion_panel("Data Export",
                                                                 shiny::fluidRow(
                                                                   shiny::column(12,
                                                                                 actionButton("data_export", "Data Export"))
                                                                 )
                                          )
                                          
                                        ),
                                        bslib::layout_column_wrap(
                                          bslib::layout_columns(
                                            bslib::card(bslib::card_header("Tutorial Video"),
                                                        bslib::card_body(
                                                          uiOutput("rendTutorialHelpText"),
                                                          uiOutput("video")
                                                          )
                                            ),
                                            bslib::card(bslib::card_header("Helpful Information"),
                                                        bslib::card_body(shiny::htmlOutput("tutorial_text")),
                                                        full_screen = TRUE
                                            ),
                                            col_widths = c(8,4)
                                          ),
                                          height = "850px",
                                          max_height = "850px"
                                        )
                                      )
                                    ),
                                    value = "tutorial"
                    )
                    #shiny::tabPanel("Tutorial", 
                    #                bslib::page_fillable(
                    #                  bslib::layout_sidebar(
                    #                    sidebar = bslib::accordion(
                    #                      bslib::accordion_panel("Data Input",
                    #                                             shiny::column(10,
                    #                                                           actionButton("load_example", "Load Example")),
                    #                                             shiny::column(10,
                    #                                                           actionButton("log_matrix", "Log Dataset")),
                    #                                             shiny::column(10,
                    #                                                           actionButton("select_species", "Select Species")),
                    #                                             shiny::column(10, 
                    #                                                           actionButton("input_matrix", "Input Matrix")),
                    #                                             shiny::column(10,
                    #                                                           actionButton("input_meta", "Input Meta")),
                    #                                             shiny::column(10, 
                    #                                                           actionButton("select_batch", "Select Batch Variable"))
                    #                      ),
                    #                      bslib::accordion_panel("Batch Correction",
                    #                                             shiny::fluidRow(
                    #                                               shiny::column(10,
                    #                                                             actionButton("load_example", "Load Example")),
                    #                                               shiny::column(10,
                    #                                                             actionButton("test", "TEST"))
                    #                                             )
                    #                      ),
                    #                      bslib::accordion_panel("Batch Evaluation",
                    #                                             shiny::fluidRow(
                    #                                               shiny::column(10,
                    #                                                             actionButton("load_example", "Load Example")),
                    #                                               shiny::column(10,
                    #                                                             actionButton("test", "TEST"))
                    #                                             )
                    #                      ),
                    #                      bslib::accordion_panel("Figure Adjustment",
                    #                                             shiny::fluidRow(
                    #                                               shiny::column(10,
                    #                                                             actionButton("load_example", "Load Example")),
                    #                                               shiny::column(10,
                    #                                                             actionButton("test", "TEST"))
                    #                                             )
                    #                      ),
                    #                      bslib::accordion_panel("Data Export",
                    #                                             shiny::fluidRow(
                    #                                               shiny::column(10,
                    #                                                             actionButton("load_example", "Load Example")),
                    #                                               shiny::column(10,
                    #                                                             actionButton("test", "TEST"))
                    #                                             )
                    #                      )
                    #                    ),
                    #                    uiOutput("video")
                    #                  )
                    #                ),
                    #                value = "tutorial"
                    #)
  )



# Server -----------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Homepage -------------------------------------------------------------------
  
  #bslib::bs_themer()
  # Homepage text
  output$homepage_text <- renderText(homepage_tutorial_text_list[["Homepage"]]$text)
  #output$homepage_text <- renderText(homepage_text$text)
  # Homepage homepage
  output$homepage <- renderImage({
    list(src = "www/BatchFLEX_Figure_1_JD_V8_20240813.png",
         contentType = 'image/png',
         width = 600,
         height = 800)
  }, deleteFile = FALSE)
  
  # Homepage links
  output$shiny_app <- renderImage({
    list(src = "www/shiny_app.png",
         contentType = 'image/png',
         width = 180,
         height = 180)
  }, deleteFile = FALSE)
  
  output$shiny_github <- renderImage({
    list(src = "www/shiny_github_nocopyright.png",
         contentType = 'image/jpg',
         width = 180,
         height = 180)
  }, deleteFile = FALSE)
  
  output$function_github <- renderImage({
    list(src = "www/function_github_nocopyright.png",
         contentType = 'image/png',
         width = 180,
         height = 180)
  }, deleteFile = FALSE)
  
  output$shiny_mergeqc <- renderImage({
    list(src = "www/Mergeqc_app.png",
         contentType = 'image/png',
         width = 180,
         height = 180)
  }, deleteFile = FALSE)
  
  output$mergeqc_github <- renderImage({
    list(src = "www/Mergeqc_github.png",
         contentType = 'image/png',
         width = 180,
         height = 180)
  }, deleteFile = FALSE)
  
  output$journal <- renderImage({
    list(src = "www/bioinfo_40_7cover.jpeg",
         contentType = 'image/jpeg',
         width = 180,
         height = 200)
  }, deleteFile = FALSE)
  
  observeEvent(input$citation, {
    text <- paste0("Test Citation Text")
    session$sendCustomMessage("txt", text)
  })
  
  observeEvent(input$get_started, {
    updateTabsetPanel(session, "navbar_id", selected = "step_1")
  })
  
  observeEvent(input$need_tutorial, {
    updateTabsetPanel(session, "navbar_id", selected = "tutorial")
  })
  
  # Tutorial -------------------------------------------------------------------
  TutorialPlayed <- reactiveVal(0)
  output$rendTutorialHelpText <- renderUI({
    
    if (TutorialPlayed() == 0) {
      h3(HTML("<b>← Click on a topic to the left to view a tutorial video</b>"), 
         style="text-align:center")
    } else {
      NULL
    }
    
  })
  
  # Tutorial Data Input
  observeEvent(input$load_microarray_example, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/ylQkEUiehyg?playlist=ylQkEUiehyg&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["load_microarray_example"]]$text)
  })
  observeEvent(input$load_rnacounts_example, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/4JaukZLc2QY?playlist=4JaukZLc2QY&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["load_rnacounts_example"]]$text)
  })
  observeEvent(input$input_data, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/9Zwn2HBGfjM?playlist=9Zwn2HBGfjM&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["input_data"]]$text)
  })
  observeEvent(input$merge_data, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/24ZZKPyzVas?playlist=24ZZKPyzVas&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["merge_data"]]$text)
  })
  observeEvent(input$log_data, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/UToZcttTFzo?playlist=UToZcttTFzo&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["log_data"]]$text)
  })
  observeEvent(input$select_species, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/Awt3Ei-GX_U?playlist=Awt3Ei-GX_U&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["select_species"]]$text)
  })
  observeEvent(input$select_batch, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/ZXtPsu95bdI?playlist=ZXtPsu95bdI&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["select_batch"]]$text)
  })
  # Set up Comparison
  observeEvent(input$evaluate_batch, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/8R2UGGIYfU0?playlist=8R2UGGIYfU0&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["evaluate_batch"]]$text)
  })
  observeEvent(input$evaluate_correction, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/5oihHITsd6w?playlist=5oihHITsd6w&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["evaluate_correction"]]$text)
  })
  
  # Tutorial Data Correction
  observeEvent(input$batch_correction_limma, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/EwJDa5N5am8?playlist=EwJDa5N5am8&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_limma"]]$text)
  })
  observeEvent(input$batch_correction_combat, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/NgiI8I-2Uz0?playlist=NgiI8I-2Uz0&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_combat"]]$text)
  })
  observeEvent(input$batch_correction_combatseq, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/29nmaPexgPU?playlist=29nmaPexgPU&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_combatseq"]]$text)
  })
  observeEvent(input$batch_correction_qnorm, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/RpZg5m8_z1g?playlist=RpZg5m8_z1g&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_qnorm"]]$text)
  })
  observeEvent(input$batch_correction_meancentering, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/il0Sg-uiPO8?playlist=il0Sg-uiPO8&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_meancentering"]]$text)
  })
  observeEvent(input$batch_correction_harman, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/YdU6p2HBpRs?playlist=YdU6p2HBpRs&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_harman"]]$text)
  })
  observeEvent(input$batch_correction_ruvg_lognorm, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/XeOy-fjsGIo?playlist=XeOy-fjsGIo&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_ruvg_lognorm"]]$text)
  })
  observeEvent(input$batch_correction_ruvg_counts, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/nrpcFMC6EQM?playlist=nrpcFMC6EQM&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_ruvg_counts"]]$text)
  })
  observeEvent(input$batch_correction_sva, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/jcijaq-ypRo?playlist=jcijaq-ypRo&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_correction_sva"]]$text)
  })
  # Tutorial Batch Evaluation
  observeEvent(input$batch_evaluation_pca, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/9Tlnqp1fGi4?playlist=9Tlnqp1fGi4&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_pca"]]$text)
  })
  observeEvent(input$batch_evaluation_mcpca, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/qNsfXBRm-5A?playlist=qNsfXBRm-5A&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_mcpca"]]$text)
  })
  observeEvent(input$batch_evaluation_pcadetails, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/_RMXgpK_dP8?playlist=_RMXgpK_dP8&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_pcadetails"]]$text)
  })
  observeEvent(input$batch_evaluation_umap, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/VNqotbPLH0g?playlist=VNqotbPLH0g&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_umap"]]$text)
  })
  observeEvent(input$batch_evaluation_cluster, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/WiQHMXlkxZY?playlist=WiQHMXlkxZY&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_cluster"]]$text)
  })
  observeEvent(input$batch_evaluation_heatmap, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/ulvA0YljtLc?playlist=ulvA0YljtLc&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_heatmap"]]$text)
  })
  observeEvent(input$batch_evaluation_diversity, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/bTV6mDgC-us?playlist=bTV6mDgC-us&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_diversity"]]$text)
  })
  observeEvent(input$batch_evaluation_rle, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/EKuNIle2kq4?playlist=EKuNIle2kq4&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_rle"]]$text)
  })
  observeEvent(input$batch_evaluation_ev, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/Z3Gd3Gti9Bs?playlist=Z3Gd3Gti9Bs&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_ev"]]$text)
  })
  observeEvent(input$batch_evaluation_pvca, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/aQ-UI5_Qvi8?playlist=aQ-UI5_Qvi8&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_pvca"]]$text)
  })
  observeEvent(input$batch_evaluation_sva, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/XKCf6GcEdtw?playlist=XKCf6GcEdtw&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_sva"]]$text)
  })
  observeEvent(input$batch_evaluation_box, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/pOrrC6gWHfg?playlist=pOrrC6gWHfg&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["batch_evaluation_box"]]$text)
  })
  # Tutorial Figure Settings
  observeEvent(input$figure_settings, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/1j5YBZNJGTc?playlist=1j5YBZNJGTc&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["figure_settings"]]$text)
  })
  
  # Tutorial Data Export
  observeEvent(input$data_export, {
    TutorialPlayed(1)
    output$video <- renderUI({
      HTML('<iframe width="100%" height="100%" src="//www.youtube.com/embed/BTZWA8m4xBg?playlist=BTZWA8m4xBg&&loop=1;rel=0&autoplay=1&mute=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
    })
    output$tutorial_text <- renderText(homepage_tutorial_text_list[["data_export"]]$text)
  })
  
  #observeEvent(input$load_example, {
  #  output$video <- renderUI({
  #    HTML('<iframe width="1200" height="600" src="//www.youtube.com/embed/wE4fhPzw2TY?playlist=wE4fhPzw2TY&&loop=1;rel=0&autoplay=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
  #  })
  #})
  #observeEvent(input$log_matrix, {
  #  output$video <- renderUI({
  #    HTML('<iframe width="1200" height="600" src="//www.youtube.com/embed/wE4fhPzw2TY?playlist=wE4fhPzw2TY&&loop=1;rel=0&autoplay=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
  #  })
  #})
  #observeEvent(input$select_species, {
  #  output$video <- renderUI({
  #    HTML('<iframe width="1200" height="600" src="//www.youtube.com/embed/wE4fhPzw2TY?playlist=wE4fhPzw2TY&&loop=1;rel=0&autoplay=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
  #  })
  #})
  #observeEvent(input$input_matrix, {
  #  output$video <- renderUI({
  #    HTML('<iframe width="1200" height="600" src="//www.youtube.com/embed/wE4fhPzw2TY?playlist=wE4fhPzw2TY&&loop=1;rel=0&autoplay=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
  #  })
  #})
  #observeEvent(input$input_meta, {
  #  output$video <- renderUI({
  #    HTML('<iframe width="1200" height="600" src="//www.youtube.com/embed/wE4fhPzw2TY?playlist=wE4fhPzw2TY&&loop=1;rel=0&autoplay=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
  #  })
  #})
  #observeEvent(input$select_batch, {
  #  output$video <- renderUI({
  #    HTML('<iframe width="1200" height="600" src="//www.youtube.com/embed/wE4fhPzw2TY?playlist=wE4fhPzw2TY&&loop=1;rel=0&autoplay=1&controls=0&showinfo=0" frameborder="0" allowfullscreen></iframe>')
  #  })
  #})
  #
  ## Homepage tutorial batch correction
  #
  ## Homepage tutorial batch evaluation
  #
  ## Homepage tutorial adjust figures
  #
  ## Homepage tutorial export
  
  # Reactive Val Start --------------------------------------------------------
  #matrix_file_raw <- reactiveVal()
  #meta_file_raw <- reactive
  matrix_raw <- reactiveVal()
  matrix_raw_ForOut <- reactiveVal()
  meta_raw <- reactiveVal()
  uncorrected_matrix <- reactiveVal()
  user_provided_batch_info <- reactiveVal()
  aligned_meta_file <- reactiveVal()
  boxplot_meta_file <- reactiveVal()
  uncorr_boxplot_meta_file <- reactiveVal()
  corr_boxplot_meta_file <- reactiveVal()
  RawCountCheck <- reactiveVal()
  files_to_download <- shiny::reactiveValues()
  MM_react <- reactiveVal(FALSE)
  FileCheckAlerts_react <- reactiveVal()
  FileCheckAlerts_reactR <- reactiveVal()
  FileCheckAlerts_reactL <- reactiveVal()
  
  observe({
    RawCountCheck(input$RawCountInput)
  })
  
  observe({
    if (RawCountCheck()) {
      if (isTruthy(input$RawMatTabs)) {
        if (input$RawMatTabs == 1) {
          if ((!input$batch_correction_method %in% c("Uncorrected","ComBatseq","RUVg")) | (!input$batch_correction_method2 %in% c("Uncorrected","ComBatseq","RUVg"))) {
            showModal(modalDialog(
              title = "Batch correction method applied in matrix transformation tab",
              paste("The batch correction method you selected will be applied to this dataset and can be visualized in the following tab.
                  The current tab is solely for viewing un-normalized RNASeq counts."),
              easyClose = TRUE,
              footer = paste("Click outside of box to exit message")
            ))
          }
        }
      }
    }
  })
  
  
  # Batch Correction Inputs ---------------------------------------------------
  output$rendbatch_correction_method <- renderUI({
    if (RawCountCheck()) {
      shiny::selectInput("batch_correction_method",
                         "Batch Correction Method 1:",
                         c("Uncorrected", "Limma", "ComBat","ComBatseq (RNAseq Counts Only)" = "ComBatseq", "Quantile Normalization", "Mean Centering", "Harman", "RUVg","SVA"),
                         selected = "Uncorrected")
    } else {
      shiny::selectInput("batch_correction_method",
                         "Batch Correction Method 1:",
                         c("Uncorrected", "Limma", "ComBat", "Quantile Normalization", "Mean Centering", "Harman", "RUVg","SVA"))
    }
  })
  
  output$rendbatch_correction_method2 <- renderUI({
    if (RawCountCheck()) {
      shiny::selectInput("batch_correction_method2",
                         "Batch Correction Method 2:",
                         c("Uncorrected", "Limma", "ComBat","ComBatseq (RNAseq Counts Only)" = "ComBatseq", "Quantile Normalization", "Mean Centering", "Harman", "RUVg","SVA"),
                         selected = "ComBatseq")
    } else {
      shiny::selectInput("batch_correction_method2",
                         "Batch Correction Method 2:",
                         c("Uncorrected", "Limma", "ComBat", "Quantile Normalization", "Mean Centering", "Harman", "RUVg","SVA"))
    }
  })
  
  observeEvent(input$RNAcountsMethod,{
    rawMethod <- input$RNAcountsMethod
    if (rawMethod == "RUVg") {
      updateSelectInput(session,"batch_correction_method2",selected = "RUVg")
    } else if (rawMethod == "ComBatseq") {
      updateSelectInput(session,"batch_correction_method2",selected = "ComBatseq")}
  })
  
  output$rendQuantNorm <- renderUI({
    req(input$batch_correction_method)
    if (!any(c("Uncorrected","Quantile Normalization") %in% input$batch_correction_method)) {
      #if (RawCountCheck()) {
      #req(input$RawCountNorm)
      #if (input$RawCountNorm != "none") {
      checkboxInput("QuantNorm","Quantile Normalize",value = F)
      #}
      #} else {
      #  checkboxInput("QuantNorm","Quantile Normalize",value = F)
      #}
    }
    
  })
  
  output$rendQuantNorm2 <- renderUI({
    req(input$batch_correction_method2)
    if (!any(c("Uncorrected","Quantile Normalization") %in% input$batch_correction_method2)) {
      #if (RawCountCheck()) {
      #req(input$RawCountNorm)
      #  if (input$RawCountNorm != "none") {
      checkboxInput("QuantNorm2","Quantile Normalize",value = F)
      #  }
      #} else {
      #  checkboxInput("QuantNorm2","Quantile Normalize",value = F)
      #}
    }
    
  })
  
  # Data Input ----------------------------------------------------------------
  
  DataLoded <- reactiveVal(0)
  output$rendDataLodedHelpText <- renderUI({
    
    if (DataLoded() == 0) {
      h3(HTML("<b>← Upload user data or load an example dataset to get started</b>"), 
         style="margin-top:80px")
    } else {
      NULL
    }
    
  })
  
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
  output$rendbatch_names <- shiny::renderUI({
    if(shiny::isTruthy(input$batch_name_delim)){
      if(input$batch_info == "Yes"){
        shiny::textInput("batch_names",
                         "Please Provide Batch Names Separated by Commas")
      }
    }
  })
  
  shiny::observeEvent(input$UseExpData, {
    DataLoded(1)
    withProgress(message = "Processing", value = 0, {
      incProgress(0.5, detail = "Loading Micro Array Example data")
      uncorrected_matrix_read <- readr::read_delim(Example_Matrix_File,
                                                   delim = '\t',
                                                   name_repair = "minimal",
                                                   na = c("", "NA", "N/A"))
      meta_read <- readr::read_delim(Example_Meta_File,
                                     delim = '\t',
                                     name_repair = "minimal",
                                     na = c("", "NA", "N/A"))
      incProgress(0.5, detail = "Complete!")
      matrix_raw(uncorrected_matrix_read)
      meta_raw(meta_read)
      MM_react(TRUE)
      updateRadioButtons(session,"HumanOrMouse",selected = "Mouse")
      updateSelectInput(session,"Init_Batch_Select",selected = "Study")
    })
  })
  shiny::observeEvent(input$UseRawExpData, {
    DataLoded(1)
    withProgress(message = "Processing", value = 0, {
      incProgress(0.5, detail = "Loading Micro Array Example data")
      uncorrected_matrix_read <- readr::read_delim(Example_MatrixCounts_File,
                                                   delim = '\t',
                                                   name_repair = "minimal",
                                                   na = c("", "NA", "N/A"))
      matrix_raw(uncorrected_matrix_read)
      meta_read <- readr::read_delim(Example_MetaCounts_File,
                                     delim = '\t',
                                     name_repair = "minimal",
                                     na = c("", "NA", "N/A"))
      incProgress(0.5, detail = "Complete!")
      RawCountCheck(TRUE)
      meta_raw(meta_read)
      MM_react(FALSE)
      updateRadioButtons(session,"HumanOrMouse",selected = "Human")
      updateSelectInput(session,"Init_Batch_Select",selected = "study")
    })
  })
  
  
  observe({
    if (isTruthy(input$uncorrected_matrix_input$datapath)) {
      DataLoded(1)
      withProgress(message = "Processing", value = 0, {
        incProgress(0.5, detail = "Loading User Matrix")
        uncorrected_matrix_read <- readr::read_delim(input$uncorrected_matrix_input$datapath,
                                                     delim = input$matrix_delim,
                                                     name_repair = "minimal",
                                                     na = c("", "NA", "N/A"))
        incProgress(0.5, detail = "Complete!")
        matrix_raw(uncorrected_matrix_read)
      })
    }
    if (isTruthy(input$user_provided_batch_info$datapath)) {
      DataLoded(1)
      withProgress(message = "Processing", value = 0, {
        incProgress(0.5, detail = "Loading User Meta Data")
        meta_read <- readr::read_delim(input$user_provided_batch_info$datapath,
                                       delim = input$meta_delim,
                                       name_repair = "minimal",
                                       na = c("", "NA", "N/A"))
        incProgress(0.5, detail = "Complete!")
        meta_raw(meta_read)
      })
    }
    #updateRadioButtons(session,"HumanOrMouse",selected = "Human")
  })
  
  meta_file_delim <- shiny::reactive({
    print(input$batch_name_delim)
  })
  meta_file_names <- shiny::reactive({
    req(input$batch_names)
    print(unlist(strsplit(input$batch_names, ",")))
  })
  observe({
    if (input$batch_info == "Yes"){
      req(meta_file_delim())
      req(meta_file_names())
      meta_names <- colnames(matrix_raw())
      meta_names <- as.data.frame(meta_names[-1])
      meta_names <- meta_names %>% dplyr::rename( "Original Sample Name" = "meta_names[-1]")
      meta_names <- cbind(
        meta_names,
        tidyr::separate_wider_delim(
          meta_names,
          cols = 1,
          delim = meta_file_delim(),
          names = meta_file_names(),
          too_few = "align_start",
          too_many = "merge"
        )
      )
      meta_names
      meta_raw(meta_names)
    } else if (input$batch_info == "No"){
      meta_raw()
    }
  })
  
  observe({
    
    if (isTruthy(matrix_raw()) & isTruthy(meta_raw())) {
      
      set.seed(input$SeedSet)
      uncorrected_matrix_read <- matrix_raw()
      user_provided_batch_info <- meta_raw()
      
      SpecDetect <- detect_species(uncorrected_matrix_read[,1])
      if (SpecDetect == "mouse") {
        updateRadioButtons(session,"HumanOrMouse",NULL,c("Human","Mouse"), selected = "Mouse")
        MM_react(TRUE)
      } else {
        updateRadioButtons(session,"HumanOrMouse",NULL,c("Human","Mouse"), selected = "Human")
        MM_react(FALSE)
      }
      
      colnames(uncorrected_matrix_read)[1] <- "Genes"
      if (TRUE %in% duplicated(uncorrected_matrix_read[, 1])) {
        data_dup <- uncorrected_matrix_read %>% dplyr::group_by(Genes) %>% dplyr::filter(n() > 1) %>% as.data.frame()
        data_nodup <- uncorrected_matrix_read %>% dplyr::group_by(Genes) %>% dplyr::filter(n() == 1) %>% as.data.frame()
        data_dup_summ <- data_dup %>%
          dplyr::group_by(Genes) %>%
          dplyr::summarise_all(max) %>%
          as.data.frame()
        uncorrected_matrix_nodupes <- rbind(data_nodup,data_dup_summ)
      } else {
        uncorrected_matrix_nodupes <- uncorrected_matrix_read
      }
      num_test <- apply(uncorrected_matrix_nodupes[,-1],2, is.numeric)
      if (all(num_test) == FALSE) {
        uncorrected_matrix_nodupes[,-1] <- mutate_all(uncorrected_matrix_nodupes[,-1], function(x) as.numeric(as.character(x)))
      }
      
      uncorrected_matrix_nodupes <- as.data.frame(uncorrected_matrix_nodupes)
      user_provided_batch_info <- as.data.frame(user_provided_batch_info)
      sampSame <- intersect(colnames(uncorrected_matrix_nodupes)[-1],user_provided_batch_info[,1])
      uncorrected_matrix_nodupes <- uncorrected_matrix_nodupes[,c("Genes",sampSame)]
      user_provided_batch_info <- user_provided_batch_info[which(user_provided_batch_info[,1] %in% sampSame),]
      
      mat <- uncorrected_matrix_nodupes
      mat <- as.matrix(mat[,-1])
      if (all(mat == floor(mat))) {
        updateCheckboxInput(session,"RawCountInput", value = T)
        RawCountCheck(TRUE)
      } else {
        updateCheckboxInput(session,"RawCountInput", value = F)
        RawCountCheck(FALSE)
      }
      
      #uncorrected_matrix_nodupes$filter <- apply(
      #  uncorrected_matrix_nodupes[,-1],
      #  1,
      #  function(x) ExprFilter2(x, 1, 0.05)
      #)
      #uncorrected_matrix_filtered <- uncorrected_matrix_nodupes[which(uncorrected_matrix_nodupes$filter == TRUE),]
      #uncorrected_matrix_filtered <- uncorrected_matrix_filtered[,-ncol(uncorrected_matrix_filtered)]
      
      uncorrected_matrix_filtered <- uncorrected_matrix_nodupes
      rownames(uncorrected_matrix_filtered) <- uncorrected_matrix_filtered[,1]
      uncorrected_numeric_matrix <- as.matrix(uncorrected_matrix_filtered[,-1])
      
      uncorrected_numeric_matrix <- uncorrected_numeric_matrix[which(rowSums(uncorrected_numeric_matrix, na.rm = T) > 0),]
      
      ## Normalize raw counts
      if (RawCountCheck()) {
        matrix_raw_ForOut(uncorrected_numeric_matrix)
      }
      ## Log counts
      #if (input$Log_Choice){
      #  withProgress(message = "Processing", value = 0, {
      #    incProgress(0.5, detail = "Logging Matrix")
      #    uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix + 1)
      #    incProgress(0.5, detail = "Complete!")
      #  })
      #}
      
      # align sample names in meta file with matrix, make sure same order
      aligned_meta_file_filtered <- user_provided_batch_info[match(colnames(uncorrected_numeric_matrix), user_provided_batch_info[,1]),]
      
      #print(head(aligned_meta_file_filtered))
      
      uncorrected_matrix(uncorrected_numeric_matrix)
      user_provided_batch_info(aligned_meta_file_filtered)
      aligned_meta_file(aligned_meta_file_filtered)
      boxplot_meta_file(aligned_meta_file_filtered)
      uncorr_boxplot_meta_file(aligned_meta_file_filtered)
      corr_boxplot_meta_file(aligned_meta_file_filtered)
      
    }
    
  })
  observe({
    updateNumericInput(session,"SVAvarNum",value = nrow(uncorrected_matrix()))
  })
  #output$rendRawCountNorm <- renderUI({
  #  if (RawCountCheck()) {
  #    selectInput("RawCountNorm","Normalize Raw Counts By:",
  #                c("No Normalization" = "none","TMM" = "TMM","Upper Quartile" = "upperquartile"), selected = "TMM")
  #  }
  #})
  observe({
    if (RawCountCheck()) {
      #req(input$RawCountNorm)
      #if (input$RawCountNorm == "none") {
      updateCheckboxInput(session,"Log_Choice", value = F)
      #}
    }
  })
  
  output$rendMatrixHeader <- renderUI({
    req(uncorrected_matrix())
    h3("Matrix Preview")
  })
  output$rendExprHead <- renderUI({
    req(uncorrected_matrix())
    radioButtons("ExprHead",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  output$uncorrected_matrix_output_input <- DT::renderDataTable({
    if (RawCountCheck()) {
      req(matrix_raw_ForOut())
      req(input$ExprHead)
      expr <- as.data.frame(matrix_raw_ForOut())
      expr <- cbind(data.frame(Genes = rownames(expr)),
                    expr)
    } else {
      req(uncorrected_matrix())
      req(input$ExprHead)
      expr <- as.data.frame(uncorrected_matrix())
      expr <- cbind(data.frame(Genes = rownames(expr)),
                    expr)
    }
    if (input$ExprHead == "View table head") {
      expr <- head(expr,c(100,100))
    }
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  output$rendMetaHeader <- renderUI({
    req(user_provided_batch_info())
    h3("Meta Preview")
  })
  output$rendClinHead <- renderUI({
    req(user_provided_batch_info())
    radioButtons("ClinHead",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  output$meta_file <- DT::renderDataTable({
    req(user_provided_batch_info())
    req(input$ClinHead)
    meta <- user_provided_batch_info()
    if (input$ClinHead == "View table head") {
      meta <- head(meta,c(100,100))
    }
    DT::datatable(meta,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
    
  })
  batch_names_from_meta <- shiny::reactive({
    batch_names_from_meta <- colnames(user_provided_batch_info())
    batch_names_from_meta
  })
  
  
  # Batch Criteria Selection --------------------------------------------------
  
  observe({
    req(input$batch_correction_method)
    if (input$batch_correction_method == "Quantile Normalization") {
      updateCheckboxInput(session,"QuantNorm",value = FALSE)
    }
  })
  
  # Initial batch selection
  output$rendInit_Batch_Select <- shiny::renderUI({
    req(user_provided_batch_info())
    batch_choices <- batch_names_from_meta()[-1]
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
      BatchSelected <- c("Select", batch_names_from_meta()[-1])[2]
    }
    shiny::selectInput("batch1_choices_limma",
                       "Batch 1 Variable",
                       c("Select", batch_names_from_meta()[-1]),
                       selected = BatchSelected)
  })
  output$batch2_selection_limma <- shiny::renderUI({
    shiny::selectInput("batch2_choices_limma",
                       "Batch 2 Variable",
                       c("Select", batch_names_from_meta()[-1]))
  })
  output$covariate_choices_limma <- shiny::renderUI({
    shiny::selectInput("covariate_choices_limma",
                       "Select Any Covariates for Correction",
                       batch_names_from_meta()[-1],
                       multiple = TRUE
    )
  })
  # Selection of batch criteria for correction.
  output$batch1_selection_limma2 <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c("Select", batch_names_from_meta()[-1])[2]
    }
    shiny::selectInput("batch1_choices_limma2",
                       "Batch 1 Variable",
                       c("Select", batch_names_from_meta()[-1]),
                       selected = BatchSelected)
  })
  output$batch2_selection_limma2 <- shiny::renderUI({
    shiny::selectInput("batch2_choices_limma2",
                       "Batch 2 Variable",
                       c("Select", batch_names_from_meta()[-1]))
  })
  output$covariate_choices_limma2 <- shiny::renderUI({
    shiny::selectInput("covariate_choices_limma2",
                       "Select Any Covariates for Correction",
                       c(batch_names_from_meta()[-1]),
                       multiple = TRUE
    )
  })
  output$batch1_selection_ComBatseq <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    shiny::selectInput("batch1_choices_ComBatseq",
                       "Select Variable for Batch",
                       c(batch_names_from_meta()[-1]),
                       selected = BatchSelected)
  })
  output$covariate_choices_ComBatseq <- shiny::renderUI({
    shiny::selectInput("covariate_choices_ComBatseq",
                       "Select Biological Variables",
                       c(batch_names_from_meta()[-1]),
                       multiple = TRUE
    )
  })
  #output$rendbatch1_choices_ComBatseq_raw <- shiny::renderUI({
  #  if (!is.null(InitialBatchSelected)) {
  #    BatchSelected <- InitialBatchSelected()
  #  } else {
  #    BatchSelected <- c(batch_names_from_meta()[-1])[1]
  #  }
  #  shiny::selectInput("batch1_choices_ComBatseq_raw",
  #                     "Select Variable for Batch",
  #                     c(batch_names_from_meta()[-1]),
  #                     selected = BatchSelected)
  #})
  #output$rendcovariate_choices_ComBatseq_raw <- shiny::renderUI({
  #  shiny::selectInput("covariate_choices_ComBatseq_raw",
  #                     "Select Biological Variables",
  #                     c(batch_names_from_meta()[-1]),
  #                     multiple = TRUE
  #  )
  #})
  output$batch1_selection_ComBatseq2 <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    shiny::selectInput("batch1_choices_ComBatseq2",
                       "Select Variable for Batch",
                       c(batch_names_from_meta()[-1]),
                       selected = BatchSelected)
  })
  output$covariate_choices_ComBatseq2 <- shiny::renderUI({
    shiny::selectInput("covariate_choices_ComBatseq2",
                       "Select Biological Variables",
                       c(batch_names_from_meta()[-1]),
                       multiple = TRUE
    )
  })
  output$batch_selection_ComBat <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    shiny::selectInput(
      "batch_selection_ComBat",
      "Select Batch",
      c(batch_names_from_meta()[-1]),
      selected = BatchSelected
    )
  })
  output$batch_selection_ComBat2 <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    shiny::selectInput(
      "batch_selection_ComBat2",
      "Select Batch",
      c(batch_names_from_meta()[-1]),
      selected = BatchSelected
    )
  })
  output$batch_selection_mean_centering <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    shiny::selectInput(
      "batch1_choices_mean_centering",
      "Select Batch for Mean Centering",
      c(batch_names_from_meta()[-1]),
      selected = BatchSelected
    )
  })
  output$batch_selection_mean_centering2 <- shiny::renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    shiny::selectInput(
      "batch1_choices_mean_centering2",
      "Select Batch for Mean Centering",
      c(batch_names_from_meta()[-1]),
      selected = BatchSelected
    )
  })
  output$batch_selection_harman <- renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    selectInput(
      "batch1_choices_harman",
      "Select Batch",
      c(batch_names_from_meta()[-1]),
      selected = BatchSelected
    )
  })
  output$treatment_selection_harman <- renderUI({
    treatOpt <- c(batch_names_from_meta()[-1])
    harmBatch <- input$batch1_choices_harman
    selectInput(
      "treatment_selection_harman",
      "Variable of Interest",
      c(batch_names_from_meta()[-1]),
      treatOpt[which(treatOpt != harmBatch)][1]
    )
  })
  output$batch_selection_harman2 <- renderUI({
    if (!is.null(InitialBatchSelected)) {
      BatchSelected <- InitialBatchSelected()
    } else {
      BatchSelected <- c(batch_names_from_meta()[-1])[1]
    }
    selectInput(
      "batch1_choices_harman2",
      "Select Batch",
      c(batch_names_from_meta()[-1]),
      selected = BatchSelected
    )
  })
  output$treatment_selection_harman2 <- renderUI({
    treatOpt <- c(batch_names_from_meta()[-1])
    harmBatch <- input$batch1_choices_harman2
    selectInput(
      "treatment_selection_harman2",
      "Variable of Interest",
      c(batch_names_from_meta()[-1]),
      treatOpt[which(treatOpt != harmBatch)][1]
    )
  })
  output$SVA_variable_of_interest_bc <- renderUI({
    VarChoices <- c(batch_names_from_meta()[-1])
    VarSelect <- VarChoices[which(VarChoices != InitialBatchSelected())]
    selectInput(
      "SVA_variable_of_interest_bc",
      "Variable of Interest",
      VarChoices,
      selected = VarSelect[1]
    )
  })
  output$SVA_variable_of_interest_bc2 <- renderUI({
    VarChoices <- c(batch_names_from_meta()[-1])
    VarSelect <- VarChoices[which(VarChoices != InitialBatchSelected())]
    selectInput(
      "SVA_variable_of_interest_bc2",
      "Variable of Interest",
      VarChoices,
      selected = VarSelect[1]
    )
  })
  
  
  # Allowing select input to provide NULL input into the functions
  batch1_choices1 <- shiny::reactive({
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
  batch1_choices2 <- shiny::reactive({
    if (input$batch_correction_method2 == "Limma"){
      if (isTruthy(input$batch1_choices_limma2)) {
        if (input$batch1_choices_limma2 == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_limma2
        }
      }
    }else if (input$batch_correction_method2 == "ComBat"){
      if (isTruthy(input$batch_selection_ComBat2)) {
        if (input$batch_selection_ComBat2 == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch_selection_ComBat2
        }
      }
    }else if (input$batch_correction_method2 == "ComBatseq"){
      if (isTruthy(input$batch1_choices_ComBatseq2)) {
        if (input$batch1_choices_ComBatseq2 == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_ComBatseq2
        }
      }
    }else if (input$batch_correction_method2 == "Mean Centering"){
      if (isTruthy(input$batch1_choices_mean_centering2)) {
        if (input$batch1_choices_mean_centering2 == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_mean_centering2
        }
      }
    }else if (input$batch_correction_method2 == "Harman"){
      if (isTruthy(input$batch1_choices_harman2)) {
        if (input$batch1_choices_harman2 == "Select"){
          batch1_choices <- NULL
        } else{
          batch1_choices <- input$batch1_choices_harman2
        }
      }
    }
  })
  
  batch1_choices <- reactive({
    
    batch_select <- unique(c(batch1_choices1(),batch1_choices2()))
    batch_select[which(isTruthy(batch_select))][1]
    
    if (isTruthy(batch_select)) {
      batch <- batch_select
      batch
    } else {
      if (isTruthy(InitialBatchSelected())) {
        batch <- InitialBatchSelected()
        batch
      } else {
        batch <- NULL
        batch
      }
    }
    
  })
  
  #batch1_choices <- reactive({
  #  batch_select <- unique(c(batch1_choices1(),batch1_choices2()))
  #  batch_select[which(isTruthy(batch_select))][1]
  #})
  batch2_choices1 <- shiny::reactive({
    if (isTruthy(input$batch2_choices_limma)) {
      if (input$batch2_choices_limma == "Select"){
        batch2_choices <- NULL
      } else{
        batch2_choices <- input$batch2_choices_limma
      }
    }
  })
  batch2_choices2 <- shiny::reactive({
    if (isTruthy(input$batch2_choices_limma2)) {
      if (input$batch2_choices_limma2 == "Select"){
        batch2_choices <- NULL
      } else{
        batch2_choices <- input$batch2_choices_limma2
      }
    }
  })
  batch2_choices <- reactive({
    batch_select <- unique(c(batch2_choices(),batch2_choices2()))
    batch_select[which(isTruthy(batch_select))][1]
  })
  harman_treatment_input1 <- reactive({
    if (isTruthy(input$treatment_selection_harman)) {
      if (input$treatment_selection_harman == "Select"){
        NULL
      }else {
        harman_treatment_input <- input$treatment_selection_harman
      }
    }
  })
  harman_treatment_input2 <- reactive({
    if (isTruthy(input$treatment_selection_harman2)) {
      if (input$treatment_selection_harman2 == "Select"){
        NULL
      }else {
        harman_treatment_input <- input$treatment_selection_harman2
      }
    }
  })
  harman_treatment_input <- reactive({
    batch_select <- unique(c(harman_treatment_input1(),harman_treatment_input2()))
    batch_select[which(isTruthy(batch_select))][1]
  })
  SVA_variable_of_interest_bc1 <- reactive({
    if (isTruthy(input$SVA_variable_of_interest_bc)) {
      if (input$SVA_variable_of_interest_bc == "Select"){
        NULL
      }else {
        SVA_variable_of_interest_bc2 <- input$SVA_variable_of_interest_bc
      }
    }
  })
  SVA_variable_of_interest_bc2 <- reactive({
    if (isTruthy(input$SVA_variable_of_interest_bc2)) {
      if (input$SVA_variable_of_interest_bc2 == "Select"){
        NULL
      }else {
        SVA_variable_of_interest_bc2 <- input$SVA_variable_of_interest_bc2
      }
    }
  })
  SVA_variable_of_interest_bc <- reactive({
    batch_select <- unique(c(SVA_variable_of_interest_bc1(),SVA_variable_of_interest_bc2()))
    batch_select[which(isTruthy(batch_select))][1]
  })
  
  # Model Matrix ------------------------------------------------
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
  
  model_matrix2 <- shiny::reactive({
    if (input$batch_correction_method2 == "Limma"){
      if (is.null(input$covariate_choices_limma2)){
        model_matrix <- NULL
      } else {
        total_covariates <- paste0(input$covariate_choices_limma2,collapse = "+")
        model_matrix <- stats::model.matrix(reformulate(total_covariates), data = as.data.frame(aligned_meta_file()))
      }
    } else if (input$batch_correction_method2 == "ComBatseq"){
      if (is.null(input$covariate_choices_ComBatseq2)){
        model_matrix <- NULL
      } else {
        counter <- 0
        for (covariate in 1:length(input$covariate_choices_ComBatseq2)){
          variable_object <- input$covariate_choices_ComBatseq2[covariate]
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
  RUVg_housekeeping1 <- reactive({
    if (input$HumanOrMouse == "Human") {
      if(input$RUVg_housekeeping_selection1 == "Eisenberg"){
        RUVg_housekeeping <- Eisenberg_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection1 == "Lin500"){
        RUVg_housekeeping <- Lin500_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection1 == "HSIAO"){
        RUVg_housekeeping <- HSIAO_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection1 == "UserInput"){
        RUVg_housekeeping_input <- read.delim(input$RUVg_user_control_genes1$datapath, sep = '\t')
        RUVg_housekeeping <- RUVg_housekeeping_input[,1]
        RUVg_housekeeping
      }
    } else {
      if(input$RUVg_housekeeping_selection1 == "Eisenberg"){
        RUVg_housekeeping <- Eisenberg_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection1 == "Lin500"){
        RUVg_housekeeping <- Lin500_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection1 == "HSIAO"){
        RUVg_housekeeping <- HSIAO_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection1 == "UserInput"){
        RUVg_housekeeping_input <- read.delim(input$RUVg_user_control_genes1$datapath, sep = '\t')
        RUVg_housekeeping <- RUVg_housekeeping_input[,1]
        RUVg_housekeeping
      }
    }
    
  })
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
  RUVg_housekeeping2 <- reactive({
    if (input$HumanOrMouse == "Human") {
      if(input$RUVg_housekeeping_selection2 == "Eisenberg"){
        RUVg_housekeeping <- Eisenberg_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection2 == "Lin500"){
        RUVg_housekeeping <- Lin500_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection2 == "HSIAO"){
        RUVg_housekeeping <- HSIAO_hkg[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection2 == "UserInput"){
        RUVg_housekeeping_input <- read.delim(input$RUVg_user_control_genes2$datapath, sep = '\t')
        RUVg_housekeeping <- RUVg_housekeeping_input[,1]
        RUVg_housekeeping
      }
    } else {
      if(input$RUVg_housekeeping_selection2 == "Eisenberg"){
        RUVg_housekeeping <- Eisenberg_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection2 == "Lin500"){
        RUVg_housekeeping <- Lin500_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection2 == "HSIAO"){
        RUVg_housekeeping <- HSIAO_hkg_mm[,1]
        RUVg_housekeeping
      }else if(input$RUVg_housekeeping_selection2 == "UserInput"){
        RUVg_housekeeping_input <- read.delim(input$RUVg_user_control_genes2$datapath, sep = '\t')
        RUVg_housekeeping <- RUVg_housekeeping_input[,1]
        RUVg_housekeeping
      }
    }
    
  })
  
  # Alerts ---------------------------------------------------------------------
  
  observe({
    
    req(input$batch_correction_method)
    req(input$batch_correction_method2)
    FileCheckAlerts_list <- c()
    FileCheckAlerts_listR <- c()
    FileCheckAlerts_listL <- c()
    if (isTruthy(input$uncorrected_matrix_input$datapath)) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                paste("Input matrix uploaded:",input$uncorrected_matrix_input$datapath,"."))
    }
    if (isTruthy(input$user_provided_batch_info$datapath)) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                paste("Input meta data uploaded:",input$user_provided_batch_info$datapath,"."))
    }
    if (input$UseExpData > 0) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                paste("Micro array example data loaded."))
    }
    if (input$UseRawExpData > 0) {
      FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                paste("RNASeq count example data loaded."))
    }
    
    FileCheckAlerts_list <- c(FileCheckAlerts_list,
                              paste("---"),
                              paste("Log data for matrix developed by method 1 (left)."))
    FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                              paste("Log data for matrix developed by method 1."))
    if (RawCountCheck()) {
      req(input$RawCountNorm)
      if (input$batch_correction_method %in% c("ComBatseq","RUVg")){
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix",input$batch_correction_method,"corrected."),
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                  paste("Input matrix",input$batch_correction_method,"corrected."),
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        if (isTruthy(input$Log_Choice_Raw)) {
          if (input$Log_Choice_Raw) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix log2 transformed."))
            FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                      paste("Input matrix log2 transformed."))
          }
        }
        if (isTruthy(input$QuantNorm)) {
          if (input$QuantNorm) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix quantile normalized."))
            FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                      paste("Input matrix quantile normalized."))
          }
        }
      } else {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        if (isTruthy(input$Log_Choice_Raw)) {
          if (input$Log_Choice_Raw) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix log2 transformed."))
            FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                       paste("Input matrix log2 transformed."))
          }
        }
        if (isTruthy(input$QuantNorm)) {
          if (input$QuantNorm & (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization"))) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix quantile normalized."))
            FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                       paste("Input matrix quantile normalized."))
          }
        }
        if (input$batch_correction_method == "Uncorrected") {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix", input$batch_correction_method))
          FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                    paste("Input matrix", input$batch_correction_method))
        } else {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix", input$batch_correction_method,"corrected."))
          FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                    paste("Input matrix", input$batch_correction_method,"corrected."))
        }
      }
    } else {
      if (isTruthy(input$Log_Choice)) {
        if (input$Log_Choice) {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix log2 transformed."))
          FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                    paste("Input matrix log2 transformed."))
        }
      }
      if (isTruthy(input$QuantNorm)) {
        if (input$QuantNorm & (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization"))) {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix quantile normalized."))
          FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                    paste("Input matrix quantile normalized."))
        }
      }
      if (input$batch_correction_method == "Uncorrected") {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix", input$batch_correction_method))
        FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                  paste("Input matrix", input$batch_correction_method))
      } else {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix", input$batch_correction_method,"corrected."))
        FileCheckAlerts_listL <- c(FileCheckAlerts_listL,
                                  paste("Input matrix", input$batch_correction_method,"corrected."))
      }
    }
    
    FileCheckAlerts_list <- c(FileCheckAlerts_list,
                              paste("---"),
                              paste("Log data for matrix developed by method 2 (right)."))
    FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                              paste("Log data for matrix developed by method 2."))
    if (RawCountCheck()) {
      req(input$batch_correction_method2)
      req(input$RawCountNorm)
      if (input$batch_correction_method2 %in% c("ComBatseq","RUVg")){
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix",input$batch_correction_method2,"corrected."),
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                  paste("Input matrix",input$batch_correction_method2,"corrected."),
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        if (isTruthy(input$Log_Choice_Raw)) {
          if (input$Log_Choice_Raw) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix log2 transformed."))
            FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                      paste("Input matrix log2 transformed."))
          }
        }
        if (isTruthy(input$QuantNorm2)) {
          if (input$QuantNorm2 & (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization"))) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix quantile normalized."))
            FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                      paste("Input matrix quantile normalized."))
          }
        }
      } else {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                  paste("Input matrix",input$RawCountNorm,"normalized."))
        if (isTruthy(input$Log_Choice_Raw)) {
          if (input$Log_Choice_Raw) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix log2 transformed."))
          }
          FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                    paste("Input matrix log2 transformed."))
        }
        if (isTruthy(input$QuantNorm2)) {
          if (input$QuantNorm2 & (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization"))) {
            FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                      paste("Input matrix quantile normalized."))
            FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                      paste("Input matrix quantile normalized."))
          }
        }
        if (input$batch_correction_method2 == "Uncorrected") {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix", input$batch_correction_method2))
          FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                    paste("Input matrix", input$batch_correction_method2))
        } else {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix", input$batch_correction_method2,"corrected."))
          FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                    paste("Input matrix", input$batch_correction_method2,"corrected."))
        }
      }
    } else {
      if (isTruthy(input$Log_Choice)) {
        if (input$Log_Choice) {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix log2 transformed."))
          FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                    paste("Input matrix log2 transformed."))
        }
      }
      if (isTruthy(input$QuantNorm2)) {
        if (input$QuantNorm2 & (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization"))) {
          FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                    paste("Input matrix quantile normalized."))
          FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                    paste("Input matrix quantile normalized."))
        }
      }
      if (input$batch_correction_method2 == "Uncorrected") {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix", input$batch_correction_method2))
        FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                  paste("Input matrix", input$batch_correction_method2))
      } else {
        FileCheckAlerts_list <- c(FileCheckAlerts_list,
                                  paste("Input matrix", input$batch_correction_method2,"corrected."))
        FileCheckAlerts_listR <- c(FileCheckAlerts_listR,
                                  paste("Input matrix", input$batch_correction_method2,"corrected."))
      }
    }
    FileCheckAlerts_react(FileCheckAlerts_list)
    FileCheckAlerts_reactR(FileCheckAlerts_listR)
    FileCheckAlerts_reactL(FileCheckAlerts_listL)
    
  })
  output$download_logs <- downloadHandler(
    filename = function() {
      paste("BatchFLEX_Processing_Notes_", format(Sys.time(),format = "%Y%m%d_%H%M"), ".txt", sep = "")
    },
    content = function(file) {
      df <- FileCheckAlerts_react()
      writeLines(df, file)
    }
  )
  output$FileCheckAlertsR <- renderPrint({
    
    req(FileCheckAlerts_reactR())
    text <- paste(FileCheckAlerts_reactR(), collapse = "\n")
    cat(text)
    
  })
  output$FileCheckAlertsR2 <- renderPrint({
    
    req(FileCheckAlerts_reactR())
    text <- paste(FileCheckAlerts_reactR(), collapse = "\n")
    cat(text)
    
  })
  output$FileCheckAlertsL <- renderPrint({
    
    req(FileCheckAlerts_reactL())
    text <- paste(FileCheckAlerts_reactL(), collapse = "\n")
    cat(text)
    
  })
  output$FileCheckAlertsL2 <- renderPrint({
    
    req(FileCheckAlerts_reactL())
    text <- paste(FileCheckAlerts_reactL(), collapse = "\n")
    cat(text)
    
  })
  
  # Correction -----------------------------------------------------------------
  matrix1_Raw <- shiny::reactive({
    
    uncorrected_numeric_matrix <- matrix_raw_ForOut()
    aligned_meta_file <- aligned_meta_file()
    
    batch_correction <- uncorrected_numeric_matrix
    
    if (exists("batch_correction")) {
      if (isTruthy(batch_correction)) {
        file_name <- paste("RNAseqCounts_InputMatrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
        #file_name <- paste(gsub(" ","",matrix1Title_react()),"_RNAseqCounts_matrix", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
        df <- as.data.frame(batch_correction)
        df <- cbind(data.frame(Genes = rownames(df)),
                    df)
        readr::write_tsv(df, file.path(temp_directory, file_name))
        batch_correction
      }
    }
    
  })
  
  ComBatseq_react1 <- reactive({
    
    uncorrected_numeric_matrix <- matrix_raw_ForOut()
    aligned_meta_file <- aligned_meta_file()
    req(matrix_raw_ForOut())
    req(input$batch_correction_method)
    req(input$batch1_choices_ComBatseq)
    if (is.null(batch1_choices1())) {
      batchChoice <- input$batch1_choices_ComBatseq
    } else {
      batchChoice <- batch1_choices1()
    }
    withProgress(message = "Processing", value = 0, {
      incProgress(0.5, detail = "Running ComBatseq Batch Correction")
      combatseq_corrected <- sva::ComBat_seq(
        as.matrix(uncorrected_numeric_matrix),
        batch = c(unlist(aligned_meta_file[,batchChoice])),
        covar_mod = model_matrix()
      )
      batch_correction <- combatseq_corrected
      incProgress(0.5, detail = "Complete!")
      batch_correction
    })
    
  })
  ComBatseq_react2 <- reactive({
    
    uncorrected_numeric_matrix <- matrix_raw_ForOut()
    aligned_meta_file <- aligned_meta_file()
    req(matrix_raw_ForOut)
    req(input$batch_correction_method2)
    req(input$batch1_choices_ComBatseq2)
    if (is.null(batch1_choices2())) {
      batchChoice <- input$batch1_choices_ComBatseq2
    } else {
      batchChoice <- batch1_choices2()
    }
    withProgress(message = "Processing", value = 0, {
      incProgress(0.5, detail = "Running ComBatseq Batch Correction")
      combatseq_corrected <- sva::ComBat_seq(
        as.matrix(uncorrected_numeric_matrix),
        batch = c(unlist(aligned_meta_file[,batchChoice])),
        covar_mod = model_matrix2()
      )
      batch_correction <- combatseq_corrected
      incProgress(0.5, detail = "Complete!")
      batch_correction
    })
    
  })
  
  matrix2_Raw <- shiny::reactive({
    
    uncorrected_numeric_matrix <- matrix_raw_ForOut()
    aligned_meta_file <- aligned_meta_file()
    
    #if (input$batch_correction_method2 == "ComBatseq") {
    #  batch_correction <- ComBatseq_react2()
    #  if (exists("batch_correction")) {
    #    if (isTruthy(batch_correction)) {
    #      file_name <- paste("RNAseqCounts_CombatSeqCorrected_matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    #      df <- as.data.frame(batch_correction)
    #      df <- cbind(data.frame(Genes = rownames(df)),
    #                  df)
    #      readr::write_tsv(df, file.path(temp_directory, file_name))
    #      batch_correction
    #    }
    #  }
    #} else if (input$batch_correction_method2 == "RUVg") {
    if (input$batch_correction_method2 == "RUVg") {
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
            round = T,
            tolerance = input$RUVg_tolerance,
            isLog = F
          )
          RUVg_correction_matrix <- as.data.frame(RUVg_correction$normalizedCounts)
          batch_correction <- RUVg_correction_matrix
          incProgress(0.5, detail = "Complete!")
        })
      }
      if (exists("batch_correction")) {
        if (isTruthy(batch_correction)) {
          file_name <- paste("RNAseqCounts_RUVgCorrected_matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
          df <- as.data.frame(batch_correction)
          df <- cbind(data.frame(Genes = rownames(df)),
                      df)
          readr::write_tsv(df, file.path(temp_directory, file_name))
          batch_correction
        }
      }
    } else if (input$batch_correction_method2 == "Uncorrected") {
      batch_correction <- matrix_raw_ForOut()
      if (exists("batch_correction")) {
        if (isTruthy(batch_correction)) {
          file_name <- paste("RNAseqCounts_InputMatrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
          df <- as.data.frame(batch_correction)
          df <- cbind(data.frame(Genes = rownames(df)),
                      df)
          readr::write_tsv(df, file.path(temp_directory, file_name))
          batch_correction
        }
      }
    } else {
      batch_correction <- ComBatseq_react2()
      if (exists("batch_correction")) {
        if (isTruthy(batch_correction)) {
          file_name <- paste("RNAseqCounts_CombatSeqCorrected_matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
          df <- as.data.frame(batch_correction)
          df <- cbind(data.frame(Genes = rownames(df)),
                      df)
          readr::write_tsv(df, file.path(temp_directory, file_name))
          batch_correction
        }
      }
    }
    
  })
  
  
  
  matrix1 <- shiny::reactive({
    
    uncorrected_numeric_matrix <- uncorrected_matrix()
    aligned_meta_file <- aligned_meta_file()
    req(input$batch_correction_method)
    
    # Raw input data
    if (RawCountCheck()) {
      uncorrected_numeric_matrix <- matrix_raw_ForOut()
      # ComBatseq method
      if (input$batch_correction_method == "ComBatseq") {
        batch_correction <- ComBatseq_react1()
        # Normalize
        if (isTruthy(batch_correction)) {
          #if (RawCountCheck()) {
            req(input$RawCountNorm)
            mat <- batch_correction
            mat_dgeList <- DGEList(counts = as.matrix(mat))
            withProgress(message = "Processing", value = 0, {
              incProgress(0.25, detail = "Normalizing Factors")
              mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
              incProgress(0.25, detail = "CPM")
              batch_correction <- edgeR::cpm(mat_dgeList_Norm)
              incProgress(0.75, detail = "Complete!")
            })
          #}
          # Log 
          if (input$Log_Choice_Raw) {
            batch_correction <- log2(batch_correction+1)
          }
          # Quantile Normalize counts
          if (isTruthy(input$QuantNorm)) {
            if (input$QuantNorm){
              withProgress(message = "Processing", value = 0, {
                incProgress(0.5, detail = "Quantile Normalization")
                uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(batch_correction), keep.names = T)
                batch_correction <- as.data.frame(uncorrected_numeric_matrix_proc)
                incProgress(0.5, detail = "Complete!")
              })
            }
          }
        }
        # RUVg Method
      } else if (input$batch_correction_method == "RUVg") {
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
              round = T,
              tolerance = input$RUVg_tolerance,
              isLog = T
            )
            RUVg_correction_matrix <- as.matrix(RUVg_correction$normalizedCounts)
            RUVg_correction_matrix[which(RUVg_correction_matrix < 0)] <- 0
            batch_correction <- RUVg_correction_matrix
            incProgress(0.5, detail = "Complete!")
          })
        }
        # Normalize
        if (isTruthy(batch_correction)) {
          #if (RawCountCheck()) {
          #req(input$RawCountNorm)
          mat <- batch_correction
          mat_dgeList <- DGEList(counts = as.matrix(mat))
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Normalizing Factors")
            mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
            incProgress(0.25, detail = "CPM")
            batch_correction <- edgeR::cpm(mat_dgeList_Norm)
            incProgress(0.75, detail = "Complete!")
          })
          #}
          # Log 
          if (input$Log_Choice_Raw) {
            batch_correction <- log2(batch_correction+1)
          }
          # Quantile Normalize counts
          if (isTruthy(input$QuantNorm)) {
            if (input$QuantNorm){
              withProgress(message = "Processing", value = 0, {
                incProgress(0.5, detail = "Quantile Normalization")
                uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(batch_correction), keep.names = T)
                batch_correction <- as.data.frame(uncorrected_numeric_matrix_proc)
                incProgress(0.5, detail = "Complete!")
              })
            }
          }
        }
        # Normalize first
      } else {
        #if (isTruthy(batch_correction)) {
          #if (RawCountCheck()) {
            req(input$RawCountNorm)
            mat <- uncorrected_numeric_matrix
            mat_dgeList <- DGEList(counts = as.matrix(mat))
            withProgress(message = "Processing", value = 0, {
              incProgress(0.25, detail = "Normalizing Factors")
              mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
              incProgress(0.25, detail = "CPM")
              uncorrected_numeric_matrix <- edgeR::cpm(mat_dgeList_Norm)
              incProgress(0.75, detail = "Complete!")
            })
          #}
          # Log 
          if (input$Log_Choice_Raw) {
            uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix+1)
          }
          # Quantile Normalize counts
          if (isTruthy(input$QuantNorm)) {
            if (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization")) {
              if (input$QuantNorm){
                withProgress(message = "Processing", value = 0, {
                  incProgress(0.5, detail = "Quantile Normalization")
                  uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
                  uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
                  incProgress(0.5, detail = "Complete!")
                })
              }
            }
          } # Batch Correct
            if(input$batch_correction_method == "Uncorrected"){
              batch_correction <- uncorrected_numeric_matrix
            } else if(input$batch_correction_method == "Limma"){
              if (isTruthy(input$batch1_choices_limma) & isTruthy(input$batch2_choices_limma)) {
                withProgress(message = "Processing", value = 0, {
                  incProgress(0.5, detail = "Running Limma batch correction")
                  batch_correction <- limma::removeBatchEffect(
                    uncorrected_numeric_matrix,
                    batch = c(unlist(aligned_meta_file[,batch1_choices1()])),
                    batch2 = c(unlist(aligned_meta_file[,batch2_choices1()])),
                    covariates = model_matrix()
                  )
                  incProgress(0.5, detail = "Complete!")
                })
              }
            }else if(input$batch_correction_method == "ComBat"){
              if (isTruthy(batch1_choices1())) {
                withProgress(message = "Processing", value = 0, {
                  incProgress(0.25, detail = "Running ComBat batch correction")
                  batch_combat <- c(unlist(aligned_meta_file[,batch1_choices1()]))
                  incProgress(0.25, detail = "Running ComBat batch correction")
                  modcombat <-  stats::model.matrix(~1, data = as.data.frame(aligned_meta_file))
                  incProgress(0.25, detail = "Running ComBat batch correction")
                  batch_correction <- sva::ComBat(
                    dat = uncorrected_numeric_matrix,
                    batch = batch_combat,
                    mod = modcombat,
                    par.prior = input$combat_parametric
                  )
                  incProgress(0.25, detail = "Complete!")
                })
              }
            }else if (input$batch_correction_method == "Quantile Normalization") {
              withProgress(message = "Processing", value = 0, {
                incProgress(0.5, detail = "Quantile Normalization")
                uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
                uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
                batch_correction <- uncorrected_numeric_matrix
                incProgress(0.5, detail = "Complete!")
              })
              batch_correction
            }else if(input$batch_correction_method == "Mean Centering"){
              if (isTruthy(batch1_choices1())) {
                withProgress(message = "Processing", value = 0, {
                  incProgress(0.25, detail = "Running Mean Centering Batch Correction")
                  mean_centering_batch = c(unlist(aligned_meta_file[,batch1_choices1()]))
                  mean_centering_data = t(uncorrected_numeric_matrix)
                  incProgress(0.25, detail = "Running Mean Centering Batch Correction")
                  mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
                  mean_center_correction <- as.data.frame(t(mean_center$xadj))
                  batch_correction <- mean_center_correction
                  incProgress(0.5, detail = "Complete!")
                })
              }
            } else if(input$batch_correction_method == "Harman"){
              req(harman_treatment_input())
              req(batch1_choices1())
              withProgress(message = "Processing", value = 0, {
                incProgress(0.25, detail = "Running Harman Batch Correction")
                harman_correction_PCA <- Harman::harman(
                  uncorrected_numeric_matrix,
                  expt = aligned_meta_file[,harman_treatment_input()],
                  batch = aligned_meta_file[,batch1_choices1()],
                  limit = input$HarmanCL1,
                  printInfo = T,
                  randseed = input$SeedSet
                )
                incProgress(0.25, detail = "Reconstructing Corrected Data")
                harman_correction <- Harman::reconstructData(harman_correction_PCA)
                batch_correction <- harman_correction
                incProgress(0.5, detail = "Complete!")
              })
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
                    round = F,
                    tolerance = input$RUVg_tolerance,
                    isLog = T
                  )
                  RUVg_correction_matrix <- as.data.frame(RUVg_correction$normalizedCounts)
                  batch_correction <- RUVg_correction_matrix
                  incProgress(0.5, detail = "Complete!")
                })
              }
            }else if(input$batch_correction_method == "SVA"){
              if (isTruthy(SVA_variable_of_interest_bc1())) {
                withProgress(message = "Processing", value = 0, {
                  incProgress(0.25, detail = "Number of surrogate variables")
                  expression_data <-  as.matrix(uncorrected_numeric_matrix)
                  mod <-  model.matrix(reformulate(SVA_variable_of_interest_bc1()), data = aligned_meta_file)
                  mod0 <- model.matrix(~1, data = aligned_meta_file)
                  n.sv <- sva::num.sv(expression_data, mod, method = input$svaMethod_bc)
                  incProgress(0.25, detail = "Running SVA")
                  svobj <- sva::sva(expression_data, mod , mod0, n.sv = n.sv)
                  incProgress(0.25, detail = "Running Frozen SVA for batch correction")
                  fsvaobj <- sva::fsva(expression_data, mod, svobj, expression_data)
                  batch_correction <- fsvaobj$db
                  incProgress(0.25, detail = "Complete!")
                })
              }
            }
        #}
      }
    } else {
      if (input$Log_Choice) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Logging Matrix")
          uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix+1)
          incProgress(0.5, detail = "Complete!")
        })
      }
      ## Quantile Normalize counts
      if (isTruthy(input$QuantNorm)) {
        if (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization")) {
          if (input$QuantNorm){
            withProgress(message = "Processing", value = 0, {
              incProgress(0.5, detail = "Quantile Normalization")
              uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
              uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
              incProgress(0.5, detail = "Complete!")
            })
          }
        }
      }
      ## Correct
      if(input$batch_correction_method == "Uncorrected"){
        batch_correction <- uncorrected_numeric_matrix
      } else if(input$batch_correction_method == "Limma"){
        if (isTruthy(input$batch1_choices_limma) & isTruthy(input$batch2_choices_limma)) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Running Limma batch correction")
            batch_correction <- limma::removeBatchEffect(
              uncorrected_numeric_matrix,
              batch = c(unlist(aligned_meta_file[,batch1_choices1()])),
              batch2 = c(unlist(aligned_meta_file[,batch2_choices1()])),
              covariates = model_matrix()
            )
            incProgress(0.5, detail = "Complete!")
          })
        }
      }else if(input$batch_correction_method == "ComBat"){
        if (isTruthy(batch1_choices1())) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Running ComBat batch correction")
            batch_combat <- c(unlist(aligned_meta_file[,batch1_choices1()]))
            incProgress(0.25, detail = "Running ComBat batch correction")
            modcombat <-  stats::model.matrix(~1, data = as.data.frame(aligned_meta_file))
            incProgress(0.25, detail = "Running ComBat batch correction")
            batch_correction <- sva::ComBat(
              dat = uncorrected_numeric_matrix,
              batch = batch_combat,
              mod = modcombat,
              par.prior = input$combat_parametric
            )
            incProgress(0.25, detail = "Complete!")
          })
        }
      }else if (input$batch_correction_method == "Quantile Normalization") {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Quantile Normalization")
          uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
          uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
          batch_correction <- uncorrected_numeric_matrix
          incProgress(0.5, detail = "Complete!")
        })
        batch_correction
      }else if(input$batch_correction_method == "Mean Centering"){
        if (isTruthy(batch1_choices1())) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Running Mean Centering Batch Correction")
            mean_centering_batch = c(unlist(aligned_meta_file[,batch1_choices1()]))
            mean_centering_data = t(uncorrected_numeric_matrix)
            incProgress(0.25, detail = "Running Mean Centering Batch Correction")
            mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
            mean_center_correction <- as.data.frame(t(mean_center$xadj))
            batch_correction <- mean_center_correction
            incProgress(0.5, detail = "Complete!")
          })
        }
      } else if(input$batch_correction_method == "Harman"){
        req(harman_treatment_input())
        req(batch1_choices1())
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Running Harman Batch Correction")
          harman_correction_PCA <- Harman::harman(
            uncorrected_numeric_matrix,
            expt = aligned_meta_file[,harman_treatment_input()],
            batch = aligned_meta_file[,batch1_choices1()],
            limit = input$HarmanCL1,
            printInfo = T,
            randseed = input$SeedSet
          )
          incProgress(0.25, detail = "Reconstructing Corrected Data")
          harman_correction <- Harman::reconstructData(harman_correction_PCA)
          batch_correction <- harman_correction
          incProgress(0.5, detail = "Complete!")
        })
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
              round = F,
              tolerance = input$RUVg_tolerance,
              isLog = T
            )
            RUVg_correction_matrix <- as.data.frame(RUVg_correction$normalizedCounts)
            batch_correction <- RUVg_correction_matrix
            incProgress(0.5, detail = "Complete!")
          })
        }
      }else if(input$batch_correction_method == "SVA"){
        if (isTruthy(SVA_variable_of_interest_bc1())) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Number of surrogate variables")
            expression_data <-  as.matrix(uncorrected_numeric_matrix)
            mod <-  model.matrix(reformulate(SVA_variable_of_interest_bc1()), data = aligned_meta_file)
            mod0 <- model.matrix(~1, data = aligned_meta_file)
            n.sv <- sva::num.sv(expression_data, mod, method = input$svaMethod_bc)
            incProgress(0.25, detail = "Running SVA")
            svobj <- sva::sva(expression_data, mod , mod0, n.sv = n.sv)
            incProgress(0.25, detail = "Running Frozen SVA for batch correction")
            fsvaobj <- sva::fsva(expression_data, mod, svobj, expression_data)
            batch_correction <- fsvaobj$db
            incProgress(0.25, detail = "Complete!")
          })
        }
      }
    }
    
    
    if (exists("batch_correction")) {
      if (isTruthy(batch_correction)) {
        file_name <- paste(matrix1DlndTitle_react(),"_matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
        df <- as.data.frame(batch_correction)
        df <- cbind(data.frame(Genes = rownames(df)),
                    df)
        readr::write_tsv(df, file.path(temp_directory, file_name))
        batch_correction
      }
    }
    
  })
  
  matrix2 <- shiny::reactive({
    
    uncorrected_numeric_matrix <- uncorrected_matrix()
    aligned_meta_file <- aligned_meta_file()
    req(input$batch_correction_method2)
    
    # Raw input data
    if (RawCountCheck()) {
      uncorrected_numeric_matrix <- matrix_raw_ForOut()
      # ComBatseq method
      if (input$batch_correction_method2 == "ComBatseq") {
        batch_correction <- ComBatseq_react2()
        #}
        # Normalize
        if (isTruthy(batch_correction)) {
          #if (RawCountCheck()) {
          req(input$RawCountNorm)
          mat <- batch_correction
          mat_dgeList <- DGEList(counts = as.matrix(mat))
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Normalizing Factors")
            mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
            incProgress(0.25, detail = "CPM")
            batch_correction <- edgeR::cpm(mat_dgeList_Norm)
            incProgress(0.75, detail = "Complete!")
          })
          #}
          # Log 
          if (input$Log_Choice_Raw) {
            batch_correction <- log2(batch_correction+1)
          }
          # Quantile Normalize counts
          if (isTruthy(input$QuantNorm2)) {
            if (input$QuantNorm2){
              withProgress(message = "Processing", value = 0, {
                incProgress(0.5, detail = "Quantile Normalization")
                uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(batch_correction), keep.names = T)
                batch_correction <- as.data.frame(uncorrected_numeric_matrix_proc)
                incProgress(0.5, detail = "Complete!")
              })
            }
          }
        }
        # RUVg Method
      } else if (input$batch_correction_method2 == "RUVg") {
        if(input$RUVg_housekeeping_selection2 == "UserInput"){
          req(input$RUVg_user_control_genes2)
        }
        if (any(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Running RUVg Batch Correction")
            RUVg_housekeeping_genes <- RUVg_housekeeping2()[which(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))]
            RUVg_matrix <- uncorrected_numeric_matrix
            RUVg_correction <- RUVSeq::RUVg(
              as.matrix(RUVg_matrix),
              cIdx = RUVg_housekeeping_genes,
              k = input$RUVg_estimate_factors2,
              drop = input$RUVg_drop_factors2,
              center = input$RUVg_mean_centered2,
              round = T,
              tolerance = input$RUVg_tolerance2,
              isLog = T
            )
            RUVg_correction_matrix <- as.matrix(RUVg_correction$normalizedCounts)
            RUVg_correction_matrix[which(RUVg_correction_matrix < 0)] <- 0
            batch_correction <- RUVg_correction_matrix
            incProgress(0.5, detail = "Complete!")
          })
        }
        # Normalize
        if (isTruthy(batch_correction)) {
          #if (RawCountCheck()) {
          req(input$RawCountNorm)
          mat <- batch_correction
          mat_dgeList <- DGEList(counts = as.matrix(mat))
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Normalizing Factors")
            mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
            incProgress(0.25, detail = "CPM")
            batch_correction <- edgeR::cpm(mat_dgeList_Norm)
            incProgress(0.75, detail = "Complete!")
          })
          #}
          # Log 
          if (input$Log_Choice_Raw) {
            batch_correction <- log2(batch_correction+1)
          }
          # Quantile Normalize counts
          if (isTruthy(input$QuantNorm2)) {
            if (input$QuantNorm2){
              withProgress(message = "Processing", value = 0, {
                incProgress(0.5, detail = "Quantile Normalization")
                uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(batch_correction), keep.names = T)
                batch_correction <- as.data.frame(uncorrected_numeric_matrix_proc)
                incProgress(0.5, detail = "Complete!")
              })
            }
          }
        }
        # Normalize first
      } else {
        #if (isTruthy(batch_correction)) {
        #if (RawCountCheck()) {
        #req(input$RawCountNorm)
        mat <- uncorrected_numeric_matrix
        mat_dgeList <- DGEList(counts = as.matrix(mat))
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Normalizing Factors")
          mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
          incProgress(0.25, detail = "CPM")
          uncorrected_numeric_matrix <- edgeR::cpm(mat_dgeList_Norm)
          incProgress(0.75, detail = "Complete!")
        })
        #}
        # Log 
        if (input$Log_Choice_Raw) {
          uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix+1)
        }
        # Quantile Normalize counts
        if (isTruthy(input$QuantNorm2)) {
          if (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization")) {
            if (input$QuantNorm2){
              withProgress(message = "Processing", value = 0, {
                incProgress(0.5, detail = "Quantile Normalization")
                uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
                uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
                incProgress(0.5, detail = "Complete!")
              })
            }
          }
        } # Batch Correct
        if(input$batch_correction_method2 == "Uncorrected"){
          batch_correction <- uncorrected_numeric_matrix
        } else if(input$batch_correction_method2 == "Limma"){
          if (isTruthy(input$batch1_choices_limma2) & isTruthy(input$batch2_choices_limma2)) {
            withProgress(message = "Processing", value = 0, {
              incProgress(0.5, detail = "Running Limma batch correction")
              batch_correction <- limma::removeBatchEffect(
                uncorrected_numeric_matrix,
                batch = c(unlist(aligned_meta_file[,batch1_choices2()])),
                batch2 = c(unlist(aligned_meta_file[,batch2_choices2()])),
                covariates = model_matrix2()
              )
              incProgress(0.5, detail = "Complete!")
            })
          }
        }else if(input$batch_correction_method2 == "ComBat"){
          if (isTruthy(batch1_choices2())) {
            withProgress(message = "Processing", value = 0, {
              incProgress(0.25, detail = "Running ComBat batch correction")
              batch_combat <- c(unlist(aligned_meta_file[,batch1_choices2()]))
              incProgress(0.25, detail = "Running ComBat batch correction")
              modcombat <-  stats::model.matrix(~1, data = as.data.frame(aligned_meta_file))
              incProgress(0.25, detail = "Running ComBat batch correction")
              batch_correction <- sva::ComBat(
                dat = uncorrected_numeric_matrix,
                batch = batch_combat,
                mod = modcombat,
                par.prior = input$combat_parametric2
              )
              incProgress(0.25, detail = "Complete!")
            })
          }
        }else if (input$batch_correction_method2 == "Quantile Normalization") {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Quantile Normalization")
            uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
            uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
            batch_correction <- uncorrected_numeric_matrix
            incProgress(0.5, detail = "Complete!")
          })
          batch_correction
        }else if(input$batch_correction_method2 == "Mean Centering"){
          if (isTruthy(batch1_choices2())) {
            withProgress(message = "Processing", value = 0, {
              incProgress(0.25, detail = "Running Mean Centering Batch Correction")
              mean_centering_batch = c(unlist(aligned_meta_file[,batch1_choices2()]))
              mean_centering_data = t(uncorrected_numeric_matrix)
              incProgress(0.25, detail = "Running Mean Centering Batch Correction")
              mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
              mean_center_correction <- as.data.frame(t(mean_center$xadj))
              batch_correction <- mean_center_correction
              incProgress(0.5, detail = "Complete!")
            })
          }
        } else if(input$batch_correction_method2 == "Harman"){
          req(harman_treatment_input())
          req(batch1_choices2())
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Running Harman Batch Correction")
            harman_correction_PCA <- Harman::harman(
              uncorrected_numeric_matrix,
              expt = aligned_meta_file[,harman_treatment_input2()],
              batch = aligned_meta_file[,batch1_choices2()],
              limit = input$HarmanCL2,
              printInfo = T,
              randseed = input$SeedSet
            )
            incProgress(0.25, detail = "Reconstructing Corrected Data")
            harman_correction <- Harman::reconstructData(harman_correction_PCA)
            batch_correction <- harman_correction
            incProgress(0.5, detail = "Complete!")
          })
        }else if(input$batch_correction_method2 == "RUVg"){
          if(input$RUVg_housekeeping_selection2 == "UserInput"){
            req(input$RUVg_user_control_genes2)
          }
          if (any(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))) {
            withProgress(message = "Processing", value = 0, {
              incProgress(0.5, detail = "Running RUVg Batch Correction")
              RUVg_housekeeping_genes <- RUVg_housekeeping2()[which(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))]
              RUVg_matrix <- uncorrected_numeric_matrix
              RUVg_correction <- RUVSeq::RUVg(
                as.matrix(RUVg_matrix),
                cIdx = RUVg_housekeeping_genes,
                k = input$RUVg_estimate_factors2,
                drop = input$RUVg_drop_factors2,
                center = input$RUVg_mean_centered2,
                round = F,
                tolerance = input$RUVg_tolerance2,
                isLog = T
              )
              RUVg_correction_matrix <- as.data.frame(RUVg_correction$normalizedCounts)
              batch_correction <- RUVg_correction_matrix
              incProgress(0.5, detail = "Complete!")
            })
          }
        }else if(input$batch_correction_method2 == "SVA"){
          if (isTruthy(SVA_variable_of_interest_bc2())) {
            withProgress(message = "Processing", value = 0, {
              incProgress(0.25, detail = "Number of surrogate variables")
              expression_data <-  as.matrix(uncorrected_numeric_matrix)
              mod <-  model.matrix(reformulate(SVA_variable_of_interest_bc2()), data = aligned_meta_file)
              mod0 <- model.matrix(~1, data = aligned_meta_file)
              n.sv <- sva::num.sv(expression_data, mod, method = input$svaMethod_bc2)
              incProgress(0.25, detail = "Running SVA")
              svobj <- sva::sva(expression_data, mod , mod0, n.sv = n.sv)
              incProgress(0.25, detail = "Running Frozen SVA for batch correction")
              fsvaobj <- sva::fsva(expression_data, mod, svobj, expression_data)
              batch_correction <- fsvaobj$db
              incProgress(0.25, detail = "Complete!")
            })
          }
        }
        #}
      }
    } else {
      if (input$Log_Choice) {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Logging Matrix")
          uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix+1)
          incProgress(0.5, detail = "Complete!")
        })
      }
      ## Quantile Normalize counts
      if (isTruthy(input$QuantNorm2)) {
        if (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization")) {
          if (input$QuantNorm2){
            withProgress(message = "Processing", value = 0, {
              incProgress(0.5, detail = "Quantile Normalization")
              uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
              uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
              incProgress(0.5, detail = "Complete!")
            })
          }
        }
      }
      ## Correct
      if(input$batch_correction_method2 == "Uncorrected"){
        batch_correction <- uncorrected_numeric_matrix
      } else if(input$batch_correction_method2 == "Limma"){
        if (isTruthy(input$batch1_choices_limma2) & isTruthy(input$batch2_choices_limma2)) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Running Limma batch correction")
            batch_correction <- limma::removeBatchEffect(
              uncorrected_numeric_matrix,
              batch = c(unlist(aligned_meta_file[,batch1_choices2()])),
              batch2 = c(unlist(aligned_meta_file[,batch2_choices2()])),
              covariates = model_matrix2()
            )
            incProgress(0.5, detail = "Complete!")
          })
        }
      }else if(input$batch_correction_method2 == "ComBat"){
        if (isTruthy(batch1_choices2())) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Running ComBat batch correction")
            batch_combat <- c(unlist(aligned_meta_file[,batch1_choices2()]))
            incProgress(0.25, detail = "Running ComBat batch correction")
            modcombat <-  stats::model.matrix(~1, data = as.data.frame(aligned_meta_file))
            incProgress(0.25, detail = "Running ComBat batch correction")
            batch_correction <- sva::ComBat(
              dat = uncorrected_numeric_matrix,
              batch = batch_combat,
              mod = modcombat,
              par.prior = input$combat_parametric2
            )
            incProgress(0.25, detail = "Complete!")
          })
        }
      }else if (input$batch_correction_method2 == "Quantile Normalization") {
        withProgress(message = "Processing", value = 0, {
          incProgress(0.5, detail = "Quantile Normalization")
          uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
          uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
          batch_correction <- uncorrected_numeric_matrix
          incProgress(0.5, detail = "Complete!")
        })
        batch_correction
      }else if(input$batch_correction_method2 == "Mean Centering"){
        if (isTruthy(batch1_choices2())) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Running Mean Centering Batch Correction")
            mean_centering_batch = c(unlist(aligned_meta_file[,batch1_choices2()]))
            mean_centering_data = t(uncorrected_numeric_matrix)
            incProgress(0.25, detail = "Running Mean Centering Batch Correction")
            mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
            mean_center_correction <- as.data.frame(t(mean_center$xadj))
            batch_correction <- mean_center_correction
            incProgress(0.5, detail = "Complete!")
          })
        }
      } else if(input$batch_correction_method2 == "Harman"){
        req(harman_treatment_input2())
        req(batch1_choices2())
        withProgress(message = "Processing", value = 0, {
          incProgress(0.25, detail = "Running Harman Batch Correction")
          harman_correction_PCA <- Harman::harman(
            uncorrected_numeric_matrix,
            expt = aligned_meta_file[,harman_treatment_input2()],
            batch = aligned_meta_file[,batch1_choices2()],
            limit = input$HarmanCL2,
            printInfo = T,
            randseed = input$SeedSet
          )
          incProgress(0.25, detail = "Reconstructing Corrected Data")
          harman_correction <- Harman::reconstructData(harman_correction_PCA)
          batch_correction <- harman_correction
          incProgress(0.5, detail = "Complete!")
        })
      }else if(input$batch_correction_method2 == "RUVg"){
        if(input$RUVg_housekeeping_selection2 == "UserInput"){
          req(input$RUVg_user_control_genes2)
        }
        if (any(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.5, detail = "Running RUVg Batch Correction")
            RUVg_housekeeping_genes <- RUVg_housekeeping2()[which(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))]
            RUVg_matrix <- uncorrected_numeric_matrix
            RUVg_correction <- RUVSeq::RUVg(
              as.matrix(RUVg_matrix),
              cIdx = RUVg_housekeeping_genes,
              k = input$RUVg_estimate_factors2,
              drop = input$RUVg_drop_factors2,
              center = input$RUVg_mean_centered2,
              round = F,
              tolerance = input$RUVg_tolerance2,
              isLog = T
            )
            RUVg_correction_matrix <- as.data.frame(RUVg_correction$normalizedCounts)
            batch_correction <- RUVg_correction_matrix
            incProgress(0.5, detail = "Complete!")
          })
        }
      }else if(input$batch_correction_method2 == "SVA"){
        if (isTruthy(SVA_variable_of_interest_bc2())) {
          withProgress(message = "Processing", value = 0, {
            incProgress(0.25, detail = "Number of surrogate variables")
            expression_data <-  as.matrix(uncorrected_numeric_matrix)
            mod <-  model.matrix(reformulate(SVA_variable_of_interest_bc2()), data = aligned_meta_file)
            mod0 <- model.matrix(~1, data = aligned_meta_file)
            n.sv <- sva::num.sv(expression_data, mod, method = input$svaMethod_bc2)
            incProgress(0.25, detail = "Running SVA")
            svobj <- sva::sva(expression_data, mod , mod0, n.sv = n.sv)
            incProgress(0.25, detail = "Running Frozen SVA for batch correction")
            fsvaobj <- sva::fsva(expression_data, mod, svobj, expression_data)
            batch_correction <- fsvaobj$db
            incProgress(0.25, detail = "Complete!")
          })
        }
      }
    }
    
    
    if (exists("batch_correction")) {
      if (isTruthy(batch_correction)) {
        file_name <- paste(matrix2DlndTitle_react(),"_matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
        df <- as.data.frame(batch_correction)
        df <- cbind(data.frame(Genes = rownames(df)),
                    df)
        readr::write_tsv(df, file.path(temp_directory, file_name))
        batch_correction
      }
    }
    
  })
  
  #matrix2 <- shiny::reactive({
  #  
  #  #req(batch1_choices2())
  #  req(input$batch_correction_method2)
  #  uncorrected_numeric_matrix <- uncorrected_matrix()
  #  aligned_meta_file <- aligned_meta_file()
  #  
  #  if (input$batch_correction_method2 == "ComBatseq") {
  #    uncorrected_numeric_matrix <- matrix_raw_ForOut()
  #    
  #    if (is.null(input$covariate_choices_ComBatseq_raw) & is.null(input$covariate_choices_ComBatseq2)) {
  #      covar_raw <- 1
  #      covar_nor <- 1
  #    } else {
  #      covar_raw <- input$covariate_choices_ComBatseq_raw
  #      covar_nor <- input$covariate_choices_ComBatseq2
  #    }
  #    if (isTruthy(batch1_choices2())) {
  #      if ((covar_raw == covar_nor) & (input$batch1_choices_ComBatseq_raw == batch1_choices2())) {
  #        batch_correction <- matrix2_Raw()
  #      } else {
  #        ## ComBatSeq
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.5, detail = "Running ComBatseq Batch Correction")
  #          combatseq_corrected <- sva::ComBat_seq(
  #            as.matrix(uncorrected_numeric_matrix),
  #            batch = c(unlist(aligned_meta_file[,batch1_choices2()])),
  #            covar_mod = model_matrix()
  #          )
  #          batch_correction <- combatseq_corrected
  #          incProgress(0.5, detail = "Complete!")
  #        })
  #      }
  #    } else {
  #      batch_correction <- NULL
  #    }
  #    
  #    ## Normalize
  #    if (isTruthy(batch_correction)) {
  #      if (RawCountCheck()) {
  #        req(input$RawCountNorm)
  #        mat <- batch_correction
  #        mat_dgeList <- DGEList(counts = as.matrix(mat))
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.25, detail = "Normalizing Factors")
  #          mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
  #          incProgress(0.25, detail = "CPM")
  #          batch_correction <- edgeR::cpm(mat_dgeList_Norm)
  #          incProgress(0.75, detail = "Complete!")
  #        })
  #      }
  #      ## Log 
  #      if (input$Log_Choice_Raw) {
  #        batch_correction <- log2(batch_correction+1)
  #      }
  #      ## Quantile Normalize counts
  #      if (isTruthy(input$QuantNorm2)) {
  #        if (input$QuantNorm2){
  #          withProgress(message = "Processing", value = 0, {
  #            incProgress(0.5, detail = "Quantile Normalization")
  #            uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(batch_correction), keep.names = T)
  #            batch_correction <- as.data.frame(uncorrected_numeric_matrix_proc)
  #            incProgress(0.5, detail = "Complete!")
  #          })
  #        }
  #      }
  #    }
  #    
  #  } else {
  #    ## Normalize raw counts
  #    if (RawCountCheck()) {
  #      req(input$RawCountNorm)
  #      mat <- uncorrected_numeric_matrix
  #      mat_dgeList <- DGEList(counts = as.matrix(mat))
  #      withProgress(message = "Processing", value = 0, {
  #        incProgress(0.25, detail = "Normalizing Factors")
  #        mat_dgeList_Norm <- edgeR::calcNormFactors(mat_dgeList, method = input$RawCountNorm)
  #        incProgress(0.25, detail = "CPM")
  #        uncorrected_numeric_matrix <- edgeR::cpm(mat_dgeList_Norm)
  #        incProgress(0.75, detail = "Complete!")
  #      })
  #    }
  #    if (input$Log_Choice_Raw) {
  #      uncorrected_numeric_matrix <- log2(uncorrected_numeric_matrix+1)
  #    }
  #    ## Quantile Normalize counts
  #    if (isTruthy(input$QuantNorm2)) {
  #      if (input$QuantNorm2){
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.5, detail = "Quantile Normalization")
  #          uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
  #          uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
  #          incProgress(0.5, detail = "Complete!")
  #        })
  #      }
  #    }
  #    ## Correct
  #    if(input$batch_correction_method2 == "Uncorrected"){
  #      batch_correction <- uncorrected_numeric_matrix
  #    } else if(input$batch_correction_method2 == "Limma"){
  #      if (isTruthy(input$batch1_choices_limma2) & isTruthy(input$batch2_choices_limma2)) {
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.5, detail = "Running Limma batch correction")
  #          batch_correction <- limma::removeBatchEffect(
  #            uncorrected_numeric_matrix,
  #            batch = c(unlist(aligned_meta_file[,batch1_choices2()])),
  #            batch2 = c(unlist(aligned_meta_file[,batch2_choices2()])),
  #            covariates = model_matrix()
  #          )
  #          incProgress(0.5, detail = "Complete!")
  #        })
  #        
  #      }
  #    }else if(input$batch_correction_method2 == "ComBat"){
  #      if (isTruthy(batch1_choices2())) {
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.25, detail = "Running ComBat batch correction")
  #          batch_combat <- c(unlist(aligned_meta_file[,batch1_choices2()]))
  #          incProgress(0.25, detail = "Running ComBat batch correction")
  #          modcombat <-  stats::model.matrix(~1, data = as.data.frame(aligned_meta_file))
  #          incProgress(0.25, detail = "Running ComBat batch correction")
  #          batch_correction <- sva::ComBat(
  #            dat = uncorrected_numeric_matrix,
  #            batch = batch_combat,
  #            mod = modcombat,
  #            par.prior = input$combat_parametric2
  #          )
  #          incProgress(0.25, detail = "Complete!")
  #        })
  #        
  #      }
  #    }else if (input$batch_correction_method2 == "Quantile Normalization") {
  #      withProgress(message = "Processing", value = 0, {
  #        incProgress(0.5, detail = "Quantile Normalization")
  #        uncorrected_numeric_matrix_proc <- preprocessCore::normalize.quantiles(as.matrix(uncorrected_numeric_matrix), keep.names = T)
  #        uncorrected_numeric_matrix <- as.data.frame(uncorrected_numeric_matrix_proc)
  #        batch_correction <- uncorrected_numeric_matrix
  #        incProgress(0.5, detail = "Complete!")
  #      })
  #      
  #    }else if(input$batch_correction_method2 == "Mean Centering"){
  #      if (isTruthy(batch1_choices2())) {
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.25, detail = "Running Mean Centering Batch Correction")
  #          mean_centering_batch = c(unlist(aligned_meta_file[,batch1_choices2()]))
  #          mean_centering_data = t(uncorrected_numeric_matrix)
  #          incProgress(0.25, detail = "Running Mean Centering Batch Correction")
  #          mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
  #          mean_center_correction <- as.data.frame(t(mean_center$xadj))
  #          batch_correction <- mean_center_correction
  #          incProgress(0.5, detail = "Complete!")
  #        })
  #        
  #      }
  #    }else if(input$batch_correction_method2 == "ComBatseq"){
  #      if (isTruthy(batch1_choices2())) {
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.5, detail = "Running ComBatseq Batch Correction")
  #          combatseq_corrected <- sva::ComBat_seq(
  #            as.matrix(uncorrected_numeric_matrix),
  #            batch = c(unlist(aligned_meta_file[,batch1_choices2()])),
  #            covar_mod = model_matrix2()
  #          )
  #          batch_correction <- combatseq_corrected
  #          incProgress(0.5, detail = "Complete!")
  #        })
  #        
  #      }
  #    }else if(input$batch_correction_method2 == "Harman"){
  #      req(harman_treatment_input2())
  #      req(batch1_choices2())
  #      withProgress(message = "Processing", value = 0, {
  #        incProgress(0.25, detail = "Running Harman Batch Correction")
  #        harman_correction_PCA <- Harman::harman(
  #          uncorrected_numeric_matrix,
  #          expt = aligned_meta_file[,harman_treatment_input2()],
  #          batch = aligned_meta_file[,batch1_choices2()],
  #          limit = input$HarmanCL2,
  #          printInfo = T,
  #          randseed = input$SeedSet
  #        )
  #        incProgress(0.25, detail = "Reconstructing Corrected Data")
  #        harman_correction <- Harman::reconstructData(harman_correction_PCA)
  #        batch_correction <- harman_correction
  #        incProgress(0.5, detail = "Complete!")
  #      })
  #      
  #    }else if(input$batch_correction_method2 == "RUVg"){
  #      if(input$RUVg_housekeeping_selection2 == "UserInput"){
  #        req(input$RUVg_user_control_genes2)
  #      }
  #      if (any(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))) {
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.5, detail = "Running RUVg Batch Correction")
  #          RUVg_housekeeping_genes <- RUVg_housekeeping2()[which(RUVg_housekeeping2() %in% rownames(uncorrected_numeric_matrix))]
  #          RUVg_matrix <- uncorrected_numeric_matrix
  #          RUVg_correction <- RUVSeq::RUVg(
  #            as.matrix(RUVg_matrix),
  #            cIdx = RUVg_housekeeping_genes,
  #            k = input$RUVg_estimate_factors2,
  #            drop = input$RUVg_drop_factors2,
  #            center = input$RUVg_mean_centered2,
  #            round = F,
  #            tolerance = input$RUVg_tolerance2,
  #            isLog = T
  #          )
  #          RUVg_correction_matrix <- as.data.frame(RUVg_correction$normalizedCounts)
  #          batch_correction <- RUVg_correction_matrix
  #          incProgress(0.5, detail = "Complete!")
  #        })
  #        
  #      }
  #    }else if(input$batch_correction_method2 == "SVA"){
  #      if (isTruthy(SVA_variable_of_interest_bc2())) {
  #        withProgress(message = "Processing", value = 0, {
  #          incProgress(0.25, detail = "Number of surrogate variables")
  #          expression_data <-  as.matrix(uncorrected_numeric_matrix)
  #          mod <-  model.matrix(reformulate(SVA_variable_of_interest_bc2()), data = aligned_meta_file)
  #          mod0 <- model.matrix(~1, data = aligned_meta_file)
  #          n.sv <- sva::num.sv(expression_data, mod, method = input$svaMethod_bc2)
  #          incProgress(0.25, detail = "Running SVA")
  #          svobj <- sva::sva(expression_data, mod , mod0, n.sv = n.sv)
  #          incProgress(0.25, detail = "Running Frozen SVA for batch correction")
  #          fsvaobj <- sva::fsva(expression_data, mod, svobj, expression_data)
  #          batch_correction <- fsvaobj$db
  #          incProgress(0.25, detail = "Complete!")
  #        })
  #        
  #      }
  #    }
  #  }
  #  
  #  
  #  if (exists("batch_correction")) {
  #    if (isTruthy(batch_correction)) {
  #      file_name <- paste(gsub(" ","_",matrix2Title_react()),"_matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
  #      df <- as.data.frame(batch_correction)
  #      df <- cbind(data.frame(Genes = rownames(df)),
  #                  df)
  #      readr::write_tsv(df, file.path(temp_directory, file_name))
  #      batch_correction
  #    }
  #  }
  #  
  #})
  
  
  
  output$Matrix1title <- renderUI({
    req(matrix1Title_react())
    h3(matrix1Title_react())
  })
  output$rendMatrix1Head <- renderUI({
    req(matrix1Title_react())
    radioButtons("Matrix1Head",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  output$Matrix2title <- renderUI({
    req(matrix2Title_react())
    h3(matrix2Title_react())
  })
  output$rendMatrix2Head <- renderUI({
    req(matrix2Title_react())
    radioButtons("Matrix2Head",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  
  output$Matrix1titleRaw <- renderUI({
    req(matrix1_Raw())
    h3("Input RNAseq Count Matrix")
    #h3(paste(matrix1Title_react()),"RNAseq count Matrix")
  })
  output$rendMatrix1HeadRaw <- renderUI({
    req(matrix1_Raw())
    radioButtons("Matrix1HeadRaw",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  output$Matrix2titleRaw <- renderUI({
    req(matrix2_Raw())
    req(input$batch_correction_method2)
    if (input$batch_correction_method2 == "RUVg") {
      h3(paste(input$batch_correction_method2,"Corrected RNAseq Count Matrix"))
    }  else if (input$batch_correction_method2 == "Uncorrected") {
      h3(paste("Input RNAseq Count Matrix"))
    } else {
      h3(paste("ComBatseq Corrected RNAseq Count Matrix"))
    }
    #h3(paste(matrix2Title_react()),"RNAseq count Matrix")
  })
  output$rendMatrix2HeadRaw <- renderUI({
    req(matrix2_Raw())
    radioButtons("Matrix2HeadRaw",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  
  output$Matrix1titleZip <- renderUI({
    req(matrix1Title_react())
    if (input$RawCountInput) {
      if (matrix1Title_react() == "Uncorrected") {
        h3(paste(input$RawCountNorm,"Normalized Matrix"))
      } else if (matrix1Title_react() == "ComBatseq"){
        h3(paste(matrix1Title_react(),input$RawCountNorm,"Normalized Matrix"))
      } else {
        h3(paste(input$RawCountNorm,"Normalized",matrix1Title_react(),"Matrix"))
      }
    } else {
      h3(paste(matrix1Title_react()),"Matrix")
    }
  })
  output$rendMatrix1HeadZip <- renderUI({
    req(matrix1Title_react())
    radioButtons("Matrix1HeadZip",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  output$Matrix2titleZip <- renderUI({
    req(matrix2Title_react())
    if (input$RawCountInput) {
      if (matrix2Title_react() == "Uncorrected") {
        h3(paste(input$RawCountNorm,"Normalized Matrix"))
      } else if (matrix2Title_react() == "ComBatseq  Corrected"){
        h3(paste(matrix2Title_react(),input$RawCountNorm,"Normalized Matrix"))
      } else {
        h3(paste(input$RawCountNorm,"Normalized",matrix2Title_react(),"Matrix"))
      }
    } else {
      h3(paste(matrix2Title_react()),"Matrix")
    }
  })
  output$rendMatrix2HeadZip <- renderUI({
    req(matrix2Title_react())
    radioButtons("Matrix2HeadZip",NULL, choices = c("View table head","View entire table"), inline = T)
  })
  
  output$uncorrected_matrix_output <- DT::renderDataTable({
    req(matrix1())
    req(input$Matrix1Head)
    expr <- as.data.frame(matrix1())
    expr <- cbind(data.frame(Genes = rownames(expr)),
                  expr)
    if (input$Matrix1Head == "View table head") {
      expr <- head(expr,c(100,100))
    }
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  output$uncorrected_matrix_outputRaw <- DT::renderDataTable({
    req(matrix1_Raw())
    req(input$Matrix1HeadRaw)
    expr <- as.data.frame(matrix1_Raw())
    expr <- cbind(data.frame(Genes = rownames(expr)),
                  expr)
    if (input$Matrix1HeadRaw == "View table head") {
      expr <- head(expr,c(100,100))
    }
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  output$uncorrected_matrix_outputZip <- DT::renderDataTable({
    req(matrix1())
    req(input$Matrix1HeadZip)
    expr <- as.data.frame(matrix1())
    expr <- cbind(data.frame(Genes = rownames(expr)),
                  expr)
    if (input$Matrix1HeadZip == "View table head") {
      expr <- head(expr,c(100,100))
    }
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  
  output$corrected_matrix <- DT::renderDataTable({
    req(matrix2())
    req(input$Matrix2Head)
    expr <- as.data.frame(matrix2())
    expr <- cbind(data.frame(Genes = rownames(expr)),
                  expr)
    if (input$Matrix2Head == "View table head") {
      expr <- head(expr,c(100,100))
    }
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  output$corrected_matrixRaw <- DT::renderDataTable({
    req(matrix2_Raw())
    req(input$Matrix2HeadRaw)
    expr <- as.data.frame(matrix2_Raw())
    expr <- cbind(data.frame(Genes = rownames(expr)),
                  expr)
    if (input$Matrix2HeadRaw == "View table head") {
      expr <- head(expr,c(100,100))
    }
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  output$corrected_matrixZip <- DT::renderDataTable({
    req(matrix2())
    req(input$Matrix2HeadZip)
    expr <- as.data.frame(matrix2())
    expr <- cbind(data.frame(Genes = rownames(expr)),
                  expr)
    if (input$Matrix2HeadZip == "View table head") {
      expr <- head(expr,c(100,100))
    }
    DT::datatable(expr,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
  })
  matrix1DlndTitle_react <- reactive({
    
    req(input$batch_correction_method)
    logqn <- NULL
    qn <- NULL
    if (isTruthy(input$QuantNorm)) {
      if (input$QuantNorm) {
        qn <- paste("QuantNorm")
      }
    }
    if (RawCountCheck()) {
      if (input$Log_Choice_Raw) {
        if (isTruthy(input$QuantNorm)) {
          if (input$QuantNorm & (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("Log2Trans_QuantNorm")
          } else {
            logqn <- paste("Log2Trans")
          }
        } else {
          logqn <- paste("Log2Trans")
        }
      } else {
        logqn <- qn
      }
    } else {
      if (input$Log_Choice) {
        if (isTruthy(input$QuantNorm)) {
          if (input$QuantNorm & (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("Log2Trans_QuantNorm")
          } else {
            logqn <- paste("Log2Trans")
          }
        } else {
          logqn <- paste("Log2Trans")
        }
      } else {
        logqn <- qn
      }
    }
    
    if (input$RawCountInput) {
      if (input$batch_correction_method == "Uncorrected") {
        paste("RNAseqCounts",input$RawCountNorm,logqn,sep = "_")
      } else if (input$batch_correction_method %in% c("ComBatseq","RUVg")){
        paste("RNAseqCounts", input$batch_correction_method,input$RawCountNorm,logqn,sep = "_")
      } else {
        paste("RNAseqCounts",input$RawCountNorm,logqn,input$batch_correction_method,sep = "_")
      }
    } else {
      paste("InputData",logqn,input$batch_correction_method,sep = "_")
    }
    
  })
  matrix1Title_react <- reactive({
    req(input$batch_correction_method)
    logqn <- NULL
    qn <- NULL
    if (isTruthy(input$QuantNorm)) {
      if (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization")) {
        if (input$QuantNorm) {
          qn <- paste("→ QuantNorm")
        }
      }
    }
    if (RawCountCheck()) {
      if (input$Log_Choice_Raw) {
        if (isTruthy(input$QuantNorm)) {
          if (input$QuantNorm & (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("→ Log2Trans → QuantNorm")
          } else {
            logqn <- paste("→ Log2Trans")
          }
        } else {
          logqn <- paste("→ Log2Trans")
        }
      } else {
        logqn <- qn
      }
    } else {
      if (input$Log_Choice) {
        if (isTruthy(input$QuantNorm)) {
          if (input$QuantNorm & (!input$batch_correction_method %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("Log2Trans → QuantNorm")
          } else {
            logqn <- paste("Log2Trans")
          }
        } else {
          logqn <- paste("Log2Trans")
        }
      } else {
        logqn <- qn
      }
    }
    
    if (input$RawCountInput) {
      if (input$batch_correction_method == "Uncorrected") {
        paste("RNAseq Counts",input$RawCountNorm,"Normalized",logqn)
      } else if (input$batch_correction_method %in% c("ComBatseq","RUVg")){
        paste("RNAseq Counts", input$batch_correction_method ,"Corrected →",input$RawCountNorm,"Normaliztion",logqn)
      } else {
        paste("RNAseq Counts",input$RawCountNorm,"Normalized",logqn,"→",input$batch_correction_method,"Corrected")
      }
    } else {
      if (input$batch_correction_method == "Uncorrected") {
        paste("Input Data",logqn,input$batch_correction_method)
      } else {
        paste("Input Data",logqn,"→",input$batch_correction_method,"Corrected")
      }
    }
    
    
  })
  matrix2Title_react <- reactive({
    req(input$batch_correction_method2)
    
    logqn <- NULL
    qn <- NULL
    if (isTruthy(input$QuantNorm2)) {
      if (input$QuantNorm2) {
        qn <- paste("→ QuantNorm")
      }
    }
    if (RawCountCheck()) {
      if (input$Log_Choice_Raw) {
        if (isTruthy(input$QuantNorm2)) {
          if (input$QuantNorm2 & (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("→ Log2Trans → QuantNorm")
          } else {
            logqn <- paste("→ Log2Trans")
          }
        } else {
          logqn <- paste("→ Log2Trans")
        }
      } else {
        logqn <- qn
      }
    } else {
      if (input$Log_Choice) {
        if (isTruthy(input$QuantNorm2)) {
          if (input$QuantNorm2 & (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("Log2Trans → QuantNorm")
          } else {
            logqn <- paste("Log2Trans")
          }
        } else {
          logqn <- paste("Log2Trans")
        }
      } else {
        logqn <- qn
      }
    }
    
    if (input$RawCountInput) {
      if (input$batch_correction_method2 == "Uncorrected") {
        paste("RNAseq Counts",input$RawCountNorm,"Normalized",logqn)
      } else if (input$batch_correction_method2 %in% c("ComBatseq","RUVg")){
        paste("RNAseq Counts", input$batch_correction_method2 ,"Corrected →",input$RawCountNorm,"Normaliztion",logqn)
      } else {
        paste("RNAseq Counts",input$RawCountNorm,"Normalized",logqn,"→",input$batch_correction_method2,"Corrected")
      }
    } else {
      if (input$batch_correction_method2 == "Uncorrected") {
        paste("Input Data",logqn,input$batch_correction_method2)
      } else {
        paste("Input Data",logqn,"→",input$batch_correction_method2,"Corrected")
      }
    }
  })
  matrix2DlndTitle_react <- reactive({
    
    req(input$batch_correction_method2)
    logqn <- NULL
    qn <- NULL
    if (isTruthy(input$QuantNorm2)) {
      if (input$QuantNorm2) {
        qn <- paste("QuantNorm")
      }
    }
    if (RawCountCheck()) {
      if (input$Log_Choice_Raw) {
        if (isTruthy(input$QuantNorm2)) {
          if (input$QuantNorm2 & (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("Log2Trans_QuantNorm")
          } else {
            logqn <- paste("Log2Trans")
          }
        } else {
          logqn <- paste("Log2Trans")
        }
      } else {
        logqn <- qn
      }
    } else {
      if (input$Log_Choice) {
        if (isTruthy(input$QuantNorm2)) {
          if (input$QuantNorm2 & (!input$batch_correction_method2 %in% c("Uncorrected","Quantile Normalization"))) {
            logqn <- paste("Log2Trans_QuantNorm")
          } else {
            logqn <- paste("Log2Trans")
          }
        } else {
          logqn <- paste("Log2Trans")
        }
      } else {
        logqn <- qn
      }
    }
    
    if (input$RawCountInput) {
      if (input$batch_correction_method2 == "Uncorrected") {
        paste("RNAseqCounts",input$RawCountNorm,logqn,sep = "_")
      } else if (input$batch_correction_method2 %in% c("ComBatseq","RUVg")){
        paste("RNAseqCounts", input$batch_correction_method2,input$RawCountNorm,logqn,sep = "_")
      } else {
        paste("RNAseqCounts",input$RawCountNorm,logqn,input$batch_correction_method2,sep = "_")
      }
    } else {
      paste("InputData",logqn,input$batch_correction_method2,sep = "_")
    }
    
  })
  
  # Matrix ---------------------------------------------------------------------
  
  #observe({
  #  if (input$RawCountInput) {
  #    updateSelectInput(session,"batch_correction_method2", selected = "ComBatseq")
  #  }
  #})
  
  output$rendMatrixTabs <- renderUI({
    
    if (input$RawCountInput) {
      #if (input$RawCountNorm != "none") {
      tabsetPanel(id = "RawMatTabs",
                  tabPanel("RNAseq Count Matrix Preview",
                           shiny::fluidRow(
                             shiny::column(6,
                                           p(),
                                           shiny::uiOutput("Matrix1titleRaw"),
                                           uiOutput("rendMatrix1HeadRaw"),
                                           DT::dataTableOutput("uncorrected_matrix_outputRaw"),
                                           p(),
                                           shiny::downloadButton("dnldsave_uncorrected_matrixRaw","Dowload Single Table")
                             ),
                             shiny::column(6,
                                           p(),
                                           shiny::uiOutput("Matrix2titleRaw"),
                                           uiOutput("rendMatrix2HeadRaw"),
                                           DT::dataTableOutput("corrected_matrixRaw"),
                                           p(),
                                           shiny::downloadButton("dnldsave_corrected_matrixRaw","Dowload Single Table")
                             )
                           ),
                           value = 1
                  ),
                  tabPanel("RNAseq Count Matrix Transformed for Correction",
                           fluidRow(
                             column(6,
                                    verbatimTextOutput("FileCheckAlertsL")
                                    ),
                             column(6,
                                    verbatimTextOutput("FileCheckAlertsR")
                                    )
                           ),
                           #tabPanel(paste0(input$RawCountNorm," Transformed Matrix"),
                           shiny::fluidRow(
                             shiny::column(6,
                                           p(),
                                           shiny::uiOutput("Matrix1title"),
                                           uiOutput("rendMatrix1Head"),
                                           DT::dataTableOutput("uncorrected_matrix_output"),
                                           p(),
                                           shiny::downloadButton("dnldsave_uncorrected_matrix","Dowload Single Table")
                             ),
                             shiny::column(6,
                                           p(),
                                           shiny::uiOutput("Matrix2title"),
                                           uiOutput("rendMatrix2Head"),
                                           DT::dataTableOutput("corrected_matrix"),
                                           p(),
                                           shiny::downloadButton("dnldsave_corrected_matrix","Dowload Single Table")
                             )
                           ),
                           value = 2
                  )
      )
    } else {
      shiny::fluidRow(
        shiny::column(6,
                      p(),
                      verbatimTextOutput("FileCheckAlertsL2"),
                      shiny::uiOutput("Matrix1title"),
                      uiOutput("rendMatrix1Head"),
                      DT::dataTableOutput("uncorrected_matrix_output"),
                      p(),
                      shiny::downloadButton("dnldsave_uncorrected_matrix","Dowload Single Table")
        ),
        shiny::column(6,
                      p(),
                      verbatimTextOutput("FileCheckAlertsR2"),
                      shiny::uiOutput("Matrix2title"),
                      uiOutput("rendMatrix2Head"),
                      DT::dataTableOutput("corrected_matrix"),
                      p(),
                      shiny::downloadButton("dnldsave_corrected_matrix","Dowload Single Table")
        )
      )
    }
    
  })
  
  
  # PCA ------------------------------------------------------------------------
  
  output$PCA1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"PCA")
  })
  output$PCA2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"PCA")
  })
  
  # Generating PCA plot
  output$rendbatch_choices_PCA <- shiny::renderUI({
    shiny::selectInput("batch_choices_PCA",
                       "Group Color",
                       choices = c("Select", batch_names_from_meta()[-1]),
                       selected = batch1_choices())
  })
  
  output$rendPCAhover <- shiny::renderUI({
    batch_choice <- input$batch_choices_PCA
    shiny::selectInput("PCAhover",
                       "Information to Display on Hover:",
                       choices = unlist(strsplit(batch_names_from_meta(), ",")),
                       selected = c(batch_names_from_meta()[1],batch_choice),
                       multiple = T)
  })
  
  ## Matrix 1 -----------------------------------------------------------------
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
  
  PCA_data <- shiny::reactive({
    req(matrix1())
    PCA_data <- cbind((t(matrix1())), aligned_meta_file())
    PCA_data
  })
  
  uncorrected_PCA_react <- reactive({
    req(matrix1())
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      #withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "PCA Object")
      if (input$PCA_type == "Cluster Annotation"){
        mat <- matrix1()
        row.names(mat) <- NULL
        mat <- as.data.frame(mat)
        mat <- mat[, sapply(mat, var) != 0]
        mat <- as.matrix(mat)
        pca <- cluster::pam(as.data.frame(t(mat)),input$cluster_number)
      } else if (input$PCA_type == "Meta Annotation"){
        mat <- matrix1()
        mat_check <- apply(mat, 1, function(a) length(unique(a))==1)
        mat_check_F <- names(mat_check[which(mat_check==FALSE)])
        mat2 <- mat[mat_check_F,]
        pca <- stats::prcomp(as.data.frame(t(mat2)), scale. = TRUE)
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
        frame = input$FrameClusters,
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
    
    if (input$PCA_type == "Cluster Annotation"){
      p <- ggplot2::autoplot(
        uncorrected_PCA_react(),
        frame = input$FrameClusters,
        frame.type = 'norm')
      ply <- ggplotly(p)
      for (i in 1:length(ply$x$data)){
        if (!is.null(ply$x$data[[i]]$name)){
          ply$x$data[[i]]$name =  gsub("\\(","",str_split(ply$x$data[[i]]$name,",")[[1]][1])
          if (ply$x$data[[i]]$hoveron == "fills") {
            ply$x$data[[i]]$showlegend = FALSE
          }
        }
      }
      ply
    } else if (input$PCA_type == "Meta Annotation"){
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
    }
    
  })
  
  output$uncorrected_PCA <- renderPlotly({
    uncorrected_PCA_plot_react()
  })
  
  ## Matrix 2 -----------------------------------------------------------------
  
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
    req(matrix2())
    PCA_data_corrected <- cbind((t(matrix2())), aligned_meta_file())
    PCA_data_corrected
  })
  
  corrected_PCA_react <- shiny::reactive({
    req(matrix2())
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      #withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "PCA Object")
      if (input$PCA_type == "Cluster Annotation"){
        mat <- matrix2()
        row.names(mat) <- NULL
        mat <- as.data.frame(mat)
        mat <- mat[, sapply(mat, var) != 0]
        mat <- as.matrix(mat)
        k <- input$cluster_number
        pca <- cluster::pam(as.data.frame(t(mat)),input$cluster_number)
        #pca <- cluster::pam(as.data.frame(t(batch_correction())),input$cluster_number)
      } else if (input$PCA_type == "Meta Annotation"){
        mat <- matrix2()
        mat_check <- apply(mat, 1, function(a) length(unique(a))==1)
        mat_check_F <- names(mat_check[which(mat_check==FALSE)])
        mat2 <- mat[mat_check_F,]
        pca <- stats::prcomp(as.data.frame(t(mat2)), scale. = TRUE)
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
        frame = input$FrameClusters,
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
    
    if (input$PCA_type == "Cluster Annotation"){
      p <- ggplot2::autoplot(
        corrected_PCA_react(),
        frame = input$FrameClusters,
        frame.type = 'norm')
      ply <- ggplotly(p)
      for (i in 1:length(ply$x$data)){
        if (!is.null(ply$x$data[[i]]$name)){
          ply$x$data[[i]]$name =  gsub("\\(","",str_split(ply$x$data[[i]]$name,",")[[1]][1])
          if (ply$x$data[[i]]$hoveron == "fills") {
            ply$x$data[[i]]$showlegend = FALSE
          }
        }
      }
      ply
    } else if (input$PCA_type == "Meta Annotation"){
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
    }
    
    
  })
  output$corrected_PCA <- renderPlotly({
    corrected_PCA_plot_react()
  })
  
  
  # PCA MC --------------------------------------------------------------------
  
  output$PCAMC1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"PCA Multiple Components")
  })
  output$PCAMC2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"PCA Multiple Components")
  })
  output$PCA_mc_color_choice <- shiny::renderUI({
    shiny::selectInput("PCA_mc_color_choice",
                       "Group color",
                       c("Select", unlist(strsplit(batch_names_from_meta(), ","))))
  })
  
  ## Matrix 1 -----------------------------------------------------------------
  
  # Generating multiple PC plot
  
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
    req(matrix1())
    mat_PCA_mc <- matrix1()
    names(mat_PCA_mc) <- NULL
    uncorrected_PCA_mc_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_PCA_mc)),
      colData = aligned_meta_file(),
      rowData = rownames(matrix1())
    )
    SummarizedExperiment::assay(uncorrected_PCA_mc_SCE, "logcounts") <- SingleCellExperiment::counts(uncorrected_PCA_mc_SCE)
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
  
  ## Matrix 2 -----------------------------------------------------------------
  
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
    req(matrix2())
    mat_PCA_mc <- matrix2()
    names(mat_PCA_mc) <- NULL
    corrected_PCA_mc_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_PCA_mc)),
      colData = aligned_meta_file(),
      rowData = rownames(matrix2())
    )
    SummarizedExperiment::assay(corrected_PCA_mc_SCE, "logcounts") <- SingleCellExperiment::counts(corrected_PCA_mc_SCE)
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
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
  
  
  # PCA Details ---------------------------------------------------------------
  
  output$PCADet1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"Corrected PCA Details")
  })
  output$PCADet2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"Corrected PCA Details")
  })
  
  ## Matrix 1 -----------------------------------------------------------------
  
  # Uncorrected PCA details plots
  uncorrected_PCA_details <- shiny::reactive({
    #req(matrix1())
    if (is.matrix(matrix1())) {
      withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
        incProgress(0.5, detail = "Running FactoMineR PCA")
        pc_obj <- FactoMineR::PCA(t(matrix1()), graph = F)
        incProgress(0.5, detail = "Complete!")
      })
      pc_obj
    }
  })
  uncorrected_PCA_details2 <- shiny::reactive({
    req(matrix1())
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      incProgress(0.5, detail = "Running stats prcomp PCA")
      pc_obj <- stats::prcomp(t(matrix1()))
      incProgress(0.5, detail = "Complete!")
    })
    pc_obj
  })
  uncorrected_scree_plot_react <- shiny::reactive({
    req(uncorrected_PCA_details())
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
    req(uncorrected_scree_plot_react())
    uncorrected_scree_plot_react()
  })
  output$PCA_factors_choices <- shiny::renderUI({
    shiny::selectInput("PCA_factors_choices",
                       "Select Factor for PCA Details",
                       c(batch_names_from_meta()[-1]))
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
    req(uncorrected_PCA_details())
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
    req(uncorrected_PCA_details())
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
  
  ## Matrix 2 -----------------------------------------------------------------
  
  #Corrected PCA Details plots
  corrected_PCA_details <- shiny::reactive({
    if (is.matrix(matrix2())) {
      withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
        incProgress(0.5, detail = "Running FactoMineR PCA")
        pc_obj <- FactoMineR::PCA(as.data.frame(t(matrix2())), graph = F)
        incProgress(0.5, detail = "Complete!")
      })
      pc_obj
    }
  })
  corrected_PCA_details2 <- shiny::reactive({
    req(matrix2())
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "Running stats prcomp PCA")
      pc_obj <- stats::prcomp(as.data.frame(t(matrix2())))
      incProgress(0.5, detail = "Complete!")
    })
    pc_obj
  })
  corrected_scree_plot_react <- shiny::reactive({
    req(corrected_PCA_details())
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
    req(corrected_scree_plot_react())
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
    req(corrected_PCA_details())
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
    req(corrected_PCA_details())
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
  
  # UMAP ----------------------------------------------------------------------
  
  output$UMAP1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"UMAP")
  })
  output$UMAP2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"UMAP")
  })
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
    mat <- matrix1()
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
  UMAP_Feature_Choices <- reactive({
    
    req(aligned_meta_file())
    Features <- NULL
    FeatCat <- input$UMAPFeatureCategory
    if (FeatCat == "Matrix Features") {
      Features <- rownames(matrix1())
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
    } else {
      Features <- NULL
      Features
    }
    
  })
  
  shiny::observe({
    
    updateSelectizeInput(session = session, inputId = "UMAPFeatSelection",
                         choices = UMAP_Feature_Choices(),
                         server = T)
    
  })
  
  UMAPGeneSetTableBack_react <- reactive({
    
    GeneSetTable_sub <- geneset_df[which(geneset_df[,1] == input$UMAPGeneSetCat),]
    GeneSetTable_sub <- GeneSetTable_sub[,-c(1,2)]
    colnames(GeneSetTable_sub) <- c("Gene Set Category","Gene Set")
    GeneSetTable_sub
    
  })
  
  output$GeneSetTableUIUMAP <- DT::renderDataTable({
    
    GeneSetTable_sub <- UMAPGeneSetTableBack_react()
    DT::datatable(GeneSetTable_sub,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  ## Matrix 1 -----------------------------------------------------------------
  
  UMAP_PCA_Proj_uncorr_Samples <- reactive({
    
    if (input$uncorrected_panel == "pca_main") {
      if (input$PCA_main_pan == "umap") {
        withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
          incProgress(0.5, detail = "Running PCA")
          mat_uncorr <- matrix1()
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
  
  UMAP_uncorr_anno_df <- shiny::reactive({
    
    FeatCat <- input$UMAPFeatureCategory
    Feature <- input$UMAPFeatSelection
    meta <- aligned_meta_file()
    NameCol <- colnames(meta)[1]
    if (FeatCat == "Matrix Features") {
      mat <- matrix1()
      featdf <- featdf[Feature,]
      #featdf <- mat[which(mat[,1] == Feature),]
      #rownames(featdf) <- featdf[,1]
      #featdf <- featdf[,-1]
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
      gs_name <- as.character(UMAPGeneSetTableBack_react()[input$GeneSetTableUIUMAP_rows_selected,ncol(UMAPGeneSetTableBack_react())])
      if (isTruthy(gs_name)) {
        gs <- geneset[gs_name]
        Feature <- gs_name
        mat <- Mat_for_ssGSEA_uncorr()
        #rownames(mat) <- mat[,1]
        #mat <- mat[,-1]
        mat <- as.matrix(mat)
        withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
          incProgress(0.5, detail = "Running ssGSEA")
          ssGSEA_param <- GSVA::ssgseaParam(mat,gs)
          ssGSEA <- GSVA::gsva(ssGSEA_param)
          incProgress(0.5, detail = "Complete!")
        })
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA[,NameCol] <- rownames(ssGSEA)
        meta <- merge(meta,ssGSEA)
      }
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
        withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
    
    #req(input$UMAPFeatSelection)
    FeatCat <- input$UMAPFeatureCategory
    plot_df <- uncorrected_umap_coord()
    rownames(plot_df) <- plot_df[,1]
    UMAPdotSize <- 2
    if (FeatCat == "Gene Set Pathways") {
      umap_annoCol <- as.character(UMAPGeneSetTableBack_react()[input$GeneSetTableUIUMAP_rows_selected,ncol(UMAPGeneSetTableBack_react())])
    } else {
      umap_annoCol <- input$UMAPFeatSelection
    }
    req(umap_annoCol)
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
  
  ## Matrix 2 -----------------------------------------------------------------
  UMAP_ImmDeconv_corr_react <- reactive({
    
    req(input$UMAPImmuneDeconvMethods)
    deconvMethod <- input$UMAPImmuneDeconvMethods
    mat <- matrix2()
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
        withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
          incProgress(0.5, detail = "Running PCA")
          mat_uncorr <- matrix2()
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
      mat <- matrix2()
      featdf <- mat[Feature,]
      #featdf <- mat[which(mat[,1] == Feature),]
      #rownames(featdf) <- featdf[,1]
      #featdf <- featdf[,-1]
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
    }else if (FeatCat == "Gene Set Pathways") {
      gs_name <- as.character(UMAPGeneSetTableBack_react()[input$GeneSetTableUIUMAP_rows_selected,ncol(UMAPGeneSetTableBack_react())])
      if (isTruthy(gs_name)) {
        gs <- geneset[gs_name]
        Feature <- gs_name
        mat <- Mat_for_ssGSEA_corr()
        #rownames(mat) <- mat[,1]
        #mat <- mat[,-1]
        mat <- as.matrix(mat)
        withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
          incProgress(0.5, detail = "Running ssGSEA")
          ssGSEA_param <- GSVA::ssgseaParam(mat,gs)
          ssGSEA <- GSVA::gsva(ssGSEA_param)
          incProgress(0.5, detail = "Complete!")
        })
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA[,NameCol] <- rownames(ssGSEA)
        meta <- merge(meta,ssGSEA)
      }
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
        withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
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
    
    #req(input$UMAPFeatSelection)
    FeatCat <- input$UMAPFeatureCategory
    plot_df <- corrected_umap_coord()
    rownames(plot_df) <- plot_df[,1]
    UMAPdotSize <- 2
    if (FeatCat == "Gene Set Pathways") {
      umap_annoCol <- as.character(UMAPGeneSetTableBack_react()[input$GeneSetTableUIUMAP_rows_selected,ncol(UMAPGeneSetTableBack_react())])
    } else {
      umap_annoCol <- input$UMAPFeatSelection
    }
    req(umap_annoCol)
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
  
  # Cluster -------------------------------------------------------------------
  
  output$Cluster1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"Cluster Plots")
  })
  output$Cluster2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"Cluster Plots")
  })
  
  observe({
    updateSelectInput(session,"ClusterMethodHeat",selected = input$ClusterMethod)
  })
  observe({
    updateSelectInput(session,"ClusterMethod",selected = input$ClusterMethodHeat)
  })
  
  ## Matrix 1 -----------------------------------------------------------------
  
  
  # uncorrected cluster analysis
  cluster_mv_features_uncorr_matrix <- reactive({
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$cluster_n_MV_features
      var_type <- input$VarianceMeasure
      mat <- matrix1()
      featColName <- colnames(mat)
      #rownames(mat) <- mat[,1]
      #mat <- mat[,-1]
      mad <- NULL
      var <- NULL
      cv <- NULL
      if (!input$Log_Choice) {
        mat <- log2(mat + 1)
      }
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(mat, 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(mat, 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(mat, 1, cv))
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = topN)
        out <- cbind(names(cv), cv[names(cv)], mat[names(cv),])
        colnames(out) <- c("Gene", "CV", colnames(mat))
        dataset <- mat[names(cv),]
      }
      
      
      zdataset <- t(apply(dataset, 1, scale))
      colnames(zdataset) <- names(dataset)
      dataset <- as.data.frame(zdataset)
      colnames(dataset) <- featColName
      
      #dataset[,featColName] <- rownames(dataset)
      #dataset <- dataset %>% dplyr::relocate(any_of(featColName))
      incProgress(0.5, detail = "Complete!")
    })
    dataset
    
  })
  uncorrected_elbow_analysis <- reactive({
    req(cluster_mv_features_uncorr_matrix())
    if (ncol(cluster_mv_features_uncorr_matrix()) <= 10) {
      k.max_var <- ncol(cluster_mv_features_uncorr_matrix())-1
    } else {
      k.max_var <- 10
    }
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      incProgress(0.5, detail = "Running Elbow Analysis")
      uncorrected_elbow_analysis <- fviz_nbclust(x = t(cluster_mv_features_uncorr_matrix()), kmeans, method = "wss",verbose = T, k.max = k.max_var)
      incProgress(0.5, detail = "Complete!")
    })
    uncorrected_elbow_analysis
  })
  output$uncorrected_elbow_plot <- renderPlot({
    uncorrected_elbow_analysis()
  })
  uncorrected_silhouette_analysis <- reactive({
    req(cluster_mv_features_uncorr_matrix())
    if (ncol(cluster_mv_features_uncorr_matrix()) <= 10) {
      k.max_var <- ncol(cluster_mv_features_uncorr_matrix())-1
    } else {
      k.max_var <- 10
    }
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      incProgress(0.5, detail = "Running Silhouette Analysis")
      uncorrected_silhouette_analysis <- fviz_nbclust(x = t(cluster_mv_features_uncorr_matrix()), kmeans, method = "silhouette",verbose = T, k.max = k.max_var)
      incProgress(0.5, detail = "Complete!")
    })
    uncorrected_silhouette_analysis
  })
  output$uncorrected_silhouette_plot <- renderPlot({
    uncorrected_silhouette_analysis()
  })
  uncorrected_dunn_index_analysis <- reactive({
    req(cluster_mv_features_uncorr_matrix())
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      incProgress(0.5, detail = "Running Dunn Index Analysis")
      if (ncol(cluster_mv_features_uncorr_matrix()) <= 10) {
        #k.max_var <- ncol(cluster_mv_features_uncorr_matrix())-1
        dunn_k <- c(2:ncol(cluster_mv_features_uncorr_matrix())-1)
      } else {
        dunn_k <- c(2:10)
      }
      #dunn_k <- c(2:10)
      dunnin <- c()
      for (i in dunn_k){
        dunnin[i] <- dunn(
          distance = dist(t(cluster_mv_features_uncorr_matrix())),
          clusters = kmeans(t(cluster_mv_features_uncorr_matrix()), i)$cluster
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
        scale_x_continuous(breaks=seq(1, max(uncorrected_dunn_index_analysis[,"cluster_number"]), 1))+
        theme(axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14))
      incProgress(0.5, detail = "Complete!")
    })
    p
    
  })
  output$uncorrected_dunn_index_plot <- renderPlot({
    uncorrected_dunn_index_analysis()
  })
  
  ## Matrix 2 -----------------------------------------------------------------
  
  cluster_mv_features_corr_matrix <- reactive({
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$cluster_n_MV_features
      var_type <- input$VarianceMeasure
      mat <- matrix2()
      featcolnames <- colnames(mat)
      #mat <- mat[,-1]
      
      mad <- NULL
      var <- NULL
      cv <- NULL
      if (!input$Log_Choice) {
        mat <- log2(mat + 1)
      }
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(mat, 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(mat, 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(mat, 1, cv))
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = topN)
        out <- cbind(names(cv), cv[names(cv)], mat[names(cv),])
        colnames(out) <- c("Gene", "CV", colnames(mat))
        dataset <- mat[names(cv),]
      }
      
      zdataset <- t(apply(dataset, 1, scale))
      colnames(zdataset) <- names(dataset)
      dataset <- as.data.frame(zdataset)
      colnames(dataset) <- featcolnames
      incProgress(0.5, detail = "Complete!")
    })
    
    dataset
    
  })
  corrected_elbow_analysis <- reactive({
    req(cluster_mv_features_corr_matrix())
    if (ncol(cluster_mv_features_uncorr_matrix()) <= 10) {
      k.max_var <- ncol(cluster_mv_features_uncorr_matrix())-1
    } else {
      k.max_var <- 10
    }
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "Running Elbow Analysis")
      corrected_elbow_analysis <- fviz_nbclust(x = t(cluster_mv_features_corr_matrix()), kmeans, method = "wss",verbose = T, k.max = k.max_var)
      incProgress(0.5, detail = "Complete!")
    })
    corrected_elbow_analysis
  })
  output$corrected_elbow_plot <- renderPlot({
    corrected_elbow_analysis()
  })
  corrected_silhouette_analysis <- reactive({
    req(cluster_mv_features_corr_matrix())
    if (ncol(cluster_mv_features_uncorr_matrix()) <= 10) {
      k.max_var <- ncol(cluster_mv_features_uncorr_matrix())-1
    } else {
      k.max_var <- 10
    }
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "Running Silhouette Analysis")
      corrected_silhouette_analysis <- fviz_nbclust(x = t(cluster_mv_features_corr_matrix()), kmeans, method = "silhouette",verbose = T, k.max = k.max_var)
      incProgress(0.5, detail = "Complete!")
    })
    corrected_silhouette_analysis
  })
  output$corrected_silhouette_plot <- renderPlot({
    corrected_silhouette_analysis()
  })
  corrected_dunn_index_analysis <- reactive({
    req(cluster_mv_features_corr_matrix())
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "Running Dunn Index Analysis")
      if (ncol(cluster_mv_features_uncorr_matrix()) <= 10) {
        #k.max_var <- ncol(cluster_mv_features_uncorr_matrix())-1
        dunn_k <- c(2:ncol(cluster_mv_features_uncorr_matrix())-1)
      } else {
        dunn_k <- c(2:10)
      }
      #dunn_k <- c(2:10)
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
        scale_x_continuous(breaks=seq(1, max(corrected_dunn_index_analysis[,"cluster_number"]), 1))+
        theme(axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14))
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$corrected_dunn_index_plot <- renderPlot({
    corrected_dunn_index_analysis()
  })
  
  # Heatmap -------------------------------------------------------------------
  
  output$Heatmap1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"Heatmap")
  })
  output$Heatmap2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"Heatmap")
  })
  
  output$rendHeatmapAnnoSel <- renderUI({
    
    shiny::selectInput("HeatmapAnnoSel",
                       "Column Annoation:",
                       c(batch_names_from_meta()[-1]),
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
  ## Matrix 1 -----------------------------------------------------------------
  
  # Generating heatmap of uncorrected data
  uncorrected_heatmap <- shiny::reactive({
    req(cluster_mv_features_uncorr_matrix())
    #req(heat_colAnn())
    
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      incProgress(0.5, detail = "Generating Heatmap")
      
      colAnn <- heat_colAnn()
      HeatRowNames <- ifelse("Turn on Row Names" %in% input$HeatRowColNames,TRUE,FALSE)
      HeatColNames <- ifelse("Turn on Column Names" %in% input$HeatRowColNames,TRUE,FALSE)
      clusterMethod <- input$ClusterMethodHeat
      
      uncorrected_matrix_heatmap <- as.matrix(cluster_mv_features_uncorr_matrix())
      uncorrected_matrix_heatmap_cols <- colnames(uncorrected_matrix_heatmap)
      #uncorrected_matrix_heatmap <- as.matrix(uncorrected_numeric_matrix()[,-1])
      #uncorrected_matrix_heatmap <- as.matrix(log2(uncorrected_numeric_matrix()[,-1] + 1))
      #uncorrected_matrix_scaled <- t(apply(uncorrected_matrix_heatmap, 1, scale))
      uncorrected_matrix_scaled <- uncorrected_matrix_heatmap
      colnames(uncorrected_matrix_scaled) <- uncorrected_matrix_heatmap_cols
      p <- suppressMessages(ComplexHeatmap::Heatmap(uncorrected_matrix_scaled, top_annotation = colAnn,
                                                    show_row_names = HeatRowNames, show_column_names = HeatColNames,
                                                    heatmap_legend_param = list(title = "Expression"),
                                                    clustering_method_columns = clusterMethod))
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$uncorrected_heatmap <- shiny::renderPlot({
    uncorrected_heatmap()
  })
  
  
  ## Matrix 2 -----------------------------------------------------------------
  # Generating heatmap of uncorrected data
  corrected_heatmap <- shiny::reactive({
    req(cluster_mv_features_corr_matrix())
    #req(heat_colAnn())
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "Generating Heatmap")
      colAnn <- heat_colAnn()
      HeatRowNames <- ifelse("Turn on Row Names" %in% input$HeatRowColNames,TRUE,FALSE)
      HeatColNames <- ifelse("Turn on Column Names" %in% input$HeatRowColNames,TRUE,FALSE)
      clusterMethod <- input$ClusterMethodHeat
      
      corrected_matrix_heatmap <- as.matrix(cluster_mv_features_corr_matrix())
      corrected_matrix_heatmap_cols <- colnames(corrected_matrix_heatmap)
      #corrected_matrix_scaled <- t(apply(corrected_matrix_heatmap, 1, scale))
      corrected_matrix_scaled <- corrected_matrix_heatmap
      colnames(corrected_matrix_scaled) <- corrected_matrix_heatmap_cols
      p <- suppressMessages(ComplexHeatmap::Heatmap(corrected_matrix_scaled, top_annotation = colAnn,
                                                    show_row_names = HeatRowNames, show_column_names = HeatColNames,
                                                    heatmap_legend_param = list(title = "Expression"),
                                                    clustering_method_columns = clusterMethod))
      incProgress(0.5, detail = "Complete!")
    })
    p
  })
  output$corrected_heatmap <- shiny::renderPlot({
    corrected_heatmap()
  })
  
  
  # Diversity Eval ------------------------------------------------------------
  
  output$DivEval1title <- renderUI({
    req(matrix1Title_react())
    h4(paste(matrix1Title_react()),input$ClusterMethod,"Cluster Bar Plot Colored by",input$BarPFillCol)
  })
  output$DivEval2title <- renderUI({
    req(matrix2Title_react())
    h4(paste(matrix2Title_react()),input$ClusterMethod,"Cluster Bar Plot Colored by",input$BarPFillCol)
  })
  output$DivEvalHEAvg1title <- renderUI({
    req(matrix1Title_react())
    h4(paste(matrix1Title_react()),"Average Heterogeneity and Eveness")
  })
  output$DivEvalHEAvg2title <- renderUI({
    req(matrix2Title_react())
    h4(paste(matrix2Title_react()),"Average Heterogeneity and Eveness")
  })
  output$DivEvalHE1title <- renderUI({
    req(matrix1Title_react())
    h4(paste(matrix1Title_react()),"Heterogeneity and Eveness")
  })
  output$DivEvalHE2title <- renderUI({
    req(matrix2Title_react())
    h4(paste(matrix2Title_react()),"Heterogeneity and Eveness")
  })
  
  #output$rendErrorMessageUncorr <- renderUI({
  #  if (!isTruthy(batch1_choices())) {
  #    #if (input$batch_correction_method == "Uncorrected") {
  #    span(textOutput("ErrorMessageUncorr"), style="color:red")
  #  }
  #})
  #output$rendErrorMessageCorr <- renderUI({
  #  if (!isTruthy(batch1_choices())) {
  #    #if (input$batch_correction_method == "Uncorrected") {
  #    span(textOutput("ErrorMessageCorr"), style="color:red")
  #  }
  #})
  #output$ErrorMessageUncorr <- renderText({"Select batch correction method and batch to continue"})
  #output$ErrorMessageCorr <- renderText({"Select batch correction method and batch to continue"})
  
  output$rendBarPFillCol <- renderUI({
    
    shiny::selectInput("BarPFillCol",
                       "Group Color",
                       choices = batch_names_from_meta()[-1],
                       selected = batch1_choices())
    
  })
  ClusterInfoTab_react <- reactive({
    
    dataset_uncorr <- cluster_mv_features_uncorr_matrix()
    dataset_corr <- cluster_mv_features_corr_matrix()
    meta <- aligned_meta_file()
    method <- input$ClusterMethod
    NumClusters <- input$NumOfCluster
    #corrMethod <- gsub(" ","_",input$batch_correction_method)
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    corrMethod2 <- gsub(" ","_",matrix2Title_react())
    
    corrMethod2 <- ifelse(corrMethod1==corrMethod2,paste0(corrMethod2,"_2"),corrMethod2)
    
    hclust_uncorr <- hclust(dist(t(dataset_uncorr)), method = method)
    hclust_uncorr_cut <- sort(cutree(hclust_uncorr, k=NumClusters))
    hclust_uncorr_cut_df <- data.frame(names(hclust_uncorr_cut),
                                       hclust_uncorr_cut)
    colnames(hclust_uncorr_cut_df) <- c(colnames(meta)[1],
                                        paste0(method,"_Cluster_",corrMethod1))
    meta <- merge(meta,hclust_uncorr_cut_df,all = T, sort = F)
    
    #rownames(dataset_corr) <- dataset_corr[,1]
    #dataset_corr <- dataset_corr[,-1]
    hclust_corr <- hclust(dist(t(dataset_corr)), method = method)
    hclust_corr_cut <- sort(cutree(hclust_corr, k=NumClusters))
    hclust_corr_cut_df <- data.frame(names(hclust_corr_cut),
                                     hclust_corr_cut)
    colnames(hclust_corr_cut_df) <- c(colnames(meta)[1],
                                      paste0(method,"_Cluster_",corrMethod2))
    meta <- merge(meta,hclust_corr_cut_df,all = T, sort = F)
    boxplot_meta_file(meta)
    meta
    
  })
  
  DiversityEval_batch <- reactive({
    
    batch <- input$BarPFillCol
    batch
    
  })
  ## Matrix 1 -----------------------------------------------------------------
  
  uncorr_bar_plot_react <- reactive({
    
    req(input$BarPFillCol)
    plot_df <- ClusterInfoTab_react()
    batch <- input$BarPFillCol
    #batch <- batch1_choices()
    method <- input$ClusterMethod
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    cluster_col <- paste0(method,"_Cluster_",corrMethod1)
    FillChoice <- input$barPfill
    plot_df[,batch] <- as.factor(plot_df[,batch])
    tickFont <- input$barPAxisTkSize
    titleFont <- input$barPAxisTtSize
    barp <- ggplot(plot_df, aes(fill = !!sym(batch), x = !!sym(cluster_col)))
    if (FillChoice) {
      barp <- barp + geom_bar(position = "fill") +
        theme_minimal()
    } else {
      barp <- barp + geom_bar() +
        theme_minimal()
    }
    barp <- barp +
      xlab(paste(method,"Cluster"))
    barp <- barp + theme(axis.text.x = element_text(size = tickFont),
                         axis.title.x = element_text(size = titleFont),
                         axis.text.y = element_text(size = tickFont),
                         axis.title.y = element_text(size = titleFont))
    barp
    
  })
  
  output$uncorr_bar_plot <- renderPlot({
    
    barp <- uncorr_bar_plot_react()
    barp
    
  })
  
  output$ClusterInfoTab <- DT::renderDataTable({
    
    meta <- ClusterInfoTab_react()
    method <- input$ClusterMethod
    BarPCol <- input$BarPFillCol
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    corrMethod2 <- gsub(" ","_",matrix2Title_react())
    corrMethod2 <- ifelse(corrMethod1==corrMethod2,paste0(corrMethod2,"_2"),corrMethod2)
    #corrMethod <- gsub(" ","_",input$batch_correction_method)
    meta <- meta %>%
      dplyr::relocate(any_of(c(paste0(method,"_Cluster_",corrMethod1),paste0(method,"_Cluster_",corrMethod2),BarPCol)),
                      .after = colnames(meta)[1])
    DT::datatable(meta,
                  options = list(lengthMenu = c(5,10,20,100,1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
    
  })
  
  uncorr_het_react <- reactive({
    
    #req(!isTruthy(batch1_choices()))
    req(DiversityEval_batch())
    meta <- ClusterInfoTab_react()
    method <- input$ClusterMethod
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    uncorr_col <- paste0(method,"_Cluster_",corrMethod1)
    batch <- DiversityEval_batch()
    
    uncorr_cluster_batch_freq <- as.data.frame.matrix(table(meta[,uncorr_col], meta[,batch]))
    
    uncorrected_heterogeneity_shannon = heterogeneity(uncorr_cluster_batch_freq, method = "shannon") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    uncorrected_heterogeneity_berger = heterogeneity(uncorr_cluster_batch_freq, method = "berger") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    uncorrected_heterogeneity_boone = heterogeneity(uncorr_cluster_batch_freq, method = "boone") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    uncorrected_heterogeneity_brillouin = heterogeneity(uncorr_cluster_batch_freq, method = "brillouin") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    uncorrected_heterogeneity_mcintosh = heterogeneity(uncorr_cluster_batch_freq, method = "mcintosh") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    
    heterogeneity_matrix_uncorr = as.data.frame(cbind(uncorrected_heterogeneity_shannon, uncorrected_heterogeneity_berger,
                                                      uncorrected_heterogeneity_boone, uncorrected_heterogeneity_brillouin,
                                                      uncorrected_heterogeneity_mcintosh))
    heterogeneity_matrix_uncorr$Cluster <- paste0(method,"_Cluster_",1:nrow(heterogeneity_matrix_uncorr))
    heterogeneity_matrix_uncorr <- heterogeneity_matrix_uncorr %>%
      dplyr::relocate(Cluster)
    colnames(heterogeneity_matrix_uncorr) <- gsub("uncorrected",corrMethod1,colnames(heterogeneity_matrix_uncorr))
    heterogeneity_matrix_uncorr
    
  })
  uncorr_avg_het_react <- reactive({
    
    #req(!isTruthy(DiversityEval_batch()))
    req(DiversityEval_batch())
    heterogeneity_matrix_uncorr <- uncorr_het_react()[,-1]
    avg_heterogeneity_vec <- apply(heterogeneity_matrix_uncorr,2,function(x) mean(x,na.rm = T))
    avg_heterogeneity_df <- data.frame(Heterogeneity_Method = names(avg_heterogeneity_vec),
                                       Average_Heterogeneity = unname(avg_heterogeneity_vec))
    avg_heterogeneity_df
    
  })
  output$uncorr_avg_het <- DT::renderDataTable({
    
    df <- uncorr_avg_het_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = 2, digits = 4)
    
  })
  output$uncorr_het <- DT::renderDataTable({
    
    df <- uncorr_het_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(#lengthMenu = c(5,10,20,100,1000),
                    #pageLength = 5,
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
  })
  
  uncorr_evn_react <- reactive({
    
    #req(!isTruthy(DiversityEval_batch()))
    req(DiversityEval_batch())
    meta <- ClusterInfoTab_react()
    method <- input$ClusterMethod
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    uncorr_col <- paste0(method,"_Cluster_",corrMethod1)
    #uncorr_col <- paste0(method,"_Cluster_Uncorrected")
    batch <- DiversityEval_batch()
    
    uncorr_cluster_batch_freq <- as.data.frame.matrix(table(meta[,uncorr_col], meta[,batch]))
    
    uncorrected_evenness_shannon = evenness(uncorr_cluster_batch_freq, method="shannon")
    uncorrected_evenness_brillouin = evenness(uncorr_cluster_batch_freq, method="brillouin")
    uncorrected_evenness_mcintosh = evenness(uncorr_cluster_batch_freq, method="mcintosh")
    uncorrected_evenness_simpson = evenness(uncorr_cluster_batch_freq, method="simpson")
    
    evenness_matrix_uncorr = as.data.frame(cbind(uncorrected_evenness_shannon, uncorrected_evenness_brillouin,
                                                 uncorrected_evenness_mcintosh, uncorrected_evenness_simpson))
    evenness_matrix_uncorr$Cluster <- paste0(method,"_Cluster_",1:nrow(evenness_matrix_uncorr))
    evenness_matrix_uncorr <- evenness_matrix_uncorr %>%
      dplyr::relocate(Cluster)
    colnames(evenness_matrix_uncorr) <- gsub("uncorrected",corrMethod1,colnames(evenness_matrix_uncorr))
    evenness_matrix_uncorr
    
  })
  
  uncorr_avg_evn_react <- reactive({
    
    #req(!isTruthy(DiversityEval_batch()))
    req(DiversityEval_batch())
    evenness_matrix_uncorr <- uncorr_evn_react()[,-1]
    avg_evenness_vec <- apply(evenness_matrix_uncorr,2,function(x) mean(x,na.rm = T))
    avg_evenness_df <- data.frame(Evenness_Method = names(avg_evenness_vec),
                                  Average_Evenness = unname(avg_evenness_vec))
    avg_evenness_df
    
  })
  output$uncorr_avg_evn <- DT::renderDataTable({
    
    df <- uncorr_avg_evn_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = 2, digits = 4)
    
  })
  output$uncorr_evn <- DT::renderDataTable({
    
    df <- uncorr_evn_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
  })
  
  
  ## Matrix 2 -----------------------------------------------------------------
  
  corr_bar_plot_react <- reactive({
    
    req(input$BarPFillCol)
    plot_df <- ClusterInfoTab_react()
    batch <- input$BarPFillCol
    #batch <- DiversityEval_batch()
    method <- input$ClusterMethod
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    corrMethod2 <- gsub(" ","_",matrix2Title_react())
    corrMethod2 <- ifelse(corrMethod1==corrMethod2,paste0(corrMethod2,"_2"),corrMethod2)
    cluster_col <- paste0(method,"_Cluster_",corrMethod2)
    #corrMethod <- gsub(" ","_",input$batch_correction_method)
    #cluster_col <- paste0(method,"_Cluster_",corrMethod,"_Corrected")
    FillChoice <- input$barPfill
    plot_df[,batch] <- as.factor(plot_df[,batch])
    tickFont <- input$barPAxisTkSize
    titleFont <- input$barPAxisTtSize
    barp <- ggplot(plot_df, aes(fill = !!sym(batch), x = !!sym(cluster_col)))
    if (FillChoice) {
      barp <- barp + geom_bar(position = "fill") +
        theme_minimal()
    } else {
      barp <- barp + geom_bar() +
        theme_minimal()
    }
    barp <- barp +
      xlab(paste(method,"Cluster"))
    barp <- barp + theme(axis.text.x = element_text(size = tickFont),
                         axis.title.x = element_text(size = titleFont),
                         axis.text.y = element_text(size = tickFont),
                         axis.title.y = element_text(size = titleFont))
    barp
    
  })
  
  output$corr_bar_plot <- renderPlot({
    
    barp <- corr_bar_plot_react()
    barp
    
  })
  corr_het_react <- reactive({
    
    ##req(!isTruthy(DiversityEval_batch()))
    req(DiversityEval_batch())
    meta <- ClusterInfoTab_react()
    method <- input$ClusterMethod
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    corrMethod2 <- gsub(" ","_",matrix2Title_react())
    corrMethod2 <- ifelse(corrMethod1==corrMethod2,paste0(corrMethod2,"_2"),corrMethod2)
    corr_col <- paste0(method,"_Cluster_",corrMethod2)
    #corrMethod <- gsub(" ","_",input$batch_correction_method)
    #corr_col <- paste0(method,"_Cluster_",corrMethod,"_Corrected")
    batch <- DiversityEval_batch()
    
    corr_cluster_batch_freq <- as.data.frame.matrix(table(meta[,corr_col], meta[,batch]))
    
    corrected_heterogeneity_shannon = heterogeneity(corr_cluster_batch_freq, method = "shannon") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    corrected_heterogeneity_berger = heterogeneity(corr_cluster_batch_freq, method = "berger") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    corrected_heterogeneity_boone = heterogeneity(corr_cluster_batch_freq, method = "boone") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    corrected_heterogeneity_brillouin = heterogeneity(corr_cluster_batch_freq, method = "brillouin") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    corrected_heterogeneity_mcintosh = heterogeneity(corr_cluster_batch_freq, method = "mcintosh") # other options include c("berger", "boone", "brillouin", "mcintosh", "shannon", "simpson")
    
    heterogeneity_matrix_corr = as.data.frame(cbind(corrected_heterogeneity_shannon, corrected_heterogeneity_berger,
                                                    corrected_heterogeneity_boone, corrected_heterogeneity_brillouin,
                                                    corrected_heterogeneity_mcintosh))
    heterogeneity_matrix_corr$Cluster <- paste0(method,"_Cluster_",1:nrow(heterogeneity_matrix_corr))
    heterogeneity_matrix_corr <- heterogeneity_matrix_corr %>%
      dplyr::relocate(Cluster)
    colnames(heterogeneity_matrix_corr) <- gsub("corrected",corrMethod2,colnames(heterogeneity_matrix_corr))
    heterogeneity_matrix_corr
    
  })
  corr_avg_het_react <- reactive({
    
    ##req(!isTruthy(DiversityEval_batch()))
    req(DiversityEval_batch())
    heterogeneity_matrix_corr <- corr_het_react()[,-1]
    avg_heterogeneity_vec <- apply(heterogeneity_matrix_corr,2,function(x) mean(x,na.rm = T))
    avg_heterogeneity_df <- data.frame(Heterogeneity_Method = names(avg_heterogeneity_vec),
                                       Average_Heterogeneity = unname(avg_heterogeneity_vec))
    avg_heterogeneity_df
    
  })
  output$corr_avg_het <- DT::renderDataTable({
    
    df <- corr_avg_het_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = 2, digits = 4)
    
  })
  output$corr_het <- DT::renderDataTable({
    
    df <- corr_het_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
  })
  
  corr_evn_react <- reactive({
    
    #req(!isTruthy(DiversityEval_batch()))
    req(DiversityEval_batch())
    meta <- ClusterInfoTab_react()
    method <- input$ClusterMethod
    corrMethod1 <- gsub(" ","_",matrix1Title_react())
    corrMethod2 <- gsub(" ","_",matrix2Title_react())
    corrMethod2 <- ifelse(corrMethod1==corrMethod2,paste0(corrMethod2,"_2"),corrMethod2)
    corr_col <- paste0(method,"_Cluster_",corrMethod2)
    #corrMethod <- gsub(" ","_",input$batch_correction_method)
    #corr_col <- paste0(method,"_Cluster_",corrMethod,"_Corrected")
    batch <- DiversityEval_batch()
    
    corr_cluster_batch_freq <- as.data.frame.matrix(table(meta[,corr_col], meta[,batch]))
    
    corrected_evenness_shannon = evenness(corr_cluster_batch_freq, method="shannon")
    corrected_evenness_brillouin = evenness(corr_cluster_batch_freq, method="brillouin")
    corrected_evenness_mcintosh = evenness(corr_cluster_batch_freq, method="mcintosh")
    corrected_evenness_simpson = evenness(corr_cluster_batch_freq, method="simpson")
    
    evenness_matrix_corr = as.data.frame(cbind(corrected_evenness_shannon, corrected_evenness_brillouin,
                                               corrected_evenness_mcintosh, corrected_evenness_simpson))
    evenness_matrix_corr$Cluster <- paste0(method,"_Cluster_",1:nrow(evenness_matrix_corr))
    evenness_matrix_corr <- evenness_matrix_corr %>%
      dplyr::relocate(Cluster)
    colnames(evenness_matrix_corr) <- gsub("corrected",corrMethod2,colnames(evenness_matrix_corr))
    evenness_matrix_corr
    
  })
  
  corr_avg_evn_react <- reactive({
    
    #req(!isTruthy(DiversityEval_batch()))
    req(DiversityEval_batch())
    evenness_matrix_corr <- corr_evn_react()[,-1]
    avg_evenness_vec <- apply(evenness_matrix_corr,2,function(x) mean(x,na.rm = T))
    avg_evenness_df <- data.frame(Evenness_Method = names(avg_evenness_vec),
                                  Average_Evenness = unname(avg_evenness_vec))
    avg_evenness_df
    
  })
  output$corr_avg_evn <- DT::renderDataTable({
    
    df <- corr_avg_evn_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = 2, digits = 4)
    
  })
  output$corr_evn <- DT::renderDataTable({
    
    df <- corr_evn_react()
    colnames(df) <- gsub("_"," ",colnames(df))
    df[,1] <- gsub("_"," ",df[,1])
    DT::datatable(df,
                  options = list(
                    dom = '',
                    scrollX = T),
                  rownames = F) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
  })
  
  # RLE -----------------------------------------------------------------------
  output$RLE1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"RLE Plot")
  })
  output$RLE2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"RLE Plot")
  })
  
  # Generating uncorrected RLE plot
  output$batch_choices_RLE <- shiny::renderUI({
    shiny::selectInput("batch_choices_RLE",
                       "Color by:",
                       #c("Select", unlist(strsplit(batch_names_from_meta(), ","))),
                       batch_names_from_meta()[-1],
                       selected = batch1_choices())
  })
  uncorrected_batch_choices_RLE2 <- shiny::reactive({
    if (isTruthy(input$batch_choices_RLE)) {
      if (input$batch_choices_RLE == "Select"){
        uncorrected_batch_choices_RLE2 <- NULL
      }else {
        uncorrected_batch_choices_RLE2 <- input$batch_choices_RLE
      }
      uncorrected_batch_choices_RLE2
    }
  })
  
  ## Matrix 1 -----------------------------------------------------------------
  
  RLE_Obj_Uncorr <- reactive({
    
    req(uncorrected_batch_choices_RLE2())
    mat_RLE <- matrix1()
    #rownames(mat_RLE) <- uncorrected_numeric_matrix()[,1]
    meta <- aligned_meta_file()
    batch.1 <- uncorrected_batch_choices_RLE2()
    rle_mat <- mat_RLE
    rle_mat_meta <- cbind(t(rle_mat), meta)
    rle_mat_meta_ordered <- rle_mat_meta[order(rle_mat_meta[[batch.1]]),]
    rle_mat_ordered <- t(rle_mat_meta_ordered[,!colnames(rle_mat_meta_ordered) %in% colnames(meta)])
    rle_meta <- rle_mat_meta_ordered[,colnames(rle_mat_meta_ordered) %in% colnames(meta)]
    
    uncorrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(rle_mat_ordered)),
      colData = rle_meta,
      rowData = rownames(mat_RLE)
    )
    uncorrected_RLE_SCE
    
  })
  uncorrected_RLE <- shiny::reactive({
    rle_obj <- RLE_Obj_Uncorr()
    req(rle_obj)
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
  
  
  ## Matrix 2 -----------------------------------------------------------------
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
    
    req(uncorrected_batch_choices_RLE2())
    mat_RLE <- matrix2()
    #mat_RLE <- batch_correction()
    meta <- aligned_meta_file()
    batch.1 <- uncorrected_batch_choices_RLE2()
    rle_mat <- mat_RLE
    rle_mat_meta <- cbind(t(rle_mat), meta)
    rle_mat_meta_ordered <- rle_mat_meta[order(rle_mat_meta[[batch.1]]),]
    rle_mat_ordered <- t(rle_mat_meta_ordered[,!colnames(rle_mat_meta_ordered) %in% colnames(meta)])
    rle_meta <- rle_mat_meta_ordered[,colnames(rle_mat_meta_ordered) %in% colnames(meta)]
    
    corrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(rle_mat_ordered)),
      colData = rle_meta,
      rowData = rownames(mat_RLE)
    )
    corrected_RLE_SCE
    
  })
  corrected_RLE <- shiny::reactive({
    rle_obj <- RLE_Obj_Corr()
    req(rle_obj)
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
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
  
  # Exp Var -------------------------------------------------------------------
  output$ExpVar1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"Explanatory Variables Plot")
  })
  output$ExpVar2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"Explanatory Variables Plot")
  })
  
  # Generating EV plot
  output$variable_choices_EV <- shiny::renderUI({
    shiny::selectInput("variable_choices_EV",
                       "Select Variables to Plot",
                       c(batch_names_from_meta()[-1]),
                       selected = head(c(batch_names_from_meta()[-1]), n = 3),
                       #c(unlist(strsplit(batch_names_from_meta(), ","))),
                       multiple = TRUE
    )
  })
  output$rendExpPlotLog <- renderUI({
    
    if (input$ExpVar_Plots == "exp_plot") {
      checkboxInput("ExpPlotLog","Log Scale X-Axis", value = T)
    }
    
  })
  
  ## Matrix 1 -----------------------------------------------------------------
  
  uncorrected_EV_df <- shiny::reactive({
    req(input$variable_choices_EV)
    my_colors <- metafolio::gg_color_hue(length(input$variable_choices_EV))
    mat_EV <- matrix1()
    names(mat_EV) <- NULL
    uncorrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_EV)),
      colData = aligned_meta_file(),
      rowData = rownames(mat_EV)
    )
    SummarizedExperiment::assay(uncorrected_EV_SCE, "logcounts") <- SingleCellExperiment::counts(uncorrected_EV_SCE)
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
  
  ## Matrix 2 -----------------------------------------------------------------
  
  corrected_EV_df <- shiny::reactive({
    req(input$variable_choices_EV)
    mat_EV <- matrix2()
    names(mat_EV) <- NULL
    my_colors <- metafolio::gg_color_hue(length(input$variable_choices_EV))
    corrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(mat_EV)),
      colData = aligned_meta_file(),
      rowData = rownames(mat_EV)
    )
    SummarizedExperiment::assay(corrected_EV_SCE, "logcounts") <- SingleCellExperiment::counts(corrected_EV_SCE)
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
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
  
  
  # PVCA ----------------------------------------------------------------------
  output$PVCA1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"PVCA Plot")
  })
  output$PVCA2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"PVCA Plot")
  })
  output$rendpvcaPct <- renderUI({
    
    if (input$ExpVar_Plots == "pvca") {
      numericInput("pvcaPct","% Threshold", value = 0.8, min = 0, max = 1)
    }
    
  })
  
  output$rendErrorMessagePVCAUncorr <- renderUI({
    if (length(input$variable_choices_EV) < 2) {
      span(textOutput("ErrorMessagePVCAUncorr_PVCA"), style="color:red")
    }
  })
  output$rendErrorMessagePVCACorr <- renderUI({
    if (length(input$variable_choices_EV) < 2) {
      span(textOutput("ErrorMessagePVCACorr_PVCA"), style="color:red")
    }
  })
  output$ErrorMessagePVCAUncorr <- renderText({"Select two or more vairables to continue"})
  output$ErrorMessagePVCACorr <- renderText({"Select two or more vairables to continue"})
  
  ## Matrix 1 -----------------------------------------------------------------
  
  pvca_mv_features_uncorr_matrix <- reactive({
    
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$pvcacluster_n_MV_features
      var_type <- input$pvcaVarianceMeasure
      mat <- matrix1()
      featColName <- colnames(mat)[1]
      #rownames(mat) <- mat[,1]
      #mat <- mat[,-1]
      
      mad <- NULL
      var <- NULL
      cv <- NULL
      if (!input$Log_Choice) {
        mat <- log2(mat + 1)
      }
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(mat, 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(mat, 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(mat, 1, cv))
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
      if (length(vars) >= 2) {
        withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
  
  
  ## Matrix 2 -----------------------------------------------------------------
  
  pvca_mv_features_corr_matrix <- reactive({
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      incProgress(0.5, detail = "Calculating Most Variable Features")
      topN <- input$pvcacluster_n_MV_features
      var_type <- input$pvcaVarianceMeasure
      mat <- matrix2()
      #rownames(mat) <- mat[,1]
      #mat <- mat[,-1]
      
      mad <- NULL
      var <- NULL
      cv <- NULL
      if (!input$Log_Choice) {
        mat <- log2(mat + 1)
      }
      if (var_type == "MAD"){
        mad <- suppressWarnings(apply(mat, 1, mad))
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = topN)
        out <- cbind(names(mad), mad[names(mad)], mat[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(mat))
        dataset <- mat[names(mad),]
      }
      else if (var_type == "VAR"){
        var <- suppressWarnings(apply(mat, 1, var))
        var <- sort(var, decreasing = T)
        var <- head(var, n = topN)
        out <- cbind(names(var), var[names(var)], mat[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(mat))
        dataset <- mat[names(var),]
      }
      else if (var_type == "CV"){
        cv <- suppressWarnings(apply(mat, 1, cv))
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
      if (length(vars) >= 2) {
        withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
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
  
  # SVA -----------------------------------------------------------------------
  output$SVA1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"SVA Analysis")
  })
  output$SVA2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"SVA Analysis")
  })
  output$uncorrected_SVA_variable_of_interest <- renderUI({
    req(aligned_meta_file())
    VarChoices <- c(batch_names_from_meta()[-1])
    meta <- aligned_meta_file()
    if (isTruthy(SVA_variable_of_interest_bc2())) {
      VarSelected <- SVA_variable_of_interest_bc2()
    } else {
      VarSelected <- VarChoices[which(!VarChoices %in% c(input$batch1_choices,input$batch2_choices))]
    }
    WorkingVars <- sapply(colnames(meta), function(x) {
      (length(table(meta[,x])) > 1) & (length(table(meta[,x])) < length(meta[,x]) | all(is.numeric(meta[,x])))
    })
    FullVars <- sapply(colnames(meta), function(x) {
      (length(meta[which(is.na(meta[,x])),x]) == 0)
    })
    FullVars <- names(FullVars[which(FullVars == TRUE)])
    WorkingVars <- names(WorkingVars[which(WorkingVars == TRUE)])
    WorkingVars <- intersect(WorkingVars,FullVars)
    VarChoices <- intersect(WorkingVars,VarChoices)
    selectInput(
      "uncorrected_SVA_variable_of_interest",
      "Select Variable of Interest",
      VarChoices,
      selected = VarSelected
    )
  })
  
  observe({
    
    VarChoices <- c(batch_names_from_meta()[-1])
    if (isTruthy(input$uncorrected_SVA_variable_of_interest)) {
      VarSelected <- input$uncorrected_SVA_variable_of_interest
      updateSelectizeInput(session,"SVAcolor", choices = VarChoices, selected = VarSelected, server = T)
    } 
    
  })
  
  output$rendSVAhover <- shiny::renderUI({
    batch_choice1 <- batch1_choices1()
    batch_choice2 <- batch1_choices2()
    batch_choice <- unique(c(batch_choice1,batch_choice2))
    VOI <- uncorrected_SVA_variable_of_interest2()
    shiny::selectInput("SVAhover",
                       "Information to Display on Hover:",
                       choices = unlist(strsplit(batch_names_from_meta(), ",")),
                       selected = c(batch_names_from_meta()[1],VOI,batch_choice),
                       multiple = T)
  })
  
  observeEvent(input$CalcSurVar, {
    output$rendSurVarEstimates <- renderUI({
      
      shiny::fluidRow(
        p(),
        tags$style(type='text/css', '#uncorrected_SVA_nsv_print {white-space: pre-wrap;}'),
        tags$style(type='text/css', '#corrected_SVA_nsv_print {white-space: pre-wrap;}'),
        column(12,
               h4(HTML("<b>Number of Estimated Surrogate Variables</b>"), 
                  style="text-align:center"),
               fluidRow(
                 column(6,
                        h5(HTML(paste("Matrix 1:",uncorrected_SVA_nsv())), 
                           style="text-align:center")
                        ),
                 column(6,
                        h5(HTML(paste("Matrix 2:",corrected_SVA_nsv())), 
                           style="text-align:center")
                        )
               ),
               hr(),
               shiny::h4("Download Surrogate Variables"),
               fluidRow(
                 column(6,
                        actionButton("save_SVA_surrogate_variables", "Add to Zip File Export")
                 ),
                 column(6,
                        shiny::downloadButton("dnldsave_SVA_surrogate_variables","Dowload Meta with SVA")
                 )
               ))
      )
      
    })
  })
  
  ## Matrix 1 -----------------------------------------------------------------
  
  SVA1_pass <- reactiveVal(1)
  SV1_error <- reactiveVal()
  SVA2_pass <- reactiveVal(1)
  SV2_error <- reactiveVal()
  uncorrected_SVA_variable_of_interest2 <- reactive({
    if (isTruthy(input$uncorrected_SVA_variable_of_interest)) {
      if (input$uncorrected_SVA_variable_of_interest == "Select"){
        NULL
      }else {
        uncorrected_SVA_variable_of_interest2 <- input$uncorrected_SVA_variable_of_interest
      }
    }
  })
  uncorrected_SVA_nsv <- eventReactive(input$CalcSurVar, {
    req(matrix1())
    if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      Feature <- uncorrected_SVA_variable_of_interest2()
      meta <- aligned_meta_file()
      colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
      withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
        incProgress(0.5, detail = "Number of surrogate variables")
        uncorrected_mod <- model.matrix(reformulate("Feature"), data = meta)
        uncorrected_mod_null <- model.matrix(~1, data = meta)
        uncorrected_SVA_nsv <- sva::num.sv(matrix1(), uncorrected_mod, method=svaMethod, vfilter = svaVarNum)
        incProgress(0.5, detail = "Complete!")
      })
      uncorrected_SVA_nsv
    }
  })
  #uncorrected_SVA_nsv <- reactive({
  #  req(matrix1())
  #  if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
  #    svaMethod <- input$svaMethod
  #    svaVarNum <- input$SVAvarNum
  #    Feature <- uncorrected_SVA_variable_of_interest2()
  #    meta <- aligned_meta_file()
  #    colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
  #    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
  #      incProgress(0.5, detail = "Number of surrogate variables")
  #      uncorrected_mod <- model.matrix(reformulate("Feature"), data = meta)
  #      uncorrected_mod_null <- model.matrix(~1, data = meta)
  #      uncorrected_SVA_nsv <- sva::num.sv(matrix1(), uncorrected_mod, method=svaMethod, vfilter = svaVarNum)
  #      incProgress(0.5, detail = "Complete!")
  #    })
  #    uncorrected_SVA_nsv
  #  }
  #})
  #observe({
  #  n.sv <- uncorrected_SVA_nsv()
  #  updateNumericInput(session,"NSV1", value = n.sv)
  #})
  observe({
    n.sv <- input$NSV1
    updateSelectInput(session,"SVAxAxis1", choices = sapply(seq(n.sv),function(x) paste0("SVA_",x)), selected = "SVA_1")
    updateSelectInput(session,"SVAyAxis1", choices = sapply(seq(n.sv),function(x) paste0("SVA_",x)), selected = "SVA_2")
  })
  observe({
    if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      #req(uncorrected_SVA_nsv())
      req(input$NSV1)
      req(matrix1())
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      Feature <- uncorrected_SVA_variable_of_interest2()
      meta <- aligned_meta_file()
      colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
      uncorrected_mod <- model.matrix(reformulate("Feature"), data = meta)
      uncorrected_mod_null <- model.matrix(~1, data = meta)
      uncorrected_SVA_matrix <- as.matrix(matrix1())
      res <- try(sva::sva(uncorrected_SVA_matrix, uncorrected_mod, uncorrected_mod_null, n.sv = input$NSV1,
                          numSVmethod = svaMethod, vfilter = svaVarNum))
      if (inherits(res, "try-error")) {
        SV1_error(sub("\n","",as.character(res)))
        SVA1_pass(2)
      } else {
        SVA1_pass(3)
      }
    }
  })
  output$rendSVAtextError1 <- renderUI({
    req(SVA1_pass())
    #req(SV1_error())
    if (SVA1_pass() == 2) {
      span(textOutput(paste("SVAtextError1")), style = "color:red")
      #span(textOutput(paste("SVAtextError1",SV1_error(), sep = '\n')), style = "color:red")
    }
  })
  output$SVAtextError1 <- renderText({
    "Parameters require adjustment"
    #paste0(c("Parameters require adjustment",SV1_error()),collapse = '\n')
  })
  uncorrected_SVA_object <- reactive({
    if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      #req(uncorrected_SVA_nsv())
      req(input$NSV1)
      req(matrix1())
      if (SVA1_pass() == 3) {
        svaMethod <- input$svaMethod
        svaVarNum <- input$SVAvarNum
        Feature <- uncorrected_SVA_variable_of_interest2()
        meta <- aligned_meta_file()
        colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
        withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
          incProgress(0.5, detail = "Running SVA")
          uncorrected_mod <- model.matrix(reformulate("Feature"), data = meta)
          uncorrected_mod_null <- model.matrix(~1, data = meta)
          uncorrected_SVA_matrix <- as.matrix(matrix1())
          uncorrected_SVA_object <- sva::sva(uncorrected_SVA_matrix, uncorrected_mod, uncorrected_mod_null, n.sv = input$NSV1,
                                             numSVmethod = svaMethod, vfilter = svaVarNum)
          save(list = ls(), file = "SVA_Env.RData",envir = environment())
          incProgress(0.5, detail = "Complete!")
        })
        #df <- as.data.frame(uncorrected_SVA_object$sv)
        #if (ncol(df) > 0) {
        #  colnames(df) <- paste0("SVA_",matrix1DlndTitle_react(),"_Surrogate_Vars_",seq(ncol(df)))
        #  df <- cbind(aligned_meta_file(),df)
        #  uncorr_boxplot_meta_file(df)
          uncorrected_SVA_object
        #}
      }
    }
  })
  observe({
    df <- as.data.frame(uncorrected_SVA_object()$sv)
    if (ncol(df) > 0) {
      colnames(df) <- paste0("SVA_",matrix1DlndTitle_react(),"_SV_",seq(ncol(df)))
      df <- cbind(aligned_meta_file(),df)
      uncorr_boxplot_meta_file(df)
    }
  })
  toListen <- reactive({
    list(uncorr_boxplot_meta_file(),corr_boxplot_meta_file())
  })
  observe({
    #main_meta <- aligned_meta_file()
    req(boxplot_meta_file())
    main_meta <- boxplot_meta_file()
    uncorr_meta <- uncorr_boxplot_meta_file()
    corr_meta <- corr_boxplot_meta_file()
    #main_meta_new <- merge(main_meta,uncorr_meta[,which(!colnames(main_meta)[-1] %in% colnames(uncorr_meta))], all.x = T)
    #main_meta_new <- merge(main_meta_new,corr_meta[,which(!colnames(main_meta)[-1] %in% colnames(corr_meta))], all.x = T)
    #main_meta_new <- merge(main_meta,uncorr_meta, by.x = colnames(main_meta)[1], by.y = colnames(uncorr_meta)[1], all.x = T)
    #main_meta_new <- merge(main_meta_new,corr_meta, by.x = colnames(main_meta)[1], by.y = colnames(corr_meta)[1], all.x = T)
    #if (any(!colnames(uncorr_meta) %in% colnames(main_meta))) {
    #  main_meta_new <- merge(main_meta,uncorr_meta, by.x = colnames(main_meta)[1], by.y = colnames(uncorr_meta)[1], all.x = T)
    #}
    #if (any(!colnames(corr_meta) %in% colnames(main_meta))) {
    #  main_meta_new <- merge(main_meta_new,corr_meta, by.x = colnames(main_meta)[1], by.y = colnames(corr_meta)[1], all.x = T)
    #}
    #boxplot_meta_file(main_meta_new)
    
    #print("main_meta")
    #print(head(main_meta))
    #print("uncorr_meta")
    #print(head(uncorr_meta))
    #print("corr_meta")
    #print(head(corr_meta))
    main_meta_new <- merge(main_meta,uncorr_meta, all.x = T)
    #print("main_meta_new1")
    #print(head(main_meta_new))
    main_meta_new <- merge(main_meta_new,corr_meta, all.x = T)
    #print("main_meta_new2")
    #print(head(main_meta_new))
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
  
  #uncorrected_SVA_probability_ggplot <- reactive({
  #  if (isTruthy(uncorrected_SVA_probability_df())) {
  #    ggplot(uncorrected_SVA_probability_df(), aes(x = probability_association_of_each_gene,
  #                                                 fill = variable_type)) +
  #      geom_density(alpha = 0.5)
  #  }
  #})
  #
  #output$uncorrected_SVA_probability_density <- renderPlot({
  #  req(uncorrected_SVA_probability_ggplot())
  #  uncorrected_SVA_probability_ggplot()
  #})
  output$uncorrected_SVA_nsv_print <- renderPrint({
    if (isTruthy(uncorrected_SVA_nsv())) {
      print(paste("Number of Estimated Surrogate Variables in Matrix 1:", uncorrected_SVA_nsv()))
      #print(paste("The Number of Estimated Surrogate Variables are:", input$NSV1))
    }
  })
  
  uncorrected_SVA_scatter_df <- reactive({
    
    req(uncorrected_SVA_object())
    df <- uncorrected_SVA_object()$sv
    if (input$SVAxAxis1 != input$SVAyAxis1) {
      colnames(df) <- sapply(seq(ncol(df)),function(x) paste0("SVA_",x))
      rownames(df) <- colnames(matrix1())
      meta <- aligned_meta_file()
      batch_choice <- uncorrected_SVA_variable_of_interest2()
      #pca <- uncorrected_PCA_react()
      #batch_choice <- input$batch_choices_PCA
      hover_choice <- input$SVAhover
      svaColor <- input$SVAcolor
      metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice,svaColor))
      metaCols <- metaCols[which(metaCols %in% colnames(meta))]
      SVA_x <- input$SVAxAxis1
      SVA_y <- input$SVAyAxis1
      uncorrected_SVA_points <- as.data.frame(df[,c(SVA_x,SVA_y)])
      uncorrected_SVA_plot_df <- merge(df,meta[,metaCols,drop = F],by.x = 0, by.y = colnames(meta)[1])
      colnames(uncorrected_SVA_plot_df)[1] <- colnames(uncorrected_SVA_plot_df)[1]
      #PCX <- gsub("^PC","",PC_x)
      #PCY <- gsub("^PC","",PC_y)
      #uncorrected_PCA_plot_df[,PC_x] <- as.numeric(uncorrected_PCA_plot_df[,PC_x]) / (pca$sdev[as.numeric(PCX)] * sqrt(nrow(uncorrected_PCA_plot_df)))
      #uncorrected_PCA_plot_df[,PC_y] <- as.numeric(uncorrected_PCA_plot_df[,PC_y]) / (pca$sdev[as.numeric(PCY)] * sqrt(nrow(uncorrected_PCA_plot_df)))
      uncorrected_SVA_plot_df$text <- NA
      for (r in seq(nrow(uncorrected_SVA_plot_df))) {
        text_vec <- c()
        for (c in metaCols) {
          text_vec <- c(text_vec,
                        paste0("\\<br\\>\\<b\\>", c, ":\\</b\\> ",uncorrected_SVA_plot_df[r,c]))
        }
        text = paste(text_vec,collapse = '')
        text = gsub("\\\\", '', text)
        text = paste0(text,"<extra></extra>") #removes outside of box text
        uncorrected_SVA_plot_df$text[r] = text
      }
      uncorrected_SVA_plot_df
    }
    
    
  })
  
  uncorrected_SVA_scatter_react <- reactive({
    
    req(uncorrected_SVA_scatter_df)
    plot_df <- uncorrected_SVA_scatter_df()
    
    meta <- aligned_meta_file()
    #plot_df <- uncorrected_PCA_plot_df()
    batch_choice <- uncorrected_SVA_variable_of_interest2()
    hover_choice <- input$SVAhover
    svaColor <- input$SVAcolor
    metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice))
    SVA_x <- input$SVAxAxis1
    SVA_y <- input$SVAyAxis1
    colnames(plot_df)[which(colnames(plot_df) == SVA_x)] <- "X"
    colnames(plot_df)[which(colnames(plot_df) == SVA_y)] <- "Y"
    dotSize <- input$SVAdotSize
    fontSize <- input$SVAfontSize
    
    if (isTruthy(batch_choice) & isTruthy(svaColor)) {
      #print(svaColor)
      #if (batch_choice != "Select") {
        colnames(plot_df)[which(colnames(plot_df) == svaColor)] <- "ColorBatch"
        p4 <- plot_ly() %>%
          add_trace(
            data = plot_df,
            x = ~X,
            y = ~Y,
            type = "scatter",
            mode = "markers",
            color = ~ColorBatch,
            legendgroup=svaColor,
            marker=list(size=dotSize,symbol = 'circle'),
            hovertemplate = ~text
          ) %>%
          layout(xaxis = list(zeroline = FALSE,title = SVA_x),
                 yaxis = list(zeroline = FALSE,title = SVA_y),
                 font = list(size = fontSize)) %>%
          config(
            toImageButtonOptions = list(
              format = "svg"
            )
          )
        p4
      #} else {
      #  p4 <- plot_ly() %>%
      #    add_trace(
      #      data = plot_df,
      #      x = ~X,
      #      y = ~Y,
      #      type = "scatter",
      #      mode = "markers",
      #      marker = list(color = "black"),
      #      marker=list(size=dotSize,symbol = 'circle'),
      #      hovertemplate = ~text
      #    ) %>%
      #    layout(xaxis = list(zeroline = FALSE,title = SVA_x),
      #           yaxis = list(zeroline = FALSE,title = SVA_y),
      #           font = list(size = fontSize)) %>%
      #    config(
      #      toImageButtonOptions = list(
      #        format = "svg"
      #      )
      #    )
      #  p4
      #}
    }
    
  })
  
  output$uncorrected_SVA_scatter <- renderPlotly({
    
    req(uncorrected_SVA_scatter_react())
    p <- uncorrected_SVA_scatter_react()
    p
    
  })
  output$uncorrected_SVA_scatter_df <- DT::renderDataTable({
    
    req(uncorrected_SVA_scatter_df())
    df <- uncorrected_SVA_scatter_df()
    df <- df[,-ncol(df)]
    colnames(df)[1] <- "SampleName"
    DT::datatable(df,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F) %>%
      formatRound(columns = grep("^SVA_",colnames(df),value = T), digits = 4)
    
  })
  
  ## Matrix 2 -----------------------------------------------------------------
  
  corrected_SVA_variable_of_interest2 <- reactive({
    if (isTruthy(input$uncorrected_SVA_variable_of_interest)) {
      if (input$uncorrected_SVA_variable_of_interest == "Select"){
        NULL
      }else {
        corrected_SVA_variable_of_interest2 <- input$uncorrected_SVA_variable_of_interest
      }
    }
  })
  corrected_SVA_nsv <- eventReactive(input$CalcSurVar, {
    if (isTruthy(corrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      req(matrix2())
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      Feature <- corrected_SVA_variable_of_interest2()
      meta <- aligned_meta_file()
      colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
      #withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
      #incProgress(0.5, detail = "Number of surrogate variables")
      corrected_mod <- model.matrix(reformulate("Feature"), data = meta)
      corrected_mod_null <- model.matrix(~1, data = meta)
      corrected_SVA_nsv <- sva::num.sv(matrix2(), corrected_mod, method = svaMethod, vfilter = svaVarNum)
      #incProgress(0.5, detail = "Complete!")
      #})
      corrected_SVA_nsv
    }
  })
  #corrected_SVA_nsv <- reactive({
  #  if (isTruthy(corrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
  #    req(matrix2())
  #    svaMethod <- input$svaMethod
  #    svaVarNum <- input$SVAvarNum
  #    Feature <- corrected_SVA_variable_of_interest2()
  #    meta <- aligned_meta_file()
  #    colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
  #    #withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
  #      #incProgress(0.5, detail = "Number of surrogate variables")
  #      corrected_mod <- model.matrix(reformulate("Feature"), data = meta)
  #      corrected_mod_null <- model.matrix(~1, data = meta)
  #      corrected_SVA_nsv <- sva::num.sv(matrix2(), corrected_mod, method = svaMethod, vfilter = svaVarNum)
  #      #incProgress(0.5, detail = "Complete!")
  #    #})
  #    corrected_SVA_nsv
  #  }
  #})
  #observe({
  #  n.sv <- corrected_SVA_nsv()
  #  updateNumericInput(session,"NSV2", value = n.sv)
  #})
  observe({
    n.sv <- input$NSV2
    updateSelectInput(session,"SVAxAxis2", choices = sapply(seq(n.sv),function(x) paste0("SVA_",x)), selected = "SVA_1")
    updateSelectInput(session,"SVAyAxis2", choices = sapply(seq(n.sv),function(x) paste0("SVA_",x)), selected = "SVA_2")
  })
  observe({
    if (isTruthy(corrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      #req(corrected_SVA_nsv())
      req(input$NSV2)
      req(matrix2())
      svaMethod <- input$svaMethod
      svaVarNum <- input$SVAvarNum
      Feature <- corrected_SVA_variable_of_interest2()
      meta <- aligned_meta_file()
      colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
      corrected_mod <- model.matrix(reformulate("Feature"), data = meta)
      corrected_mod_null <- model.matrix(~1, data = meta)
      corrected_SVA_matrix <- as.matrix(matrix2())
      res <- try(sva::sva(corrected_SVA_matrix, corrected_mod, corrected_mod_null, n.sv = input$NSV2,
                          numSVmethod = svaMethod, vfilter = svaVarNum))
      if (inherits(res, "try-error")) {
        SV2_error(sub("\n","",as.character(res)))
        SVA2_pass(2)
      } else {
        SVA2_pass(3)
      }
    }
  })
  output$rendSVAtextError2 <- renderUI({
    req(SVA2_pass())
    if (SVA2_pass() == 2) {
      span(textOutput("SVAtextError2"), style = "color:red")
    }
  })
  output$SVAtextError2 <- renderText({
    "Parameters require adjustment"
  })
  corrected_SVA_object <- reactive({
    if (isTruthy(uncorrected_SVA_variable_of_interest2()) & isTruthy(aligned_meta_file())) {
      #req(corrected_SVA_nsv())
      req(input$NSV2)
      req(SVA2_pass())
      if (SVA2_pass() == 3) {
        svaMethod <- input$svaMethod
        svaVarNum <- input$SVAvarNum
        Feature <- corrected_SVA_variable_of_interest2()
        meta <- aligned_meta_file()
        colnames(meta)[which(colnames(meta) == Feature)] <- "Feature"
        withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
          incProgress(0.5, detail = "Running SVA")
          corrected_mod <- model.matrix(reformulate("Feature"), data = meta)
          corrected_mod_null <- model.matrix(~1, data = meta)
          corrected_SVA_matrix <- as.matrix(matrix2())
          corrected_SVA_object <- sva::sva(corrected_SVA_matrix, corrected_mod, corrected_mod_null, n.sv = input$NSV2,
                                           numSVmethod = svaMethod, vfilter = svaVarNum)
          incProgress(0.5, detail = "Complete!")
        })
        #df <- as.data.frame(corrected_SVA_object$sv)
        #if (ncol(df) > 0) {
        #  colnames(df) <- paste0("SVA_",matrix2DlndTitle_react(),"_Surrogate_Vars_",seq(ncol(df)))
        #  #colnames(df) <- paste0("SVA_",input$batch_correction_method,"_Corrected_Surrogate_Vars_",seq(ncol(df)))
        #  df <- cbind(aligned_meta_file(),df)
        #  corr_boxplot_meta_file(df)
          corrected_SVA_object
        #}
      }
    }
  })
  observe({
    df <- as.data.frame(corrected_SVA_object()$sv)
    if (ncol(df) > 0) {
      colnames(df) <- paste0("SVA_",matrix2DlndTitle_react(),"_SV_",seq(ncol(df)))
      df <- cbind(aligned_meta_file(),df)
      corr_boxplot_meta_file(df)
    }
  })
  corrected_SVA_probability_df <- reactive({
    if (isTruthy(corrected_SVA_object())) {
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
    if (isTruthy(corrected_SVA_probability_df())) {
      ggplot(corrected_SVA_probability_df(), aes(x = probability_association_of_each_gene,fill = variable_type)) +
        geom_density(alpha = 0.5)
    }
  })
  
  output$corrected_SVA_probability_density <- renderPlot({
    req(corrected_SVA_probability_ggplot())
    corrected_SVA_probability_ggplot()
  })
  
  output$corrected_SVA_nsv_print <- renderPrint({
    if (isTruthy(corrected_SVA_nsv())) {
      print(paste("Estimated Surrogate Variables in Matrix 2:", corrected_SVA_nsv()))
      #print(paste("The Number of Estimated Surrogate Variables are:", input$NSV2))
    }
  })
  
  corrected_SVA_scatter_df <- reactive({
    
    req(corrected_SVA_object())
    df <- corrected_SVA_object()$sv
    if (input$SVAxAxis2 != input$SVAyAxis2) {
      colnames(df) <- sapply(seq(ncol(df)),function(x) paste0("SVA_",x))
      rownames(df) <- colnames(matrix2())
      meta <- aligned_meta_file()
      batch_choice <- uncorrected_SVA_variable_of_interest2()
      #pca <- uncorrected_PCA_react()
      #batch_choice <- input$batch_choices_PCA
      hover_choice <- input$SVAhover
      svaColor <- input$SVAcolor
      metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice,svaColor))
      metaCols <- metaCols[which(metaCols %in% colnames(meta))]
      SVA_x <- input$SVAxAxis2
      SVA_y <- input$SVAyAxis2
      corrected_SVA_points <- as.data.frame(df[,c(SVA_x,SVA_y)])
      corrected_SVA_plot_df <- merge(df,meta[,metaCols,drop = F],by.x = 0, by.y = colnames(meta)[1])
      colnames(corrected_SVA_plot_df)[1] <- colnames(corrected_SVA_plot_df)[1]
      #PCX <- gsub("^PC","",PC_x)
      #PCY <- gsub("^PC","",PC_y)
      #uncorrected_PCA_plot_df[,PC_x] <- as.numeric(uncorrected_PCA_plot_df[,PC_x]) / (pca$sdev[as.numeric(PCX)] * sqrt(nrow(uncorrected_PCA_plot_df)))
      #uncorrected_PCA_plot_df[,PC_y] <- as.numeric(uncorrected_PCA_plot_df[,PC_y]) / (pca$sdev[as.numeric(PCY)] * sqrt(nrow(uncorrected_PCA_plot_df)))
      corrected_SVA_plot_df$text <- NA
      for (r in seq(nrow(corrected_SVA_plot_df))) {
        text_vec <- c()
        for (c in metaCols) {
          text_vec <- c(text_vec,
                        paste0("\\<br\\>\\<b\\>", c, ":\\</b\\> ",corrected_SVA_plot_df[r,c]))
        }
        text = paste(text_vec,collapse = '')
        text = gsub("\\\\", '', text)
        text = paste0(text,"<extra></extra>") #removes outside of box text
        corrected_SVA_plot_df$text[r] = text
      }
      corrected_SVA_plot_df
    }
    
    
  })
  
  corrected_SVA_scatter_react <- reactive({
    
    req(corrected_SVA_scatter_df)
    plot_df <- corrected_SVA_scatter_df()
    
    meta <- aligned_meta_file()
    #plot_df <- uncorrected_PCA_plot_df()
    batch_choice <- uncorrected_SVA_variable_of_interest2()
    hover_choice <- input$SVAhover
    metaCols <- unique(c(colnames(meta)[1],batch_choice,hover_choice))
    SVA_x <- input$SVAxAxis2
    SVA_y <- input$SVAyAxis2
    svaColor <- input$SVAcolor
    colnames(plot_df)[which(colnames(plot_df) == SVA_x)] <- "X"
    colnames(plot_df)[which(colnames(plot_df) == SVA_y)] <- "Y"
    dotSize <- input$SVAdotSize
    fontSize <- input$SVAfontSize
    
    if (isTruthy(batch_choice) & isTruthy(svaColor)) {
      #if (batch_choice != "Select") {
        colnames(plot_df)[which(colnames(plot_df) == svaColor)] <- "ColorBatch"
        p4 <- plot_ly() %>%
          add_trace(
            data = plot_df,
            x = ~X,
            y = ~Y,
            type = "scatter",
            mode = "markers",
            color = ~ColorBatch,
            legendgroup=svaColor,
            marker=list(size=dotSize,symbol = 'circle'),
            hovertemplate = ~text
          ) %>%
          layout(xaxis = list(zeroline = FALSE,title = SVA_x),
                 yaxis = list(zeroline = FALSE,title = SVA_y),
                 font = list(size = fontSize)) %>%
          config(
            toImageButtonOptions = list(
              format = "svg"
            )
          )
        p4
      #} else {
      #  p4 <- plot_ly() %>%
      #    add_trace(
      #      data = plot_df,
      #      x = ~X,
      #      y = ~Y,
      #      type = "scatter",
      #      mode = "markers",
      #      marker = list(color = "black"),
      #      marker=list(size=dotSize,symbol = 'circle'),
      #      hovertemplate = ~text
      #    ) %>%
      #    layout(xaxis = list(zeroline = FALSE,title = SVA_x),
      #           yaxis = list(zeroline = FALSE,title = SVA_y),
      #           font = list(size = fontSize)) %>%
      #    config(
      #      toImageButtonOptions = list(
      #        format = "svg"
      #      )
      #    )
      #  p4
      #}
    }
    
  })
  
  output$corrected_SVA_scatter <- renderPlotly({
    
    req(corrected_SVA_scatter_react())
    p <- corrected_SVA_scatter_react()
    p
    
  })
  output$corrected_SVA_scatter_df <- DT::renderDataTable({
    
    req(corrected_SVA_scatter_df())
    df <- corrected_SVA_scatter_df()
    df <- df[,-ncol(df)]
    colnames(df)[1] <- "SampleName"
    DT::datatable(df,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F) %>%
      formatRound(columns = grep("^SVA_",colnames(df),value = T), digits = 4)
    
  })
  
  # Box Plot ------------------------------------------------------------------
  output$BP1title <- renderUI({
    req(matrix1Title_react())
    h3(paste(matrix1Title_react()),"Box Plot")
  })
  output$BP2title <- renderUI({
    req(matrix2Title_react())
    h3(paste(matrix2Title_react()),"Box Plot")
  })
  output$rendBPsampSubset <- renderUI({
    
    #FeatChoices <- c("Select All Samples",colnames(aligned_meta_file())[-1])
    FeatChoices <- c("Select All Samples",colnames(boxplot_meta_file())[-1])
    selectInput("BPsampSubset","Subset Samples By:",choices = FeatChoices, selected = FeatChoices[1])
    
  })
  
  output$rendBPsampCriteria <- renderUI({
    
    req(input$BPsampSubset)
    subSelect <- input$BPsampSubset
    if (subSelect != "Select All Samples") {
      #sampCrit <- unique(aligned_meta_file()[,subSelect])
      sampCrit <- unique(boxplot_meta_file()[,subSelect])
      selectInput("BPsampCrit","Sample Criteria:",choices = sampCrit)
    }
    
  })
  
  output$rendBPgroupCriteria <- renderUI({
    
    #GroupChoices <- colnames(aligned_meta_file())[-1]
    GroupChoices <- colnames(boxplot_meta_file())[-1]
    if (isTruthy(input$BPsampSubset)) {
      GroupChoices <- GroupChoices[which(GroupChoices!=input$BPsampSubset)]
    }
    selectInput("BPgroupCriteria","Grouping Criteria:",choices = GroupChoices, selected = GroupChoices[1])
    
  })
  
  output$rendBPgroupSelection <- renderUI({
    
    req(input$BPgroupCriteria)
    req(boxplot_meta_file())
    #meta <- aligned_meta_file()
    meta <- boxplot_meta_file()
    groupCrit <- input$BPgroupCriteria
    if (input$BPremoveSingles == T) {
      tab <- table(meta[,groupCrit])
      meta <- meta[meta[,groupCrit] %in% names(tab[tab >= 2]), ]
    }
    GroupSelec <- unique(meta[,groupCrit])
    selectInput("BPgroupSelection","Select Groups:",choices = GroupSelec, selected = GroupSelec, multiple = T)
    
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
    req(uncorrected_matrix())
    if (FeatCat == "Matrix Features") {
      Features <- rownames(uncorrected_matrix())
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
  
  ## Matrix 1 -----------------------------------------------------------------
  
  PCA_Proj_uncorr_Samples <- reactive({
    req(matrix1())
    if (input$uncorrected_panel == "box") {
      withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
        incProgress(0.5, detail = "Running PCA")
        mat_uncorr <- matrix1()
        #rownames(mat_uncorr) <- mat_uncorr[,1]
        #mat_uncorr <- mat_uncorr[,-1]
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
  
  ImmDeconv_uncorr_react <- reactive({
    
    req(input$ImmuneDeconvMethods)
    req(matrix1())
    deconvMethod <- input$ImmuneDeconvMethods
    mat <- matrix1()
    #rownames(mat) <- mat[,1]
    #mat <- mat[,-1]
    withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
    
    req(matrix1())
    mat <- matrix1()
    if (input$HumanOrMouse == "Mouse") {
      mat_conv <- MouseToHuman(mat,MM_HS_Conv,gene_col = FALSE)
      mat_conv
    } else {
      mat <- matrix1()
      mat
    }
    
  })
  
  output$renduncorrected_Box_plot <- renderUI({
    plotHeight <- input$BPplotHeight
    plotWidth <- input$BPplotWidth
    shinyjqui::jqui_resizable(shiny::plotOutput("uncorrected_Box_plot",height = plotHeight, width = plotWidth))
  })
  
  CohortBPPlot_df_react <- reactive({
    
    req(input$BPsampSubset)
    req(matrix1())
    req(input$BPFeatureCategory)
    req(input$BPFeatSelection)
    req(input$BPgroupCriteria)
    req(boxplot_meta_file())
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
      mat <- as.data.frame(matrix1())
      featdf <- mat[Feature,]
      #rownames(featdf) <- featdf[,1]
      #featdf <- featdf[,-1]
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
        #rownames(mat) <- mat[,1]
        #mat <- mat[,-1]
        mat <- as.matrix(mat)
        withProgress(message = paste("Processing",matrix1Title_react()), value = 0, {
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
        #meta
      }
      if (Feature %in% colnames(meta)) {
        meta <- meta %>%
          group_by(!!sym(groupCrit)) %>%
          mutate(Outlier = ifelse(is_outlier(!!sym(Feature)),TRUE,FALSE)) %>%
          as.data.frame()
        meta
      }
    }
  })
  
  
  
  CohortBPPlot_react <- reactive({
    
    req(CohortBPPlot_df_react())
    req(input$BPFeatureCategory)
    req(input$BPFeatSelection)
    req(input$BPgroupCriteria)
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
    hjust_orient <- 0.5
    vjust_orient <- 0.5
    axis_orient <- as.numeric(input$BPxAxisOrient)  # X-axis label orientation
    if (axis_orient == 90) {                          # Adjust hjust if orientation is 0
      hjust_orient <- 1                                # Initial hjust
      vjust_orient <- 0.5
    } else if (axis_orient == 45) {                        # Adjust hjust if orientation is 0
      hjust_orient <- 1
      vjust_orient <- 1
    }
    BPorder <- input$BPplotXaxOrder
    BPGroupSelect <- input$BPgroupSelection
    
    if (isTruthy(Feature)) {
      plotdf <- plotdf_full[,c(NameCol,groupCrit,Feature)]
      plotdf <- plotdf[which(plotdf[,groupCrit] %in% BPGroupSelect),]
      colnames(plotdf) <- c("SampleName","Group","Feature")
      plotdf[,"Group"] <- gsub("[_-]"," ",plotdf[,"Group"])
      plotdf[,"Group"] <- ifelse(nchar(plotdf[,"Group"])>15,str_wrap(gsub("[_-]"," ",plotdf[,"Group"]), width = 15),plotdf[,"Group"])
      
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
      barp <- barp + theme(axis.text.x = element_text(size = Xaxis_font,angle = axis_orient, hjust = hjust_orient, vjust = vjust_orient),
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
  
  output$uncorrected_Box_plot_df <- DT::renderDataTable({
    
    req(CohortBPPlot_df_react())
    df <- CohortBPPlot_df_react()
    DT::datatable(df,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
    
  })
  
  
  ## Matrix 2 -----------------------------------------------------------------
  
  PCA_Proj_corr_Samples <- reactive({
    req(matrix2())
    if (input$uncorrected_panel == "box") {
      withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
        incProgress(0.5, detail = "Running PCA")
        mat_corr <- matrix2()
        #rownames(mat_corr) <- mat_corr[,1]
        #mat_corr <- mat_corr[,-1]
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
    req(matrix2())
    deconvMethod <- input$ImmuneDeconvMethods
    mat <- matrix2()
    #rownames(mat) <- mat[,1]
    #mat <- mat[,-1]
    withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
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
    
    req(matrix2())
    mat <- matrix2()
    if (input$HumanOrMouse == "Mouse") {
      mat_conv <- MouseToHuman(mat,MM_HS_Conv,gene_col = FALSE)
      mat_conv
    } else {
      mat <- matrix2()
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
    req(matrix2())
    mat <- matrix2()
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
      mat <- as.data.frame(matrix2())
      featdf <- mat[Feature,]
      #rownames(featdf) <- featdf[,1]
      #featdf <- featdf[,-1]
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
        #rownames(mat) <- mat[,1]
        #mat <- mat[,-1]
        mat <- as.matrix(mat)
        withProgress(message = paste("Processing",matrix2Title_react()), value = 0, {
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
        
        #meta
      }
      if (Feature %in% colnames(meta)) {
        meta <- meta %>%
          group_by(!!sym(groupCrit)) %>%
          mutate(Outlier = ifelse(is_outlier(!!sym(Feature)),TRUE,FALSE)) %>%
          as.data.frame()
        meta
      }
    }
    
  })
  
  CohortBPPlot_corr_react <- reactive({
    
    req(CohortBPPlot_corr_df_react())
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
    hjust_orient <- 0.5
    vjust_orient <- 0.5
    axis_orient <- as.numeric(input$BPxAxisOrient)  # X-axis label orientation
    if (axis_orient == 90) {                          # Adjust hjust if orientation is 0
      hjust_orient <- 1                                # Initial hjust
      vjust_orient <- 0.5
    } else if (axis_orient == 45) {                        # Adjust hjust if orientation is 0
      hjust_orient <- 1
      vjust_orient <- 1
    }
    BPorder <- input$BPplotXaxOrder
    BPGroupSelect <- input$BPgroupSelection
    
    if (isTruthy(Feature)) {
      plotdf <- plotdf_full[,c(NameCol,groupCrit,Feature)]
      plotdf <- plotdf[which(plotdf[,groupCrit] %in% BPGroupSelect),]
      colnames(plotdf) <- c("SampleName","Group","Feature")
      plotdf[,"Group"] <- gsub("[_-]"," ",plotdf[,"Group"])
      plotdf[,"Group"] <- ifelse(nchar(plotdf[,"Group"])>15,str_wrap(gsub("[_-]"," ",plotdf[,"Group"]), width = 15),plotdf[,"Group"])
      
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
      barp <- barp + theme(axis.text.x = element_text(size = Xaxis_font, angle = axis_orient, hjust = hjust_orient, vjust = vjust_orient),
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
  
  output$corrected_Box_plot_df <- DT::renderDataTable({
    
    req(CohortBPPlot_corr_df_react())
    df <- CohortBPPlot_corr_df_react()
    DT::datatable(df,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  rownames = F)
    
  })
  
  
  # File Download -------------------------------------------------------------
  
  ## Matrix 1 -------------------------------------------------------------
  # Saving files to ZIP
  #dir.create(temp_directory)
  #observe({
  ##shiny::observeEvent(input$save_uncorrected_matrix, {
  #  #withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
  #    file_name <- paste(gsub(" ","",matrix1Title_react()),"_matrix", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
  #    df <- as.data.frame(matrix1())
  #    df <- cbind(data.frame(Genes = rownames(df)),
  #                df)
  #    readr::write_tsv(df, file.path(temp_directory, file_name))
  #    print(paste(matrix1Title_react(),"matrix Ready for Zip"))
  #    #incProgress(1, detail = "Complete!")
  #  #})
  #  
  #  #print("Uncorrected_matrix Ready for Zip")
  #})
  #shiny::observeEvent(input$save_uncorrected_matrixZip, {
  #  withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
  #    file_name <- paste(gsub(" ","",matrix1Title_react()),"_matrix", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
  #    df <- as.data.frame(matrix1())
  #    df <- cbind(data.frame(Genes = rownames(df)),
  #                df)
  #    readr::write_tsv(df, file.path(temp_directory, file_name))
  #    print(paste(matrix1Title_react(),"matrix Ready for Zip"))
  #    incProgress(1, detail = "Complete!")
  #  })
  #  
  #  #print("Uncorrected_matrix Ready for Zip")
  #})
  output$dnldsave_uncorrected_matrix <- downloadHandler(
    filename = function() {
      paste0(matrix1DlndTitle_react(),"_Matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- as.data.frame(matrix1())
      df <- cbind(data.frame(Genes = rownames(df)),
                  df)
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  output$dnldsave_uncorrected_matrixRaw <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_Input_RNASeqCounts_Matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- as.data.frame(matrix1())
      df <- cbind(data.frame(Genes = rownames(df)),
                  df)
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  #output$dnldsave_uncorrected_matrixZip <- downloadHandler(
  #  filename = function() {
  #    paste0("BatchFlex_",gsub(" ","",matrix1Title_react()),"_Matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
  #  },
  #  content = function(file) {
  #    df <- as.data.frame(matrix1())
  #    df <- cbind(data.frame(Genes = rownames(df)),
  #                df)
  #    write.table(df,file,sep = '\t', row.names = F)
  #  }
  #)
  shiny::observeEvent(input$save_uncorrected_PCA_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_PCA_plot",input$save_uncorrected_PCA_plot, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_PCA_plot_forDlnd_react(),
                      height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
      #print(paste(matrix1DlndTitle_react(),"PCA_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_PCA_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_PCA","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_PCA_plot_forDlnd_react()
      ggplot2::ggsave(file,p, height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_PCA_mc_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_PCA_mc_plot",input$save_uncorrected_PCA_mc_plot, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_PCA_multiple_components(),
                      height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
      #print(paste(matrix1DlndTitle_react(),"PCA_mc_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_PCA_mc_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_PCA_MultipleComponents","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_PCA_multiple_components()
      ggplot2::ggsave(file,p, height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_scree_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_Scree_plot",input$save_uncorrected_scree_plot, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_scree_plot_react(),
                      height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
      #print(paste(matrix1DlndTitle_react(),"Scree_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_scree_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Scree_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_scree_plot_react()
      ggplot2::ggsave(file,p, height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_PCA_components, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_PC_components",input$save_uncorrected_PCA_components, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(as.data.frame(uncorrected_PCA_details2()$x), file.path(temp_directory, file_name))
      #print(paste(matrix1DlndTitle_react(),"PCA_details Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnld <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_PC_Components","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- save_uncorrected_PCA_components()
      write.table(as.data.frame(df$x),file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorrected_contribution_table, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_contribution_table",input$save_uncorrected_contribution_table, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(uncorrected_PCA_individuals(), file.path(temp_directory, file_name))
      #print(paste(matrix1DlndTitle_react(),"contribution_table Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_contribution_table <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Contribution_Table","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- uncorrected_PCA_individuals()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorrected_contribution_counts, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_contribution_counts",input$save_uncorrected_contribution_counts, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(uncorrected_contribution_counts(), file.path(temp_directory, file_name))
      #print(paste(matrix1DlndTitle_react(),"contribution_counts Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_contribution_counts <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Contribution_Counts","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- uncorrected_contribution_counts()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorrected_UMAP, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_UMAP",input$save_uncorrected_UMAP, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_UMAP_react(),
                      height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
      #print(paste(matrix1DlndTitle_react(),"UMAP Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_UMAP <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_UMAP","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_UMAP_react()
      ggplot2::ggsave(file,p, height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
    }
  )
  observeEvent(input$save_uncorrected_elbow_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_elbow_plot",input$save_uncorrected_elbow_plot, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_elbow_analysis(),
                      height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
      #print(paste(matrix1DlndTitle_react(),"elbow_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_elbow_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Elbow_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_elbow_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_uncorrected_silhouette_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_silhouette_plot",input$save_uncorrected_silhouette_plot, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_silhouette_analysis(),
                      height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
      #print(paste(matrix1DlndTitle_react(),"silhouette_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_silhouette_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Silhouette_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_silhouette_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_uncorrected_dunn_index_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_dunn_index_plot",input$save_uncorrected_dunn_index_plot, "_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_dunn_index_analysis(),
                      height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
      #print(paste(matrix1DlndTitle_react(),"dunn_index_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_dunn_index_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Dunn_Index_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_dunn_index_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_heatmap, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(temp_directory, "/", gsub(" ","",matrix1DlndTitle_react()),"_heatmap",input$save_uncorrected_heatmap,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      svg(filename = file_name, height = input$heatmapHeight, width = input$heatmapWidth)
      ComplexHeatmap::draw(uncorrected_heatmap())
      dev.off()
      #print(paste(matrix1DlndTitle_react(),"heatmap Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  
  output$dnldsave_uncorrected_heatmap <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Heatmap","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      svg(filename = file, height = input$heatmapHeight, width = input$heatmapWidth)
      ComplexHeatmap::draw(uncorrected_heatmap())
      dev.off()
    }
  )
  shiny::observeEvent(input$save_ClusterInfoTab, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste("Meta_Cluster_Info","_",format(Sys.time(),format = "%Y%m%d_%H%M"),input$save_ClusterInfoTab,".tsv", sep = "")
      readr::write_tsv(boxplot_meta_file(), file.path(temp_directory, file_name))
      #print("Cluster_Info Ready for Zip")
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_ClusterInfoTab <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Meta_Cluster_Info","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- boxplot_meta_file()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  
  observeEvent(input$save_uncorrected_barplot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_Bar_plot",input$save_uncorrected_barplot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggsave(filename = file_name, path = temp_directory, plot = uncorr_bar_plot_react(),
             height = input$barPHeight, width = input$barPWidth, units = input$barPUnits)
      #print(paste(matrix1DlndTitle_react(),"Bar_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_barplot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Bar_plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorr_bar_plot_react()
      ggplot2::ggsave(file,p, height = input$barPHeight, width = input$barPWidth, units = input$barPUnits)
    }
  )
  
  shiny::observeEvent(input$save_uncorr_avg_het_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_Average_Heterogeneity_",input$save_uncorr_avg_het_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(uncorr_avg_het_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix1DlndTitle_react()," Average_Heterogeneity Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorr_avg_het_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix1DlndTitle_react()),"_Average_Heterogeneity_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- uncorr_avg_het_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorr_het_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_Heterogeneity_",input$save_uncorr_het_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(uncorr_het_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix1DlndTitle_react()," Heterogeneity Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorr_het_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix1DlndTitle_react()),"_Heterogeneity_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- uncorr_het_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  
  shiny::observeEvent(input$save_uncorr_avg_evn_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_Average_Eveness_",input$save_uncorr_avg_evn_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(uncorr_avg_evn_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix1DlndTitle_react()," Average_Eveness Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorr_avg_evn_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix1DlndTitle_react()),"_Average_Eveness_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- uncorr_avg_evn_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_uncorr_evn_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_Eveness_",input$save_uncorr_evn_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(uncorr_evn_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix1DlndTitle_react()," Eveness Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorr_evn_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix1DlndTitle_react()),"_Eveness_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- uncorr_evn_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  
  shiny::observeEvent(input$save_uncorrected_RLE_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_RLE_plot",input$save_uncorrected_RLE_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_RLE(),
                      height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
      #print(paste(matrix1DlndTitle_react(),"RLE_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_RLE_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_RLE_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_RLE()
      ggplot2::ggsave(file,p, height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_EV_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_EV_plot",input$save_uncorrected_EV_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = uncorrected_EV(),
                      height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
      #print(paste(matrix1DlndTitle_react(),"EV_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_EV_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_EV_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_EV()
      ggplot2::ggsave(file,p, height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_pvca_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_PVCA_Plot",input$save_uncorrected_pvca_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = pvca_uncorr_react(),
                      height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
      #print(paste(matrix1DlndTitle_react(),"PVCA_Plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_pvca_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_PVCA_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- pvca_uncorr_react()
      ggplot2::ggsave(file,p, height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
    }
  )
  #observeEvent(input$save_uncorrected_SVA_probability_density, {
  #  withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
  #    file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_SVA_probability_density",input$save_uncorrected_SVA_probability_density,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
  #    ggsave(filename = file_name, path = temp_directory, plot = uncorrected_SVA_probability_ggplot(),
  #           height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
  #    #print(paste(matrix1DlndTitle_react(),"SVA_probability_density Ready for Zip"))
  #    incProgress(1, detail = "Complete!")
  #  })
  #})
  #output$dnldsave_uncorrected_SVA_probability_density <- downloadHandler(
  #  filename = function() {
  #    paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_SVA_Probability_Density","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
  #  },
  #  content = function(file) {
  #    p <- uncorrected_SVA_probability_ggplot()
  #    ggplot2::ggsave(file,p, height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
  #  }
  #)
  observeEvent(input$save_uncorrected_SVA_probability_density, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_SVA_ScatterPlot",input$save_uncorrected_SVA_probability_density,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggsave(filename = file_name, path = temp_directory, plot = uncorrected_SVA_scatter_react(),
             height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
      #print(paste(matrix1DlndTitle_react(),"SVA_probability_density Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_SVA_probability_density <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_SVA_ScatterPlot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- uncorrected_SVA_scatter_react()
      ggplot2::ggsave(file,p, height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
    }
  )
  observeEvent(input$save_uncorrected_SVA_scatter_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_SurrogateVars_",input$save_uncorrected_SVA_scatter_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      df <- uncorrected_SVA_scatter_df()
      df <- df[,-ncol(df)]
      colnames(df)[1] <- "SampleName"
      write.table(df,file.path(temp_directory, file_name), sep = '\t', row.names = F)
      #print(paste(matrix1DlndTitle_react(),"SVA_probability_density Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_SVA_scatter_df <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_SurrogateVars_","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".txt")
    },
    content = function(file) {
      df <- uncorrected_SVA_scatter_df()
      df <- df[,-ncol(df)]
      colnames(df)[1] <- "SampleName"
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  observeEvent(input$save_SVA_surrogate_variables, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste("Meta_withSV_",input$save_SVA_surrogate_variables,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(boxplot_meta_file(), file.path(temp_directory, file_name))
      #print("surrogate_variables Ready for Zip")
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_SVA_surrogate_variables <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_","Meta_withSV","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- boxplot_meta_file()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  observeEvent(input$save_uncorrected_Box_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_Box_plot",input$save_uncorrected_Box_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggsave(filename = file_name, path = temp_directory, plot = CohortBPPlot_react(),
             height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
      #print(paste(matrix1DlndTitle_react(),"Box_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_Box_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_Box_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- CohortBPPlot_react()
      ggplot2::ggsave(file,p, height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
    }
  )
  shiny::observeEvent(input$save_uncorrected_Box_plot_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix1DlndTitle_react()),"_BoxPlot_data",input$save_uncorrected_Box_plot_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(CohortBPPlot_df_react(), file.path(temp_directory, file_name))
      #print(paste(matrix1DlndTitle_react(),"BoxPlot_data Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_uncorrected_Box_plot_df <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix1DlndTitle_react()),"_BoxPlot_data","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- CohortBPPlot_df_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  ## Matrix 2 ---------------------------------------------------------------
  #shiny::observeEvent(input$save_corrected_matrix, {
  #  withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
  #    file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_matrix", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
  #    df <- as.data.frame(matrix2())
  #    df <- cbind(data.frame(Genes = rownames(df)),
  #                df)
  #    readr::write_tsv(df, file.path(temp_directory, file_name))
  #    print(paste(matrix2DlndTitle_react(),"matrix Ready for Zip"))
  #    incProgress(1, detail = "Complete!")
  #  })
  #})
  #shiny::observeEvent(input$save_corrected_matrixZip, {
  #  withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
  #    file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_matrix", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
  #    df <- as.data.frame(matrix2())
  #    df <- cbind(data.frame(Genes = rownames(df)),
  #                df)
  #    readr::write_tsv(df, file.path(temp_directory, file_name))
  #    print(paste(matrix2DlndTitle_react(),"matrix Ready for Zip"))
  #    incProgress(1, detail = "Complete!")
  #  })
  #})
  output$dnldsave_corrected_matrix <- downloadHandler(
    filename = function() {
      paste0(matrix2DlndTitle_react(),"_Matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- as.data.frame(matrix2())
      df <- cbind(data.frame(Genes = rownames(df)),
                  df)
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  output$dnldsave_corrected_matrixRaw <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_RNASeqCounts_",input$batch_correction_method2,"_Corrected_Matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- as.data.frame(matrix2())
      df <- cbind(data.frame(Genes = rownames(df)),
                  df)
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  #output$dnldsave_corrected_matrixZip <- downloadHandler(
  #  filename = function() {
  #    paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Matrix","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
  #  },
  #  content = function(file) {
  #    df <- as.data.frame(matrix2())
  #    df <- cbind(data.frame(Genes = rownames(df)),
  #                df)
  #    write.table(df,file,sep = '\t', row.names = F)
  #  }
  #)
  shiny::observeEvent(input$save_corrected_PCA_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_PCA_plot",input$save_corrected_PCA_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_PCA_plot_forDlnd_react(),
                      height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
      #print(paste(matrix2DlndTitle_react(),"PCA_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_PCA_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_PCA_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_PCA_plot_forDlnd_react()
      ggplot2::ggsave(file,p, height = input$pcaHeight, width = input$pcaWidth, units = input$pcaUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_PCA_mc_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_PCA_mc_plot",input$save_corrected_PCA_mc_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_PCA_multiple_components(),
                      height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
      #print(paste(matrix2DlndTitle_react(),"PCA_mc_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_PCA_mc_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_PCA_Multiple_Components_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_PCA_multiple_components()
      ggplot2::ggsave(file,p, height = input$pca_mcHeight, width = input$pca_mcWidth, units = input$pca_mcUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_scree_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_Scree_plot",input$save_corrected_scree_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_scree_plot_react(),
                      height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
      #print(paste(matrix2DlndTitle_react(),"Scree_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_scree_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Scree_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_scree_plot_react()
      ggplot2::ggsave(file,p, height = input$pca_dtHeight, width = input$pca_dtWidth, units = input$pca_dtUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_PCA_components, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_PC_components",input$save_corrected_PCA_components,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(as.data.frame(corrected_PCA_details2()$x), file.path(temp_directory, file_name))
      #print(paste(matrix2DlndTitle_react(),"PCA_details Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_PCA_components <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_PC_Components","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- corrected_PCA_details2()
      write.table(as.data.frame(df$x),file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corrected_contribution_table, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_contribution_table",input$save_corrected_contribution_table,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(corrected_PCA_individuals(), file.path(temp_directory, file_name))
      #print(paste(matrix2DlndTitle_react(),"contribution_table Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_contribution_table <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Contribution_Table","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- corrected_PCA_individuals()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corrected_contribution_counts, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_contribution_counts",input$save_corrected_contribution_counts,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(corrected_contribution_counts(), file.path(temp_directory, file_name))
      #print(paste(matrix2DlndTitle_react(),"contribution_counts Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_contribution_counts <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Contribution_Counts","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- corrected_contribution_counts()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corrected_UMAP, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_UMAP",input$save_corrected_UMAP,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_UMAP_react(),
                      height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
      #print(paste(matrix2DlndTitle_react(),"UMAP Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_UMAP <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_UMAP","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_UMAP_react()
      ggplot2::ggsave(file,p, height = input$umapHeight, width = input$umapWidth, units = input$umapUnits)
    }
  )
  observeEvent(input$save_corrected_elbow_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_elbow_plot",input$save_corrected_elbow_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_elbow_analysis(),
                      height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
      #print(paste(matrix2DlndTitle_react(),"elbow_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_elbow_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Elbow_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_elbow_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_corrected_silhouette_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_silhouette_plot",input$save_corrected_silhouette_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_silhouette_analysis(),
                      height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
      #print(paste(matrix2DlndTitle_react(),"silhouette_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_silhouette_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Silhouette_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_silhouette_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  observeEvent(input$save_corrected_dunn_index_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_dunn_index_plot",input$save_corrected_dunn_index_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_dunn_index_analysis(),
                      height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
      #print(paste(matrix2DlndTitle_react(),"dunn_index_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_dunn_index_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Dunn_Index_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_dunn_index_analysis()
      ggplot2::ggsave(file,p, height = input$cluster_mainHeight, width = input$cluster_mainWidth, units = input$cluster_mainUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_heatmap, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(temp_directory,"/",gsub(" ","",matrix2DlndTitle_react()),"_heatmap",input$save_corrected_heatmap,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      svg(filename = file_name, height = input$heatmapHeight, width = input$heatmapWidth)
      ComplexHeatmap::draw(corrected_heatmap())
      dev.off()
      #print(paste(matrix2DlndTitle_react(),"heatmap Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_heatmap <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Heatmap","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      svg(filename = file, height = input$heatmapHeight, width = input$heatmapWidth)
      ComplexHeatmap::draw(corrected_heatmap())
      dev.off()
    }
  )
  
  observeEvent(input$save_corrected_barplot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_Bar_plot",input$save_corrected_barplot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggsave(filename = file_name, path = temp_directory, plot = corr_bar_plot_react(),
             height = input$barPHeight, width = input$barPWidth, units = input$barPUnits)
      #print(paste(matrix2DlndTitle_react(),"Bar_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_barplot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Bar_plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corr_bar_plot_react()
      ggplot2::ggsave(file,p, height = input$barPHeight, width = input$barPWidth, units = input$barPUnits)
    }
  )
  
  shiny::observeEvent(input$save_corr_avg_het_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_Average_Heterogeneity_",input$save_corr_avg_het_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(corr_avg_het_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix2DlndTitle_react()," Average_Heterogeneity Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corr_avg_het_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix2DlndTitle_react()),"_Average_Heterogeneity_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- corr_avg_het_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corr_het_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_Heterogeneity_",input$save_corr_het_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(corr_het_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix2DlndTitle_react()," Heterogeneity Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corr_het_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix2DlndTitle_react()),"_Heterogeneity_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- corr_het_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  
  shiny::observeEvent(input$save_corr_avg_evn_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_Average_Eveness_",input$save_corr_avg_evn_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(corr_avg_evn_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix2DlndTitle_react()," Average_Eveness Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corr_avg_evn_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix2DlndTitle_react()),"_Average_Eveness_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- corr_avg_evn_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  shiny::observeEvent(input$save_corr_evn_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_Eveness_",input$save_corr_evn_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(corr_evn_react(), file.path(temp_directory, file_name))
      #print(paste0(matrix2DlndTitle_react()," Eveness Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corr_evn_df <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",matrix2DlndTitle_react()),"_Eveness_", format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
    },
    content = function(file) {
      df <- corr_evn_react()
      write.table(df,file,sep = '\t', row.names = F)
    }
  )
  
  shiny::observeEvent(input$save_corrected_RLE_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_RLE_plot",input$save_corrected_RLE_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_RLE(),
                      height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
      #print(paste(matrix2DlndTitle_react(),"RLE_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_RLE_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_RLE_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_RLE()
      ggplot2::ggsave(file,p, height = input$rleHeight, width = input$rleWidth, units = input$rleUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_EV_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_EV_plot",input$save_corrected_EV_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_EV(),
                      height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
      #print(paste(matrix2DlndTitle_react(),"EV_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_EV_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_EV_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_EV()
      ggplot2::ggsave(file,p, height = input$exp_varHeight, width = input$exp_varWidth, units = input$exp_varUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_pvca_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_PVCA_Plot",input$save_corrected_pvca_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = pvca_corr_react(),
                      height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
      #print(paste(matrix2DlndTitle_react(),"PVCA_Plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_pvca_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_PVCA_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- pvca_corr_react()
      ggplot2::ggsave(file,p, height = input$pvcaHeight, width = input$pvcaWidth, units = input$pvcaUnits)
    }
  )
  #observeEvent(input$save_corrected_SVA_probability_density, {
  #  withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
  #    file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_SVA_probability_density",input$save_corrected_SVA_probability_density,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
  #    ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_SVA_probability_ggplot(),
  #                    height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
  #    #print(paste(matrix2DlndTitle_react(),"SVA_probability_density Ready for Zip"))
  #    incProgress(1, detail = "Complete!")
  #  })
  #})
  #output$dnldsave_corrected_SVA_probability_density <- downloadHandler(
  #  filename = function() {
  #    paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_SVA_Probability_Density","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
  #  },
  #  content = function(file) {
  #    p <- corrected_SVA_probability_ggplot()
  #    ggplot2::ggsave(file,p, height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
  #  }
  #)
  observeEvent(input$save_corrected_SVA_probability_density, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_SVA_ScatterPlot",input$save_corrected_SVA_probability_density,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = corrected_SVA_scatter_react(),
                      height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
      #print(paste(matrix2DlndTitle_react(),"SVA_probability_density Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_SVA_probability_density <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_SVA_ScatterPlot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- corrected_SVA_scatter_react()
      ggplot2::ggsave(file,p, height = input$svaHeight, width = input$svaWidth, units = input$svaUnits)
    }
  )
  observeEvent(input$save_corrected_SVA_scatter_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_SurrogateVars_",input$save_corrected_SVA_scatter_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      df <- corrected_SVA_scatter_df()
      df <- df[,-ncol(df)]
      colnames(df)[1] <- "SampleName"
      write.table(df,file.path(temp_directory, file_name), sep = '\t', row.names = F)
      #print(paste(matrix1DlndTitle_react(),"SVA_probability_density Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_SVA_scatter_df <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_SurrogateVars_","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".txt")
    },
    content = function(file) {
      df <- corrected_SVA_scatter_df()
      df <- df[,-ncol(df)]
      colnames(df)[1] <- "SampleName"
      write.table(df,file, sep = '\t', row.names = F)
    }
  )
  observeEvent(input$save_corrected_Box_plot, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_Box_plot",input$save_corrected_Box_plot,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg", sep = "")
      ggplot2::ggsave(filename = file_name, path = temp_directory, plot = CohortBPPlot_corr_react(),
                      height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
      #print(paste(matrix2DlndTitle_react(),"Box_plot Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_Box_plot <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_Box_Plot","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".svg")
    },
    content = function(file) {
      p <- CohortBPPlot_corr_react()
      ggplot2::ggsave(file,p, height = input$BPHeight, width = input$BPWidth, units = input$BPUnits)
    }
  )
  shiny::observeEvent(input$save_corrected_Box_plot_df, {
    withProgress(message = paste("Adding to Step-3 zip file"), value = 0, {
      file_name <- paste(gsub(" ","",matrix2DlndTitle_react()),"_BoxPlot_data",input$save_corrected_Box_plot_df,"_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv", sep = "")
      readr::write_tsv(CohortBPPlot_corr_df_react(), file.path(temp_directory, file_name))
      #print(paste(matrix2DlndTitle_react(),"BoxPlot_data Ready for Zip"))
      incProgress(1, detail = "Complete!")
    })
  })
  output$dnldsave_corrected_Box_plot_df <- downloadHandler(
    filename = function() {
      paste0("BatchFlex_",gsub(" ","",matrix2DlndTitle_react()),"_BoxPlot_data","_",format(Sys.time(),format = "%Y%m%d_%H%M"),".tsv")
    },
    content = function(file) {
      df <- CohortBPPlot_corr_df_react()
      write.table(df,file,sep = '\t', row.names = F)
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
    #print(my_files())
    updatePickerInput(session, "select_save_files", choices = my_files(),
                      selected = my_files())
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
      paste(input$file_name, "_", format(Sys.time(),format = "%Y%m%d_%H%M"), ".zip", sep = "")
      #if (input$file_name == "BatchFlex" | !isTruthy(input$file_name)) {
      #  paste(input$file_name, "_", format(Sys.time(),format = "%Y%m%d_%H%M"), ".zip", sep = "")
      #} else {
      #  paste(input$file_name, "_", format(Sys.time(),format = "%Y%m%d"), ".zip", sep = "")
      #}
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