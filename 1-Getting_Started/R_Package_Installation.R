####----BatchFlex Shiny App R Package Installation----####
# author: Alyssa Obermayer (alyssa.obermayer@moffitt.org)

# Documentation
#' R script for installation of all R packages needed to run the BatchFlex Shiny App of tools

##--This Will install immunedeconv package--##
##---R Version 4.1 or greater is required---##

R.major <- as.numeric(R.Version()$major)
R.minor <- as.numeric(R.Version()$minor)

if (R.major >= 4 & R.minor >= 1) {
  
  install.packages("remotes")
  library(remotes)
  remotes::install_github("omnideconv/immunedeconv")
  library(immunedeconv)
  
}


##--Other R packages--##

## Check if packages are installed
c("shiny","shinyjqui","limma","bapred","DT","shinycssloaders","metafolio","corrplot",
  "readr","dplyr","tidyverse","zip","glue","factoextra","FactoMineR","magick","ggfortify",
  "ggplotify","draw","clValid","shinyWidgets","tidyr","umap","plotly","statVisual",
  "svglite","tabula","folio","bslib","readtext","shinythemes","uuid","stringr")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))


##--Bioconductor Packages--##

## Check if Bioconductor specific packages are installed
bioCpacks <- c("GSVA","sva","preprocessCore","scater","affyPLM",
               "Harman","RUVSeq","ComplexHeatmap","pvca","edgeR")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))



