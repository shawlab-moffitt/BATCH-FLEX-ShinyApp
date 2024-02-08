
##------- BATCH-FLEX-Shiny Packages Installation ------- ##


# List of base packages
base_packages <- c(
  "base",
  "grid",
  "stats4",
  "stats",
  "graphics",
  "grDevices",
  "utils",
  "datasets",
  "methods"
)

# List of packages loaded after running the script
loaded_packages <- c(
  "pvca",
  "ComplexHeatmap",
  "GSVA",
  "RUVSeq",
  "EDSeq",
  "ShortRead",
  "GenomicAlignments",
  "Rsamtools",
  "Biostrings",
  "XVector",
  "Harman",
  "scater",
  "scuttle",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "GenomicRanges",
  "GenomeInfoDb",
  "IRanges",
  "S4Vectors",
  "MatrixGenerics",
  "matrixStats",
  "statVisual",
  "plotly",
  "umap",
  "edgeR",
  "EPIC",
  "shinyWidgets",
  "clValid",
  "draw",
  "ggplotify",
  "ggfortify",
  "magick",
  "FactoMineR",
  "factoextra",
  "glue",
  "zip",
  "lubridate",
  "forcats",
  "stringr",
  "purrr",
  "tidyr",
  "tibble",
  "ggplot2",
  "tidyverse",
  "dplyr",
  "readr",
  "corrplot",
  "metafolio",
  "shinycssloaders",
  "DT",
  "bapred",
  "affyPLM",
  "preprocessCore",
  "gcrma",
  "affy",
  "Biobase",
  "BiocGenerics",
  "sva",
  "BiocParallel",
  "genefilter",
  "mgcv",
  "nlme",
  "MASS",
  "lme4",
  "glmnet",
  "Matrix",
  "limma",
  "shinyjqui",
  "shiny",
  "svglite",
  "tabula",
  "folio"
)

# List of Bioconductor packages
bioc_packages <- c('sva', 'affyPLM', 'scater', 'Harman', 'RUVSeq', 'ComplexHeatmap', 'pvca', 'bapred', 'statVisual')

# List of GitHub packages
github_packages <- c('omnideconv/immunedeconv')

# Remove packages that will be installed from Bioconductor or GitHub
loaded_packages <- loaded_packages[!(loaded_packages %in% bioc_packages | loaded_packages %in% github_packages)]

# Install base packages if not installed
install_base_packages <- base_packages[!(base_packages %in% rownames(installed.packages()))]
if (length(install_base_packages) > 0) {
  install.packages(install_base_packages)
}

# Install loaded packages if not installed
install_loaded_packages <- loaded_packages[!(loaded_packages %in% rownames(installed.packages()))]
if (length(install_loaded_packages) > 0) {
  install.packages(install_loaded_packages)
}

# Check if Bioconductor Manager is installed, if not, install it
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
installed_bioc_packages <- bioc_packages %in% rownames(installed.packages())
if (any(!installed_bioc_packages)) {
  BiocManager::install(bioc_packages[!installed_bioc_packages], ask = FALSE)
}

# Check if 'remotes' package is installed, if not, install it
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# Install packages from GitHub
remotes::install_github(github_packages)

