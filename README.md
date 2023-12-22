# BATCH-FLEX
BatchFLEX: Comprehensive software for assessing and correcting batch effects, comparing batch correction methods, and exporting corrected matrices for downstream analysis

# Introduction
BatchFLEX is an intuitive Shiny App that can be viewed as a web interface through shiny.io or accessed locally by downloading and installing an app.R or docker container, or by installing an R-package. BatchFLEX is designed with a clean user interface and easy to follow steps to prompt users to input a data matrix and meta file, to select and implement the desired method of batch correction, to assess the impact of the batch correction method using side-by-side comparisons of pre and post corrected graphs and statistics, and finally to export a corrected matrix and accompanying diagnostics for downstream analysis. The intuitive user interface allows for hands-on evaluation and selection of the most optimal batch correction method for any dataset and for users with any background. BatchFLEX is the most comprehensive batch correction and evaluation tool and can easily be implemented into any pipeline using the R-package wrapper. Notable features of BatchFLEX include implementation of a wide variety of batch correction methods such as ComBat, ComBatSeq, Harman, LIMMA, RUVg, and Mean centering, incorporation of multiple different types of evaluation methods such as PCA, cluster analysis, gene-variable association analysis, scree plots, SVA, uMAP, RLE, boxplots, and heatmaps, and inclusion of a simple and comprehensive export feature for use of corrected matrices in downstream analysis and for record keeping to justify the correction method selection for publications. 

An overview of the BatchFLEX suite of tools can be found on our GitHub page (https://github.com/shawlab-moffitt/BATCHFLEX), which includes source code, example data, and an installation guide. Additionally, BatchFLEX can be accessed through shiny.io at (https://shawlab-moffitt.shinyapps.io/batch_flex/). 

# Installation
The BatchFLEX suite can be downloaded (cloned) and installed through the GitHub repository. The downloaded file can be unzipped to a destination folder, which should be set as the working directory or file path. Of note, some of the example files (e.g., gene set files) use relative paths, so the program may fail to identify the file if a working directory is not properly set. The suite was developed in R version 4.1. 
•	Install BatchFLEX suite GitHub repository 
** git clone https://github.com/shawlab-moffitt/BATCHFLEX

*	Download and unzip repository https://github.com/shawlab-moffitt/BATCHFLEX
**	Set working directory to BatchFLEX folder
**	Install required R packages
**	Suite of tools was built on R version 4.2
**	R script for package installation is provided in the “1-Getting_Started” folder




