# Use the rocker/r-ver image for R 4.3.1
FROM rocker/r-ver:4.3.2

# Maintainer information
LABEL maintainer="Tim Shaw <timothy.shaw@moffitt.org>"

# Install necessary system libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    libudunits2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    libssh2-1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libhdf5-dev \
    libfftw3-dev \
    libopenblas-dev \
    libatlas-base-dev \
    libmagick++-dev \
    libudunits2-dev \
    libgdal-dev \
    libgsl0-dev \
    libv8-dev \
    libzstd-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
  && rm -rf /var/lib/apt/lists/*

# Set library paths for R packages
ENV LD_LIBRARY_PATH=/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/lib:/usr/local/lib/R/lib


# Install required R packages
RUN R -e "install.packages(c('shiny','shinyjqui','limma','bapred','DT','shinycssloaders','metafolio','corrplot','readr','dplyr','tidyverse','zip','glue','factoextra','FactoMineR',
'magick','ggfortify','ggplotify','draw','clValid','shinyWidgets','tidyr','umap','plotly','statVisual','svglite','tabula','folio','bslib','readtext','shinythemes','uuid','stringr'))"


# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" \
  && R -e "BiocManager::install(c('sva', 'affyPLM', 'scater', 'Harman', 'RUVSeq', 'ComplexHeatmap', 'pvca', 'bapred', 'statVisual','GSVA','edgeR'), ask = FALSE)"

# Install 'remotes' package and use it to install 'omnideconv/immunedeconv' from GitHub
RUN R -e "install.packages('remotes')" \   
  && R -e "remotes::install_github('omnideconv/immunedeconv')"

# Install preprocessCore with specific configuration
RUN R -e "BiocManager::install('preprocessCore', configure.args = c(preprocessCore = '--disable-threading'), force = TRUE, update = TRUE, type = 'source')"

# Set the Shiny port and host
RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))"

# Create a user for running the app
RUN addgroup --system app && adduser --system --ingroup app app

# Set the working directory
WORKDIR /app

# Copy the R script files into the container at /app/
COPY app/www /app/www/
COPY app/R /app/R/
COPY app/Mouse_To_Human_Gene_Conversion /app/Mouse_To_Human_Gene_Conversion/
COPY app/GeneSet_Data /app/GeneSet_Data/
COPY app/Example_Data /app/Example_Data/
COPY app/default_housekeeping /app/default_housekeeping/
COPY app/Batchflex_Homepage /app/Batchflex_Homepage/
COPY app/app.R /app/

# Set user permissions
RUN chown app:app -R /app

# Expose the port that Shiny will run on
EXPOSE 3838

# Run the Shiny app when the container launches
USER app
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]