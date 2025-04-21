# Use the official R Shiny image as the base
FROM rocker/shiny:4.2.0

# Set the working directory inside the container
WORKDIR /srv/shiny-server

# Install system dependencies (required for some R packages like curl, libcurl, etc.)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R dependencies from CRAN and Bioconductor
# These libraries are used in your web app
RUN R -e "install.packages(c('shinydashboard', 'shinydashboardPlus', 'DT', 'dplyr', 'readr', 'shinyjs', 'plotly', 'ggplot2', 'ggrepel', 'EnhancedVolcano', 'pheatmap', 'RColorBrewer', 'shinyWidgets', 'shinyalert', 'tibble'))" && \
    R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('DESeq2', 'limma', 'biomaRt', 'clusterProfiler', 'org.Hs.eg.db', 'org.Mm.eg.db'))"

# Copy the entire app directory (including all subfolders like R/, www/, images/) into the container
COPY . .

# Expose the port Shiny uses to communicate with the outside world
EXPOSE 3838

# Make sure the entrypoint is to run the Shiny app (app.R)
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app.R', host = '0.0.0.0', port = 3838)"]
