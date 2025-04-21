# RNAlytics/R/global_config.R

# Load required R packages
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(DESeq2)
library(limma)
library(dplyr)
library(readr)
library(shinyjs)
library(plotly)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(shinyWidgets)
library(BiocManager)
library(tibble)
library(shinyalert)
library(biomaRt)  
library(clusterProfiler)  
library(org.Hs.eg.db)  
library(org.Mm.eg.db)  
options(repos = BiocManager::repositories())

# Global constants
APP_VERSION <- "v1.0.0"
DEVELOPER <- "Jash Trivedi, Azarian Lab, University of Central Florida"
RELEASE_DATE <- "April 02, 2025"
