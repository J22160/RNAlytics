# RNAlytics/R/ui_components.R
# UI components for each tab in the Shiny app

source("R/utilities.R")

homeUI <- function() {
  fluidPage(
    tags$head(
      tags$style(HTML("
        /* Add custom global animation and styles */
        @keyframes fadeInUp {
          0% {
            opacity: 0;
            transform: translateY(30px);
          }
          100% {
            opacity: 1;
            transform: translateY(0);
          }
        }

        @keyframes fadeInText {
          0% {
            opacity: 0;
            transform: translateY(30px);
          }
          100% {
            opacity: 1;
            transform: translateY(0);
          }
        }

        .animated-card {
          animation: fadeInUp 1s ease-in-out;
          transition: transform 0.3s ease, box-shadow 0.3s ease;
          background-color: #ffffff;
          border-radius: 12px;
          box-shadow: 0 4px 8px rgba(0,0,0,0.1);
          padding: 15px;
          margin-bottom: 20px;
          text-align: center;
        }

        .animated-card:hover {
          transform: translateY(-6px) scale(1.05);
          box-shadow: 0 12px 24px rgba(0,0,0,0.15);
        }

        .feature-title {
          font-size: 22px;
          color: #e74c3c;
          margin-top: 10px;
        }

        .feature-text {
          font-size: 14px;
          color: #555555;
          margin-top: 10px;
        }

        h1, h2, h3 {
          font-family: 'Segoe UI', sans-serif;
        }

        .intro-text {
          font-size: 18px;
          color: #333;
          line-height: 1.6;
        }

        /* Style button with transition */
        .btn-primary {
          background-color: #3498db;
          border-color: #3498db;
          color: #ffffff;
          font-size: 16px;
          padding: 12px 25px;
          border-radius: 5px;
          transition: background-color 0.3s ease, transform 0.2s ease;
        }

        .btn-primary:hover {
          background-color: #5dade2;
          border-color: #5dade2;
          transform: translateY(-3px);
        }

        h1 {
          font-size: 98px;
          color: #2c3e50;
          margin-bottom: 15px;
          animation: fadeInText 4s ease-out;
        }

        h2 {
          font-size: 28px;
          color: #2980b9;
          margin-top: 20px;
          margin-bottom: 10px;
        }

        h3 {
          font-size: 24px;
          color: #34495e;
          margin-bottom: 10px;
        }

        /* Background color */
        body {
          background-color: #f7f1e1; /* light peach */
        }

        /* Responsive design */
        @media (max-width: 768px) {
          .feature-box {
            width: 45%;
          }
        }

        @media (max-width: 480px) {
          .feature-box {
            width: 100%;
            flex-direction: column;
            text-align: center;
          }

          .feature-image {
            margin-right: 0;
            margin-bottom: 10px;
          }

          h1 {
            font-size: 76px;
          }

          h2 {
            font-size: 24px;
          }

          h3 {
            font-size: 20px;
          }

          p {
            font-size: 16px;
          }

          .btn-primary {
            font-size: 14px;
            padding: 8px 16px;
          }
        }
      "))
    ),
    
    fluidRow(
      box(
        width = 12,
        status = "primary",
        solidHeader = FALSE,
        
        tags$h1("Welcome to RNAlytics", style = "color: #073667; text-align: center; font-size: 76px; animation: fadeInText 3s ease-out;"),
        tags$p("RNAlytics is a powerful open-source web app that makes RNA-seq data analysis fast, intuitive, and code-free. Designed for researchers and biologists, it streamlines the workflow from raw counts to differential expression and enrichment analysis.",
          class = "intro-text", style = "text-align: center;"
        ),
        
        br(),
        tags$h2("Why RNAlytics?", style = "color: #2980b9; font-size: 28px;"),
        tags$p(
          "RNA-seq analysis often involves complex steps—data cleaning, statistical testing, and biological interpretation...",
          class = "intro-text"
        ),
        tags$ul(
          tags$li("User-Friendly: No coding required—just upload your data and start analyzing."),
          tags$li("Comprehensive: Covers preprocessing, differential analysis, visualization, and gene/pathway enrichment."),
          tags$li("Flexible: Supports multiple species and statistical methods."),
          style = "font-size: 16px;"
        ),
        
        br(),
        tags$h2("Explore the Power of RNAlytics", style = "color: #2980b9; font-size: 28px;"),
        
        # Feature Grid Row 1
        tags$div(
          style = "display: flex; justify-content: space-around; flex-wrap: wrap;",
          
          tags$div(
            class = "animated-card", style = "width: 30%; text-align: center;",
            tags$img(src = "images/data_preprocessing.png", width = "80px", height = "80px"),
            tags$div(class = "feature-title", "Data Preprocessing"),
            tags$p("Filter and normalize your raw gene expression data with ease using DESeq2-based methods.", class = "feature-text")
          ),
          
          tags$div(
            class = "animated-card", style = "width: 30%; text-align: center;",
            tags$img(src = "images/ensembl_conversion.png", width = "80px", height = "80px"),
            tags$div(class = "feature-title", "Ensembl Conversion"),
            tags$p("Convert Ensembl IDs to gene symbols for seamless analysis across human and mouse datasets.", class = "feature-text")
          ),
          
          tags$div(
            class = "animated-card", style = "width: 30%; text-align: center;",
            tags$img(src = "images/gene_expression.png", width = "80px", height = "80px"),
            tags$div(class = "feature-title", "Differential Expression"),
            tags$p("Perform robust statistical analysis with DESeq2, limma, or t-tests, tailored to your needs.", class = "feature-text")
          )
        ),
        
        # Feature Grid Row 2
        tags$div(
          style = "display: flex; justify-content: space-around; flex-wrap: wrap;",
          
          tags$div(
            class = "animated-card", style = "width: 30%; text-align: center;",
            tags$img(src = "images/visualization.png", width = "80px", height = "80px"),
            tags$div(class = "feature-title", "Visualization"),
            tags$p("Generate interactive PCA, volcano, and heatmap plots to explore your data visually.", class = "feature-text")
          ),
          
          tags$div(
            class = "animated-card", style = "width: 30%; text-align: center;",
            tags$img(src = "images/kegg_enrichment.png", width = "80px", height = "80px"),
            tags$div(class = "feature-title", "KEGG Enrichment"),
            tags$p("Discover enriched KEGG pathways with detailed results and visualizations.", class = "feature-text")
          ),
          
          tags$div(
            class = "animated-card", style = "width: 30%; text-align: center;",
            tags$img(src = "images/go_enrichment.png", width = "80px", height = "80px"),
            tags$div(class = "feature-title", "GO Enrichment"),
            tags$p("Analyze Gene Ontology terms (BP, MF, CC) with detailed results and visualizations.", class = "feature-text")
          )
        ),
        
        br(),
        tags$h2("Key Benefits", style = "color: #2980b9; font-size: 28px;"),
        tags$ul(
          tags$li("Speed: Process large datasets quickly with optimized algorithms."),
          tags$li("Interactivity: Real-time plot customization and exploration."),
          tags$li("Reproducibility: Download results and plots for reports or publications."),
          tags$li("Support: Active development and community feedback welcomed."),
          style = "font-size: 16px;"
        ),
        
        br(), br(),
        tags$p(
          paste("Version:", APP_VERSION, "| Released:", RELEASE_DATE, "| Developed by:", DEVELOPER),
          style = "font-size: 14px; text-align: center;"
        ),
        tags$p(
          "Stay tuned for updates with more features and enhanced functionality!",
          style = "font-size: 14px; text-align: center;"
        )
      )
    )
  )
}

# Preprocessing tab UI
preprocessingUI <- function() {
  fluidRow(
    box(
      title = "Data Preprocessing",
      width = 3,
      solidHeader = TRUE,
      status = "primary",
      collapsible = TRUE,
      fileInput("rawFile", "Upload Raw Gene Expression Data (.txt)", accept = ".txt"),
      numericInput("minCount", "Minimum Count Threshold", value = 5, min = 0),
      checkboxInput("normalize", "Normalize Data (DESeq2)", value = TRUE),
      actionButton("preprocessData", "Preprocess Data", class = "btn-primary", icon = icon("filter")),
      br(), br(),
      downloadButton("downloadPreprocessed", "Download Preprocessed Data", disabled = TRUE)
    ),
    box(
      title = "Preprocessed Data Preview",
      width = 9,
      solidHeader = TRUE,
      status = "success",
      collapsible = TRUE,
      DTOutput("preprocessedTable")
    )
  )
}

# Ensembl Conversion tab UI
ensemblConversionUI <- function() {
  fluidRow(
    box(
      title = "Ensembl ID to Gene Symbol Conversion",
      width = 3,
      solidHeader = TRUE,
      status = "primary",
      collapsible = TRUE,
      fileInput("ensemblFile", "Upload Data with Ensembl IDs (.txt)", accept = ".txt"),
      selectInput("species", "Select Species", choices = c("Human" = "hsapiens_gene_ensembl", "Mouse" = "mmusculus_gene_ensembl"), selected = "hsapiens_gene_ensembl"),
      actionButton("convertEnsembl", "Convert IDs", class = "btn-primary", icon = icon("exchange-alt")),
      br(), br(),
      downloadButton("downloadConverted", "Download Converted Data", disabled = TRUE)
    ),
    box(
      title = "Converted Data Preview",
      width = 9,
      solidHeader = TRUE,
      status = "success",
      collapsible = TRUE,
      DTOutput("convertedTable")
    )
  )
}

# Analysis tab UI
analysisUI <- function() {
  fluidRow(
    box(
      title = "Input & Settings",
      width = 4,
      solidHeader = TRUE,
      status = "primary",
      collapsible = TRUE,
      fileInput("file", "Upload Gene Expression Data (.txt)", accept = ".txt"),
      fileInput("annotationFile", "Upload Annotation File (.txt)", accept = ".txt"),
      uiOutput("groupDefinitionsUI"),
      uiOutput("sampleAssignmentsUI"),
      uiOutput("comparisonSelectionUI"),
      selectInput("statTest", "Statistical Test", choices = c("DESeq2", "limma"), selected = "DESeq2"),
      numericInput("logFCThresh", "Log2 Fold Change Threshold", value = 2, min = 0, step = 0.1),
      numericInput("padjThresh", "Adjusted P-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
      actionButton("runAnalysis", "Run Analysis", class = "btn-primary", icon = icon("cogs")),
      br(), br(),
      downloadButton("downloadFullResults", "Full Results", disabled = TRUE),
      downloadButton("downloadFilteredResults", "Filtered Results", disabled = TRUE)
    ),
    box(
      title = "Analysis Results",
      width = 8,
      solidHeader = TRUE,
      status = "success",
      collapsible = TRUE,
      numericInput("filterValue", "Filter Value for log2FoldChange", value = 0, step = 0.1),
      selectInput("filterType", "Filter Type", choices = c("Greater Than or Equal To" = "gte", "Less Than or Equal To" = "lte"), selected = "gte"),
      DTOutput("resultsTable"),
      verbatimTextOutput("errorLog")
    )
  )
}

# Visualization tab UI
visualizationUI <- function() {
  fluidRow(
    box(
      width = 12,
      title = "Interactive Visualizations",
      status = "info",
      solidHeader = TRUE,
      collapsible = TRUE,
      tabsetPanel(
        id = "vizTabs",
        tabPanel("PCA Plot", fluidRow(column(9, plotOutput("pcaPlot", height = "800px", width = "1000px")), column(3, pcaCustomUI()))),
        tabPanel("Volcano Plot", fluidRow(column(9, plotOutput("volcanoPlot", height = "800px", width = "1000px")), column(3, volcanoCustomUI()))),
        tabPanel("Heatmap", fluidRow(column(9, plotOutput("heatmapPlot", height = "800px", width = "1000px")), column(3, heatmapCustomUI()))),
        tabPanel("Gene Count Plot", fluidRow(column(9, plotOutput("geneCountPlot", height = "600px")), column(3, geneCountCustomUI())))
      )
    )
  )
}


#PCA cutomization panel
pcaCustomUI <- function() {
  box(
    width = NULL,
    status = "warning",
    collapsible = TRUE,
    title = "PCA Customization",
    textInput("pcaTitle", "Title", "PCA of Normalized Gene Expression"),
    textInput("pcaXLabel", "X-axis Label", "PC1"),
    textInput("pcaYLabel", "Y-axis Label", "PC2"),
    pickerInput(
      "pcaColors", 
      "Condition Colors", 
      choices = c("Red", "Blue", "Green", "Purple", "Orange", "Yellow", "Cyan", "Magenta"), 
      selected = c("Red", "Blue", "Green"),  # Default to 3 colors
      multiple = TRUE,
      options = list(`min-options` = 2)  # Enforce at least 2 colors
    ),
    numericInput("pcaPointSize", "Point Size", value = 3, min = 1, step = 0.5),
    numericInput("pcaLabelSize", "Label Size", value = 3, min = 1, step = 0.5),
    numericInput("pcaFontSize", "Font Size", value = 14, min = 8, step = 1),
    downloadButton("downloadPCAPlot", "Download PCA", disabled = TRUE)
  )
}

# Volcano customization panel
volcanoCustomUI <- function() {
  box(
    width = NULL,
    status = "warning",
    collapsible = TRUE,
    title = "Volcano Customization",
    textInput("volcanoTitle", "Title", "Volcano Plot of Differential Expression"),
    textInput("volcanoSubtitle", "Subtitle", "Significant Genes Highlighted"),
    textInput("volcanoXLabel", "X-axis Label", "Log2 Fold Change"),
    textInput("volcanoYLabel", "Y-axis Label", "-log10 Adjusted P-value"),
    pickerInput("volcanoColors", "Colors (Not Sig, Sig, sig+LFC, sig+pval)", choices = c("Grey", "Red", "Blue", "Green", "Orange", "Purple", "Yellow"), selected = c("Grey", "Red", "Blue", "Green"), multiple = TRUE),
    numericInput("volcanoPointSize", "Point Size", value = 2, min = 1, step = 0.5),
    numericInput("volcanoLabelSize", "Label Size", value = 3, min = 1, step = 0.5),
    selectInput("volcanoLegendPos", "Legend Position", choices = c("right", "left", "top", "bottom"), selected = "right"),
    numericInput("volcanoFontSize", "Font Size", value = 14, min = 8, step = 1),
    selectInput("volcanoTopGenes", "Label Top Genes", choices = c("None" = 0, "Top 10" = 10, "Top 20" = 20, "Top 30" = 30), selected = 0),
    downloadButton("downloadVolcanoPlot", "Download Volcano", disabled = TRUE)
  )
}

# Heatmap customization panel
heatmapCustomUI <- function() {
  box(
    width = NULL,
    status = "warning",
    collapsible = TRUE,
    title = "Heatmap Customization",
    textInput("heatmapTitle", "Title", "Heatmap of Top DE Genes"),
    pickerInput("heatmapPalette", "Color Palette", choices = c("RdBu", "YlOrRd", "Greens", "Blues"), selected = "RdBu"),
    numericInput("heatmapFontSize", "Font Size", value = 14, min = 8, step = 1),
    selectInput("heatmapClusterMethod", "Clustering Method", choices = c("euclidean", "manhattan", "correlation"), selected = "euclidean"),
    checkboxInput("heatmapShowRownames", "Show Row Names", value = TRUE),
    downloadButton("downloadHeatmapPlot", "Download Heatmap", disabled = TRUE)
  )
}

geneCountCustomUI <- function() {
  box(
    width = NULL,
    status = "warning",
    collapsible = TRUE,
    title = "Gene Count Plot Customization",
    selectizeInput("geneCountGene", "Select Gene", choices = NULL),
    textInput("geneCountTitle", "Title", "Gene Expression Across Samples"),
    pickerInput(
      "geneCountColors", 
      "Condition Colors", 
      choices = c("Red", "Blue", "Green", "Purple", "Orange", "Yellow", "Cyan", "Magenta"), 
      selected = c("Red", "Blue", "Green"),  # Default to 3 colors
      multiple = TRUE,
      options = list(`min-options` = 2)
    ),
    numericInput("geneCountPointSize", "Point Size", value = 3, min = 1, step = 0.5),
    numericInput("geneCountFontSize", "Font Size", value = 14, min = 8, step = 1),
    downloadButton("downloadGeneCountPlot", "Download Gene Count Plot", disabled = TRUE)
  )
}

# KEGG Pathway Enrichment tab UI
pathwayEnrichmentUI <- function() {
  fluidRow(
    box(
      title = "KEGG Pathway Enrichment Analysis",
      width = 4,
      solidHeader = TRUE,
      status = "primary",
      collapsible = TRUE,
      selectInput("enrichSource", "Enrichment Source", choices = c("KEGG"), selected = "KEGG"),
      selectInput("enrichSpecies", "Select Species", choices = c("Human" = "hsa", "Mouse" = "mmu"), selected = "hsa"),
      numericInput("enrichPval", "P-value Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      actionButton("runEnrichment", "Run Enrichment", class = "btn-primary", icon = icon("cogs")),
      br(), br(),
      downloadButton("downloadEnrichResults", "Download Results", disabled = TRUE),
      downloadButton("downloadEnrichDotPlot", "Download Dot Plot", disabled = TRUE)
    ),
    box(
      title = "KEGG Enrichment Results",
      width = 8,
      solidHeader = TRUE,
      status = "success",
      collapsible = TRUE,
      DTOutput("enrichTable"),
      plotOutput("enrichDotPlot", height = "500px", width = "800px")
    )
  )
}

# GO Enrichment tab UI
goEnrichmentUI <- function() {
  fluidRow(
    box(
      title = "Gene Ontology Enrichment Analysis",
      width = 4,
      solidHeader = TRUE,
      status = "primary",
      collapsible = TRUE,
      selectInput("goSpecies", "Select Species", choices = c("Human" = "hsa", "Mouse" = "mmu"), selected = "hsa"),
      selectInput("goOntology", "Ontology", choices = c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC"), selected = "BP"),
      numericInput("goPval", "P-value Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      actionButton("runGOEnrichment", "Run GO Enrichment", class = "btn-primary", icon = icon("cogs")),
      br(), br(),
      downloadButton("downloadGOResults", "Download Results", disabled = TRUE),
      downloadButton("downloadGODotPlot", "Download Dot Plot", disabled = TRUE)
    ),
    box(
      title = "GO Enrichment Results",
      width = 8,
      solidHeader = TRUE,
      status = "success",
      collapsible = TRUE,
      DTOutput("goTable"),
      plotOutput("goDotPlot", height = "500px", width = "800px")
    )
  )
}

# Help tab UI
helpUI <- function() {
  fluidRow(
    box(
      width = 12,
      title = "User Guide",
      status = "info",
      solidHeader = TRUE,
      collapsible = TRUE,
      tags$h3("How to Use RNAlytics"),
      tags$p("1. Home: Learn about the app."),
      tags$p("2. Preprocessing: Clean and normalize your data."),
      tags$p("3. Ensembl Conversion: Convert Ensembl IDs to gene symbols."),
      tags$p("4. Analysis: Upload data, add annotations, and run tests."),
      tags$p("5. Visualization: Explore plots."),
      tags$p("6. KEGG Enrichment: Analyze KEGG pathways of DEGs."),
      tags$p("7. GO Enrichment: Analyze Gene Ontology terms of DEGs."),
      tags$p("For support, contact jashtrivedi221@gmail.com.")
    )
  )
}

# Main app UI with resource path for images
appUI <- function() {
  # Add resource path for images folder
  addResourcePath("images", "images")
  
  dashboardPage(
    skin = "blue-light",
    dashboardHeader(title = tags$span("RNAlytics", style = "font-size: 32px; font-weight: bold; color: #34495e;"), titleWidth = 300),
    dashboardSidebar(
      width = 250,
      sidebarMenu(
        id = "tabs",
        menuItem("Home", tabName = "home", icon = icon("house", lib = "font-awesome")),
        menuItem("Preprocessing", tabName = "preprocessing", icon = icon("filter", lib = "font-awesome")),
        menuItem("Ensembl Conversion", tabName = "ensembl_conversion", icon = icon("exchange-alt", lib = "font-awesome")),
        menuItem("Analysis", tabName = "analysis", icon = icon("vial", lib = "font-awesome")),
        menuItem("Visualization", tabName = "visualization", icon = icon("chart-line", lib = "font-awesome")),
        menuItem("KEGG Enrichment", tabName = "pathway_enrichment", icon = icon("sitemap", lib = "font-awesome")),
        menuItem("GO Enrichment", tabName = "go_enrichment", icon = icon("sitemap", lib = "font-awesome")),
        menuItem("Help", tabName = "help", icon = icon("info-circle", lib = "font-awesome"))
      ),
      tags$hr(),
      tags$div(
        style = "padding: 10px; color: #000000; font-size: 14px;",
        paste("Version:", APP_VERSION),
        br(),
        DEVELOPER
      ),
      tags$hr(),
      selectInput("comparisonHistory", "Comparison History", choices = c("None" = ""), selected = "", width = "100%")
    ),
    dashboardBody(
      useShinyjs(),
      useShinyalert(),
      tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
      tabItems(
        tabItem(tabName = "home", homeUI()),
        tabItem(tabName = "preprocessing", preprocessingUI()),
        tabItem(tabName = "ensembl_conversion", ensemblConversionUI()),
        tabItem(tabName = "analysis", analysisUI()),
        tabItem(tabName = "visualization", visualizationUI()),
        tabItem(tabName = "pathway_enrichment", pathwayEnrichmentUI()),
        tabItem(tabName = "go_enrichment", goEnrichmentUI()),
        tabItem(tabName = "help", helpUI())
      )
    )
  )
}