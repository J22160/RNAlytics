# RNAlytics/R/server_logic.R
# Main server-side logic


source("R/data_processing.R")
source("R/utilities.R")
source("R/plotting.R")

appServer <- function(input, output, session) {
  # Reactive values to store app state
  rv <- reactiveValues(
    rawData = NULL,
    preprocessedData = NULL,
    convertedData = NULL,
    exprData = NULL,
    annotationData = NULL,
    groups = NULL,
    sampleAssignments = NULL,
    comparison = NULL,
    fullResults = NULL,
    filteredResults = NULL,
    normCounts = NULL,
    enrichResultObj = NULL,
    enrichResults = NULL,
    goResultObj = NULL,
    goResults = NULL,
    error = NULL,
    history = list()
  )
    
  # Preprocessing tab logic
  observeEvent(input$preprocessData, {
    preprocessData(input, rv)
  })
  
  output$preprocessedTable <- renderDT({
    req(rv$preprocessedData)
    datatable(
      rv$preprocessedData,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = TRUE
    )
  })
  
  output$downloadPreprocessed <- downloadHandler(
    filename = "Preprocessed_Data.csv",
    content = function(file) { write.csv(rv$preprocessedData, file, row.names = TRUE) }
  )
  
  # Ensembl Conversion tab logic
  observeEvent(input$convertEnsembl, {
    convertEnsemblIDs(input, rv, session)
  })
  
  output$convertedTable <- renderDT({
    req(rv$convertedData)
    datatable(rv$convertedData, options = list(pageLength = 10, scrollX= TRUE), rownames = TRUE)
  })
  
  output$downloadConverted <- downloadHandler(
    filename = "Converted_Data.csv",
    content = function(file) { write.csv(rv$convertedData, file, row.names = TRUE) }
  )
  
  # Analysis tab logic
  observeEvent(input$file, {
    loadExpressionData(input, rv)
  })
  
  observeEvent(input$annotationFile, {
    loadAnnotationData(input, rv)
  })
  
  output$groupDefinitionsUI <- renderUI({
    req(rv$exprData)
    tagList(
      numericInput("numGroups", "Number of Groups/Conditions", value = 2, min = 2, max = ncol(rv$exprData)),
      uiOutput("groupNamesUI")
    )
  })
  
  output$groupNamesUI <- renderUI({
    req(input$numGroups)
    lapply(1:input$numGroups, function(i) {
      textInput(paste0("groupName_", i), paste("Group", i, "Name"), value = paste0("Group", i))
    })
  })
  
  output$sampleAssignmentsUI <- renderUI({
    req(rv$exprData, input$numGroups)
    sampleCols <- colnames(rv$exprData)
    groupNames <- sapply(1:input$numGroups, function(i) input[[paste0("groupName_", i)]], USE.NAMES = FALSE)
    if (any(is.null(groupNames) | groupNames == "")) return(NULL)
    rv$groups <- groupNames
    lapply(seq_along(groupNames), function(i) {
      selectizeInput(paste0("samplesGroup_", i), paste("Samples for", groupNames[i]), choices = sampleCols, multiple = TRUE)
    })
  })
  
  output$comparisonSelectionUI <- renderUI({
    req(rv$groups)
    groupNames <- rv$groups
    tagList(
      selectInput("comparisonBase", "Base Condition", choices = groupNames),
      selectInput("comparisonContrast", "Contrast Condition", choices = groupNames)
    )
  })
  
  observe({
    req(rv$exprData, input$numGroups)
    groupNames <- sapply(1:input$numGroups, function(i) input[[paste0("groupName_", i)]], USE.NAMES = FALSE)
    if (any(is.null(groupNames) | groupNames == "")) return()
    assignments <- lapply(1:input$numGroups, function(i) input[[paste0("samplesGroup_", i)]])
    names(assignments) <- groupNames
    rv$sampleAssignments <- assignments
    req(input$comparisonBase, input$comparisonContrast)
    rv$comparison <- c(input$comparisonBase, input$comparisonContrast)
  })
  
  observeEvent(input$runAnalysis, {
    runDifferentialAnalysis(input, rv, session)
  })
  
  observeEvent(input$comparisonHistory, {
    req(input$comparisonHistory != "")
    rv$fullResults <- rv$history[[input$comparisonHistory]]$fullResults
    rv$filteredResults <- rv$history[[input$comparisonHistory]]$filteredResults
    shinyjs::enable("downloadFullResults")
    shinyjs::enable("downloadFilteredResults")
    shinyjs::enable("downloadPCAPlot")
    shinyjs::enable("downloadVolcanoPlot")
    shinyjs::enable("downloadHeatmapPlot")
    shinyjs::enable("downloadGeneCountPlot")
  })
  
  output$resultsTable <- renderDT({
    req(rv$fullResults)
    filtered_data <- rv$fullResults
    if (input$filterType == "gte") {
      filtered_data <- filtered_data %>% dplyr::filter(log2FoldChange >= input$filterValue)
    } else if (input$filterType == "lte") {
      filtered_data <- filtered_data %>% dplyr::filter(log2FoldChange <= input$filterValue)
    }
    datatable(filtered_data, filter = "top", options = list(pageLength = 10, autoWidth = TRUE, searchHighlight = TRUE, dom = 'lfrtip'), rownames = FALSE) %>%
      formatRound("log2FoldChange", digits = 4) %>%
      formatSignif("padj", digits = 5)
  })
  
  output$errorLog <- renderPrint({
    if (!is.null(rv$error)) rv$error else "No errors detected."
  })
  
  output$downloadFullResults <- downloadHandler(
    filename = function() { paste("Full_DE_Results_", input$statTest, "_", Sys.Date(), ".csv", sep = "") },
    content = function(file) { write.csv(rv$fullResults, file, row.names = FALSE) }
  )
  
  output$downloadFilteredResults <- downloadHandler(
    filename = function() { paste("Filtered_DE_Results_", input$statTest, "_", Sys.Date(), ".csv", sep = "") },
    content = function(file) { write.csv(rv$filteredResults, file, row.names = FALSE) }
  )
  
  # Visualization tab logic
  observe({
    req(rv$normCounts)
    updateSelectizeInput(session, "geneCountGene", choices = rownames(rv$normCounts))
  })

  observe({
    req(rv$sampleAssignments)
    n_groups <- length(rv$sampleAssignments)
    updatePickerInput(
      session,
      "pcaColors",
      selected = head(c("Red", "Blue", "Green", "Purple", "Orange", "Yellow"), n_groups),
      choices = c("Red", "Blue", "Green", "Purple", "Orange", "Yellow", "Cyan", "Magenta"),
      options = list(`min-options` = n_groups)
  )
  })

  observe({
  req(rv$exprData, input$numGroups)
  groupNames <- sapply(1:input$numGroups, function(i) input[[paste0("groupName_", i)]], USE.NAMES = FALSE)
  if (any(is.null(groupNames) | groupNames == "")) return()
  assignments <- lapply(1:input$numGroups, function(i) input[[paste0("samplesGroup_", i)]])
  names(assignments) <- groupNames
  rv$sampleAssignments <- assignments
  req(input$comparisonBase, input$comparisonContrast)
  rv$comparison <- c(input$comparisonBase, input$comparisonContrast)
  })

observe({
  req(rv$sampleAssignments)
  n_groups <- length(rv$sampleAssignments)
  updatePickerInput(
    session,
    "pcaColors",
    selected = head(c("Red", "Blue", "Green", "Purple", "Orange", "Yellow"), n_groups),
    choices = c("Red", "Blue", "Green", "Purple", "Orange", "Yellow", "Cyan", "Magenta"),
    options = list(`min-options` = n_groups)
  )
  updatePickerInput(
    session,
    "geneCountColors",
    selected = head(c("Red", "Blue", "Green", "Purple", "Orange", "Yellow"), n_groups),
    choices = c("Red", "Blue", "Green", "Purple", "Orange", "Yellow", "Cyan", "Magenta"),
    options = list(`min-options` = n_groups)
  )
  })

  output$pcaPlot <- renderPlot({
    req(rv$normCounts, rv$sampleAssignments)
    generatePCAPlot(input, rv)
  })

  output$heatmapPlot <- renderPlot({
    req(rv$normCounts, rv$fullResults, rv$sampleAssignments)
    generateHeatmapPlot(input, rv)
  })
  
  output$volcanoPlot <- renderPlot({
    req(rv$fullResults)
    generateVolcanoPlot(input, rv)
  })
  
  output$geneCountPlot <- renderPlot({
    req(rv$normCounts, rv$sampleAssignments, input$geneCountGene)
    generateGeneCountPlot(input, rv)
  })

  output$downloadPCAPlot <- downloadHandler(
    filename = function() { paste("PCA_Plot_", Sys.Date(), ".pdf", sep = "") },
    content = function(file) { savePCAPlot(input, rv, file) }
  )
  
  output$downloadVolcanoPlot <- downloadHandler(
    filename = function() { paste("Volcano_Plot_", Sys.Date(), ".pdf", sep = "") },
    content = function(file) { saveVolcanoPlot(input, rv, file) }
  )
  
  output$downloadHeatmapPlot <- downloadHandler(
    filename = function() { paste("Heatmap_", Sys.Date(), ".pdf", sep = "") },
    content = function(file) { saveHeatmapPlot(input, rv, file) }
  )
  
  output$downloadGeneCountPlot <- downloadHandler(
    filename = function() { paste("GeneCount_", input$geneCountGene, "_", Sys.Date(), ".pdf", sep = "") },
    content = function(file) { saveGeneCountPlot(input, rv, file) }
  )

  # KEGG Pathway Enrichment tab logic
  observeEvent(input$runEnrichment, {
    runPathwayEnrichment(input, rv, session)
  })
  
  output$enrichTable <- renderDT({
    req(rv$enrichResults)
    datatable(
      rv$enrichResults,
      options = list(pageLength = 10, scrollX= TRUE),
      rownames = FALSE
    ) %>%
      formatSignif(c("pvalue", "p.adjust", "qvalue"), digits = 3)
  })
  
  output$enrichDotPlot <- renderPlot({
    req(rv$enrichResultObj)
    species_name <- ifelse(input$enrichSpecies == "hsa", "Human", "Mouse")
    dotplot(rv$enrichResultObj, showCategory = 10) +
      ggtitle(paste("KEGG Pathway Enrichment -", species_name)) +
      theme_minimal()
  })
  
  output$downloadEnrichResults <- downloadHandler(
    filename = function() { paste("KEGG_Enrichment_Results_", input$enrichSpecies, "_", Sys.Date(), ".csv", sep = "") },
    content = function(file) { write.csv(rv$enrichResults, file, row.names = FALSE) }
  )
  
  output$downloadEnrichDotPlot <- downloadHandler(
    filename = function() { paste("KEGG_DotPlot_", input$enrichSpecies, "_", Sys.Date(), ".pdf", sep = "") },
    content = function(file) {
      species_name <- ifelse(input$enrichSpecies == "hsa", "Human", "Mouse")
      p <- dotplot(rv$enrichResultObj, showCategory = 10) +
        ggtitle(paste("KEGG Pathway Enrichment -", species_name)) +
        theme_minimal()
      ggsave(file, plot = p, width = 8, height = 5, dpi = 800)
    }
  )
  
  # GO Enrichment tab logic
  observeEvent(input$runGOEnrichment, {
    runGOEnrichment(input, rv, session)
  })
  
  output$goTable <- renderDT({
    req(rv$goResults)
    datatable(
      rv$goResults,
      options = list(pageLength = 10, autoWidth = TRUE, scrollX= TRUE),
      rownames = FALSE
    ) %>%
      formatSignif(c("pvalue", "p.adjust", "qvalue"), digits = 3)
  })
  
  output$goDotPlot <- renderPlot({
    req(rv$goResultObj)
    species_name <- ifelse(input$goSpecies == "hsa", "Human", "Mouse")
    ontology_name <- switch(input$goOntology,
                            "BP" = "Biological Process",
                            "MF" = "Molecular Function",
                            "CC" = "Cellular Component")
    dotplot(rv$goResultObj, showCategory = 10) +
      ggtitle(paste("GO Enrichment (", ontology_name, ") -", species_name)) +
      theme_minimal()
  })
  
  output$downloadGOResults <- downloadHandler(
    filename = function() { paste("GO_Enrichment_Results_", input$goSpecies, "_", input$goOntology, "_", Sys.Date(), ".csv", sep = "") },
    content = function(file) { write.csv(rv$goResults, file, row.names = FALSE) }
  )
  
  output$downloadGODotPlot <- downloadHandler(
    filename = function() { paste("GO_DotPlot_", input$goSpecies, "_", input$goOntology, "_", Sys.Date(), ".pdf", sep = "") },
    content = function(file) {
      species_name <- ifelse(input$goSpecies == "hsa", "Human", "Mouse")
      ontology_name <- switch(input$goOntology,
                              "BP" = "Biological Process",
                              "MF" = "Molecular Function",
                              "CC" = "Cellular Component")
      p <- dotplot(rv$goResultObj, showCategory = 10) +
        ggtitle(paste("GO Enrichment (", ontology_name, ") -", species_name)) +
        theme_minimal()
      ggsave(file, plot = p, width = 8, height = 5, dpi = 800)
    }
  )
}