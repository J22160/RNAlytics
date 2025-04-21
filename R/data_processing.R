# RNAlytics/R/data_processing.R
# Functions for data loading, preprocessing, and analysis

source("R/utilities.R")

# Preprocess raw data with duplicate row name handling
preprocessData <- function(input, rv) {
  req(input$rawFile)
  tryCatch({
    # Read the file and handle duplicates
    df <- read_tsv(input$rawFile$datapath, col_names = TRUE)
    id_col <- names(df)[1]  # Assume first column is the identifier
    
    # Check for duplicate row names
    if (any(duplicated(df[[id_col]]))) {
      df <- df %>%
        group_by(!!sym(id_col)) %>%
        summarise(across(everything(), sum, na.rm = TRUE)) %>%
        as.data.frame()
      showNotification("Duplicate row names detected in raw data. Aggregated by summing counts.", type = "warning")
    }
    
    # Convert to row names after aggregation
    df <- df %>% column_to_rownames(var = id_col)
    validateData(df)
    
    # Apply minimum count threshold
    df <- df[rowSums(df) >= input$minCount, ]
    
    # Normalize if requested
    if (input$normalize) {
      colData <- data.frame(sample = colnames(df), condition = "unknown")
      dds <- DESeqDataSetFromMatrix(countData = round(df), colData = colData, design = ~ 1)
      dds <- estimateSizeFactors(dds)
      df <- as.data.frame(counts(dds, normalized = TRUE))
    }
    
    rv$preprocessedData <- df
    rv$exprData <- df
    shinyjs::enable("downloadPreprocessed")
    showNotification("Data preprocessed successfully! Check the preview table.", type = "success")
  }, error = function(e) {
    rv$error <- paste("Preprocessing error:", e$message)
    rv$preprocessedData <- NULL
    shinyjs::disable("downloadPreprocessed")
    showNotification(paste("Preprocessing failed:", e$message), type = "error")
    print(paste("Error details:", e$message))
  })
}

# Optimized Ensembl ID conversion with duplicate handling
convertEnsemblIDs <- function(input, rv, session) {
  req(input$ensemblFile)
  withProgress(message = "Converting Ensembl IDs...", value = 0, {
    tryCatch({
      incProgress(0.1, detail = "Reading input file")
      df <- read_tsv(input$ensemblFile$datapath, col_names = TRUE)
      id_col <- names(df)[1]  # First column as Ensembl IDs
      
      # Aggregate duplicate Ensembl IDs
      if (any(duplicated(df[[id_col]]))) {
        df <- df %>%
          group_by(!!sym(id_col)) %>%
          summarise(across(everything(), sum, na.rm = TRUE)) %>%
          as.data.frame()
        showNotification("Duplicate Ensembl IDs detected. Aggregated by summing counts.", type = "warning")
      }
      
      df <- df %>% column_to_rownames(var = id_col)
      validateData(df)
      
      incProgress(0.3, detail = "Connecting to Ensembl database")
      ensembl <- useMart("ensembl", dataset = input$species)
      
      # Batch processing to reduce query load
      ensembl_ids <- unique(rownames(df))
      batch_size <- 10000
      batches <- split(ensembl_ids, ceiling(seq_along(ensembl_ids) / batch_size))
      gene_symbols <- data.frame()
      
      incProgress(0.5, detail = "Fetching gene symbols in batches")
      for (i in seq_along(batches)) {
        batch_ids <- batches[[i]]
        batch_symbols <- getBM(
          attributes = c("ensembl_gene_id", "external_gene_name"),
          filters = "ensembl_gene_id",
          values = batch_ids,
          mart = ensembl
        )
        gene_symbols <- rbind(gene_symbols, batch_symbols)
        incProgress(0.5 + 0.4 * (i / length(batches)), detail = paste("Processed batch", i, "of", length(batches)))
      }
      
      incProgress(0.9, detail = "Processing results")
      gene_symbols <- gene_symbols %>%
        group_by(ensembl_gene_id) %>%
        summarise(external_gene_name = first(na.omit(external_gene_name))) %>%
        filter(!is.na(ensembl_gene_id))
      
      converted <- df %>%
        rownames_to_column("ensembl_gene_id") %>%
        left_join(gene_symbols, by = "ensembl_gene_id") %>%
        mutate(Gene = ifelse(is.na(external_gene_name) | external_gene_name == "", ensembl_gene_id, external_gene_name)) %>%
        dplyr::select(-ensembl_gene_id, -external_gene_name) %>%
        group_by(Gene) %>%
        summarise(across(everything(), sum, na.rm = TRUE)) %>%
        column_to_rownames("Gene")
      
      incProgress(1.0, detail = "Complete")
      rv$convertedData <- converted
      rv$exprData <- converted
      shinyjs::enable("downloadConverted")
      showNotification("Ensembl IDs converted successfully!", type = "success")
    }, error = function(e) {
      rv$error <- paste("Conversion error:", e$message)
      rv$convertedData <- NULL
      shinyjs::disable("downloadConverted")
      showNotification(paste("Conversion failed:", e$message), type = "error")
    })
  })
}

# Load expression data with duplicate handling
loadExpressionData <- function(input, rv) {
  req(input$file)
  tryCatch({
    df <- read_tsv(input$file$datapath, col_names = TRUE)
    id_col <- names(df)[1]  # First column as gene IDs
    
    # Aggregate duplicate row names
    if (any(duplicated(df[[id_col]]))) {
      df <- df %>%
        group_by(!!sym(id_col)) %>%
        summarise(across(everything(), sum, na.rm = TRUE)) %>%
        as.data.frame()
      showNotification("Duplicate row names detected in expression data. Aggregated by summing counts.", type = "warning")
    }
    
    df <- df %>% column_to_rownames(var = id_col)
    validateData(df)
    rv$exprData <- df
    rv$error <- NULL
    shinyjs::enable("runAnalysis")
    showNotification("Data uploaded successfully!", type = "success")
  }, error = function(e) {
    rv$error <- paste("File upload error:", e$message)
    rv$exprData <- NULL
    shinyjs::disable("runAnalysis")
    shinyjs::disable("downloadFullResults")
    shinyjs::disable("downloadFilteredResults")
    showNotification(paste("Upload failed:", e$message), type = "error")
  })
}

# Load annotation data (unchanged, no row names involved)
loadAnnotationData <- function(input, rv) {
  req(input$annotationFile)
  tryCatch({
    annot <- read_tsv(input$annotationFile$datapath, col_names = TRUE)
    if (!"Gene" %in% colnames(annot)) stop("Annotation file must contain a 'Gene' column.")
    rv$annotationData <- annot
    showNotification("Annotations uploaded successfully!", type = "success")
  }, error = function(e) {
    rv$annotationData <- NULL
    showNotification(paste("Annotation upload failed:", e$message), type = "error")
  })
}

# Run differential expression analysis with enhanced variance filtering
runDifferentialAnalysis <- function(input, rv, session) {
  req(rv$exprData, rv$sampleAssignments, rv$comparison, input$statTest, input$logFCThresh, input$padjThresh)
  
  withProgress(message = "Running Analysis...", value = 0, {
    tryCatch({
      incProgress(0.2, detail = "Preparing data")
      
      colData <- data.frame(sample = colnames(rv$exprData), condition = NA_character_)
      for (group in names(rv$sampleAssignments)) {
        colData$condition[colData$sample %in% rv$sampleAssignments[[group]]] <- group
      }
      colData <- colData[complete.cases(colData), ]
      colData$condition <- factor(colData$condition)
      
      countData <- rv$exprData[, colData$sample]
      countData <- countData[rowSums(countData) > 1, ]
      
      # Early check for sufficient samples
      if (ncol(countData) < 2) {
        stop("Insufficient samples (need at least 2) for differential analysis.")
      }
      
      incProgress(0.5, detail = "Performing statistical test")
      
      if (input$statTest == "DESeq2") {
        countData <- round(countData)
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
        dds <- DESeq(dds, fitType = "local")
        res <- results(dds, contrast = c("condition", rv$comparison[2], rv$comparison[1])) %>%
          as.data.frame() %>%
          tibble::as_tibble(rownames = "Gene") %>%
          dplyr::select(Gene, log2FoldChange, padj) %>%
          dplyr::mutate(padj = ifelse(is.na(padj), 1, padj))
        norm_counts <- log2(counts(dds, normalized = TRUE) + 1)
        # Enhanced variance filtering
        variances <- apply(norm_counts, 2, var, na.rm = TRUE)
        valid_cols <- variances > 0
        if (!all(valid_cols)) {
          norm_counts <- norm_counts[, valid_cols, drop = FALSE]
          colData <- colData[valid_cols, , drop = FALSE]
          showNotification("Removed samples with zero variance from normalized counts.", type = "warning")
        }
        # Check if enough samples remain
        if (ncol(norm_counts) < 2) {
          stop("After filtering, insufficient variable samples (need at least 2) for analysis.")
        }
        rv$normCounts <- norm_counts
      } else if (input$statTest == "limma") {
        design <- model.matrix(~ 0 + condition, data = colData)
        colnames(design) <- levels(colData$condition)
        fit <- lmFit(log2(countData + 1), design)
        contrast.matrix <- makeContrasts(contrasts = paste(rv$comparison[2], "-", rv$comparison[1]), levels = design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH") %>%
          tibble::as_tibble(rownames = "Gene") %>%
          dplyr::select(Gene, logFC, adj.P.Val) %>%
          dplyr::rename(log2FoldChange = logFC, padj = adj.P.Val)
        norm_counts <- log2(countData + 1)
        # Enhanced variance filtering
        variances <- apply(norm_counts, 2, var, na.rm = TRUE)
        valid_cols <- variances > 0
        if (!all(valid_cols)) {
          norm_counts <- norm_counts[, valid_cols, drop = FALSE]
          colData <- colData[valid_cols, , drop = FALSE]
          showNotification("Removed samples with zero variance from normalized counts.", type = "warning")
        }
        # Check if enough samples remain
        if (ncol(norm_counts) < 2) {
          stop("After filtering, insufficient variable samples (need at least 2) for analysis.")
        }
        rv$normCounts <- norm_counts
      } else if (input$statTest == "t-test") {
        group1 <- colData$sample[colData$condition == rv$comparison[1]]
        group2 <- colData$sample[colData$condition == rv$comparison[2]]
        res <- apply(countData, 1, function(x) {
          t.test(x[group1], x[group2], var.equal = FALSE)
        }) %>%
          lapply(function(x) c(log2FC = log2(mean(x$estimate[2]) / mean(x$estimate[1])), pval = x$p.value)) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          tibble::as_tibble(rownames = "Gene") %>%
          dplyr::mutate(padj = p.adjust(pval, method = "BH")) %>%
          dplyr::select(Gene, log2FoldChange = log2FC, padj)
        norm_counts <- log2(countData + 1)
        # Enhanced variance filtering
        variances <- apply(norm_counts, 2, var, na.rm = TRUE)
        valid_cols <- variances > 0
        if (!all(valid_cols)) {
          norm_counts <- norm_counts[, valid_cols, drop = FALSE]
          colData <- colData[valid_cols, , drop = FALSE]
          showNotification("Removed samples with zero variance from normalized counts.", type = "warning")
        }
        # Check if enough samples remain
        if (ncol(norm_counts) < 2) {
          stop("After filtering, insufficient variable samples (need at least 2) for analysis.")
        }
        rv$normCounts <- norm_counts
      }
      
      if (!is.null(rv$annotationData)) {
        res <- res %>% left_join(rv$annotationData, by = "Gene")
      }
      
      rv$fullResults <- res %>% dplyr::arrange(padj)
      rv$filteredResults <- rv$fullResults %>% dplyr::filter(abs(log2FoldChange) >= input$logFCThresh & padj <= input$padjThresh)
      
      comp_name <- paste(rv$comparison[2], "vs", rv$comparison[1])
      rv$history[[comp_name]] <- list(fullResults = rv$fullResults, filteredResults = rv$filteredResults)
      updateSelectInput(session, "comparisonHistory", choices = c("None" = "", names(rv$history)))
      
      rv$error <- NULL
      shinyjs::enable("downloadFullResults")
      shinyjs::enable("downloadFilteredResults")
      shinyjs::enable("downloadPCAPlot")
      shinyjs::enable("downloadVolcanoPlot")
      shinyjs::enable("downloadHeatmapPlot")
      
      incProgress(1.0, detail = "Complete")
      showNotification("Analysis complete! Check the Visualization tab.", type = "success")
    }, error = function(e) {
      rv$error <- paste("Analysis error:", e$message)
      rv$fullResults <- NULL
      rv$filteredResults <- NULL
      rv$normCounts <- NULL
      shinyjs::disable("downloadFullResults")
      shinyjs::disable("downloadFilteredResults")
      shinyjs::disable("downloadPCAPlot")
      shinyjs::disable("downloadVolcanoPlot")
      shinyjs::disable("downloadHeatmapPlot")
      showNotification(paste("Analysis failed:", e$message), type = "error")
    })
  })
}

# Run pathway enrichment analysis (KEGG) with fixed syntax
runPathwayEnrichment <- function(input, rv, session) {
  req(rv$filteredResults)
  withProgress(message = "Running KEGG Enrichment...", value = 0, {
    tryCatch({
      incProgress(0.2, detail = "Preparing gene list")
      degs <- rv$filteredResults$Gene
      annotation_db <- switch(input$enrichSpecies,
                              "hsa" = org.Hs.eg.db,
                              "mmu" = org.Mm.eg.db)
      entrez_ids <- mapIds(annotation_db, keys = degs, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
      entrez_ids <- na.omit(entrez_ids)
      
      incProgress(0.5, detail = "Performing KEGG enrichment")
      enrich_result <- enrichKEGG(
        gene = entrez_ids,
        organism = input$enrichSpecies,
        pvalueCutoff = input$enrichPval,
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500
      )
      
      incProgress(1.0, detail = "Complete")
      rv$enrichResultObj <- enrich_result
      rv$enrichResults <- as.data.frame(enrich_result)
      shinyjs::enable("downloadEnrichResults")
      shinyjs::enable("downloadEnrichDotPlot")
      showNotification("KEGG enrichment completed successfully!", type = "success")
    }, error = function(e) {
      rv$enrichResultObj <- NULL
      rv$enrichResults <- NULL
      shinyjs::disable("downloadEnrichResults")
      shinyjs::disable("downloadEnrichDotPlot")
      showNotification(paste("KEGG enrichment failed:", e$message), type = "error")
    })
  })
}

# Run Gene Ontology enrichment analysis
runGOEnrichment <- function(input, rv, session) {
  req(rv$filteredResults)
  withProgress(message = "Running GO Enrichment...", value = 0, {
    tryCatch({
      incProgress(0.2, detail = "Preparing gene list")
      degs <- rv$filteredResults$Gene
      annotation_db <- switch(input$goSpecies,
                              "hsa" = org.Hs.eg.db,
                              "mmu" = org.Mm.eg.db)
      entrez_ids <- mapIds(annotation_db, keys = degs, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
      entrez_ids <- na.omit(entrez_ids)
      
      incProgress(0.5, detail = "Performing GO enrichment")
      go_result <- enrichGO(
        gene = entrez_ids,
        OrgDb = annotation_db,
        ont = input$goOntology,
        pvalueCutoff = input$goPval,
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      )
      
      incProgress(1.0, detail = "Complete")
      rv$goResultObj <- go_result
      rv$goResults <- as.data.frame(go_result)
      shinyjs::enable("downloadGOResults")
      shinyjs::enable("downloadGODotPlot")
      showNotification("GO enrichment completed successfully!", type = "success")
    }, error = function(e) {
      rv$goResultObj <- NULL
      rv$goResults <- NULL
      shinyjs::disable("downloadGOResults")
      shinyjs::disable("downloadGODotPlot")
      showNotification(paste("GO enrichment failed:", e$message), type = "error")
    })
  })
}