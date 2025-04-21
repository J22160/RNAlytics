# RNAlytics/R/plotting.R
# Plotting functions for visualizations

library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)
library(RColorBrewer)

# Generate PCA plot (static with ggplot2)
generatePCAPlot <- function(input, rv) {
  req(rv$normCounts, rv$sampleAssignments)
  
  pca_res <- prcomp(t(rv$normCounts), scale. = TRUE)
  pca_data <- as.data.frame(pca_res$x)
  pca_data$sample <- rownames(pca_data)
  pca_data$condition <- NA_character_
  
  for (group in names(rv$sampleAssignments)) {
    pca_data$condition[pca_data$sample %in% rv$sampleAssignments[[group]]] <- group
  }
  
  percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)[1:2]
  
  n_conditions <- length(unique(pca_data$condition))
  n_colors <- length(input$pcaColors)
  if (n_colors < n_conditions) {
    extra_colors <- brewer.pal(max(3, n_conditions - n_colors), "Set1")
    pca_colors <- c(input$pcaColors, extra_colors[1:(n_conditions - n_colors)])
  } else {
    pca_colors <- input$pcaColors[1:n_conditions]
  }
  
  ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample)) +
    geom_point(size = input$pcaPointSize) +
    geom_text(vjust = -1, size = input$pcaLabelSize) +
    labs(
      title = input$pcaTitle,
      x = paste0(input$pcaXLabel, " (", percentVar[1], "% variance)"),
      y = paste0(input$pcaYLabel, " (", percentVar[2], "% variance)")
    ) +
    scale_color_manual(values = pca_colors) +
    theme_minimal(base_size = input$pcaFontSize) +
    theme(legend.title = element_blank())
}

# Save PCA plot
savePCAPlot <- function(input, rv, file) {
  req(rv$normCounts, rv$sampleAssignments)
  
  pca_res <- prcomp(t(rv$normCounts), scale. = TRUE)
  pca_data <- as.data.frame(pca_res$x)
  pca_data$sample <- rownames(pca_data)
  pca_data$condition <- NA_character_
  
  for (group in names(rv$sampleAssignments)) {
    pca_data$condition[pca_data$sample %in% rv$sampleAssignments[[group]]] <- group
  }
  
  percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)[1:2]
  
  n_conditions <- length(unique(pca_data$condition))
  n_colors <- length(input$pcaColors)
  if (n_colors < n_conditions) {
    extra_colors <- brewer.pal(max(3, n_conditions - n_colors), "Set1")
    pca_colors <- c(input$pcaColors, extra_colors[1:(n_conditions - n_colors)])
  } else {
    pca_colors <- input$pcaColors[1:n_conditions]
  }
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample)) +
    geom_point(size = input$pcaPointSize) +
    geom_text(vjust = -1, size = input$pcaLabelSize) +
    labs(
      title = input$pcaTitle,
      x = paste0(input$pcaXLabel, " (", percentVar[1], "% variance)"),
      y = paste0(input$pcaYLabel, " (", percentVar[2], "% variance)")
    ) +
    scale_color_manual(values = pca_colors) +
    theme_minimal(base_size = input$pcaFontSize) +
    theme(legend.title = element_blank())
  
  ggsave(file, plot = p, width = 16, height = 7, dpi = 1200)
}

# Generate Enhanced Volcano plot
generateVolcanoPlot <- function(input, rv) {
  req(rv$fullResults)
  data <- rv$fullResults
  
  top_n <- as.numeric(input$volcanoTopGenes)
  if (top_n > 0) {
    top_genes <- data %>%
      arrange(padj) %>%
      head(top_n) %>%
      pull(Gene)
    data$label <- ifelse(data$Gene %in% top_genes, data$Gene, NA)
  } else {
    data$label <- NA
  }
  
  EnhancedVolcano(
    toptable = data,
    lab = data$label,
    x = "log2FoldChange",
    y = "padj",
    title = input$volcanoTitle,
    subtitle = input$volcanoSubtitle,
    caption = paste("Thresholds: log2FC =", input$logFCThresh, ", padj =", input$padjThresh),
    xlab = input$volcanoXLabel,
    ylab = input$volcanoYLabel,
    pCutoff = input$padjThresh,
    FCcutoff = input$logFCThresh,
    pointSize = input$volcanoPointSize,
    labSize = input$volcanoLabelSize,
    col = input$volcanoColors[1:4],
    colAlpha = 0.6,
    legendPosition = input$volcanoLegendPos,
    legendLabels = c("Not Sig.", "Log2 FC", "P-adj", "Sig. + LFC"),
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    maxoverlaps = top_n,
    max.overlaps = top_n
  ) +
    theme_minimal(base_size = input$volcanoFontSize)
}

# Save Enhanced Volcano plot
saveVolcanoPlot <- function(input, rv, file) {
  p <- generateVolcanoPlot(input, rv)
  ggsave(file, plot = p, width = 10, height = 10, dpi = 1200)
}

# Generate Heatmap plot
generateHeatmapPlot <- function(input, rv) {
  req(rv$normCounts, rv$fullResults, rv$sampleAssignments)
  
  top_genes <- rv$fullResults %>% 
    arrange(padj) %>% 
    head(50) %>% 
    pull(Gene)
  heatmap_data <- rv$normCounts[top_genes, , drop = FALSE]
  
  annotation_col <- data.frame(
    Sample = colnames(heatmap_data),
    Condition = NA_character_
  )
  for (group in names(rv$sampleAssignments)) {
    annotation_col$Condition[annotation_col$Sample %in% rv$sampleAssignments[[group]]] <- group
  }
  rownames(annotation_col) <- annotation_col$Sample
  annotation_col <- annotation_col[, "Condition", drop = FALSE]
  
  condition_levels <- unique(annotation_col$Condition)
  ann_colors <- list(Condition = setNames(input$pcaColors[1:length(condition_levels)], condition_levels))
  
  pheatmap(
    heatmap_data,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = input$heatmapClusterMethod,
    clustering_distance_cols = input$heatmapClusterMethod,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = input$heatmapShowRownames,
    show_colnames = TRUE,
    treeheight_row = 50,
    treeheight_col = 50,
    main = input$heatmapTitle,
    color = colorRampPalette(rev(brewer.pal(11, input$heatmapPalette)))(100),
    fontsize = input$heatmapFontSize
  )
}

# Save Heatmap plot
saveHeatmapPlot <- function(input, rv, file) {
  req(rv$normCounts, rv$fullResults, rv$sampleAssignments)
  
  top_genes <- rv$fullResults %>% 
    arrange(padj) %>% 
    head(50) %>% 
    pull(Gene)
  heatmap_data <- rv$normCounts[top_genes, , drop = FALSE]
  
  annotation_col <- data.frame(
    Sample = colnames(heatmap_data),
    Condition = NA_character_
  )
  for (group in names(rv$sampleAssignments)) {
    annotation_col$Condition[annotation_col$Sample %in% rv$sampleAssignments[[group]]] <- group
  }
  rownames(annotation_col) <- annotation_col$Sample
  annotation_col <- annotation_col[, "Condition", drop = FALSE]
  
  condition_levels <- unique(annotation_col$Condition)
  ann_colors <- list(Condition = setNames(input$pcaColors[1:length(condition_levels)], condition_levels))
  
  pheatmap(
    heatmap_data,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = input$heatmapClusterMethod,
    clustering_distance_cols = input$heatmapClusterMethod,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = input$heatmapShowRownames,
    show_colnames = TRUE,
    treeheight_row = 50,
    treeheight_col = 50,
    main = input$heatmapTitle,
    color = colorRampPalette(rev(brewer.pal(11, input$heatmapPalette)))(100),
    fontsize = input$heatmapFontSize,
    filename = file,
    width = 15,
    height = 18,
    dpi = 1200
  )
}

# Generate Gene Count Plot (updated with color validation)
generateGeneCountPlot <- function(input, rv) {
  req(rv$normCounts, rv$sampleAssignments, input$geneCountGene)
  
  gene_data <- data.frame(
    Sample = colnames(rv$normCounts),
    Expression = rv$normCounts[input$geneCountGene, ]
  )
  
  gene_data$Condition <- NA_character_
  for (group in names(rv$sampleAssignments)) {
    gene_data$Condition[gene_data$Sample %in% rv$sampleAssignments[[group]]] <- group
  }
  gene_data$Condition <- factor(gene_data$Condition, levels = names(rv$sampleAssignments))
  
  # Ensure enough colors for all conditions
  n_conditions <- length(unique(gene_data$Condition))
  n_colors <- length(input$geneCountColors)
  if (n_colors < n_conditions) {
    extra_colors <- brewer.pal(max(3, n_conditions - n_colors), "Set1")
    gene_colors <- c(input$geneCountColors, extra_colors[1:(n_conditions - n_colors)])
  } else {
    gene_colors <- input$geneCountColors[1:n_conditions]
  }
  
  ggplot(gene_data, aes(x = Sample, y = Expression, color = Condition)) +
    geom_point(size = input$geneCountPointSize) +
    labs(
      title = paste(input$geneCountTitle, "-", input$geneCountGene),
      x = "Samples",
      y = "Normalized Count (log2)"
    ) +
    scale_color_manual(values = gene_colors) +
    theme_minimal(base_size = input$geneCountFontSize) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Save Gene Count Plot (updated with color validation)
saveGeneCountPlot <- function(input, rv, file) {
  req(rv$normCounts, rv$sampleAssignments, input$geneCountGene)
  
  gene_data <- data.frame(
    Sample = colnames(rv$normCounts),
    Expression = rv$normCounts[input$geneCountGene, ]
  )
  
  gene_data$Condition <- NA_character_
  for (group in names(rv$sampleAssignments)) {
    gene_data$Condition[gene_data$Sample %in% rv$sampleAssignments[[group]]] <- group
  }
  gene_data$Condition <- factor(gene_data$Condition, levels = names(rv$sampleAssignments))
  
  n_conditions <- length(unique(gene_data$Condition))
  n_colors <- length(input$geneCountColors)
  if (n_colors < n_conditions) {
    extra_colors <- brewer.pal(max(3, n_conditions - n_colors), "Set1")
    gene_colors <- c(input$geneCountColors, extra_colors[1:(n_conditions - n_colors)])
  } else {
    gene_colors <- input$geneCountColors[1:n_conditions]
  }
  
  p <- ggplot(gene_data, aes(x = Sample, y = Expression, color = Condition)) +
    geom_point(size = input$geneCountPointSize) +
    labs(
      title = paste(input$geneCountTitle, "-", input$geneCountGene),
      x = "Samples",
      y = "Normalized Count (log2)"
    ) +
    scale_color_manual(values = gene_colors) +
    theme_minimal(base_size = input$geneCountFontSize) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file, plot = p, width = 12, height = 6, dpi = 1200)
}