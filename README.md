
# RNAlytics: A Shiny-Based RNA-Seq Analysis Web Application

**Version**: 1.0.0  
**Last Updated**: April 10, 2025  
**Author**: Jash Trivedi  
**License**: MIT

`RNAlytics` is a powerful, interactive web application built with R and Shiny, designed for comprehensive RNA-Seq data analysis and visualization. It enables researchers, bioinformaticians, and data scientists to perform differential expression analysis, conduct pathway and functional enrichment, and generate publication-quality visualizations—all within an intuitive interface. With a modular architecture leveraging industry-standard R packages, `RNAlytics` balances usability and advanced functionality, reflecting decades of software engineering best practices.

[Run RNAlytics Online](https://rnalytics.shinyapps.io/RNAlytics_app/) 

---

## Table of Contents
1. [Overview](#overview)
2. [Features](#features)
   - [Differential Expression Analysis](#differential-expression-analysis)
   - [PCA Plot](#pca-plot)
   - [Volcano Plot](#volcano-plot)
   - [Heatmap](#heatmap)
   - [Gene Count Plot](#gene-count-plot)
   - [KEGG Pathway Enrichment](#kegg-pathway-enrichment)
   - [Gene Ontology (GO) Analysis](#gene-ontology-go-analysis)
3. [Use Cases](#use-cases)
4. [Installation](#installation)
   - [Prerequisites](#prerequisites)
   - [Local Deployment with Shiny Server](#local-deployment-with-shiny-server)
   - [Troubleshooting Deployment](#troubleshooting-deployment)
5. [Usage](#usage)
6. [Support](#support)
7. [Acknowledgments](#acknowledgments)

---

## Overview

`RNAlytics` streamlines RNA-Seq workflows by integrating data processing, statistical analysis, functional enrichment, and visualization into a single platform. Built with R 4.4.3 and Shiny, it leverages packages like `DESeq2`, `ggplot2`, `pheatmap`, `EnhancedVolcano` and `clusterProfiler` to deliver high-performance analytics. The app’s modular design—split across `data_processing.R`, `plotting.R`, `server_logic.R`, `ui_components.R`, and `utilities.R`—ensures maintainability and extensibility.

---

## Features

### Differential Expression Analysis
- **Description**: Performs differential expression analysis using `DESeq2` or `Limma` on user-uploaded RNA-Seq count data. Users define sample groups and comparisons for downstream analysis and visualization.
- **Inputs**: Count matrix (.txt), group definitions, comparison base/contrast.
- **Outputs**: Table with log2 fold changes and adjusted p-values; downloadable as CSV.

### PCA Plot
- **Description**: Visualizes sample variance via a static PCA plot based on normalized counts. Points are colored by condition, with customizable labels and aesthetics.
- **Inputs**: Title, axis labels, point/label/font sizes, condition colors.
- **Outputs**: Scatter plot with variance percentages, downloadable as PDF (16x7 inches, 1200 DPI).

### Volcano Plot
- **Description**: Highlights differentially expressed genes (DEGs) with a volcano plot, based on log2 fold change and adjusted p-value. Top genes are labeled dynamically.
- **Inputs**: Title, subtitle, axis labels, thresholds, point/label sizes, colors, top gene count.

### Heatmap
- **Description**: Displays a clustered heatmap of the top 50 DEGs, with Z-score scaling and dendrograms. Samples are annotated by condition.
- **Inputs**: Title, clustering method, row name visibility, font size, color palette.
- **Outputs**: Static heatmap, downloadable as PDF (15x18 inches, 1200 DPI).

### Gene Count Plot
- **Description**: Plots normalized expression (log2 counts) of a user-specified gene across samples, colored by condition, for detailed gene-level analysis.
- **Inputs**: Gene name (dropdown), title, point/font sizes, condition colors.
- **Outputs**: Scatter plot with rotated x-axis labels, downloadable as PDF (12x6 inches, 1200 DPI).


### KEGG Pathway Enrichment
- **Description**: Performs Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway enrichment analysis on DEGs, identifying overrepresented biological pathways. Results are stored and visualized as tables or plots.
- **Inputs**: DEG list from differential analysis, organism type (e.g. human or mouse), p-value cutoff.
- **Outputs**: Table of enriched pathways with p-values, gene counts, and pathway IDs; downloadable as CSV. Optional bar or dot plot visualization.


### Gene Ontology (GO) Analysis
- **Description**: Conducts Gene Ontology enrichment analysis on DEGs, categorizing genes into Biological Process (BP), Molecular Function (MF), and Cellular Component (CC) terms. Results enhance functional interpretation.
- **Inputs**: DEG list, ontology type (BP/MF/CC), organism type, p-value cutoff.
- **Outputs**: Table of enriched GO terms with p-values, gene ratios, and term descriptions; downloadable as CSV. Optional visualization (e.g., bar plot or network).

---

## Use Cases

1. **Exploratory Analysis**:
   - **Scenario**: Assess sample clustering and gene expression variability.
   - **Features**: PCA Plot, Gene Count Plot.
   - **Outcome**: Detects outliers and validates experimental design.

2. **Differential Expression Studies**:
   - **Scenario**: Compare gene expression across conditions.
   - **Features**: Differential Expression Analysis, Volcano Plot, Heatmap.
   - **Outcome**: Identifies and visualizes significant DEGs.

3. **Pathway and Functional Insights**:
   - **Scenario**: Investigate biological implications of DEGs in a disease model.
   - **Features**: KEGG Pathway Enrichment, GO Analysis.
   - **Outcome**: Reveals enriched pathways (e.g., metabolism) and GO terms (e.g., immune response), aiding hypothesis generation.

4. **Candidate Gene Validation**:
   - **Scenario**: Examine expression of a specific gene across samples.
   - **Features**: Gene Count Plot.
   - **Outcome**: Confirms gene behavior for targeted studies.

5. **Comprehensive Genomic Workflow**:
   - **Scenario**: Analyze RNA-Seq data end-to-end for a publication.
   - **Features**: All features combined.
   - **Outcome**: Produces statistical results, visualizations, and functional annotations in one pipeline.

---


## Support
- **Email**: jashtrivedi221@gmail.com
- **Logs**: Share `/var/log/shiny-server.log` for deployment support.

---

## Acknowledgments
This web application was developed as part of an effort to streamline RNA-Seq data analysis for researchers and students. I would like to thank the open-source R and Bioconductor communities for their powerful tools, and the developers of packages like DESeq2, clusterProfiler, and EnhancedVolcano for making high-quality bioinformatics accessible. Special thanks to my mentors and peers for their feedback during development and testing.