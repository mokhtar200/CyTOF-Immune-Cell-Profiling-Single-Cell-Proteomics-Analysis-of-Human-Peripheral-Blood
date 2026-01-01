# CyTOF-Immune-Cell-Profiling-Single-Cell-Proteomics-Analysis-of-Human-Peripheral-Blood

## Project Overview

This project performs a comprehensive single-cell proteomics analysis on human peripheral blood samples using CyTOF (mass cytometry) data from the Levine_32dim.fcs dataset.  
The goal is to identify immune cell populations, annotate clusters, and investigate differential marker expression patterns across clusters.

Key steps in the analysis:

1. Data Loading and Preprocessing  
   - Load FCS files.  
   - Arcsinh transformation of protein expression values.  
   - Remove non-informative channels (Time, event_number, DNA).  
   - Handle missing or infinite values.

2. Clustering
   - Perform clustering using FlowSOM.  
   - Reduce the number of clusters to avoid excessive fragmentation.  

3. SingleCellExperiment Object
   - Store expression matrix and cluster assignments in a `SingleCellExperiment` object.  

4. Dimensionality Reduction
   - Principal Component Analysis (PCA).  
   - UMAP visualization.  
   - Plot UMAP colored by clusters and annotated cell types.

5. Cluster Annotation
   - Annotate clusters into canonical immune cell types (T cells, B cells, NK cells, Monocytes, Myeloid).  
   - Reduce "Unknown" clusters using marker-based rules.

6. Differential Marker Analysis
   - Compute mean marker expression per cluster.  
   - Generate heatmaps, boxplots, and volcano plots for each marker.  


