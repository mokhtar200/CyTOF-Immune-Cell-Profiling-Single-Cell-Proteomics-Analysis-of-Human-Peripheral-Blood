# -------------------------------
# 0. Load required packages
# -------------------------------
library(flowCore)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(FlowSOM)
library(ggrepel)

# -------------------------------
# 1. Load FCS file
# -------------------------------
fcs_file <- "Levine_32dim.fcs"
fcs <- read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)

# -------------------------------
# 2. Extract and clean expression matrix
# -------------------------------
expr <- exprs(fcs)

# arcsinh transform (CyTOF standard)
expr <- asinh(expr / 5)

# Replace Inf / NaN values with 0
expr[!is.finite(expr)] <- 0

# Remove non-protein channels (Time, event_number, DNA)
bad_channels <- grep("Time|event|DNA", colnames(expr), ignore.case = TRUE)
if(length(bad_channels) > 0) expr <- expr[, -bad_channels]

expr_all <- expr  # keep for FlowSOM

cat("Expression matrix dimensions: ", dim(expr_all), "\n")

# -------------------------------
# 3. Clustering using FlowSOM
# -------------------------------
set.seed(123)
fsom <- FlowSOM::FlowSOM(
  as.matrix(expr_all),
  compensate = FALSE,
  transform = FALSE,
  scale = FALSE,
  colsToUse = 1:ncol(expr_all),
  nClus = 30   # reduced clusters to avoid fragmentation
)

cluster_id <- fsom$metaclustering  # 1:30
length(cluster_id)  # should match number of cells

# -------------------------------
# 4. Create SingleCellExperiment object
# -------------------------------
sce <- SingleCellExperiment(
  assays = list(exprs = t(expr_all))
)
sce$cluster_id <- factor(cluster_id)  # assign FlowSOM clusters

# -------------------------------
# 5. Dimensionality reduction
# -------------------------------
# PCA
sce <- runPCA(sce, exprs_values = "exprs")

# UMAP
sce <- runUMAP(sce, dimred = "PCA")

# Plot UMAP colored by cluster
plotUMAP(sce, colour_by = "cluster_id")

# -------------------------------
# 6. Cluster annotation (minimize Unknown)
# -------------------------------
# Compute cluster mean expression
expr_df <- as.data.frame(expr_all)
expr_df$cluster <- sce$cluster_id

cluster_means <- expr_df %>%
  group_by(cluster) %>%
  summarise(across(.cols = where(is.numeric), .fns = mean))

cluster_matrix <- as.matrix(cluster_means[,-1])
rownames(cluster_matrix) <- cluster_means$cluster

# Annotation rules based on canonical markers
cluster_labels <- rep("Unknown", nrow(cluster_matrix))

cluster_labels[cluster_matrix[,"CD3"] > 0.3 & cluster_matrix[,"CD4"] > 0.2] <- "CD4 T cell"
cluster_labels[cluster_matrix[,"CD3"] > 0.3 & cluster_matrix[,"CD8"] > 0.2] <- "CD8 T cell"
cluster_labels[cluster_matrix[,"CD19"] > 0.2] <- "B cell"
cluster_labels[cluster_matrix[,"CD56"] > 0.2] <- "NK cell"
cluster_labels[cluster_matrix[,"CD14"] > 0.2] <- "Monocyte"
cluster_labels[cluster_matrix[,"CD33"] > 0.2] <- "Myeloid"

# Map back to single cells
sce$cell_type <- cluster_labels[as.numeric(sce$cluster_id)]

# UMAP colored by annotated cell type
plotUMAP(sce, colour_by = "cell_type")

# -------------------------------
# 7. Differential marker analysis
# -------------------------------
markers <- colnames(expr_all)

# Mean expression per cluster
cluster_means <- expr_df %>%
  group_by(cluster) %>%
  summarise(across(.cols = where(is.numeric), .fns = mean))

cluster_matrix <- as.matrix(cluster_means[,-1])
rownames(cluster_matrix) <- cluster_means$cluster

# Heatmap of all markers
pheatmap(cluster_matrix,
         scale="row",
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         main="Differential Marker Expression Heatmap")

# Boxplots per marker
for(marker in markers){
  p <- ggplot(expr_df, aes(x=cluster, y=.data[[marker]], fill=cluster)) +
    geom_boxplot(outlier.size=0.5) +
    theme_minimal() +
    labs(title=paste("Expression of", marker, "per cluster"), y=marker) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
  
  ggsave(filename = paste0("Boxplot_", marker, ".png"), plot=p, width=10, height=6)
}

# Volcano-like differential expression per marker
volcano_df <- data.frame()
for(marker in markers){
  for(cl in unique(expr_df$cluster)){
    expr_cl <- expr_df[expr_df$cluster==cl, marker]
    expr_rest <- expr_df[expr_df$cluster!=cl, marker]
    
    log2FC <- log2(mean(expr_cl)+1e-5) - log2(mean(expr_rest)+1e-5)
    pval <- t.test(expr_cl, expr_rest)$p.value
    
    volcano_df <- rbind(volcano_df,
                        data.frame(cluster=cl, marker=marker, log2FC=log2FC, pval=pval))
  }
}

# Example volcano plot
marker_example <- "CD3"
df_plot <- volcano_df[volcano_df$marker==marker_example,]

ggplot(df_plot, aes(x=log2FC, y=-log10(pval), label=cluster)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 100) +
  theme_minimal() +
  labs(title=paste("Differential Expression Volcano for", marker_example),
       x="log2 Fold Change", y="-log10(p-value)")


