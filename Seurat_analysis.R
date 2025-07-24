##############################################
# R Script for cMoP Single-cell RNA-seq Analysis
# Annotated and Reviewed: Libraries, Comments
##############################################

# Load required libraries
library(Seurat)
library(BuenColors)
library(DoubletFinder)
library(DropletUtils)
library(SCP)
library(SeuratWrappers)
library(SingleCellExperiment)
library(SoupX)
library(UCell)
library(biomaRt)
library(cluster)
library(clustree)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(scCustomize)
library(scDblFinder)
library(viridis)
library(patchwork)
library(Matrix)
library(DESeq2)
library(pheatmap)
library(WriteXLS)
library(EnhancedVolcano)

# Set Seurat options
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 3e+12)

# Source custom helper scripts
source("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Helper_functions.R")
source("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/cellcycle.R")

# Load Seurat object
cMoP = readRDS("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered/merged_cMoP_seurat.rds")
cMoP <- JoinLayers(cMoP)

# Add quality control and gene content metrics
cMoP <- Add_Cell_QC_Metrics(cMoP, species ="human")
cMoP <- AddTFPercentToMetadata(cMoP)
cMoP <- AddLncRNApercentToMetadata(cMoP)

# Extract donor and sample identifiers
cMoP$donor <- gsub(".*cMoP_", "", cMoP$Sample_Name)
cMoP$cells <- gsub("_HD.*", "", cMoP$Sample_Name)

# Normalize, reduce dimensions and cluster
cMoP <- SCTransform(cMoP, vars.to.regress = NULL, verbose = FALSE)
cMoP <- RunPCA(cMoP, verbose = FALSE)
cMoP <- RunUMAP(cMoP, dims = 1:30)
cMoP <- FindNeighbors(cMoP, dims = 1:30)

# Cluster over a resolution sweep
resolutions <- seq(0.1, 0.8, by = 0.1)
for (res in resolutions) {
  cMoP <- FindClusters(cMoP, resolution =res, algorithm = 4, random.seed=123)
}

# Visualize cluster resolution tree
clustree(cMoP, prefix = "SCT_snn_res.", show_axis = TRUE)

# Generate UMAP plots for each resolution
res <- sort(grep("^SCT_snn_res\\.", colnames(cMoP@meta.data), value = TRUE))
plot_list <- list()
for (i in seq_along(res)) {
  plot_list[[i]] <- DimPlot_scCustom(cMoP, reduction = "umap", label = TRUE, label.size = 6,
                                     repel = TRUE, pt.size = 0.01, group.by = res[[i]]) & NoLegend()
}
gridExtra::grid.arrange(grobs = plot_list, ncol = 4)

# Marker analysis
Idents(cMoP) <- "SCT_snn_res.0.4"

# Cell cycle and QC plots
p1 <- CellDimPlot(cMoP, group.by = "cells", reduction = "UMAP", theme_use = "theme_blank", palcolor = c("#9EB8FC", "#D600D6"))
p2 <- CellDimPlot(cMoP, group.by = "cells", reduction = "UMAP", theme_use = "theme_blank", split.by = "cells",
                  palcolor = c("#9EB8FC", "#D600D6"), bg_color = "white")

# Save Figure 1C
pdf("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Figure1C_A.pdf", width = 5, height = 5, useDingb = FALSE)
print(p1)
dev.off()

pdf("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Figure1C_B.pdf", width = 5, height = 5, useDingb = FALSE)
print(p2)
dev.off()

# Run pseudotime inference with Monocle3
cMoP <- RunMonocle3(srt = cMoP, reduction = "umap", assay = "SCT", clusters = "SCT_snn_res.0.4")

# Compare pseudotime across groups
cMoP$cells <- factor(cMoP$cells, levels = c("CLEC12A_cMoP", "CD115_cMoP"))
v2 <- VlnPlot_scCustom(cMoP, features = "Monocle3_Pseudotime", group.by = "cells", plot_boxplot = TRUE)
v2 + ggpubr::stat_compare_means()

pdf("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/FigureSE.pdf", width = 4, height = 5, useDingb = FALSE)
print(v2)
dev.off()

##############################################
# Pseudo-bulk DE analysis between conditions
##############################################
cMoP = readRDS("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered/merged_cMoP_seurat.rds")
cMoP <- JoinLayers(cMoP)
cMoP$donor <- gsub(".*cMoP_", "", cMoP$Sample_Name)
cMoP$cells <- gsub("_HD.*", "", cMoP$Sample_Name)
cMoP$group_id <- paste(cMoP$donor, cMoP$cells, sep = "_")

# Generate pseudobulk matrix
counts <- GetAssayData(cMoP, assay = "RNA", slot = "counts")
meta <- cMoP@meta.data
group_ids <- unique(meta$group_id)

pseudobulk_counts <- sapply(group_ids, function(g) {
  cell_names <- rownames(meta)[meta$group_id == g]
  if (length(cell_names) == 1) {
    counts[, cell_names]
  } else {
    Matrix::rowSums(counts[, cell_names])
  }
})

pseudobulk_counts <- as.matrix(pseudobulk_counts)
colnames(pseudobulk_counts) <- group_ids
rownames(pseudobulk_counts) <- rownames(counts)

group_df <- data.frame(
  group_id = colnames(pseudobulk_counts),
  donor = sapply(strsplit(colnames(pseudobulk_counts), "_"), `[`, 1)
)
group_df$celltype <- gsub("^[^_]+_", "", group_df$group_id)
group_df$donor <- factor(group_df$donor)
group_df$celltype <- factor(group_df$celltype)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(pseudobulk_counts),
  colData = group_df,
  design = ~ donor + celltype
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("celltype", "CD115_cMoP", "CLEC12A_cMoP"), alpha = 0.05)
summary(res)

# Filter significant results
res_sig <- subset(res, padj <= 0.05)
res_sig[order(res_sig$log2FoldChange), ]
table(abs(res_sig$log2FoldChange) > 0.5)

# Export results
WriteXLS::WriteXLS(as.data.frame(res_sig), ExcelFileName = "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Kawamura_vs_Calzetti.xls", row.names = TRUE)

# Apply fold-change filters
res_sig_1.5 <- subset(res, padj <= 0.05 & abs(res$log2FoldChange) > log2(1.5))

# Plot volcano
p4 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      FCcutoff = log2(1.5),
                      pCutoff = 0.05,
                      xlim = c(-1.8, 1.8),
                      ylim = c(-0.2, 9.5),
                      parseLabels = TRUE,
                      legendLabels = c('Not sig.', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
                      legendPosition = 'bottom',
                      pointSize = 1)
p4

pdf("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/FigureS1C.pdf", width = 8, height = 8, useDingb = FALSE)
print(p4)
dev.off()

# Correlation matrix plot
rld <- rlog(dds)
cor_mat <- cor(assay(rld), method = "pearson")
p3 <- pheatmap(cor_mat,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
               breaks = seq(0.990, 1, length.out = 101),
               main = "CLEC12A_cMoP and CD115_cMoP - Pearson Correlation")

pdf("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Figure1D.pdf", width = 5, height = 5, useDingb = FALSE)
print(p3)
dev.off()

