########################################################################################################################
# cMoP PNAS Figure Generation Script — Final Polished Version
# Author: Francisco Bianchetto
# Description:
#   Generates publication-quality figures (Figure 1C and Figure 2A–D)
#   from integrated and velocity-processed Seurat objects.
########################################################################################################################

# --------------------------------------------
# Setup
# --------------------------------------------
options(Seurat.object.assay.version = "v3")  # Ensure Seurat v3 object compatibility
library(Seurat)
library(RColorBrewer)
library(SCP)

# Output directory
outdir <- "/home/patgen/OneDrive/cMoP/PNAS/Plots"

# Color palettes
pal_cells     <- c("#D600D6", "#9EB8FC")  # Figure 1C
pal_clusters  <- c("#7F3C8D", "#E69F00", "#0072B2", "#009E73", "#8B0000")  # Figures 2A–D

# --------------------------------------------
# Figure 1C — UMAP by cell type (Harmony reduction)
# --------------------------------------------

# Load integrated cMoP dataset
cMoP <- readRDS("/home/patgen/OneDrive/cMoP/PNAS/R_script/cMoP_integrated.rds")

# UMAP colored by "cells" (merged view)
p1 <- CellDimPlot(
  cMoP,
  group.by   = "cells",
  reduction  = "umap.harmony",
  theme_use  = "theme_blank",
  palcolor   = pal_cells
)

# Split UMAP by "cells" (two panels)
p2 <- CellDimPlot(
  cMoP,
  group.by   = "cells",
  reduction  = "umap.harmony",
  theme_use  = "theme_blank",
  split.by   = "cells",
  palcolor   = pal_cells,
  bg_color   = "white"
)

# Save plots
pdf(file.path(outdir, "Figure_1C_A_UMAP_merged.pdf"), width = 5, height = 5, useDingb = FALSE)
print(p1)
dev.off()

pdf(file.path(outdir, "Figure_1C_B_UMAP_split.pdf"), width = 8, height = 8, useDingb = FALSE)
print(p2)
dev.off()

# --------------------------------------------
# Figure 2A — Cluster UMAP (Harmony reduction)
# --------------------------------------------

# Load integrated object with cluster assignments
seurat_clusters <- readRDS("/home/patgen/OneDrive/cMoP/PNAS/R_script/integrated_seurat_harmony_selected_setting_clusters.rds")

Idents(seurat_clusters) <- "SCT_snn_res.0.6"
seurat_clusters$clusters <- factor(paste0("c", seurat_clusters$SCT_snn_res.0.6),
                                   levels = c("c1", "c2", "c3", "c4", "c5"))

pdf(file.path(outdir, "Figure_2A_UMAP_clusters.pdf"), width = 4, height = 5, useDingb = FALSE)
CellDimPlot(
  seurat_clusters,
  reduction     = "umap.harmony",
  group.by      = "clusters",
  theme_use     = "theme_blank",
  label_insitu  = FALSE,
  label         = FALSE,
  palcolor      = pal_clusters,
  pt.size       = 0.1
)
dev.off()

# --------------------------------------------
# Figure 2B — DotPlot of selected marker genes
# --------------------------------------------

# Marker genes for cluster comparison
genes <- c(
  "PRTN3","ELANE","AZU1","MS4A3","MPO","HGF","MYB","F13A1","NUCB2","TFRC",
  "FKBP4","NME1","HSPH1","ABCE1","HSPE1","RCC1","PAICS","DKC1","GNL3","NOLC1",
  "S100A10","HLA-DQA1","CIITA","CD86","HLA-DOA","HLA-DMB","HLA-DQB1","S100A6",
  "S100A4","CD74","ISG15","IFIT3","MX1","IFI44L","MX2","IFIT1","OAS2","RSAD2",
  "IFIT2","OAS3","VCAN","S100A8","NCF2","S100A9","CYBB","CCR2","FGR","CCR1",
  "CTSH","LYZ"
)

Idents(seurat_clusters) <- "clusters"

pdf(file.path(outdir, "Figure_2B_DotPlot_DEGs.pdf"), width = 5, height = 10, useDingb = FALSE)
DotPlot(seurat_clusters, features = genes, dot.scale = 4, cluster.idents = FALSE) +
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)) +
  theme(
    axis.title   = element_blank(),
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 1, size = 10),
    axis.text.y  = element_text(color = "grey20", size = 12),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
dev.off()

# --------------------------------------------
# Figure 2C — RNA Velocity UMAP
# --------------------------------------------

# Load Seurat object with RNA velocity results
seurat_velocity <- readRDS("/home/patgen/OneDrive/cMoP/PNAS/R_script/seurat_velocity_computed.rds")

seurat_velocity$clusters <- factor(paste0("c", seurat_velocity$SCT_snn_res.0.6),
                                   levels = c("c1", "c2", "c3", "c4", "c5"))

# Background UMAP for velocity overlay
pdf(file.path(outdir, "Figure_2C_UMAP_velocity_background.pdf"), width = 4, height = 5, useDingb = FALSE)
CellDimPlot(
  seurat_velocity,
  reduction     = "umap.harmony",
  group.by      = "clusters",
  theme_use     = "theme_blank",
  label_insitu  = FALSE,
  label         = FALSE,
  palcolor      = pal_clusters,
  pt.size       = 0.1,
  pt.alpha      = 0.2
)
dev.off()

# Velocity streamlines
velocito <- VelocityPlot(
  srt                 = seurat_velocity,
  velocity            = "dynamical",
  reduction           = "umap.harmony",
  dims                = c(1, 2),
  plot_type           = "stream",
  streamline_L        = 5,
  streamline_minL     = 1,
  streamline_res      = 1,
  streamline_n        = 10,
  streamline_width    = c(0, 1),
  streamline_alpha    = 1,
  streamline_palette  = "RdYlBu",
  streamline_bg_color = "black",
  streamline_bg_stroke= 0.5,
  cutoff_perc         = 5,
  density             = 1,
  smooth              = 1,
  scale               = 1,
  group_by            = "clusters",
  group_palette       = pal_clusters,
  aspect.ratio        = 1,
  title               = "Cell velocity",
  legend.position     = "right",
  theme_use           = "theme_scp",
  seed                = 11
)

pdf(file.path(outdir, "Figure_2C_UMAP_velocity_vectors.pdf"), width = 4, height = 5, useDingb = FALSE)
velocito
dev.off()

# --------------------------------------------
# Figure 2D — Dynamical pseudotime distribution
# --------------------------------------------

# Reorder clusters according to pseudotime progression
seurat_velocity$clusters <- factor(
  seurat_velocity$SCT_snn_res.0.6,
  levels = c("c1", "c2", "c4", "c5", "c3")
)

# Plot pseudotime distribution across clusters
pdf(file.path(outdir, "Figure_2D_Pseudotime.pdf"), width = 4, height = 5, useDingb = FALSE)
FeatureStatPlot(
  srt          = seurat_velocity,
  group.by     = "clusters",
  stat.by      = "dynamical_pseudotime",
  add_box      = TRUE,
  comparisons  = list(c("c1","c2"), c("c2","c4"), c("c4","c5"), c("c5","c3")),
  palcolor     = pal_clusters,
  ylab         = "Dynamical pseudotime",
  title        = "RNA Velocity: Dynamical Pseudotime"
) + NoLegend()
dev.off()

