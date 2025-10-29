###################################################################################################
## 1. — Remove ambient (environmental) RNA contamination with SoupX
###################################################################################################

library(Seurat)
library(SoupX)
library(DropletUtils)

# Define input/output directories
filt_dir   <- "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Default_setting/sid_3439/sid-3439_RSEC_MolsPerCell_MEX/"
raw_dir    <- "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Default_setting/sid_3439/RSEC_MolsPerCell_Unfiltered_MEX/"
output_dir <- "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered/soupX_filt"

# Read filtered (cells) and raw (droplets) matrices
toc <- Read10X(filt_dir)   # table of counts (filtered)
tod <- Read10X(raw_dir)    # table of droplets (raw)

# Ensure gene consistency between matrices
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]

# Create a Seurat object for quick clustering (used by SoupX)
srat <- SCTransform(CreateSeuratObject(toc), verbose = FALSE) |>
  RunPCA(verbose = FALSE) |>
  RunUMAP(dims = 1:30, verbose = FALSE) |>
  FindNeighbors(dims = 1:30, verbose = FALSE) |>
  FindClusters(verbose = TRUE)

# Initialize SoupChannel and add cluster and UMAP info
soup <- SoupChannel(tod, toc)
soup <- setClusters(soup, setNames(srat$seurat_clusters, Cells(srat)))
soup <- setDR(soup, Embeddings(srat, "umap"))

# Estimate and adjust for ambient RNA contamination
soup <- autoEstCont(soup)
adj  <- adjustCounts(soup, roundToInt = TRUE)

# Write SoupX-corrected count matrix to disk
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
DropletUtils:::write10xCounts(output_dir, adj)

###################################################################################################
## 2. Load SoupX-filtered data & metadata
###################################################################################################
setwd("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered")
sc <- load_sample_data("3439")

sc$donor <- gsub(".*_", "", sc$Sample_Name)
sc$cells <- gsub("_[^_]*$", "", sc$Sample_Name)
sc <- subset(sc, subset = !(Sample_Name %in% c("Multiplet", "Undetermined")))

###################################################################################################
## 3. Subset, QC, normalization
###################################################################################################
cMoP <- subset(sc, subset = cells %in% c("cMoP_FC", "cMoP_Kawa"))
cMoP$cells <- gsub("cMoP_FC", "CD115_cMoP", cMoP$cells)
cMoP$cells <- gsub("cMoP_Kawa", "CLEC12A_cMoP", cMoP$cells)
cMoP$cells <- factor(cMoP$cells, levels = c("CLEC12A_cMoP", "CD115_cMoP"))

cMoP <- Add_Cell_QC_Metrics(cMoP, species = "human", add_cell_cycle = TRUE, overwrite = TRUE)
cMoP <- AddTFPercentToMetadata(cMoP)
cMoP <- AddLncRNApercentToMetadata(cMoP)
cMoP <- subset(cMoP, subset = percent_mito < 20)
cMoP$CC.Difference <- cMoP$S.Score - cMoP$G2M.Score

Vln_list <- lapply(c("percent_lncRNA", "percent_TF"), 
                   \(f) VlnPlot_scCustom(cMoP, features = f, plot_boxplot = TRUE, group.by = "Sample_Name") + NoLegend())
wrap_plots(Vln_list, ncol = 2)
FeatureScatter_scCustom(cMoP, feature1 = "percent_lncRNA", feature2 = "percent_TF")

###################################################################################################
## 4. Normalization, dimensionality reduction, clustering
###################################################################################################
DefaultAssay(cMoP) <- "RNA"
cMoP <- SCTransform(cMoP, vars.to.regress = c("S.Score","G2M.Score"), verbose = FALSE)
cMoP <- RunPCA(cMoP) |> FindNeighbors(dims = 1:50) |> 
  FindClusters(resolution = 0.1, cluster.name = "unintegrated_clusters", algorithm = 4) |>
  RunUMAP(dims = 1:50, reduction.name = "umap.unintegrated")

DimPlot_scCustom(cMoP, group.by = "unintegrated_clusters", split.by = "donor")

###################################################################################################
## 5. Harmony integration & reclustering
###################################################################################################
cMoP <- IntegrateLayers(
  object = cMoP, method = HarmonyIntegration, normalization.method = "SCT",
  new.reduction = "integrated.sct.harmony", reduction = "pca",
  dims = 1:50, group.by = "donor"
)

cMoP <- FindNeighbors(cMoP, dims = 1:50, reduction = "integrated.sct.harmony") |>
  RunUMAP(dims = 1:50, reduction = "integrated.sct.harmony", reduction.name = "umap.harmony")

resolutions <- seq(0.1, 0.8, by = 0.1)
for (r in resolutions) cMoP <- FindClusters(cMoP, resolution = r, algorithm = 4, random.seed = 123)

plots <- lapply(grep("^SCT_snn_res\\.", colnames(cMoP@meta.data), value = TRUE),
                \(r) DimPlot_scCustom(cMoP, reduction = "umap.harmony", label = TRUE, repel = TRUE, pt.size = 0.01, group.by = r) & NoLegend())
gridExtra::grid.arrange(grobs = plots, ncol = 4)

###################################################################################################
## 6. Marker discovery & visualization
###################################################################################################
Idents(cMoP) <- "SCT_snn_res.0.1"
cMoP <- PrepSCTFindMarkers(cMoP)
markers <- FindAllMarkers(cMoP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- markers %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>% slice_head(n = 15)
DotPlot(cMoP, features = unique(top_markers$gene), dot.scale = 4, cluster.idents = FALSE) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(100, "RdBu"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

###################################################################################################
## 6. Save object
###################################################################################################
saveRDS(cMoP, file = "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/cMoP_raw.rds")

