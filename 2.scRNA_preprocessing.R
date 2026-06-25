library(Seurat)

###################################################################################################
## 1. Load SoupX-filtered data & metadata
###################################################################################################
source("/home/patgen/OneDrive/cMoP/Cell_Reports/Toward-a-more-precise-immunophenotypic-definition-of-human-cMoPs/Helper_functions.R")
sc <- readRDS(file = "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Default_setting/sid_3439/sid-3439_Seurat_SoupX_corrected.rds")

sc$donor <- gsub(".*_", "", sc$Sample_Name)
sc$cells <- gsub("_[^_]*$", "", sc$Sample_Name)
sc <- subset(sc, subset = !(Sample_Name %in% c("Multiplet", "Undetermined")))

###################################################################################################
## 2. Subset, QC, normalization
###################################################################################################
cMoP <- subset(sc, subset = cells %in% c("cMoP_FC", "cMoP_Kawa"))
cMoP$cells <- gsub("cMoP_FC", "CD115_cMoP", cMoP$cells)
cMoP$cells <- gsub("cMoP_Kawa", "CLEC12A_cMoP", cMoP$cells)
cMoP$cells <- factor(cMoP$cells, levels = c("CLEC12A_cMoP", "CD115_cMoP"))
cMoP <- subset(cMoP,subset = DoubletFinder =="Singlet")

cMoP <- Add_Cell_QC_Metrics(cMoP, species = "human", add_cell_cycle = TRUE, overwrite = TRUE)
cMoP <- AddTFPercentToMetadata(cMoP)
cMoP <- AddLncRNApercentToMetadata(cMoP)
cMoP$CC.Difference <- cMoP$S.Score - cMoP$G2M.Score

Vln_list <- lapply(c("percent_lncRNA", "percent_TF"), 
                   \(f) VlnPlot_scCustom(cMoP, features = f, plot_boxplot = TRUE, group.by = "Sample_Name") + NoLegend())
wrap_plots(Vln_list, ncol = 2)
VlnPlot(cMoP,features = c("percent_mito","nFeature_RNA"),group.by = "donor")
cMoP <- filter_genes(cMoP)
cMoP <- subset(cMoP, subset = nFeature_RNA > 500 & percent_mito < 20)
###################################################################################################
## 3. Normalization, dimensionality reduction, clustering
###################################################################################################
DefaultAssay(cMoP) <- "RNA"
cMoP <- JoinLayers(cMoP)
cMoP[["RNA"]] <- split(cMoP[["RNA"]], f = cMoP$donor)
cMoP <- SCTransform(cMoP, vars.to.regress = c("S.Score","G2M.Score"), verbose = FALSE)
cMoP <- RunPCA(cMoP) |> FindNeighbors(dims = 1:50) |> 
  FindClusters(resolution = 0.1, cluster.name = "unintegrated_clusters", algorithm = 4) |>
  RunUMAP(dims = 1:50, reduction.name = "umap.unintegrated")

DimPlot_scCustom(cMoP, group.by = "unintegrated_clusters", split.by = "donor")

###################################################################################################
## 4. Harmony integration & reclustering
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
## 5. Marker discovery & visualization
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

