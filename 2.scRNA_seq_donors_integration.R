###################################################################################################
## Filter out cluster 3, debris cells , setting 1
###################################################################################################
# Load raw object
cMoP <- readRDS("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/cMoP_raw.rds")
# Remove stressed cluster (cluster 3)
cMoP <- subset(cMoP, idents = 3,invert=TRUE) ### exclude cluster 3 (stressed cells)

# Reset assay and clean metadata
DefaultAssay(cMoP) <- "RNA"
meta <- cMoP@meta.data[, !grepl("SCT_snn_res", colnames(cMoP@meta.data))]
cMoP <- JoinLayers(cMoP)
counts = LayerData(cMoP, assay = "RNA", layer = "counts")

# Recreate object with clean metadata
cMoP <- CreateSeuratObject(counts = counts,meta.data = meta)
# Normalize and split by donor (for integration)
DefaultAssay(cMoP) <- "RNA"
cMoP <- JoinLayers(cMoP)
cMoP <- NormalizeData(cMoP)
cMoP[["RNA"]] <- split(cMoP[["RNA"]], f = cMoP$donor)

# SCTransform with cell cycle regression
vars.to.regress = c("S.Score","G2M.Score")
cMoP <- SCTransform(cMoP, vars.to.regress = vars.to.regress, verbose = FALSE)

# PCA + clustering + UMAP (pre-integration)
cMoP <- RunPCA(cMoP, verbose = FALSE)
cMoP <- FindNeighbors(cMoP, dims = 1:50, reduction = "pca")
cMoP <- FindClusters(cMoP, resolution = 0.1, cluster.name = "unintegrated.clusters", algorithm = 4)
cMoP <- RunUMAP(cMoP, dims = 1:50, reduction.name = "umap.unintegrated")
DimPlot_scCustom(cMoP,group.by = "unintegrated.clusters",reduction = "umap.unintegrated",split.by = "donor")
FeatureDimPlot(cMoP,features = c("S.Score","G2M.Score"),reduction = "umap.unintegrated")
CellStatPlot(srt =cMoP, stat.by = "donor", group.by = "unintegrated_clusters", label = TRUE)

# ------------------------
# Harmony integration step
# ------------------------
cMoP <- IntegrateLayers(
  object = cMoP,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  verbose = FALSE,
  new.reduction = "integrated.sct.harmony",
  reduction = "pca",
  dims = 1:50,
  group.by = "donor"
)

# Run neighbors and UMAP on Harmony
cMoP <- FindNeighbors(cMoP, dims = 1:50, reduction = "integrated.sct.harmony")
cMoP <- RunUMAP(cMoP, dims = 1:50, reduction = "integrated.sct.harmony", reduction.name = "umap.harmony")

# Clustering at multiple resolutions
resolutions <- seq(0.1, 0.8, by = 0.1)

for (res in resolutions) {
  cMoP <- FindClusters(cMoP, resolution = res, algorithm = 4, random.seed = 123)
  colnames(cMoP@meta.data)[ncol(cMoP@meta.data)] <- paste0("SCT_snn_res.", res)
}


res <- sort(grep("^SCT_snn_res\\.", colnames(cMoP@meta.data), value = TRUE))
plot_list <- list()
for(i in 1:length(res)){
  plot_list[[i]] <- DimPlot_scCustom(cMoP, reduction = "umap.harmony",label=TRUE,label.size = 6,repel=TRUE,pt.size = 0.01,group.by = res[[i]]) & NoLegend()
}

gridExtra::grid.arrange(grobs = plot_list, ncol=4)

saveRDS(cMoP,file="/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/cMoP_integrated.rds")
