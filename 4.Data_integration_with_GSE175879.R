###################################################################################################
## MERGE AND CLEAN SEURAT OBJECTS (donor-specific, SoupX-filtered, DoubletFinder-processed)
###################################################################################################
### DATA from GSE175879
# Load donor-specific Seurat objects (already filtered by SoupX and DoubletFinder)
donor_list <- readRDS(
  file = "/home/patgen/working_dir/Data_analysis/SingleCells_NCPs_cMoPs/donor_list_soupx_doubletFinder.rds")

# Keep only cMoP singlets from each donor
donor_list <- lapply(donor_list, function(x) {
  subset(x, subset = Sample_Name == "cMoP" & DoubletFinder == "Singlet")
})

# Annotate donor IDs manually
donor_list[[1]]$donor <- "4"
donor_list[[2]]$donor <- "5"
donor_list[[3]]$donor <- "6"

# Annotate sequencing runs
donor_list[[1]]$run <- "1"
donor_list[[2]]$run <- "2"
donor_list[[3]]$run <- "3"

###################################################################################################
## MERGE DONORS FROM RUNS 1–3
###################################################################################################

seurat_runs123 <- merge(donor_list[[1]],y = c(donor_list[[2]], donor_list[[3]]),add.cell.ids = c("run1", "run2", "run3"))

###################################################################################################
## LOAD AND CLEAN THE INTEGRATED cMoP OBJECT (RUN4)
###################################################################################################

cMoP <- readRDS(file = "/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/cMoP_integrated.rds")

# Visualize integrated object before subsetting
DimPlot_scCustom(cMoP, reduction = "umap.harmony")

# Subset to specific population
cMoP <- subset(cMoP, subset = cells == "CD115_cMoP")

# Reset assay and clean metadata
DefaultAssay(cMoP) <- "RNA"

# Remove SCT clustering resolutions and retain first 7 metadata columns
meta <- cMoP@meta.data[, !grepl("SCT_snn_res", colnames(cMoP@meta.data))]
meta <- meta[, 1:7]

# Join layers (ensures consistent assay structure)
cMoP <- JoinLayers(cMoP)

# Extract raw RNA counts
counts <- LayerData(cMoP, assay = "RNA", layer = "counts")

###################################################################################################
## REBUILD SEURAT OBJECT FOR RUN4 WITH CLEAN METADATA
###################################################################################################

seurat_run3 <- CreateSeuratObject(counts = counts,meta.data = meta)

# Add QC metrics and cell cycle scoring
seurat_run3 <- scCustomize::Add_Cell_QC_Metrics(seurat_run3,species = "human",add_cell_cycle = TRUE,overwrite = TRUE,top_pct_name = "TopGene")

# Annotate metadata
seurat_run3$run <- "run4"
colnames(seurat_run3) <- paste("run4", colnames(seurat_run3), sep = "_")

# Harmonize donor and run annotations with previous objects
seurat_run3$donor <- paste0("donor", seurat_run3$donor)
seurat_runs123$donor <- paste0("donor", seurat_runs123$donor)

# Add population labels and format run IDs
seurat_runs123$cells <- "CD115_cMoP"
seurat_runs123$run <- paste0("run", seurat_runs123$run)

###################################################################################################
## MERGE RUN4 WITH RUNS 1–3 AND FINALIZE METADATA
###################################################################################################

merged_seurat <- merge(seurat_runs123,y = seurat_run3,add.cell.ids = NULL)

DefaultAssay(merged_seurat) <- "RNA"

# Join layers again to ensure full assay compatibility
merged_seurat <- JoinLayers(merged_seurat)

# Create batch variable (run + donor)
merged_seurat$batch <- paste(merged_seurat$run, merged_seurat$donor, sep = "_")

# Remove redundant suffixes
merged_seurat$batch <- gsub("_cMoP", "", merged_seurat$batch)

# Save merged Seurat object
saveRDS(merged_seurat, file = "/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/merged_seurat.rds")

###################################################################################################
## 4. SCTransform NORMALIZATION
###################################################################################################

filtered_seurat = readRDS("/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/merged_seurat.rds")

# Split object by batch for SCTransform normalization
seurat_list <- SplitObject(filtered_seurat, split.by = "batch")

# Run SCTransform on each batch (regressing out technical variables)
seurat_list <- lapply(seurat_list, function(x) {
  x <- SCTransform(
    x,
    method = "glmGamPoi",
    vars.to.regress = c("percent_mito", "S.Score", "G2M.Score", "percent_ribo"),
    verbose = TRUE
  )
  return(x)
})

###################################################################################################
## 5. INTEGRATION PREPARATION
###################################################################################################

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

# Prepare SCT integration
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

# Run PCA on each object
seurat_list <- lapply(seurat_list, function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

###################################################################################################
## 6. HARMONY INTEGRATION
###################################################################################################

# Merge objects before Harmony integration
integrated_seurat <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)])

# Set default assay to SCT
DefaultAssay(integrated_seurat) <- "SCT"

# Run PCA and Harmony integration
integrated_seurat <- RunPCA(integrated_seurat, features = features, verbose = FALSE)
integrated_seurat <- RunHarmony(integrated_seurat,group.by.vars = c("run", "donor"))

###################################################################################################
## 7. DIMENSIONALITY REDUCTION AND CLUSTERING
###################################################################################################

# Run UMAP using Harmony embeddings
integrated_seurat <- RunUMAP(integrated_seurat,
  reduction = "harmony",dims = 1:30,
  reduction.name = "umap_harmony"
)

# Find neighbors and clusters
integrated_seurat <- FindNeighbors(
  integrated_seurat,
  reduction = "harmony",
  dims = 1:30
)

# Explore multiple clustering resolutions
resolutions <- seq(0.1, 0.8, by = 0.1)
for (res in resolutions) {
integrated_seurat <- FindClusters(integrated_seurat,resolution = res,algorithm = 4)
}

# Visualize all resolutions
res <- sort(grep("^SCT_snn_res\\.", colnames(integrated_seurat@meta.data), value = TRUE))
plot_list <- list()
for (i in seq_along(res)) {
  plot_list[[i]] <- DimPlot_scCustom(
    integrated_seurat,
    reduction = "umap_harmony",
    label = TRUE,label.size = 6,repel = TRUE,pt.size = 0.01,group.by = res[[i]]) & NoLegend()
}
gridExtra::grid.arrange(grobs = plot_list, ncol = 4)

# Save integrated object
saveRDS(integrated_seurat,file = "/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/integrated_seurat_harmony_selected_setting.rds")

###################################################################################################
## 8. CLUSTER VISUALIZATION AND MARKER DISCOVERY
###################################################################################################
library(scCustomize)

# Load integrated object
integrated_seurat <- readRDS(
  file = "/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/integrated_seurat_harmony_selected_setting.rds")

# Visualize selected resolution
DimPlot_scCustom(integrated_seurat,reduction = "umap_harmony",group.by = c("SCT_snn_res.0.6"),label = TRUE,label.size = 8) + NoLegend()

# Set cluster identity and prepare for marker detection
Idents(integrated_seurat) <- "SCT_snn_res.0.6"
integrated_seurat <- PrepSCTFindMarkers(integrated_seurat)

# Find cluster markers
all_markers <- FindAllMarkers(integrated_seurat,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25,test.use = "wilcox",min.diff.pct = -Inf)

# Extract top markers per cluster
top_markers <- all_markers %>%group_by(cluster) %>% filter(avg_log2FC > 1) %>%slice_head(n = 20) %>% ungroup()

# Visualize top markers
DotPlot(integrated_seurat,
  features = unique(top_markers$gene),dot.scale = 4,cluster.idents = FALSE) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 100, name = "RdBu"))) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(color = "grey20", size = 12),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

