# Load required libraries
library(Seurat)       
 library(ggplot2)       
library(FNN)          
library(igraph)        # 
library(SCP)           # For RNA velocity and trajectory analysis
source("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Helper_functions.R")  # Custom helper functions
#########################################################################################################################
####### RNA velocity analysis

# Load the integrated Seurat object
integrated_seurat <- readRDS(file = "/home/patgen/OneDrive/cMoP/PNAS/R_script/integrated_seurat_harmony_selected_setting.rds")
seurat_obj <- integrated_seurat

# Ensure compatibility with Seurat v3 assays
options(Seurat.object.assay.version = "v3")
assay_names <- names(seurat_obj@assays)
for(assay_name in assay_names) {
  if(inherits(seurat_obj[[assay_name]], "Assay5")) {
    cat("Converting", assay_name, "assay\n")
    seurat_obj[[assay_name]] <- as(seurat_obj[[assay_name]], "Assay")
  }
}
DefaultAssay(seurat_obj) <- "RNA"

# Filter genes using a custom helper function
seurat_obj <- filter_genes(seurat_obj)

# Function to retain only protein-coding genes
keep_protein_coding <- function(seurat_obj, species = "hsapiens") {
  library(biomaRt)
  mart <- useMart("ensembl", dataset = paste0(species, "_gene_ensembl"))
  genes <- rownames(seurat_obj)
  
  attr <- ifelse(species == "hsapiens", "hgnc_symbol", "mgi_symbol")
  annot <- getBM(attributes = c(attr, "gene_biotype"),
                 filters = attr,
                 values = genes,
                 mart = mart)
  
  keep <- annot[[attr]][annot$gene_biotype == "protein_coding"]
  subset(seurat_obj, features = intersect(rownames(seurat_obj), keep))
}
# Keep only protein-coding genes
seurat_obj <- keep_protein_coding(seurat_obj)

# Set working directory for velocity analysis
setwd("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Cell/Velocity")

# Load required libraries again for safety
library(Seurat)
library(SeuratDisk)
library(SCP)

#### Load loom files for each donor/run and process

# --- Donor 1 ---
loom_data <- Connect(filename = "/home/patgen/working_dir/Data_analysis/Velocity/Velocity_input/velocyto/Combined_BM-d1_4XDRQ.loom", mode = "r")
print("Available datasets:")
print(names(loom_data))

# Extract spliced and unspliced matrices
spliced_raw <- loom_data[["layers/spliced"]]
unspliced_raw <- loom_data[["layers/unspliced"]]

# Print dimensions
print(paste("Spliced dimensions:", dim(spliced_raw)[1], "x", dim(spliced_raw)[2]))

# Convert to full matrices and transpose to match Seurat format (cells x genes)
spliced_matrix <- t(spliced_raw[,])
unspliced_matrix <- t(unspliced_raw[,])

# Extract cell and gene names
if("col_attrs" %in% names(loom_data) && "CellID" %in% names(loom_data[["col_attrs"]])) {
  cell_names <- loom_data[["col_attrs/CellID"]][]
} else {
  cell_names <- paste0("Cell_", 1:ncol(spliced_raw))
}

if("row_attrs" %in% names(loom_data) && "Gene" %in% names(loom_data[["row_attrs"]])) {
  gene_names <- loom_data[["row_attrs/Gene"]][]
} else {
  gene_names <- paste0("Gene_", 1:nrow(spliced_raw))
}

# Close loom connection
loom_data$close_all()

# Set proper row and column names
rownames(spliced_matrix)   <- gene_names
rownames(unspliced_matrix) <- gene_names

# Clean up column names and add run-specific prefix
clean_names <- gsub("x", "", gsub("Combined_BM-d1_4XDRQ:", "", cell_names))
clean_names <- paste('run1', clean_names, sep="_")
colnames(spliced_matrix)   <- clean_names
colnames(unspliced_matrix) <- clean_names

# Match the genes and cells to main Seurat object
genes <- rownames(seurat_obj)
cells <- colnames(seurat_obj)
common_cells <- intersect(cells, clean_names)
unspliced_matrix <- unspliced_matrix[genes, match(common_cells, clean_names)]
spliced_matrix   <- spliced_matrix[genes, match(common_cells, clean_names)]

# Create Seurat object for velocity
velocity_seurat_d1 <- CreateSeuratObject(counts = spliced_matrix, project = "velocity")
velocity_seurat_d1[["unspliced"]] <- CreateAssayObject(counts = unspliced_matrix)

# Subset and reorder cells to match main Seurat object
common_cells <- intersect(colnames(seurat_obj), colnames(velocity_seurat_d1))
ordered_cells <- colnames(seurat_obj)[colnames(seurat_obj) %in% common_cells]
seurat_obj_subset_d1 <- subset(seurat_obj, cells = ordered_cells)
velocity_seurat_d1 <- subset(velocity_seurat_d1, cells = ordered_cells)

# Verify matching order
all(colnames(seurat_obj_subset_d1) == colnames(velocity_seurat_d1))  # TRUE
all(rownames(seurat_obj_subset_d1) == rownames(velocity_seurat_d1))   # TRUE

# Add velocity assays to main object
seurat_obj_subset_d1[["spliced"]] <- velocity_seurat_d1[["RNA"]]
seurat_obj_subset_d1[["unspliced"]] <- velocity_seurat_d1[["unspliced"]]

#### run 2
loom_data <- Connect(filename = "/home/patgen/working_dir/Data_analysis/Velocity/Velocity_input/velocyto/Combined_BM-d2_ICMTB.loom", mode = "r")
# First, let's check what's available
print("Available datasets:")
print(names(loom_data))

# Extract the matrices
spliced_raw <- loom_data[["layers/spliced"]]
unspliced_raw <- loom_data[["layers/unspliced"]]

# Get dimensions
print(paste("Spliced dimensions:", dim(spliced_raw)[1], "x", dim(spliced_raw)[2]))

# Extract full matrices
spliced_matrix <- spliced_raw[,]
unspliced_matrix <- unspliced_raw[,]

# Transpose (loom format is genes x cells, Seurat expects cells x genes for input)
spliced_matrix <- t(spliced_matrix)
unspliced_matrix <- t(unspliced_matrix)

# Get cell names if available
if("col_attrs" %in% names(loom_data) && "CellID" %in% names(loom_data[["col_attrs"]])) {
  cell_names <- loom_data[["col_attrs/CellID"]][]
} else {
  # Create cell names
  cell_names <- paste0("Cell_", 1:ncol(spliced_raw))
}

# Get gene names if available
if("row_attrs" %in% names(loom_data) && "Gene" %in% names(loom_data[["row_attrs"]])) {
  gene_names <- loom_data[["row_attrs/Gene"]][]
} else {
  # Create gene names
  gene_names <- paste0("Gene_", 1:nrow(spliced_raw))
}

# Close the connection
loom_data$close_all()

# Set proper row and column names
rownames(spliced_matrix)   <- gene_names
rownames(unspliced_matrix) <- gene_names

# Clean up column names
clean_names <- gsub("x", "", gsub("Combined_BM-d2_ICMTB:", "", cell_names))
clean_names <- paste('run2',clean_names,sep="_")
colnames(spliced_matrix)   <- clean_names
colnames(unspliced_matrix) <- clean_names

genes = rownames(seurat_obj);length(genes)
cells = colnames(seurat_obj);length(cells)

common_cells <- intersect(cells, clean_names)
unspliced_matrix <- unspliced_matrix[genes, match(common_cells, clean_names)]
spliced_matrix   <- spliced_matrix[genes, match(common_cells, clean_names)]

# Now create the Seurat object
velocity_seurat_d2 <- CreateSeuratObject(counts = spliced_matrix, project = "velocity")
velocity_seurat_d2[["unspliced"]] <- CreateAssayObject(counts = unspliced_matrix)

# 1. Get intersection of cell names
common_cells <- intersect(colnames(seurat_obj), colnames(velocity_seurat_d2))

# 2. Define the order (using seurat_obj as reference)
ordered_cells <- colnames(seurat_obj)[colnames(seurat_obj) %in% common_cells]

# 3. Subset *and* reorder explicitly
seurat_obj_subset_d2 <- subset(seurat_obj, cells = ordered_cells)
velocity_seurat_d2 <- subset(velocity_seurat_d2, cells = ordered_cells)

# 4. Verify order
all(colnames(seurat_obj_subset_d2) == colnames(velocity_seurat_d2))
# Should return TRUE
all(rownames(seurat_obj_subset_d2)==rownames(velocity_seurat_d2))
# Add velocity assays to your main object
seurat_obj_subset_d2[["spliced"]] <- velocity_seurat_d2[["RNA"]]
seurat_obj_subset_d2[["unspliced"]] <- velocity_seurat_d2[["unspliced"]]

#### run 3
loom_data <- Connect(filename = "/home/patgen/working_dir/Data_analysis/Velocity/Velocity_input/velocyto/Combined_BM-d3_58B2Y.loom", mode = "r")
# First, let's check what's available
print("Available datasets:")
print(names(loom_data))

# Extract the matrices
spliced_raw <- loom_data[["layers/spliced"]]
unspliced_raw <- loom_data[["layers/unspliced"]]

# Get dimensions
print(paste("Spliced dimensions:", dim(spliced_raw)[1], "x", dim(spliced_raw)[2]))

# Extract full matrices
spliced_matrix <- spliced_raw[,]
unspliced_matrix <- unspliced_raw[,]

# Transpose (loom format is genes x cells, Seurat expects cells x genes for input)
spliced_matrix <- t(spliced_matrix)
unspliced_matrix <- t(unspliced_matrix)

# Get cell names if available
if("col_attrs" %in% names(loom_data) && "CellID" %in% names(loom_data[["col_attrs"]])) {
  cell_names <- loom_data[["col_attrs/CellID"]][]
} else {
  # Create cell names
  cell_names <- paste0("Cell_", 1:ncol(spliced_raw))
}

# Get gene names if available
if("row_attrs" %in% names(loom_data) && "Gene" %in% names(loom_data[["row_attrs"]])) {
  gene_names <- loom_data[["row_attrs/Gene"]][]
} else {
  # Create gene names
  gene_names <- paste0("Gene_", 1:nrow(spliced_raw))
}

# Close the connection
loom_data$close_all()

# Set proper row and column names
rownames(spliced_matrix)   <- gene_names
rownames(unspliced_matrix) <- gene_names

# Clean up column names
clean_names <- gsub("x", "", gsub("Combined_BM-d3_58B2Y:", "", cell_names))
clean_names <- paste('run3',clean_names,sep="_")
colnames(spliced_matrix)   <- clean_names
colnames(unspliced_matrix) <- clean_names

genes = rownames(seurat_obj);length(genes)
cells = colnames(seurat_obj);length(cells)

common_cells <- intersect(cells, clean_names)
unspliced_matrix <- unspliced_matrix[genes, match(common_cells, clean_names)]
spliced_matrix   <- spliced_matrix[genes, match(common_cells, clean_names)]

# Now create the Seurat object
velocity_seurat_d3 <- CreateSeuratObject(counts = spliced_matrix, project = "velocity")
velocity_seurat_d3[["unspliced"]] <- CreateAssayObject(counts = unspliced_matrix)

# 1. Get intersection of cell names
common_cells <- intersect(colnames(seurat_obj), colnames(velocity_seurat_d3))

# 2. Define the order (using seurat_obj as reference)
ordered_cells <- colnames(seurat_obj)[colnames(seurat_obj) %in% common_cells]

# 3. Subset *and* reorder explicitly
seurat_obj_subset_d3 <- subset(seurat_obj, cells = ordered_cells)
velocity_seurat_d3 <- subset(velocity_seurat_d3, cells = ordered_cells)

# 4. Verify order
all(colnames(seurat_obj_subset_d3) == colnames(velocity_seurat_d3))
# Should return TRUE
all(rownames(seurat_obj_subset_d3)==rownames(velocity_seurat_d3))
# Add velocity assays to your main object
seurat_obj_subset_d3[["spliced"]] <- velocity_seurat_d3[["RNA"]]
seurat_obj_subset_d3[["unspliced"]] <- velocity_seurat_d3[["unspliced"]]

### run 4
loom_data <- Connect(filename = "/home/patgen/working_dir/Data_analysis/Velocity/Velocity_input/velocyto/Combined_sid-3439_IUZBW.loom", mode = "r")
# First, let's check what's available
print("Available datasets:")
print(names(loom_data))

# Extract the matrices
spliced_raw <- loom_data[["layers/spliced"]]
unspliced_raw <- loom_data[["layers/unspliced"]]

# Get dimensions
print(paste("Spliced dimensions:", dim(spliced_raw)[1], "x", dim(spliced_raw)[2]))

# Extract full matrices
spliced_matrix <- spliced_raw[,]
unspliced_matrix <- unspliced_raw[,]

# Transpose (loom format is genes x cells, Seurat expects cells x genes for input)
spliced_matrix <- t(spliced_matrix)
unspliced_matrix <- t(unspliced_matrix)

# Get cell names if available
if("col_attrs" %in% names(loom_data) && "CellID" %in% names(loom_data[["col_attrs"]])) {
  cell_names <- loom_data[["col_attrs/CellID"]][]
} else {
  # Create cell names
  cell_names <- paste0("Cell_", 1:ncol(spliced_raw))
}

# Get gene names if available
if("row_attrs" %in% names(loom_data) && "Gene" %in% names(loom_data[["row_attrs"]])) {
  gene_names <- loom_data[["row_attrs/Gene"]][]
} else {
  # Create gene names
  gene_names <- paste0("Gene_", 1:nrow(spliced_raw))
}

# Close the connection
loom_data$close_all()

# Set proper row and column names
rownames(spliced_matrix)   <- gene_names
rownames(unspliced_matrix) <- gene_names

# Clean up column names
clean_names <- gsub("x", "", gsub("Combined_sid-3439_IUZBW:", "", cell_names))
clean_names <- paste('run4',clean_names,sep="_")

colnames(spliced_matrix)   <- clean_names
colnames(unspliced_matrix) <- clean_names

# colnames(seurat_obj) <- sub(".*_", "", colnames(seurat_obj))
genes = rownames(seurat_obj);length(genes)
cells = colnames(seurat_obj);length(cells)
common_cells <- intersect(cells, clean_names)
unspliced_matrix <- unspliced_matrix[genes, match(common_cells, clean_names)]
spliced_matrix   <- spliced_matrix[genes, match(common_cells, clean_names)]
# Now create the Seurat object
velocity_seurat_d4 <- CreateSeuratObject(counts = spliced_matrix, project = "velocity")
velocity_seurat_d4[["unspliced"]] <- CreateAssayObject(counts = unspliced_matrix)

# 1. Get intersection of cell names
common_cells <- intersect(colnames(seurat_obj), colnames(velocity_seurat_d4))

# 2. Define the order (using seurat_obj as reference)
ordered_cells <- colnames(seurat_obj)[colnames(seurat_obj) %in% common_cells]

# 3. Subset *and* reorder explicitly
seurat_obj_subset_d4 <- subset(seurat_obj, cells = ordered_cells)
velocity_seurat_d4 <- subset(velocity_seurat_d4, cells = ordered_cells)

# 4. Verify order
all(colnames(seurat_obj_subset_d4) == colnames(velocity_seurat_d4))
# Should return TRUE
all(rownames(seurat_obj_subset_d4)==rownames(velocity_seurat_d4))
# Add velocity assays to your main object
seurat_obj_subset_d4[["spliced"]] <- velocity_seurat_d4[["RNA"]]
seurat_obj_subset_d4[["unspliced"]] <- velocity_seurat_d4[["unspliced"]]

#### Merge velocities across runs
velocity_seurat_merged <- merge(velocity_seurat_d1, c(velocity_seurat_d2, velocity_seurat_d4), add.cell.ids = NULL)
print(velocity_seurat_merged)
names(velocity_seurat_merged@assays)

# Subset and reorder merged velocity object to match main Seurat object
common_cells <- intersect(colnames(seurat_obj), colnames(velocity_seurat_merged))
ordered_cells <- colnames(seurat_obj)[colnames(seurat_obj) %in% common_cells]
seurat_obj_subset <- subset(seurat_obj, cells = ordered_cells)
velocity_seurat_merged <- subset(velocity_seurat_merged, cells = ordered_cells)

# Verify matching order
all(colnames(seurat_obj_subset) == colnames(velocity_seurat_merged))  # TRUE
all(rownames(seurat_obj_subset) == rownames(velocity_seurat_merged))   # TRUE

# Add merged velocity assays to main object
seurat_obj_subset[["spliced"]] <- velocity_seurat_merged[["RNA"]]
seurat_obj_subset[["unspliced"]] <- velocity_seurat_merged[["unspliced"]]

# Plot UMAP of the integrated Seurat object
DimPlot_scCustom(seurat_obj_subset, reduction = "umap_harmony")

# Save and reload Seurat object
saveRDS(seurat_obj_subset, file = "seurat_obj_subset.rds")
seurat_obj_subset = readRDS("seurat_obj_subset.rds")
setwd("/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Cell/Velocity")

# Run SCVELO analysis
library(SCP)
options(Seurat.object.assay.version = "v3")
res = "SCT_snn_res.0.6"
reduction = "umap_harmony"
palcolor = c("#7F3C8D", "#E69F00", "#0072B2", "#009E73", "#8B0000")
dirpath = "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Cell/Velocity/Plots_white"
res = "SCT_snn_res.0.6"
seurat_velocity_computed <- RunSCVELO(
  srt = seurat_obj_subset, group_by = res,mode ="dynamical",magic_impute = TRUE,palcolor="#FFFFFF00",dpi = 1200,dirpath=dirpath,
  linear_reduction = "pca", nonlinear_reduction = reduction,n_pcs=30,calculate_velocity_genes=FALSE,save = TRUE)
saveRDS(seurat_velocity_computed, file = "seurat_velocity_computed.rds")
