# Load necessary libraries
library(Seurat)      # For creating and handling Seurat objects
library(Matrix)      # For working with sparse matrices
library(dplyr)       # For data manipulation
library(stringr)     # For string operations like pattern matching

# === Set working directory containing filtered matrix output ===
data_dir <- "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered/output_mtx/annotated_barcodes"  
setwd(data_dir)

# === Step 1: Identify all matrix files ===
# Find all .mtx.gz files (one per sample)
matrix_files <- list.files(pattern = "_matrix\\.mtx\\.gz$")

# === Step 2: Create a list to hold Seurat objects ===
seurat_list <- list()

# Loop over each matrix file and process it
for (file in matrix_files) {
  # Extract sample name by removing suffix from filename
  sample_name <- str_remove(file, "_matrix\\.mtx\\.gz$")
  
  message("🔄 Processing: ", sample_name)  # Log progress
  
  # === Define input file paths ===
  matrix_path   <- gzfile(paste0(sample_name, "_matrix.mtx.gz"), "rt")         # Expression matrix
  barcodes_path <- paste0(sample_name, "_barcodes.csv.gz")                     # Try CSV first
  if (!file.exists(barcodes_path)) {
    barcodes_path <- paste0(sample_name, "_barcodes.tsv.gz")                   # Fallback to TSV
  }
  features_path <- gzfile(paste0(sample_name, "_features.tsv.gz"), "rt")       # Feature/gene info
  
  # === Load sparse count matrix ===
  counts <- readMM(matrix_path)  # Matrix Market format
  
  # === Load barcodes (cell IDs) ===
  barcodes <- read.delim(gzfile(barcodes_path), header = TRUE, stringsAsFactors = FALSE)
  colnames(counts) <- barcodes$V1  # Assign barcodes to columns of matrix
  
  # === Load features (gene names) ===
  features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- features$V1  # Assign gene IDs to rows of matrix
  
  # === Create Seurat object for the current sample ===
  seu <- CreateSeuratObject(counts = counts, project = sample_name)
  
  # Annotate metadata with sample ID
  seu$orig.ident <- sample_name
  seu$Sample_Name <- sample_name
  
  # Store in the list
  seurat_list[[sample_name]] <- seu
}

# === Step 3: Merge all Seurat objects into one ===
merged_seurat <- merge(
  seurat_list[[1]],                # First object
  y = seurat_list[-1],            # All other objects
  add.cell.ids = names(seurat_list),  # Prefix cell names with sample ID
  project = "Merged_cMoP"
)

# === Step 4: Save merged object to disk ===
saveRDS(
  merged_seurat, 
  file = "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered/merged_cMoP_seurat.rds"
)

# === Final message ===
print("✅ Merging complete.")
print(merged_seurat)
