library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)

# === Set data directory ===
data_dir <- "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered/output_mtx/annotated_barcodes"  
setwd(data_dir)

# === Step 1: Find all matrix files ===
matrix_files <- list.files(pattern = "_matrix\\.mtx\\.gz$")

# === Step 2: Initialize list for Seurat objects ===
seurat_list <- list()

for (file in matrix_files) {
  # Extract sample name prefix
  sample_name <- str_remove(file, "_matrix\\.mtx\\.gz$")
  
  message("🔄 Processing: ", sample_name)
  
  # Define file paths
  matrix_path   <- gzfile(paste0(sample_name, "_matrix.mtx.gz"), "rt")
  barcodes_path <- paste0(sample_name, "_barcodes.csv.gz")
  if (!file.exists(barcodes_path)) {
    barcodes_path <- paste0(sample_name, "_barcodes.tsv.gz")
  }
  features_path <- gzfile(paste0(sample_name, "_features.tsv.gz"), "rt")
  
  # === Read matrix ===
  counts <- readMM(matrix_path)
  
  # === Read barcodes ===
  barcodes <- read.delim(gzfile(barcodes_path), header = TRUE, stringsAsFactors = FALSE)
  colnames(counts) <- barcodes$V1
  
  # === Read features ===
  features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- features$V1
  
  # === Create Seurat object ===
  seu <- CreateSeuratObject(counts = counts, project = sample_name)
  
  # Add sample ID metadata
  seu$orig.ident <- sample_name
  seu$Sample_Name <- sample_name
  
  # Add to list
  seurat_list[[sample_name]] <- seu
}

# === Step 3: Merge all into one Seurat object ===
merged_seurat <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "Merged_cMoP"
)

# === Optional: Save merged object ===
saveRDS(merged_seurat, file = "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/R_rhapsody2.3/Soupx_filtered/merged_cMoP_seurat.rds")

# Done!
print("✅ Merging complete.")
print(merged_seurat)
