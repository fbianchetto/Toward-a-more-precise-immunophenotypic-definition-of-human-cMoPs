#!/usr/bin/env Rscript

# SoupX Decontamination for Multiple BD Rhapsody Samples
# This script removes ambient RNA contamination from single-cell data

# Load required libraries
library(SoupX)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
options(future.globals.maxSize = 2000 * 1024^2)
# Set base working directory
base_dir <- "/home/patgen/working_dir/Data_analysis/SingleCells_cMoP/Default_setting"

# Define samples to process
samples <- c("sid_3439")

# DoubletFinder parameters
doublet_rate <- 0.075  # 7.5% expected doublet rate (adjust based on your experiment)
pca_dims <- 1:30       # PCA dimensions to use

# ============================================================================
# DoubletFinder function
# ============================================================================

run_doubletfinder <- function(seurat_obj, doublet_rate = 0.075, dims = 1:30) {
  
  cat("Running DoubletFinder...\n")
  
  # Parameter sweep
  cat("  Performing parameter sweep...\n")
  sweep.res <- paramSweep(seurat_obj, PCs = dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Pick optimal pK
  pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  cat(sprintf("  Optimal pK: %.3f\n", pK))
  
  # Calculate expected number of doublets
  nExp <- round(doublet_rate * ncol(seurat_obj))
  cat(sprintf("  Expected doublets: %d (%.1f%%)\n", nExp, doublet_rate * 100))
  
  # Run DoubletFinder
  cat("  Running DoubletFinder classification...\n")
  seurat_obj <- doubletFinder(seurat_obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp, sct = TRUE)
  
  # Identify and store classification column
  df_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(df_col) == 1) {
    seurat_obj$DoubletFinder <- seurat_obj[[df_col]][,1]
    
    # Print summary
    doublet_summary <- table(seurat_obj$DoubletFinder)
    cat("\n  DoubletFinder Results:\n")
    print(doublet_summary)
    cat(sprintf("  Doublet rate: %.2f%%\n", 
                100 * doublet_summary["Doublet"] / sum(doublet_summary)))
  }
  
  return(seurat_obj)
}

# ============================================================================
# Function to process each sample
# ============================================================================

process_sample <- function(sample_name, base_dir) {
  
  cat("\n")
  cat("========================================================================\n")
  cat(sprintf("Processing sample: %s\n", sample_name))
  cat("========================================================================\n\n")
  
  # Set working directory for this sample
  sample_dir <- file.path(base_dir, sample_name)
  setwd(sample_dir)
  
  # Define file prefix (replace underscore with hyphen for filenames)
  file_prefix <- gsub("_", "-", sample_name)
  
  # ------------------------------------------------------------------------
  # 1. Load Data
  # ------------------------------------------------------------------------
  
  cat("Loading BD Rhapsody data...\n")
  
  # Define directories
  filtered_dir <- paste0(file_prefix, "_RSEC_MolsPerCell_MEX")
  unfiltered_dir <- paste0(file_prefix, "_RSEC_MolsPerCell_Unfiltered_MEX")
  seurat_file <- paste0(file_prefix, "_Seurat.rds")
  
  # Check if files exist
  if (!dir.exists(filtered_dir)) {
    cat(sprintf("ERROR: Filtered directory not found: %s\n", filtered_dir))
    return(NULL)
  }
  if (!dir.exists(unfiltered_dir)) {
    cat(sprintf("ERROR: Unfiltered directory not found: %s\n", unfiltered_dir))
    return(NULL)
  }
  if (!file.exists(seurat_file)) {
    cat(sprintf("ERROR: Seurat file not found: %s\n", seurat_file))
    return(NULL)
  }
  
  # Load matrices
  toc <- Read10X(filtered_dir)  # Filtered (cells)
  tod <- Read10X(unfiltered_dir)  # Unfiltered (cells + empty droplets)
  
  # Load existing Seurat object
  seurat_obj <- readRDS(seurat_file)
  
  cat("Updating Seurat object...\n")
  seurat_obj <- UpdateSeuratObject(seurat_obj)
  
  # Filter out Undetermined and Multiplet
  cat("Filtering out Undetermined and Multiplet cells...\n")
  seurat_obj <- subset(seurat_obj, 
                       subset = Sample_Name %in% c("Undetermined", "Multiplet"), 
                       invert = TRUE)
  
  cat(sprintf("Cells after filtering: %d\n", ncol(seurat_obj)))
  
  seurat_obj <- scCustomize::Add_Cell_QC_Metrics(
    seurat_obj, species = "human", add_cell_cycle = TRUE, 
    overwrite = TRUE, top_pct_name = "TopGene")
  seurat_obj <- subset(seurat_obj,subset = percent_mito < 25 )
  cat(sprintf("Cells after filtering by MT-: %d\n", ncol(seurat_obj)))
  # Process with SCTransform and dimensionality reduction
  cat("Running SCTransform and dimensionality reduction...\n")
  seurat_obj <- SCTransform(seurat_obj) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30)
  
  # ------------------------------------------------------------------------
  # 2.5 Run DoubletFinder (before SoupX)
  # ------------------------------------------------------------------------
  
  seurat_obj <- run_doubletfinder(seurat_obj, 
                                  doublet_rate = doublet_rate, 
                                  dims = pca_dims)
  
  # ------------------------------------------------------------------------
  # 2. Prepare SoupChannel Object
  # ------------------------------------------------------------------------
  
  cat("Creating SoupChannel object...\n")
  
  # Get cell barcodes from filtered Seurat object
  cells <- colnames(seurat_obj)
  
  # Find common genes between filtered and unfiltered matrices
  common_genes <- intersect(rownames(toc), rownames(tod))
  cat(sprintf("Common genes: %d\n", length(common_genes)))
  
  # Subset matrices to common genes and cells
  tod <- tod[common_genes, ]
  toc <- toc[common_genes, cells]
  
  # Verify gene order matches
  if (!all(rownames(tod) == rownames(toc))) {
    stop("Gene order mismatch between filtered and unfiltered matrices!")
  }
  
  # Check for and handle NAs in matrices
  cat("Checking for NA/NaN values...\n")
  
  # Check filtered matrix (toc)
  na_count_toc <- sum(is.na(toc@x))
  if (na_count_toc > 0) {
    cat(sprintf("  Warning: Found %d NA/NaN values in filtered matrix. Replacing with 0.\n", na_count_toc))
    toc@x[is.na(toc@x)] <- 0
  } else {
    cat("  Filtered matrix: No NA/NaN values found.\n")
  }
  
  # Check unfiltered matrix (tod)
  na_count_tod <- sum(is.na(tod@x))
  if (na_count_tod > 0) {
    cat(sprintf("  Warning: Found %d NA/NaN values in unfiltered matrix. Replacing with 0.\n", na_count_tod))
    tod@x[is.na(tod@x)] <- 0
  } else {
    cat("  Unfiltered matrix: No NA/NaN values found.\n")
  }
  
  # Check for infinite values
  inf_count_toc <- sum(is.infinite(toc@x))
  inf_count_tod <- sum(is.infinite(tod@x))
  
  if (inf_count_toc > 0) {
    cat(sprintf("  Warning: Found %d infinite values in filtered matrix. Replacing with 0.\n", inf_count_toc))
    toc@x[is.infinite(toc@x)] <- 0
  }
  
  if (inf_count_tod > 0) {
    cat(sprintf("  Warning: Found %d infinite values in unfiltered matrix. Replacing with 0.\n", inf_count_tod))
    tod@x[is.infinite(tod@x)] <- 0
  }
  
  # Verify matrices are valid sparse matrices
  if (!inherits(toc, "dgCMatrix") || !inherits(tod, "dgCMatrix")) {
    cat("  Converting matrices to dgCMatrix format...\n")
    toc <- as(toc, "dgCMatrix")
    tod <- as(tod, "dgCMatrix")
  }
  
  cat("  Matrix validation complete.\n")
  
  # Create SoupChannel
  sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
  
  # ------------------------------------------------------------------------
  # 3. Add Clustering Information
  # ------------------------------------------------------------------------
  
  cat("Adding clustering information...\n")
  cat("Available metadata columns:\n")
  print(colnames(seurat_obj@meta.data))
  
  # Find cluster information
  clusters <- NULL
  cluster_candidates <- c("seurat_clusters", "SCT_snn_res.0.8", "RNA_snn_res.1", 
                          "RNA_snn_res.0.5", "wsnn_res.0.8", "Cell_Type_Experimental")
  
  for (col in cluster_candidates) {
    if (col %in% colnames(seurat_obj@meta.data)) {
      clusters <- seurat_obj@meta.data[[col]]
      cat(sprintf("Using clustering from: %s\n", col))
      break
    }
  }
  
  # Fallback: search for any cluster column
  if (is.null(clusters)) {
    cluster_cols <- grep("cluster|res\\.|type", 
                         colnames(seurat_obj@meta.data), 
                         value = TRUE, ignore.case = TRUE)
    if (length(cluster_cols) > 0) {
      clusters <- seurat_obj@meta.data[[cluster_cols[1]]]
      cat(sprintf("Using clustering from: %s\n", cluster_cols[1]))
    }
  }
  
  # Last resort: use UMAP k-means
  if (is.null(clusters)) {
    cat("No existing clusters found. Creating clusters with k-means...\n")
    umap_embed <- Embeddings(seurat_obj, "umap")
    clusters <- kmeans(umap_embed, centers = 10)$cluster
    cat("Created 10 clusters using k-means on UMAP coordinates\n")
  }
  
  # Ensure clusters match cell barcodes
  if (length(clusters) != ncol(toc)) {
    cell_names <- colnames(toc)
    seurat_names <- colnames(seurat_obj)
    matching_idx <- match(cell_names, seurat_names)
    clusters <- clusters[matching_idx]
  }
  
  # Add clusters to SoupChannel
  sc <- setClusters(sc, setNames(as.character(clusters), colnames(toc)))
  
  # Add UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, "umap")
  sc <- setDR(sc, umap_coords[colnames(toc), ])
  
  # ------------------------------------------------------------------------
  # 4. Estimate Contamination
  # ------------------------------------------------------------------------
  
  cat("Estimating contamination fraction...\n")
  
  # Automatic estimation
  sc <- autoEstCont(sc, doPlot = TRUE)
  
  # Save contamination plot
  ggsave(paste0(file_prefix, "_SoupX_contamination_estimation.pdf"), 
         width = 10, height = 8)
  
  cat(sprintf("Estimated contamination fraction: %.3f\n", 
              sc$metaData$rho[1]))
  
  # ------------------------------------------------------------------------
  # 5. Visualize Soup Profile
  # ------------------------------------------------------------------------
  
  soup_profile <- data.frame(
    gene = rownames(sc$soupProfile),
    expression = sc$soupProfile[, 1]
  )
  soup_profile <- soup_profile[order(-soup_profile$expression), ]
  
  cat("\nTop 20 genes in ambient RNA profile:\n")
  print(head(soup_profile, 20))
  
  # Plot top contaminating genes (excluding mitochondrial genes)
  # Filter out MT genes for better visualization
  non_mt_genes <- soup_profile$gene[!grepl("^MT-", soup_profile$gene)]
  
  if (length(non_mt_genes) >= 3) {
    top_genes <- non_mt_genes[1:3]
  } else {
    top_genes <- soup_profile$gene[1:min(3, nrow(soup_profile))]
  }
  
  # Check if genes exist in the count matrix
  valid_genes <- top_genes[top_genes %in% rownames(sc$toc)]
  
  if (length(valid_genes) > 0) {
    pdf(paste0(file_prefix, "_SoupX_marker_maps.pdf"), width = 12, height = 10)
    tryCatch({
      plotMarkerMap(sc, valid_genes)
    }, error = function(e) {
      cat(sprintf("Warning: Could not create marker maps: %s\n", e$message))
      plot.new()
      text(0.5, 0.5, "Marker map generation failed", cex = 1.5)
    })
    dev.off()
  } else {
    cat("Warning: No valid genes found for marker map plotting\n")
  }
  
  # ------------------------------------------------------------------------
  # 6. Adjust Counts
  # ------------------------------------------------------------------------
  
  cat("\nAdjusting counts to remove contamination...\n")
  out <- adjustCounts(sc, roundToInt = TRUE)
  
  # ------------------------------------------------------------------------
  # 7. Save Corrected Data
  # ------------------------------------------------------------------------
  
  cat("Saving decontaminated data...\n")
  
  # Create output directory
  output_dir <- paste0(file_prefix, "_SoupX_corrected")
  dir.create(output_dir, showWarnings = FALSE)
  
  # Save as MEX format
  writeMM(out, file.path(output_dir, "matrix.mtx"))
  writeLines(colnames(out), file.path(output_dir, "barcodes.tsv"))
  write.table(data.frame(gene_id = rownames(out), 
                         gene_name = rownames(out)),
              file.path(output_dir, "features.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Compress files
  system(sprintf("gzip %s/matrix.mtx", output_dir))
  system(sprintf("gzip %s/barcodes.tsv", output_dir))
  system(sprintf("gzip %s/features.tsv", output_dir))
  
  # ------------------------------------------------------------------------
  # 8. Create Corrected Seurat Object
  # ------------------------------------------------------------------------
  
  cat("Creating corrected Seurat object...\n")
  
  seurat_corrected <- CreateSeuratObject(
    counts = out,
    project = paste0(sample_name, "_SoupX"),
    meta.data = seurat_obj@meta.data
  )
  
  # seurat_corrected <- CreateSeuratObject(
  #   counts = out,
  #   project = paste0(sample_name, "_SoupX")
  # )
  
  # Add DoubletFinder results visualization
  if ("DoubletFinder" %in% colnames(seurat_obj@meta.data)) {
    cat("Creating DoubletFinder visualization...\n")
    
    # Transfer UMAP coordinates for visualization
    seurat_corrected[["umap"]] <- CreateDimReducObject(
      embeddings = Embeddings(seurat_obj, "umap"),
      key = "UMAP_",
      assay = DefaultAssay(seurat_corrected)
    )
    
    # Create DoubletFinder plot
    p1 <- DimPlot(seurat_corrected, reduction = "umap", group.by = "DoubletFinder",
                  cols = c("Singlet" = "blue", "Doublet" = "red")) +
      ggtitle(paste0(sample_name, " - DoubletFinder Classification"))
    
    ggsave(paste0(file_prefix, "_DoubletFinder_UMAP.pdf"), 
           plot = p1, width = 10, height = 8)
  }
  
  # Save corrected Seurat object (with DoubletFinder annotations)
  saveRDS(seurat_corrected, paste0(file_prefix, "_Seurat_SoupX_corrected.rds"))
  
  # ------------------------------------------------------------------------
  # 9. Summary Statistics
  # ------------------------------------------------------------------------
  
  total_before <- sum(toc)
  total_after <- sum(out)
  reduction <- (total_before - total_after) / total_before * 100
  
  # DoubletFinder stats
  if ("DoubletFinder" %in% colnames(seurat_obj@meta.data)) {
    doublet_counts <- table(seurat_obj$DoubletFinder)
    n_doublets <- as.numeric(doublet_counts["Doublet"])
    n_singlets <- as.numeric(doublet_counts["Singlet"])
    doublet_pct <- 100 * n_doublets / (n_doublets + n_singlets)
  } else {
    n_doublets <- NA
    n_singlets <- NA
    doublet_pct <- NA
  }
  
  summary_stats <- data.frame(
    Sample = sample_name,
    Metric = c("Total UMIs before SoupX", 
               "Total UMIs after SoupX", 
               "UMIs removed (%)",
               "Mean contamination fraction", 
               "Number of cells",
               "Singlets (DoubletFinder)",
               "Doublets (DoubletFinder)",
               "Doublet rate (%)"),
    Value = c(total_before, 
              total_after, 
              sprintf("%.2f%%", reduction),
              sprintf("%.3f", mean(sc$metaData$rho)),
              ncol(out),
              n_singlets,
              n_doublets,
              sprintf("%.2f%%", doublet_pct))
  )
  
  write.csv(summary_stats, 
            paste0(file_prefix, "_SoupX_DoubletFinder_summary.csv"), 
            row.names = FALSE)
  
  cat("\n")
  print(summary_stats)
  
  cat(sprintf("\n=== %s processing complete! ===\n", sample_name))
  
  return(summary_stats)
}

# ============================================================================
# Main execution: Process all samples
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("SoupX Batch Processing for BD Rhapsody Samples\n")
cat("========================================================================\n")
cat(sprintf("Samples to process: %s\n", paste(samples, collapse = ", ")))
cat(sprintf("Base directory: %s\n", base_dir))
cat("\n")

# Store results
all_results <- list()

# Process each sample
for (sample in samples) {
  tryCatch({
    result <- process_sample(sample, base_dir)
    all_results[[sample]] <- result
  }, error = function(e) {
    cat(sprintf("\nERROR processing %s: %s\n", sample, e$message))
    all_results[[sample]] <- NULL
  })
}

# ============================================================================
# Combined Summary
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("COMBINED SUMMARY\n")
cat("========================================================================\n\n")

# Combine all results
if (length(all_results) > 0) {
  combined_summary <- do.call(rbind, all_results)
  print(combined_summary)
  
  # Save combined summary
  setwd(base_dir)
  write.csv(combined_summary, "All_Samples_SoupX_Summary.csv", row.names = FALSE)
  cat("\nCombined summary saved to: All_Samples_SoupX_Summary.csv\n")
}

cat("\n=== All samples processed! ===\n")

