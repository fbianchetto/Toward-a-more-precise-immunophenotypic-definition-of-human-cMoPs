##################################################################################################
### ATAC-seq Pipeline: Nucleosome-Free Regions (NFR) Differential Accessibility Analysis
###
### Pipeline overview:
###   - Input:  Consensus peaks from MACS3 (Chromatin accessibility profiling by ATAC-seq)
###   - Counts: Generated via createIterativeOverlapPeakSet.R + macs3_counts_table.sh
###   - Cells:  HSC_MPP, cMoP, NCP1
##################################################################################################

# ==============================================================================
# 0. SETUP: LIBRARIES, PATHS, HELPERS
# ==============================================================================

setwd("/home/patgen/working_dir/Data_analysis/ATAC_seq_RUN31/NFRs/")
source("/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/ATAC_seq_functions.R")

# --- Paths ---
path_plots   <- "/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/plots/"
path_objects <- "/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/"

# --- Libraries ---
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(readxl)
library(purrr)
library(cluster)
library(factoextra)

# Genomics
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(ChIPpeakAnno)
library(ChIPseeker)
library(biomaRt)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Differential accessibility
library(DESeq2)
library(pcaExplorer)

# Visualization
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(wesanderson)

# ==============================================================================
# 1. DATA IMPORT & PRE-PROCESSING
# ==============================================================================

# --- Raw counts ---
counts       <- read.delim("hsc_cMoPs_NCP1.paired", comment.char = "#")
counts$regions <- paste(counts$Chr, counts$Start, counts$End, sep = "-")

exclude <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "regions")
df <- counts[, !(colnames(counts) %in% exclude)]
rownames(df) <- counts$Geneid

# --- Sample metadata ---
cells   <- c("HSC_MPP", "cMoP", "NCP1")
coldata <- read.csv("SampleTable.csv", sep = ",") %>%
  mutate_if(is.character, as.factor) %>%
  mutate_if(is.integer, as.factor)

coldata$Cells <- gsub("HSC", "HSC_MPP", coldata$Cells)
coldata$Cells  <- factor(coldata$Cells, levels = cells)
rownames(coldata) <- coldata$Bam_files

# --- Align count matrix to metadata ---
df <- df[, rownames(coldata)]
stopifnot(all(rownames(coldata) == colnames(df)))

# ==============================================================================
# 2. GENOMIC COORDINATES: CANONICAL CHROMOSOMES ONLY
# ==============================================================================

prec <- readRDS("/home/patgen/working_dir/Data_analysis/ATAC_seq_RUN31/mergingPeaks/hsc_cMoPs_NCP1.fwp.filter.non_overlapping.rds")
prec_canonical <- keepStandardChromosomes(prec, pruning.mode = "tidy")

# ==============================================================================
# 3. DESeq2 OBJECT CONSTRUCTION
# ==============================================================================

dds <- DESeqDataSetFromMatrix(countData = df, colData = coldata, design = ~ Exp + Cells)

# Attach genomic metadata; rename rows and columns for clarity
meta        <- as.data.frame(prec_canonical)
meta$region <- paste(meta$seqnames, paste(meta$start, meta$end, sep = "-"), sep = ":")
mcols(dds)  <- prec_canonical
rownames(dds) <- meta$region

colnames(dds) <- gsub("_NFR.bam", "", colnames(dds))
colnames(dds) <- gsub("CD38",     "HSC_MPP", colnames(dds))

# ==============================================================================
# 4. FILTERING & DISPERSION ESTIMATION
# ==============================================================================

# Keep rows with >= 10 normalised counts in >= 2 samples
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 2
dds <- dds[idx, ]

# Compare dispersion fit types (diagnostic only)
par_fit <- estimateDispersions(dds, fitType = "parametric")
loc_fit <- estimateDispersions(dds, fitType = "local")
plotDispEsts(par_fit, main = "Dispersion: parametric")
plotDispEsts(loc_fit, main = "Dispersion: local")

# ==============================================================================
# 5. DIFFERENTIAL ACCESSIBILITY ANALYSIS
# ==============================================================================

# Wald test for pairwise comparisons
dds <- DESeq(dds, fitType = "local")

# LRT for ANOVA-style multi-group comparison (Cell type effect)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ Exp, fitType = "local")
res      <- results(dds_lrt)

# ==============================================================================
# 6. QUALITY CONTROL PLOTS
# ==============================================================================

# Cook's distance (outlier detection)
plot_cooks(dds, colData(dds))

# Variance-stabilised counts
vstMat     <- assay(varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local"))
notAllZero <- rowMeans(vstMat) > min(vstMat)

# RLE: normalization effectiveness across Raw / Norm / Filtered / VST counts
plots <- plot_countNorm(dds, colData(dds))
a <- plots[[1]]; b <- plots[[2]]; c <- plots[[3]]; d <- plots[[4]]; e <- plots[[5]]

filt_norm_counts <- dcast(plots[[6]][, c("Var2", "Var1", "value")], Var2 ~ Var1)
rownames(filt_norm_counts) <- filt_norm_counts$Var2
filt_norm_counts <- filt_norm_counts[, -1]

grid.arrange(
  a + theme(legend.position = "none"),
  b + theme(legend.position = "none"),
  d + theme(legend.position = "none"),
  c + theme(legend.position = "none"),
  bottom = e, ncol = 2, nrow = 2,
  left = textGrob("log2(counts)", gp = gpar(fontsize = 21, fontfamily = "Helvetica"), rot = 90)
)

rle_plots <- plot_countRLE(dds, sample_info = colData(dds), filt_norm_counts)
grid.arrange(
  rle_plots[[1]], rle_plots[[2]], rle_plots[[3]], rle_plots[[4]],
  ncol = 2, nrow = 2, bottom = rle_plots[[5]],
  left = textGrob("Relative Log Expression", rot = 90, gp = gpar(fontsize = 18))
)

# ==============================================================================
# 7. EXPLORATORY ANALYSIS: PCA & P-VALUE DISTRIBUTIONS
# ==============================================================================

sigvar_peaks <- rownames(res)[!is.na(res$padj) & res$padj < 0.01]
rld      <- rlog(dds_lrt, blind = TRUE, fitType = "local")
rld_sign <- rld[sigvar_peaks, ]

# PCA of significant regions
pcaplot(rld_sign, intgroup = "Cells")

# P-value histogram (should be uniform with peak near 0)
hist(res$pvalue[res$baseMean > 1], breaks = 30, col = "grey50", border = "white",
     main = "P-value Distribution", xlab = "p-value")

# Independent filtering: features passing adjusted p-value thresholds
use <- res$baseMean > metadata(res)$filterThreshold
h1  <- hist(res$padj[!use], breaks = 0:50/50, plot = FALSE)
h2  <- hist(res$padj[ use], breaks = 0:50/50, plot = FALSE)
colors_filt <- c(`do not pass` = "khaki", `pass` = "powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colors_filt, space = 0, ylab = "frequency")
text(x = c(0, length(h1$counts)), y = 0, label = c(0, 1), adj = c(0.5, 1.7), xpd = NA)
legend("topright", fill = rev(colors_filt), legend = rev(names(colors_filt)))

# ==============================================================================
# 8. SIGNIFICANT REGIONS (DACs): FILTERING & ANNOTATION
# ==============================================================================

# --- Filter: padj <= 0.01 → 24,131 regions ---
res_sig <- subset(res, padj <= 0.01)

hist(res_sig$padj[res_sig$baseMean > 1],
     breaks = seq(0, 0.01, 0.001), col = "grey50", border = "white",
     main = "P-adj Distribution (Significant Regions)", xlab = "Adjusted P-value")

de     <- rownames(res_sig)
select <- which(rownames(dds) %in% de)

# --- Attach LRT statistics to GRanges ---
DAC_regions <- mcols(dds[de, ])$X
DAC_regions$baseMean      <- res_sig$baseMean
DAC_regions$log2FoldChange <- res_sig$log2FoldChange
DAC_regions$padj           <- res_sig$padj
names(DAC_regions) <- rownames(res_sig)

# --- Annotate to nearest TSS using Ensembl Release 100 (GRCh38) ---
gtf       <- import("/home/patgen/working_dir/Data_analysis/ATAC_seq_RUN31/Homo_sapiens.GRCh38.100_chr.gtf")
genes_gtf <- gtf[gtf$type == "gene"]
genes_gtf$ensembl_gene_id <- genes_gtf$gene_id
names(genes_gtf) <- genes_gtf$gene_id

DAC_regions <- annotatePeakInBatch(
  DAC_regions,
  featureType            = "TSS",
  output                 = "nearestLocation",
  PeakLocForDistance     = "start",
  FeatureLocForDistance  = "TSS",
  multiple               = FALSE,
  AnnotationData         = genes_gtf
)

# Append HGNC symbols and biotypes via Ensembl April-2020 archive
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                host = "https://apr2020.archive.ensembl.org/")

DAC_regions <- addGeneIDs(
  annotatedPeak  = DAC_regions,
  mart           = mart,
  feature_id_type = "ensembl_gene_id",
  IDs2Add        = c("hgnc_symbol", "gene_biotype")
)

# Standardise names as chr:start-end; remove duplicates from annotation merge
names(DAC_regions) <- paste(
  seqnames(DAC_regions),
  paste(start(DAC_regions), end(DAC_regions), sep = "-"),
  sep = ":"
)

DAC_regions_no_duplicated <- DAC_regions[!duplicated(names(DAC_regions)), ]
saveRDS(DAC_regions_no_duplicated, file = paste0(path_objects, "DAC_regions_no_duplicated.rds"))
stopifnot(all(names(DAC_regions_no_duplicated) == rownames(res_sig)))

# ==============================================================================
# 9. PCA & BATCH CORRECTION
# ==============================================================================

# --- Pre-correction PCA ---
pcaplot(rld[de, ], intgroup = "Cells", ntop = length(de), ellipse = TRUE)
pcaplot(rld[de, ], intgroup = "Cells", ntop = length(de), ellipse = TRUE, pcX = 2, pcY = 3)

# --- Batch effect removal (limma) ---
# Remove Exp (batch) while preserving Cells (biological) signal
mod <- model.matrix(~ Cells, colData(rld))
mat <- limma::removeBatchEffect(assay(rld), batch = rld$Exp, design = mod)

rld_batch       <- rld
assay(rld_batch) <- mat

# --- Post-correction PCA ---
pcaplot(rld_batch[de, ], intgroup = "Cells", ntop = length(de), ellipse = FALSE) +
  ggtitle("Batch Corrected (PC1 vs PC2)", paste(length(de), "regions"))
pcaplot(rld_batch[de, ], intgroup = "Cells", ntop = length(de), ellipse = TRUE, pcX = 1, pcY = 3) +
  ggtitle("Batch Corrected (PC1 vs PC3)", paste(length(de), "regions"))

# --- Export PCA (Figure 4A) ---
pdf(file = paste0(path_plots, "Figure4A.pdf"), width = 10, height = 10)
pcaplot(rld_batch[de, ], intgroup = "Cells", ntop = length(de), ellipse = FALSE) +
  ggtitle("Batch Corrected PCA", paste(length(de), "Significant Regions"))
dev.off()

# ==============================================================================
# 10. HIERARCHICAL CLUSTERING (Figure 4B)
# ==============================================================================

suppressWarnings(
  d_hc <- stats::dist(t(assay(rld_batch)), method = "euclidean")
)
hc    <- hclust(d_hc, method = "ward.D2")
dend  <- as.dendrogram(hc)
dend2 <- dendextend::seriate_dendrogram(dend, d_hc, method = "OLO")

nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
par(mar = c(2, 1, 1, 3))
plot(dend2, xlab = "Distance (Height)", nodePar = nodePar, horiz = TRUE)

pdf(file = paste0(path_plots, "Figure4B.pdf"), width = 10, height = 15)
plot(dend2, xlab = "Distance (Height)", nodePar = nodePar, horiz = TRUE)
dev.off()

# ==============================================================================
# 11. K-MEANS CLUSTERING (Figure 4C)
# ==============================================================================

# --- Min-Max scaling per region (pattern-based clustering) ---
df_assay <- assay(rld_batch)[de, ]
df_assay <- t(apply(df_assay, 1, function(x) (x - min(x)) / (max(x) - min(x))))

target_samples <- c("HSC_MPP_SID23", "HSC_MPP_SID33", "HSC_MPP_SID41",
                    "cMop_SID22",    "cMop_SID32",    "cMop_SID40",
                    "NCP1_2_SID1",   "NCP1_2_SID11",  "NCP1_2_SID39")
x <- df_assay[, target_samples]

# --- K-means (k = 6) ---
set.seed(123)
km <- 6
res.kmDEG <- eclust(x, FUNcluster = "kmeans", nstart = 25, k = km, nboot = 500, seed = 123)
saveRDS(res.kmDEG, file = paste0(path_objects, "res.kmDEG.rds"))

# --- Silhouette validation ---
fviz_silhouette(res.kmDEG)
sil           <- res.kmDEG$silinfo$widths[, 1:3]
neg_sil_index <- which(sil[, "sil_width"] < 0)
message("Regions with negative silhouette: ", length(neg_sil_index))

# --- Cluster ordering via centroid dendrogram ---
res.kmDEG <- readRDS(file = paste0(path_objects, "res.kmDEG.rds"))
pheatmap::pheatmap(res.kmDEG$centers, cluster_cols = FALSE)

d_km  <- dist(res.kmDEG$centers, method = "euclidean")
hc_km <- hclust(d_km, method = "ward.D2")
dend_km  <- as.dendrogram(hc_km)
dend_km2 <- dendextend::seriate_dendrogram(dend_km, d_km, method = "OLO")
plot(dend_km2)

kmeans_levels <- c(6, 2, 3, 1, 4, 5)
split <- factor(res.kmDEG$cluster, levels = kmeans_levels)

# Cluster rename map (biological ordering m1–m6)
cluster_map <- c("6" = "m1", "2" = "m2", "3" = "m3", "1" = "m4", "4" = "m5", "5" = "m6")

# ==============================================================================
# 12. HEATMAP (Figure 4C)
# ==============================================================================

mat           <- as.matrix(res.kmDEG$data)
stopifnot(all(rownames(mat) == names(DAC_regions_no_duplicated)))
rownames(mat) <- paste(names(DAC_regions_no_duplicated), DAC_regions_no_duplicated$hgnc_symbol, sep = "_")

f1 <- colorRamp2(seq(0, 1, length = 3), c("blue", "white", "red"))
colors_set <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(6, "Set3"))

ha <- HeatmapAnnotation(
  celltype = c(rep("HSC_MPP", 3), rep("cMoP", 3), rep("NCP1", 3))
)

names(split) <- rownames(mat)

reorder.hmap <- Heatmap(
  mat, col = f1, name = "Relative Accessibility", split = split,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  height = unit(150, "mm"), width = unit(50, "mm"),
  column_names_gp = gpar(fontsize = 8), top_annotation = ha
)

cluster_hmap <- Heatmap(
  split,
  col             = structure(colors_set, names = as.character(1:km)),
  show_row_names  = FALSE, cluster_rows = TRUE,
  show_column_names = FALSE, show_heatmap_legend = FALSE,
  name = "cluster", width = unit(2, "mm"),
  show_row_dend = FALSE, split = split, height = unit(150, "mm")
)

# Add gene label annotations from curated list
df_genes        <- read_excel("/home/patgen/OneDrive/ATAC_seq_NCPs_RUN31/pdf_files/DAC_regions interenting genes.xlsx", sheet = 2)
interesting_loci <- paste(paste(df_genes$seqnames, paste(df_genes$start, df_genes$end, sep = "-"), sep = ":"),
                          df_genes$hgnc_symbol, sep = "_")
at   <- which(rownames(mat) %in% interesting_loci)
anno <- rowAnnotation(link = anno_mark(at = at, labels = rownames(mat)[at],
                                       labels_gp = gpar(fontsize = 6), padding = 1))

p_heatmap <- cluster_hmap + reorder.hmap + cluster_hmap + anno

pdf(file = paste0(path_plots, "Figure_4C.pdf"), width = 10, height = 10)
p_heatmap
dev.off()

# ==============================================================================
# 13. CLUSTER METADATA & BED FILE EXPORT
# ==============================================================================

DAC_regions_no_duplicated$cluster <- res.kmDEG$cluster

df_dac           <- as.data.frame(DAC_regions_no_duplicated)
df_dac$log2_score   <- log2(df_dac$score)
df_dac$log_padj     <- -log10(df_dac$padj)
df_dac$log_baseMean <- log2(df_dac$baseMean)
df_dac$cluster      <- factor(df_dac$cluster, levels = kmeans_levels)

# Standardise peak name prefixes
df_dac$name <- gsub("HSC",    "HSC/MPP_peak",  df_dac$name)
df_dac$name <- gsub("NCP1_2", "NCP1_2_peak",   df_dac$name)
df_dac$name <- gsub("cMoP",   "cMoP_peak",     df_dac$name)

selected_columns <- c("seqnames", "start", "end", "name", "score", "strand",
                      "insideFeature", "fromOverlappingOrNearest", "feature", "hgnc_symbol")

# Split peaks by cluster and export BED files
peaks_by_cluster <- setNames(
  lapply(1:6, function(i) subset(df_dac, cluster == i)[, selected_columns]),
  paste0("k", 1:6)
)

beds <- list(
  m1 = peaks_by_cluster$k6, m2 = peaks_by_cluster$k2, m3 = peaks_by_cluster$k3,
  m4 = peaks_by_cluster$k1, m5 = peaks_by_cluster$k4, m6 = peaks_by_cluster$k5,
  all_DACs = df_dac[, selected_columns]
)

setwd("/home/patgen/OneDrive/cMoP/Cell_Reports/R_script/DACRs")
for (n in names(beds)) {
  write.table(beds[[n]], file = paste0(n, "_annotated_peaks.bed"),
              row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
}

# ==============================================================================
# 14. ACCESSIBILITY BOXPLOTS BY CLUSTER (Figure S4A)
# ==============================================================================

mat_cluster           <- as.data.frame(cbind(res.kmDEG$data, cluster = res.kmDEG$cluster))
colnames(mat_cluster) <- c(rep("HSCs/MPPs", 3), rep("cMoPs", 3), rep("NCP1/2s", 3), "cluster")

df_mat         <- reshape2::melt(mat_cluster, id = "cluster", variable.name = "samples", value.name = "zscore")
df_mat$cluster <- cluster_map[as.character(df_mat$cluster)]
df_mat$cluster <- factor(df_mat$cluster, levels = paste0("m", 1:6))

o <- ggplot(df_mat, aes(x = samples, y = zscore, group = samples)) +
  geom_boxplot() +
  facet_grid(cols = vars(cluster), scales = "free")

cairo_pdf(file = paste0(path_plots, "Figure_S4A.pdf"), width = 10, height = 10)
print(o)
dev.off()

# ==============================================================================
# 15. GREAT ENRICHMENT ANALYSIS
# ==============================================================================

library(rGREAT)

# Helper: convert cluster data.frame to GRanges
df2grange <- function(df) {
  GRanges(
    seqnames = droplevels(df$seqnames),
    ranges   = IRanges(start = df$start, end = df$end),
    strand   = df$strand,
    name     = df$name,
    score    = df$score
  )
}

# Helper: filter GREAT results to significant, enriched terms
filtersGO <- function(GO) {
  subset(GO, Binom_Raw_PValue < 0.05 & Binom_Adjp_BH < 0.05 &
           Hyper_Adjp_BH < 0.05 & Binom_Fold_Enrichment > 2)
}

# Build GRanges list and submit GREAT jobs per cluster
cluster_grs  <- lapply(peaks_by_cluster, df2grange)
job_list     <- lapply(cluster_grs, function(gr) {
  submitGreatJob(gr, bg = NULL, species = "hg38", rule = "basalPlusExt",
                 request_interval = 10, version = "4.0")
})

# Retrieve and filter Biological Process results
job_list_BP <- lapply(names(job_list), function(i) {
  filtersGO(getEnrichmentTable(job_list[[i]], ontology = "GO Biological Process"))
})
names(job_list_BP) <- names(job_list)

# ==============================================================================
# 16. GO DOTPLOT (Figure S4_3C)
# ==============================================================================

top <- 10
top_combined <- map_dfr(names(job_list_BP), function(cl) {
  job_list_BP[[cl]][seq_len(top), ] %>%
    dplyr::filter(!is.na(name)) %>%
    mutate(cluster = cl, log_p = -log10(Binom_Raw_PValue))
})

top_combined$name <- factor(
  top_combined$name,
  levels = top_combined %>%
    group_by(name) %>%
    summarize(avg_log_p = mean(log_p), .groups = "drop") %>%
    arrange(avg_log_p) %>%
    pull(name)
)
top_combined$cluster <- cluster_map[as.character(top_combined$cluster)]
top_combined$cluster <- factor(top_combined$cluster, levels = paste0("m", 1:6))

BP_plot <- ggplot(top_combined, aes(x = cluster, y = name,
                                    color = log_p, size = Binom_Fold_Enrichment)) +
  geom_point() +
  labs(title  = "Top 10 Enriched GO Terms: Biological Process",
       x      = "Cluster", y = "GO Term",
       color  = expression(-log[10](p)),
       size   = "Binomial Fold Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), plot.title = element_text(hjust = 0.5))

pdf(file = paste0(path_plots, "FigureS2_3C.pdf"), width = 10, height = 10)
BP_plot
dev.off()

# ==============================================================================
# 17. GENE–REGION ASSOCIATIONS FROM GREAT
# ==============================================================================

# Extract associated genes per cluster, ordered by significance
GeneAssociations_list <- lapply(names(job_list), function(i) {
  ga <- getRegionGeneAssociations(job_list[[i]])
  names(ga) <- paste(seqnames(ga), paste(start(ga), end(ga), sep = "-"), sep = ":")
  idx        <- match(names(ga), names(DAC_regions_no_duplicated))
  ga$padj    <- DAC_regions_no_duplicated$padj[idx]
  ga[order(mcols(ga)$padj), ]
})
names(GeneAssociations_list) <- names(job_list)

# Pull gene lists per cluster
great_list <- lapply(GeneAssociations_list, function(ga) ga$annotated_genes@unlistData)

# Update gene symbols (cached to reduce API calls)
cache <- new.env(parent = emptyenv())
update_one_cached <- function(g) {
  if (exists(g, cache)) return(cache[[g]])
  res <- tryCatch(UpdateSymbolList(g), error = function(e) NA_character_)
  cache[[g]] <- res
  res
}

great_list_updated <- lapply(names(great_list), function(i) {
  message("Updating symbols: ", i)
  genes <- unique(as.character(great_list[[i]]))
  genes <- genes[!grepl("ENSG", genes) & !is.na(genes) & genes != ""]
  UpdateSymbolList(symbols = genes)
})
names(great_list_updated) <- names(great_list)

# ==============================================================================
# 18. SINGLE-CELL PROJECTION (Figure 4D)
# ==============================================================================

sc <- readRDS(paste0(path_objects, "NCPs_cMoPs.rds"))
DefaultAssay(sc) <- "RNA"
sc <- NormalizeData(sc)

# Harmonise cell-type labels
sc$Sample_Name <- gsub("NM",   "NCP",  sc$Sample_Name)
sc$Sample_Name <- gsub("cMOP", "cMoP", sc$Sample_Name)
sc$Sample_Name <- factor(sc$Sample_Name, levels = c("cMoP", "NCP1", "NCP2", "NCP3", "NCP4"))

# --- UMAP coloured by cell type (Figure 3D_1) ---
cell_colors <- c("#3B4CC0", "#B40426", "#1B9E77", "#E69F00", "#7F3C8D")
p_umap <- DimPlot(sc, reduction = "umap_learn50PC", label = TRUE,
                  group.by = "Sample_Name", cols = cell_colors)

pdf(file = paste0(path_plots, "Figure_4D_1.pdf"), width = 10, height = 10)
p_umap
dev.off()

# --- UCell module scores: project cluster-associated genes onto single cells ---
# Filter gene lists to genes present in the scRNA object
sc_genes <- rownames(sc)
great_list_sc <- lapply(great_list_updated, function(g) {
  g[g %in% sc_genes]
})
great_list_top500 <- lapply(great_list_sc, head, 500)
names(great_list_top500) <- names(great_list_updated)

sc <- UCell::AddModuleScore_UCell(sc, features = great_list_top500, name = "_cluster")

# Rename module score columns to m1–m6
sc_col_map <- c("great_k6_cluster" = "m1", "great_k2_cluster" = "m2",
                "great_k3_cluster" = "m3", "great_k1_cluster" = "m4",
                "great_k4_cluster" = "m5", "great_k5_cluster" = "m6")
cols <- colnames(sc@meta.data)
cols[cols %in% names(sc_col_map)] <- sc_col_map[cols[cols %in% names(sc_col_map)]]
colnames(sc@meta.data) <- cols

# --- Feature UMAPs for m4, m5, m6 (Figure 3D_2) ---
p4 <- SCP::FeatureDimPlot(sc, features = "m4", reduction = "umap_learn50PC",
                          theme_use = "theme_blank", ncol = 2, pt.size = 0.1, upper_cutoff = 0.26)
p5 <- SCP::FeatureDimPlot(sc, features = "m5", reduction = "umap_learn50PC",
                          theme_use = "theme_blank", ncol = 2, pt.size = 0.1, upper_cutoff = 0.26)
p6 <- SCP::FeatureDimPlot(sc, features = "m6", reduction = "umap_learn50PC",
                          theme_use = "theme_blank", ncol = 2, pt.size = 0.1, upper_cutoff = 0.285)

p_feature <- cowplot::plot_grid(p4, p5, p6, ncol = 2, align = "v")

cairo_pdf(file = paste0(path_plots, "Figure_4D_2.pdf"), width = 10, height = 10)
print(p_feature)
dev.off()
