# ==============================================================================
# HELPER FUNCTIONS FOR ATAC-seq DATA VISUALIZATION
# ==============================================================================
library(wesanderson)
pal <- wes_palette(50, name = "Zissou1", type = "continuous")
#' Custom Notebook Theme for GGPlot2
#' @description Sets a consistent, high-resolution aesthetic for plots in notebooks.
#' @param base_size Font size for text elements.
#' @param base_family Font family for text elements.
theme_notebook <- function(base_size = 18, base_family = "helvetica") {
  theme_set(theme_minimal(base_size = base_size, base_family = base_family)) +
    theme(
      plot.title     = element_text(face = "bold", size = 20, hjust = 0.5),
      axis.title     = element_text(face = "bold", size = rel(1)),
      axis.title.y   = element_text(angle = 90, vjust = 2, size = 20),
      axis.title.x   = element_text(vjust = -0.2, size = 20),
      axis.text      = element_text(size = 20),
      axis.line      = element_line(colour = "black"),
      axis.ticks     = element_line(),
      legend.key     = element_rect(colour = NA),
      legend.key.size = unit(0.5, "cm"),
      legend.margin  = unit(0.5, "cm"),
      legend.text    = element_text(size = 14),
      legend.title   = element_text(size = 16)
    )
}

# Define global palette
pal <- wes_palette(50, name = "Zissou1", type = "continuous")

# ------------------------------------------------------------------------------
# QC PLOTS
# ------------------------------------------------------------------------------

#' Plot Cook's Distance
#' @description Generates boxplots of Cook's distance to identify potential outliers in DESeq2 data.
#' @param dds A DESeqDataSet object.
#' @param sample_info Dataframe containing metadata for samples.
plot_cooks <- function(dds, sample_info) {
  require(reshape2)
  
  # Extract Cook's distance and transform to log10
  cooks <- melt(as.data.frame(log10(assays(dds)[["cooks"]])))
  cooks <- merge(cooks, sample_info, by.x = "variable", by.y = 0, all = T)
  
  p <- ggplot(cooks, aes(y = value, x = variable, colour = Sample)) +
    geom_boxplot() +
    theme_notebook() +
    theme(axis.text.x = element_blank()) +
    labs(y = "Cooks Distance (log10)", x = "Samples") +
    scale_colour_manual(values = pal)
  
  return(p)
}

#' Relative Log Expression (RLE) Comparison
#' @description Compares Raw, Normalized, Filtered, and VST counts using RLE plots.
#' @param dds DESeqDataSet object.
#' @param sample_info Metadata for samples.
#' @param filt_norm_counts Matrix of filtered normalized counts.
plot_countRLE <- function(dds, sample_info, filt_norm_counts) {
  
  # Internal helper to calculate RLE dataframes
  calc_rle <- function(count_mat, meta) {
    medians <- apply(count_mat, 1, median)
    rle <- log2(count_mat / medians)
    df <- reshape2::melt(t(rle))
    merge(df, meta, by.x = "Var1", by.y = 0, all = T)
  }
  
  # 1. Raw counts
  p1 <- ggplot(calc_rle(counts(dds), sample_info), aes(y = value, x = Var1, colour = Cells)) +
    geom_boxplot() + theme_notebook() + theme(axis.text.x = element_blank()) +
    geom_hline(yintercept = 0, lty = "dashed") + labs(title = "Raw counts", y = "RLE", x = "")
  
  # 2. Normalized counts
  p2 <- ggplot(calc_rle(counts(dds, normalized = TRUE), sample_info), aes(y = value, x = Var1, colour = Cells)) +
    geom_boxplot() + theme_notebook() + theme(axis.text.x = element_blank()) +
    geom_hline(yintercept = 0, lty = "dashed") + labs(title = "Normalized counts", y = "", x = "")
  
  # 3. Filtered normalized counts
  p3 <- ggplot(calc_rle(filt_norm_counts, sample_info), aes(y = value, x = Var1, colour = Cells)) +
    geom_boxplot() + theme_notebook() + theme(axis.text.x = element_blank()) +
    geom_hline(yintercept = 0, lty = "dashed") + labs(title = "Filtered counts", y = "", x = "")
  
  # 4. VST counts
  p4 <- ggplot(calc_rle(vstMat, sample_info), aes(y = value, x = Var1, colour = Cells)) +
    geom_boxplot() + theme_notebook() + theme(axis.text.x = element_blank()) +
    geom_hline(yintercept = 0, lty = "dashed") + labs(title = "VST counts", y = "", x = "")
  
  # Legend Extraction
  get_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    return(tmp$grobs[[leg]])
  }
  
  key <- get_legend(p4 + theme(legend.direction = "horizontal"))
  return(list(p1, p2, p3, p4, key))
}

library(reshape2)
plot_countNorm <- function(dds, sample_info){
  
  # raw counts
  raw <- reshape2::melt(t(counts(dds)))
  raw <- merge(raw, sample_info, by.x="Var1", by.y=0, all=T)
  Palette <- wes_palette("Zissou1", 50, type = "continuous")
  p1 <- ggplot(raw, aes(y=log2(value), x=Var1, colour=Cells)) + 
    geom_boxplot() + 
    theme_notebook() +
    theme(axis.text.x=element_blank()) +
    labs(title="Raw counts", x="", y="") + 
    scale_colour_manual(values=Palette)
  
  # norm counts = counts(dds) / sizeFactors(dds)
  norm <- reshape2::melt(t(counts(dds, normalized=T)))
  norm <- merge(norm, sample_info, by.x="Var1", by.y=0, all=T)
  
  p2 <- ggplot(norm, aes(y=log2(value), x=Var1, colour=Cells)) + 
    geom_boxplot() + 
    theme_notebook() +
    theme(axis.text.x=element_blank()) +
    labs(title="Normalised counts", y="", x="") + 
    scale_colour_manual(values=Palette)
  
  # VST counts
  colnames(vstMat) <- dds@colData@rownames
  data <-reshape2::melt(vstMat[notAllZero,])
  data <- merge(data, sample_info, by.x="Var2", by.y=0, all=T)
  
  p3 <- ggplot(data, aes(y=value, x=Var2, colour=Cells)) + 
    geom_boxplot()  + 
    theme_notebook() +
    theme(axis.text.x=element_blank()) +
    labs(title="VST counts", y="", x="") + 
    scale_colour_manual(values=Palette) 
  
  # normalised counts with excluded counts (padj == NA) removed
  # counts excluded by: Cooks cutoff & independent filtering
  filtered_genes <- subset(res, res@listData$padj != "NA")@"rownames" # get sig DE gene names
  filt_norm_counts <- subset(reshape2::melt(t(counts(dds,normalized=T))), Var2 %in% filtered_genes)
  
  filt_norm_counts_info <- merge(filt_norm_counts, sample_info, by.x="Var1", by.y=0, all=T)
  
  p4 <- ggplot(filt_norm_counts_info, aes(y=log2(value), x=Var1, colour=Cells)) + 
    geom_boxplot() + 
    theme_notebook() +
    theme(axis.text.x=element_blank()) +
    labs(title="Filtered normalised counts", x="", y="") + 
    scale_colour_manual(values=Palette) +
    theme(legend.direction="horizontal")
  
  get_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)}
  
  key <- get_legend(p4)
  
  plots = list(p1, p2, p3, p4, key, filt_norm_counts)
  
  return(plots)
}

# ------------------------------------------------------------------------------
# ENRICHMENT ANALYSIS PLOTS
# ------------------------------------------------------------------------------

#' Custom DotPlot for ClusterProfiler
#' @description Creates a refined DotPlot for enrichment results with custom styling.
#' @param clusterProfilerOutput Result object from enrichGO or compareCluster.
#' @param title String for plot title.
#' @param showCategory Number of top categories to display per cluster.
MyDotPlot <- function(clusterProfilerOutput = NULL, title = NULL, showCategory = 10, decreasing = FALSE) {
  require(dplyr)
  
  df <- fortify(clusterProfilerOutput)
  df$Log <- -log(df$p.adjust, 10)
  
  # Select top categories per cluster
  top_df <- df %>% group_by(Cluster) %>% top_n(n = showCategory, wt = Log)
  
  plot_theme <- theme(
    axis.title.x = element_text(size = 8),
    axis.text.x  = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.y  = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title   = element_text(hjust = 0.5, size = 12),
    legend.key.size = unit(10, 'mm'),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8)
  )
  
  q <- ggplot(top_df, aes(x = Cluster, y = Description, size = GeneRatio, color = p.adjust)) +
    geom_point() +
    scale_color_continuous(low = "red", high = "blue", name = "p.adjust", guide = guide_colorbar(reverse = TRUE)) +
    plot_theme + xlab(NULL) + ylab(NULL) + ggtitle(title)
  
  return(q)
}