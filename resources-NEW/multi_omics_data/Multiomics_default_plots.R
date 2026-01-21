# Default plots for Multiomics (generate manually and save to PNG)
# These plots are shown on the default page

# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

# load the data
fetal_heart <- readRDS("fetal_heart_multiomics_annotated.rds")

# UMAP demonstrate different cell populations based on RNA expression
p1 <- DimPlot(fetal_heart, reduction = "umap.rna", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("RNA") + labs(x = 'UMAP_1',y = 'UMAP_2')

# UMAP demonstrate different cell populations based on ATAC data
p2 <- DimPlot(fetal_heart, reduction = "umap.atac.peaks",label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC") + labs(x = 'UMAP_1',y = 'UMAP_2')

# dotplot demonstrate the marker expression of each cell population
DefaultAssay(fetal_heart) <-'RNA'
# Set levels - use tryCatch to handle if some levels don't exist
tryCatch({
  levels(fetal_heart) <- c('VCM','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','ACM')
}, error = function(e) {
  # If setting levels fails, just use existing levels
  cat("Note: Could not set all desired levels, using existing levels\n")
})
p3 <- DotPlot(fetal_heart,features = c("MYL2","MYH11","MYBPC1","HCN4","HBG1","STMN2","CD163","CD96","PDGFRA","PAX1","WT1","CDH5","NPPA")) + RotatedAxis() + labs(x = '',y = '')

# Combine plots: 2 on top, 1 on bottom (bottom plot same width as one top plot)
combined <- (p1 + p2) / (p3 + plot_spacer())

# Save to PNG in the same directory
ggsave("Multiomics_default_plots.png", plot = combined, width = 20, height = 15, device = "png")
cat("Saved default plots to Multiomics_default_plots.png\n")
