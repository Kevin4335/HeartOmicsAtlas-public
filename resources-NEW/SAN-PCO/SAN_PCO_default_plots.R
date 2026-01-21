# Default plots for SAN-PCO (generate manually and save to PNG)
# These plots are shown on the default page

# load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)

# load the data
SAN_PCO <- readRDS("SAN-PCO_annotated.rds")

levels(SAN_PCO) <- c("ACM","Endothelial","Epicardial","Epithelial","FB",
                     "Neuronal","Neural_Crest","Proliferating","SAN")

# UMAP demonstrate different cell populations based on RNA expression
# show in the default page
p1 <- UMAPPlot(SAN_PCO,label=TRUE,label.size = 4, repel = TRUE) + labs(x = "UMAP_1",y = "UMAP_2")

# dotplot demonstrate the marker expression of each cell population
# show in the default page
p2 <- DotPlot(SAN_PCO,features = c("NPPA","CDH5","WT1","EPCAM","COL1A1","STMN2","SOX2","TOP2A","SHOX2")) + RotatedAxis() + labs(x = '',y = '')

# Combine plots: side by side
combined <- p1 + p2

# Save to PNG in the same directory
ggsave("SAN_PCO_default_plots.png", plot = combined, width = 20, height = 8, device = "png")
cat("Saved default plots to SAN_PCO_default_plots.png\n")
