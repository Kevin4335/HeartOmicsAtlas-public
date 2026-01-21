# Default plots for ACM_VCM_SAN (generate manually and save to PNG)
# These plots are shown on the default page

# load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)

# load the data
Combine <- readRDS("ACM_VCM_SAN_annotated.rds")

levels(Combine) <- c("ACM","Endothelial","Epicardial","Epithelial","FB_1","FB_2",
                     "Neuronal","Neural_Crest","Proliferating","SAN","VCM")

# UMAP demonstrate different cell populations based on RNA expression
# show in the default page
p1 <- UMAPPlot(Combine,label=TRUE,label.size =5,repel = TRUE) + labs(x = 'UMAP_1', y = 'UMAP_2')

p2 <- UMAPPlot(Combine,label=TRUE,split.by = "orig.ident",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

# dotplot demonstrate the marker expression of each cell population
# show in the default page
p3 <- DotPlot(Combine,features = c("NPPA","CDH5","WT1","EPCAM","POSTN","PDZRN4","STMN2","SOX2","TOP2A","SHOX2","MYL2")) + RotatedAxis()

# Combine plots: 2 on top, 1 on bottom (bottom plot same width as one top plot)
combined <- (p1 + p2) / (p3 + plot_spacer())

# Save to PNG in the same directory
ggsave("ACM_VCM_SAN_default_plots.png", plot = combined, width = 20, height = 15, device = "png")
cat("Saved default plots to ACM_VCM_SAN_default_plots.png\n")
