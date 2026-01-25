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
SAN_PCO <- RenameIdents(SAN_PCO, 'FB' = 'Fibroblast')
SAN_PCO$celltype <- Idents(SAN_PCO)
levels(SAN_PCO) <- c("ACM","Endothelial","Epicardial","Epithelial","Fibroblast",
                     "Neuronal","Neural_Crest","Proliferating","SAN")

# UMAP demonstrate different cell populations based on RNA expression
# first plot to demonstrate in default page
p1 <- UMAPPlot(SAN_PCO,label=TRUE,label.size = 6, repel = TRUE) + labs(x = 'UMAP_1',y = 'UMAP_2') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# dotplot demonstrate the marker expression of each cell population
# second plot to demonstrate in default page
p2 <- DotPlot(SAN_PCO,features = c("NPPA","CDH5","WT1","EPCAM","COL1A1","STMN2","SOX2","TOP2A","SHOX2")) +
  labs(x='',y='') +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  RotatedAxis()

# Combine plots: side by side
combined <- p1 + p2

# Save to PNG in the same directory
ggsave("SAN_PCO_default_plots.png", plot = combined, width = 20, height = 8, device = "png")
cat("Saved default plots to SAN_PCO_default_plots.png\n")
