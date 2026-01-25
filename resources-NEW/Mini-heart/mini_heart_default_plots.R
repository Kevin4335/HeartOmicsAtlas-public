# Default plots for Mini-heart (generate manually and save to PNG)
# These plots are shown on the default page

# load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)

# load the data
heart <- readRDS("mini_heart_annotated.rds")
heart <- RenameIdents(heart, 'FB' = 'Fibroblast')
heart$celltype <- Idents(heart)
levels(heart) <- c("ACM","Endothelial","Epicardial","Epithelial","Fibroblast",
                   "Neuronal","Proliferating","SAN","VCM")

# UMAP demonstrate different cell populations based on RNA expression
# first plot to demonstrate in default page
p1 <- UMAPPlot(heart,label=TRUE,label.size =6,repel = TRUE) + labs(x = 'UMAP_1',y = 'UMAP_2') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# dotplot demonstrate the marker expression of each cell population
# second plot to demonstrate in default page
p2 <- DotPlot(heart,features = c("NPPA","CDH5","WT1","EPCAM","DCN","STMN2","TOP2A","SHOX2","HEY2")) +
  labs(x='',y='') +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  RotatedAxis()

# Combine plots: side by side
combined <- p1 + p2

# Save to PNG in the same directory
ggsave("mini_heart_default_plots.png", plot = combined, width = 20, height = 8, device = "png")
cat("Saved default plots to mini_heart_default_plots.png\n")
