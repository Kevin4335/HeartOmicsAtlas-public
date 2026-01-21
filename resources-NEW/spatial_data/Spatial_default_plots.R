# Default plots for Spatial (generate manually and save to PNG)
# These plots are shown on the default page

# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

# load the data
fetal_heart <- readRDS("fetal_heart_0103_16um_annotated.rds")

# Update Seurat object to ensure compatibility with spatial plotting
fetal_heart <- UpdateSeuratObject(fetal_heart)

# UMAP demonstrate different cell populations
p1 <- DimPlot(fetal_heart,reduction = "umap.016um",label = TRUE, repel = TRUE) + labs(x = "UMAP1",y = "UMAP2")

# Spatial plot demonstrate different cell populations
mycols <- c(
  "orange",
  "darkgreen",
  "#F781BF", 
  "#FF7F00", 
  "#984EA3",
  "#AEC6CF",
  "#E41A1C",
  "#A65628" 
)
p2 <- SpatialDimPlot(fetal_heart,label = TRUE,label.size = 3,repel = TRUE) + scale_fill_manual(values = mycols) + theme_void()

# Spatial plot to show the marker expression of each cluster using Dotplot
p3 <- DotPlot(fetal_heart,features = c("MYL7","VWF","DCN","SPP1","PLP1","HBG1","SHOX2","MYH11")) + RotatedAxis() + labs(x = '',y='')

# Combine plots: 2 on top, 1 on bottom (bottom plot same width as one top plot)
combined <- (p1 + p2) / (p3 + plot_spacer())

# Save to PNG in the same directory
ggsave("Spatial_default_plots.png", plot = combined, width = 20, height = 15, device = "png")
cat("Saved default plots to Spatial_default_plots.png\n")
