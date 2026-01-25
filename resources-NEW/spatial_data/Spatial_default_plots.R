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
fetal_heart <- RenameIdents(fetal_heart, 'EC' = 'Endothelial','FB' = 'Fibroblast')
fetal_heart$celltype <- Idents(fetal_heart)
levels(fetal_heart) <- c("ACM","Endothelial","Fibroblast","Macrophages","Neuronal",'RBC',"SAN","SMC")

# Update Seurat object - required for SpatialDimPlot to align with tissue_hires_image.png
tryCatch({
  fetal_heart <- UpdateSeuratObject(fetal_heart)
  cat("Seurat object updated successfully\n")
}, error = function(e) {
  cat("Warning: UpdateSeuratObject failed:", conditionMessage(e), "\n")
  cat("SpatialDimPlot may fail; if so, ensure ggplot2 is 3.5.x (not 4.x)\n")
})

# UMAP demonstrate different cell populations
# first plot to demonstrate in default page
p1 <- DimPlot(fetal_heart,reduction = "umap.016um",label = TRUE,label.size = 6, repel = TRUE) +
  labs(x = "UMAP1",y = "UMAP2") +
  theme(legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position = "none")

# Spatial plot demonstrate different cell populations
# second plot to demonstrate in default page
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
p2 <- SpatialDimPlot(fetal_heart,label = TRUE,label.size = 5,repel = TRUE) + scale_fill_manual(values = mycols) +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.position = "none")

# Spatial plot to show the marker expression of each cluster using Dotplot
# third plot to demonstrate in default page
p3 <- DotPlot(fetal_heart,features = c("MYL7","VWF","DCN","SPP1","PLP1","HBG1","SHOX2","MYH11")) +
  labs(x='',y='') +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  RotatedAxis()

# Combine plots: 2 on top, 1 on bottom (bottom plot same width as one top plot)
combined <- (p1 + p2) / (p3 + plot_spacer())

# Save to PNG in the same directory
ggsave("Spatial_default_plots.png", plot = combined, width = 20, height = 15, device = "png")
cat("Saved default plots to Spatial_default_plots.png\n")
