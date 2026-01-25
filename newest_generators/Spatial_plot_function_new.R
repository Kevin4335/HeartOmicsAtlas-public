# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

#### example genes to list in the search panel
#### SHOX2, MYH11, MYL7, PLP1

# load the data
fetal_heart <- readRDS("fetal_heart_0103_16um_annotated.rds")
fetal_heart <- RenameIdents(fetal_heart, 'EC' = 'Endothelial','FB' = 'Fibroblast')
fetal_heart$celltype <- Idents(fetal_heart)
levels(fetal_heart) <- c("ACM","Endothelial","Fibroblast","Macrophages","Neuronal",'RBC',"SAN","SMC")

# UMAP demonstrate different cell populations
# first plot to demonstrate in default page
DimPlot(fetal_heart,reduction = "umap.016um",label = TRUE,label.size = 6, repel = TRUE) + 
  labs(x = "UMAP1",y = "UMAP2") + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
  

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
SpatialDimPlot(fetal_heart,label = TRUE,label.size = 5,repel = TRUE) + scale_fill_manual(values = mycols) + 
  guides(fill = guide_legend(override.aes = list(size = 10))) + 
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 20))

# Spatial plot to show the marker expression of each cluster using Dotplot
# third plot to demonstrate in default page
DotPlot(fetal_heart,features = c("MYL7","VWF","DCN","SPP1","PLP1","HBG1","SHOX2","MYH11")) + 
  labs(x='',y='') + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
  RotatedAxis()

# Spatial plot showing the expression of indicated gene
# first plot to show when search for specific gene
SpatialFeaturePlot(fetal_heart,features = c("SHOX2"),alpha = c(0.1,1),pt.size.factor = 5) + 
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 18,vjust = 0.8))

# UMAP plot showing the expression of indicated gene
# second plot to show when search for specific gene
FeaturePlot(fetal_heart,features = c("SHOX2")) + labs(x='UMAP1',y='UMAP2')+
  theme(legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# Violin plot showing the expression of indicated gene
# third plot to show when search for specific gene
VlnPlot(fetal_heart,features = c("SHOX2"),pt.size = 0) + labs(x='') + 
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# Dot plot showing the expression of indicated gene
# fourth plot to show when search for specific gene
DotPlot(fetal_heart,features = c("SHOX2")) + labs(x = '',y = '')+ 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
  RotatedAxis()
