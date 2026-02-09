# load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)

#### example genes to list in the search panel
#### SHOX2, NPPA, MYL7, CDH5

# load the data
heart <- readRDS("mini_heart_annotated.rds")
heart <- RenameIdents(heart, 'FB' = 'Fibroblast')
heart$celltype <- Idents(heart)
levels(heart) <- c("ACM","Endothelial","Epicardial","Epithelial","Fibroblast",
                   "Neuronal","Proliferating","SAN","VCM")

# UMAP demonstrate different cell populations based on RNA expression
# first plot to demonstrate in default page
UMAPPlot(heart,label=TRUE,label.size =6,repel = TRUE) + labs(x = 'UMAP_1',y = 'UMAP_2') + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# dotplot demonstrate the marker expression of each cell population
# second plot to demonstrate in default page
DotPlot(heart,features = c("NPPA","CDH5","WT1","EPCAM","DCN","STMN2","TOP2A","SHOX2","HEY2")) + 
  labs(x='',y='') + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
  RotatedAxis()

# Violin plot showing the expression of indicated gene
# first plot to show when search for specific gene
VlnPlot(heart,features = "SHOX2",pt.size = 0) + labs(x = '') +
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# UMAP to demonstrate the expression of selected gene
# second plot to show when search for specific gene
FeaturePlot(heart,features = "SHOX2") + 
  labs(x = 'UMAP_1',y = "UMAP_2") +
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# Dotplot to demonstrate the expression of selected gene
# third plot to show when search for specific gene
DotPlot(heart,features = "SHOX2")  + labs(x = '',y = '')+ 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
  RotatedAxis()

