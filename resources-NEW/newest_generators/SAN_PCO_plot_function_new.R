# load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)

#### example genes to list in the search panel
#### SHOX2, NPPA, CDH5, STMN2

# load the data
SAN_PCO <- readRDS("SAN-PCO_annotated.rds")
SAN_PCO <- RenameIdents(SAN_PCO, 'FB' = 'Fibroblast')
SAN_PCO$celltype <- Idents(SAN_PCO)
levels(SAN_PCO) <- c("ACM","Endothelial","Epicardial","Epithelial","Fibroblast",
                     "Neuronal","Neural_Crest","Proliferating","SAN")

# UMAP demonstrate different cell populations based on RNA expression
# first plot to demonstrate in default page
UMAPPlot(SAN_PCO,label=TRUE,label.size = 6, repel = TRUE) + labs(x = 'UMAP_1',y = 'UMAP_2') + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# dotplot demonstrate the marker expression of each cell population
# second plot to demonstrate in default page
DotPlot(SAN_PCO,features = c("NPPA","CDH5","WT1","EPCAM","COL1A1","STMN2","SOX2","TOP2A","SHOX2")) +
  labs(x='',y='') + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
  RotatedAxis()

# Violin plot showing the expression of indicated gene
# first plot to show when search for specific gene
VlnPlot(SAN_PCO,features = "STMN2",pt.size = 0) + labs(x = '') +
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# UMAP to demonstrate the expression of selected gene
# second plot to show when search for specific gene
FeaturePlot(SAN_PCO,features = "STMN2") + 
  labs(x = 'UMAP_1',y = "UMAP_2") +
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# Dotplot to demonstrate the expression of selected gene
# third plot to show when search for specific gene
DotPlot(SAN_PCO,features = "STMN2") + labs(x = '',y = '')+ 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
  RotatedAxis()

