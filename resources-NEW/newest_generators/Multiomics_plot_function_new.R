# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

#### example genes to list in the search panel
#### TBX5, NPPA, MYL2, SYN3
# load the data
fetal_heart <- readRDS("fetal_heart_multiomics_annotated.rds")
levels(fetal_heart) <- c('ACM','Endothelial','Epicardial','Epithelial','Fibroblast','Lymphoid','Myeloid','Neuronal_1','Neuronal_2','RBC','SAN','SMC','VCM')

# UMAP demonstrate different cell populations based on RNA expression
# first plot to demonstrate in default page
DimPlot(fetal_heart, reduction = "umap.rna", label = TRUE, label.size = 6, repel = TRUE) + 
  ggtitle("RNA") + labs(x = 'UMAP_1',y = 'UMAP_2') + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(legend.text = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# UMAP demonstrate different cell populations based on ATAC data
# second plot to demonstrate in default page
DimPlot(fetal_heart, reduction = "umap.atac.peaks",label = TRUE, label.size = 6, repel = TRUE) +
  ggtitle("ATAC") + labs(x = 'UMAP_1',y = 'UMAP_2')  + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(legend.text = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# UMAP demonstrate different cell populations based on WNN 
# third plot to demonstrate in default page
DimPlot(fetal_heart, reduction = "wnn.umap",label = TRUE, label.size = 6, repel = TRUE) +
  ggtitle("WNN") + labs(x = 'UMAP_1',y = 'UMAP_2')  + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(legend.text = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# dotplot demonstrate the marker expression of each cell population
# fourth plot to demonstrate in default page
DefaultAssay(fetal_heart) <-'RNA'

DotPlot(fetal_heart,features = c("NPPA","CDH5","WT1","PAX1","PDGFRA","CD96","CD163","SYN3","NRXN1","HBG1","HCN1","MYH11","MYL2"),group.by = "celltype") + 
  labs(x='',y='') + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) + 
  RotatedAxis()

# UMAP plot showing the expression of indicated gene
# first plot to show when search for specific gene
DefaultAssay(fetal_heart) <-'RNA'
FeaturePlot(fetal_heart,features = c("TBX5")) + labs(x = 'UMAP_1',y = 'UMAP_2') + 
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# Violin plot showing the expression of indicated gene
# second plot to show when search for specific gene
DefaultAssay(fetal_heart) <-'RNA'
VlnPlot(fetal_heart,features = c("TBX5"),pt.size = 0) + labs(x = '') + 
  theme(plot.title = element_text(size = 20,hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# IGV coverage plot showing the peaks at the regulatory region of  indicated gene
# third plot to show when search for specific gene
DefaultAssay(fetal_heart) <-'ATAC'
frags <- Fragments(fetal_heart[['ATAC']])
Fragments(fetal_heart[['ATAC']]) <- NULL
newpath1 <- "3655/atac_fragments.tsv.gz"
newpath2 <- "3675/atac_fragments.tsv.gz"
newpath3 <- "7668/atac_fragments.tsv.gz"
newpath4 <- "3688/atac_fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]],new.path = newpath1)
frags[[2]] <- UpdatePath(frags[[2]],new.path = newpath2)
frags[[3]] <- UpdatePath(frags[[3]],new.path = newpath3)
frags[[4]] <- UpdatePath(frags[[4]],new.path = newpath4)
Fragments(fetal_heart[['ATAC']]) <- frags
DefaultAssay(fetal_heart) <-'peaks'
Fragments(fetal_heart[['peaks']]) <- Fragments(fetal_heart[['ATAC']])
CoveragePlot(fetal_heart, region = "TBX5",features = "TBX5",annotation = TRUE,peaks = FALSE)


