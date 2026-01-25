# Default plots for Multiomics (generate manually and save to PNG)
# These plots are shown on the default page

# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

# load the data
fetal_heart <- readRDS("fetal_heart_multiomics_annotated.rds")
levels(fetal_heart) <- c('ACM','Endothelial','Epicardial','Epithelial','Fibroblast','Lymphoid','Myeloid','Neuronal_1','Neuronal_2','RBC','SAN','SMC','VCM')

# UMAP demonstrate different cell populations based on RNA expression
# first plot to demonstrate in default page
p1 <- DimPlot(fetal_heart, reduction = "umap.rna", label = TRUE, label.size = 6, repel = TRUE) +
  ggtitle("RNA") + labs(x = 'UMAP_1',y = 'UMAP_2') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_fixed()

# UMAP demonstrate different cell populations based on ATAC data
# second plot to demonstrate in default page
p2 <- DimPlot(fetal_heart, reduction = "umap.atac.peaks",label = TRUE, label.size = 6, repel = TRUE) +
  ggtitle("ATAC") + labs(x = 'UMAP_1',y = 'UMAP_2') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_fixed()

# UMAP demonstrate different cell populations based on WNN
# third plot to demonstrate in default page
p3 <- DimPlot(fetal_heart, reduction = "wnn.umap",label = TRUE, label.size = 6, repel = TRUE) +
  ggtitle("WNN") + labs(x = 'UMAP_1',y = 'UMAP_2') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  coord_fixed()

# dotplot demonstrate the marker expression of each cell population
# fourth plot to demonstrate in default page
DefaultAssay(fetal_heart) <-'RNA'
p4 <- DotPlot(fetal_heart,features = c("NPPA","CDH5","WT1","PAX1","PDGFRA","CD96","CD163","SYN3","NRXN1","HBG1","HCN1","MYH11","MYL2"),group.by = "celltype") +
  labs(x='',y='') +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  RotatedAxis()

# Combine plots: 2x2 grid
combined <- (p1 + p2) / (p3 + p4)

# Save to PNG (20x20 so each panel is 10x10, comparable to Spatial/ACM_VCM_SAN which use width 20)
ggsave("Multiomics_default_plots.png", plot = combined, width = 20, height = 20, device = "png")
cat("Saved default plots to Multiomics_default_plots.png\n")
