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
Combine <- RenameIdents(Combine, 'FB_1' = 'Fibroblast_1','FB_2' = 'Fibroblast_2')
Combine$celltype <- Idents(Combine)
levels(Combine) <- c("ACM","Endothelial","Epicardial","Epithelial","Fibroblast_1","Fibroblast_2",
                     "Neuronal","Neural_Crest","Proliferating","SAN","VCM")

# UMAP demonstrate different cell populations based on RNA expression
# first plot to demonstrate in default page
p1 <- UMAPPlot(Combine,label=TRUE,label.size =6,repel = TRUE) + labs(x = 'UMAP_1',y = 'UMAP_2') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# UMAP demonstrate different cell populations based on RNA expression with splitted demonstration
# second plot to demonstrate in default page
p2 <- UMAPPlot(Combine,label=TRUE,split.by = "orig.ident",label.size =5,repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# dotplot demonstrate the marker expression of each cell population
# third plot to demonstrate in default page
p3 <- DotPlot(Combine,features = c("NPPA","CDH5","WT1","EPCAM","POSTN","PDZRN4","STMN2","SOX2","TOP2A","SHOX2","MYL2")) +
  labs(x='',y='') +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  RotatedAxis()

# Combine plots: split UMAP (3 panels) full first row; single UMAP + DotPlot on bottom
combined <- p2 / (p1 + p3)

# Save to PNG in the same directory
ggsave("ACM_VCM_SAN_default_plots.png", plot = combined, width = 20, height = 15, device = "png")
cat("Saved default plots to ACM_VCM_SAN_default_plots.png\n")
