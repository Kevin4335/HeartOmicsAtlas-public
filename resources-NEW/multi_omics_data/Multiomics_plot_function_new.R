# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(jsonlite)
library(httpuv)

# load the data
fetal_heart <- readRDS("fetal_heart_multiomics_annotated.rds")

# UMAP demonstrate different cell populations based on RNA expression
DimPlot(fetal_heart, reduction = "umap.rna", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("RNA") + labs(x = 'UMAP_1',y = 'UMAP_2')

# UMAP demonstrate different cell populations based on ATAC data
DimPlot(fetal_heart, reduction = "umap.atac.peaks",label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC") + labs(x = 'UMAP_1',y = 'UMAP_2')

# dotplot demonstrate the marker expression of each cell population
DefaultAssay(fetal_heart) <-'RNA'
levels(fetal_heart) <- c('VCM','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','ACM')
DotPlot(fetal_heart,features = c("MYL2","MYH11","MYBPC1","HCN4","HBG1","STMN2","CD163","CD96","PDGFRA","PAX1","WT1","CDH5","NPPA")) + RotatedAxis() + labs(x = '',y = '')

# Setup fragments for ATAC coverage plots (done once)
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

# Function to generate MultiOmics plots for a gene
multiOmics <- function(gene, png_path) {
  tryCatch({
    # UMAP plot showing the expression of indicated gene
    DefaultAssay(fetal_heart) <-'RNA'
    p3 <- FeaturePlot(fetal_heart,features = gene) + labs(x = 'UMAP_1',y = 'UMAP_2')
    
    # Violin plot showing the expression of indicated gene
    DefaultAssay(fetal_heart) <-'RNA'
    p4 <- VlnPlot(fetal_heart,features = gene,pt.size = 0) + labs(x = '')
    
    # IGV coverage plot showing the peaks at the regulatory region of indicated gene
    DefaultAssay(fetal_heart) <-'peaks'
    p5 <- CoveragePlot(fetal_heart, region = gene,features = gene,annotation = TRUE,peaks = FALSE, pextend.upstream = 2000, extend.downstream = 2000)
    
    combined_plot <- p3 + p4 + p5 + plot_layout(ncol=2)
    ggsave(png_path, plot = combined_plot, width = 15, height = 10, device = "png")
    cat("Saved MultiOmics image to", png_path, "\n")
  }, error = function(e) {
    cat("Error generating MultiOmics plot:", e$message, "\n")
  })
}

# Define HTTP app
app <- list(
  call = function(req) {
    url <- req$PATH_INFO

    if (url == "/genes") {
      gene_list <- rownames(fetal_heart[["RNA"]])
      response_body <- toJSON(gene_list, auto_unbox = TRUE)
      return(list(
        status = 200L,
        headers = list(
          'Access-Control-Allow-Origin' = '*',
          'Content-Type' = 'application/json',
          'Content-Length' = as.character(nchar(response_body))
        ),
        body = response_body
      ))
    }

    if (grepl("^/genes/", url)) {
      gene_name <- sub("^/genes/", "", url)

      if (!(gene_name %in% rownames(fetal_heart[["RNA"]]))) {
        return(list(
          status = 404L,
          headers = list(
            'Access-Control-Allow-Origin' = '*',
            'Content-Type' = 'text/plain'
          ),
          body = paste("Gene not found:", gene_name)
        ))
      }

      png_file <- tempfile(fileext = ".png")
      multiOmics(gene_name, png_file)
      png_data <- readBin(png_file, what = "raw", n = file.info(png_file)$size)

      return(list(
        status = 200L,
        headers = list(
          'Access-Control-Allow-Origin' = '*',
          'Content-Type' = 'image/png',
          'Content-Length' = as.character(length(png_data))
        ),
        body = png_data
      ))
    }

    return(list(
      status = 400L,
      headers = list(
        'Access-Control-Allow-Origin' = '*',
        'Content-Type' = 'text/plain'
      ),
      body = "Invalid request"
    ))
  }
)

# Start server on port 9026
server <- startServer("0.0.0.0", 9026, app)
cat("MultiOmics Server started on http://localhost:9026\n")

while (TRUE) {
  service()
  Sys.sleep(0.001)
}