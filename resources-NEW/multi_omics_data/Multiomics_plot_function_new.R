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

# Function to generate MultiOmics plots for a gene
multiOmics <- function(gene, png_path) {
  tryCatch({
    # UMAP plot showing the expression of indicated gene
    DefaultAssay(fetal_heart) <-'RNA'
    p3 <- FeaturePlot(fetal_heart,features = gene) + labs(x = 'UMAP_1',y = 'UMAP_2')
    
    # Violin plot showing the expression of indicated gene
    DefaultAssay(fetal_heart) <-'RNA'
    p4 <- VlnPlot(fetal_heart,features = gene,pt.size = 0) + labs(x = '')
    
    # Setup fragments for ATAC coverage plots
    # Use direct path assignment (like old version) to avoid "Cells already present" error
    DefaultAssay(fetal_heart) <-'peaks'
    newpath1 <- "3655/atac_fragments.tsv.gz"
    newpath2 <- "3675/atac_fragments.tsv.gz"
    newpath3 <- "7668/atac_fragments.tsv.gz"
    newpath4 <- "3688/atac_fragments.tsv.gz"
    
    # Directly modify fragment paths in the object (like old version)
    # This avoids the "Cells already present" error by not reassigning fragments
    if (length(fetal_heart@assays[["peaks"]]@fragments) >= 4) {
      fetal_heart@assays[["peaks"]]@fragments[[1]]@path <- newpath1
      fetal_heart@assays[["peaks"]]@fragments[[2]]@path <- newpath2
      fetal_heart@assays[["peaks"]]@fragments[[3]]@path <- newpath3
      fetal_heart@assays[["peaks"]]@fragments[[4]]@path <- newpath4
    }
    
    # IGV coverage plot showing the peaks at the regulatory region of indicated gene
    p5 <- CoveragePlot(fetal_heart, region = gene,features = gene,annotation = TRUE,peaks = FALSE)
    
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
