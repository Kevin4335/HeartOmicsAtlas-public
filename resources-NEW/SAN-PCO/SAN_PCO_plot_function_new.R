# load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(httpuv)
library(jsonlite)

# load the data
SAN_PCO <- readRDS("SAN-PCO_annotated.rds")

levels(SAN_PCO) <- c("ACM","Endothelial","Epicardial","Epithelial","FB",
                     "Neuronal","Neural_Crest","Proliferating","SAN")

# UMAP demonstrate different cell populations based on RNA expression
# show in the default page
UMAPPlot(SAN_PCO,label=TRUE,label.size = 4, repel = TRUE) + labs(x = "UMAP_1",y = "UMAP_2")

# dotplot demonstrate the marker expression of each cell population
# show in the default page
DotPlot(SAN_PCO,features = c("NPPA","CDH5","WT1","EPCAM","COL1A1","STMN2","SOX2","TOP2A","SHOX2")) + RotatedAxis() + labs(x = '',y = '')

# Violin plot showing the expression of indicated gene
# show after type in the search gene 
VlnPlot(SAN_PCO,features = "STMN2",pt.size = 0) + labs(x = "")

# UMAP to demonstrate the expression of selected gene
# show after type in the search gene 
FeaturePlot(SAN_PCO,features = "STMN2") + labs(x = "UMAP_1",y = "UMAP_2")

# Dotplot to demonstrate the expression of selected gene
# show after type in the search gene 
DotPlot(SAN_PCO,features = "STMN2") + labs(x = "",y = "")

# Plot generator for SAN_PCO data
sanpcoOmics <- function(genes, png_path) {
  tryCatch({
    p2 <- DotPlot(SAN_PCO,features = genes) + labs(x = "",y = "")
    p3 <- FeaturePlot(SAN_PCO,features = genes) + labs(x = "UMAP_1",y = "UMAP_2")
    p4 <- VlnPlot(SAN_PCO,features = genes,pt.size = 0) + labs(x = "")

    combined <- p2 + p3 + p4 + plot_layout(ncol = 2)
    ggsave(png_path, plot = combined, width = 15, height = 10, device = "png")
    cat("Saved image to", png_path, "\n")
  }, error = function(e) {
    cat("Error generating SAN_PCO plot:", e$message, "\n")
  })
}

# Decode hex to string (optional, not used unless extended)
hex_to_string <- function(hex_str) {
  hex_split <- strsplit(hex_str, "(?<=..)", perl = TRUE)[[1]]
  raw_vec <- as.raw(as.hexmode(hex_split))
  result_str <- rawToChar(raw_vec)
  return(result_str)
}

# HTTP server definition
app <- list(
  call = function(req) {
    url <- req$PATH_INFO

    if (url == "/genes") {
      gene_list <- rownames(SAN_PCO)
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

      if (!(gene_name %in% rownames(SAN_PCO))) {
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
      sanpcoOmics(gene_name, png_file)
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

    response_body <- "Invalid endpoint"
    return(list(
      status = 404L,
      headers = list(
        'Content-Type' = 'text/plain'
      ),
      body = response_body
    ))
  }
)

# Start server on port 9028 (adjust if needed)
server <- startServer("0.0.0.0", 9028, app)
cat("SAN_PCO Omics server started at http://localhost:9028\n")

# Keep it running
while (TRUE) {
  service()
  Sys.sleep(0.001)
}
