library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(httpuv)
library(jsonlite)

# Load Sinoid dataset
Sinoid <- readRDS("Sinoid_annotated.rds")
levels(Sinoid) <- c("SAN","Proliferating","Neuronal","Epithelial","Epicardial","Endothelial")

# Function to generate plots and save to PNG
sinoidOmics <- function(genes, png_path) {
  tryCatch({
    p2 <- DotPlot(Sinoid,features = genes) + labs(x = 'UMAP_1', y = 'UMAP_2')
    p3 <- FeaturePlot(Sinoid,features = genes) + labs(x = 'UMAP_1',y = "UMAP_2")
    p4 <- VlnPlot(Sinoid,features = genes,pt.size = 0.5) + labs(x = '')

    combined <- p2 + p3 + p4 + plot_layout(ncol = 2)
    ggsave(png_path, plot = combined, width = 15, height = 10, device = "png")
    cat("Saved image to", png_path, "\n")
  }, error = function(e) {
    cat("Error generating Sinoid plot:", e$message, "\n")
  })
}

# Decode hex to string (if using hex-based JSON API later)
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

    # Endpoint for full gene list
    if (url == "/genes") {
      gene_list <- rownames(Sinoid)
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

    # Endpoint for single gene PNG
    if (grepl("^/genes/", url)) {
      gene_name <- sub("^/genes/", "", url)

      if (!(gene_name %in% rownames(Sinoid))) {
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
      sinoidOmics(gene_name, png_file)

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

    # Optional: fallback handler
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

# Start server on port 9026 (or whatever’s free)
server <- startServer("0.0.0.0", 9027, app)
cat("SinoidOmics server started at http://localhost:9027\n")

# Keep server running
while (TRUE) {
  service()
  Sys.sleep(0.001)
}
