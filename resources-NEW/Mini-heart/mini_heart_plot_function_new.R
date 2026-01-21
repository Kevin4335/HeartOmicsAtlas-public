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
heart <- readRDS("mini_heart_annotated.rds")

levels(heart) <- c("ACM","Endothelial","Epicardial","Epithelial","FB",
                   "Neuronal","Proliferating","SAN","VCM")

# Plot generator for Mini-heart data
miniHeartOmics <- function(genes, png_path) {
  tryCatch({
    # Violin plot showing the expression of indicated gene
    # show after type in the search gene 
    p4 <- VlnPlot(heart,features = genes,pt.size = 0) + labs(x = '')
    
    # UMAP to demonstrate the expression of selected gene
    # show after type in the search gene 
    p3 <- FeaturePlot(heart,features = genes) + labs(x = 'UMAP_1',y = "UMAP_2")
    
    # Dotplot to demonstrate the expression of selected gene
    # show after type in the search gene 
    p2 <- DotPlot(heart,features = genes) + labs(x = 'UMAP_1', y = 'UMAP_2')

    combined <- p2 + p3 + p4 + plot_layout(ncol = 2)
    ggsave(png_path, plot = combined, width = 15, height = 10, device = "png")
    cat("Saved image to", png_path, "\n")
  }, error = function(e) {
    cat("Error generating Mini-heart plot:", e$message, "\n")
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
      gene_list <- rownames(heart)
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

      if (!(gene_name %in% rownames(heart))) {
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
      miniHeartOmics(gene_name, png_file)
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

# Start server on port 9029 (new port)
server <- startServer("0.0.0.0", 9029, app)
cat("Mini-heart Omics server started at http://localhost:9029\n")

# Keep it running
while (TRUE) {
  service()
  Sys.sleep(0.001)
}
