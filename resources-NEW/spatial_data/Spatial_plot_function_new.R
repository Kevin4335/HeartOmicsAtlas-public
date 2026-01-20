# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(CellChat)
library(httpuv)
library(jsonlite)

# load the data
fetal_heart <- readRDS("fetal_heart_0103_16um_annotated.rds")

# UMAP demonstrate different cell populations
DimPlot(fetal_heart,reduction = "umap.016um",label = TRUE, repel = TRUE) + labs(x = "UMAP1",y = "UMAP2")

# Spatial plot demonstrate different cell populations
mycols <- c(
  "orange",
  "darkgreen",
  "#F781BF", 
  "#FF7F00", 
  "#984EA3",
  "#AEC6CF",
  "#E41A1C",
  "#A65628" 
)
SpatialDimPlot(fetal_heart,label = TRUE,label.size = 3,repel = TRUE) + scale_fill_manual(values = mycols) + theme_void()

# Spatial plot to show the marker expression of each cluster using Dotplot
DotPlot(fetal_heart,features = c("MYL7","VWF","DCN","SPP1","PLP1","HBG1","SHOX2","MYH11")) + RotatedAxis() + labs(x = '',y='')

#  Spatial plot showing the expression of indicated gene
SpatialFeaturePlot(fetal_heart,features = c("SHOX2"),alpha = c(0.1,1),pt.size.factor = 5)

# UMAP plot showing the expression of indicated gene
FeaturePlot(fetal_heart,features = c("SHOX2")) + labs(x='UMAP1',y='UMAP2')

# Violin plot showing the expression of indicated gene
VlnPlot(fetal_heart,features = c("SHOX2"),pt.size = 0) + labs(x='')

# Dot plot showing the expression of indicated gene
DotPlot(fetal_heart,features = c("SHOX2")) + RotatedAxis() + labs(x = '',y = '')

# Function to generate plots for a given gene and save to PNG
spatialOmics <- function(genes, png_path) {
  tryCatch({
    # Spatial plot showing the expression of indicated gene
    p3 <- SpatialFeaturePlot(fetal_heart,features = genes,alpha = c(0.1,1),pt.size.factor = 5)
    
    # UMAP plot showing the expression of indicated gene
    p4 <- FeaturePlot(fetal_heart,features = genes) + labs(x='UMAP1',y='UMAP2')
    
    # Violin plot showing the expression of indicated gene
    p5 <- VlnPlot(fetal_heart,features = genes,pt.size = 0) + labs(x='')
    
    # Dot plot showing the expression of indicated gene
    p6 <- DotPlot(fetal_heart,features = genes) + RotatedAxis() + labs(x = '',y = '')
    
    combined_plot <- p3 + p4 + p5 + p6 + plot_layout(ncol = 2)
    ggsave(png_path, plot = combined_plot, width = 15, height = 10, device = "png")
    cat("Saved image to", png_path, "\n")
  }, error = function(e) {
    cat("Error generating plot:", e$message, "\n")
  })
}

# Utility: Decode URL-encoded hex to JSON string
hex_to_string <- function(hex_str) {
  hex_split <- strsplit(hex_str, "(?<=..)", perl = TRUE)[[1]]
  raw_vec <- as.raw(as.hexmode(hex_split))
  result_str <- rawToChar(raw_vec)
  return(result_str)
}

# Define HTTP app
app <- list(
  call = function(req) {
    url <- req$PATH_INFO

    # Serve full gene list at /genes
    if (url == "/genes") {
      gene_list <- rownames(fetal_heart)  # ALL genes
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

    # Serve gene plot at /genes/GENE_NAME
    if (grepl("^/genes/", url)) {
      gene_name <- sub("^/genes/", "", url)

      # Validate gene exists
      if (!(gene_name %in% rownames(fetal_heart))) {
        return(list(
          status = 404L,
          headers = list(
            'Access-Control-Allow-Origin' = '*',
            'Content-Type' = 'text/plain'
          ),
          body = paste("Gene not found:", gene_name)
        ))
      }

      # Generate temp file path for PNG
      png_file <- tempfile(fileext = ".png")
      spatialOmics(gene_name, png_file)

      # Read PNG file as raw binary
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

    # Otherwise, treat as hex JSON fallback
    json_data <- hex_to_string(substr(url, 2, nchar(url)))
    data <- fromJSON(json_data)

    if (data$f == 25) {
      cat("spatialOmics\n")
      spatialOmics(data$p1, data$p2)
    }

    response_body <- "finished"
    return(list(
      status = 200L,
      headers = list(
        'Access-Control-Allow-Origin' = '*',
        'Content-Type' = 'text/plain',
        'Content-Length' = as.character(nchar(response_body))
      ),
      body = response_body
    ))
  }
)

# Start server
server <- startServer("0.0.0.0", 9025, app)
cat("Server started on http://localhost:9025\n")

# Keep it running
while (TRUE) {
  service()
  Sys.sleep(0.001)
}
