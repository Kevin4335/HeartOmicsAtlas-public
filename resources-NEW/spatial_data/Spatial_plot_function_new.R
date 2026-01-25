# load the required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(httpuv)
library(jsonlite)

#### example genes to list in the search panel
#### SHOX2, MYH11, MYL7, PLP1

# load the data
fetal_heart <- readRDS("fetal_heart_0103_16um_annotated.rds")
fetal_heart <- RenameIdents(fetal_heart, 'EC' = 'Endothelial','FB' = 'Fibroblast')
fetal_heart$celltype <- Idents(fetal_heart)
levels(fetal_heart) <- c("ACM","Endothelial","Fibroblast","Macrophages","Neuronal",'RBC',"SAN","SMC")

# Update Seurat object - required for SpatialFeaturePlot to align with tissue_hires_image.png
tryCatch({
  fetal_heart <- UpdateSeuratObject(fetal_heart)
  cat("Seurat object updated successfully\n")
}, error = function(e) {
  cat("Warning: UpdateSeuratObject failed:", conditionMessage(e), "\n")
})

# Function to generate plots for a given gene and save to PNG
# Spatial plot showing the expression of indicated gene - first plot when search for specific gene
# UMAP plot showing the expression of indicated gene - second plot when search for specific gene
# Violin plot showing the expression of indicated gene - third plot when search for specific gene
# Dot plot showing the expression of indicated gene - fourth plot when search for specific gene
spatialOmics <- function(genes, png_path) {
  tryCatch({
    # Spatial plot showing the expression of indicated gene
    p3 <- SpatialFeaturePlot(fetal_heart,features = c(genes),alpha = c(0.1,1),pt.size.factor = 5) +
      theme(legend.text = element_text(size = 15),
            legend.title = element_text(size = 18,vjust = 0.8))

    # UMAP plot showing the expression of indicated gene
    p4 <- FeaturePlot(fetal_heart,features = c(genes)) + labs(x='UMAP1',y='UMAP2')+
      theme(legend.text = element_text(size = 16),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))

    # Violin plot showing the expression of indicated gene
    p5 <- VlnPlot(fetal_heart,features = c(genes),pt.size = 0) + labs(x='') +
      theme(plot.title = element_text(size = 20,hjust = 0.5),
            legend.text = element_text(size = 16),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))

    # Dot plot showing the expression of indicated gene
    p6 <- DotPlot(fetal_heart,features = c(genes)) + labs(x = '',y = '')+
      theme(axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16)) +
      RotatedAxis()

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
      gene_list <- rownames(fetal_heart)
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

      png_file <- tempfile(fileext = ".png")
      spatialOmics(gene_name, png_file)
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

while (TRUE) {
  service()
  Sys.sleep(0.001)
}
