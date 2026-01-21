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
Combine <- readRDS("ACM_VCM_SAN_annotated.rds")

levels(Combine) <- c("ACM","Endothelial","Epicardial","Epithelial","FB_1","FB_2",
                     "Neuronal","Neural_Crest","Proliferating","SAN","VCM")

# UMAP demonstrate different cell populations based on RNA expression
# show in the default page
UMAPPlot(Combine,label=TRUE,label.size =5,repel = TRUE) + labs(x = 'UMAP_1', y = 'UMAP_2')

UMAPPlot(Combine,label=TRUE,split.by = "orig.ident",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

# dotplot demonstrate the marker expression of each cell population
# show in the default page
DotPlot(Combine,features = c("NPPA","CDH5","WT1","EPCAM","POSTN","PDZRN4","STMN2","SOX2","TOP2A","SHOX2","MYL2")) + RotatedAxis()

# Violin plot showing the expression of indicated gene
# show after type in the search gene 
VlnPlot(Combine,features = "SHOX2",pt.size = 0) + labs(x = '')

# UMAP to demonstrate the expression of selected gene
# show after type in the search gene 
FeaturePlot(Combine,features = "SHOX2") + labs(x = 'UMAP_1',y = "UMAP_2")

# Dotplot to demonstrate the expression of selected gene
# show after type in the search gene 
DotPlot(Combine,features = "SHOX2") + labs(x = 'UMAP_1', y = 'UMAP_2')

# Plot generator for ACM_VCM_SAN data
acmvcmOmics <- function(genes, png_path) {
  tryCatch({
    p2 <- DotPlot(Combine,features = genes) + RotatedAxis() + labs(x = '',y = '')
    p3 <- FeaturePlot(Combine,features = genes) + labs(x = 'UMAP_1',y = "UMAP_2")
    p4 <- VlnPlot(Combine,features = genes,pt.size = 0) + labs(x = '')

    combined <- p2 + p3 + p4 + plot_layout(ncol = 2)
    ggsave(png_path, plot = combined, width = 15, height = 10, device = "png")
    cat("Saved image to", png_path, "\n")
  }, error = function(e) {
    cat("Error generating ACM_VCM_SAN plot:", e$message, "\n")
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
      gene_list <- rownames(Combine)
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

      if (!(gene_name %in% rownames(Combine))) {
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
      acmvcmOmics(gene_name, png_file)
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

# Start server on port 9027
server <- startServer("0.0.0.0", 9027, app)
cat("ACM_VCM_SAN Omics server started at http://localhost:9027\n")

# Keep it running
while (TRUE) {
  service()
  Sys.sleep(0.001)
}
