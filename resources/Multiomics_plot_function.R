# Load necessary packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(jsonlite)
library(httpuv)

# Load MultiOmics data
fetal_heart_multi <- readRDS("fetal_heart_multiomics_annotated.rds")
levels(fetal_heart_multi) <- c('VCM','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','ACM')
fragment_dirs <- c("3655", "3675", "7668", "3688")
fragment_paths <- paste0("../fragments/", fragment_dirs, "/atac_fragments.tsv.gz")

for (i in seq_along(fragment_paths)) {
  fetal_heart_multi@assays[["peaks"]]@fragments[[i]]@path <- fragment_paths[[i]]
}
# Function to generate MultiOmics plots for a gene
multiOmics <- function(gene, png_path) {
  tryCatch({
    DefaultAssay(fetal_heart_multi) <- 'RNA'
    p3 <- FeaturePlot(fetal_heart_multi,features = gene) + labs(x = 'UMAP_1',y = 'UMAP_2')
    p4 <- VlnPlot(fetal_heart_multi,features = gene,pt.size = 0) + labs(x = '')

    DefaultAssay(fetal_heart_multi) <- 'peaks'
    levels(fetal_heart_multi) <- c('VCM','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','ACM')

    p5 <- CoveragePlot(fetal_heart_multi, region = gene,features = gene,annotation = TRUE,peaks = FALSE, pextend.upstream = 2000, extend.downstream = 2000)


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
      gene_list <- rownames(fetal_heart_multi[["RNA"]])
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

      if (!(gene_name %in% rownames(fetal_heart_multi[["RNA"]]))) {
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
