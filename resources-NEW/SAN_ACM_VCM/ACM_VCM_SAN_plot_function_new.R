# load the required packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(httpuv)
library(jsonlite)

#### example genes to list in the search panel
#### SHOX2, NPPA, MYL7, STMN2

# load the data
Combine <- readRDS("ACM_VCM_SAN_annotated.rds")
Combine <- RenameIdents(Combine, 'FB_1' = 'Fibroblast_1','FB_2' = 'Fibroblast_2')
Combine$celltype <- Idents(Combine)
levels(Combine) <- c("ACM","Endothelial","Epicardial","Epithelial","Fibroblast_1","Fibroblast_2",
                     "Neuronal","Neural_Crest","Proliferating","SAN","VCM")
# Origin for Sinoid / ACO / VCO (from revised script)
Combine$origin <- ifelse(Combine$orig.ident == 'ACM', 'ACO',
                         ifelse(Combine$orig.ident == 'VCM', 'VCO', 'Sinoid'))
Combine$origin <- factor(Combine$origin, levels = c('ACO', 'VCO', 'Sinoid'))

# Plot generator: subtype "sinoid", "aco", or "vco" filters to that origin; 4-panel layout from revised script
acmvcmOmics <- function(genes, png_path, subtype = "sinoid") {
  tryCatch({
    # Map frontend subtype to origin level
    origin_level <- switch(subtype, "aco" = "ACO", "vco" = "VCO", "Sinoid")
    if (is.null(origin_level)) origin_level <- "Sinoid"
    sub <- subset(Combine, subset = origin == origin_level)

    # 1) Violin (subset)
    p1 <- VlnPlot(sub, features = genes, pt.size = 0) + labs(x = '') +
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            legend.text = element_text(size = 16),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))

    # 2) Feature plot (subset)
    p2 <- FeaturePlot(sub, features = genes) +
      labs(x = 'UMAP_1', y = "UMAP_2") +
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            legend.text = element_text(size = 16),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))

    # 3) Cell type UMAP for selected type only (subset) - no longer split by all three origins
    p3 <- DimPlot(sub, label = TRUE, label.size = 5, repel = TRUE) +
      labs(x = 'UMAP_1', y = "UMAP_2", title = paste0("Cell types (", origin_level, ")")) +
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            legend.text = element_text(size = 16),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))

    # 4) Dotplot (subset)
    p4 <- DotPlot(sub, features = genes) + labs(x = '', y = '') +
      theme(axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16)) +
      RotatedAxis()

    combined <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
    ggsave(png_path, plot = combined, width = 20, height = 14, device = "png")
    cat("Saved image to", png_path, "\n")
  }, error = function(e) {
    cat("Error generating plot:", conditionMessage(e), "\n")
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
      # Remove any query string from gene_name (e.g. "SHOX2" from "SHOX2?subtype=aco")
      gene_name <- strsplit(gene_name, "?", fixed = TRUE)[[1]][1]

      # Parse subtype from query string (sinoid, aco, vco)
      subtype <- "sinoid"
      qs <- req$QUERY_STRING
      if (!is.null(qs) && nchar(qs) > 0) {
        qs <- sub("^[?]", "", qs)  # strip leading ? if present
        pairs <- strsplit(qs, "&", fixed = TRUE)[[1]]
        for (p in pairs) {
          kv <- strsplit(p, "=", fixed = TRUE)[[1]]
          if (length(kv) >= 2 && trimws(kv[1]) == "subtype") {
            subtype <- trimws(utils::URLdecode(kv[2]))
            if (!subtype %in% c("sinoid", "aco", "vco")) subtype <- "sinoid"
            break
          }
        }
      }

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
      acmvcmOmics(gene_name, png_file, subtype)
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

while (TRUE) {
  service()
  Sys.sleep(0.001)
}
