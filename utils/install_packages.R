# install_packages.R

cran_packages <- c(
  "Seurat",
  "patchwork",
  "ggplot2",
  "httpuv",
  "jsonlite",
  "dplyr",
  "EnhancedVolcano",
  "RColorBrewer",
  "Signac"
)

# Install CRAN packages
install.packages("BiocManager", repos = "http://cran.rstudio.com/")
install.packages(cran_packages, repos = "http://cran.rstudio.com/")

# Install Bioconductor packages
BiocManager::install(c("EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38"), ask = FALSE)

# Verify all packages installed
required <- c(cran_packages, "BiocManager", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38")
missing <- setdiff(required, rownames(installed.packages()))

if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
} else {
  cat("✅ All packages installed successfully\n")
}
