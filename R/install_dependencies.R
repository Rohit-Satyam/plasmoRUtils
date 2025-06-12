#' Install package dependencies
#'
#' Installs missing dependencies from both CRAN and Bioconductor
#' @export
install_dependencies <- function() {
  # Check if BiocManager is installed, install if not
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager",dependencies = TRUE)
  }

  # Get package dependencies from DESCRIPTION
  deps <- unique(
    c(renv::dependencies(quiet = TRUE)$Package,
      unlist(tools::package_dependencies("plasmoRUtils", recursive = TRUE))
    ))

    needed <- deps[!deps %in% installed.packages()[,"Package"]]

    if (length(needed) == 0) {
      message("All dependencies are already installed")
      return(invisible())
    }

    # Bioconductor packages (hardcoded but is there an automatic way to do this?)
    bioc_pkgs <- c("rmarkdown","pRoloc","knitr","BiocStyle","DESeq2","styler","utils","IRanges","BiocGenerics","rtracklayer","scuttle","txdbmaker","topGO","drawProteins","GenomicFeatures","biomaRt","AnnotationForge","Biostrings","GenomeInfoDb","SingleCellExperiment","SingleR","NOISeq","GenomicRanges","BSgenome")
    bioc_needed <- intersect(needed, bioc_pkgs)
    cran_needed <- setdiff(needed, bioc_pkgs)

    # Install CRAN packages first
    if (length(cran_needed) > 0) {
      message("Installing CRAN packages: ", paste(cran_needed, collapse = ", "))
      install.packages(cran_needed, dependencies = TRUE)
    }

    # Install Bioconductor packages
    if (length(bioc_needed) > 0) {
      message("Installing Bioconductor packages: ", paste(bioc_needed, collapse = ", "))
      BiocManager::install(bioc_needed)
    }

    message("All dependencies installed successfully")
    }
