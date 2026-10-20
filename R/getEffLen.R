#' Get effective gene lengths form TPM calculation
#'
#' This function provides effective length (sum of lengths of exons) of the genes for calculating TPM values. Calculates the union of all annotated exons across every transcript associated with a gene.
#'
#' @import dplyr
#' @export
#'
#' @param gtf Provide path or URL to GTF file.
#' @param format Format of the feature file i.e. "gtf" or "gff3". Default "gff3".
#'
#' @return df numeric. This function returns a data frame with 2 columns: "GeneID", "Length".
#' @examples
#' \dontrun{
#'
#' getEffLen("/data/PlasmoDB-67_Pfalciparum3D7.gtf")
#' }
getEffLen <- function(gtf = NULL, format = "gff3") {
  pkgs <- c("GenomicFeatures", "IRanges", "BiocGenerics")

  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(
      "Required package(s) missing: ",
      paste(missing_pkgs, collapse = ", "),
      ". Please install them using BiocManager::install().",
      call. = FALSE
    )
  }

  if (length(gtf) != 1L || is.na(gtf) || !nzchar(gtf)) stop("Provide one GFF3/GTF path")


  txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = format)
  exonic <- GenomicFeatures::exonsBy(txdb, by = "gene")
  red_exonic <- IRanges::reduce(exonic)

  gene_lengths <- red_exonic %>%
    BiocGenerics::width() %>%
    vapply(sum, numeric(1)) %>%
    tibble::enframe(name = "GeneID", value = "Length")

  return(gene_lengths)
  ## source: https://www.biostars.org/p/185665/
}

# styler:::style_active_file()
