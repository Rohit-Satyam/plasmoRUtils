#' Get effective gene lengths form TPM calculation
#'
#' This function provides effective length (sum of lengths of exons) of the genes for calculating TPM values.
#'
#' @import dplyr
#' @import GenomicFeatures
#' @importFrom BiocGenerics width
#' @importFrom IRanges reduce
#' @export
#'
#' @param gtf Provide path or URL to GTF file.
#' @param format Format of the feature file i.e. "gtf" or "gff3". Default "gff3".
#'
#' @return df numeric. This function returns a data frame with 2 columns: "GeneID", "Length".
#' @examples
#' \dontrun{
#' baseurl <- "https://plasmodb.org/common/downloads/release-68/"
#' getEffLen(paste0(baseurl,
#' "Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff"))
#'
#' OR
#'
#' getEffLen("/data/PlasmoDB-67_Pfalciparum3D7.gtf")
#' }
getEffLen <- function(gtf = NULL, format = "gff3") {
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = format)
  exonic <- GenomicFeatures::exonsBy(txdb, by = "gene")
  red_exonic <- IRanges::reduce(exonic)

  gene_lengths <- red_exonic %>%
    BiocGenerics::width() %>%
    vapply(sum, numeric(1)) %>%
    tibble::enframe(name = "GeneID", value = "Length")

  return(gene_lengths)
  ## source: https://www.biostars.org/p/185665/
}

#styler:::style_active_file()
