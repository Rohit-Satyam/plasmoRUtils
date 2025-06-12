#' Calculate TPM values from count data
#'
#' This function provides ability to compute quick TPM values.
#'
#' @import dplyr tibble
#' @export
#'
#' @param counts Count matrix containing raw counts.
#' @param featureLength Effective length of genes generated from \code{getEfflen}.
#'
#' @return df numeric. This function returns a dataframe of TPM normalized counts and final column containing feature length.
#' @examples
#' \dontrun{
#' ## Effective length of the gene
#' gene_info <- data.frame(GeneID = c("gene1", "gene2", "gene3"), Length = c(1000, 1500, 2000))
#' count_matrix <- matrix(c(10, 20, 30, 40, 50, 60),
#'   nrow = 3, ncol = 2,
#'   dimnames = list(c("gene1", "gene2", "gene3"), c("sample1", "sample2"))
#' )
#' test <- easyTPM(count_matrix, gene_info)
#' }
#'
easyTPM <- function(counts, featureLength) {
  effLen <- featureLength %>% dplyr::filter(GeneID %in% rownames(counts))

  rpk <- counts / (effLen$Length / 1000)
  scaling_factor <- colSums(rpk)

  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6
  tpm_df <- as.data.frame(tpm)

  tpm_df <- tpm_df %>%
    tibble::rownames_to_column("GeneID") %>%
    dplyr::left_join(effLen, by = "GeneID") %>%
    tibble::column_to_rownames("GeneID")

  return(tpm_df)
}

#styler:::style_active_file()
