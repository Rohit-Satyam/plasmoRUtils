#' Calculate transcripts per million (TPM) from counts
#'
#' Normalise a raw count matrix using matched gene lengths. Gene lengths are aligned by gene identifier before calculating reads per kilobase and scaling each sample to one million TPM.
#'
#' @import dplyr tibble
#' @export
#'
#' @param counts Count matrix containing raw counts.
#' @param featureLength Effective length of genes generated from [getEffLen()].
#'
#' @return df numeric. This function returns a dataframe of TPM normalized counts and final column containing feature length.
#' @examples
#' \dontrun{
#' ## Effective length of the gene
#' gene_info <- data.frame(GeneID = c("gene1", "gene2", "gene3"), Length = c(1000, 1500, 2000))
#' count_matrix <- matrix(c(10, 20, 30, 40, 50, 60),
#'   nrow = 3, ncol = 2,
#'   dimnames = list(c("gene2", "gene1", "gene3"), c("sample1", "sample2"))
#' )
#' test <- easyTPM(count_matrix, gene_info)
#' }
#'
easyTPM <- function(counts, featureLength) {
  ## Exception handling
  if (is.null(dim(counts)) || length(dim(counts)) != 2L || !is.numeric(counts)) {
    stop("counts must be a numeric matrix or matrix-like object")
  }
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) ||
    anyDuplicated(rownames(counts)) || any(!nzchar(rownames(counts)))) {
    stop("counts must have unique, nonempty gene row names")
  }
  if (!all(c("GeneID", "Length") %in% names(featureLength))) {
    stop("featureLength must contain GeneID and Length")
  }
  ids <- as.character(featureLength$GeneID)
  lens <- as.numeric(featureLength$Length)
  if (anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) {
    stop("featureLength has missing or duplicate gene identifiers")
  }
  idx <- match(rownames(counts), ids) # aligns lengths, not just membership
  if (anyNA(idx)) stop("No length for: ", paste(rownames(counts)[is.na(idx)], collapse = ", "))
  lens <- lens[idx] ## fixed ordering of gene lengths
  if (any(!is.finite(lens) | lens <= 0)) stop("Lengths must be positive and finite")
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("Counts must be finite and nonnegative")
  }

  rpk <- counts / (lens / 1000) ## fixed ordering of gene lengths
  scaling_factor <- colSums(rpk)

  if (any(!is.finite(scaling_factor) | scaling_factor == 0)) {
    stop("At least one sample has zero/invalid total RPK")
  }

  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6
  tpm_df <- as.data.frame(tpm)

  tpm_df <- tpm_df %>%
    tibble::rownames_to_column("GeneID") %>%
    dplyr::left_join(featureLength, by = "GeneID") %>%
    tibble::column_to_rownames("GeneID")

  return(tpm_df)
}

# styler:::style_active_file()
# styler::style_dir("R", filetype = "R")