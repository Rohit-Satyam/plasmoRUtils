#' Quickly plot PhenoPlasm summary tables
#'
#' This function generates Disruptability and Mutant Phenotype tables in R, mirroring the style of Phenoplasm visualizations.
#'
#' @import mgsub gt
#' @importFrom plyr ldply
#' @importFrom reshape2 dcast
#' @importFrom stringr str_c str_trim
#' @importFrom stats setNames
#' @export
#'
#' @param file Path of Phenotype.txt file obtained from PhenoPlasm database or a data frame.
#' @param skip Number of lines to skip in the file. Default 2.
#'
#' @return A gt table plot.
#' @examples
#' \dontrun{
#'
#' ## Read the table generated from Phenoplasm advance search and pass the resulting data frame
#' df <- read.csv("phenotype.txt", skip = 2, sep = "\t") %>%
#' dplyr::select(-3, -4) %>%
#' dplyr::rename_with(~ gsub("Sprozoite", "Sporozoite", .x))
#' easyPhplplottbl(df)
#'
#' ## Pass the file path directly
#' easyPhplplottbl("phenotype.txt")
#'
#' ## Load example data frame
#' data(pf3d7PhplTable)
#'
#' easyPhplplottbl(pf3d7PhplTable)
#' }
#'

easyPhplplottbl <- function(file, skip=2) {

  # Define icon mappings
  gene <- data.frame(
    abbr = c("V", "R", "U", "D", "A", "E", "I", "T", "C", "S"),
    html = c(
      '<span style="color:#00AA00 !important;">&#10004;</span>',
      '<span style="color:darkred;">&#10060;</span>',
      '<span style="color:#00AA00 !important;">&#9989;</span>',
      '<span style="color:darkred;">&#10071;</span>',
      '<span style="color:darkred;">&#10071;</span>',
      '<span style="color:darkred;">&#11145;</span>',
      '<span style="color:darkred;">&#10228;</span>',
      '<span style="color:darkred;">&#128202;</span>',
      '<span style="color:darkred;">&#9728;</span>',
      '<span style="color:darkred;">&#128997;</span>'
    )
  )

  ortho <- data.frame(
    abbr = c("V", "R", "U", "D", "A", "E", "I", "T", "C", "S"),
    html = c(
      '<span style="color:darkgreen; opacity:0.4;">&#10004;</span>',
      '<span style="color:darkred; opacity:0.4;">&#10060;</span>',
      '<span style="color:darkgreen; opacity:0.4;">&#9989;</span>',
      '<span style="color:darkred; opacity:0.4;">&#10071;</span>',
      '<span style="color:darkred; opacity:0.4;">&#10071;</span>',
      '<span style="color:darkred; opacity:0.4;">&#11145;</span>',
      '<span style="color:darkred; opacity:0.4;">&#10228;</span>',
      '<span style="color:darkred; opacity:0.4;">&#128202;</span>',
      '<span style="color:darkred; opacity:0.4;">&#9728;</span>',
      '<span style="color:darkred; opacity:0.4;">&#128997;</span>'
    )
  )

  # Read and preprocess data
  if(is.character(file)){
    ph <- read.csv(file, skip = skip, sep = "\t") %>%
      dplyr::select(-3, -4) %>%
      dplyr::rename_with(~ gsub("Sprozoite", "Sporozoite", .x))
  } else {
    ph <- file
  }


  # Replace abbreviations with HTML icons
  gene_cols <- grep("Gene", colnames(ph), value = TRUE)
  ortho_cols <- grep("Ortho", colnames(ph), value = TRUE)

  ph <- ph %>%
    dplyr::mutate(across(all_of(gene_cols[2:length(gene_cols)]),
                  ~ mgsub::mgsub(.x, gene$abbr, gene$html))) %>%
    dplyr::mutate(across(all_of(ortho_cols),
                  ~ mgsub::mgsub(.x, ortho$abbr, ortho$html)))

  # Process columns
  col_unique <- colnames(ph[, 3:ncol(ph)]) %>%
    mgsub::mgsub(c("Gene", "Ortholog"), c("", "")) %>%
    unique()
  # Replace NAs
  ph[is.na(ph)] <- ""
  result_list <- lapply(col_unique, function(x) {
    gene_col <- paste0("Gene", x)
    ortholog_col <- paste0("Ortholog", x)

    ph %>%
      dplyr::select(Gene, all_of(c(gene_col, ortholog_col))) %>%
      dplyr::mutate(Gene_Ortholog = stringr::str_c(.[[2]], .[[3]], sep = " ") %>% stringr::str_trim()) %>%
      dplyr::select(Gene, Gene_Ortholog)
  }) %>% stats::setNames(col_unique)

  # Combine and reshape
  temp <- plyr::ldply(result_list)
  temp2 <- reshape2::dcast(temp, Gene ~ .id, value.var = "Gene_Ortholog")

  # Convert to HTML after dcast
  temp2 <- temp2 %>%
    mutate(across(all_of(col_unique), ~ lapply(.x, gt::html)))
  tbl <- gt::gt(temp2) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    )
  return(tbl)
}

