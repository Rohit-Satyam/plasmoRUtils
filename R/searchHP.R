#' Fetch Protein-protein interaction for given gene IDs from Hit Predict database
#'
#' This function searches the Hitpredict database to retrieve the Experimental Protein-Protein Interaction data.
#'
#' @import dplyr
#' @importFrom biomaRt useEnsemblGenomes getBM
#' @importFrom S4Vectors merge
#' @importFrom readr read_lines
#' @export
#'
#' @param geneID Single gene ID.
#' @param taxid Taxon ID of the organism of interest. Default: 36329.
#' @param uniprotToGID To convert Uniprot ID to gene ID. Set TRUE for Plasmodium geneIDs only.
#'
#' @return A data frame, containing 11 columns: "Interaction", "Interactor", "Name", "Experiments", "Category", "Method.Score", "Annotation.Score", "Interaction.Score", "Confidence", "QueryID", "ensembl_gene_id".
#' @examples
#' \dontrun{
#' test <- searchHP("PF3D7_0418300")
#'
#' ## To use it for other organism, turn off uniprotToGID and provide taxid of the organism
#' test <- searchHP("BRCA1",taxid = "9606" , uniprotToGID = FALSE)
#' }
#'

searchHP <- function(geneID, taxid="36329", uniprotToGID=TRUE) {
  url <- glue::glue("http://www.hitpredict.org/proteins.php?Value={geneID}&Species={taxid}")

  webpage <- rvest::read_html(url)


  htp_link <- webpage %>%
    rvest::html_nodes("a") %>%
    rvest::html_attr("href") %>%
    stringr::str_subset("./htp_int") %>%
    unique()

  if (length(htp_link)!=0) {
    alldf <- purrr::map(htp_link, function(x){
      url <- paste0("http://www.hitpredict.org/", stringr::str_replace(x, "./htp_int", "htp_int_txt"))
      skip_lines_start <- 3
      skip_lines_end <- 1
      total_lines <- readr::read_lines(url) %>% length()
      nrows <- total_lines - skip_lines_start - skip_lines_end - 1
      message(glue::glue("\033[0;32mPPI found for: {geneID}\033[0m\n"))
      data <- utils::read.table(url, header = TRUE, sep = "\t", skip = skip_lines_start, nrows = nrows) %>%
        dplyr::mutate(QueryID = geneID)
    })

    data <- dplyr::bind_rows(alldf)

    if (uniprotToGID) {
      pfa_ensembl <- biomaRt::useEnsemblGenomes(biomart = "protists_mart", dataset = "pfalciparum_eg_gene")
      converted <- biomaRt::getBM(attributes = c("uniprotsptrembl", "ensembl_gene_id"),
                         filters = "uniprotsptrembl",
                         values = data$Interactor,
                         mart = pfa_ensembl)

      data <- S4Vectors::merge(data,converted, all.x=TRUE, by.x="Interactor", by.y="uniprotsptrembl")
    }

    return(data)
  } else {
    message("No interaction found in HitPredict Database")
  }
}
