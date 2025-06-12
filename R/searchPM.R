#' Fetch articles from PubMed database
#'
#' This function searches the PubMed corpus for the articles that contains your Gene ID of interest.
#'
#' @import dplyr
#' @import easyPubMed
#' @export
#'
#' @param geneID Character vector of Gene IDs.
#' @param org Scientific name of the organism. Default "Plasmodium falciparum"
#' @param query String of user defined custom queries. If you wish to pass your own combination od terms use this argument.
#' @param from To define the start year for querying the articles.
#' @param to To define the end year for querying the articles.
#' @param verbose Disable to turn off the messages printed by the function.
#' @return A data frame, containing 9 columns: "pmid"    "doi"     "title"   "year"    "month"   "day"     "jabbrv"  "journal" "GeneID" .
#' @examples
#' \dontrun{
#' test <- searchPM(geneID = c("PF3D7_0420300","PF3D7_0621000"))
#' }
#'
searchPM <- function(geneID, org = "Plasmodium falciparum", query = NULL, from = 2010, to = 2025, verbose = TRUE) {
  results <- purrr::map(geneID, function(gid) {
    search_term <- if (is.null(query)) {
      paste0('("', org, '"[All Fields] AND "', gid, '"[tiab:~0]) AND (', from, '[PDAT] : ', to, '[PDAT])')
    } else {
      query
    }

    my_entrez_id <- easyPubMed::get_pubmed_ids(search_term)
    if (as.numeric(my_entrez_id$Count) > 0) {
      xml <- easyPubMed::articles_to_list(easyPubMed::fetch_pubmed_data(my_entrez_id, format = "xml"))
      final_df <- purrr::map_dfr(xml, article_to_df, max_chars = -1, getAuthors = FALSE)
      partialdf <- final_df %>% dplyr::select(1, 2, 3, 5, 6, 7, 8, 9) %>% dplyr::mutate(GeneID = gid)
      if (verbose) cat(paste("PubMed Query used for", gid, "was: \n", my_entrez_id$QueryTranslation, "\n"))
      return(partialdf)
    } else {
      message("No results found for gene ID: ", gid)
      return(NULL)
    }
  })

  combined <- dplyr::bind_rows(results) %>%
    dplyr::distinct(.data$title, .keep_all = TRUE) %>%
    dplyr::mutate(title = stringr::str_trim(gsub("<.*?>", " ", .data$title)))

  return(combined)
}
