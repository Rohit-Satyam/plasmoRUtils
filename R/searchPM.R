#' Fetch articles from PubMed database
#'
#' This function searches the PubMed corpus for the articles that contains your Gene ID of interest.
#'
#' @import dplyr purrr stringr
#' @import easyPubMed
#' @export
#'
#' @param geneID Character vector of Gene IDs.
#' @param org Scientific name of the organism. Default \emph{Plasmodium falciparum}
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

searchPM <- function(geneID,
                     org = "Plasmodium falciparum",
                     query = NULL,
                     from = 2010,
                     to = 2026,
                     verbose = TRUE) {
  results <- purrr::map(geneID, function(gid) {
    search_term <- if (is.null(query)) {
      paste0(
        '("', org, '"[All Fields] AND "',
        gid,
        '"[tiab:~0]) AND (',
        from,
        '[PDAT] : ',
        to,
        '[PDAT])'
      )
    } else {
      query
    }
    
    my_entrez_id <- easyPubMed::epm_query(
      query_string = search_term,
      verbose = FALSE
    )
    
    meta <- easyPubMed::get_epm_meta(my_entrez_id)
    
    record_count <- meta$count %||%
      meta$exp_count %||%
      length(easyPubMed::get_epm_uilist(my_entrez_id))
    
    if (as.numeric(record_count) > 0) {
      epm_obj <- easyPubMed::epm_fetch(
        my_entrez_id,
        format = "xml",
        verbose = FALSE
      )
      
      epm_obj <- easyPubMed::epm_parse(
        epm_obj,
        compact_output = TRUE,
        include_abstract = TRUE,
        max_authors = 1,
        verbose = FALSE
      )
      
      final_df <- easyPubMed::get_epm_data(epm_obj)
      
      ## Keep original positional selection as much as possible,
      ## but avoid selecting positions that do not exist
      keep_pos <- intersect(c(1, 2, 3, 5, 6, 7, 8, 9), seq_along(final_df))
      
      partialdf <- final_df %>%
        dplyr::select(dplyr::all_of(names(final_df)[keep_pos]))
      
      ## Ensure a usable title column exists
      title_col <- names(final_df)[
        stringr::str_detect(names(final_df), stringr::regex("title", ignore_case = TRUE))
      ][1]
      
      if (!is.na(title_col) && !"title" %in% names(partialdf)) {
        partialdf$title <- final_df[[title_col]]
      }
      
      if (!"title" %in% names(partialdf)) {
        stop(
          "No title-like column found in PubMed results. Available columns are: ",
          paste(names(final_df), collapse = ", "),
          call. = FALSE
        )
      }
      
      partialdf <- partialdf %>%
        dplyr::mutate(GeneID = gid)
      
      if (verbose) {
        query_translation <- meta$query_translation %||%
          meta$QueryTranslation %||%
          search_term
        
        cat(
          paste(
            "PubMed Query used for",
            gid,
            "was: \n",
            query_translation,
            "\n"
          )
        )
      }
      
      return(partialdf)
    } else {
      message("No results found for gene ID: ", gid)
      return(NULL)
    }
  })
  
  combined <- dplyr::bind_rows(results)
  
  if (nrow(combined) == 0) {
    return(combined)
  }
  
  combined <- combined %>%
    dplyr::distinct(.data$title, .keep_all = TRUE) %>%
    dplyr::mutate(
      title = stringr::str_trim(gsub("<.*?>", " ", .data$title))
    )
  
  return(combined)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    y
  } else {
    x
  }
}
