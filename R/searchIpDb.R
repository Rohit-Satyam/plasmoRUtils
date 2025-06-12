#' Fetch orthologs from InParanoiDb
#'
#' This function retrieves the Plasmodium gene orthologs from InParanoiDB 9 database across all life forms
#'
#' @importFrom ggsci scale_fill_nejm
#' @import dplyr purrr
#' @import biomaRt
#' @importFrom glue glue
#' @import rvest
#' @export
#'
#' @param geneID Gene ID of Plasmodium falciparum.
#' @param ... Additional arguments that can be passed to \code{toGeneid} function.
#' @examples
#' \dontrun{
#' ids <- c("PF3D7_0807800", "PF3D7_1023900")
#' df <- searchIpDb(ids)
#' }
#'
searchIpDb <- function(geneID, ...) {
  pfa_ensembl <- toGeneid(geneID, from = "ensembl", ...)
  converted <- pfa_ensembl %>%
    dplyr::select(c(`Gene ID`, `UniProt ID(s)`)) %>%
    tidyr::separate_rows(`UniProt ID(s)`, sep = ",") %>%
    dplyr::rename(inputid = `Gene ID`, uniprotid = `UniProt ID(s)`)

  failed <- dplyr::setdiff(geneID, converted$inputid)
  if (length(failed) > 0) message(glue::glue("Following genes failed to convert to Uniprot ID: {failed}\n"))

  if (length(converted$inputid) == 0) {
    message("No matching input IDs found. Exiting function.")
    return(NULL)
  }

  res.list <- purrr::map(converted$uniprotid, function(protid) {
    url <- glue::glue("https://inparanoidb.sbc.su.se/orthologs/{protid}&1/")
    webpage <- tryCatch(
      {
        rvest::read_html(url)
      },
      error = function(e) {
        message("Error reading URL: ", url, "; Error message: ", e$message)
        return(NULL)
      }
    )

    if (is.null(webpage)) {
      message("No results were fetched\n")
      return(NULL)
    } else {
      df <- webpage %>%
        rvest::html_nodes("table") %>%
        rvest::html_table(fill = TRUE) %>%
        plyr::ldply(.) %>%
        dplyr::mutate(queryid = converted[converted$uniprotid == protid[1], ]$inputid) %>%
        dplyr::select(where(~ !all(is.na(.))))
      cat(paste("success", protid, "\n"))
      return(df)
    }
  }) %>%
    purrr::keep(~ !is.null(.x)) %>%
    purrr::map(~ mutate(.x, `Seed Score info_outline` = as.character(`Seed Score info_outline`))) %>%
    dplyr::bind_rows()
}

# styler:::style_active_file()
