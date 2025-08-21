#' Fetch orthologs from InParanoiDb
#'
#' This function retrieves the orthologs from [InParanoiDB 9](https://inparanoidb.sbc.su.se/) database for 640 species given set of Ensembl Gene IDs or uniprot IDs usinf database API. 
#' 
#' To view list of species covered by InParanoiDB 9, use \code{listipdb()} function.   
#'
#' @import dplyr purrr
#' @importFrom tidyr separate_rows
#' @importFrom readr read_csv
#' @importFrom glue glue
#' @export
#'
#' @param geneID Gene ID of \emph{Plasmodium falciparum} or VEupathDB enlisted organisms that are also covered by InParanoiDB. If providing uniprot ID, set idtype="uniprot" to prevent ID conversion.
#' @param idtype Set this to "uniprot" if using uniprot IDs directly and GenID to uniprot ID conversion is not required.
#' @param ... Additional arguments that can be passed to \code{toGeneid()} function. This comes in very handy if you are working with parasite gene IDs other than \emph{Plasmodium}.
#' 
#' @return A data frame, containing 10 columns.
#'  
#' @seealso
#'  \code{\link[plasmoRUtils]{listipdb}}
#'  
#' @examples
#' \dontrun{
#' df <-  searchIpDb( c("PF3D7_0807800", "PF3D7_1023900"))
#' df <- searchIpDb( c("C5LD32", "A5KAC7"),idtype = "uniprot")
#' }
#'
searchIpDb <- function(geneID, ..., idtype = "ensembl") {
  if (idtype == "ensembl") {
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
  } else {
    converted <- data.frame(uniprotid = geneID, inputid = geneID)
  }

  res.list <- purrr::map(converted$uniprotid, function(protid) {
    url <- glue::glue("https://inparanoidb.sbc.su.se/download/proteinorthologs/{protid}&all&csv")
    webpage <- tryCatch(
      {
        readr::read_csv(url, progress = FALSE, show_col_types = FALSE)
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
        dplyr::mutate(queryid = converted[converted$uniprotid == protid[1], ]$inputid)
      cat(paste("success", protid, "\n"))
      return(df)
    }
  }) %>%
    purrr::keep(~ !is.null(.x)) %>%
    # purrr::map(~ mutate(.x, `Seed Score info_outline` = as.character(`Seed Score info_outline`))) %>%
    dplyr::bind_rows()
}

#styler:::style_active_file()
