#' Fetch protein domains from TED database
#'
#' A convenience function to quickly access The Encyclopedia of Domains (TED) database and fetch domain boundary information for given Uniprot IDs. For information about the column names the users are requested to refer to the TED database at https://ted.cathdb.info/.
#'
#' @import dplyr purrr
#' @importFrom httr GET content
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @param uniprotid A character vector of uniprot IDs.
#' @param returnCATHdesc Logical. Set this on to get the description of the CATH ID from CATH database.
#'
#' @return A data frame, domain boundaries and other information provided by TED. For details visit the TED database.
#' @examples
#' \dontrun{
#' df <- searchTedConsensus(
#' c("Q7K6A1","Q8IAP8","C0H4D0","C6KT90","Q8IBJ7"),
#' returnCATHdesc=FALSE)
#' }
#'



searchTedConsensus <- function(uniprotid="", returnCATHdesc=TRUE) {
  urls <- paste0("https://ted.cathdb.info/api/v1/uniprot/summary/", trimws(uniprotid))
  res <- urls %>%
    purrr::map(httr::GET) %>%
    purrr::map(~ httr::content(.x, "text", encoding = "UTF-8")) %>%
    purrr::map(jsonlite::fromJSON) %>%
    purrr::map_dfr(~if (.x$count > 0) .x$data else NULL)

  ## removing the last column which is list
  res <- res[, -ncol(res)]
  if(returnCATHdesc){
    desc <- lapply(res$cath_label, function(x){
      if(x!="-"){
        page <- rvest::read_html(paste0("https://www.cathdb.info//version/latest/cathnode/",x))
        page  %>% html_elements("h2") %>% html_text() %>% .[1]
      } else{
        return(NULL)
      }

    }) %>% as.character()
    res$cath_label_desc <- desc
  }
  absentIDs <- setdiff(uniprotid, unique(res$uniprot_acc))

  if (length(absentIDs) > 0) {
    message(glue::glue("Following Uniprot IDs returned no results: {paste(absentIDs, collapse = ' ')}\n"))
  }
  return(res)
}



