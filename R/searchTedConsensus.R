#' Fetch protein domains from TED database
#'
#' A convenience function to quickly access The Encyclopedia of Domains (TED) database and fetch domain boundary information for given Uniprot IDs. For information about the column names the users are requested to refer to the TED database at https://ted.cathdb.info/.
#'
#' @import dplyr purrr httr2
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @param uniprotid A character vector of uniprot IDs.
#' @param returnCATHdesc Logical. Set this on to get the description of the CATH ID from CATH database.
#' @param timeout Timeout time in seconds
#' @param max_tries Maximum retries allowed.
#'
#' @return A data frame, domain boundaries and other information provided by TED. For details visit the TED database.
#' @examples
#' \dontrun{
#' df <- searchTedConsensus(
#' c("Q7K6A1","Q8IAP8","C0H4D0","C6KT90","Q8IBJ7"),
#' returnCATHdesc=FALSE)
#' }
#'

searchTedConsensus <- function(uniprotid = "",
                               returnCATHdesc = TRUE,
                               timeout = 60,
                               max_tries = 3) {
  fetch_response <- function(url) {
    httr2::request(url) |>
      httr2::req_user_agent("plasmoRUtils searchTedConsensus") |>
      httr2::req_timeout(timeout) |>
      httr2::req_retry(
        max_tries = max_tries,
        retry_on_failure = TRUE,
        is_transient = function(resp) {
          httr2::resp_status(resp) %in% c(408, 429, 500, 502, 503, 504)
        }
      ) |>
      httr2::req_error(is_error = function(resp) FALSE) |>
      httr2::req_perform()
  }
  
  urls <- paste0(
    "https://ted.cathdb.info/api/v1/uniprot/summary/",
    trimws(uniprotid)
  )
  
  urls <- urls %>%
    purrr::map(fetch_response) %>%
    purrr::map2_chr(urls, ~ {
      code <- httr2::resp_status(.x)
      if (code < 400 && code != 404) .y else NA_character_
    }) %>%
    purrr::discard(is.na)
  
  res <- urls %>%
    purrr::map(fetch_response) %>%
    purrr::map(~ httr2::resp_body_string(.x)) %>%
    purrr::map(jsonlite::fromJSON) %>%
    purrr::map_dfr(~ if (.x$count > 0) .x$data else NULL)
  
  ## removing the last column which is list
  res <- res[, -ncol(res)]
  
  if (returnCATHdesc) {
    desc <- lapply(res$cath_label, function(x) {
      if (x != "-") {
        page_url <- paste0(
          "https://www.cathdb.info//version/latest/cathnode/",
          x
        )
        
        page_txt <- fetch_response(page_url) %>%
          httr2::resp_body_string()
        
        page <- rvest::read_html(page_txt)
        
        page %>%
          rvest::html_elements("h2") %>%
          rvest::html_text() %>%
          .[1]
      } else {
        return(NULL)
      }
    }) %>% as.character()
    
    res$cath_label_desc <- desc
  }
  
  absentIDs <- setdiff(uniprotid, unique(res$uniprot_acc))
  
  if (length(absentIDs) > 0) {
    message(
      paste0(
        "Following Uniprot IDs returned no results: ",
        paste(absentIDs, collapse = " "),
        "\n"
      )
    )
  }
  
  return(res)
}