#' Get pre-configured tables from OrthoMCL-DB
#'
#' A convenience function to quickly fetch preconfigured table of from OrthoMCL-DB. This function can be used to map old OG IDs to new and vice-versa given a vector of OG IDs. If no IDs are provided, entire data tables fetched from the database will be returned, which can be manipulated or processed further as per user need.
#'
#' @import httr2 jsonlite dplyr
#' @export
#'
#' @param og_ids Orthogroup IDs.
#' @param customField Preconfigured tables to be downloaded from OrthoMCL. 
#' @param timeout Timeout time in seconds
#' @param max_tries Maximum retries allowed.
#' 
#' \itemize{
#' \item \emph{"TaxonCounts"}: Phlyetic Distribution of Proteins
#' \item \emph{"PhyleticDistributionCounts"}: Phyletic Distribution clade abundance
#' \item \emph{"previousGroups"}: All previous groups."Default"
#' \item \emph{"GroupStat"}: Group Statistics
#' \item \emph{"SimilarGroup"}: Similar Groups
#' \item \emph{"EcNumber"}: Summary of EC Numbers
#' \item \emph{"PFams"}: Summary of Pfam domains
#' \item \emph{"ProteinPFams"}: PFam Architecture of Each Protein
#' \item \emph{"KeywordFrequency"}: Keyword Frequency
#'}
#'
#' @return A data frame.
#' @examples
#' \dontrun{
#' df <- getPreconfiguredTableOrthomcl("OG3_10277")
#' df <- getPreconfiguredTableOrthomcl(customField="SimilarGroup")
#' }
#'

getPreconfiguredTableOrthomcl <- function(og_ids = NULL,
                                          customField = "previousGroups",
                                          timeout = 60,
                                          max_tries = 3) {
  base_url <- paste0(
    "https://orthomcl.org/orthomcl/",
    "service/record-types/group/searches/AllGroups/reports/tableTabular"
  )
  
  reportConfig <- jsonlite::toJSON(
    list(
      tables = list(customField),
      includeHeader = TRUE,
      attachmentType = "csv"
    ),
    auto_unbox = TRUE
  )
  
  resp <- httr2::request(base_url) |>
    httr2::req_url_query(
      reportConfig = as.character(reportConfig)
    ) |>
    httr2::req_user_agent("plasmoRUtils getPreconfiguredTableOrthomcl") |>
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
  
  txt <- httr2::resp_body_string(resp)
  
  if (httr2::resp_status(resp) >= 400) {
    stop(
      "OrthoMCL request failed with HTTP ",
      httr2::resp_status(resp),
      "\n\nServer response:\n",
      txt,
      call. = FALSE
    )
  }
  
  res <- readr::read_csv(
    I(txt),
    progress = FALSE,
    show_col_types = FALSE
  )
  
  if (
    is.null(og_ids) ||
    length(og_ids) == 0 ||
    all(is.na(og_ids)) ||
    all(trimws(as.character(og_ids)) == "")
  ) {
    return(res)
  }
  
  og_ids <- trimws(as.character(og_ids))
  og_ids <- og_ids[!is.na(og_ids) & og_ids != ""]
  
  hits <- res |>
    dplyr::filter(
      dplyr::if_any(
        dplyr::everything(),
        ~ as.character(.x) %in% og_ids
      )
    )
  
  return(hits)
}
