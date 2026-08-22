#' Get pre-configured tables from OrthoMCL-DB
#'
#' A convenience function to quickly fetch preconfigured tables from
#' OrthoMCL-DB. This function can be used to map old OG IDs to new and
#' vice-versa given a vector of OG IDs. If no IDs are provided, the entire
#' data table fetched from the database is returned and a copy of it is cached to
#' reduce the next query time.
#'
#' @import httr2 jsonlite dplyr
#' @export
#'
#' @param og_ids Orthogroup IDs.
#' @param customField Preconfigured table to be downloaded from OrthoMCL.
#' @param api_key Character. VEuPathDB API key used for authentication.
#'   By default, the value is obtained from the `VEUPATHDB_API_KEY`
#'   environment variable.
#' @param timeout Timeout time in seconds.
#' @param max_tries Maximum retries allowed.
#' @param refresh Logical. If `TRUE`, ignores any locally cached copy and
#'   downloads the requested table again from OrthoMCL. Defaults to `FALSE`.
#'   
#' \itemize{
#' \item \emph{"TaxonCounts"}: Phyletic Distribution of Proteins
#' \item \emph{"PhyleticDistributionCounts"}: Phyletic Distribution clade abundance
#' \item \emph{"previousGroups"}: All previous groups. "Default"
#' \item \emph{"GroupStat"}: Group Statistics
#' \item \emph{"SimilarGroup"}: Similar Groups
#' \item \emph{"EcNumber"}: Summary of EC Numbers
#' \item \emph{"PFams"}: Summary of Pfam domains
#' \item \emph{"ProteinPFams"}: PFam Architecture of Each Protein
#' \item \emph{"KeywordFrequency"}: Keyword Frequency
#' }
#'
#' @return A data frame.
#'
#' @examples
#' \dontrun{
#' df <- getPreconfiguredTableOrthomcl("OG3_10277")
#' df <- getPreconfiguredTableOrthomcl(customField = "SimilarGroup")
#' }
#'

getPreconfiguredTableOrthomcl <- function(
    og_ids = NULL,
    customField = "previousGroups",
    api_key = Sys.getenv("VEUPATHDB_API_KEY"),
    timeout = 60,
    max_tries = 3,
    refresh = FALSE
) {
  
  if (!nzchar(api_key)) {
    stop("VEuPathDB API key is missing.", call. = FALSE)
  }
  
  # Cache location
  cache_dir <- tools::R_user_dir(
    "plasmoRUtils",
    which = "cache"
  )
  
  dir.create(
    cache_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  cache_file <- file.path(
    cache_dir,
    paste0(
      "orthomcl_",
      customField,
      ".rds"
    )
  )
  
  # Use cached table unless refresh = TRUE
  if (file.exists(cache_file) && !refresh) {
    
    res <- readRDS(cache_file)
    
  } else {
    
    base_url <- paste0(
      "https://orthomcl.org/orthomcl/",
      "service/record-types/group/searches/",
      "AllGroups/reports/tableTabular"
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
      httr2::req_auth_bearer_token(api_key) |>
      httr2::req_url_query(
        reportConfig = as.character(reportConfig)
      ) |>
      httr2::req_user_agent(
        "plasmoRUtils getPreconfiguredTableOrthomcl"
      ) |>
      httr2::req_timeout(timeout) |>
      httr2::req_retry(
        max_tries = max_tries,
        retry_on_failure = TRUE,
        is_transient = function(resp) {
          httr2::resp_status(resp) %in%
            c(408, 429, 500, 502, 503, 504)
        }
      ) |>
      httr2::req_error(
        is_error = function(resp) FALSE
      ) |>
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
    
    # Save for subsequent calls
    saveRDS(res, cache_file)
  }
  
  # Return complete table if no OG IDs supplied
  if (
    is.null(og_ids) ||
    length(og_ids) == 0 ||
    all(is.na(og_ids)) ||
    all(trimws(as.character(og_ids)) == "")
  ) {
    return(res)
  }
  
  og_ids <- trimws(as.character(og_ids))
  og_ids <- og_ids[
    !is.na(og_ids) & og_ids != ""
  ]
  
  res |>
    dplyr::filter(
      dplyr::if_any(
        dplyr::everything(),
        ~ as.character(.x) %in% og_ids
      )
    )
}