
#' Retrieve gene annotations using a mixed list of IDs
#'
#' Retrieves gene annotations across all databases integrated into
#' VEuPathDB, without requiring organism or database information.
#'
#' Uses the VEuPathDB GeneByLocusTag search.
#'
#' @import httr2 readr jsonlite
#' @export
#'
#' @param ids Character vector of gene identifiers.
#' @param customFields Optional character vector of attributes.
#'   If NULL, default mapping attributes are returned.
#'   See getTableAttributes() for available attributes.
#' @param api_key VEuPathDB API key.
#'
#' @return A data frame containing the input IDs, matched gene IDs,
#'   organism information and requested annotations.
#'
#' @examples
#' \dontrun{
#' df <- getTableVeupathdb(
#'   ids = c("PF3D7_1372300", "PBANKA_1355600","A0A143ZZR8","W7JSP8" )
#' )
#' }

getTableVeupathdb <- function(
    ids,
    customFields = NULL,
    api_key = Sys.getenv("VEUPATHDB_API_KEY")
) {
  
  # Validate API key
  if (!is.character(api_key) ||
      length(api_key) != 1L ||
      is.na(api_key) ||
      !nzchar(trimws(api_key))) {
    stop("VEuPathDB API key is missing.", call. = FALSE)
  }
  
  # Validate input IDs
  if (!is.character(ids) || anyNA(ids)) {
    stop("ids must be a character vector without NA values.")
  }
  
  ids <- unique(trimws(ids))
  ids <- ids[nzchar(ids)]
  
  if (!length(ids)) {
    stop("No valid input IDs supplied.")
  }
  
  # Integrated VEuPathDB endpoint
  base_url <- "https://veupathdb.org/veupathdb/service"
  
  # Preserve the session across requests
  cookie_file <- tempfile()
  on.exit(unlink(cookie_file), add = TRUE)
  
  # STEP 1: Create an input-ID dataset
  response <- httr2::request(
    paste0(base_url, "/users/current/datasets")
  ) |>
    httr2::req_auth_bearer_token(api_key) |>
    httr2::req_cookie_preserve(cookie_file) |>
    httr2::req_body_json(
      list(
        sourceType = "idList",
        sourceContent = list(
          ids = as.list(ids)
        )
      ),
      auto_unbox = TRUE
    ) |>
    httr2::req_perform()
  
  dataset <- httr2::resp_body_json(response)
  
  dataset_id <- dataset$id
  
  if (is.null(dataset_id) || length(dataset_id) != 1L ||
      is.na(dataset_id)) {
    stop("VEuPathDB did not return a valid dataset ID.")
  }
  
  # STEP 2: Configure attributes
  if (is.null(customFields)) {
    
    fields <- c(
      "primary_key",
      "source_id",
      "dataset_order",
      "organism",
      "gene_product",
      "input_id",
      "uniprot_ids",
      "gene_entrez_id",
      "gene_previous_ids",
      "gene_name",
      "strand"
    )
    
  } else {
    
    if (!is.character(customFields) ||
        !length(customFields) ||
        anyNA(customFields) ||
        any(!nzchar(trimws(customFields)))) {
      stop("customFields must be a non-empty character vector.")
    }
    
    fields <- unique(c(
      "primary_key",
      "input_id",
      "dataset_order",
      "organism",
      customFields
    ))
  }
  
  # Build report configuration
  report_config <- jsonlite::toJSON(
    list(
      attributes = fields,
      includeHeader = TRUE,
      attachmentType = "plain"
    ),
    auto_unbox = TRUE
  )
  
  # STEP 3: Retrieve the annotation table
  url <- paste0(
    base_url,
    "/record-types/transcript/searches/",
    "GeneByLocusTag/reports/attributesTabular"
  )
  
  response <- httr2::request(url) |>
    httr2::req_auth_bearer_token(api_key) |>
    httr2::req_cookie_preserve(cookie_file) |>
    httr2::req_headers(Accept = "text/plain") |>
    httr2::req_url_query(
      ds_gene_ids = as.character(dataset_id),
      reportConfig = as.character(report_config)
    ) |>
    httr2::req_perform()
  
  # STEP 4: Parse the returned TSV
  out <- readr::read_tsv(
    I(httr2::resp_body_string(response)),
    col_types = readr::cols(
      .default = readr::col_character()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  
  # Remove the prefix from previous IDs, when requested
  if ("Previous ID(s)" %in% names(out)) {
    out$`Previous ID(s)` <- sub(
      "^Previous IDs: ",
      "",
      out$`Previous ID(s)`
    )
  }
  
  return(out)
}
