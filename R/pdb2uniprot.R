#' Quick mapping of PDB IDs to Uiprot IDs
#'
#' A convenience function to quickly convert the PDB IDs to Uniprot IDs. If a protein is multimeric, corresponding Uniprot IDs are returned for them. This function uses PDBe API.
#'
#' @import httr2 jsonlite dplyr tidyr purrr
#' @export
#'
#' @param pdbid A single PDB Id.
#' @param timeout Timeout time in seconds
#' @param max_tries Maximum retries allowed.
#'
#'
#' @return A data frame, Uniprot IDs, Chain IDs and start and end coordinates of the chains.
#' @examples
#' \dontrun{
#' df <- pdb2uniprot("9FIA")
#' }
#'

pdb2uniprot <- function(pdbid,
                        timeout = 60,
                        max_tries = 3) {
  url <- paste0(
    "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/",
    trimws(pdbid)
  )
  
  response <- httr2::request(url) |>
    httr2::req_user_agent("plasmoRUtils pdb2uniprot") |>
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
  
  data <- suppressMessages(httr2::resp_body_string(response))
  
  if (httr2::resp_status(response) >= 400) {
    stop(
      "PDBe request failed with HTTP ",
      httr2::resp_status(response),
      "\n\nServer response:\n",
      data,
      call. = FALSE
    )
  }
  
  json_data <- jsonlite::fromJSON(data)
  attributes <- json_data[[1]][[1]]
  
  if (length(attributes) > 1) {
    df <- attributes %>%
      purrr::map_dfr(~ {
        nested_cols <- .x$mappings %>%
          dplyr::select(where(~ is.data.frame(.x) || (is.list(.x) && all(purrr::map_lgl(.x, is.data.frame))))) %>%
          names()
        
        .x$mappings %>%
          tidyr::unnest_wider(any_of(nested_cols), names_sep = ".")
      }, .id = "attribute")
    
    ## Correcting chain names mistaken as NA
    struct_cols <- grep("struct_asym_id", colnames(df), value = TRUE)
    
    # Replace NA with "NA" in all matching columns
    df <- df %>%
      dplyr::mutate(across(all_of(struct_cols), ~ tidyr::replace_na(., "NA")))
    
    df$query <- pdbid
    
    return(df)
  } else {
    df <- purrr::map_dfr(attributes, ~ as.data.frame(.x$mappings), .id = "attribute")
    
    ## Identifying nested columns and flattening them
    nested_cols <- df %>%
      dplyr::select(where(~ is.data.frame(.x) || (is.list(.x) && all(purrr::map_lgl(.x, is.data.frame))))) %>%
      names()
    
    ## Flattening the columns
    df <- df %>%
      tidyr::unnest_wider(any_of(nested_cols), names_sep = ".")
    
    ## Correcting chain names mistaken as NA
    struct_cols <- grep("struct_asym_id", colnames(df), value = TRUE)
    
    # Replace NA with "NA" in all matching columns
    df <- df %>%
      dplyr::mutate(across(all_of(struct_cols), ~ tidyr::replace_na(., "NA")))
    
    df$query <- pdbid
    
    return(df)
  }
}
