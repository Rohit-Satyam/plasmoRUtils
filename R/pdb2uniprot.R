#' Quick mapping of PDB IDs to Uiprot IDs
#'
#' A convenience function to quickly convert the PDB IDs to Uniprot IDs. If a protein is multimeric, corresponding Uniprot IDs are returned for them. This function uses PDBe API.
#'
#' @import httr jsonlite dplyr
#' @importFrom purrr map_dfr
#' @export
#'
#' @param pdbid A single PDB Id.
#'
#'
#' @return A data frame, Uniprot IDs, Chain IDs and start and end coordinates of the chains.
#' @examples
#' \dontrun{
#' df <- pdb2uniprot("9FIA")
#' }
#'

pdb2uniprot <- function(pdbid) {
  url <- paste0("https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/", trimws(pdbid))

  response <- httr::GET(url)
  data <- suppressMessages(httr::content(response, "text"))
  json_data <- jsonlite::fromJSON(data)
  attributes <- json_data[[1]][[1]]

  if(length(attributes)>1){
    df <- attributes %>%
      map_dfr(~ {
        nested_cols <- .x$mappings %>%
          dplyr::select(where(~ is.data.frame(.x) || (is.list(.x) && all(purrr::map_lgl(.x, is.data.frame)) ))) %>%
          names()
        .x$mappings %>%
          unnest_wider(any_of(nested_cols), names_sep = ".")
      },.id="attribute")

    ## Correcting chain names mistaken as NA
    struct_cols <- grep("struct_asym_id", colnames(df), value = TRUE)

    # Replace NA with "NA" in all matching columns
    df <- df %>%
      dplyr::mutate(across(all_of(struct_cols), ~ tidyr::replace_na(., "NA")))
    df$query <- pdbid

    return(df)
  } else {
    df <-purrr::map_dfr(attributes, ~ as.data.frame(.x$mappings), .id = "attribute")

    ## Identifying nested columns and flattening them
    nested_cols <- df %>%
      dplyr::select(where(~ is.data.frame(.x) || (is.list(.x) && all(purrr::map_lgl(.x, is.data.frame)) ))) %>%
      names()

    ## Flattening the columns
    df <- df %>%
      unnest_wider(any_of(nested_cols), names_sep = ".")

    ## Correcting chain names mistaken as NA
    struct_cols <- grep("struct_asym_id", colnames(df), value = TRUE)

    # Replace NA with "NA" in all matching columns
    df <- df %>%
      dplyr::mutate(across(all_of(struct_cols), ~ tidyr::replace_na(., "NA")))
    df$query <- pdbid

    return(df)

  }
}
