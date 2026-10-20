
#' Convert identifiers from one type to another
#'
#' Converts between current gene IDs, previous gene IDs,
#' UniProt IDs and Entrez IDs using getTableVeupathdb(). All possible combinations are allowed.
#'
#' Supports mixed-organism identifiers without requiring
#' organism or database information.
#'
#' All available gene-identifier associations are retained,
#' including one-to-many and many-to-many mappings.
#' 
#' With release 71 of VEuPathDB, only users with subscription and API key can access the database programmatically. Use `usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using `VEUPATHDB_API_KEY` variable.
#'
#' @import dplyr tidyr
#' @export
#'
#' @param inputid Character vector of identifiers to convert.Can be vector of mixed identifiers.
#' @param from Input identifier type. One of "ensembl", "uniprot",
#'   "entrez" or "old".
#' @param to Desired output identifier type. One of "ensembl",
#'   "uniprot", "entrez" or "old". Defaults to "ensembl".
#' @param api_key VEuPathDB API key.
#'
#' @return A data frame containing Input ID, Output ID, Gene ID
#'   and Organism. Multiple mappings are returned as separate rows.
#'
#' @examples
#' \dontrun{
#'
#' # Gene IDs to UniProt IDs
#' df <- toGeneid(
#'   c("PF3D7_1372300", "PBANKA_1355600"),
#'   from = "ensembl",
#'   to = "uniprot"
#' )
#'
#' # UniProt IDs to gene IDs
#' df <- toGeneid(
#'   c("A0A509APC6", "Q75UY1"),
#'   from = "uniprot",
#'   to = "ensembl"
#' )
#'
#' # Previous IDs to UniProt IDs
#' df <- toGeneid(
#'   "PBANKA_135560",
#'   from = "old",
#'   to = "uniprot"
#' )
#' }

toGeneid <- function(
    inputid,
    from = "",
    to = "ensembl",
    api_key = Sys.getenv("VEUPATHDB_API_KEY")
) {
  
  types <- c("ensembl", "uniprot", "entrez", "old")
  
  from <- match.arg(from, types)
  to <- match.arg(to, types)
  
  if (!is.character(inputid) || anyNA(inputid)) {
    stop(
      "inputid must be a character vector without NA values.",
      call. = FALSE
    )
  }
  
  inputid <- unique(trimws(inputid))
  inputid <- inputid[nzchar(inputid)]
  
  if (!length(inputid)) {
    stop("No valid input IDs supplied.", call. = FALSE)
  }
  
  # Retrieve matching records from integrated VEuPathDB
  fetchtable <- getTableVeupathdb(
    ids = inputid,
    api_key = api_key
  )
  
  # Define the corresponding output columns
  id_cols <- c(
    ensembl = "Gene ID",
    uniprot = "UniProt ID(s)",
    entrez  = "Entrez Gene ID",
    old     = "Previous ID(s)"
  )
  
  # Consistent empty-result structure
  empty <- tibble::tibble(
    `Input ID` = character(),
    `Output ID` = character(),
    `Gene ID` = character(),
    Organism = character()
  )
  
  if (nrow(fetchtable) == 0L) {
    message(
      "No matching VEuPathDB records found for the supplied IDs."
    )
    return(empty)
  }
  
  # Check required attributes
  required <- unique(c(
    "Gene ID",
    "Organism",
    unname(id_cols[c(from, to)])
  ))
  
  missing_cols <- setdiff(required, names(fetchtable))
  
  if (length(missing_cols)) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Convert an identifier column into one ID per row
  extract_ids <- function(column, type, value_name) {
    
    ref <- tibble::tibble(
      gene_id = as.character(fetchtable[["Gene ID"]]),
      organism = as.character(fetchtable[["Organism"]]),
      value = as.character(fetchtable[[column]])
    )
    
    # Remove prefix from historical gene IDs
    if (type == "old") {
      ref$value <- sub(
        "^Previous IDs:\\s*",
        "",
        ref$value
      )
    }
    
    # Split multiple identifiers
    ref <- tidyr::separate_rows(
      ref,
      value,
      sep = "[,;]"
    )
    
    # Remove whitespace and invalid identifiers
    ref$value <- trimws(ref$value)
    
    ref <- ref %>%
      dplyr::filter(
        !is.na(.data$value),
        !.data$value %in% c("", "N/A", "NA", "NULL"),
        !is.na(.data$gene_id),
        nzchar(.data$gene_id)
      ) %>%
      dplyr::distinct()
    
    # Rename the identifier column
    names(ref)[names(ref) == "value"] <- value_name
    
    return(ref)
  }
  
  # Extract source identifiers
  source <- extract_ids(
    column = unname(id_cols[from]),
    type = from,
    value_name = "input_id"
  )
  
  # Retain only the IDs requested by the user
  source <- source %>%
    dplyr::filter(.data$input_id %in% inputid)
  
  # Extract destination identifiers
  target <- extract_ids(
    column = unname(id_cols[to]),
    type = to,
    value_name = "output_id"
  )
  
  # 7. Map source IDs to destination IDs
  # Gene ID + Organism serves as the common mapping key.
  # This also combines annotations across transcript rows
  # belonging to the same gene.
  df <- source %>%
    dplyr::inner_join(
      target,
      by = c("gene_id", "organism"),
      relationship = "many-to-many"
    ) %>%
    dplyr::transmute(
      `Input ID` = .data$input_id,
      `Output ID` = .data$output_id,
      `Gene ID` = .data$gene_id,
      Organism = .data$organism
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(
      match(.data[["Input ID"]], inputid),
      .data[["Gene ID"]],
      .data[["Output ID"]]
    )
  
  # Report unresolved input IDs
  unresolved <- setdiff(
    inputid,
    unique(source$input_id)
  )
  
  if (length(unresolved)) {
    message(
      "Following IDs could not be resolved: ",
      paste(unresolved, collapse = ", ")
    )
  }
  
  # Report IDs having no annotation of the requested target type
  no_target <- setdiff(
    unique(source$input_id),
    unique(df[["Input ID"]])
  )
  
  if (length(no_target)) {
    message(
      "No ", to, " mapping available for: ",
      paste(no_target, collapse = ", ")
    )
  }
  
  return(df)
}
