#' Browse available VEuPathDB databases attributes
#'
#' Retrieves metadata describing the attributes available for a specified
#' VEuPathDB record type. The returned attribute names can be supplied to
#' the `customFields` argument of [getTable()]. 
#' 
#' With release 71 of VEuPathDB, only users with subscription and API key can access the database programmatically. Use `usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using `VEUPATHDB_API_KEY` variable.
#'
#' @import dplyr httr2
#' @param db Character. Name of the VEuPathDB database, such as
#'   `"plasmodb"`, `"toxodb"`, `"fungidb"` etc.
#' @param record_type Character. VEuPathDB record type for which attributes
#'   should be retrieved. Defaults to `"transcript"`.
#' @param api_key Character. VEuPathDB API key used for authentication.
#'   By default, the value is obtained from the `VEUPATHDB_API_KEY`
#'   environment variable.
#'
#' @return A data frame containing metadata for the available attributes,
#'   including the internal attribute name, display name, data type,
#'   reporting availability, and help text where available.
#'
#' @details
#' VEuPathDB record types can contain many reportable attributes.
#' `getAttributes()` provides a convenient way to browse these attributes
#' without manually consulting the VEuPathDB web service.
#'
#' The `name` column contains the attribute identifiers expected by the
#' `customFields` argument of [getTable()], while `displayName` provides
#' human-readable labels.
#'
#' @examples
#' \dontrun{
#' # Browse transcript attributes available in PlasmoDB
#' attrs <- getTableAttributes(
#'   db = "plasmodb",
#'   record_type = "transcript", api_key=Sys.getenv("VEUPATHDB_API_KEY")
#' )
#'
#' head(attrs)
#'
#' # Search for protein-related attributes
#' subset(
#'   attrs,
#'   grepl("protein", displayName, ignore.case = TRUE)
#' )
#'
#' # Use selected attribute names with getTable()
#' genes <- getTableAttributes(
#'   org = "Plasmodium falciparum 3D7",
#'   db = "plasmodb",
#'   customFields = c(
#'     "primary_key",
#'     "gene_product",
#'     "protein_length"
#'   )
#' )
#' }
#'
#' @seealso [getTable()] for retrieving data using selected attributes.
#' @family VEuPathDB data retrieval
#'
#' @export
#' 


getTableAttributes <- function(
    db = "plasmodb",
    record_type = "transcript",
    api_key = Sys.getenv("VEUPATHDB_API_KEY")
) {
  
  db_short <- c(
    toxodb = "toxo",
    plasmodb = "plasmo",
    hostdb = "hostdb",
    amoebadb = "amoeba",
    cryptodb = "cryptodb",
    fungidb = "fungidb",
    giardiadb = "giardiadb",
    microsporidiadb = "micro",
    piroplasmadb = "piro",
    trichdb = "trichdb",
    tritrypdb = "tritrypdb"
  )[tolower(db)]
  
  url <- paste0(
    "https://", db, ".org/",
    db_short,
    "/service/record-types/",
    record_type
  )
  
  resp <- httr2::request(url) |>
    httr2::req_auth_bearer_token(api_key) |>
    httr2::req_perform()
  
  metadata <- httr2::resp_body_json(
    resp,
    simplifyVector = TRUE
  )
  
  metadata$attributes |>
    dplyr::select(
      name,
      displayName,
      columnDataType,
      isDisplayable,
      isSortable,
      isInReport,
      help
    )
}