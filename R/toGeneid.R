#' Convert Other IDs to Ensembl gene IDs
#'
#' A convenience function to quickly convert the Uniprot or Entrez Ids to Ensembl gene IDs,  using VEuPathDB specialized databases.It also provides description and gene symbol for input Ids.
#'
#' @import dplyr tidyr
#' @importFrom glue glue
#' @export
#'
#' @param inputid A character vector of IDs. Can be Ensembl, Uniprot, Entrez or or old Pf ids.
#' @param from To describe the type of Input ID. Possible values: "old", "uniprot". "entrez", "ensembl"
#' @param to To describle the type of output ID desired. Possible values: "emsembl".
#' @param org Organism to which the IDs belongs to.Possible values: toxodb, plasmodb, hostdb, amoebadb, cryptodb, fungidb, giardiadb, microsporidiadb, piroplasmadb, trichdb, tritrypdb.
#' @param db Database in which organism is present.
#' @param ... Additional arguments that can be passed to the \code{getTable} function.
#'
#' @return A data frame, containing Gene IDs, gene description and gene Symbols and more.
#' @examples
#' \dontrun{
#' df <- toGeneid(
#' c("PF3D7_0420300", "PF3D7_0621000"),
#'       from="ensembl")
#' }
#'

toGeneid <- function (inputid, from = "", to = "", org="Plasmodium falciparum 3D7", db="plasmodb",...)
{
  fetchtable <- getTable(org=org,db=db,...)
  get_ref <- function(col_name) {
    fetchtable %>%
      dplyr::select(`Gene ID`, `Product Description`, `Gene Name or Symbol`, all_of(col_name)) %>%
      tidyr::separate_rows(all_of(col_name),sep = "[,;]") %>%
      dplyr::filter(all_of(col_name) !="N/A") %>%
      dplyr::distinct()
  }
  check_failed <- function(ref, col_name) {
    failed <- unique(inputid[!inputid %in% ref[[col_name]]])
    if (length(failed) > 0) {
      message(glue::glue("Following genes failed to convert: {paste(failed, collapse = ' ')}\n"))
    }
    ref %>% dplyr::filter(.data[[col_name]] %in% inputid)
  }

  if (from == "old" && to == "ensembl") {
    ref <- get_ref("Previous ID(s)")
    df <- check_failed(ref, "Previous ID(s)")
  } else if (from == "uniprot" && to == "ensembl") {
    ref <- get_ref("UniProt ID(s)")
    df <- check_failed(ref, "UniProt ID(s)")
  } else if (from == "entrez" && to == "ensembl") {
    ref <- get_ref("Entrez Gene ID")
    df <- check_failed(ref, "Entrez Gene ID")
  } else if (from == "ensembl") {
    ref <- fetchtable %>% dplyr::distinct()
    df <- check_failed(ref, "Gene ID")
  } else {
    message("Provide all necessary arguments. For help use ?toGeneid")
    return(NULL)
  }
  return(df)
}

