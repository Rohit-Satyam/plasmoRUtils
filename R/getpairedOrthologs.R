#' Fetch paired orthologs
#'
#' This function provides ability to query InParanoiDB 9 and OrthoMCL 7 to get paired orthologs between two species of interest.
#'
#' @import dplyr
#' @importFrom glue glue
#' @import readr
#' @importFrom jsonlite toJSON
#' @importFrom tidyr separate_rows
#' @export
#'
#' @param from ID of query organism for which orthologs are to be fetched. To view organisms indexed in the database and their respective IDs use \code{listipdb()} or \code{listOrthomcl()}.
#' @param to ID of target organism against which will be queried for orthologs. To view organisms indexed in the database and their respective IDs use \code{listipdb()} or \code{listOrthomcl()}.
#' @param db Define database to be queried. Possible values: "orthomcl","ipdb".
#' @param customFields Additional field to be fetched from OrthoMCL 7. "primary_key" and "target_id" are mandatory fields. Additional popular fields include: "group_name","product","source_id","num_core","num_peripheral","length","sequence","taxon_name","abbreviation","core_peripheral","ec_numbers","pfam_domains" etc.For more fields refer to OrthoMCL REST query builder.
#' @param transform Logical. In case of InParanoiDB, transform collapses the orthologs by Group IDs so that you have unique rows. For OrthoMCL, it performs the opposite function and separates the rows so that you have only one ID in each column. When set to FALSE, the function returns the ortholog query results in the raw form provided by the database.
#'
#'
#' @return df This function returns a dataframe of orthologs and additional requested fields.
#' @examples
#' \dontrun{
#' df <- getpairedOrthologs(from=1742855, to=1747281, db="orthomcl", transform = FALSE)
#' df <- getpairedOrthologs(from=1742855, to=1747281, db="orthomcl", 
#' transform = FALSE, customFields = c("primary_key","target_id","group_name"))
#' df <- getpairedOrthologs(from=1742855, to=1747281, db="orthomcl", transform = TRUE)
#' df <- getpairedOrthologs(from=423536,to=36329, transform = FALSE)
#' df <- getpairedOrthologs(from=423536,to=36329, transform = TRUE)
#' }
#'

getpairedOrthologs <- function(from, to, db = c("ipdb", "orthomcl"),
                               customFields = NULL, transform = TRUE) {
  db <- match.arg(db)
  
  # helper: collapse rows within a group as comma-separated unique values
  collapse_by_group <- function(x) {
    dplyr::group_by(x, `Group-id`) %>%
      dplyr::summarise(
        dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ",")),
        .groups = "drop"
      )
  }
  
  if (db == "ipdb") {
    url <- glue::glue("https://inparanoidb.sbc.su.se/download/sqltable/{from}&{to}&prot")
    df  <- readr::read_tsv(url, col_names = FALSE, show_col_types = FALSE, progress = FALSE)
    colnames(df) <- c("Group-id", "Bitscore", "Species", "Inparalog-score",
                      "Protein-name", "Seed-score")
    df$Species <- sub("\\.fa$", "", df$Species)
    
    if (!transform) return(df)
    
    # taxonomy name map from user-supplied helper
    tax <- listipdb()  # expects columns: #NCBITaxID, Species_name
    
    # split -> label species nicely -> collapse -> prefix -> join
    query  <- df %>%
      dplyr::filter(Species == from) %>%
      dplyr::mutate(Species = tax$Species_name[match(from, tax$`#NCBITaxID`)]) %>%
      collapse_by_group() %>%
      dplyr::rename_with(~ paste0("query_", .x))
    
    target <- df %>%
      dplyr::filter(Species == to) %>%
      dplyr::mutate(Species = tax$Species_name[match(to, tax$`#NCBITaxID`)]) %>%
      collapse_by_group() %>%
      dplyr::rename_with(~ paste0("target_", .x))
    
    return(dplyr::full_join(query, target,
                            by = c("query_Group-id" = "target_Group-id")))
  }
  
  # ---- orthomcl branch ----
  attrs  <- if (is.null(customFields)) c("primary_key", "target_id") else customFields
  report <- jsonlite::toJSON(
    list(attributes = attrs, includeHeader = TRUE, attachmentType = "plain"),
    auto_unbox = TRUE
  )
  u <- paste0(
    "https://orthomcl.org/orthomcl/service/record-types/sequence/",
    "searches/BySharedOrtholog/reports/attributesTabular",
    "?query_organism_type_ahead=", from,
    "&target_organism_type_ahead=", to,
    "&reportConfig=", utils::URLencode(report, reserved = TRUE)
  )
  
  df <- readr::read_tsv(u, show_col_types = FALSE, progress = FALSE)
  if (transform && "Target ID" %in% names(df)) {
    df <- tidyr::separate_rows(df, `Target ID`, sep = ",")
  }
  df
}
