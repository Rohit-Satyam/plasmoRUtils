#' Get tables with custom fields
#'
#' A convenience function to quickly fetch table of Gene IDs, Protein IDs, Gene Symbols, Annotations and many more columns from database of your choice such as PlasmoDB, ToxoDB, PiroplasmaDB among other VEuPathDB pathogen databases. 
#' 
#' With release 71 of VEuPathDB, only users with subscription and API key can access the database programmatically. Use `usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using `VEUPATHDB_API_KEY` variable.
#'
#' @import dplyr
#' @export
#'
#' @param org Full name of organism of interest as specified in VEuPathDB. To find the exact name of the organism, use `listVeupathdb` function.
#' @param db Character Name of the database in which the organism is present. These can be one of the following: "toxodb","plasmodb","hostdb","amoebadb","cryptodb","fungidb","giardiadb","microsporidiadb","piroplasmadb","trichdb","tritrypdb".
#' @param customFields A vector of custom fields desired to be fetched. "primary_key" is mandatory field. Use [getTableAttributes()] to browse the attributes available for the
#'   selected database and record type. If `NULL`, a default set of
#'   attributes is returned.
#'
#' @param api_key VEuPathDB API key.
#' 
#' @seealso [getTableAttributes()] for browsing the attributes available
#'   for a given record type.
#'   
#' @return A data frame, containing "Gene ID", "Product Description", "Gene Strand", "Gene Name or Symbol", "Previous ID(s)", "Entrez Gene ID", "UniProt ID(s)", "Protein Length", "TM Domains" and "SignalP Peptide" for all the genes present in the organism of interest.
#' @examples
#' \dontrun{
#' 
#' df <- getTable(org="Plasmodium falciparum 3D7", 
#' db="plasmodb",api_key = Sys.getenv("VEUPATHDB_API_KEY"))
#' 
#' df <- getTable(org="Plasmodium falciparum 3D7", 
#' db="plasmodb", api_key = Sys.getenv("VEUPATHDB_API_KEY"), 
#' customFields = c("primary_key","protein_sequence"))

#' }
#'

getTable <- function(
    org,
    db = "toxodb",
    customFields = NULL,
    api_key = Sys.getenv("VEUPATHDB_API_KEY")
) {
  
  if (!nzchar(api_key)) {
    stop("VEuPathDB API key is missing.", call. = FALSE)
  }
  
  part1 <- paste0(
    "service/record-types/transcript/searches/",
    "GenesByTaxon/reports/attributesTabular?organism=%5B%22"
  )
  
  query <- utils::URLencode(trimws(org))
  
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
  
  if (is.na(db_short)) {
    stop("Unsupported database: ", db, call. = FALSE)
  }
  
  
  if (is.null(customFields)) {
    
    part2 <- paste0(
      "%22%5D&reportConfig={",
      "%22attributes%22:[",
      "%22primary_key%22,",
      "%22gene_product%22,",
      "%22strand%22,",
      "%22gene_name%22,",
      "%22gene_previous_ids%22,",
      "%22gene_entrez_id%22,",
      "%22uniprot_ids%22,",
      "%22protein_length%22,",
      "%22tm_count%22,",
      "%22signalp_60_probability%22",
      "],",
      "%22includeHeader%22:true,",
      "%22attachmentType%22:%22plain%22}"
    )
    
    url <- paste0(
      "https://", db, ".org/",
      db_short, "/",
      part1, query, part2
    )
    
    out <- .read_authenticated_tsv(url,api_key=api_key)
    
    if ("Previous ID(s)" %in% names(out)) {
      out <- out %>%
        dplyr::mutate(
          `Previous ID(s)` = stringr::str_remove(
            `Previous ID(s)`,
            "Previous IDs: "
          )
        )
    }
    
    out
    
  } else {
    
    encodeit <- paste0(
      "%22",
      customFields,
      "%22",
      collapse = ","
    )
    
    part2 <- paste0(
      "%22%5D&reportConfig={%22attributes%22:[",
      encodeit,
      "],%22includeHeader%22:true,",
      "%22attachmentType%22:%22plain%22}"
    )
    
    url <- paste0(
      "https://", db, ".org/",
      db_short, "/",
      part1, query, part2
    )
    
    .read_authenticated_tsv(url,api_key=api_key)
  }
}
