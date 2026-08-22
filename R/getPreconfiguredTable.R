#' Get pre-configured tables
#'
#' A convenience function to quickly fetch preconfigured table of Signal Peptide ranges, Pathways, Pubmed entries related to genes, Annotations and etc from database of your choice such as PlasmoDB, ToxoDB, PiroplasmaDB among other VEuPathDB pathogen databases.
#'
#' @import dplyr tidyr
#' @export
#'
#' @param org Full name of organism of interest as specified in VEuPathDB. To find the exact name of the organism, use `listVeupathdb` function.
#' @param db Character Name of the database in which the organism is present. These can be one of the following: "toxodb","plasmodb","hostdb","amoebadb","cryptodb","fungidb","giardiadb","microsporidiadb","piroplasmadb","trichdb","tritrypdb".
#' @param customField Preconfigured table that you wish the fetch. Pass only one value at a time from the following: "GeneModelDump",   "GeneTranscripts",   "Alias",   "GeneLinkouts",   "GeneLocation",   "PubMed",   "OrthologsLite",   "LowComplexity",   "PdbSimilarities",   "3dPreds",   "AlphaFoldLinkouts",   "ProteinProperties",   "InterPro",   "SignalP",   "TMHMM",   "ECNumbers",   "ECNumbersInferred",   "protein_length",   "chromosome",   "location_text",   "sequence_id",   "gene_ortholog_number",   "gene_orthomcl_name",   "gene_paralog_number",   "MetabolicPathwaysMPMP",   "MetabolicPathways",   "CompoundsMetabolicPathways",   "Y2hInteractions",   "MassSpecDownload",   "MassSpecMod",   "Epitopes" etc. You can view list of 58 preconfigured tables across VEuPathDB databases using `getPreconfiguredTable(listtables = TRUE)`
#' @param listtables To get list of preconfigured tables that can be supplied as an argument to `customField`.
#' @param api_key Character. VEuPathDB API key used for authentication.
#'   By default, the value is obtained from the `VEUPATHDB_API_KEY`
#'   environment variable.

#' @return A data frame.
#' 
#' @details
#' With release 71 of VEuPathDB, only users with subscription and API key can access the database programmatically. Use `usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using `VEUPATHDB_API_KEY` variable.
#' 
#' @examples
#' \dontrun{
#' getPreconfiguredTable(listtables = TRUE)
#' df <- getPreconfiguredTable(org = "Plasmodium falciparum 3D7",
#'      db = "plasmodb",customField = "Y2hInteractions")
#' }
#'

getPreconfiguredTable <- function(
    org = NULL,
    db = "plasmodb",
    customField = "Y2hInteractions",
    listtables = FALSE,
    api_key = Sys.getenv("VEUPATHDB_API_KEY")
) {
  
  # Return available preconfigured tables
  if (listtables) {
    
    url <- paste0(
      "https://raw.githubusercontent.com/",
      "VEuPathDB/ApiCommonModel/master/",
      "Model/lib/wdk/ontology/individuals.txt"
    )
    
    x <- readLines(url, warn = FALSE)
    
    # Split each ontology line into tab-separated fields
    fields <- strsplit(x, "\t", fixed = TRUE)
    
    # Keep actual GeneRecord table definitions that are downloadable
    keep <- vapply(
      fields,
      function(z) {
        length(z) >= 6 &&
          z[4] == "GeneRecordClasses.GeneRecordClass" &&
          z[5] == "table" &&
          "download" %in% z
      },
      logical(1)
    )
    
    fields <- fields[keep]
    
    # Field 6 contains the actual table name
    table_name <- vapply(
      fields,
      function(z) z[6],
      character(1)
    )
    
    return(sort(unique(table_name)))
  }
  
  # org is required when fetching a table
  if (is.null(org)) {
    stop(
      "`org` must be provided when listtables = FALSE.",
      call. = FALSE
    )
  }
  
  # URL components
  part1 <- paste0(
    "service/record-types/transcript/searches/",
    "GenesByTaxon/reports/tableTabular?organism=%5B%22"
  )
  
  query <- utils::URLencode(org)
  
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
  
  part2 <- paste0(
    "%22%5D&reportConfig=%7B%22tables%22%3A%5B%22",
    customField,
    "%22%5D%2C%22includeHeader%22%3Atrue",
    "%2C%22attachmentType%22%3A%22plain%22%7D"
  )
  
  url <- paste0(
    "https://", db, ".org/",
    db_short, "/",
    part1, query, part2
  )
  
  out <- .read_authenticated_tsv(
    url,
    api_key = api_key
  )
  
  return(out)
}

