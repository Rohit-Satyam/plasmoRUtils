#' List genomes and metadata in VEuPathDB
#'
#' A convenience function to quickly fetch table of genomes and their associated metadata from VEupathDB.
#'
#' @import readr jsonlite
#' @export
#'
#' @param customFields A vector of custom fields desired to be fetched. "primary_key" is mandatory field. Other fields can be supplied and can be chosen from (but are not limited to): `"annotation_source"`, `"annotation_version"`, `"arraygenecount"`, `"chipchipgenecount"`, `"chromosomeCount"`, `"codinggenecount"`, `"communitycount"`, `"contigCount"`, `"ecnumbercount"`, `"estcount"`, `"genecount"`, `"genecount_number"`, `"genome_source"`, `"genome_version"`, `"gocount"`, `"is_in_apollo"`, `"is_reference_strain"`, `"megabps"`, `"ncbi_tax_id"`, `"ncbi_taxon_url"`, `"organism"`, `"organism_full"`, `"orthologcount"`, `"othergenecount"`, `"popsetcount"`, `"project_id"`, `"proteomicscount"`, `"pseudogenecount"`, `"rnaseqcount"`, `"rtpcrcount"`, `"snpcount"`, `"species"`, `"species_ncbi_tax_id"`, `"species_ncbi_taxon_url"`, `"supercontigCount"`, `"tfbscount"`, `"URLcdsFasta"`, `"URLGenomeFasta"`, `"URLgff"`, `"URLproteinFasta"`, `"URLtranscriptFasta"`. For more fields, refer to the VEuPathDB Documentation
#' @param api_key Character. VEuPathDB API key used for authentication.
#'   By default, the value is obtained from the `VEUPATHDB_API_KEY`
#'   environment variable.
#' @param timeout Timeout time in seconds. Defaults to 60.
#' @param max_tries Maximum number of request attempts. Defaults to 3.
#' 
#' @return A data frame containing information about genomes present in VEuPathDB and their attributes.
#' @examples
#' \dontrun{
#' df <- listVeupathdb()
#' df <- listVeupathdb(customFields=c("species", "project_id"))
#' }
#'

listVeupathdb <- function(
    customFields = NULL,
    api_key = Sys.getenv("VEUPATHDB_API_KEY"),
    timeout = 60,
    max_tries = 3
) {
  
  if (!nzchar(api_key)) {
    stop(
      "VEuPathDB API key is missing. ",
      "Set the VEUPATHDB_API_KEY environment variable.",
      call. = FALSE
    )
  }
  
  defaultFields <- c(
    "primary_key",
    "species",
    "URLGenomeFasta",
    "URLcdsFasta",
    "URLtranscriptFasta",
    "URLproteinFasta",
    "project_id",
    "genecount",
    "URLgff",
    "genome_source",
    "annotation_source"
  )
  
  if (is.null(customFields)) {
    fields <- defaultFields
  } else {
    # primary_key is required
    fields <- unique(c("primary_key", customFields))
  }
  
  report_config <- jsonlite::toJSON(
    list(
      attributes = fields,
      includeHeader = TRUE,
      attachmentType = "plain"
    ),
    auto_unbox = TRUE
  )
  
  base_url <- paste0(
    "https://veupathdb.org/veupathdb/service/",
    "record-types/organism/searches/",
    "GenomeDataTypes/reports/attributesTabular"
  )
  
  resp <- httr2::request(base_url) |>
    httr2::req_auth_bearer_token(api_key) |>
    httr2::req_url_query(
      reportConfig = as.character(report_config)
    ) |>
    httr2::req_user_agent("plasmoRUtils listVeupathdb") |>
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
      "VEuPathDB request failed with HTTP ",
      httr2::resp_status(resp),
      "\n\nServer response:\n",
      txt,
      call. = FALSE
    )
  }
  
  readr::read_tsv(
    I(txt),
    progress = FALSE,
    show_col_types = FALSE
  )
}