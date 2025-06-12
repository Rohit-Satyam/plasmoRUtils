#' List genomes and metadata in VEuPathDB
#'
#' A convenience function to quickly fetch table of genomes and their associated metadata from VEupathDB.
#'
#' @import readr
#' @importFrom glue glue
#' @export
#'
#' @param customFields A vector of custom fields desired to be fetched. "primary_key" is mandatory field. Other fields can be supplied and can be chosen from (but are not limited to): "annotation_source", "annotation_version", "arraygenecount", "chipchipgenecount", "chromosomeCount", "codinggenecount", "communitycount", "contigCount", "ecnumbercount", "estcount", "genecount", "genecount_number", "genome_source", "genome_version", "gocount", "is_in_apollo", "is_reference_strain", "megabps", "ncbi_tax_id", "ncbi_taxon_url", "organism", "organism_full", "orthologcount", "othergenecount", "popsetcount", "project_id", "proteomicscount", "pseudogenecount", "rnaseqcount", "rtpcrcount", "snpcount", "species", "species_ncbi_tax_id", "species_ncbi_taxon_url", "supercontigCount", "tfbscount", "URLcdsFasta", "URLGenomeFasta", "URLgff", "URLproteinFasta", "URLtranscriptFasta". For more fields, refer to the VEuPathDB Documentation
#'
#' @return A data frame containing information about genomes present in VEuPathDB and their attributes.
#' @examples
#' \dontrun{
#' df <- listVeupathdb()
#' df <- listVeupathdb(customFields=c("species", "project_id"))
#' }
#'


listVeupathdb <- function(customFields=NULL){

  if(is.null(customFields)){
    df <- readr::read_tsv(
      utils::URLencode(
        'https://veupathdb.org/veupathdb/service/record-types/organism/searches/GenomeDataTypes/reports/attributesTabular?reportConfig={"attributes": ["primary_key","species","URLGenomeFasta","URLcdsFasta","URLtranscriptFasta","URLproteinFasta","project_id","genecount","URLgff","project_id","URLgff","project_id","genome_source","annotation_source"],"includeHeader":true,"attachmentType":"plain"}'
        ),progress = FALSE,show_col_types = FALSE
      )
    return(df)


  } else{
    # Collapse into comma-separated string of quoted field names
    field_str <- glue_collapse(glue('"{customFields}"'), sep = ",")

    # Build full URL with correct JSON query syntax
    url <- glue(
      'https://veupathdb.org/veupathdb/service/record-types/organism/searches/GenomeDataTypes/reports/attributesTabular?reportConfig={{"attributes": [{field_str}],"includeHeader":true,"attachmentType":"plain"}}'
    )

    df <- readr::read_tsv(utils::URLencode(url),progress = FALSE,show_col_types = FALSE)
    return(df)
  }



}

