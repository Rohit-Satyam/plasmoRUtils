#' Get tables with custom fields
#'
#' A convenience function to quickly fetch table of Gene IDs, Protein IDs, Gene Symbols, Annotations and many more columns from database of your choice such as PlasmoDB, ToxoDB, PiroplasmaDB among other VEuPathDB pathogen databases.
#'
#' @import dplyr tidyr
#' @importFrom glue glue
#' @export
#'
#' @param org Full name of organism of interest as specified in VEuPathDB. To find the exact name of the organism, use `listVeupathdb` function.
#' @param db Character Name of the database in which the organism is present. These can be one of the following: "toxodb","plasmodb","hostdb","amoebadb","cryptodb","fungidb","giardiadb","microsporidiadb","piroplasmadb","trichdb","tritrypdb".
#' @param customFields A vector of custom fields desired to be fetched. "primary_key" is mandatory field. Other fields can be supplied and can be chosen from (but are not limited to): "organism",   "gene_location_text",   "gene_product",   "gene_type",  "exon_count",   "gene_exon_count",   "gene_transcript_count",   "three_prime_utr_length",   "five_prime_utr_length",   "strand",   "is_pseudo",   "transcript_length",   "is_deprecated",   "gene_name",   "gene_source_id",   "transcript_product",   "protein_length",   "chromosome",   "location_text",   "sequence_id",   "gene_ortholog_number",   "gene_orthomcl_name",   "gene_paralog_number",   "cds_length",   "molecular_weight",   "isoelectric_point",   "tm_count",   "signalp_peptide",   "predicted_go_id_component",   "predicted_go_component",   "predicted_go_id_function",   "predicted_go_function",   "predicted_go_id_process",   "predicted_go_process",   "annotated_go_id_component",   "annotated_go_component",   "annotated_go_id_function",   "annotated_go_function",   "annotated_go_id_process",   "annotated_go_process",   "ec_numbers",   "ec_numbers_derived"
#'
#' @return A data frame, containing "Gene ID", "Product Description", "Gene Strand", "Gene Name or Symbol", "Previous ID(s)", "Entrez Gene ID", "UniProt ID(s)", "Protein Length", "TM Domains" and "SignalP Peptide" for all the genes present in the organism of interest.
#' @examples
#' \dontrun{
#' df <- getTable(org="Plasmodium falciparum 3D7", db="plasmodb")
#' }
#'

getTable <- function(org,db="toxodb",customFields=NULL){
  # URL components
  part1 <- "service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular?organism=%5B%22"

  query <-  utils::URLencode(trimws(org))
  db_short <- c(toxodb = "toxo", plasmodb = "plasmo", hostdb = "hostdb", amoebadb = "amoeba", cryptodb = "cryptodb",
                fungidb = "fungidb", giardiadb = "giardiadb", microsporidiadb = "micro", piroplasmadb = "piro",
                trichdb = "trichdb", tritrypdb = "tritrypdb")[tolower(db)]

  if(is.null(customFields)){
    part2 <- "%22%5D&reportConfig={%22attributes%22:[%22primary_key%22,%22gene_product%22,%22strand%22,%22gene_name%22,%22gene_previous_ids%22,%22gene_entrez_id%22,%22uniprot_ids%22,%22protein_length%22,%22tm_count%22,%22signalp_peptide%22],%22includeHeader%22:true,%22attachmentType%22:%22plain%22}"

    url <- paste0("https://", db, ".org/", db_short, "/", part1, query, part2)
    readr::read_tsv(url,show_col_types = FALSE,progress = FALSE) %>%
      dplyr::mutate(`Previous ID(s)` = stringr::str_remove(`Previous ID(s)`, "Previous IDs: "))

  } else {
    encodeit <- paste0('%22', customFields, '%22', collapse = ',')
    part2 <- paste0("%22%5D&reportConfig={%22attributes%22:[",encodeit,"],%22includeHeader%22:true,%22attachmentType%22:%22plain%22}")

    url <- paste0("https://", db, ".org/", db_short, "/", part1, query, part2)
    readr::read_tsv(url,progress = FALSE,show_col_types = FALSE)
  }


}

