#' Get pre-configured tables
#'
#' A convenience function to quickly fetch preconfigured table of Signal Peptide ranges, Pathways, Pubmed entries related to genes, Annotations and etc from database of your choice such as PlasmoDB, ToxoDB, PiroplasmaDB among other VEuPathDB pathogen databases.
#'
#' @import dplyr tidyr
#' @importFrom glue glue
#' @export
#'
#' @param org Full name of organism of interest as specified in VEuPathDB. To find the exact name of the organism, use `listVeupathdb` function.
#' @param db Character Name of the database in which the organism is present. These can be one of the following: "toxodb","plasmodb","hostdb","amoebadb","cryptodb","fungidb","giardiadb","microsporidiadb","piroplasmadb","trichdb","tritrypdb".
#' @param customField Preconfigured table that you wish the fetch. Pass only one value at a time from the following: "GeneModelDump",   "GeneTranscripts",   "Alias",   "GeneLinkouts",   "GeneLocation",   "PubMed",   "OrthologsLite",   "LowComplexity",   "PdbSimilarities",   "3dPreds",   "AlphaFoldLinkouts",   "ProteinProperties",   "InterPro",   "SignalP",   "TMHMM",   "ECNumbers",   "ECNumbersInferred",   "protein_length",   "chromosome",   "location_text",   "sequence_id",   "gene_ortholog_number",   "gene_orthomcl_name",   "gene_paralog_number",   "MetabolicPathwaysMPMP",   "MetabolicPathways",   "CompoundsMetabolicPathways",   "Y2hInteractions",   "MassSpecDownload",   "MassSpecMod",   "Epitopes" etc.
#'
#' @return A data frame.
#' @examples
#' \dontrun{
#' df <- getPreconfiguredTable(org = "Plasmodium falciparum 3D7",
#'      db = "plasmodb",customField = "Y2hInteractions")
#' }
#'

getPreconfiguredTable <- function(org,db="plasmodb",customField="Y2hInteractions"){
  # URL components
  part1 <- "service/record-types/transcript/searches/GenesByTaxon/reports/tableTabular?organism=%5B%22"

  query <-  utils::URLencode(org)
  db_short <- c(toxodb = "toxo", plasmodb = "plasmo", hostdb = "hostdb", amoebadb = "amoeba", cryptodb = "cryptodb",
                fungidb = "fungidb", giardiadb = "giardiadb", microsporidiadb = "micro", piroplasmadb = "piro",
                trichdb = "trichdb", tritrypdb = "tritrypdb")[tolower(db)]

  part2 <- paste0("%22%5D&reportConfig=%7B%22tables%22%3A%5B%22",customField,"%22%5D%2C%22includeHeader%22%3Atrue%2C%22attachmentType%22%3A%22plain%22%7D")

  url <- paste0("https://", db, ".org/", db_short, "/", part1, query, part2)
  readr::read_tsv(url,progress = F,show_col_types = FALSE)
}
