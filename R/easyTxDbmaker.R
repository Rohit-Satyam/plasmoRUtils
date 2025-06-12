
#' Quickly make TxDb object
#'
#' A convenience function to make OrgDb packages quickly to be used with other GO enrichment packages such as ClusterProfiler.
#'
#' @import GenomeInfoDb
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_split
#' @import txdbmaker rtracklayer
#' @export
#'
#' @param gff A link to GFF file from VEuPathDB or path to GFF file if the file is present locally.
#' @param fasta Genome FASTA file obtained from VEuPathDB.
#' @param abbr Abbreviation of the organism.
#' @param org Name of the organism including genus and species.
#' @param taxid NCBI Taxonomy ID. This can be obtained from https://www.ncbi.nlm.nih.gov/taxonomy.
#' @param db Name of the database with release information.
#'
#' @return A tar.gz file that can be installed as a package and can be used with GO enrichment tools such as ClusterProfiler.
#' @examples
#' \dontrun{
#'  txdb<-easyTxDbmaker(
#'  gff="https://toxodb.org/common/downloads/release-68/TgondiiME49/gff/data/ToxoDB-68_TgondiiME49.gff",
#'  fasta="https://toxodb.org/common/downloads/release-68/TgondiiME49/fasta/data/ToxoDB-68_TgondiiME49_Genome.fasta",
#'  abbr="TgondiiME49",
#'  taxid=508771,org = "Toxoplasma gondii ME49",
#'  db = "ToxoDB release 68")
#'
#'
#' }
#'


easyTxDbmaker <- function(gff, fasta, abbr="TgondiiME49",org = "Toxoplasma gondii",taxid=508771, db = "ToxoDB release 68"){

  metadata <- data.frame(name=c("GFF/GTF","FASTA"), value=c(gff,fasta))
  fasta_file <- Biostrings::readDNAStringSet(fasta,format = "fasta")

  # Calculate lengths and return seqinfo object
  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = names(fasta_file) %>% stringr::str_split(pattern = "[ ]",simplify = TRUE) %>% .[,1],
    seqlengths = width(fasta_file)
  )


  txdb <- txdbmaker::makeTxDbFromGFF(gff,chrominfo = seqinfo,metadata = metadata, organism = org,taxonomyId = taxid,dataSource = db)

  genome(txdb) <- abbr

  return(txdb)


}


