#' Get NOISeq::readData ready annotations
#'
#' A convenience function to make a data frame containing biological annotations required by NOISeq to run NOISeq::readData function using custom GTF/GFF file and FASTA.
#'
#' @importFrom stringr str_split
#' @importFrom tools file_path_sans_ext
#' @import dplyr rtracklayer AnnotationForge
#' @export
#'
#' @param gff A link to GFF file from VEuPathDB or path to GTF file produced using AGAT.
#' @param fasta A link or path to the genome fasta file.
#' @param name Name of the organism.
#' @param select Type of features to be selected. By default we select "protein_coding_gene","ncRNA_gene"and "pseudogene" as they cover all the genes in the VEuPathDB annotation files. If using AGAT formatted GTF file, using "gene" would be sufficient as AGAT put all the genes types under gene tag.
#' @param geneidcol Use the tag that refers to the gene IDs in the GTF/GFF file.
#' @param genetype Use the tag that refers to the gene subtypes such as "ebi_biotype".
#'
#' @return A dataframe containing annotations per gene such as GC content, gene description, gene-length, gene start, gene end coordinates and chromosome information.
#' @examples
#' \dontrun{
#'  df <- easyNOISeqAnnot(
#'  gff="https://toxodb.org/common/downloads/release-68/EpraecoxHoughton/gff/data/ToxoDB-68_EpraecoxHoughton.gff",
#'  fasta = "https://toxodb.org/common/downloads/release-68/EpraecoxHoughton/fasta/data/ToxoDB-68_EpraecoxHoughton_Genome.fasta")
#'
#'
#' }
#'


easyNOISeqAnnot <- function(gff,fasta,name="Tgondii", select=c("protein_coding_gene","ncRNA_gene","pseudogene"), geneidcol="ID",genetype="ebi_biotype"){

  GTF <- rtracklayer::import.gff(gff, genome=name,feature.type=select)
  gene_ids <- S4Vectors::mcols(GTF)[[geneidcol]]
  gene_biotypes <- S4Vectors::mcols(GTF)[[genetype]]
  fasta_file <- Biostrings::readDNAStringSet(fasta,format = "fasta")

  names(fasta_file) <- names(fasta_file) %>% stringr::str_split(pattern = "[ ]",simplify = TRUE) %>% .[,1]

  seqs <- Biostrings::getSeq(fasta_file, GTF)
  gc_counts <- Biostrings::letterFrequency(seqs, letters = c("G", "C"), as.prob = FALSE)
  lengths <- BiocGenerics::width(GTF)
  gc_content <- rowSums(gc_counts) / lengths

  ## reading the FASTA reference
  gene_info <- data.frame(
    gene_id = gene_ids,
    gene_biotype = gene_biotypes,
    gc=gc_content,
    desc=GTF$description,
    length=lengths,
    starts=BiocGenerics::start(GTF),
    ends=BiocGenerics::end(GTF),
    chr=GenomeInfoDb::seqnames(GTF),
    stringsAsFactors = FALSE)

  rownames(gene_info) <- gene_info$gene_id

  ## Removing special characters if AGAT processed GTF file is used
  gene_info$desc<-gsub("[+]"," ", gene_info$desc)
  return(gene_info)
}
