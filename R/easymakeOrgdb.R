#' Create Org.db package quickly
#'
#' A convenience function to make OrgDb packages quickly to be used with other GO enrichment packages such as ClusterProfiler.
#'
#' @importFrom stringr str_split
#' @importFrom tools file_path_sans_ext
#' @importFrom readr read_tsv
#' @import dplyr rtracklayer AnnotationForge
#' @export
#'
#' @param gff A link to GFF file from VEuPathDB or path to GFF file if the file is present locally.
#' @param gaf Gene Ontology file obtained from VEuPathDB.
#' @param out.dir Output directory where the package will be saved.
#' @param taxid Taxonomy ID of the organism. This can be obtained from https://www.ncbi.nlm.nih.gov/taxonomy.
#' @param genus Genus of the organism. This will be used to construct the name of the package.
#' @param sp Species with or without strain information.
#' @param version Version of the package if you choose to maintain and share the package.
#' @param verbose Display the messages when running makeOrgPackage function.
#' @param maintainer Email Id of the package builder. Default "John doe <johndoe@gmail.com>"
#'
#' @return A tar.gz file that can be installed as a package and can be used with GO enrichment tools such as ClusterProfiler.
#' @examples
#' \dontrun{
#'  easymakeOrgdb(
#'   gff="https://plasmodb.org/common/downloads/release-68/PbergheiANKA/gff/data/PlasmoDB-68_PbergheiANKA.gff",
#'   gaf="https://plasmodb.org/common/downloads/release-68/PbergheiANKA/gaf/PlasmoDB-68_PbergheiANKA_Curated_GO.gaf.gz",
#'   out.dir=".", taxid=5823,genus="Plasmodium",
#'   sp="bergheiANKA",
#'   version=0.1,
#'   verbose = FALSE,
#'   maintainer="John doe <johndoe@gmail.com>")
#'
#'
#' }
#'


easymakeOrgdb <- function(gff="https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff",gaf="https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gaf/PlasmoDB-68_Pfalciparum3D7_Curated_GO.gaf.gz",out.dir=".", taxid=36329,genus="Plasmodium",sp="falciparum3D7",version=0.1, verbose = FALSE, maintainer="John doe <johndoe@gmail.com>"){


  gff_data <- rtracklayer::import.gff(gff, format="gff", genome=tools::file_path_sans_ext(basename(gff)) %>%
  stringr::str_split(., pattern = "_",n = 3,simplify = TRUE) %>%
                                        .[,2]) %>%
    .[.$type == "protein_coding_gene", ]

  gene_2_GO <- readr::read_tsv(gaf, comment = "!", col_names = FALSE,trim_ws = TRUE, progress = FALSE,show_col_types = FALSE)

  colnames(gene_2_GO) <- paste0("col",seq_len(ncol(gene_2_GO)))
  ## Subset the gff_data
  #gff_data <- gff_data[gff_data$ID %in% unique(gene_2_GO$X2),]

  go_mappings <- data.frame(
    GID = gene_2_GO$col2,    # Gene ID
    GO = gene_2_GO$col5,          # GO Term
    EVIDENCE = gene_2_GO$col7,
    stringsAsFactors = FALSE
  )

  gene_info <- data.frame(
    GID = gff_data$ID,            # Gene ID from GFF
    SYMBOL = gff_data$Name,       # Gene symbol
    CHR = as.character(GenomeInfoDb::seqnames(gff_data)),           # Chromosome
    START = BiocGenerics::start(gff_data),            # Start position
    END = BiocGenerics::end(gff_data),                # End position
    STRAND = as.character(BiocGenerics::strand(gff_data)),
    stringsAsFactors = FALSE
  )

  ## Since some gene Symbols are absent so substitute with the GIds
  gene_info[is.na(gene_info[,2]),]$SYMBOL <- gene_info[is.na(gene_info[,2]),]$GID

  library(AnnotationForge)
  orgpath <- paste0(out.dir,"/","org.",substr(genus, 1, 1),sp,".eg.db")
  if (file.exists(orgpath)) {
    message(orgpath, " already exists, deleting it.")
    ret <- unlink(orgpath, recursive = TRUE)
  }

  AnnotationForge::makeOrgPackage(chromosome=unique(gene_info),
                 go=unique(go_mappings),
                 version=as.character(version),
                 maintainer= maintainer,
                 author=maintainer,
                 outputDir = out.dir,
                 tax_id=as.character(taxid),
                 genus=genus,
                 species=sp,
                 goTable="go",verbose = verbose)
  system(paste0("tar -czvf ", paste0(basename(orgpath),".",as.character(version),".tar.gz "),orgpath))
  system(paste("rm -rf ",orgpath))
}
