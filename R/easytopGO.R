#' Performing quick ORA analysis
#'
#' A convenience function to quickly perform GO Term Enrichment analysis using TopGO.The results then can be plotted using \code{easyGOPlot}
#'
#' @import dplyr tibble
#' @import biomaRt topGO
#' @importFrom readr read_tsv
#' @export
#'
#' @param geneID A vector of named p-values. The names should be gene IDs.
#' @param bkggset A character vector of gene IDs that should be used as background.
#' @param gaf A URL or path to the .gaf file obtained from PlasmoDB should you choose not to use biomaRt. When using this argument, set useBiomart=FALSE and useGAF=TRUE.
#' @param useBiomart Logical To enable usage of BiomaRt to fetch GO terms. Default: TRUE.
#' @param useGAF Logical To enable usage of custom .gaf file to fetch GO terms. Default: FALSE
#' @param mart Argument to specify the mart for BiomaRt functions. Default: "protists_mart".
#' @param gset Argument to specify the geneset to be used by BiomaRt functions. Default: "pfalciparum_eg_gene"
#' @param algo Argument to specify which algorithm to be used for enrichment by topGO. For possible options use \code{topGO::whichAlgorithms()}
#' @param stats Argument to specify the statistical test to be used for enrichment by topGO. For possible values use \code{topGO::whichTests()}. Default: "ks".
#' @param category Specify the category of Over-representation analysis such as "BP" for Biological Process, "MF" for Molecular Function and "CC" for Cellular Component Enrichment.
#' @param fdr logical. Perform multiple testing correction testing. Default (FALSE)
#' @param correction Method to be used to calculate adjusted p-value. Possible values: ""holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr". Read section 6.2 in topGO documentation before performing correction. Correction when using elim and weight is usually not recommended.
#'
#' @return A dataframe of enriched terms with GO description and genes filteres by uncorrected p-values.
#' @examples
#' \dontrun{
#' ## making gene list from DESEq2
#' geneList <- subset(res, regulate=="Up") %>% .$padj
#' names(geneList) <- subset(res, regulate=="Up") %>% .$Geneid
#'
#' ## background genes will be the genes tested for differential expression
#' background.gset <- res$Geneid
#' baseurl <- "https://plasmodb.org/common/downloads/Current_Release/"
#' url<-paste0(baseurl,"Pfalciparum3D7/gaf/PlasmoDB-68_Pfalciparum3D7_GO.gaf.gz")
#' gores<-easytopGO(geneID = geneList,useGAF = TRUE,useBiomart = FALSE,gaf=url,
#' bkggset = background.gset, category = "BP", stats = "ks")
#' }
#'

easytopGO <- function(geneID, bkggset="",gaf="",useBiomart=TRUE,useGAF=FALSE, mart = "protists_mart", gset = "pfalciparum_eg_gene",algo="weight01",stats="ks",category="BP", fdr=FALSE, correction="BY"){
  require(topGO)
  if(useBiomart){
    db <- biomaRt::useMart(biomart = mart, host = "https://protists.ensembl.org", dataset = gset)
    go_ids <- biomaRt::getBM(attributes = c('go_id', 'ensembl_gene_id', 'namespace_1003'), filters = 'ensembl_gene_id', values = bkggset, mart = db) %>%
      dplyr::filter(go_id != "")
    gene_2_GO <- utils::unstack(go_ids[, c(1, 2)])

    goEnrichment <- .usetopGO(stats=stats,category=category,geneID=geneID,gene_2_GO=gene_2_GO, algo=algo,fdr=fdr,correction=correction)
    return(goEnrichment)
  }

  if(all(useGAF && gaf!="")){
    gene_2_GO <- readr::read_tsv(gaf, comment = "!", col_names = FALSE, trim_ws = TRUE, progress = FALSE,show_col_types = FALSE) %>%
      dplyr::select(X2, X5) %>%
      dplyr::distinct() %>%
      dplyr::group_by(X2) %>%
      dplyr::summarise(PFID_list = list(X5)) %>%
      tibble::deframe() %>%
      base::Filter(function(x) length(x) > 0, .)

    goEnrichment <- .usetopGO(stats=stats,category=category,geneID=geneID,gene_2_GO=gene_2_GO,algo=algo,fdr=fdr,correction=correction)
    return(goEnrichment)
  }
}

utils::globalVariables(c("go_id","ensembl_gene_id","namespace_1003","X2","X5"))
