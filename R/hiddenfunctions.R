#' plasmoRUtils
#'
#' Hidden function to fetch the expandable tables from MPMP database
#'
#' @import glue rvest
#' @param url MPMP url.
#' @keywords internal

.clickable <- function(url) {
  df <- rvest::read_html(url) %>%
    rvest::html_elements("table.table-bordered.table-hover") %>%
    rvest::html_table()

  if (length(df) == 0) {
    message(glue::glue("No table found for the given URL: {url}\n"))
  } else {
    message(glue::glue("\033[0;32mfetched successfully: {url}\033[0m\n"))
    return(df[[1]])
  }
}

#' plasmoRUtils
#'
#' Hidden function to download static tables from MPMP database
#'
#' @import stringr rvest dplyr
#' @importFrom janitor row_to_names
#' @param url MPMP url.
#'
#' @keywords internal
.nonclickable <- function (url){
  df <- rvest::read_html(url) %>% rvest::html_elements("tr") %>%
    rvest::html_text2() %>% as.data.frame() %>% dplyr::pull(1)
  if (length(df) < 2) {
    message(glue::glue("No table found for the given URL: {url}\n"))
  }
  else {
    pattern <- "PF3D7.*|^MAL*|PF.*|.pre-tRNA-|RNAzID|U5RNA|ORF|Surf|.*t000.*|.*m000.*|PF.*|3D7surf.*|Pf[0-9]+|[0-9]+\\.t00[0-9]+|1396.pre-trna-gly-1|1981.m00174|MaL13P1.80|pf14_0741"
    df <- df[grep("Annotation|PfID|PlasmoDB|Gene ID|NTac|Protein ID|description|Gene ID",
                  df, ignore.case = TRUE):length(df)] %>% stringr::str_replace_all("\t+", ":") %>%
      stringr::str_replace_all("^:+|:+$", "") %>%
      trimws() %>%
      purrr::discard(~. == "") %>%
      #grep("PF3D7*|MAL*|PC*|*tRNA*|ORF|RNA|pf|Surf|*t000*|*m000*", ., ignore.case = TRUE, value = TRUE) %>%
      grep(pattern, ., ignore.case = TRUE, value = TRUE) %>%
      stringr::str_split(":", simplify = TRUE) %>% as.data.frame()

    ## getting only PFIds
    df <- df[,which(sapply(df, function(column) any(grepl(pattern, column))))]
    ## The the returned object has more than one column itneeds to be processed
    df <- if(is.data.frame(df)){df[sapply(df, function(column) grepl(pattern, column))]} else{df}
    df <- df[!grepl("Annotation|PfID|PlasmoDB|Gene ID|NTac|Protein ID|description|Gene ID",df)]

    ##remove longer strings:
    df <- df[!(nchar(df) >= 25)]
    #colnames(df) <- c("PfID", "Annotation")
    message(glue::glue("\033[0;32mfetched successfully: {url}\033[0m\n"))
    return(unique(df))
  }
}

#' Fastest IfElse function: https://github.com/ICJIA/r-user-group/issues/11
#'
#' @keywords internal

.fast_ifelse2 <- function(test, yes, no) {
  stopifnot(identical(class(yes), class(no)))

  out <- rep(NA, length(test))
  out[test] <- yes
  out[!test] <- no
  class(out) <- class(yes)

  out
}

#' plasmoRUtils
#'
#' Hidden function called by \code{easytopGO} to run topGO.
#'
#' @import stringr topGO dplyr
#' @importFrom pRoloc goIdToTerm
#' @importFrom methods new
#'
#' @keywords internal
.usetopGO <- function(stats=stats,category=category,geneID=geneID,gene_2_GO=gene_2_GO,algo=algo,fdr=fdr,correction=correction){
  GOdata <- methods::new('topGOdata', ontology = category, allGenes = geneID, annot = topGO::annFUN.gene2GO, gene2GO = gene_2_GO, geneSel = function(allScore) allScore < 0.05)
  results <- topGO::runTest(GOdata, algorithm = algo, statistic = stats)

  goEnrichment <- topGO::GenTable(GOdata, stats = results, orderBy = 6, topNodes = length(topGO::usedGO(GOdata)))
  goEnrichment[["stats"]] <- as.numeric(goEnrichment[["stats"]] )
  ## Idea of FDR borrowed from: https://github.com/federicomarini/pcaExplorer/issues/5 and mosdef::run_topGO()
  if (fdr) {
    goEnrichment[["padj"]] <- p.adjust(goEnrichment[["stats"]], method = correction)
  }

  if (nrow(goEnrichment) == 0) return(message("No enriched term found within significant threshold"))



  get.genes <- plyr::ldply(topGO::genesInTerm(GOdata, goEnrichment$GO.ID), rbind) %>%
    tidyr::unite("col", 2:ncol(.), sep = ",", na.rm = TRUE)
  goEnrichment$associated_genes <- get.genes$col

  ## get full term for incomplete go description
  incomplete <- goEnrichment$GO.ID[grep("*\\.\\.\\.$",goEnrichment$Term)]
  goEnrichment$Term[grep("*\\.\\.\\.$",goEnrichment$Term)] <- pRoloc::goIdToTerm(incomplete, names = TRUE, keepNA = TRUE)

  goEnrichment$Term <- goEnrichment$Term %>%
    factor(levels = goEnrichment$Term) %>%
    # paste(goEnrichment$GO.ID, ., sep = ", ") %>%
    factor(levels = rev(.))
  goEnrichment <- goEnrichment %>% dplyr::filter(stats < 0.05)

  colnames(goEnrichment)[colnames(goEnrichment) == "stats"] <- stats

  return(goEnrichment)
}

cleanFun <- function(htmlString) {
  return(gsub("<.*?>", "", htmlString))
}

convert_last_letter <- function(str) {
  sub("([A-Za-z])$", "\\L\\1", str, perl = TRUE)
}

utils::globalVariables(".")
