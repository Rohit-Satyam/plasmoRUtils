#' plasmoRUtils
#'
#' Hidden function to fetch the expandable tables from MPMP database
#'
#' @import rvest
#' @param url MPMP url.
#' @keywords internal

.clickable <- function(url) {
  df <- rvest::read_html(url) %>%
    rvest::html_elements("table.table-bordered.table-hover") %>%
    rvest::html_table()

  if (length(df) == 0) {
    message(paste0("No table found for the given URL: ", url, "\n"))
  } else {
    message(paste0("\033[0;32mfetched successfully: ", url, "\033[0m\n"))
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
    message(paste0("No table found for the given URL: ", url, "\n"))
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
    message(paste0("\033[0;32mfetched successfully: ", url, "\033[0m\n"))
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
#' Hidden function called by [easytopGO()] to run topGO.
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

#' plasmoRUtils
#'
#' Hidden function for searchRS
#'
#' @import httr2 jsonlite
#'
#' @keywords internal
.rs_api_search_page <- function(query,
                                offset = 0,
                                article_type = "Research Article",
                                status = "all") {
  resp <- httr2::request("https://www.researchsquare.com/api/search") |>
    httr2::req_url_query(
      articleType = article_type,
      offset = offset,
      status = status,
      unified = query
    ) |>
    httr2::req_user_agent("plasmoRUtils searchRS") |>
    httr2::req_timeout(60) |>
    httr2::req_retry(
      max_tries = 3,
      retry_on_failure = TRUE,
      is_transient = function(resp) {
        httr2::resp_status(resp) %in% c(408, 429, 500, 502, 503, 504)
      }
    ) |>
    httr2::req_perform()
  
  txt <- httr2::resp_body_string(resp)
  out <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  out$result
}

#' plasmoRUtils
#'
#' Hidden function for searchRS 
#'
#'
#' @keywords internal
.escape_regex <- function(x) {
  vapply(strsplit(x, "", fixed = TRUE), function(chars) {
    meta <- chars %in% c(
      "\\", ".", "|", "(", ")", "[", "]", "{", "}",
      "^", "$", "*", "+", "?"
    )
    
    paste0(ifelse(meta, paste0("\\", chars), chars), collapse = "")
  }, character(1))
}


#' plasmoRUtils
#'
#' Hidden function for searchRS
#'
#' @import purrr
#'
#' @keywords internal
.make_flexible_regex <- function(terms) {
  terms <- unique(terms[!is.na(terms) & nzchar(terms)])
  
  if (length(terms) == 0) {
    return("(?!)")
  }
  
  term_regex <- purrr::map_chr(terms, function(x) {
    parts <- unlist(strsplit(x, "[-_[:space:]]+"))
    parts <- .escape_regex(parts)
    
    paste(parts, collapse = "[-_[:space:]]*")
  })
  
  paste0(
    "(?i)(?<![A-Za-z0-9])(?:",
    paste(term_regex, collapse = "|"),
    ")(?![A-Za-z0-9])"
  )
}

#' plasmoRUtils
#'
#' Hidden function for searchRS
#'
#' @import stringr httr2 rvest
#'
#' @keywords internal
.rs_article_text <- function(url_path) {
  url <- ifelse(
    stringr::str_starts(url_path, "https?://"),
    url_path,
    paste0("https://www.researchsquare.com", url_path)
  )
  
  resp <- httr2::request(url) |>
    httr2::req_user_agent("plasmoRUtils searchRS") |>
    httr2::req_timeout(60) |>
    httr2::req_retry(
      max_tries = 3,
      retry_on_failure = TRUE,
      is_transient = function(resp) {
        httr2::resp_status(resp) %in% c(408, 429, 500, 502, 503, 504)
      }
    ) |>
    httr2::req_perform()
  
  txt <- httr2::resp_body_string(resp)
  
  html <- rvest::read_html(txt)
  body <- rvest::html_element(html, "body")
  
  rvest::html_text2(body)
}

#' plasmoRUtils
#'
#' Hidden function to get tables from VEuPathDB databases with authentication
#'
#' @import httr2 readr
#'
#' @keywords internal
#' 
.read_authenticated_tsv <- function(url,api_key) {
  
  resp <- httr2::request(url) %>%
    httr2::req_auth_bearer_token(api_key) %>%
    httr2::req_headers(
      Accept = "text/plain"
    ) %>%
    httr2::req_perform()
  
  txt <- httr2::resp_body_string(resp)
  
  readr::read_tsv(
    I(txt),
    show_col_types = FALSE,
    progress = FALSE
  )
}
