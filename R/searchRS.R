#' Fetch articles from Research Square database
#'
#' This function searches the Research Square corpus for the articles that contains your Gene ID/ Gene name of interest.
#'
#' @import dplyr
#' @import httr2
#' @import rvest
#' @import stringr
#' @import tibble
#' @importFrom jsonlite fromJSON
#' @importFrom purrr map_chr
#' @export
#'
#' @param geneID Gene ID such as PF3D7_XXXX or Gene Symbol.
#' @param org Scientific name of the organism. Eg. \emph{Plasmodium falciparum}
#' @param article_type The type of Article to return. Possible values: "Research Article", "Systematic Review", "Method Article", "Short Report", "Case Report", "Data Note", "Video".
#' @param sleep_sec Sleep time between subsequent queries.
#' @param gene_aliases Gene aliases, short form or alternate names.
#' @param org_aliases Common name or aliases to filter irrelevant articles.
#' @param max_pages Maximum number of pages to parse.
#' @return A data frame, containing 10 columns such as Research square ID, Title, Authors, Date, Status, Title of the journal, Article type, URLs, Gene ID present (logical), Organism name present (logical).
#' 
#' @examples
#' \dontrun{
#' 
#' hits <- searchRS(
#' geneID = "PfMEDA1",
#' org = "Plasmodium falciparum",
#' max_pages = 10
#' )
#' 
#' hits <- searchRS(
#' geneID = "PfAP2-I",
#' org = "Plasmodium falciparum",
#' gene_aliases = c(
#'   "Pf AP2-I",
#'   "pfap2-i",
#'   "PF3D7_1007700",
#'   "Apetala 2 Invasion",
#'   "AP2-I transcription factor"
#' ),
#' org_aliases = c(
#'   "P. falciparum",
#'   "3D7",
#'   "malaria parasite"
#' ),
#' max_pages = 10
#' )
#' }
#'



searchRS <- function(geneID,
                     org,
                     gene_aliases = character(),
                     org_aliases = character(),
                     max_pages = 5,
                     sleep_sec = 0.2,
                     status = "all",
                     article_type = "Research Article") {
  query <- paste(geneID, org)
  
  first <- .rs_api_search_page(
    query = query,
    offset = 0,
    article_type = article_type,
    status = status
  )
  
  page_limit <- first$limit %||% 10
  
  offsets <- seq(0, by = page_limit, length.out = max_pages)
  
  raw_hits <- purrr::map_dfr(offsets, function(off) {
    Sys.sleep(sleep_sec)
    
    page <- .rs_api_search_page(
      query = query,
      offset = off,
      article_type = article_type,
      status = status
    )
    
    if (is.null(page$data) || length(page$data) == 0) {
      return(tibble::tibble())
    }
    
    tibble::as_tibble(page$data)
  })
  
  if (nrow(raw_hits) == 0) {
    return(raw_hits)
  }
  
  raw_hits <- dplyr::distinct(
    raw_hits,
    article_identity,
    doi_version,
    .keep_all = TRUE
  )
  
  gene_rx <- .make_flexible_regex(c(geneID, gene_aliases))
  org_rx  <- .make_flexible_regex(c(org, org_aliases))
  
  out <- dplyr::mutate(
    raw_hits,
    full_url = paste0("https://www.researchsquare.com", .data$url),
    article_text = purrr::map_chr(
      .data$url,
      ~ tryCatch(
        .rs_article_text(.x),
        error = function(e) NA_character_
      )
    ),
    search_blob = paste(.data$title, .data$article_text, sep = "\n\n"),
    has_gene = stringr::str_detect(.data$search_blob, gene_rx),
    has_org = stringr::str_detect(.data$search_blob, org_rx),
    keep = .data$has_gene & .data$has_org
  )
  
  out <- dplyr::filter(out, .data$keep)
  
  dplyr::select(
    out,
    article_identity,
    title,
    authors,
    posted_at,
    status,
    journal_title,
    article_type,
    full_url,
    has_gene,
    has_org
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    y
  } else {
    x
  }
}
