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
#'   geneID = "PfMEDA1",
#'   org = "Plasmodium falciparum",
#'   max_pages = 10
#' )
#'
#' hits <- searchRS(
#'   geneID = "PfAP2-I",
#'   org = "Plasmodium falciparum",
#'   gene_aliases = c(
#'     "Pf AP2-I",
#'     "pfap2-i",
#'     "PF3D7_1007700",
#'     "Apetala 2 Invasion",
#'     "AP2-I transcription factor"
#'   ),
#'   org_aliases = c(
#'     "P. falciparum",
#'     "3D7",
#'     "malaria parasite"
#'   ),
#'   max_pages = 10
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
  # -------------------------------------------------------
  # 1. Validate input arguments
  # -------------------------------------------------------

  if (!is.numeric(max_pages) ||
    length(max_pages) != 1L ||
    is.na(max_pages) ||
    !is.finite(max_pages) ||
    max_pages < 1 ||
    max_pages != floor(max_pages)) {
    stop("max_pages must be a positive integer.")
  }

  if (!is.numeric(sleep_sec) ||
    length(sleep_sec) != 1L ||
    is.na(sleep_sec) ||
    !is.finite(sleep_sec) ||
    sleep_sec < 0) {
    stop("sleep_sec must be a non-negative number.")
  }

  max_pages <- as.integer(max_pages)

  # -------------------------------------------------------
  # 2. Construct search query
  # -------------------------------------------------------

  query <- paste(geneID, org)

  # -------------------------------------------------------
  # 3. Fetch the first page only once
  # -------------------------------------------------------

  first <- .rs_api_search_page(
    query = query,
    offset = 0,
    article_type = article_type,
    status = status
  )

  # Determine the API pagination limit.
  # Fall back to 10 if the API returns an invalid value.

  page_limit <- suppressWarnings(
    as.integer(first$limit %||% 10L)
  )

  if (length(page_limit) != 1L ||
    is.na(page_limit) ||
    page_limit < 1L) {
    page_limit <- 10L
  }

  # Retain results from the initial API request.

  first_hits <- if (!is.null(first$data) &&
    length(first$data) > 0L) {
    tibble::as_tibble(first$data)
  } else {
    tibble::tibble()
  }

  # -------------------------------------------------------
  # 4. Generate offsets for remaining pages
  # -------------------------------------------------------

  # Example:
  # max_pages = 5, page_limit = 10
  #
  # First request:  offset 0
  # Remaining:      10, 20, 30, 40

  offsets <- if (max_pages > 1L) {
    page_limit * seq_len(max_pages - 1L)
  } else {
    integer(0)
  }

  # -------------------------------------------------------
  # 5. Fetch subsequent pages
  # -------------------------------------------------------

  remaining_hits <- purrr::map_dfr(offsets, function(off) {
    # Delay only between API requests.
    Sys.sleep(sleep_sec)

    page <- .rs_api_search_page(
      query = query,
      offset = off,
      article_type = article_type,
      status = status
    )

    # Return an empty tibble for missing/empty responses.
    if (is.null(page$data) || length(page$data) == 0L) {
      return(tibble::tibble())
    }

    tibble::as_tibble(page$data)
  })

  # -------------------------------------------------------
  # 6. Combine first page with remaining pages
  # -------------------------------------------------------

  raw_hits <- dplyr::bind_rows(
    first_hits,
    remaining_hits
  )

  # Stop processing if no articles were retrieved.
  if (nrow(raw_hits) == 0L) {
    return(raw_hits)
  }

  # -------------------------------------------------------
  # 7. Remove duplicate articles
  # -------------------------------------------------------

  raw_hits <- dplyr::distinct(
    raw_hits,
    article_identity,
    doi_version,
    .keep_all = TRUE
  )

  # -------------------------------------------------------
  # 8. Construct flexible regex patterns
  # -------------------------------------------------------

  gene_rx <- .make_flexible_regex(
    c(geneID, gene_aliases)
  )

  org_rx <- .make_flexible_regex(
    c(org, org_aliases)
  )

  # -------------------------------------------------------
  # 9. Retrieve article text and identify matches
  # -------------------------------------------------------

  out <- dplyr::mutate(
    raw_hits,
    full_url = paste0(
      "https://www.researchsquare.com",
      .data$url
    ),
    article_text = purrr::map_chr(
      .data$url,
      ~ tryCatch(
        .rs_article_text(.x),
        error = function(e) NA_character_
      )
    ),
    search_blob = paste(
      .data$title,
      .data$article_text,
      sep = "\n\n"
    ),
    has_gene = stringr::str_detect(
      .data$search_blob,
      gene_rx
    ),
    has_org = stringr::str_detect(
      .data$search_blob,
      org_rx
    ),
    keep = .data$has_gene & .data$has_org
  )

  # -------------------------------------------------------
  # 10. Retain articles matching both gene and organism
  # -------------------------------------------------------

  out <- dplyr::filter(
    out,
    .data$keep
  )

  # -------------------------------------------------------
  # 11. Return selected columns
  # -------------------------------------------------------

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


# ---------------------------------------------------------
# NULL-coalescing helper
# ---------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) ||
    length(x) == 0L ||
    (length(x) == 1L && is.na(x))) {
    y
  } else {
    x
  }
}
