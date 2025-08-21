#' Fetch articles for given gene IDs from Google Scholar
#'
#' This function searches the Google Scholar corpus recursively for the articles that contains your Gene ID of interest.
#'
#' **Warning**: Scraping Google Scholar is against their [Terms of Service](https://scholar.google.com/intl/en/scholar/about.html). We advise users to use this function for querying few IDs (not more than 20) per day.
#' Proceeding with this function may result in your IP being blocked temporarily.
#'
#' @import dplyr
#' @import rvest
#' @import utils
#' @importFrom plyr ldply
#' @importFrom purrr map_dfr
#' @importFrom stringr str_extract
#' @importFrom polyglotr google_translate
#' @import tibble tibble
#' @export
#'
#' @param geneIDs Character vector of Gene IDs. If you want to use gene symbols, use organism name alongside to avoid articles that might have similar abbreviated word. eg. use "AP2-P AND Plasmodium".
#' @param year_start Limit the results to starting year of interest.
#' @param year_end Limit the results to end year of interest.
#' @param max_pages Maximum number of pages to scrap.
#' @param verbose Print warnings.
#' @param translate Translate the paper titles to english or desired language. Use two letter code.eg: "fr" for french, "en" for english and "es" for spanish.
#' @return A data frame, containing 5 columns: GeneID, Title of the article, Year of Publication, Url and Authors.
#' @examples
#' \dontrun{
#' ## We have a fake ID: PF3D7_0420300OR
#' res <- searchGSC(
#' geneIDs=c("PF3D7_0420300 OR MAL4P1.192 OR Q8I1N6 OR PFD0985w","PF3D7_0621000","PF3D7_0420300OR"),
#' translate = "en",
#' year_start = 2018, 
#' year_end   = 2021)
#' 
#' test <- searchGSC(geneID = c("AP2-P AND Plasmodium", "AP2-I"))
#' }
#'
searchGSC <- function(geneIDs,
                      year_start = NULL,
                      year_end = NULL,
                      max_pages = 2,
                      sleep_secs = 10,
                      verbose = TRUE, translate = NULL) {
  require(rvest)

  stopifnot(
    is.null(year_start) || is.numeric(year_start),
    is.null(year_end) || is.numeric(year_end)
  )

  # empty tibble with stable schema
  empty_tbl <- tibble::tibble(
    Query   = character(),
    Title   = character(),
    Url     = character(),
    Authors = character(),
    Year    = integer()
  )

  # Build URL and handle year specific queries. Also, if user use + sign rather than AND or other reserved characters, setting reserved = TRUE
  build_url <- function(q_raw, start_idx, ylo = NULL, yhi = NULL) {
    q <- utils::URLencode(q_raw, reserved = TRUE)
    if (is.null(ylo) && is.null(yhi)) {
      sprintf(
        "https://scholar.google.com/scholar?as_vis=1&hl=en&as_sdt=0,5&q=%s&start=%d",
        q, start_idx
      )
    } else {
      # If only one bound provided, use it for both (exact year). Setting as_vis=1 to avoid citations. Add as_rr=1 if you wish to modify this function to get review articles only
      if (is.null(ylo) && !is.null(yhi)) ylo <- yhi
      if (!is.null(ylo) && is.null(yhi)) yhi <- ylo
      sprintf(
        paste0(
          "https://scholar.google.com/scholar?",
          "as_vis=1&hl=en&as_sdt=0,5&q=%s&as_ylo=%d&as_yhi=%d&start=%d"
        ),
        q, as.integer(ylo), as.integer(yhi), start_idx
      )
    }
  }

  # Fetch & parse one page. Always returns a tibble.
  fetch_page <- function(query, start_idx, ylo = NULL, yhi = NULL) {
    url <- build_url(query, start_idx, ylo, yhi)

    page <- tryCatch(
      rvest::read_html(url),
      error = function(e) {
        if (verbose) {
          message(sprintf(
            "Read error at start=%d for '%s': %s",
            start_idx, query, conditionMessage(e)
          ))
        }
        return(NULL)
      }
    )
    if (is.null(page)) {
      return(empty_tbl)
    }

    # Result blocks (ignore sidebars/other)
    results <- rvest::html_nodes(page, ".gs_r .gs_ri")
    if (length(results) == 0) {
      return(empty_tbl)
    }

    # Titles: prefer anchor text; fall back to container text
    titles <- results %>%
      rvest::html_node(".gs_rt a") %>%
      rvest::html_text2()
    # For items without a link fill from .gs_rt container
    missing <- which(is.na(titles) | titles == "")
    if (length(missing)) {
      titles_all <- results %>%
        rvest::html_node(".gs_rt") %>%
        rvest::html_text2()
      titles[missing] <- titles_all[missing]
    }

    links <- results %>%
      rvest::html_node(".gs_rt a") %>%
      rvest::html_attr("href")

    meta <- results %>%
      rvest::html_node(".gs_a") %>%
      rvest::html_text2() ## Better leave the author string preprocessing.

    years <- stringr::str_extract(meta, "\\b\\d{4}\\b") %>% as.integer()

    tibble::tibble(
      Query   = rep(query, length(titles)),
      Title   = titles %||% character(length(titles)),
      Url     = links %||% rep(NA_character_, length(titles)),
      Authors = meta %||% rep(NA_character_, length(titles)),
      Year    = years
    )
  }

  # loop over each query and paginate until a page yields no new rows
  all <- purrr::map_dfr(geneIDs, function(gid) {
    seen_titles <- character()
    acc <- list()

    for (i in seq_len(max_pages)) {
      start_idx <- (i - 1) * 10
      dat <- fetch_page(gid, start_idx, year_start, year_end)

      # stop if page is empty
      if (nrow(dat) == 0) {
        if (i == 1L && verbose) {
          message(sprintf("No results for '%s'.", gid))
        }
        break
      }

      # keep only genuinely new titles (avoid dupes across pages)
      dat <- dat[!(dat$Title %in% seen_titles) & !is.na(dat$Title) & nzchar(dat$Title), , drop = FALSE]
      if (nrow(dat) == 0) break

      seen_titles <- c(seen_titles, dat$Title)
      acc[[length(acc) + 1L]] <- dat

      Sys.sleep(sleep_secs) # Increase the sleep time. helps avoid CAPTCHAs/rate limits
    }

    if (length(acc)) dplyr::bind_rows(acc) else empty_tbl
  })

  # If absolutely nothing came back, return empty with right columns
  if (nrow(all) == 0) {
    return(empty_tbl)
  }

  # Translate paper titles to desired language
  if (!is.null(translate)) {
    require(polyglotr)
    all$translation <- polyglotr::google_translate(all$Title, target_language = translate)
  }

  # De-dup and order; guard all steps
  all %>%
    dplyr::distinct(Title, .keep_all = TRUE) %>%
    dplyr::arrange(Query, dplyr::desc(Year))
}

#styler:::style_active_file()
