#' Fetch articles for given gene IDs from Google Scholar
#'
#' This function searches the Google Scholar corpus for the articles that contains your Gene ID of interest.
#'
#' **Warning**: Scraping Google Scholar is against their [Terms of Service](https://scholar.google.com/intl/en/scholar/about.html). We advise users to use this function for querying few IDs (not more than 20) per day.
#' Proceeding with this function may result in your IP being blocked temporarily.
#'
#' @import dplyr
#' @import rvest
#' @importFrom plyr ldply
#' @import tibble tibble
#' @export
#'
#' @param geneID Character vector of Gene IDs. If you want to use gene symbols, use organism name alongside to avoid articles that might have similar abbreviated word. eg. use "AP2-P AND Plasmodium".
#'
#' @return A data frame, containing 5 columns: GeneID, Title of the article, Year of Publication, Url and Authors.
#' @examples
#' \dontrun{
#' test <- searchGSC(c("PF3D7_0420300", "PF3D7_0621000"))
#' test <- searchGSC(geneID = c("AP2-P AND Plasmodium","AP2-I"))
#' }
#'
searchGSC <- function(geneID) {
  fetch_results <- function(gid) {
    url <- stringr::str_c("https://scholar.google.com/scholar?q=", utils::URLencode(gid))
    page <- tryCatch(rvest::read_html(url), error = function(e) return(NULL))
    if (is.null(page)) return(NULL)

    entries <- page %>% rvest::html_nodes(".gs_r")

    if (length(entries) == 0) return(NULL)

    results <- lapply(entries, function(entry) {
      # Title (works with or without <a> tag)
      title_node <- entry %>% html_node(".gs_rt")
      title <- title_node %>% html_text(trim = TRUE)%>%
        stringr::str_remove("^(\\[[^]]+\\])+\\s*") %>%
        stringr::str_trim()

      # Link (may be missing)
      link_node <- title_node %>% html_node("a")
      link <- if (!is.na(link_node)) html_attr(link_node, "href") else NA

      # Author info & Year
      author_info <- entry %>% html_node(".gs_a") %>% html_text(trim = TRUE)
      year <- stringr::str_extract(author_info, "\\b\\d{4}\\b")
      authors <- stringr::str_remove(author_info, "\\s*-\\s*.*")  # Remove everything after dash

      tibble(
        GeneID = gid,
        Title = title,
        Year = year,
        Url = link,
        Authors = authors
      )
    })

    bind_rows(results)
  }


  combined <- purrr::map_dfr(geneID, fetch_results) %>%
    .[!(is.na(.$Title)),] %>%
    dplyr::distinct(Title, .keep_all = TRUE) %>%
    dplyr::mutate(Title = stringr::str_trim(stringr::str_replace_all(Title, "<.*?>", " "))) %>% dplyr::arrange(Year)

  return(combined)
}

#styler:::style_active_file()
