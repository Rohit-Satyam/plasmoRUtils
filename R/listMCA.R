#' List data sets in MCA
#'
#' This function lists all the datasets available on Malaria Cell Atlas Database.
#'
#' @import dplyr
#' @import rvest
#' @importFrom rlang .data
#'
#' @return df A dataframe listing all the datasets available on MCA and URls that can be used with \code{easyMCA} function.
#' @export
#' @examples
#' \dontrun{
#'   df <- listMCA()
#' }
#'
listMCA <- function() {
  # Read and parse the table from the MCA website
  temp <- rvest::read_html("https://www.malariacellatlas.org/data-sets/") %>%
    rvest::html_table() %>%
    .[[1]] %>%
    dplyr::filter(stringr::str_detect(.data$`Description and links`, pattern = "Download"))

  # Extract and clean up the description and links
  temp2 <- temp$`Description and links` %>%
    stringr::str_split(pattern = "\n", simplify = TRUE) %>%
    .[, 1:2] %>%
    trimws() %>%
    as.data.frame() %>%
    dplyr::rename(Set = 1, Description = 2)

  # Combine the cleaned data with the original table
  temp <- temp %>%
    dplyr::select(-`Description and links`) %>%
    dplyr::bind_cols(temp2) %>%
    dplyr::mutate(links = paste0(
      "https://www.malariacellatlas.org",
      rvest::read_html("https://www.malariacellatlas.org/data-sets/") %>%
        rvest::html_nodes(xpath = "//a[text()='Download data']") %>%
        rvest::html_attr("href")
    ))

  # View and return the final table
  View(temp)
  return(temp)
}


