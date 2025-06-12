#' Fetch data tables from malaria.tools database
#'
#' This function retrieves data from malaria.tools and generates a dataframe containing the Stage of Parasite in which the gene is highly expressed.
#'
#' @import dplyr
#' @import stringr
#' @importFrom glue glue
#' @import rvest
#' @export
#'
#' @param geneID A character vector of Gene IDs of Plasmodium falciparum or Plasmodium berghi.
#'
#' @return A plot (or data thereof) of TPM values across multiple stages of parasite.
#' @examples
#' \dontrun{
#'   geneID <- c("PBANKA_0100600", "PBANKA_0102900", "PF3D7_0102900")
#'   ## To get condition specificity and tissue specificity data
#'   res <- searchMT(geneID = geneID)
#' }
#'
searchMT <- function(geneID) {
  data("malariatools",envir = environment())

  index <- dplyr::filter(malariatools, gene %in% geneID) %>%
    dplyr::pull(2)

  ll <- purrr::map(index, function(x) {
    url <- glue::glue("https://malaria.sbs.ntu.edu.sg/profile/view/{x}")
    webpage <- rvest::read_html(url)

    col1 <- webpage %>%
      rvest::html_elements("h1") %>%
      rvest::html_text2() %>%
      stringr::str_remove("Expression profile for ") %>%
      stringr::str_split("[ ]", n = 2, simplify = TRUE) %>%
      .[, 1]

    col2 <- webpage %>%
      rvest::html_elements("p") %>%
      rvest::html_text2() %>%
      dplyr::first() %>%
      stringr::str_remove_all('"')

    col3 <- webpage %>%
      rvest::html_elements("span") %>%
      rvest::html_text2() %>%
      stringr::str_subset("Condition specificity|Tissue specificity") %>%
      .[1] %>%
      stringr::str_remove("Condition specificity: ")

    col4 <- webpage %>%
      rvest::html_elements("span") %>%
      rvest::html_text2() %>%
      stringr::str_subset("Condition specificity|Tissue specificity") %>%
      .[2] %>%
      stringr::str_remove("Tissue specificity: ")

    tibble::tibble(geneid = col1, description = col2, `Condition Specificity` = col3, `Tissue Specificity` = col4)
  })

  df <- dplyr::bind_rows(ll)
  return(df)
}

#styler:::style_active_file()
