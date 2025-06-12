#' Fetch Kinases and Phosphatases from Kipho Database
#'
#' This function provides ability to query your gene IDs to KiPho Database.
#'
#' @import dplyr
#' @import rvest
#' @export
#'
#' @param org Abbreviation of organism of interest.
#' Plasmodium Species:
#' pb: Plasmodium berghii
#' pv: Plasmodium vivax
#' pf: Plasmodium falciparum
#' pc: Plasmodium chabaudi
#' @param type Type of protein class i.e. "kinase" or "phosphatase". Default: "kinase"
#'
#'
#' @return df This function returns a dataframe of kinases/phosphatases in Plasmodium species.
#' @examples
#' \dontrun{
#' test <- searchKipho(org="pf")
#' }
#'


searchKipho <- function(org = "pf", type = "kinase") {
  baseurl <- "https://bioinfo.icgeb.res.in/kipho/"
  urls <- list(
    kinase = list(
      pb = glue::glue("{baseurl}kinase_PBANKA.php"),
      pv = glue::glue("{baseurl}kinase_PVX.php"),
      pf = glue::glue("{baseurl}kinase_PF3D7.php"),
      pc = glue::glue("{baseurl}kinase_PCHAS.php")
    ),
    phosphatase = list(
      pb = glue::glue("{baseurl}phosphatase_PBANKA.php"),
      pv = glue::glue("{baseurl}phosphatase_PVX.php"),
      pf = glue::glue("{baseurl}phosphatase_PF3D7.php"),
      pc = glue::glue("{baseurl}phosphatase_PCHAS.php")
    )
  )

  if (type %in% c("kinase", "phosphatase")) {
    url <- urls[[type]][[org]]
    table <- rvest::read_html(url) %>%
      rvest::html_nodes("table") %>%
      rvest::html_table(fill = TRUE) %>%
      .[[1]] %>%
      dplyr::filter(!apply(., 1, function(row) all(row == ""))) %>%
      dplyr::filter(!stringr::str_detect(`Gene ID`, "Gene ID")) %>%
      dplyr::select(-tail(names(.), 5))

    return(table)
  } else {
    message("Invalid type argument. Use type = 'kinase' or type = 'phosphatase'.")
  }
}
