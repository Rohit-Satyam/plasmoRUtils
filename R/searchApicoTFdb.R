#' Fetch data from ApicoTFdb
#'
#' This function provides ability to query your gene IDs to ApicoTFDb.
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
#' pk: Plasmodium knowlesi
#' py: Plasmodium yoelii
#' pc: Plasmodium chabaudi
#'
#'By Other Apicomplexan Species
#'tg49: Toxoplasma Gondii ME49
#'tg89: Toxoplasma Gondii P89
#'cp: Cryptosporidium parvum
#'em: Eimeria maxima
#'bb: Babesia bovis
#'et: Eimeria tenella
#'nu: Neurospora caninum
#'cy: Cyclospora cayetanensis
#'
#'
#' @return df This function returns a dataframe of transcription regulators for organism of interest from ApicoTFDb.
#' @examples
#' \dontrun{
#' test <- searchApicoTFdb(org="pf")
#' }
#'


searchApicoTFdb <- function(org="pf"){
  urls <- list(
    pb = "https://bioinfo.icgeb.res.in/PtDB/htmls/berghii_new.html",
    pv = "https://bioinfo.icgeb.res.in/PtDB/htmls/vivax_new.html",
    pf = "https://bioinfo.icgeb.res.in/PtDB/htmls/pfal_new.html",
    pk = "https://bioinfo.icgeb.res.in/PtDB/htmls/knw_new.html",
    pc = "https://bioinfo.icgeb.res.in/PtDB/htmls/chas_new.html",
    py = "https://bioinfo.icgeb.res.in/PtDB/htmls/yoelii_new.html",
    tg49 = "https://bioinfo.icgeb.res.in/PtDB/htmls/t_gondii_49.html",
    tg89 = "https://bioinfo.icgeb.res.in/PtDB/htmls/t_gondii_new_89.html",
    cp = "https://bioinfo.icgeb.res.in/PtDB/htmls/c_parvum_new.html",
    em = "https://bioinfo.icgeb.res.in/PtDB/htmls/eimeria_maxima.html",
    bb = "https://bioinfo.icgeb.res.in/PtDB/htmls/babaesia.html",
    et = "https://bioinfo.icgeb.res.in/PtDB/htmls/emeria_tellena.html",
    nu = "https://bioinfo.icgeb.res.in/PtDB/htmls/neurospora.html",
    cy = "https://bioinfo.icgeb.res.in/PtDB/htmls/cyclospora.html"
  )

  url <- urls[[org]]
  webpage <- rvest::read_html(url)
  table <- webpage %>%
    rvest::html_nodes("table") %>%
    rvest::html_table(fill = TRUE) %>%
    .[[1]] %>%
    dplyr::filter(!apply(., 1, function(row) all(row == ""))) %>%
    dplyr::filter(!stringr::str_detect(`Gene ID`, "Gene ID"))

  return(table)
}
