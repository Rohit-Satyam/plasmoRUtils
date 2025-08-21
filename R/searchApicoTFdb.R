#' Fetch data from ApicoTFdb
#'
#' This function provides ability to query your gene IDs to [ApicoTFDb](https://bioinfo.icgeb.res.in/PtDB/index.html).
#'
#' @import dplyr
#' @import rvest
#' @export
#'
#' @param org Abbreviation of organism of interest.
#' @param fetch Describe the tables to be fetched. Default: "all" will fetch all TFs. To fetch TRs,CRRs,RNA-regs or Experimentally verified TF, use "trs","crrs","rnaregs" and "exptfs" respectively. When using other than "all", org argument will be ignored.
#'
#' \itemize{
#' \strong{Plasmodium Species:}
#'
#' \item pb: \emph{Plasmodium berghii}
#' \item pv: \emph{Plasmodium vivax}
#' \item pf: \emph{Plasmodium falciparum}
#' \item pk: \emph{Plasmodium knowlesi}
#' \item py: \emph{Plasmodium yoelii}
#' \item pc: \emph{Plasmodium chabaudi}
#' }
#'
#' \itemize{
#' \strong{By Other Apicomplexan Species}
#'
#' \item tg49: \emph{Toxoplasma Gondii} ME49
#' \item tg89: \emph{Toxoplasma Gondii} P89
#' \item cp: \emph{Cryptosporidium parvum}
#' \item em: \emph{Eimeria maxima}
#' \item bb: \emph{Babesia bovis}
#' \item et: \emph{Eimeria tenella}
#' \item nu: \emph{Neurospora caninum}
#' \item cy: \emph{Cyclospora cayetanensis}
#' }
#'
#' @return df This function returns a dataframe of transcription regulators for organism of interest from ApicoTFDb.
#' @examples
#' \dontrun{
#' test <- searchApicoTFdb(org = "pf")
#' }
#'
searchApicoTFdb <- function(org = "pf", fetch = "all") {
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
  if (fetch == "all") {
    url <- urls[[org]]
    webpage <- rvest::read_html(url)
    table <- webpage %>%
      rvest::html_nodes("table") %>%
      rvest::html_table(fill = TRUE) %>%
      .[[1]] %>%
      dplyr::filter(!apply(., 1, function(row) all(row == ""))) %>%
      dplyr::filter(!stringr::str_detect(`Gene ID`, "Gene ID"))

    return(table)
  } else if (fetch == "trrs") {
    url <- "https://bioinfo.icgeb.res.in/PtDB/htmls/TR.html"
    webpage <- rvest::read_html(url)
    table <- webpage %>%
      rvest::html_nodes("table") %>%
      rvest::html_table(fill = TRUE) %>%
      .[[1]] %>%
      .[!(stringr::str_detect(.$`Gene-ID`, "Gene ID")), ] %>%
      .[, -ncol(.)]
    return(table)
  } else if (fetch == "crrs") {
    url <- "https://bioinfo.icgeb.res.in/PtDB/htmls/CRR.html"
    webpage <- rvest::read_html(url)
    table <- webpage %>%
      rvest::html_nodes("table") %>%
      rvest::html_table(fill = TRUE) %>%
      .[[1]] %>%
      .[!(stringr::str_detect(.$`Gene-ID`, "Gene ID")), ] %>%
      .[, -ncol(.)]
    return(table)
  } else if (fetch == "rnaregs") {
    url <- "https://bioinfo.icgeb.res.in/PtDB/htmls/RNA.1.html"
    webpage <- rvest::read_html(url)
    table <- webpage %>%
      rvest::html_nodes("table") %>%
      rvest::html_table(fill = TRUE) %>%
      .[[1]] %>%
      .[!(stringr::str_detect(.$`Gene-ID`, "Gene ID")), ] %>%
      .[, -ncol(.)]
    return(table)
  } else if (fetch == "exptfs") {
    url <- "https://bioinfo.icgeb.res.in/PtDB/htmls/exp_new.html"
    webpage <- rvest::read_html(url)
    table <- webpage %>%
      rvest::html_nodes("table") %>%
      rvest::html_table(fill = TRUE) %>%
      .[[1]] %>%
      .[!(stringr::str_detect(.$`Gene_ID`, "Gene_ID")), ]
    return(table)
  }
}

#styler:::style_active_file()
