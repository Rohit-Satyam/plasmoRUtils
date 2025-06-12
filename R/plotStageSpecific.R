#' Making plots similar to malaria.tools
#'
#' This function retrieves data from malaria.tools and generates expression value plots (in TPM) similar to those produced by the website. Use this function to create publication-ready plots.
#'
#' @importFrom ggsci scale_fill_nejm
#' @importFrom plotly ggplotly
#' @import dplyr
#' @import stringr
#' @import ggplot2
#' @importFrom glue glue
#' @import rvest
#' @export
#'
#' @param geneID Gene ID of Plasmodium falciparum or Plasmodium berghi.
#' @param returnData Logical. Use true to return dataframe used for making plots.
#' @param plotify To make plots interactive using plotly.
#'
#' @return A plot (or data theirof) of TPM values across multiple stages of parasite.
#' @examples
#' \dontrun{
#'   geneID <- c("PBANKA_0100600", "PBANKA_0102900", "PF3D7_0102900")
#'   ## To get Plot similar to malaria.tools
#'   res <- plotTissueSpecific(geneID = "PBANKA_0100600")
#' }
#'
plotStageSpecific <- function(geneID, returnData = FALSE, plotify = FALSE) {
  data("malariatools",envir = environment())

  index <- dplyr::filter(malariatools, gene %in% geneID) %>%
    dplyr::pull(2)

  url <- ifelse(grepl("PF3D7", geneID),
                glue::glue("https://malaria.sbs.ntu.edu.sg/profile/download/plot/{index}/2"),
                glue::glue("https://malaria.sbs.ntu.edu.sg/profile/download/plot/{index}/1")
  )

  line <- readLines(url, warn = FALSE) %>%
    stringr::str_replace_all("Gametocyte \\(male\\)", "Gametocyte(male)") %>%
    stringr::str_replace_all("Gametocyte \\(female\\)", "Gametocyte(female)")

  components <- unlist(stringr::str_split(line, " "))

  df <- matrix(components[-(1:4)], ncol = 4, byrow = TRUE) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::rename_with(~ components[1:4]) %>%
    dplyr::mutate(across(c(mean, min, max), as.numeric)) %>%
    dplyr::mutate(condition = factor(condition, levels = unique(condition)))

  if (returnData) {
    return(df)
  }

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = condition, y = mean, fill = condition)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.9), width = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = min), position = ggplot2::position_dodge(width = 0.9), alpha = 0.5, shape = 20, size = 2.5, color = "#636363") +
    ggplot2::geom_point(ggplot2::aes(y = max), position = ggplot2::position_dodge(width = 0.9), alpha = 0.5, shape = 20, size = 2.5, color = "#636363") +
    ggplot2::labs(x = "", y = "TPM") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggsci::scale_fill_nejm()

  if (plotify) {
    return(plotly::ggplotly(plot))
  } else {
    return(plot)
  }
}

