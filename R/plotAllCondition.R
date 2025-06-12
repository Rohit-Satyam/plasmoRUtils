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
#' @param geneID Single Gene ID of Plasmodium falciparum or Plasmodium berghi.
#' @param returnData Logical. Use true to return dataframe used for making plots.
#' @param plotify To make plots interactive using plotly.
#'
#' @return A plot (or data theirof) of TPM values across multiple stages of parasite.
#' @examples
#' \dontrun{
#'   #'   ## To get Plot similar to malaria.tools
#'   res <- plotAllCondition(geneID = "PBANKA_0100600")
#' }
#'
plotAllCondition <- function(geneID, returnData = FALSE, plotify = FALSE) {
  data("malariatools",envir = environment())

  index <- malariatools %>%
    dplyr::filter(gene %in% geneID) %>%
    dplyr::pull(2)

  url <- glue::glue("https://malaria.sbs.ntu.edu.sg/profile/download/plot/{index}")

  df <- read.table(url, sep = "\t", header = TRUE) %>%
    mutate(group = stringr::str_split(condition, pattern = "[ :,]", n = 2, simplify = TRUE)[, 1])

  if (returnData) {
    return(df)
  }


  plot <- ggplot(df, aes(x = condition, y = mean, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    geom_point(aes(y = min), position = position_dodge(width = 0.9), alpha = 0.5, shape = 20, size = 2.5, color = "#636363") +
    geom_point(aes(y = max), position = position_dodge(width = 0.9), alpha = 0.5, shape = 20, size = 2.5, color = "#636363") +
    labs(x = "", y = "TPM") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
    ) +
    ggsci::scale_fill_nejm()

  if (plotify) {
    return(plotly::ggplotly(plot))
  } else {
    return(plot)
  }
}
