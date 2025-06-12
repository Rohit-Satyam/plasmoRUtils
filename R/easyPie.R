#' Pie chart to summarize searchMT results
#'
#' This function make a donut chart to represent distribution of input gene IDs across different stages of Malaria given a result object from \code{plotTissueSpecific} function.
#'
#' @import echarts4r
#' @import dplyr tidyr
#' @export
#'
#' @param df A dataframe obtained from plotTissueSpecific(returnData=TRUE).
#' @param col Column to plot as donut chart. Default: "Tissue Specificity"
#' @return A plot (or data thereof) of domains present in the list of gene IDs.
#' @examples
#' \dontrun{
#'   geneID <- c("PBANKA_0100600", "PBANKA_0102900", "PF3D7_0102900")
#'   ## To get Plot similar to malaria.tools
#'   res <- searchMT(geneID = geneID)
#'   res %>% easyPie()
#' }
#'
easyPie <- function(df, col="Tissue Specificity"){
  stringr::str_split(df[[col]],pattern="[ ]",n = 2,simplify = TRUE) %>%
    `colnames<-`(c("Stages","discard"))%>%
    as.data.frame() %>%
    dplyr::group_by(Stages)  %>%
    dplyr::summarise(Genes=n()) %>%
    mutate(labels=paste0(Stages,"(",Genes,")")) %>%
    e_charts(labels) |>
    e_pie(Genes, radius = c("20%", "40%"),legend = FALSE) |>
    e_labels(position = "outside", fontSize = 17,fontWeight = "bold", color = "black", lineLength = 20, lineStyle = list(color = "#333", width = 1))
}
