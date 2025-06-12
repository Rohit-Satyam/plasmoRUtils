#' Quick ORA result plotting
#'
#' A convenience function to quickly plot the results of \code{easytopGO} or ORA results obtained from VEupathDB component databases.
#'
#' @import dplyr
#' @importFrom scales label_wrap
#' @importFrom rlang .data
#' @import ggplot2 topGO
#' @export
#'
#' @param res Output of \code{easytopGO}.
#' @param title Title of the plot.
#' @param limit No of terms to plot.
#' @param desc Column name containing the GO description.If using VEuPathDB compnent database, this would be "Name" column.
#' @param genecounts Column name containing number of query genes associated with GO description. Used for bubble size. For VEupathDB this would be "Result count" column.
#' @param sortby Name of the statistics column to sort the terms by and use for plotting.
#'
#' @return A publication ready ggplot2 object.
#' @examples
#' \dontrun{
#' ## get enrichment results from easytopGO
#' gores <- easytopGO(geneID = geneList, bkggset = background.gset, stats = "ks")
#' baseurl <- "https://plasmodb.org/common/downloads/Current_Release/"
#' url <- paste0(baseurl, "Pfalciparum3D7/gaf/PlasmoDB-68_Pfalciparum3D7_GO.gaf.gz")
#' gores <- easytopGO(
#'   geneID = geneList, useGAF = TRUE, useBiomart = FALSE, gaf = url,
#'   bkggset = background.gset, category = "BP", stats = "ks"
#' )
#'
#' ## Then feed them to easyGOPlot
#' plot <- easyGOPlot(gores, title = "GO Enrichment Biological Processes")
#' }
#'
easyGOPlot <- function(res, title = "GO Biological processes", limit = 20, desc = "Term", genecounts = "Annotated", sortby = "ks") {
  ## To take care of marginal cases where max and min no. of genes are equal then take the unique values for ranges. Also limit plotting everything
  ranges <- unique(round(seq(min(res[, genecounts]), max(res[, genecounts]),
    length.out = ifelse(nrow(res) < 3, nrow(res), 3)
  )))

  ## Handling duplicated GO terms using the GO description and other columns. We see a lot of them when using VEUpathDB HAF files

  res <- res %>% .[!duplicated(.[, c(colnames(.)[2:ncol(.)])]), ]

  ## Transform P-values
  res[[sortby]] <- -log10(as.numeric(res[[sortby]]))
  col <- sortby

  ## Prevent plotting everything
  res[[desc]] <- factor(res[[desc]], levels = res[order(res[[sortby]], decreasing = FALSE), ][[desc]])
  res <- res %>%
    slice_head(n = ifelse(nrow(res) < limit, nrow(res), limit))

  ggplot2::ggplot(res, aes(
    x = .data[[desc]], y = .data[[col]],
    size = .data[[genecounts]], fill = .data[[col]]
  )) +
    ggplot2::geom_hline(
      yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
      linetype = c("dotted", "longdash", "solid"),
      colour = c("black", "black", "black"),
      size = c(0.5, 1, 1.5), alpha = 0.5
    ) +
    ggplot2::expand_limits(y = 1) +
    ggplot2::geom_point(shape = 21) +
    ggplot2::scale_size(
      range = c(3, 8), name = "No. of Genes",
      breaks = ranges,
      labels = ranges
    ) +
    ggplot2::scale_fill_continuous(low = "#BDD2E6", high = "#7d1d79") +
    ggplot2::xlab("") +
    ggplot2::ylab("Enrichment score") +
    ggplot2::labs(
      title = title,
      subtitle = paste("Top", length(res[[desc]]), "terms ordered by", col, "p-value"),
      caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001",
      fill = paste0("-log10(", col, ")")
    ) +
    ggplot2::theme_bw(base_size = 24) +
    ggplot2::theme(
      legend.position = "right",
      legend.background = ggplot2::element_rect(),
      plot.title = ggplot2::element_text(angle = 0, size = 11, face = "bold", vjust = 1),
      plot.subtitle = ggplot2::element_text(angle = 0, size = 10, face = "bold", vjust = 1),
      plot.caption = ggplot2::element_text(angle = 0, size = 10, face = "bold", vjust = 1),
      axis.text.x = ggplot2::element_text(angle = 0, size = 10, face = "bold", hjust = 1.1),
      axis.text.y = ggplot2::element_text(angle = 0, size = 10, face = "bold", vjust = 0.5),
      axis.title = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.key = ggplot2::element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.text = ggplot2::element_text(size = 11, face = "bold"),
      title = ggplot2::element_text(size = 11, face = "bold")
    ) +
    ggplot2::guides(color = guide_legend(paste("-log10(", col, ")\n"))) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = scales::label_wrap(40))
}
#styler:::style_active_file()
