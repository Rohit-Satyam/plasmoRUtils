#' Plotting domains on protein body
#'
#' This function retrieves data from malaria.tools and generates a dataframe containing the Stage of Parasite in which the gene is highly expressed.
#'
#' @import dplyr
#' @import biomaRt
#' @importFrom glue glue glue_collapse
#' @import drawProteins
#' @export
#'
#' @param geneID A character vector of Gene IDs of Plasmodium falciparum or Plasmodium berghi.Remove version from the gene ids.
#' @param mart Name of Ensembl Biomart. Default: "protists_mart"
#' @param gset Gene-set of organism. This should be changed if you are dealing with other Plasmodium species. Default: "pfalciparum_eg_gene"
#' @param input Input id type. The gene IDs from PlasmoDB are "ensembl_gene_id".
#' @param fetchid Desired output ids
#' @param returnData To return converted gene ids with domain information fetched from uniprot.
#'
#' @return A plot (or data thereof) of domains present in the list of gene IDs.
#' @examples
#' \dontrun{
#'   ## Search proper mart
#'   ## View(listDatasets(biomaRt::useEnsemblGenomes(biomart = "protists_mart")))
#'   ## To get domain information from uniprot and prepare publication ready figures
#'   plot <- plotDomain(geneID = c("PF3D7_0518900", "PF3D7_0602800", "PF3D7_0624600"))
#'
#'   ## Currently pberghei doesn't work, see issue: https://github.com/grimbough/biomaRt/issues/110
#'   plot <- plotDomain(geneID = c("PBANKA_0100600", "PBANKA_0102900"), gset = "pberghei_eg_gene")
#'
#'   ## Change plot domain colors, if desired. Say you have 15 domains in all 3 proteins combined
#'   palette <- randomcoloR::distinctColorPalette(15)
#'   plot + scale_fill_manual(values = palette)
#' }
#'
plotDomain <- function(geneID, mart = "protists_mart", gset = "pfalciparum_eg_gene", input = "ensembl_gene_id", fetchid = "uniprotsptrembl", returnData = FALSE) {
  ensembl_mart <- biomaRt::useEnsemblGenomes(biomart = mart, dataset = gset)
  converted <- biomaRt::getBM(attributes = c(input, fetchid), filters = input, values = geneID, mart = ensembl_mart) %>%
    dplyr::rename(inputid = 1, uniprotid = 2)

  ## track the failed ID
  failed <- setdiff(geneID, converted$inputid)
  if (length(failed) > 0) message(glue::glue("Following genes failed to convert to Uniprot ID and will not be plotted: {failed}\n"))

  plotdata <- drawProteins::get_features(glue::glue_collapse(converted$uniprotid, sep = " ")) %>%
    drawProteins::feature_to_dataframe() %>%
    dplyr::mutate(entryName = converted$inputid[match(accession, converted$uniprotid)])

  if (returnData) {
    return(plotdata)
  }

  # Create the base plot using drawProteins functions
  p <- plotdata %>%
    drawProteins::draw_canvas() %>%
    drawProteins::draw_chains(., plotdata, outline = "white", fill = "#9DCADF") %>%
    drawProteins::draw_domains(., plotdata, label_domains = FALSE) %>%
    drawProteins::draw_repeat(., plotdata) %>%
    drawProteins::draw_motif(., plotdata) %>%
    drawProteins::draw_phospho(., plotdata, size = 8)

  # Now add the theme and scale separately to the plot object
  p <- p +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(plotdata$end), by = round(max(plotdata$end) / 10)),
      limits = c(-1000, max(plotdata$end))
    )

  return(p)
}

#styler:::style_active_file()
