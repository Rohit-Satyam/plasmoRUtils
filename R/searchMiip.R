#' Fetch Protein-Protein interactions from MIIP database
#'
#' This function retrieves Protein-protein interaction data from [MIIP database](http://www.hpppi.iicb.res.in/pfnet/).
#'
#'
#' @import dplyr tidyr
#' @import rvest
#' @import purrr
#' @export
#'
#' @param geneID A character vector of Gene IDs of \emph{Plasmodium falciparum}.
#'
#' @return A data frame of Protein protein interaction provided by MIIP database.
#' @examples
#' \dontrun{
#'  df <- searchMiip(c("PF3D7_0807800","PF3D7_1023900"))
#' }
#'
searchMiip <- function(geneID) {
  baseurl <- "http://www.hpppi.iicb.res.in/pfnet/"
  stages <- c("gam-t1", "mer-t1", "spo-t1", "tro-t1", "rin-t1", "sch-t1")
  #data(pfPlasmodbv68)
  pfPlasmodbv68 <- getTable(org="Plasmodium falciparum 3D7", db="plasmodb")

  pfPlasmodb <- pfPlasmodbv68 %>%
    dplyr::select(`Gene ID`, `Product Description`, `Gene Name or Symbol`, `Previous ID(s)`) %>%
    tidyr::separate_rows(`Previous ID(s)`, sep = ";")
  ## Making subset of pfPlasmodb for query IDs
  inputids <- dplyr::filter(pfPlasmodb, `Gene ID` %in% geneID)

  ## Recursively accessing multiple pages of MIIP to get the PPI data and their stages
  t2 <- purrr::map(stages, function(x) {
    webpage <- rvest::read_html(glue::glue(baseurl, x, ".html"))
    table <- webpage %>%
      rvest::html_nodes("table") %>%
      rvest::html_table(fill = TRUE) %>%
      .[[1]] %>%
      janitor::row_to_names(row_number = 1) %>%
      .[, c(1, 2, 9, 10)]
    colnames(table) <- c("interactorA", "descriptionA", "interactorB", "descriptionB")
    dplyr::bind_rows(
      dplyr::filter(table, interactorA %in% inputids$`Previous ID(s)`),
      dplyr::filter(table, interactorB %in% inputids$`Previous ID(s)`)
    )
  }) %>%
    rlang::set_names(c("gametocyte", "merozoite", "sporozoites", "trophozoites", "ring", "schizont")) %>%
    purrr::keep(~ nrow(.x) > 0)

  t3 <- purrr::map_df(names(t2), ~ dplyr::mutate(t2[[.x]], stage = .x))

  ## Replacing old PF ids with new IDs and description.
  if (nrow(t3) != 0) {
    t4 <- t3 %>%
      dplyr::mutate(
        interactorA = pfPlasmodb$`Gene ID`[match(interactorA, pfPlasmodb$`Previous ID(s)`)],
        descriptionA = pfPlasmodb$`Product Description`[match(interactorA, pfPlasmodb$`Gene ID`)],
        interactorB = pfPlasmodb$`Gene ID`[match(interactorB, pfPlasmodb$`Previous ID(s)`)],
        descriptionB = pfPlasmodb$`Product Description`[match(interactorB, pfPlasmodb$`Gene ID`)]
      )
    return(t4)
  } else {
    message("None of your geneIDs were found in MIIP Database. Kindly visit the following URL for manual search: http://www.hpppi.iicb.res.in/pfnet/stage-map.html ")
  }
}

#styler:::style_active_file()
