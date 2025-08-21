#' List species metadata present in InparanoiDB 9
#'
#' A convenience function to quickly fetch table of all InParanoiDB 9 species and their taxonomy id. 
#'
#' @import utils
#' @export
#'
#' @return A data frame containing information about species present in InParanoiDB9 and their taxon ID.
#' @examples
#' \dontrun{
#' listipdb()
#' }
#'


listipdb <- function(){
  species <- read.csv("https://inparanoidb.sbc.su.se/download/specieslist", check.names = FALSE)
}
